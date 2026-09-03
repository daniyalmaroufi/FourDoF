#!/usr/bin/env python3
"""
ndi_tip_xyz.py -- publish the 3D position of a single stray marker from an
Ethernet NDI tracker (Polaris Vega / Lyra) on the topic `tip_xyz`.

Talks the ASCII Combined API directly over TCP, so it needs no NDI SDK and no
cisst/ROS driver stack.  Intended for the case where exactly one retro-reflective marker is
in the measurement volume.

A passive tool definition (.rom) must still be supplied via ~rom: the Polaris
only illuminates retro-reflective markers while a passive port handle is
enabled.  Any .rom will do -- it is a dummy and will read MISSING.

Published topic:
    tip_xyz   (geometry_msgs/PointStamped)

Parameters:
    ~ip           (str)   tracker IP address                    [required]
    ~port         (int)   API port                              [8765]
    ~rate         (float) polling rate in Hz                    [60.0]
    ~frame_id     (str)   frame_id stamped on the message       ["ndi"]
    ~units        (str)   "m" (ROS convention) or "mm"          ["m"]
    ~skip_out_of_volume (bool) drop markers flagged OOV         [True]
    ~timeout      (float) socket timeout in seconds             [2.0]
    ~rom          (str)   path to a passive tool .rom, loaded as
                          a dummy so the illuminators fire        [""]
    ~vsel         (int)   measurement volume to select, 0 = leave
                          the default (see SFLIST 03)             [0]
"""

import math
import socket

import rospy
from geometry_msgs.msg import PointStamped


class NDIClient(object):
    """Minimal ASCII-API client: connect, start tracking, poll stray markers."""

    def __init__(self, host, port=8765, timeout=2.0):
        self.host = host
        self.port = port
        self.timeout = timeout
        self.sock = None
        self.buf = b""

    # -- connection ------------------------------------------------------

    def connect(self, rom_path=None, vsel=0):
        self.close()
        self.buf = b""
        self.sock = socket.create_connection((self.host, self.port), self.timeout)
        self.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        self._expect_okay("INIT ")
        if vsel:
            self._expect_okay("VSEL %d" % vsel)
        if rom_path:
            self.load_passive_tool(rom_path)
        self._expect_okay("TSTART 80")

    def load_passive_tool(self, rom_path):
        """Allocate, load and enable one passive port handle.

        The Polaris only fires its illuminators for retro-reflective markers
        when a passive port handle is enabled.  Without this, the sensor sees
        nothing and TX reports zero stray markers.  The tool itself will read
        MISSING unless its actual geometry is in view -- that is expected and
        harmless; it exists only to turn the illuminators on.
        """
        with open(rom_path, "rb") as handle:
            rom = handle.read()
        if len(rom) > 960:
            raise IOError("%s is %d bytes; the limit is 960" % (rom_path, len(rom)))

        # allocate a port handle for a passive (wireless) tool
        reply = self._cmd("PHRQ *********1****")
        port = reply[0:2]
        if not port.isalnum():
            raise IOError("PHRQ returned %r" % reply)
        rospy.loginfo("allocated passive port handle %s", port)

        # write the definition in 64-byte chunks, padded to a multiple of 64
        padded = rom + b"\x00" * (-len(rom) % 64)
        for chunk in range(len(padded) // 64):
            block = padded[chunk * 64:(chunk + 1) * 64]
            self._expect_okay("PVWR %s%04X%s"
                              % (port, chunk * 64, block.hex().upper()))

        self._expect_okay("PINIT %s" % port)
        self._expect_okay("PENA %sD" % port)   # D = dynamic
        rospy.loginfo("loaded and enabled %s -- illuminators active", rom_path)
        return port

    def close(self):
        if self.sock is None:
            return
        try:
            self._expect_okay("TSTOP ")
        except Exception:
            pass
        try:
            self.sock.close()
        except Exception:
            pass
        self.sock = None

    # -- protocol --------------------------------------------------------

    def _cmd(self, text):
        """Send one ASCII command and return its reply, minus the trailing CR.

        The space delimiter form makes the CRC optional, so none is appended.
        """
        self.sock.sendall(text.encode("ascii") + b"\r")
        while b"\r" not in self.buf:
            chunk = self.sock.recv(4096)
            if not chunk:
                raise IOError("tracker closed the connection")
            self.buf += chunk
        reply, _, self.buf = self.buf.partition(b"\r")
        return reply.decode("ascii")

    def _expect_okay(self, text):
        """Accept OKAY, and WARNING (e.g. PINIT returns WARNING05 on firmware
        006+ when it picks a default marker wavelength -- harmless)."""
        reply = self._cmd(text)
        if reply.startswith("OKAY"):
            return reply
        if reply.startswith("WARNING"):
            rospy.logwarn("%s returned %s", text.strip(), reply[:9])
            return reply
        raise IOError("command %r returned %r" % (text.strip(), reply))

    def stray_markers(self):
        """Poll once.  Returns a list of (x_mm, y_mm, z_mm, out_of_volume)."""
        # TX reply option 1801 = transforms (0001) + report OOV (0800)
        #                        + 3D positions of stray passive markers (1000)
        r = self._cmd("TX 1801")
        i = 0
        n_handles = int(r[0:2], 16)
        i += 2
        for _ in range(n_handles):           # skip past any defined tools
            i = r.index("\n", i) + 1
        n = int(r[i:i + 2], 16)
        i += 2
        oov_chars = int(math.ceil(n / 4.0))  # 4 markers packed per hex char
        oov_bits = int(r[i:i + oov_chars], 16) if oov_chars else 0
        i += oov_chars
        markers = []
        for k in range(n):
            # each coordinate is sign + 6 digits with an implied decimal
            # point two places in, i.e. "+001234" -> 12.34 mm
            x, y, z = (int(r[i + j * 7:i + (j + 1) * 7]) / 100.0 for j in range(3))
            i += 21
            # the bitfield is padded at the front, so marker 0 is the most
            # significant of the n used bits, not the least
            markers.append((x, y, z, bool(oov_bits >> (n - 1 - k) & 1)))
        return markers


def pick_marker(markers, last_xyz):
    """Choose which marker is 'the' marker when more than one is reported.

    Reflections and stray IR sources show up as extra entries, and NDI does not
    label markers across frames, so index 0 is not guaranteed to be yours.
    Nearest-to-last-known is a cheap and usually sufficient association.
    """
    if len(markers) == 1:
        return markers[0]
    if last_xyz is None:
        rospy.logwarn_throttle(
            5.0, "%d markers visible and no previous position to match "
                 "against; using the first one" % len(markers))
        return markers[0]
    rospy.logwarn_throttle(
        5.0, "%d markers visible (reflections?); tracking the one nearest to "
             "the last known position" % len(markers))

    def dist2(m):
        return sum((m[j] - last_xyz[j]) ** 2 for j in range(3))

    return min(markers, key=dist2)


def main():
    rospy.init_node("ndi_tip_xyz")

    ip = rospy.get_param("~ip", "172.31.1.10")
    if not ip:
        rospy.logfatal("parameter ~ip is required (the tracker's IP address)")
        return
    port = int(rospy.get_param("~port", 8765))
    rate_hz = float(rospy.get_param("~rate", 60.0))
    frame_id = rospy.get_param("~frame_id", "ndi")
    units = rospy.get_param("~units", "m")
    skip_oov = bool(rospy.get_param("~skip_out_of_volume", True))
    timeout = float(rospy.get_param("~timeout", 2.0))
    rom = rospy.get_param("~rom", "/home/daniyal/ctsdr_ndi_ws/src/cisst-saw/sawNDITracker/core/share/roms/8700338.rom")
    vsel = int(rospy.get_param("~vsel", 0))

    if units not in ("m", "mm"):
        rospy.logfatal("~units must be 'm' or 'mm', got %r" % units)
        return
    scale = 0.001 if units == "m" else 1.0

    if not rom:
        rospy.logwarn("no ~rom given; on Polaris Vega the illuminators stay "
                      "off without an enabled passive port handle, so no "
                      "passive marker will ever be reported")

    pub = rospy.Publisher("tip_xyz", PointStamped, queue_size=10)
    client = NDIClient(ip, port, timeout)
    rate = rospy.Rate(rate_hz)
    last_xyz = None

    rospy.on_shutdown(client.close)

    while not rospy.is_shutdown():
        # (re)connect, retrying until the node is asked to shut down
        if client.sock is None:
            try:
                rospy.loginfo("connecting to NDI tracker at %s:%d", ip, port)
                client.connect(rom or None, vsel)
                rospy.loginfo("connected, tracking started")
            except Exception as exc:
                rospy.logerr_throttle(5.0, "connect failed: %s" % exc)
                rospy.sleep(2.0)
                continue

        try:
            markers = client.stray_markers()
        except Exception as exc:
            rospy.logerr("lost connection to tracker: %s", exc)
            client.sock = None
            continue

        n_seen = len(markers)
        if skip_oov:
            markers = [m for m in markers if not m[3]]

        if not markers:
            if n_seen:
                rospy.logwarn_throttle(
                    2.0, "%d marker(s) detected but all flagged out of volume; "
                         "aim the sensor or set _skip_out_of_volume:=false"
                         % n_seen)
            else:
                rospy.logwarn_throttle(2.0, "no marker in view")
            rate.sleep()
            continue

        x, y, z, _ = pick_marker(markers, last_xyz)
        last_xyz = (x, y, z)

        msg = PointStamped()
        msg.header.stamp = rospy.Time.now()
        msg.header.frame_id = frame_id
        msg.point.x = x * scale
        msg.point.y = y * scale
        msg.point.z = z * scale
        pub.publish(msg)

        rate.sleep()


if __name__ == "__main__":
    try:
        main()
    except rospy.ROSInterruptException:
        pass