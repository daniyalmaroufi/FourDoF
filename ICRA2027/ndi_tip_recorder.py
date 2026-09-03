#!/usr/bin/env python3
"""
ndi_record_xyz.py -- record the /tip_xyz topic to a CSV, then plot it.

    rosrun your_package ndi_record_xyz.py
    rosrun your_package ndi_record_xyz.py _experiment_name:=curve_plane_calib
    rosrun your_package ndi_record_xyz.py _out:=run1.csv _duration:=30

Stop early with Ctrl-C.  Rows are flushed as they arrive, so an interrupted
run still leaves a complete, readable file.  When recording ends the data is
plotted automatically (disable with _plot:=false).

Parameters:
    ~topic            (str)   topic to record                        ["/tip_xyz"]
    ~experiment_name  (str)   used to auto-name the CSV; ignored if
                              ~out is set; "" falls back to "tip_xyz"  [""]
    ~out              (str)   CSV path; "" auto-names by date        [""]
    ~duration         (float) seconds to record, 0 = until Ctrl-C     [0]
    ~plot             (bool)  show the plot when recording ends      [True]
    ~save_plot        (str)   also write the figure here              [""]
    ~units            (str)   axis label for the plot only            ["m"]
"""

import csv
import datetime
import os
import re
import sys

import rospy
from geometry_msgs.msg import PointStamped


class Recorder(object):
    def __init__(self, path, duration=0.0):
        self.path = path
        self.duration = duration
        self.count = 0
        self.t0 = None
        self.done = False
        self.handle = open(path, "w", newline="")
        self.writer = csv.writer(self.handle)
        self.writer.writerow(["stamp", "elapsed", "x", "y", "z", "frame_id"])
        self.handle.flush()

    def callback(self, msg):
        if self.done:
            return
        stamp = msg.header.stamp.to_sec()
        if self.t0 is None:
            self.t0 = stamp
            rospy.loginfo("first sample received, recording to %s", self.path)
        elapsed = stamp - self.t0

        self.writer.writerow(["%.6f" % stamp, "%.6f" % elapsed,
                              "%.6f" % msg.point.x,
                              "%.6f" % msg.point.y,
                              "%.6f" % msg.point.z,
                              msg.header.frame_id])
        self.handle.flush()          # survive an abrupt Ctrl-C
        self.count += 1

        if self.count % 100 == 0:
            rospy.loginfo("%d samples (%.1f s)", self.count, elapsed)

        if self.duration and elapsed >= self.duration:
            self.done = True
            rospy.loginfo("reached %.1f s, stopping", self.duration)
            rospy.signal_shutdown("duration reached")

    def close(self):
        if self.handle and not self.handle.closed:
            self.handle.close()


def _slugify(name):
    """Make a string safe to use as a filename prefix."""
    name = re.sub(r"\s+", "_", name.strip())
    return re.sub(r"[^A-Za-z0-9_-]", "", name)


def main():
    rospy.init_node("ndi_record_xyz")

    topic = rospy.get_param("~topic", "/tip_xyz")
    experiment_name = rospy.get_param("~experiment_name", "")
    out = rospy.get_param("~out", "")
    duration = float(rospy.get_param("~duration", 0.0))
    do_plot = bool(rospy.get_param("~plot", True))
    save_plot = rospy.get_param("~save_plot", "")
    units = rospy.get_param("~units", "m")

    if not out:
        prefix = _slugify(experiment_name) or "tip_xyz"
        out = datetime.datetime.now().strftime("%s_%%Y%%m%%d_%%H%%M%%S.csv" % prefix)
    out = os.path.abspath(os.path.expanduser(out))

    recorder = Recorder(out, duration)
    rospy.on_shutdown(recorder.close)
    rospy.Subscriber(topic, PointStamped, recorder.callback, queue_size=200)

    if duration:
        rospy.loginfo("recording %s for %.1f s -> %s", topic, duration, out)
    else:
        rospy.loginfo("recording %s until Ctrl-C -> %s", topic, out)
    rospy.loginfo("waiting for the first message on %s ...", topic)

    rospy.spin()
    recorder.close()

    if recorder.count == 0:
        rospy.logwarn("no messages arrived on %s -- is the tracker node "
                      "running and publishing?", topic)
        return 1

    rospy.loginfo("wrote %d samples to %s", recorder.count, out)

    if do_plot:
        try:
            import ndi_plot_xyz
        except ImportError:
            print("\nndi_plot_xyz.py not importable from here; plot with:\n"
                  "    python3 ndi_plot_xyz.py %s" % out)
            return 0
        try:
            ndi_plot_xyz.plot(out, units=units,
                              save=save_plot or None, show=True)
        except Exception as exc:
            print("\nplotting failed (%s); the CSV is fine, plot it with:\n"
                  "    python3 ndi_plot_xyz.py %s" % (exc, out))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except rospy.ROSInterruptException:
        pass