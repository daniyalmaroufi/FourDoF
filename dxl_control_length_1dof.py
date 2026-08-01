#!/usr/bin/env python
"""
CLI tool for controlling the Dynamixel-driven linear stage through ROS.
Supports both absolute and relative positioning commands in mm.
Linear Stage: 5mm pitch leadscrew, 819.2 pulses/mm
"""

import readline  # Enables arrow key command history in input()
import rospy
from std_msgs.msg import Float64MultiArray
import sys

# Linear Stage Specifications
PULSES_PER_REVOLUTION = 4096
MM_PER_REVOLUTION = 5.0
CONVERSION_FACTOR = PULSES_PER_REVOLUTION / MM_PER_REVOLUTION  # 819.2 pulses/mm

def pulses_to_mm(pulses):
    """Convert pulses to mm."""
    return pulses / CONVERSION_FACTOR


class MotorControlCLI(object):
    """Interactive CLI for controlling the linear stage via ROS (in mm)."""
    
    def __init__(self):
        rospy.init_node('motor_control_cli', anonymous=True)
        
        # Publisher for motor commands
        self.cmd_pub = rospy.Publisher('/dxl/control_action', Float64MultiArray, queue_size=10)
        
        # Subscribe to motor position feedback
        self.current_position = None  # [mm]
        self.position_sub = rospy.Subscriber('dxl/motor_position', Float64MultiArray, self.position_callback)
        
        # Wait a bit for ROS connections to establish
        rospy.sleep(0.5)
        
        # Default speed for movements [mm/s]
        self.default_speed = 5.0
        
    def position_callback(self, msg):
        """Update current position from encoder feedback (in mm)."""
        if msg.data and len(msg.data) > 0:
            # Position is published directly in mm
            self.current_position = msg.data[0]
    
    def send_relative_command(self, delta_mm, speed_mm_per_sec=None):
        """Send a relative movement command.
        
        Args:
            delta_mm: float, mm to move (positive or negative)
            speed_mm_per_sec: float, mm/second (optional, uses default if not provided)
        """
        if self.current_position is None:
            print("Error: Current position unknown. Waiting for feedback...")
            rospy.sleep(0.5)
            if self.current_position is None:
                print("Error: Still no position feedback. Make sure the motor node is running.")
                return
        
        if speed_mm_per_sec is None:
            speed_mm_per_sec = self.default_speed
        
        # Calculate relative movement
        target_position = self.current_position + delta_mm
        
        print(f"Relative move: {delta_mm:+.2f} mm from {self.current_position:.2f} mm to {target_position:.2f} mm at {speed_mm_per_sec:.2f} mm/s")
        
        msg = Float64MultiArray()
        msg.data = [target_position, speed_mm_per_sec]
        self.cmd_pub.publish(msg)
        print(f"Sent command: goal_position={target_position:.2f} mm, velocity={speed_mm_per_sec:.2f} mm/s")
    
    def send_absolute_command(self, target_mm, speed_mm_per_sec=None):
        """Send command to move to absolute position.
        
        Args:
            target_mm: float, target position in mm
            speed_mm_per_sec: float, mm/second for movement (optional)
        """
        if self.current_position is None:
            print("Error: Current position unknown. Waiting for feedback...")
            rospy.sleep(0.5)
            if self.current_position is None:
                print("Error: Still no position feedback. Make sure the motor node is running.")
                return
        
        if speed_mm_per_sec is None:
            speed_mm_per_sec = self.default_speed
        
        # Take a snapshot of current position to avoid race conditions with async subscriber
        current_pos_snapshot = self.current_position
        
        # Calculate delta
        delta_mm = target_mm - current_pos_snapshot
        
        if abs(delta_mm) < 0.1:  # Already close enough (0.1mm threshold)
            print(f"Already at target position! (current: {current_pos_snapshot:.2f} mm, target: {target_mm:.2f} mm)")
            return
        
        print(f"\nMoving from {current_pos_snapshot:.2f} mm to {target_mm:.2f} mm (delta={delta_mm:+.2f} mm)")
        print(f"Speed: {speed_mm_per_sec:.2f} mm/s")
        print()
        
        # Send command: [goal_position_mm, velocity_mm_s]
        msg = Float64MultiArray()
        msg.data = [target_mm, speed_mm_per_sec]
        self.cmd_pub.publish(msg)
    
    def get_current_position(self):
        """Get and display current linear stage position."""
        if self.current_position is not None:
            print(f"Current position: {self.current_position:.2f} mm")
            return self.current_position
        else:
            print("Position unknown (no feedback received yet)")
            return None
    
    def set_default_speed(self, speed_mm_per_sec):
        """Set the default speed for movements."""
        self.default_speed = float(speed_mm_per_sec)
        print(f"Default speed set to: {self.default_speed:.2f} mm/s")
    
    def stop(self):
        """Stop the motor by sending current position as goal."""
        if self.current_position is None:
            print("Error: Current position unknown. Cannot stop.")
            return
        
        msg = Float64MultiArray()
        msg.data = [self.current_position, 0.1]  # Stay at current position
        self.cmd_pub.publish(msg)
        print(f"Stop command sent (holding at {self.current_position:.2f} mm)")
    
    def print_help(self):
        """Print available commands."""
        help_text = """
╔════════════════════════════════════════════════════════════════════════════╗
║        DYNAMIXEL LINEAR STAGE CONTROL CLI (XL430-W250-T, 5mm pitch)        ║
╠════════════════════════════════════════════════════════════════════════════╣
║                                                                            ║
║  MOVEMENT COMMANDS:                                                        ║
║  ─────────────────                                                         ║
║                                                                            ║
║   r <delta_mm> [speed]   - Relative move (increment from current)          ║
║                            delta_mm:  distance to move (positive/negative) ║
║                            speed:     mm/s (optional, default: 5.0)        ║
║                                                                            ║
║      Example:  r 10 8     - Move forward 10 mm at 8 mm/s                   ║
║      Example:  r -5 3     - Move backward 5 mm at 3 mm/s                   ║
║                                                                            ║
║   a <position_mm> [speed] - Absolute move (move to exact position)         ║
║                            position_mm: target position (e.g., 0-500)      ║
║                            speed:       mm/s (optional, default: 5.0)      ║
║                                                                            ║
║      Example:  a 100 10   - Move to 100 mm at 10 mm/s                      ║
║      Example:  a 0 5      - Return to 0 mm at 5 mm/s                       ║
║                                                                            ║
╠════════════════════════════════════════════════════════════════════════════╣
║  UTILITY COMMANDS:                                                         ║
║  ──────────────────                                                        ║
║                                                                            ║
║   p                  - Print current position (mm)                         ║
║   speed <mm_per_sec> - Set default speed for all commands                  ║
║   s or stop          - Stop (hold at current position)                     ║
║   h or help          - Show this help                                      ║
║   q or quit or exit  - Exit CLI                                            ║
║                                                                            ║
╠════════════════════════════════════════════════════════════════════════════╣
║  LINEAR STAGE SPECIFICATIONS:                                              ║
║  ─────────────────────────────                                             ║
║                                                                            ║
║   Hardware:    Dynamixel XL430-W250-T with 5mm pitch leadscrew             ║
║   Resolution:  819.2 pulses/mm (4096 pulses/rev ÷ 5mm/rev)                 ║
║   Range:       Typically -500 to +500 mm (configurable)                    ║
║                                                                            ║
║   Max Speed:   ~500 mm/s (depends on voltage and load)                     ║
╚════════════════════════════════════════════════════════════════════════════╝

⚠  Make sure dxl_control_length_node is running before using this CLI!
"""
        print(help_text)
    
    def run(self):
        """Main interactive loop."""
        print("\n" + "="*70)
        print("  DYNAMIXEL LINEAR STAGE CONTROL CLI (5mm pitch)")
        print("="*70)
        print("\nResolution: 819.2 pulses/mm (4096 pulses/rev, 5mm/rev)")
        print("Type 'h' for help, 'q' to quit\n")
        
        # Wait for initial position
        print("Waiting for linear stage position feedback...", end="", flush=True)
        timeout = rospy.Time.now() + rospy.Duration(3.0)
        while self.current_position is None and rospy.Time.now() < timeout:
            rospy.sleep(0.1)
            print(".", end="", flush=True)
        print()
        
        if self.current_position is not None:
            print(f"✓ Connected! Current position: {self.current_position:.2f} mm\n")
        else:
            print("⚠ Warning: No position feedback. Make sure dxl_control_length_node is running!\n")
        
        while not rospy.is_shutdown():
            try:
                # Get user input
                user_input = input("stage> ").strip()
                
                if not user_input:
                    continue
                
                # Parse command
                parts = user_input.split()
                cmd = parts[0].lower()
                
                # Handle commands
                if cmd in ['q', 'quit', 'exit']:
                    print("Exiting...")
                    break
                
                elif cmd in ['h', 'help']:
                    self.print_help()
                
                elif cmd == 'p':
                    self.get_current_position()
                
                elif cmd == 's' or cmd == 'stop':
                    self.stop()
                
                elif cmd == 'speed':
                    if len(parts) < 2:
                        print("Error: Missing speed value. Usage: speed <value>")
                    else:
                        try:
                            speed = float(parts[1])
                            if speed <= 0:
                                print("Error: Speed must be positive")
                            else:
                                self.set_default_speed(speed)
                        except ValueError:
                            print("Error: Speed must be a number")
                
                elif cmd == 'r':
                    # Relative command
                    if len(parts) < 2:
                        print("Error: Missing delta. Usage: r <delta_mm> [speed]")
                    else:
                        try:
                            delta = float(parts[1])
                            speed = float(parts[2]) if len(parts) > 2 else None
                            self.send_relative_command(delta, speed)
                        except ValueError:
                            print("Error: Invalid number format")
                
                elif cmd == 'a':
                    # Absolute command
                    if len(parts) < 2:
                        print("Error: Missing target position. Usage: a <position_mm> [speed]")
                    else:
                        try:
                            target = float(parts[1])
                            speed = float(parts[2]) if len(parts) > 2 else None
                            self.send_absolute_command(target, speed)
                        except ValueError:
                            print("Error: Invalid number format")
                
                else:
                    print(f"Unknown command: '{cmd}'. Type 'h' for help.")
            
            except KeyboardInterrupt:
                print("\n\nInterrupted. Type 'q' to quit or continue entering commands.")
            except EOFError:
                print("\nExiting...")
                break
            except Exception as e:
                print(f"Error: {e}")
        
        print("CLI terminated.")


def main():
    try:
        cli = MotorControlCLI()
        cli.run()
    except rospy.ROSInterruptException:
        print("ROS interrupted")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()