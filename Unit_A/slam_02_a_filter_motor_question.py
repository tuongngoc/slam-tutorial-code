# Implement the first move model for the Lego robot.
# 02_a_filter_motor
# Claus Brenner, 31 OCT 2012
from math import sin, cos, pi, atan2
from pylab import *
from lego_robot import *

# This function takes the old (x, y, heading) pose and the motor ticks
# (ticks_left, ticks_right) and returns the new (x, y, heading).
def filter_step(old_pose, motor_ticks, ticks_to_mm, robot_width):

    # Find out if there is a turn at all.
    if motor_ticks[0] == motor_ticks[1]:
        # No turn. Just drive straight.

        x, y, theta = old_pose
        x = x + motor_ticks[0]*cos(theta)*ticks_to_mm
        y = y + motor_ticks[0]*sin(theta)*ticks_to_mm
        return (x, y, theta)

    else:
        x, y, theta = old_pose
        l, r = motor_ticks
        l = l * ticks_to_mm
        r = r * ticks_to_mm
        # Turn. Compute alpha, R, etc.
        alpha = (r-l)/robot_width
        R = l/alpha
        Rw2 = (R + robot_width/2)
        cx, cy = (x - Rw2*sin(theta),
                  y - Rw2*(-cos(theta)))
        theta = (theta + alpha) % (2*pi)

        x = cx + Rw2*sin(theta)
        y = cy + Rw2*(-cos(theta))
        return (x, y, theta)

if __name__ == '__main__':
    # Empirically derived conversion from ticks to mm.
    ticks_to_mm = 0.349

    # Measured width of the robot (wheel gauge), in mm.
    robot_width = 150.0

    # Read data.
    logfile = LegoLogfile()
    logfile.read("robot4_motors.txt")

    # Start at origin (0,0), looking along x axis (alpha = 0).
    # pose = (0.0, 0.0, 0.0)
    pose = (1850.0, 1897.0, 213.0 / 180.0 * pi)

    # Loop over all motor tick records generate filtered position list.
    filtered = []
    for ticks in logfile.motor_ticks:
        pose = filter_step(pose, ticks, ticks_to_mm, robot_width)
        filtered.append(pose)

    # Draw result.
    for pose in filtered:
        print(pose)
        plot([p[0] for p in filtered], [p[1] for p in filtered], 'bo')
    show()
