# For each cylinder in the scan, find its cartesian coordinates,
# in the scanner's coordinate system.
# Write the result to a file which contains all cylinders, for all scans.
# 03_d_find_cylinders_cartesian
# Claus Brenner, 09 NOV 2012
from lego_robot import *
from math import sin, cos

# Find the derivative in scan data, ignoring invalid measurements.
def compute_derivative(scan, min_dist):
    jumps = [ 0 ]
    for i in range(1, len(scan) - 1):
        l = scan[i-1]
        r = scan[i+1]
        if l > min_dist and r > min_dist:
            derivative = (r - l) / 2.0
            jumps.append(derivative)
        else:
            jumps.append(0)
    jumps.append(0)
    return jumps

# For each area between a left falling edge and a right rising edge,
# determine the average ray number and the average depth.
def find_cylinders(scan, scan_derivative, jump, min_dist):
    cylinder_list = []
    on_cylinder = False
    sum_ray, sum_depth, rays = 0.0, 0.0, 0

    for i in range(len(scan_derivative)):
        df = scan_derivative[i]
        if df <= -jump:
            on_cylinder = True
            sum_ray, sum_depth, rays = 0.0, 0.0, 0

        if on_cylinder:
            if abs(df) <= jump: # Practically constant
                sum_ray += i
                sum_depth += scan[i]
                rays += 1
            elif df >= jump:
                on_cylinder = False
                cylinder_list.append( (sum_ray/rays, sum_depth/rays) )


    return cylinder_list

def compute_cartesian_coordinates(cylinders, cylinder_offset):
    result = []
    for c in cylinders:
        ray, depth = c
        depth = depth + cylinder_offset

        ray_angle = LegoLogfile.beam_index_to_angle(ray)
        x = depth*cos(ray_angle)
        y = depth*sin(ray_angle)
        result.append( (x, y) ) # Replace this by your (x,y)
    return result


if __name__ == '__main__':

    minimum_valid_distance = 20.0
    depth_jump = 100.0
    cylinder_offset = 90.0

    # Read the logfile which contains all scans.
    logfile = LegoLogfile()
    logfile.read("robot4_scan.txt")

    # Write a result file containing all cylinder records.
    # Format is: D C x[in mm] y[in mm] ...
    # With zero or more points.
    # Note "D C" is also written for otherwise empty lines (no
    # cylinders in scan)
    out_file = open("cylinders.txt", "w")
    for scan in logfile.scan_data:
        # Find cylinders.
        der = compute_derivative(scan, minimum_valid_distance)
        cylinders = find_cylinders(scan, der, depth_jump,
                                   minimum_valid_distance)
        cartesian_cylinders = compute_cartesian_coordinates(cylinders,
                                                            cylinder_offset)
        # Write to file.
        output = 'D C '+' '.join("%.1f %.1f" % c for c in cartesian_cylinders)+'\n'
        out_file.write(output)

    out_file.close()
