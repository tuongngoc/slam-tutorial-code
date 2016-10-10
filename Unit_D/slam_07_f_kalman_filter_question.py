# The full Kalman filter, consisting of prediction and correction step.
#
# slam_07_f_kalman_filter
# Claus Brenner, 12.12.2012
from lego_robot import *
from math import sin, cos, pi, atan2, sqrt
from numpy import *
from slam_d_library import get_observations, write_cylinders


class ExtendedKalmanFilter:
    def __init__(self, state, covariance,
                 robot_width, scanner_displacement,
                 control_motion_factor, control_turn_factor,
                 measurement_distance_stddev, measurement_angle_stddev):
        # The state. This is the core data of the Kalman filter.
        self.state = state
        self.covariance = covariance

        # Some constants.
        self.robot_width = robot_width
        self.scanner_displacement = scanner_displacement
        self.control_motion_factor = control_motion_factor
        self.control_turn_factor = control_turn_factor
        self.measurement_distance_stddev = measurement_distance_stddev
        self.measurement_angle_stddev = measurement_angle_stddev
        self.Q = diag([measurement_distance_stddev**2,
                       measurement_angle_stddev**2])
    @staticmethod
    def g(state, control, w):
        x, y, theta = state
        l, r = control
        if r != l:
            alpha = (r - l) / w
            rad = l/alpha
            g1 = x + (rad + w/2.)*(sin(theta+alpha) - sin(theta))
            g2 = y + (rad + w/2.)*(-cos(theta+alpha) + cos(theta))
            g3 = (theta + alpha + pi) % (2*pi) - pi
        else:
            g1 = x + l * cos(theta)
            g2 = y + l * sin(theta)
            g3 = theta

        return array([g1, g2, g3])

    @staticmethod
    def dg_dstate(state, control, w):
        theta = state[2]
        l, r = control
        if r != l:

            # This is for the case r != l.
            # g has 3 components and the state has 3 components, so the
            # derivative of g with respect to all state variables is a
            # 3x3 matrix.
            alpha = (r - l) / w
            rad = l/alpha
            Rw2 = rad + w/2
            m = array([ [1, 0, Rw2*(cos(theta + alpha) - cos(theta))],
                        [0, 1, Rw2*(sin(theta + alpha) - sin(theta))],
                        [0, 0, 1]
                ])

        else:
            # This is for the special case r == l.
            m = array([ [1, 0, -l*sin(theta)],
                        [0, 1,  l*cos(theta)],
                        [0, 0, 1]
                ])
        return m

    @staticmethod
    def dg_dcontrol(state, control, w):
        theta = state[2]
        l, r = tuple(control)
        if r != l:
            alpha = (r - l) / w
            rad = l/alpha
            # This is for the case l != r.
            # Note g has 3 components and control has 2, so the result
            # will be a 3x2 (rows x columns) matrix.
            dg1dl = (w*r)*( sin(theta+alpha) - sin(theta))/(r-l)**2 \
                   -(r+l)*cos(theta+alpha)/(2*(r-l))
            dg2dl = (w*r)*(-cos(theta+alpha) + cos(theta))/(r-l)**2 \
                   -(r+l)*sin(theta+alpha)/(2*(r-l))


            dg1dr = (-w*l)*( sin(theta+alpha) - sin(theta))/(r-l)**2 \
                   +(r+l)*cos(theta+alpha)/(2*(r-l))
            dg2dr = (-w*l)*(-cos(theta+alpha) + cos(theta))/(r-l)**2 \
                   +(r+l)*sin(theta+alpha)/(2*(r-l))
        else:
            # This is for the special case l == r.
            dg1dl = 0.5*(cos(theta) + (l/w)*sin(theta))
            dg2dl = 0.5*(sin(theta) - (l/w)*cos(theta))
            dg1dr = 0.5*(-(l/w)*sin(theta) + cos(theta))
            dg2dr = 0.5*( (l/w)*cos(theta) + sin(theta))

        dg3dl = -1/w
        dg3dr =  1/w

        m = array([[dg1dl, dg1dr],
                   [dg2dl, dg2dr],
                   [dg3dl, dg3dr]])

        return m

    @staticmethod
    def get_error_ellipse(covariance):
        """Return the position covariance (which is the upper 2x2 submatrix)
           as a triple: (main_axis_angle, stddev_1, stddev_2), where
           main_axis_angle is the angle (pointing direction) of the main axis,
           along which the standard deviation is stddev_1, and stddev_2 is the
           standard deviation along the other (orthogonal) axis."""
        eigenvals, eigenvects = linalg.eig(covariance[0:2,0:2])
        angle = atan2(eigenvects[1,0], eigenvects[0,0])
        return (angle, sqrt(eigenvals[0]), sqrt(eigenvals[1]))

    def predict(self, control):
        """The prediction step of the Kalman filter."""
        # covariance' = G * covariance * GT + R
        # where R = V * (covariance in control space) * VT.
        # Covariance in control space depends on move distance.
        l, r = control

        sigL2 = (self.control_motion_factor*l)**2 \
                + (self.control_turn_factor*(r-l))**2
        sigR2 = (self.control_motion_factor*r)**2 \
                + (self.control_turn_factor*(r-l))**2

        control_cov = diag([sigL2, sigR2])
        V = self.dg_dcontrol(self.state, control, self.robot_width)
        G = self.dg_dstate(self.state, control, self.robot_width)
        self.covariance = G @ self.covariance @ G.T + V @ control_cov @ V.T

        # state' = g(state, control)
        self.state = self.g(self.state, control, self.robot_width)

    @staticmethod
    def h(state, landmark, scanner_displacement):
        """Takes a (x, y, theta) state and a (x, y) landmark, and returns the
           measurement (range, bearing)."""
        dx = landmark[0] - (state[0] + scanner_displacement * cos(state[2]))
        dy = landmark[1] - (state[1] + scanner_displacement * sin(state[2]))
        r = sqrt(dx * dx + dy * dy)
        alpha = (atan2(dy, dx) - state[2] + pi) % (2*pi) - pi

        return array([r, alpha])

    @staticmethod
    def dh_dstate(state, landmark, scanner_displacement):
        x, y, theta = state
        xm, ym = landmark
        d = scanner_displacement

        xl = x + d*cos(theta)
        yl = y + d*sin(theta)
        q = (xm - xl)**2 + (ym - yl)**2

        dx = xm - xl
        dy = ym - yl

        # Note that:
        # x y theta is state[0] state[1] state[2]
        # x_m y_m is landmark[0] landmark[1]
        # The Jacobian of h is a 2x3 matrix.
        drdx = -dx/sqrt(q)
        drdy = -dy/sqrt(q)
        drdtheta = d*(dx*sin(theta) - dy*cos(theta))/sqrt(q)

        dadx = dy/q
        dady = -dx/q
        dadtheta = - d/q * (dx*cos(theta) + dy*sin(theta)) - 1
        return array([[drdx, drdy, drdtheta],
                      [dadx, dady, dadtheta]])

    def correct(self, measurement, landmark):
        """The correction step of the Kalman filter."""
        Z = array(measurement)
        Z[1] = atan2(sin(Z[1]), cos(Z[1]))

        H = self.dh_dstate(self.state, landmark, self.scanner_displacement)
        P = self.covariance

        K = P @ H.T @ linalg.inv(H @ P @ H.T + self.Q)
        pred_z = self.h(self.state, landmark, self.scanner_displacement)
        innov = Z - pred_z
        innov[1] = atan2(sin(innov[1]), cos(innov[1]))

        self.state = self.state + K @ innov
        self.covariance = (eye(3) - K @ H) @ P

if __name__ == '__main__':
    # Robot constants.
    scanner_displacement = 30.0
    ticks_to_mm = 0.349
    robot_width = 155.0

    # Cylinder extraction and matching constants.
    minimum_valid_distance = 20.0
    depth_jump = 100.0
    cylinder_offset = 90.0
    max_cylinder_distance = 300.0

    # Filter constants.
    control_motion_factor = 0.35  # Error in motor control.
    control_turn_factor = 0.6  # Additional error due to slip when turning.
    measurement_distance_stddev = 200.0  # Distance measurement error of cylinders.
    measurement_angle_stddev = 15.0 / 180.0 * pi  # Angle measurement error.

    # Measured start position.
    initial_state = array([1850.0, 1897.0, 213.0 / 180.0 * pi])
    # Covariance at start position.
    initial_covariance = diag([100.0**2, 100.0**2, (10.0 / 180.0 * pi) ** 2])
    # Setup filter.
    kf = ExtendedKalmanFilter(initial_state, initial_covariance,
                              robot_width, scanner_displacement,
                              control_motion_factor, control_turn_factor,
                              measurement_distance_stddev,
                              measurement_angle_stddev)

    # Read data.
    logfile = LegoLogfile()
    logfile.read("robot4_motors.txt")
    logfile.read("robot4_scan.txt")
    logfile.read("robot_arena_landmarks.txt")
    reference_cylinders = [l[1:3] for l in logfile.landmarks]

    # Loop over all motor tick records and all measurements and generate
    # filtered positions and covariances.
    # This is the Kalman filter loop, with prediction and correction.
    states = []
    covariances = []
    matched_ref_cylinders = []
    for i in range(len(logfile.motor_ticks)):
        # Prediction.
        control = array(logfile.motor_ticks[i]) * ticks_to_mm
        kf.predict(control)

        # Correction.
        observations = get_observations(
            logfile.scan_data[i],
            depth_jump, minimum_valid_distance, cylinder_offset,
            kf.state, scanner_displacement,
            reference_cylinders, max_cylinder_distance)
        for j in range(len(observations)):
            kf.correct(*observations[j])

        # Log state, covariance, and matched cylinders for later output.
        states.append(kf.state)
        covariances.append(kf.covariance)
        matched_ref_cylinders.append([m[1] for m in observations])

    # Write all states, all state covariances, and matched cylinders to file.
    f = open("kalman_prediction_and_correction.txt", "w")
    for i in range(len(states)):
        # Output the center of the scanner, not the center of the robot.
        print("F %f %f %f" % \
            tuple(states[i] + [scanner_displacement * cos(states[i][2]),
                               scanner_displacement * sin(states[i][2]),
                               0.0]), file=f)
        # Convert covariance matrix to angle stddev1 stddev2 stddev-heading form
        e = ExtendedKalmanFilter.get_error_ellipse(covariances[i])
        print("E %f %f %f %f" % (e + (sqrt(covariances[i][2,2]),)), file=f)
        # Also, write matched cylinders.
        write_cylinders(f, "W C", matched_ref_cylinders[i])

    f.close()
