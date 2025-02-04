import xp
import socket
import numpy as np
import quaternion
import time
from dataclasses import dataclass
from UIElements import UIElements

@dataclass
class PoseData:
    plane_orientation: quaternion.quaternion
    local_acceleration: np.ndarray
    modified_acceleration: np.ndarray
    total_acceleration: np.ndarray
    ball_pose: quaternion.quaternion

"""
The relative anglular position, velocity, and acceleration is constant between the aircrafts
frame of reference and that of the pilot. However, the linear velocity/acceleration can differ 
significantly.  

This function provides the transform of linear velocity/acceleration to the pilots frame of reference. 
"""
def pilot_acceleration_correction(velocity, acceleration, angular_velocity, angular_acceleration, pilots_relative_position):
    # v_pilot = v_cg + omega * Δx
    transformed_velocity     = velocity + angular_velocity * pilots_relative_position

    # a_pilot = a_cg + omega * Δx + d(omega)/dt * Δx
    transformed_acceleration = acceleration + (angular_velocity + angular_acceleration) * pilots_relative_position 

    return transformed_velocity, transformed_acceleration
    

def calculate_pose(orientation, velocity, acceleration, on_ground, dt, prev_modified_acceleration, \
                   current_ball_pose, angular_velocity, angular_acceleration, pilots_relative_position):
    # Convert the plane's orientation to the correct frame:
    plane_orientation = quaternion.quaternion(orientation.w, orientation.y, -orientation.z, -orientation.x)
    plane_orientation *= np.sign(plane_orientation.w)

    # Convert the plane's acceleration to the correct frame:
    local_acceleration_xp = quaternion.rotate_vectors(plane_orientation.inverse(), acceleration)
    local_acceleration = np.array([
        -local_acceleration_xp[2],
        -local_acceleration_xp[0],
        local_acceleration_xp[1]
    ])

    # The -z on local_acceleration is limited. This is done because the NOVA would have to flip upside down
    # to handle this, and during the flipping motion the occupant experiences a lot of incorrect intermediate
    # forces. Instead, some x acceleration is bled in to give the sensation of falling / climbing.
    modified_acceleration = np.array([
        local_acceleration[0] + 0.2 * local_acceleration[2],
        local_acceleration[1],
        max(local_acceleration[2], 0)
    ])

    ## TODO: Correction, make sure use correct vels/accels.
    ## velocity, acceleration = pilot_acceleration_correction(velocity, acceleration, angular_velocity, angular_acceleration, pilots_relative_position)
    
    # A low pass filter on the acceleration to eliminate high frequency spikes. While on the ground X-Plane's
    # acceleration values are almost nonsensical, so a much smoothing factor is required.
    if on_ground:
        # smoothing factor of np.power(x, dt) means after 1 second, a ratio of `x` of the old data remains,
        # and a ratio of `1 - x` of the new data has been mixed in.
        smoothing_factor = np.power(0.1, dt)
    else:
        smoothing_factor = np.power(0.005, dt)

    modified_acceleration = prev_modified_acceleration * smoothing_factor + (1 - smoothing_factor) * modified_acceleration

    # Find where gravity is relative to the Ball
    gravity = np.array([0, 9.81, 0])
    local_gravity = quaternion.rotate_vectors(plane_orientation.inverse(), gravity)
    local_gravity = np.array([
        -local_gravity[2],
        -local_gravity[0],
        local_gravity[1]
    ])

    # The "down_target" is the direction in which you want the NOVA occupant to experience
    # a force. We're mostly limited to 1g due to the nature of the NOVA. By aligning gravity
    # with this down vector the occupant will experience 1g acceleration in the opposite direction.
    total_acceleration = modified_acceleration + local_gravity
    down_target = -total_acceleration
    down_target /= np.linalg.norm(down_target)

    # You might be tempted to find the quaternion that takes you directly
    # from global ground to the "down_target", however this results in a singularity when upside down.
    # Instead, the necessary rotation to rotate the current down to the target down is found. Then the
    # previous quaternion is adjusted with this rotation.

    # "down_target" is in a local frame of reference. In order to find where it
    # actually is on the NOVA, we have to rotate it by the current orientation.
    adjusted_down_target = quaternion.rotate_vectors(current_ball_pose, down_target)

    # Now we can rotate this target to point down:
    global_down = np.array([0, 0, -1.0])
    axis = np.cross(adjusted_down_target, global_down)

    rotation = quaternion.quaternion(
        1 + np.dot(adjusted_down_target, global_down),
        axis[0],
        axis[1],
        axis[2]
    )
    rotation = rotation.normalized()

    # This is a rotation in the global frame, so left multiply.
    new_ball_orientation = rotation * current_ball_pose
    new_ball_orientation *= np.sign(new_ball_orientation.w)

    # The global yaw axis can still be rotated around without altering our perception of down.
    # It can be used to induce a sensation of rotation (most important for yawing while upright during taxi)
    # The idea is to calculate what a global yaw would feel like locally, and then dot product that normal
    # vector with the current local angular velocity. This number will be how much global yaw radians / s gets
    # added.
    local_angular_velocity = np.radians([
        angular_velocity[0],
        -angular_velocity[1],
        -angular_velocity[2]
    ])

    global_yaw = np.array([0, 0, 1.0])
    local_yaw = quaternion.rotate_vectors(new_ball_orientation, global_yaw)
    similarity = np.dot(local_yaw, local_angular_velocity)
    yaw_rotation = quaternion.from_rotation_vector(global_yaw * similarity * dt)

    new_ball_orientation = yaw_rotation * new_ball_orientation

    return PoseData(
        plane_orientation=plane_orientation,
        local_acceleration=local_acceleration,
        modified_acceleration=modified_acceleration,
        total_acceleration=total_acceleration,
        ball_pose=new_ball_orientation
    )


class PythonInterface:
    def XPluginStart(self):
        self.Name = "MotionMixing v1.1"
        self.Sig = "motionMixing.demos.xppython3"
        self.Desc = "Eight360's motion plugin for X-Plane."

        self.UI = UIElements()
        self.myMenu = self.UI.createMenu("Quaternion Motion Driver", self.MyMenuHandlerCallback)
        xp.appendMenuItem( self.myMenu, "ON",  'ON')
        xp.appendMenuItem( self.myMenu, "OFF", 'OFF' )
        self.quaternion_toggle = True

        self.address = ("192.168.20.124", 28360)  # Change this to match with the Base's address.

        # For more data refs, see https://www.siminnovations.com/xplane/dataref/
        self.q_dataref = xp.findDataRef("sim/flightmodel/position/q")

        self.local_vx_dataref = xp.findDataRef("sim/flightmodel/position/local_vx")    # m/s
        self.local_vy_dataref = xp.findDataRef("sim/flightmodel/position/local_vy")    # m/s
        self.local_vz_dataref = xp.findDataRef("sim/flightmodel/position/local_vz")    # m/s

        self.local_ax_dataref = xp.findDataRef("sim/flightmodel/position/local_ax")    # m/s^2
        self.local_ay_dataref = xp.findDataRef("sim/flightmodel/position/local_ay")    # m/s^2
        self.local_az_dataref = xp.findDataRef("sim/flightmodel/position/local_az")    # m/s^2

        self.local_roll_rate_dataref  = xp.findDataRef("sim/flightmodel/position/P")   # deg/sec
        self.local_pitch_rate_dataref = xp.findDataRef("sim/flightmodel/position/Q")   # deg/sec
        self.local_yaw_rate_dataref   = xp.findDataRef("sim/flightmodel/position/R")   # deg/sec

        self.local_roll_accel_dataref  = xp.findDataRef("sim/flightmodel/position/P_dot")   # deg/sec^2
        self.local_pitch_accel_dataref = xp.findDataRef("sim/flightmodel/position/Q_dot")   # deg/sec^2
        self.local_yaw_accel_dataref   = xp.findDataRef("sim/flightmodel/position/R_dot")   # deg/sec^2

        self.dx_dataref = xp.findDataRef("sim/graphics/view/pilots_head_x")   # meters
        self.dy_dataref = xp.findDataRef("sim/graphics/view/pilots_head_y")   # meters
        self.dz_dataref = xp.findDataRef("sim/graphics/view/pilots_head_z")   # meters
  
        self.on_ground_dataref = xp.findDataRef("sim/flightmodel/failures/onground_any")

        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.prev_time = None
        self.dt = 1 / 60  # Initial estimate.

        self.pose_data: PoseData | None = None

        window_info = (50, 600, 300, 400, 1,
                      self.DrawWindowCallback,
                      None,
                      None,
                      None,
                      None,
                      0,
                      xp.WindowDecorationRoundRectangle,
                      xp.WindowLayerFloatingWindows,
                      None)
        self.window_id = xp.createWindowEx(window_info)
        self.UI.attachWindowToDestroy(self.window_id, 10)

        xp.registerFlightLoopCallback(self.FlightLoopCallback, 1.0, 0)
        return self.Name, self.Sig, self.Desc

    def XPluginStop(self):
        # Unregister the callback
        xp.unregisterFlightLoopCallback(self.FlightLoopCallback, 0)

        self.sock.close()

    def XPluginEnable(self):
        return 1

    def XPluginDisable(self):
        pass

    def XPluginReceiveMessage(self, inFromWho, inMessage, inParam):
        pass

    def MyMenuHandlerCallback(self, inMenuRef, inItemRef):
        if inItemRef == 'ON':
            self.quaternion_toggle = True
            self.UI.ShowMessage("Quaternion Motion Driver: Activated")
        elif inItemRef == 'OFF':
            self.quaternion_toggle = False
            self.UI.ShowMessage("Quaternion Motion Driver: Deactivated")


    def FlightLoopCallback(self, elapsedMe, elapsedSim, counter, refcon):
        if not self.quaternion_toggle:
            return -1

        current_time = time.monotonic()
        dt = (current_time - self.prev_time) if self.prev_time is not None else (1 / 60)
        self.prev_time = current_time

        temp = []
        xp.getDatavf(self.q_dataref, temp, count=4)
        plane_orientation_xp = quaternion.quaternion(*temp)
        del temp

        velocity_xp = np.array([
            xp.getDataf(self.local_vx_dataref),
            xp.getDataf(self.local_vy_dataref),
            xp.getDataf(self.local_vz_dataref)
        ])

        acceleration_xp = np.array([
            xp.getDataf(self.local_ax_dataref),
            xp.getDataf(self.local_ay_dataref),
            xp.getDataf(self.local_az_dataref)
        ])

        on_ground = xp.getDatai(self.on_ground_dataref)

        angular_velocity_xp = np.array([
            xp.getDataf(self.local_roll_rate_dataref),
            xp.getDataf(self.local_pitch_rate_dataref),
            xp.getDataf(self.local_yaw_rate_dataref)
        ])

        angular_acceleration_xp = np.array([
            xp.getDataf(self.local_roll_accel_dataref),
            xp.getDataf(self.local_pitch_accel_dataref),
            xp.getDataf(self.local_yaw_accel_dataref)
        ])

        pilots_relative_position_xp = np.array([
            xp.getDataf(self.dx_dataref),
            xp.getDataf(self.dy_dataref),
            xp.getDataf(self.dz_dataref)
        ])
        
        prev_modified_acceleration = self.pose_data.modified_acceleration if self.pose_data is not None else np.array([0, 0, 0])
        ball_pose = self.pose_data.ball_pose if self.pose_data is not None else quaternion.quaternion(1, 0, 0, 0)
        self.pose_data = calculate_pose(
            plane_orientation_xp,
            velocity_xp,
            acceleration_xp,
            on_ground,
            dt,
            prev_modified_acceleration,
            ball_pose,
            angular_velocity_xp,
            angular_acceleration_xp,
            pilots_relative_position_xp
        )

        message = f'-1,{self.pose_data.ball_pose.w},{self.pose_data.ball_pose.x},{self.pose_data.ball_pose.y},{self.pose_data.ball_pose.z}\n'

        self.sock.sendto(message.encode(), self.address)

        return -1.0  # -1.0 is a request to be called again next frame.

    def DrawWindowCallback(self, inWindowID, inRefcon):
        # First we get the location of the window passed in to us.
        (left, top, right, bottom) = xp.getWindowGeometry(inWindowID)
        """
        We now use an XPLMGraphics routine to draw a translucent dark
        rectangle that is our window's shape.
        """
        xp.drawTranslucentDarkBox(left, top, right, bottom)
        color = 1.0, 1.0, 1.0

        lines = ["Eight360 X-Plane Plugin"]

        xp.drawString(color, left + 5, top - 20, "Eight360 X-Plane Plugin", 0, xp.Font_Basic)

        if self.pose_data is not None:
            acceleration = self.pose_data.local_acceleration
            lines.append(f"Acceleration: {acceleration[0]:.2f}, {acceleration[1]:.2f}, {acceleration[2]:.2f}")

            modified_acceleration = self.pose_data.modified_acceleration
            lines.append(f"Modified Acceleration: {modified_acceleration[0]:.2f}, {modified_acceleration[1]:.2f}, {modified_acceleration[2]:.2f}")

            total_acceleration = self.pose_data.total_acceleration
            lines.append(f"Total Acceleration: {total_acceleration[0]:.2f}, {total_acceleration[1]:.2f}, {total_acceleration[2]:.2f}")

            plane_orientation = self.pose_data.plane_orientation
            lines.append(f"Plane Orientation: {plane_orientation.w:.2f}, {plane_orientation.x:.2f}, {plane_orientation.y:.2f}, {plane_orientation.z:.2f}")

            sent_quaternion = self.pose_data.ball_pose
            lines.append(f"Sent Orientation: {sent_quaternion.w:.2f}, {sent_quaternion.x:.2f}, {sent_quaternion.y:.2f}, {sent_quaternion.z:.2f}")

        for n, line in enumerate(lines):
            height = top - 20 - n * 10
            xp.drawString(color, left + 5, height, line, 0, xp.Font_Basic)

 