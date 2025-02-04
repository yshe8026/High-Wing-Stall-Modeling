from XPPython3 import xp
import numpy as np
import initalisePlugin as ip
from Autopilot import LinearAutopilots
import os

from enum import Enum

class ThrottleControls(Enum):
    EMPTY = 0
    IndicatedAirspeedHold = 1

class VerticalControls(Enum):
    EMPTY = 0
    AltitudeHold = 1
    VerticalSpeedHold = 2

class LateralControls(Enum):
    EMPTY = 0
    HeadingHold = 1
    WingLeveler = 3

class AltitudeModes(Enum):
    Pitch = 3 
    VerticalSpeed = 4
    LevelChange = 5
    AltitudeHold = 6
    TerrainFollow = 7
    Glideslope = 8
    VNAV = 9

class XpCOMM(str):
    HEADING_HOLD                = "sim/autopilot/heading_hold"
    NAVIGATION                  = "sim/autopilot/NAV"
    ALT_HOLD                    = "sim/autopilot/altitude_hold"
    ALT_ARM                     = "sim/autopilot/altitude_arm"
    VERTICAL_SPEED              = "sim/autopilot/vertical_speed"
    AUTO_THROTTLE_ON            = "sim/autopilot/autothrottle_on"
    AUTO_THROTTLE_OFF           = "sim/autopilot/autothrottle_off"
    AUTO_THROTTLE_TOGGLE        = "sim/autopilot/autothrottle_toggle"
    WING_LEVELER                = "sim/autopilot/wing_leveler"
    FLIGHT_DIRECTOR_ON          = "sim/autopilot/fdir_on"
    FLIGHT_DIRECTOR_OFF         = "sim/autopilot/fdir_off"

    def register(self):

        self.HEADING_HOLD          = xp.findCommand(self.HEADING_HOLD)
        self.NAVIGATION            = xp.findCommand(self.NAVIGATION)
        self.ALT_HOLD              = xp.findCommand(self.ALT_HOLD)
        self.ALT_ARM               = xp.findCommand(self.ALT_ARM)
        self.VERTICAL_SPEED        = xp.findCommand(self.VERTICAL_SPEED)
        self.AUTO_THROTTLE_ON      = xp.findCommand(self.AUTO_THROTTLE_ON)
        self.AUTO_THROTTLE_OFF     = xp.findCommand(self.AUTO_THROTTLE_OFF)
        self.AUTO_THROTTLE_TOGGLE  = xp.findCommand(self.AUTO_THROTTLE_TOGGLE)
        self.WING_LEVELER          = xp.findCommand(self.WING_LEVELER)
        self.FLIGHT_DIRECTOR_ON    = xp.findCommand(self.FLIGHT_DIRECTOR_ON)
        self.FLIGHT_DIRECTOR_OFF   = xp.findCommand(self.FLIGHT_DIRECTOR_OFF)

class PythonInterface:

    def __init__(self):
        self.Name = "Custom autopilot"
        self.Sig  = "flightsim.Python.custom_autopilot"
        self.Desc = "An autopilot plugin"

        self.state_DREF                 = "customPhysics/X_state" # READ ONLY
        self.control_DREF               = "customPhysics/U_state" # Read and write
        self.physics_DREF               = "customPhysics/physics_status" # Gets the status of the custom physics
        self.heading_mag_DREF           = "sim/flightmodel2/position/mag_psi" # Dref for magnetic heading
        self.altitude_target_DREF       = "sim/cockpit/autopilot/altitude" # Dref for storing altitude value 
        self.vertical_speed_target_DREF = "sim/cockpit/autopilot/vertical_velocity" #Dref for storing vertical speed values 
        self.airspeed_target_DREF       = "sim/cockpit/autopilot/airspeed" # Dref for airspeed 
        self.APisMach_DREF              = "sim/cockpit/autopilot/airspeed_is_mach" # boolean to change airspeed to mach 
        self.heading_target_DREF        = "sim/cockpit/autopilot/heading_mag" # Dref for magnetic heading 
        self.flight_director_DREF       = "sim/cockpit2/autopilot/flight_director_mode" # 0 is off, 1 is on (Flight director)
        self.roll_target_DREF           = "sim/cockpit/autopilot/flight_director_roll" # Dref for roll target
        self.HDG_yaw_damper_DREF        = "sim/cockpit/switches/yaw_damper_on" # Dref for yaw damper
        self.autopilot_state_DREF       = "sim/cockpit/autopilot/autopilot_state"
        self.auto_throttle_DREF         = "sim/cockpit2/autopilot/autothrottle_enabled" # Dref for auto throttle
        self.autopilot_ai_override_DREF = "sim/operation/override/override_plane_ai_autopilot"
        self.autopilot_override_DREF    = "sim/operation/override/override_autopilot"

        ## Panel mode setting DREFs 
        self.altitide_mode_DREF         = 'sim/cockpit2/autopilot/altitude_mode'

        self.APState                = 0
        self.heading_mag            = 0
        self.curr_altitude          = 0
        self.curr_IAS               = 0
        self.altitude_target        = 0
        self.vertical_speed_target  = 0
        self.airspeed_target        = 0
        self.APisMach               = 0
        self.heading_target         = 0
        self.roll_target            = 0
        self.auto_throttle          = 0

        self.HDG_yaw_damper         = 0

        self.physics_status         = False

        self.active_throttle_contrl = ThrottleControls.EMPTY
        self.active_vertical_contrl = VerticalControls.EMPTY
        self.active_lateral_contrl  = LateralControls.EMPTY

        # Boolean to check if autopilot is engaged. It is true if 
        # if any of the autopilot modes are on (Not EMPTY)
        self.autopilot_engaged = False 

        config_path = os.path.join('Resources', 'plugins', 'PythonPlugins', 'config.ini')
        user_defined = ip.InitalisePlugin(config_path)
        autopilot_settings = user_defined.autopilot_settings()
        controllers = user_defined.load_controller()
        if controllers is None:
            raise Exception('Failed to load controller')
        cntrl_var = controllers.get_controller_variables()
        self.controller = LinearAutopilots(cntrl_var, autopilot_settings)
        # Logging 
        self.log = user_defined.logging_on()

        self.XpComm = XpCOMM()

        self.target_error = {
            'phi_c': -1,
            'psi_c': -1,
            'u_c': -1,
            'vs_c': -1,
            'alt_c': -1,
        }

        self.X = np.zeros(12)
        self.U = np.zeros(5)
        self.U0 = np.zeros(5)

        self.autopilot_update_rate = 0.1 # seconds
        self.dt = self.autopilot_update_rate

        self.autopilot_loop = None


    def XPluginStart(self):

        # Register the DREFs
        self.state_DREF                 = xp.findDataRef(self.state_DREF)
        self.control_DREF               = xp.findDataRef(self.control_DREF)
        self.physics_DREF               = xp.findDataRef(self.physics_DREF)
        self.heading_mag_DREF           = xp.findDataRef(self.heading_mag_DREF)
        self.altitude_target_DREF       = xp.findDataRef(self.altitude_target_DREF)
        self.vertical_speed_target_DREF = xp.findDataRef(self.vertical_speed_target_DREF)
        self.airspeed_target_DREF       = xp.findDataRef(self.airspeed_target_DREF)
        self.APisMach_DREF              = xp.findDataRef(self.APisMach_DREF)
        self.heading_target_DREF        = xp.findDataRef(self.heading_target_DREF)
        self.flight_director_DREF       = xp.findDataRef(self.flight_director_DREF)
        self.roll_target_DREF           = xp.findDataRef(self.roll_target_DREF)
        self.HDG_yaw_damper_DREF        = xp.findDataRef(self.HDG_yaw_damper_DREF)
        self.autopilot_state_DREF       = xp.findDataRef(self.autopilot_state_DREF)
        self.auto_throttle_DREF         = xp.findDataRef(self.auto_throttle_DREF)

        self.autopilot_ai_override_DREF = xp.findDataRef(self.autopilot_ai_override_DREF)
        self.autopilot_override_DREF    = xp.findDataRef(self.autopilot_override_DREF)

        self.altitide_mode_DREF         = xp.findDataRef(self.altitide_mode_DREF)

        # Register the command handler
        self.XpComm.register()

        self.override_autopilot(True)

        return self.Name, self.Sig, self.Desc
    
    
    def XPluginEnable(self):

        self.autopilot_loop = xp.createFlightLoop(self.autopilot_callback, 
                                    xp.FlightLoop_Phase_AfterFlightModel)
        xp.scheduleFlightLoop(self.autopilot_loop, -1, 1)
        return 1
    
    def XPluginDisable(self):
        if self.autopilot_loop:
            xp.destroyFlightLoop(self.autopilot_loop)
            self.autopilot_loop = None
        return

    def XPluginStop(self):
        self.XPluginDisable()
        return 1
    
    def XPluginReceiveMessage(self, inFromWho, inMessage, inParam):
        pass

    def read_custom_physics(self):
        self.physics_status = xp.getDatai(self.physics_DREF) == 1
        pass 

    def read_altitude_target(self):
        """ Reads the altitude target from X-Plane """
        self.altitude_target = xp.getDataf(self.altitude_target_DREF)
        self.altitude_target *= 0.3048 # Convert to meters
        return
    
    def write_altitude_target(self, altitude):
        """ Writes the altitude target to X-Plane """
        xp.setDataf(self.altitude_target_DREF, altitude * 3.28084) # Convert to feet
        return

    def read_vertical_speed_target(self):
        """ Reads the vertical speed target from X-Plane """
        self.vertical_speed_target = xp.getDataf(self.vertical_speed_target_DREF)
        self.vertical_speed_target *= 0.3048 / 60 # Convert to meters
        return
    
    def write_vertical_speed_target(self, vertical_speed):
        """ Writes the vertical speed target to X-Plane """
        xp.setDataf(self.vertical_speed_target_DREF, vertical_speed * 196.85)
        return
    
    def read_airspeed_target(self):
        """ Reads the airspeed target from X-Plane """
        self.airspeed_target = xp.getDataf(self.airspeed_target_DREF)
        if xp.getDatai(self.APisMach_DREF) == 1:
            self.airspeed_target *= 340.29
        else:
            self.airspeed_target *= 0.514444
        return
    
    def write_airspeed_target(self, airspeed):
        """ Writes the airspeed target to X-Plane """
        xp.setDataf(self.airspeed_target_DREF, airspeed * 1.94384)
        return
    
    def read_heading_target(self):
        """ Reads the heading target from X-Plane """
        self.heading_target = xp.getDataf(self.heading_target_DREF)
        self.heading_target = np.deg2rad(self.heading_target)
        self.heading_target %= 2 * np.pi
        return
    
    def write_heading_target(self, heading):
        """ Writes the heading target to X-Plane """
        xp.setDataf(self.heading_target_DREF, np.rad2deg(heading) % 360)
        return
    
    def read_roll_target(self):
        """ Reads the roll target from X-Plane """
        self.roll_target = xp.getDataf(self.roll_target_DREF)
        self.roll_target = np.deg2rad(self.roll_target)

        # Convert to -pi to pi
        self.roll_target = (self.roll_target + np.pi) % (2 * np.pi) - np.pi
        if self.log:
            xp.log(f'Roll Target {self.roll_target}')
        return
    
    def write_roll_target(self, roll: float):
        """ Writes the roll target to X-Plane """
        xp.setDataf(self.roll_target_DREF, np.rad2deg(roll))
        return
    
    def read_HDG_yaw_damper(self):
        """ Reads the yaw damper state from X-Plane """
        self.HDG_yaw_damper = xp.getDatai(self.HDG_yaw_damper_DREF)
        return
    
    def write_HDG_yaw_damper(self, state):
        """ Writes the yaw damper state to X-Plane """
        xp.setDatai(self.HDG_yaw_damper_DREF, state)
        return
    
    def read_heading_mag(self):
        """ Reads the magnetic heading from X-Plane """
        self.heading_mag = xp.getDataf(self.heading_mag_DREF)
        self.heading_mag = np.deg2rad(self.heading_mag)
        self.heading_mag %= 2 * np.pi
        return
    
    def read_flight_director(self) -> int:
        """ Reads the flight director state from X-Plane """
        self.flight_director = int(xp.getDatai(self.flight_director_DREF))
        return self.flight_director
    
    def write_flight_director(self, state):
        """ Writes the flight director state to X-Plane """
        xp.setDatai(self.flight_director_DREF, state)
        return
    
    def read_auto_throttle(self):
        """ Reads the auto throttle state from X-Plane """
        self.auto_throttle = xp.getDatai(self.auto_throttle_DREF)
        return
    
    def write_auto_throttle(self, state):
        """ Writes the auto throttle state to X-Plane """
        xp.setDatai(self.auto_throttle_DREF, state)
        return

    def read_state(self):
        """ Reads the state vector from X-Plane """
        X = []
        xp.getDatavf(self.state_DREF, X, 0, 12)
        self.X = np.array(X)
        self.read_heading_mag()
        for i in range(6, 6 + 2):
            self.X[i] = (self.X[i] + np.pi) % (2 * np.pi) - np.pi
        
        self.X[8] = (self.heading_mag) % (2 * np.pi)

        self.heading_mag = xp.getDataf(self.heading_mag_DREF)
        self.curr_altitude = -self.X[11]
        self.curr_IAS = np.sqrt(np.power(self.X[0], 2) + np.power(self.X[2], 2)) * \
                 np.sin(self.X[7] - np.arctan(self.X[2]/self.X[0]))
        return
    
    def read_control(self):
        """ Reads the control vector from X-Plane """
        U = []
        xp.getDatavf(self.control_DREF, U, 0, 5)
        self.U = np.array(U)
        return
    
    def write_control(self, U):
        """ Writes the control vector to X-Plane """
        xp.setDatavf(self.control_DREF, U, 0, 5)
        return
    
    def override_autopilot(self, state: bool):
        """ Overrides the autopilot """
        xp.setDatai(self.autopilot_override_DREF, int(state))
        return
    
    def read_altitude_mode(self) -> int:
        return xp.getDatai(self.altitide_mode_DREF)

    def write_altitude_mode(self, mode: int):
        xp.setDatai(self.altitide_mode_DREF, mode)
        return
    
    def call_lateral_controller(self) -> np.ndarray:
        """
        Call the lateral controller and return the control inputs
        """
        assert (self.X.shape == (12,))

        UdeltaZero = np.zeros(5)

        if self.active_lateral_contrl == LateralControls.EMPTY:
            return UdeltaZero
        elif self.active_lateral_contrl == LateralControls.WingLeveler:
            if self.log:
                xp.log(f'Wing Leveler called with target {self.roll_target}')
            return self.controller.WL(self.X, self.dt, self.roll_target)
        elif self.active_lateral_contrl == LateralControls.HeadingHold:
            if self.log:
                xp.log(f'Heading hold called with target f{self.heading_target}')
            return self.controller.HDG(self.X,self.dt, self.heading_target)
        else:
            # Will never reach here
            raise NotImplementedError
        
    def call_vertical_controller(self) -> np.ndarray:
        '''
        Call the longitudinal controller and return the control inputs
        '''
        assert (self.X.shape == (12,))

        UdeltaZero = np.zeros(5)

        if self.active_vertical_contrl == VerticalControls.EMPTY:
            return UdeltaZero
        elif self.active_vertical_contrl == VerticalControls.AltitudeHold:
            if self.log:
                xp.log(f'Altitiude hold called with target {self.altitude_target}')
            return self.controller.ALT(self.X, self.dt, self.altitude_target)
        elif self.active_vertical_contrl == VerticalControls.VerticalSpeedHold:
            if self.log:    
                xp.log(f'Vertical Speed caleed with target {self.vertical_speed_target}')
            return self.controller.VS(self.X, self.dt, self.vertical_speed_target)
        else:
            # Will never reach here
            raise NotImplementedError
    def call_throttle_controller(self) -> np.ndarray:
        '''
        Call the throttle controller and return the control inputs
        '''
        assert (self.X.shape == (12,))

        xp.log('Autothrottle CALLED! ')

        UdeltaZero = np.zeros(5)

        if self.active_throttle_contrl == ThrottleControls.EMPTY:
            return UdeltaZero
        elif self.active_throttle_contrl == ThrottleControls.IndicatedAirspeedHold:
            if self.log:
                xp.log(f'IAS called with target {self.airspeed_target}')
            return self.controller.IAS(self.X, self.dt, self.airspeed_target)
        else:
            # Will never reach here
            raise NotImplementedError

    def call_yaw_damper(self) -> np.ndarray:

        assert (self.X.shape == (12,))
        UdeltaZero = np.zeros(5)

        if self.HDG_yaw_damper:
            return self.controller.Yaw_Damper(self.X, self.dt)
        else:
            return UdeltaZero


    def read_autopilot_state(self) -> bool:
        """ Reads the autopilot state from X-Plane """
        self.APState = xp.getDatai(self.autopilot_state_DREF)
        # This is for reference
        # bit_values = [
        #     (1, ThrottleControls.IndicatedAirspeedHold), # Autothrottle Speed Engage
        #     (2, LateralControls.HeadingHold),            # Heading Select Engage
        #     (4, LateralControls.WingLeveler),            # Roll Hold Engage
        #     (8, VerticalControls.EMPTY),                 # Speed-by-pitch Engage
        #     (16, VerticalControls.VerticalSpeedHold),    # V/S Engage
        #     (32, VerticalControls.AltitudeHold),         # Altitude Hold Arm
        #     (64, VerticalControls.EMPTY),                # Flight Level Change Engage
        #     (128, VerticalControls.EMPTY),               # Pitch Hold Engage
        #     (256, LateralControls.EMPTY),                # Nav Armed
        #     (512, LateralControls.EMPTY),                # Nav Engaged
        #     (1024, VerticalControls.EMPTY),              # Glideslope Armed
        #     (2048, VerticalControls.EMPTY),              # Glideslope Engaged
        #     (4096, VerticalControls.EMPTY),              # VNAV Speed Armed
        #     (8192, VerticalControls.EMPTY),              # VNAV Speed Engaged
        #     (16384, VerticalControls.AltitudeHold),      # Altitude Hold Engaged
        #     (32768, LateralControls.EMPTY),              # TO/GA (Lateral)
        #     (65536, VerticalControls.EMPTY),             # TO/GA (Vertical)
        #     (131072, VerticalControls.EMPTY),            # VNAV Path Armed
        #     (262144, VerticalControls.EMPTY),            # VNAV Path Engaged
        #     (524288, LateralControls.EMPTY),             # GPSS Engaged
        #     (1048576, LateralControls.HeadingHold),      # Heading Hold Engaged
        #     (2097152, LateralControls.EMPTY),            # Turn Rate Engaged
        #     (4194304, LateralControls.EMPTY),            # Track Engaged
        #     (8388608, VerticalControls.EMPTY)            # Flight Path Angle Engaged
        # ]

        # Split bit_values into 3 lists for each of the 3 autopilot states
        ThrottleControls_list = [
            (1, ThrottleControls.IndicatedAirspeedHold)
        ]
        LateralControls_list = [
            (2, LateralControls.HeadingHold),
            (4, LateralControls.WingLeveler),
            (256, LateralControls.EMPTY),
            (512, LateralControls.EMPTY),
            (32768, LateralControls.EMPTY),
            (524288, LateralControls.EMPTY),
            (1048576, LateralControls.HeadingHold),
            (2097152, LateralControls.EMPTY),
            (4194304, LateralControls.EMPTY),
        ]
        VerticalControls_list = [
            (8, VerticalControls.EMPTY),
            (16, VerticalControls.VerticalSpeedHold),
            (32, VerticalControls.AltitudeHold),
            (64, VerticalControls.EMPTY),
            (128, VerticalControls.EMPTY),
            (1024, VerticalControls.EMPTY),
            (2048, VerticalControls.EMPTY),
            (4096, VerticalControls.EMPTY),
            (8192, VerticalControls.EMPTY),
            (16384, VerticalControls.AltitudeHold),
            (65536, VerticalControls.EMPTY),
            (131072, VerticalControls.EMPTY),
            (262144, VerticalControls.EMPTY),
            (8388608, VerticalControls.EMPTY)
        ]

        for bit_value, control_type in ThrottleControls_list:
            if self.APState & bit_value:
                self.active_throttle_contrl = control_type
                break
        
        for bit_value, control_type in LateralControls_list:
            if self.APState & bit_value:
                self.active_lateral_contrl = control_type
                break
        for bit_value, control_type in VerticalControls_list:
            if self.APState & bit_value:
                self.active_vertical_contrl = control_type
                break

        # Check if any of the autopilot modes are on
        if self.active_throttle_contrl != ThrottleControls.EMPTY or \
           self.active_lateral_contrl != LateralControls.EMPTY or \
           self.active_vertical_contrl != VerticalControls.EMPTY:
            self.autopilot_engaged = True
        else:
            self.autopilot_engaged = False
        
        # Read Controller targets
        self.read_heading_target()
        self.read_roll_target()
        self.read_altitude_target()
        self.read_vertical_speed_target()
        self.read_airspeed_target()
        self.read_auto_throttle()
        self.read_HDG_yaw_damper()

        # Read
        return self.autopilot_engaged
    
    def read_xplane(self):
        """ Reads the state and control vectors from X-Plane """
        self.read_state()
        self.read_control()
        # Need to be enforced to zero
        self.write_roll_target(0)
        self.override_autopilot(True)
        self.read_autopilot_state()
        return

    def autopilot_callback(self, sinceLast, elapsedTime, counter, refcon):
        if self.log:
            xp.log(f"Autopilot called at time {elapsedTime}")
        
        self.dt = sinceLast
        if self.log:
            xp.log(self.get_status())
        
        self.read_custom_physics()
        if self.read_flight_director() != 2 or not self.physics_status:
            if self.log:
                xp.log(f'Controller reset called at {elapsedTime}')
            self.controller.reset_autopilot()
            self.write_auto_throttle(0)
            return self.autopilot_update_rate
        self.read_xplane()
        # THIS IS A NOT OPTIMIMAL FIX. NEED TO FIND A BETTER WAY TO DO THIS
        # Currently, auto thorottle is always enabled when the autopilot is engaged
        self.write_auto_throttle(1)

        deltaU = np.zeros(5)
        deltaU += self.call_lateral_controller()
        deltaU += self.call_vertical_controller()
        deltaU += self.call_throttle_controller()
        deltaU += self.call_yaw_damper()
        self.write_control(deltaU)

        return self.autopilot_update_rate
    
    def autopilot_mode(self):
        '''This function is used to change the autopilot mode swithing behaviour'''
        pass
        # if self.active_vertical_contrl == VerticalControls.VerticalSpeedHold:
        #     if self.vertical_speed_target >= 0:
        #         # Check the target altitude is greater than current altitude 
        #         if self.altitude_target > self.curr_altitude:
        #             self.write_altitude_mode(AltitudeModes.VerticalSpeed.value)
        #         else:
        #             self.write_altitude_target(self.curr_altitude)
        #             self.write_altitude_mode(AltitudeModes.AltitudeHold.value)
        #             xp.commandOnce(self.XpComm.ALT_HOLD)
        #     else:
        #         if self.altitude_target <= self.curr_altitude:
        #             self.write_altitude_mode(AltitudeModes.VerticalSpeed.value)
        #         else:
        #             self.write_altitude_target(self.curr_altitude)
        #             self.write_altitude_mode(AltitudeModes.AltitudeHold.value)
            
        #     # If target altitude is less than 100ft to current altitude, switch to altitude hold
        #     if np.abs(self.altitude_target - self.curr_altitude) < 100:
        #         self.write_altitude_mode(AltitudeModes.AltitudeHold.value)
        #         xp.commandOnce(self.XpComm.ALT_HOLD)
            
        # if self.active_vertical_contrl == VerticalControls.AltitudeHold:
        #     # Check the target altitude is greater than current altitude
        #     if self.altitude_target > self.curr_altitude:
        #         self.write_altitude_mode(AltitudeModes.AltitudeHold.value)
        #     else:
        #         self.write_altitude_target(self.curr_altitude)
        #         self.write_altitude_mode(AltitudeModes.VerticalSpeed.value)

        #     # If target 
            
            

                

    def get_status(self) -> str:
        '''
        Returns the status of the target state based on the current active controllers
        '''
        self.calc_error()

        status = ''

        # Print Current State
        if self.physics_status:
            status += 'Custom Physics activation detected!\n'
        else:
            status += 'Custom Physics is currently deactivated\n'
        status += 'Current State:\n'
        status += f'Velocity   : u =     {self.X[0]:5.2f} m/s\n'
        status += f'           : v =     {self.X[1]:5.2f} m/s\n'
        status += f'           : w =     {self.X[2]:5.2f} m/s\n'
        status += f'Orientation: phi   = {np.rad2deg(self.X[6]):5.2f} deg\n'
        status += f'           : theta = {np.rad2deg(self.X[7]):5.2f} deg\n'
        status += f'           : psi   = {np.rad2deg(self.X[8]):5.2f} deg\n'
        status += f'Altitude   : h =     {-self.X[11]:5.2f} m\n'
        status += '\n'
        status += '-----------------------------------------------------\n'
        status += 'Controller Status:\n'
        status += f'Active Lateral Controller     : {self.active_lateral_contrl.name}\n'
        if (self.active_lateral_contrl != LateralControls.EMPTY):
            if (self.active_lateral_contrl == LateralControls.WingLeveler):
                status += f'Current Bank Angle Error      : {np.rad2deg(self.target_error["phi_c"]):5.2f} deg\n'
            elif (self.active_lateral_contrl == LateralControls.HeadingHold):
                status += f'Current Heading Error         : {np.rad2deg(self.target_error["psi_c"]):5.2f} deg\n'
        else:
            status += f'No Active Lateral Controller\n'
        status += f'Active Vertical Controller: {self.active_vertical_contrl.name}\n'
        if (self.active_vertical_contrl != VerticalControls.EMPTY):
            if (self.active_vertical_contrl == VerticalControls.AltitudeHold):
                status += f'Current Altitude Error        : {self.target_error["alt_c"]:5.2f} m\n'
            elif (self.active_vertical_contrl == VerticalControls.VerticalSpeedHold):
                status += f'Current Vertical Speed Error  : {self.target_error["vs_c"]:5.2f} m/s\n'
        else:
            status += f'No Active Longitudinal Controller\n'
        status += f'Active Throttle Controller: {self.active_throttle_contrl.name}\n'
        if (self.active_throttle_contrl != ThrottleControls.EMPTY):
            if (self.active_throttle_contrl == ThrottleControls.IndicatedAirspeedHold):
                status += f'Current Airspeed Error        : {self.target_error["u_c"]:5.2f} m/s\n'
        else:
            status += f'No Active Throttle Controller\n'
        status += '\n'
        return status

    def calc_error(self):
        '''Function to check if the target state has been reached'''
        # Check if the lateral target has been reached
        if self.active_lateral_contrl == LateralControls.EMPTY:
            pass
        elif self.active_lateral_contrl == LateralControls.WingLeveler:
            self.target_error['phi_c'] = np.abs(self.X[6] - self.roll_target)
        elif self.active_lateral_contrl == LateralControls.HeadingHold:
            self.target_error['psi_c'] = np.abs(self.heading_mag - self.heading_target) % (2 * np.pi)
        
        # Check if the longitudinal target has been reached
        if self.active_vertical_contrl == VerticalControls.EMPTY:
            pass
        elif self.active_vertical_contrl == VerticalControls.AltitudeHold:
            self.target_error['alt_c'] = np.abs(-self.X[11] - self.altitude_target)
        elif self.active_vertical_contrl == VerticalControls.VerticalSpeedHold:
            vs = np.sqrt(np.power(self.X[0], 2) + np.power(self.X[2], 2)) * \
                 np.sin(self.X[7] - np.arctan(self.X[2]/self.X[0]))
            error = np.abs(vs - self.vertical_speed_target)
            self.target_error['vs_c'] = error
        
        if self.active_throttle_contrl == ThrottleControls.IndicatedAirspeedHold:
            error = np.abs(self.X[0] - self.airspeed_target)
            self.target_error['u_c'] = error