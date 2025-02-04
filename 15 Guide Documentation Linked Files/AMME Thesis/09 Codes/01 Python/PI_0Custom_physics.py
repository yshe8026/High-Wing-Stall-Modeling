from XPPython3 import xp
import physicsEngine.flightModel as fm
import numpy as np

import os

import UIElements as ui
import initalisePlugin as ip

from CoodinateTransform.transform import LocalToBody, BodyToLocal

class PythonInterface:

    @staticmethod
    def printState(X):
        '''
        Prints the state vector
        '''
        string = ''
        string += f'Velocity: {X[0:3]} \n'
        string += f'PQR     : {X[3:6]} \n'
        string += f'Angles  : {X[6:9]} \n'
        string += f'Position: {X[9:12]} \n'
        xp.log(string)

    """
    ======================================================================================================================
        \brief  Initalise Plug-In:
                + Set up data references
                + Set up hotkeys
                + Initalise class variables

        NOTE: The class is initalised once, when X-Plane is first booted.
              If the user wants to update any settings, one must re-initalise class:
                    1: Locate the plug-in drop down menu
                    2. XPPython3 -> Reload Script
    ======================================================================================================================
    """
    def __init__(self):
        self.Name = "Physics Engine"
        self.Sig  = "flightsim.Python.custom_physics"
        self.Desc = "A plug-in that allows the user to control the aircraft physics"

        self.positionDREF = ["sim/flightmodel/position/local_x",
                             "sim/flightmodel/position/local_y",
                             "sim/flightmodel/position/local_z"]
        self.velocityDREF = ["sim/flightmodel/position/local_vx", 
                             "sim/flightmodel/position/local_vy", 
                             "sim/flightmodel/position/local_vz"]
        self.angleDREF = ["sim/flightmodel/position/phi", 
                          "sim/flightmodel/position/theta", 
                          "sim/flightmodel/position/psi"]
        self.angular_velocityDREF = ["sim/flightmodel/position/Prad", 
                                     "sim/flightmodel/position/Qrad", 
                                     "sim/flightmodel/position/Rrad"]
        self.quaternionDREF = "sim/flightmodel/position/q"
        self.altitudeDREF = "sim/flightmodel/position/elevation"

        # self.trimDREF = [
        #     'sim/flight_controls/pitch_trim_down',
        #     'sim/flight_controls/pitch_trim_up',
        #     'sim/flight_controls/aileron_trim_right',
        #     'sim/flight_controls/aileron_trim_left'
        # ]

        ############################################################################################################
        # Richard Trim (place 1 of 5)
        self.trimDREF = [
            'sim/cockpit2/controls/elevator_trim',
            'sim/cockpit2/controls/aileron_trim',
            'sim/cockpit2/controls/rudder_trim'
        ]
        ############################################################################################################

        self.controlDREF = ['sim/flightmodel/engine/ENGN_thro' ,
                            'sim/joystick/yoke_pitch_ratio',
                            'sim/joystick/yoke_roll_ratio', 
                            'sim/joystick/yoke_heading_ratio', 
                            'sim/flightmodel/controls/flaprqst'] # Richard changed this from 'sim/flightmodel/controls/flaprat' to 'sim/flightmodel/controls/flaprqst' for the flap fix
        self._controlDescription = ['Throttle', 'Elevator', 'Aileron', 'Rudder', 'Flaps']

        self.flightDirectorDREF  = 'sim/cockpit2/autopilot/flight_director_mode'
        self.disable_controlDREF = 'sim/operation/override/override_joystick'

        # THIS CUSTOM DREF IS CURRENTLY USED BY THE AUTOPILOT PLUGIN
        self.physics_statusDREF  = 'customPhysics/physics_status'
        
        self.control_surface_aileronDREF   = 'sim/flightmodel2/wing/aileron1_deg'
        self.control_surface_elevatorDREF  = 'sim/flightmodel2/wing/elevator1_deg'
        self.control_surface_rudderDREF    = 'sim/flightmodel2/wing/rudder1_deg'
        ############################################################################################################
        # Richard Flap (place 1 of 3)
        self.control_surface_flapDREF    = 'sim/flightmodel2/wing/flap1_deg'
        ############################################################################################################
        self.override_control_surfacesDREF = 'sim/operation/override/override_control_surfaces'

        self.on_groundDREF = 'sim/flightmodel/failures/onground_any'

        # Data ref for pausing the sim
        self.pauseSimDREF = 'sim/time/paused'

        # Data ref used for activating physics
        ########### THIS COMMAND IS FOR PHYSICS ACTIVATION #############
        self.physics_command = 'customPhysics/physics_command'
        self.physics_command_toggle_lastcall = 0

        # Hotkeys 
        self.physics_hotkey = None
        self.disable_autopilot_hkey = None

        # Vector to store values from the sim
        self.Pos = np.zeros(3)
        self.Vel = np.zeros(3)
        self.PQR = np.zeros(3)
        self.Ang = np.zeros(3)
        ############################################################################################################
        # Richard Trim (place 2 of 5)
        self.Trim = np.zeros(3)
        ############################################################################################################
        self.Control = [[], 0, 0, 0, 0]
        self.q = []

        self.UStick = np.zeros(5)

        self.UStickRef = np.zeros(5)
        self.UStickRefSet = False

        self.U_autopilot = np.zeros(5)

        self.X = np.zeros(12)

        # Load flight physics 
        self.flightModel = fm.flightModel()

        config_path = os.path.join('Resources', 'plugins', 'PythonPlugins', 'config.ini')
        self.user_defined = ip.InitalisePlugin(config_path)
    
        self.flightModel.update_parameters(self.user_defined.aircraft_derivatives())
        self.mean_ctrl_limits = self.user_defined.aircraft_ctrl_limits()
        self.ctrl_limits_upper = self.user_defined.aircraft_control_upper()
        self.ctrl_limits_lower = self.user_defined.aircraft_control_lower()

        # Logging 
        self.log = self.user_defined.logging_on()

        # Flag to start/stop the flight model
        self.flight_model_on = False

        self.state_accessor = None

        self.UI = ui.UIElements()

        # X-Plane autopilot variables
        self.altitude_target = "sim/cockpit/autopilot/altitude" # Dref for storing altitude value 
        self.vertical_speed_target  = "sim/cockpit/autopilot/vertical_velocity" #Dref for storing vertical speed values 
        self.airspeed_target  = "sim/cockpit/autopilot/airspeed" # Dref for airspeed 
        self.APisMach = "sim/cockpit/autopilot/airspeed_is_mach" # boolean to change airspeed to mach 
        self.heading_target = "sim/cockpit/autopilot/heading_mag" # Dref for magnetic heading 
        self.flight_director = "sim/cockpit2/autopilot/flight_director_mode" # 0 is off, 1 is on (Flight director)
        self.autopilot_state = "sim/cockpit/autopilot/autopilot_state"

        self.panel_update_rate  = 1 # seconds
        self.panel_loop         = None
        self.autopilotDisengagePercentage = 0.3
    
    """
    ======================================================================================================================
        \brief  Samples hotkey status, enabling/disabling the custom physics module depending on hotkey status.
    ======================================================================================================================
    """
    def flight_model_hotkey(self, nonarg1=None, nonarg2 = None, nonarg3=None):
        '''
        Enables/Disables the custom physics
        '''
        current_time = xp.getElapsedTime()
        # This is to ensure that the toggle is not called multiple times 
        # in a second
        if current_time - self.physics_command_toggle_lastcall < 1:
            return 0
        else:
            self.physics_command_toggle_lastcall = current_time

        xp.log('Hotkey pressed')

        if self.getOnGroundStatus(): 
            self.flight_model_on = False
            self.overrideControlSurfaces(False)
            self.UI.ShowMessage('Cannot start physics while on ground!')
            return 2

        # Toggle the flag
        # if flight model is on, turn it off
        if self.flight_model_on: 
            # Disable the flight model
            self.flight_model_on = False
            # Enable control surfaces
            self.overrideControlSurfaces(False)
            self.UI.ShowMessage('R35 Custom Physics: Off')
            xp.log('Custom Physics Stopped')
        else:
            # Enable the flight model
            self.flight_model_on = True
            # Disable control surfaces
            self.overrideControlSurfaces(True)
            self.UI.ShowMessage('R35 Custom Physics: On')
            xp.log('Custom Physics Started')
        return 2
    

    """
    ======================================================================================================================
        \brief  XPluginStart:
                + Required by XPPython3
                + Called once by X-Plane on startup (or when plugins are re-starting as part of reload) to set data accessors
                  or datarefs.
                + New Datarefs/Accessors must be registered here.
                + You need to return three strings:
                    - self.Name, self.Sig, self.Desc
                        self.Name:  Plug-in Name
                        self.Sig:   Plug-in Signature
                        self.Desc:  Plug-in Description
                    - These are used by X-Plane's Plug-in manager
    ======================================================================================================================
    """
    def XPluginStart(self):

        for i in range(len(self.positionDREF)):
            self.positionDREF[i] = xp.findDataRef(self.positionDREF[i])

        for i in range(len(self.velocityDREF)):
            self.velocityDREF[i] = xp.findDataRef(self.velocityDREF[i])
        
        for i in range(len(self.angleDREF)):
            self.angleDREF[i] = xp.findDataRef(self.angleDREF[i])
        
        for i in range(len(self.angular_velocityDREF)):
            self.angular_velocityDREF[i] = xp.findDataRef(self.angular_velocityDREF[i])

        self.quaternionDREF = xp.findDataRef(self.quaternionDREF)
        self.altitudeDREF = xp.findDataRef(self.altitudeDREF)

        for i in range(len(self.controlDREF)):
            self.controlDREF[i] = xp.findDataRef(self.controlDREF[i])

        for i in range(len(self.trimDREF)):
            self.trimDREF[i] = xp.findDataRef(self.trimDREF[i])

        self.pauseSimDREF = xp.findDataRef(self.pauseSimDREF)

        self.flightDirectorDREF = xp.findDataRef(self.flightDirectorDREF)

        self.disable_controlDREF = xp.findDataRef(self.disable_controlDREF)
        
        # Drefs for control surfaces
        self.control_surface_aileronDREF   = xp.findDataRef(self.control_surface_aileronDREF)
        self.control_surface_elevatorDREF  = xp.findDataRef(self.control_surface_elevatorDREF)
        self.control_surface_rudderDREF    = xp.findDataRef(self.control_surface_rudderDREF)
        ############################################################################################################
        # Richard Flap (place 2 of 3)
        self.control_surface_flapDREF    = xp.findDataRef(self.control_surface_flapDREF)
        ############################################################################################################
        self.override_control_surfacesDREF = xp.findDataRef(self.override_control_surfacesDREF)

        # Dref for on ground
        self.on_groundDREF = xp.findDataRef(self.on_groundDREF)

        # Commands for starting physics
        xp.createCommand(self.physics_command, 'Start/Stop Custom Physics')
        self.physics_command = xp.findCommand(self.physics_command)
        xp.registerCommandHandler(self.physics_command, self.flight_model_hotkey, 1, 0)


        # Use F to start and stop the flight model
        self.physics_hotkey = xp.registerHotKey(xp.VK_9, xp.DownFlag, "Start/Stop Custom Physics",
                                                 self.flight_model_hotkey, 0)
        
        self.state_accessor = xp.registerDataAccessor('customPhysics/X_state', 
                                    readFloatArray=self.get_state_dataref)
        self.control_accessor = xp.registerDataAccessor('customPhysics/U_state',
                                    readFloatArray=self.get_control_dataref,
                                    writeFloatArray=self.write_control_dataref)
        
        self.physics_status = xp.registerDataAccessor(self.physics_statusDREF,
                                            readInt=self.getPhysicsStatus)
        
        # X plane autopilot datarefs
        self.altitude_target            = xp.findDataRef(self.altitude_target)
        self.vertical_speed_target      = xp.findDataRef(self.vertical_speed_target)
        self.airspeed_target            = xp.findDataRef(self.airspeed_target)
        self.APisMach                   = xp.findDataRef(self.APisMach)   
        self.heading_target             = xp.findDataRef(self.heading_target)
        self.flight_director            = xp.findDataRef(self.flight_director)
        self.autopilot_state            = xp.findDataRef(self.autopilot_state)


        return self.Name, self.Sig, self.Desc

    """
    ======================================================================================================================
        \brief  Get the current status of the physics 
    ======================================================================================================================
    """
    def getPhysicsStatus(self,readRefCon):
        return int(self.flight_model_on)
    
    def getOnGroundStatus(self):
        return xp.getDatai(self.on_groundDREF)
    
    def overrideControlSurfaces(self, override: bool = True):
        xp.setDatai(self.override_control_surfacesDREF, int(override))
        return

    
    """
    ======================================================================================================================
        \brief  Function to get state dataref
    ======================================================================================================================
    """
    def get_state_dataref(self, refCon, values, offset, count):
        if values is None:
            return len(self.X)
        values.extend(self.X[offset:offset + count])
        return min(count, len(self.X) - offset)
    
    """
    ======================================================================================================================
        \brief  Function to get control dataref
    ======================================================================================================================
    """
    def get_control_dataref(self, refCon, values, offset, count):
        if values is None:
            return len(self.U)
        U = self.mapStickInputToU(self.UStick)
        values.extend(U[offset:offset + count])
        return min(count, 5 - offset)
    
    """
    ======================================================================================================================
        \brief   Write to the control dataref. This function changes the control inputs by adding the values to the current 
                 control
    ======================================================================================================================
    """
    def write_control_dataref(self, refCon, values, offset, count):
        if not values:
            xp.log('write_control_dataref: No values provided!')
        self.U_autopilot[offset:offset + count] = values[0:count]
        return

    """
    ======================================================================================================================
        \brief  Called once by X-Plane on quit (or when plugins are exiting as part of reload)
                Return is ignored
    ======================================================================================================================
    """
    def XPluginStop(self):
        if self.flight_loop:
            xp.destroyFlightLoop(self.flight_loop)
            self.flight_loop = None

        if self.position_loop:
            xp.destroyFlightLoop(self.position_loop)
            self.position_loop = None
        
        if self.physics_hotkey:
            xp.unregisterHotKey(self.physics_hotkey)
            self.physics_hotkey = None

        xp.unregisterDataAccessor(self.state_accessor)
        self.state_accessor = None
        xp.unregisterDataAccessor(self.control_accessor)
        self.control_accessor = None
        xp.unregisterDataAccessor(self.physics_status)
        self.physics_hotkey = None

    """
    ======================================================================================================================
        \brief  Required by XPPython3 to enable plug-in when requested.  
                + Creates flight loop, position loop, and hotkey monitoring
    ======================================================================================================================
    """
    def XPluginEnable(self):
        # 
        self.position_loop = xp.createFlightLoop(self.positionCallback, 
                                xp.FlightLoop_Phase_BeforeFlightModel)
        self.flight_loop = xp.createFlightLoop(self.PositionVelocityCallback, 
                                xp.FlightLoop_Phase_AfterFlightModel)
        # self.panel_loop =  xp.createFlightLoop(self.panel_callback,
        #                         xp.FlightLoop_Phase_AfterFlightModel)
        xp.scheduleFlightLoop(self.position_loop, -1, 1)
        xp.scheduleFlightLoop(self.flight_loop, -1, 1)
        # xp.scheduleFlightLoop(self.panel_loop, -1, 1)
        return 1

    """
    ======================================================================================================================
        \brief  Called once by X-Plane, when plugin is requested to be disabled. All plugins are disabled prior to Stop.
                
                + Destroys flight loop, position loop, and hotkey monitoring,

                Return is ignored
    ======================================================================================================================
    """
    def XPluginDisable(self):
        # Destroy the flight loop
        if self.flight_loop:
            xp.destroyFlightLoop(self.flight_loop)
            self.flight_loop = None
        
        # Destroy the position loop
        if self.position_loop:
            xp.destroyFlightLoop(self.position_loop)
            self.position_loop = None

        # Destroy the hotkey
        if self.physics_hotkey:
            xp.unregisterHotKey(self.physics_hotkey)
            self.physics_hotkey = None

        xp.unregisterDataAccessor(self.state_accessor)
        self.state_accessor = None
        xp.unregisterDataAccessor(self.control_accessor)
        self.control_accessor = None
        xp.unregisterDataAccessor(self.physics_status)
        self.physics_hotkey = None
        return

    """
    ======================================================================================================================
        \brief  Pause Sim
    ======================================================================================================================
    """
    def pauseFlightIntegration(self):
        xp.setDatai(self.pauseSimDREF, 1)
        return
    
    """
    ======================================================================================================================
        \brief  Resumes flight after paused
    ======================================================================================================================
    """
    def resumeFlightIntegration(self):
        xp.setDatai(self.pauseSimDREF, 0)
        return
    
    """
    ======================================================================================================================
        \brief   Function gets the current state and control inputs of the aircraft. This is called before the flight model 
                 is integrated
    ======================================================================================================================
    """
    def positionCallback(self, sinceLast, elapsedTime, counter, refcon):
        
        # If paused, return 
        self.paused = xp.getDatai(self.pauseSimDREF)
        if self.paused == 1:
            return -1
        
        self.getState()
        self.getControl()

        # Read the flight director mode
        fd = xp.getDatai(self.flightDirectorDREF)

        if fd <= 1:
            self.U_autopilot[:] = 0
        else:
            if self.UStickRefSet == False:
                self.UStickRef = self.UStick.copy()
                self.UStickRefSet = True
            else:
                # Check if StickRef is within 10% of the current stick input
                if np.any(np.abs(self.UStickRef[1:3] - self.UStick[1:3]) \
                            > self.autopilotDisengagePercentage):
                    # Switch off the autopilot
                    self.U_autopilot[:] = 0
                    self.UStickRefSet = False
                    xp.setDatai(self.flightDirectorDREF, 1)

        # Map to correct stick inputs
        self.U = self.mapStickInputToU(self.UStick)
        self.U += self.U_autopilot
        
        
        if self.log:
            xp.log(f'Before Flight Loop Position: {self.Pos[0]}, {self.Pos[1]}, {self.Pos[2]} at t = {elapsedTime}')

        return -1     

    """
    ======================================================================================================================
        \brief  Called by X-Plane whenever a plugin message is being sent to your plugin. Messages include MSG_PLANE_LOADED, 
                MSG_ENTERED_VR, etc., as described in XPLMPlugin module.
                
                Messages may be custom inter-plugin messages, as defined by other plugins.
                
                Return is ignored
    ======================================================================================================================
    """
    def XPluginReceiveMessage(self, inFromWho, inMessage, inParam):
        pass

    """
    ======================================================================================================================
        \brief  Function to get the current state of the aircraft from the sim
    ======================================================================================================================
    """
    def getState(self) -> np.ndarray:
        '''
        Returns
        ------------
        X : State Vector (12, )
        '''
        # Get the current position
        for i in range(3):
            self.Pos[i] = xp.getDatad(self.positionDREF[i])
        # Get the current velocity
        for i in range(3):
            self.Vel[i] = xp.getDataf(self.velocityDREF[i])
        # Get the current angles
        for i in range(3):
            self.Ang[i] = np.deg2rad(xp.getDataf(self.angleDREF[i]))
        # Get the current angular velocity
        for i in range(3):
            self.PQR[i] = xp.getDataf(self.angular_velocityDREF[i])
        ############################################################################################################
        # Richard Trim (place 3 of 5)
        # Get the trim settings
        for i in range(3):
            self.Trim[i] = xp.getDataf(self.trimDREF[i])        
        ############################################################################################################
        # Get the current altitude
        self.Alt = xp.getDataf(self.altitudeDREF)

        # Convert the velocity from local to body
        self.Vel = LocalToBody(self.Vel, self.Ang)

        # Get the current quaternion
        xp.getDatavf(self.quaternionDREF, self.q)

        # Build the state vector 
        self.X = np.zeros(12)
        self.X[0:3] = self.Vel
        self.X[3:6] = self.PQR
        self.X[6:9] = self.Ang
        self.X[9:12] = self.Pos

        self.X[9:11] = 0
        self.X[11] = -self.Alt

        return self.X
    
    """
    ======================================================================================================================
        \brief  Function to set the current state of the aircraft in the sim
    ======================================================================================================================
    """
    def setState(self, X, dt):
        '''
        Parameters
        ------------
        X : State Vector (12, )
        dt : Time step in seconds

        Returns
        ------------
        None
        '''

        self.Vel = X[0:3]
        self.PQR = X[3:6]
        # Convert the velocity from body to local
        self.Vel = BodyToLocal(self.Vel, self.Ang)
        # self.Vel[0] = 41.1556; self.Vel[1] = 0; self.Vel[2] = 0; self.PQR[0] = 0; self.PQR[1] = 0; self.PQR[2] = 0
        self.Ang = X[6:9]
        # Set the current velocity
        for i in range(3):
            xp.setDataf(self.velocityDREF[i], self.Vel[i])
        # Set the current angular velocity
        for i in range(3):
            xp.setDataf(self.angular_velocityDREF[i], self.PQR[i])
        # Set Euler Angles
        for i in range(3):
            xp.setDataf(self.angleDREF[i], np.rad2deg(self.Ang[i]))

        # Update the quaternion
        self.updateQuaternion()

        # Set the current angles
        for i in range(3):
            xp.setDataf(self.angleDREF[i], np.rad2deg(self.Ang[i]))
        
        # Calculate the new position (Euler integration)
        self.Pos = self.Pos + self.Vel * dt

        # Set the new position
        for i in range(3):
            xp.setDatad(self.positionDREF[i], self.Pos[i])

        ############################################################################################################
        # Richard Trim (place 4 of 5)
        # Set the current Trim settings
        for i in range(3):
            xp.setDataf(self.trimDREF[i], self.Trim[i])
        ############################################################################################################       
        
        return 
    
    """
    ======================================================================================================================
        \brief  Function to get the current control inputs from the sim
                This returns the Joystick inputs
    ======================================================================================================================
    """
    def getControl(self):

        # Get the current engine throttle ratio and store it at self.Control[0]
        self.throttle_len = xp.getDatavf(self.controlDREF[0], self.Control[0])
        for i in range(1, 5):
            self.Control[i] = xp.getDataf(self.controlDREF[i])

        self.UStick = np.zeros(5)
        self.UStick[0] = np.mean(self.Control[0])
        self.UStick[1:] = self.Control[1:]

        return self.UStick
    
    """
    ======================================================================================================================
        \brief   Function to set the control inputs in the sim

        NOTE: THIS FUNCTION MUST BE CALLED AFTER CALLING GET CONTROL
    ======================================================================================================================
    """
    def setControl(self, throttle_stick):

        throttle = [throttle_stick[0]] * self.throttle_len
        xp.setDatavf(self.controlDREF[0], throttle)

        # Set the control surfaces
        for i in range(1, 5):
            xp.setDataf(self.controlDREF[i], throttle_stick[i])
        return 
    
    """
    ======================================================================================================================
        \brief  Sets the state of the aircraft. i.e. Sends the custom aircraft position, orientation, and velocity to 
                X-Plane
    ======================================================================================================================
    """
    def PositionVelocityCallback(self, sinceLast, elapsedTime, counter, refcon):

        # Check if the simulation is paused
        self.paused = xp.getDatai(self.pauseSimDREF)
        if self.paused == 1:
            return -1

        # Get the current state 
        X = self.X
        U = self.U

        ############################################################################################################
        # Richard Trim (place 5 of 5)
        Trim = self.Trim
        #-------------------------------------------------------------------------------------------------------
        # User Input
        # The larger the number, the more readily you can feel the effect of the trim upon pressing trim buttons
        elevator_trim_effectiveness = 0.1 # (0.1) Meaning 5.73 deg up or down trim is possible
        aileron_trim_effectiveness = 0.05 # (0.05) Meaning 2.86 deg left or right trim is possible
        rudder_trim_effectiveness = 0.05 # (0.05) Meaning 2.86 deg left or right trim is possible
        #-------------------------------------------------------------------------------------------------------
        U[1] = U[1] - elevator_trim_effectiveness * Trim[0] # elevator trim
        U[2] = U[2] - aileron_trim_effectiveness * Trim[1] # aileron trim
        U[3] = U[3] - rudder_trim_effectiveness * Trim[2] # rudder trim
        ############################################################################################################   

        ##############################################################################################################################
        # Richard Flap Experimentation (place 1 of 2)

        # U from last time step
        global global_U 
        global global_flap_is_going_down
        global global_flap_is_going_up
        global global_flap_setting_target
        global global_flap_setting_start_point
        global global_flap_setting_intermediate

        if (U[4] - global_U[4]) < np.deg2rad(-10):
            global_flap_setting_target = U[4]
            global_flap_setting_start_point = global_U[4]
            U[4] = global_flap_setting_start_point
            global_flap_setting_intermediate = global_flap_setting_start_point
            global_flap_is_going_down = True

        elif (U[4] - global_U[4]) > np.deg2rad(10):
            global_flap_setting_target = U[4]
            global_flap_setting_start_point = global_U[4]
            U[4] = global_flap_setting_start_point
            global_flap_setting_intermediate = global_flap_setting_start_point
            global_flap_is_going_up = True

        # Check whether flap is going up, going down, or idle (not mentioned)
        if (global_flap_is_going_down == True) and (global_flap_setting_intermediate > global_flap_setting_target):
            global_flap_setting_intermediate = global_flap_setting_intermediate - np.deg2rad(0.2)
            U[4] = global_flap_setting_intermediate

        elif (global_flap_is_going_up == True) and (global_flap_setting_intermediate < global_flap_setting_target):
            global_flap_setting_intermediate = global_flap_setting_intermediate + np.deg2rad(0.2)
            U[4] = global_flap_setting_intermediate

        global_U = U
        ##############################################################################################################################

        if self.log:
            # Print the current control inputs
            xp.log(f'Current Control Inputs: {U}')
            # Print the current state to the console 
            xp.log('Current State: \n')
            self.printState(X)

        Xg = np.zeros_like(X)
        dt = sinceLast; # Time step in seconds

        X_out = self.flightModel.midpoint_integration(dt, X, Xg, U)

        if self.flight_model_on == False:
            return -1

        # Set the new state
        self.setState(X_out, dt)
        self.moveControlSurfaces(U)

        if self.log:
            # Log
            xp.log('New State: \n')
            self.printState(X_out)

        return -1

    """
    ======================================================================================================================
        \brief   Updates the quaternion based on the euler angles and sends it to the sim
    ======================================================================================================================
    """
    def updateQuaternion(self):
        if self.log:
            xp.log(f'Current Quaternion: {self.q}')
        # Convert the euler angles to a quaternion
        phi = self.Ang[0] / 2
        theta = self.Ang[1] / 2
        psi = self.Ang[2] / 2
        # # Calculate the quaternion
        self.q[0] =  np.cos(psi) * np.cos(theta) * np.cos(phi) + np.sin(psi) * np.sin(theta) * np.sin(phi)
        self.q[1] =  np.cos(psi) * np.cos(theta) * np.sin(phi) - np.sin(psi) * np.sin(theta) * np.cos(phi)
        self.q[2] =  np.cos(psi) * np.sin(theta) * np.cos(phi) + np.sin(psi) * np.cos(theta) * np.sin(phi)
        self.q[3] = -np.cos(psi) * np.sin(theta) * np.sin(phi) + np.sin(psi) * np.cos(theta) * np.cos(phi)

        if self.log:
            xp.log(f'New Quaternion: {self.q}')

        # Set the quaternion to the sim
        xp.setDatavf(self.quaternionDREF, self.q)

        return 

    """
    ======================================================================================================================
        \brief  Maps strick input and throttle to Control Vector U (5,)
    ======================================================================================================================
    """
    def mapStickInputToU(self, throttle_stick : np.ndarray):
        '''
        Parameters
        ------------
        throttle_stick : Vector of X plane Control input 
        
        LIMITS : Control surface limits (5, ) (Note LIMITS[0] == 1)

        Returns
        ------------
        U : Control Vector (5, )

        '''
        
        # Map the control inputs to the correct range
        throttle_stick = throttle_stick.flatten() * self.mean_ctrl_limits
        if self.log:
            xp.log('Mapping Stick Inputs to Control Inputs')
            xp.log(f'Current Stick Inputs: {throttle_stick}')
        return throttle_stick
 
    """
    ======================================================================================================================
        \brief  Maps Control Vector U (5,) to stick Input
    ======================================================================================================================
    """
    def mapUtoStickInput(self, U: np.ndarray):
        '''
        Paramters
        -----------
        U       : Control Vector (5,)
        LIMITS  : Control Limtis for Plane (5, )

        Returns
        -----------
        throttle_stick   : Returns the stick input (5, 1)

        Example
        -----------
        throttle_stick = self.mapUtoStickInput(U)
        '''

        throttle_stick = np.clip(U / self.mean_ctrl_limits, -1, 1)
        return throttle_stick


    def moveControlSurfaces(self, control : np.ndarray):
        '''
        Function to move the control surfaces
        '''

        # limit the control surface deflection
        control = np.clip(control, self.ctrl_limits_lower, self.ctrl_limits_upper)

        elevator_len = xp.getDatavf(self.control_surface_elevatorDREF)
        elevatorDeflection = [np.rad2deg(control[1])] * elevator_len
        xp.setDatavf(self.control_surface_elevatorDREF, elevatorDeflection)

        aleron_len = xp.getDatavf(self.control_surface_aileronDREF)
        aleronDeflection = [0] * aleron_len
        for i in range(aleron_len):
            sign = -1 if i % 2 == 1 else 1
            aleronDeflection[i] = np.rad2deg(control[2]) * sign
        xp.setDatavf(self.control_surface_aileronDREF, aleronDeflection)
        
        rudder_len = xp.getDatavf(self.control_surface_rudderDREF)
        rudderDeflection = [control[3]] * rudder_len
        for i in range(rudder_len):
            sign = -1 if i % 2 == 1 else 1
            rudderDeflection[i] = np.rad2deg(control[3]) * sign * (-1)
        xp.setDatavf(self.control_surface_rudderDREF, rudderDeflection)

        ############################################################################################################
        # Richard Flap (place 3 of 3)
        flap_len = xp.getDatavf(self.control_surface_flapDREF)
        flapDeflection = [-np.rad2deg(control[4])] * flap_len
        xp.setDatavf(self.control_surface_flapDREF, flapDeflection)
        ############################################################################################################

        return




##############################################################################################################################
# Richard Flap Experimentation (place 2 of 2)
global_U = np.ones(5)
global_flap_is_going_down = False
global_flap_is_going_up = False
global_flap_setting_target = 0
global_flap_setting_target = 0
global_flap_setting_start_point = 0
global_flap_setting_intermediate = 0
##############################################################################################################################