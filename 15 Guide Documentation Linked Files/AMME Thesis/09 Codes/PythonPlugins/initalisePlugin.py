import numpy as np
import configparser
import os
import importlib.util as import_util

from Autopilot.controllerTemplate import Controllers
import FlightDerivatives as lookup

class InitalisePlugin:
    def __init__(self, config_file):

        if not os.path.exists(config_file):
            raise FileNotFoundError('Configuration File Not Pathed Correctly')
        
        configfile = configparser.ConfigParser()
        self.config_file = config_file
        try:
           configfile.read(config_file)
        except FileNotFoundError:
            raise FileNotFoundError(f"Configuration file '{config_file}' not found.")
        
        self.controller_data = None
        self.controller_file = None
        
        self.loadConfig(configfile)

        self.FlightData =  lookup.aircraft_info[self.config['Aircraft_Settings']['aircraft']]['FlightData']


    def loadConfig(self, configfile):
        
        self.config = {}

        if 'Home' not in configfile:
            raise ValueError("The 'Home' section is missing in the configuration file.")
        self.config['home'] = {}
        self.config['home']['XP_Version'] = int(configfile['Home']['XP_Version'])

        if 'Simulation_Settings' not in configfile:
            raise ValueError("The 'Simulation_Settings' section is missing in the configuration file.")
        self.config['sim_settings'] = {}
        self.config['sim_settings']['logging'] = configfile.getboolean('Simulation_Settings','logging')
        
        if 'Aircraft_Settings' not in configfile:
            raise ValueError("The 'Aircraft_Settings' section is missing in the configuration file.")
        self.config['Aircraft_Settings'] = {}
        self.config['Aircraft_Settings']['aircraft'] = configfile['Aircraft_Settings']['aircraft']

        if 'Autopilot_Settings' in configfile:
            self.config['Autopilot_Settings'] = {}
            self.config['Autopilot_Settings']['max_bank_angle'] = np.deg2rad(float(configfile['Autopilot_Settings']['max_bank_angle']))
            self.config['Autopilot_Settings']['max_pitch_rate'] = np.deg2rad(float(configfile['Autopilot_Settings']['max_pitch_rate']))
            
            controller_file_name = configfile['Autopilot_Settings']['autopilot_controller']
            # Controller file are present at StudentControllers folder
            parent_dir = os.path.dirname(self.config_file)
            self.controller_file = os.path.join(parent_dir, 'StudentControllers', controller_file_name)
            if not os.path.exists(self.controller_file):
                raise FileNotFoundError(f"Controller file '{self.controller_file }' not found.")
        
        if 'Navigation_Settings' in configfile:
            self.config['Navigation_Settings'] = {}
            self.config['Navigation_Settings']['Navigation_On'] = configfile.getboolean('Navigation_Settings', 'navigation_on')
            self.config['Navigation_Settings']['preset_path']   = configfile.getboolean('Navigation_Settings', 'preset_path')

    def load_controller(self) -> Controllers | None:
        spec = import_util.spec_from_file_location("Controllers", self.controller_file)
        if spec is None or spec.loader is None:
            raise Exception(f"File loader failed '{self.controller_file }'.")
        controller_module = import_util.module_from_spec(spec)
        if controller_module is None:
            raise Exception(f"Module loader failed '{self.controller_file }'.")
        
        spec.loader.exec_module(controller_module)
        controller = controller_module.Controllers()
        return controller
    
    def logging_on(self):
        return self.config['sim_settings']['logging']

    def autopilot_settings(self):
        return self.config['Autopilot_Settings']
    
    def aircraft_ctrl_limits(self):
        lower = self.FlightData['ControlLimits']['Lower']
        upper = self.FlightData['ControlLimits']['Upper']
        
        return np.array([upper[0], 
                         -0.5* (upper[1]-lower[1]),
                         -0.5* (upper[2]-lower[2]),
                         -0.5* (upper[3]-lower[3]),
                         lower[4]])
    
    def aircraft_control_upper(self):
        return self.FlightData['ControlLimits']['Upper']
    
    def aircraft_control_lower(self):
        return self.FlightData['ControlLimits']['Lower']

    def aircraft_derivatives(self):
        # import lookup

        # FlightData =  lookup.aircraft_info[self.config['Aircraft_Settings']['aircraft']]['FlightData']

        aircraft_data_array = [
            self.FlightData['Inertial']['g'],
            self.FlightData['Inertial']['m'],
            self.FlightData['Inertial']['Ixx'],
            self.FlightData['Inertial']['Iyy'],
            self.FlightData['Inertial']['Izz'],
            self.FlightData['Inertial']['Ixz'],
            self.FlightData['Geometric']['S'],
            self.FlightData['Geometric']['c'],
            self.FlightData['Geometric']['b'],
            self.FlightData['Propeller']['P_max'],
            self.FlightData['Propeller']['eta'],
            self.FlightData['ControlLimits']['Lower'],
            self.FlightData['ControlLimits']['Upper'],
            self.FlightData['Aero']['alpha_o'],
            self.FlightData['Aero']['Cdo'],
            self.FlightData['Aero']['k'],
            self.FlightData['Aero']['CLa'],
            self.FlightData['Aero']['CLq'],
            self.FlightData['Aero']['CLad'],
            self.FlightData['Aero']['CLde'],
            self.FlightData['Aero']['CLdf'],
            self.FlightData['Aero']['CLo'],
            self.FlightData['Aero']['Cyb'],
            self.FlightData['Aero']['Cybd'],
            self.FlightData['Aero']['Cyp'],
            self.FlightData['Aero']['Cyr'],
            self.FlightData['Aero']['Cyda'],
            self.FlightData['Aero']['Cydr'],
            self.FlightData['Aero']['Cmo'],
            self.FlightData['Aero']['Cma'],
            self.FlightData['Aero']['Cmq'],
            self.FlightData['Aero']['Cmad'],
            self.FlightData['Aero']['Cmde'],
            self.FlightData['Aero']['Cmdf'],
            self.FlightData['Aero']['Cnb'],
            self.FlightData['Aero']['Cnbd'],
            self.FlightData['Aero']['Cnp'],
            self.FlightData['Aero']['Cnr'],
            self.FlightData['Aero']['Cnda'],
            self.FlightData['Aero']['Cndr'],
            self.FlightData['Aero']['Clb'],
            self.FlightData['Aero']['Clbd'],
            self.FlightData['Aero']['Clp'],
            self.FlightData['Aero']['Clr'],
            self.FlightData['Aero']['Clda'],
            self.FlightData['Aero']['Cldr']
        ]

        return np.concatenate([np.array(arr, dtype=float).flatten() for arr in aircraft_data_array])
