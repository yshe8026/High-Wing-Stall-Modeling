"""
======================================================================================================================

  This file is property of the University of Sydney: you cannot redistribute
  it and/or modify it under the academic integrity agreement of the University.

 \file    flightModel.py
 \author  Brendan Waters <brendan.waters@sydney.edu.au>, Ankith Anil Das <ankith.anildas@sydney.edu.au>
 \date    15/07/2023

 \brief   This class encapsulates the flight physics for a custom aircraft simulation.
          To improve the users performance, the class is jit compiled through the numba library.
          
          This introduces some restrictions on numerical methods, albeit, necessary to minimise 
          the computational overhead.   

          It is recomended to refer to the numba documentation before attempting to implement
          any code: https://numba.pydata.org/numba-doc/latest/user/overview.html

\Changes  6/08/2023: Removed unwanted reshaping vector to column vector and fixed bugs in calculate_aero_forces.
          (Cx, Cy, Cz) (Cl, Cm, Cn) were calculated as 3x3 matrices instead of 3x1 vectors
======================================================================================================================
"""
import numpy as np

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# For determining relative path
# xp is run at:
# C:\Users\Richard\00 Richard Apps\Steam\steamapps\common\X-Plane 11
# The required modules are located at:
# C:\Users\Richard\00 Richard Apps\Steam\steamapps\common\X-Plane 11\Resources\plugins\PythonPlugins\physicsEngine\M01_xp_calculate_dX_flow_sep_by_dt.py
import sys
import os

# Determine the path to the directory containing the module
module_path = os.path.join(os.getcwd(), 'Resources', 'plugins', 'PythonPlugins', 'physicsEngine')

# Add the module's directory to the Python path
sys.path.append(module_path)
# Modules required by the stall model
from M01_xp_calculate_dX_flow_sep_by_dt import M01_xp_calculate_dX_flow_sep_by_dt
# from M02_update_X_flow_sep import M02_update_X_flow_sep
# from M03_update_global_time import M03_update_global_time
from M04_calculate_alpha_crit import M04_calculate_alpha_crit
from M06_calculate_alpha_dot import M06_calculate_alpha_dot
# from M05_write_variables_3_interp import M05_write_variables_3_interp
# from M07_calculate_CL import M07_calculate_CL
# from M08_calculate_Cm import M08_calculate_Cm
# from M09_calculate_CD import M09_calculate_CD
from M10_xp_calculate_CL_Cm_CD import M10_xp_calculate_CL_Cm_CD
from T02_automation import FlightData_Geometric_c
# from M11_calculate_az_ay_from_buffet_simplified import M11_calculate_az_ay_from_buffet_simplified
from M11_calculate_az_ay_from_buffet import BuffetNoiseGenerator
from scipy.interpolate import interp1d # For remote flight testing input U
from M12_nasa_earth_atmosphere_model import M12_nasa_earth_atmosphere_model
from M12_nasa_mars_atmosphere_model import M12_nasa_mars_atmosphere_model
from scipy.interpolate import interp1d # For Richard Real Time Gust Module
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

from numba import float64   
from numba.experimental import jitclass
import xp # Richard import xp for generating custom log in XPPython3Log # check out on line 455
TINY_NUMBER = 1e-12

spec = [
    ('FlightData_Inertial_g', float64),
    ('FlightData_Inertial_m', float64),
    ('FlightData_Inertial_Ixx', float64),
    ('FlightData_Inertial_Iyy', float64),
    ('FlightData_Inertial_Izz', float64),
    ('FlightData_Inertial_Ixz', float64),
    ('FlightData_Geometric_S', float64),
    ('FlightData_Geometric_c', float64),
    ('FlightData_Geometric_b', float64),
    ('FlightData_Propeller_P_max', float64),
    ('FlightData_Propeller_eta', float64),
    ('FlightData_ControlLimits_Lower', float64[:]),
    ('FlightData_ControlLimits_Upper', float64[:]),
    ('FlightData_Aero_alpha_o', float64),
    ('FlightData_Aero_Cdo', float64),
    ('FlightData_Aero_k', float64),
    ('FlightData_Aero_CLa', float64),
    ('FlightData_Aero_CLq', float64),
    ('FlightData_Aero_CLad', float64),
    ('FlightData_Aero_CLde', float64),
    ('FlightData_Aero_CLdf', float64),
    ('FlightData_Aero_CLo', float64),
    ('FlightData_Aero_Cyb', float64),
    ('FlightData_Aero_Cybd', float64),
    ('FlightData_Aero_Cyp', float64),
    ('FlightData_Aero_Cyr', float64),
    ('FlightData_Aero_Cyda', float64),
    ('FlightData_Aero_Cydr', float64),
    ('FlightData_Aero_Cmo', float64),
    ('FlightData_Aero_Cma', float64),
    ('FlightData_Aero_Cmq', float64),
    ('FlightData_Aero_Cmad', float64),
    ('FlightData_Aero_Cmde', float64),
    ('FlightData_Aero_Cmdf', float64),
    ('FlightData_Aero_Cnb', float64),
    ('FlightData_Aero_Cnbd', float64),
    ('FlightData_Aero_Cnp', float64),
    ('FlightData_Aero_Cnr', float64),
    ('FlightData_Aero_Cnda', float64),
    ('FlightData_Aero_Cndr', float64),
    ('FlightData_Aero_Clb', float64),
    ('FlightData_Aero_Clbd', float64),
    ('FlightData_Aero_Clp', float64),
    ('FlightData_Aero_Clr', float64),
    ('FlightData_Aero_Clda', float64),
    ('FlightData_Aero_Cldr', float64),
]


# @jitclass(spec) # Richard mute it to prevent intervention with standard Python
class flightModel: 
    def __init__(self) -> None:
        """
        Defines the default aircraft (PC-9) physical parameters and aerodynamic coefficients.
        NOTE: These are overwritten by aircraft_info() if config specifies specific plane
        """
        # Inertial Data
        self.FlightData_Inertial_g = 9.81
        self.FlightData_Inertial_m = 2087
        self.FlightData_Inertial_Ixx = 5066
        self.FlightData_Inertial_Iyy = 6578
        self.FlightData_Inertial_Izz = 10975
        self.FlightData_Inertial_Ixz = 203
        
        # Geometric Data
        self.FlightData_Geometric_S = 16.29
        self.FlightData_Geometric_c = 1.652
        self.FlightData_Geometric_b = 10.12
        
        # Propeller Data
        self.FlightData_Propeller_P_max = 950000
        self.FlightData_Propeller_eta = 0.8
        
        # Control Data
        DtoR = np.pi / 180
        self.FlightData_ControlLimits_Lower = np.array([0, -25*DtoR, -25*DtoR, -25*DtoR, -45*DtoR])
        self.FlightData_ControlLimits_Upper = np.array([1, 25*DtoR, 25*DtoR, 25*DtoR,  0*DtoR])
        
        # Aerodynamic Data
        self.FlightData_Aero_alpha_o = -3.0 / 57.3
        self.FlightData_Aero_Cdo = 0.020
        self.FlightData_Aero_k = 0.050
        self.FlightData_Aero_CLa = 5.827
        self.FlightData_Aero_CLq = 7.960
        self.FlightData_Aero_CLad = -1.987
        self.FlightData_Aero_CLde = 0.532
        self.FlightData_Aero_CLdf = 0
        self.FlightData_Aero_CLo = -self.FlightData_Aero_CLa * self.FlightData_Aero_alpha_o
        self.FlightData_Aero_Cyb = -0.507
        self.FlightData_Aero_Cybd = -0.0032
        self.FlightData_Aero_Cyp = -0.128
        self.FlightData_Aero_Cyr = 0.336
        self.FlightData_Aero_Cyda = 0.000
        self.FlightData_Aero_Cydr = 0.050
        self.FlightData_Aero_Cmo = 0.06
        self.FlightData_Aero_Cma = -0.802
        self.FlightData_Aero_Cmq = -17.72
        self.FlightData_Aero_Cmad = -2.210
        self.FlightData_Aero_Cmde = -1.822
        self.FlightData_Aero_Cmdf = 0
        self.FlightData_Aero_Cnb = 0.107
        self.FlightData_Aero_Cnbd = 0.0015
        self.FlightData_Aero_Cnp = -0.0226
        self.FlightData_Aero_Cnr = -0.160
        self.FlightData_Aero_Cnda = 0.0048
        self.FlightData_Aero_Cndr = -0.115
        self.FlightData_Aero_Clb = -0.0852
        self.FlightData_Aero_Clbd = -0.0004
        self.FlightData_Aero_Clp = -0.328
        self.FlightData_Aero_Clr = 0.0776
        self.FlightData_Aero_Clda = -0.164
        self.FlightData_Aero_Cldr = 0.0302      

    """
    ======================================================================================================================
     \brief   Updates the flight derivatives to match specified aircraft in config.ini
    ======================================================================================================================
    """
    def update_parameters(self, parameters: np.ndarray):
        # Inertial Data
        self.FlightData_Inertial_g   = parameters[0]
        self.FlightData_Inertial_m   = parameters[1]
        self.FlightData_Inertial_Ixx = parameters[2]
        self.FlightData_Inertial_Iyy = parameters[3]
        self.FlightData_Inertial_Izz = parameters[4]
        self.FlightData_Inertial_Ixz = parameters[5]
        
        # Geometric Data
        self.FlightData_Geometric_S = parameters[6]
        self.FlightData_Geometric_c = parameters[7]
        self.FlightData_Geometric_b = parameters[8]
        
        # Propeller Data
        self.FlightData_Propeller_P_max = parameters[9]
        self.FlightData_Propeller_eta   = parameters[10]
        
        # Control Data
        self.FlightData_ControlLimits_Lower = parameters[11:16]
        self.FlightData_ControlLimits_Upper = parameters[16:21]
        
        # Aerodynamic Data
        self.FlightData_Aero_alpha_o = parameters[21]
        self.FlightData_Aero_Cdo  = parameters[22]
        self.FlightData_Aero_k    = parameters[23]
        self.FlightData_Aero_CLa  = parameters[24]
        self.FlightData_Aero_CLq  = parameters[25]
        self.FlightData_Aero_CLad = parameters[26]
        self.FlightData_Aero_CLde = parameters[27]
        self.FlightData_Aero_CLdf = parameters[28]
        self.FlightData_Aero_CLo  = parameters[29]
        self.FlightData_Aero_Cyb  = parameters[30]
        self.FlightData_Aero_Cybd = parameters[31]
        self.FlightData_Aero_Cyp  = parameters[32]
        self.FlightData_Aero_Cyr  = parameters[33]
        self.FlightData_Aero_Cyda = parameters[34]
        self.FlightData_Aero_Cydr = parameters[35]
        self.FlightData_Aero_Cmo  = parameters[36]
        self.FlightData_Aero_Cma  = parameters[37]
        self.FlightData_Aero_Cmq  = parameters[38]
        self.FlightData_Aero_Cmad = parameters[39]
        self.FlightData_Aero_Cmde = parameters[40]
        self.FlightData_Aero_Cmdf = parameters[41]
        self.FlightData_Aero_Cnb  = parameters[42]
        self.FlightData_Aero_Cnbd = parameters[43]
        self.FlightData_Aero_Cnp  = parameters[44]
        self.FlightData_Aero_Cnr  = parameters[45]
        self.FlightData_Aero_Cnda = parameters[46]
        self.FlightData_Aero_Cndr = parameters[47]
        self.FlightData_Aero_Clb  = parameters[48]
        self.FlightData_Aero_Clbd = parameters[49]
        self.FlightData_Aero_Clp  = parameters[50]
        self.FlightData_Aero_Clr  = parameters[51]
        self.FlightData_Aero_Clda = parameters[52]
        self.FlightData_Aero_Cldr = parameters[53]

    """
    ======================================================================================================================
     \brief   Calculates the aerodynamic forces and moments
              NOTE Usage: ForceCoeff, MomentCoeff = aero4560_aero(X, Xg, Xdot, U, FlightData)
    ======================================================================================================================
    """
    def calculate_aero_forces(self, X, Xg, Xdot, U):
        # Input:
        # X         = state vector
        # Xg        = gust vector, zero if not given [x;y;z] (m/s)
        # Xdot      = state derivative vector
        # U         = control vector
        # FlightData = dictionary containing flight data
        # Output:
        # ForceCoeff = force coefficients [Cx; Cy; Cz]
        # MomentCoeff = moment coefficients [Cl; Cm; Cn]
        
        # Calculate Airpath Properties
        # Richard Note: u=X[0], v=X[1], w=X[2] are aircraft speeds relative to still air in body axes (still air is considered still relative to ground beneath)
        # Richard Note: ug=Xg[0], vg=Xg[1], wg=Xg[2] are wind gust speeds relative to hypothetical "still air" in body axes (still air is considered still relative to ground beneath)
        # Richard Note: pg=Xg[3], qg=Xg[4], rg=Xg[5] are rotational wind gust speeds relative to hypothetical "still air" in body axes (think of them as vortices in turbulent air)
        # Richard Note: For example, ug=Xg[0]= 10 kn roughly represents a 10 kn tail wind
        # Richard Note: According to this formulation, all expressions like X[0]+Xg[0] need sign flips, becoming like X[0]-Xg[0]. Currently no gust model, so not necessary.
        # Richard Note: From a system dynamics perspective, Xdot = A_aero * (X - Xg) + A_kin * X + B * u (if we think in terms of linearly simplified dynamics), this would be A_aero.
        Vt = np.sqrt((X[0]-Xg[0])**2 + (X[1]-Xg[1])**2 + (X[2]-Xg[2])**2)
        alpha = np.arctan2((X[2]-Xg[2]), (X[0]-Xg[0]))
        beta = np.arcsin((X[1]-Xg[1]) / Vt)

        # Calculate Atmosphere Properties
        #--------------------------------------------------------
        # Custom atmosphere models (place 1 of 2)
        # User input 
        # 0: original Earth atmosphere model
        # 1: NASA Earth atmosphere model
        # 2: NASA Mars atmosphere model
        which_atmosphere_model_to_use = 1
        #--------------------------------------------------------
        if which_atmosphere_model_to_use == 1: # option 1
            rho, a_sound_speed = M12_nasa_earth_atmosphere_model(-X[11])
        elif which_atmosphere_model_to_use == 2: # option 2
            rho, a_sound_speed = M12_nasa_mars_atmosphere_model(-X[11])
        else: # option 0
            rho, a_sound_speed = self.atomosphere_model(-X[11])
        #--------------------------------------------------------
        # np.seterr(all='raise')
        qbar = 0.5 * rho * Vt**2

        # Normalized Parameters
        bo2V = 0.5 * self.FlightData_Geometric_b / Vt
        co2V = 0.5 * self.FlightData_Geometric_c / Vt
        # da_a = (Xdot[2]/X[0] - X[2]*Xdot[0]/X[0]**2) * co2V
        # db_a = (Xdot[1]/X[0] - X[1]*Xdot[0]/X[0]**2) * bo2V
        da_a = (Xdot[2]/(X[0]-Xg[0]) - (X[2]-Xg[2])*Xdot[0]/(X[0]-Xg[0])**2) * co2V
        db_a = (Xdot[1]/(X[0]-Xg[0]) - (X[1]-Xg[1])*Xdot[0]/(X[0]-Xg[0])**2) * bo2V
        p_a = (X[3]-Xg[3]) * bo2V
        q_a = (X[4]-Xg[4]) * co2V
        r_a = (X[5]-Xg[5]) * bo2V

        # ===================================================================
        # Richard Additional Code Block 1 (Start)
        # ===================================================================
        # Declare that we will use the global clock
        global global_time
        # # note: initialization lasts around 2.3 seconds, or 230 iterations, during which Vt and term_5 can be nan
        # init_wait_time = 5

        # Calculate the time rate of change of X_flow_sep
        global global_alpha
        # use the global dt determined by xplane
        global dt # declare dt as a global variable so that I could print it in calculate_aero_forces() (Maybe useless now)

        # Declare the variable as global to ensure it is recognized at the global scope (time step size for the LAST timestep)
        global global_dt
        try:
            # Try to print the value of the global variable 'global_dt'
            print(global_dt)
        except NameError:
            # If 'global_dt' has not been defined, catch the NameError exception
            # Assign the value 0.01 to 'global_dt' since it was not previously defined
            global_dt = 0.01

        # dt_time_step = global_dt
        # alpha_dot = (alpha - global_alpha) / dt_time_step
        # global_alpha = alpha

        # a_1 = 22.5
        # alpha_star = 20 * (np.pi / 180)
        # # tau_1 = 11.93 * (self.FlightData_Geometric_c / Vt)
        # # tau_2 = 6.66 * (self.FlightData_Geometric_c / Vt)
        # # tau_1 = 11.93 * (self.FlightData_Geometric_c / 52.65)
        # # tau_2 = 6.66 * (self.FlightData_Geometric_c / 52.65)
        # # tau_1 = 11.93 * (self.FlightData_Geometric_c / (Vt if Vt == Vt else 52.65)) #uses the trick that nan is not equal to nan
        # # tau_2 = 6.66 * (self.FlightData_Geometric_c / (Vt if Vt == Vt else 52.65))
        # tau_1 = (1.00) * 11.93 * (self.FlightData_Geometric_c / (Vt if global_time > init_wait_time else 52.65))
        # tau_2 = (1.00) * 6.66 * (self.FlightData_Geometric_c / (Vt if global_time > init_wait_time else 52.65))
        global X_flow_sep  # Declare that we will use the global variable
        global dX_flow_sep_by_dt
        # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep
        # dX_flow_sep_by_dt = -(1/tau_1) * X_flow_sep + (1/(2*tau_1)) * (1 - np.tanh(a_1 * (alpha - alpha_star))) * (alpha - tau_2 * alpha_dot)
        
        # global dt_time_step
        # dt_time_step = 0.01


        # # old term naming convention
        # term_0 = 1/(2*tau_1)
        # term_1 = a_1 * (alpha - alpha_star)
        # # term_2 = 1 - np.tanh(a_1 * (alpha - alpha_star))
        # term_2 = 1 - np.tanh(a_1 * ((alpha if global_time > init_wait_time else 0.0873) - alpha_star))
        # # term_3 = alpha - tau_2 * alpha_dot
        # term_3 = (alpha if global_time > init_wait_time else 0.0873) - tau_2 * (alpha_dot if global_time > init_wait_time else 0)
        # term_4 = -(1 / tau_1)
        # term_5 = term_0 * term_2 * term_3

        # # new term naming convention
        # term_0 = -(1 / tau_1)

        # term_5 = a_1 * (alpha - alpha_star)

        # term_2 = 1/(2*tau_1)
        # # term_3 = 1 - np.tanh(a_1 * (alpha - alpha_star))
        # term_3 = 1 - np.tanh(a_1 * ((alpha if global_time > init_wait_time else 0.0873) - alpha_star))
        # # term_4 = alpha - tau_2 * alpha_dot
        # term_4 = (alpha if global_time > init_wait_time else 0.0873) - tau_2 * (alpha_dot if global_time > init_wait_time else 0)

        # term_1 = term_2 * term_3 * term_4


        # # record trimmed value of alpha (through flight test) (rad)
        # alpha_trimmed = 0.029276918705747886
        # # alpha_trimmed = 6.4 * (np.pi / 180)
        # ric_multiplier = 1/alpha_trimmed
        # # Enforce modification on differential equation to achieve dX_flow_sep_by_dt = 0 when X = 1 (not guarateed to be physical)
        # term_3 = ric_multiplier * term_3
        # # revise term_1 accordingly
        # term_1 = term_2 * term_3 * term_4

        # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + term_0 # ok
        # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + (-1) # ok
        # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + np.random.uniform(-0.2, 2.8) % ok
        # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + (term_5 if global_time > 3 else np.random.uniform(-0.2, 2.8))

        # # with old term naming convention
        # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + (term_5 if global_time > init_wait_time else 2.5)

        # # with new term naming convention
        # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + (term_1 if global_time > init_wait_time else 2.5)

        # # prevent unphysically large X_flow_sep due to large global_dt
        # if global_time <= init_wait_time or global_dt >= 0.25:
        #     X_flow_sep = 1
        #     dX_flow_sep_by_dt = 0

        # # make sure X_flow_sep stay between 0 and 1
        # if X_flow_sep > 1:
        #     X_flow_sep = 1
        #     dX_flow_sep_by_dt = 0

        # if X_flow_sep < 0:
        #     X_flow_sep = 0
        #     dX_flow_sep_by_dt = 0

        # dX_flow_sep_by_dt = term_4 + term_0 * term_2 * term_3 # not ok

        # X_flow_sep = X_flow_sep + dX_flow_sep_by_dt * dt_time_step  # Modify the global variable

        # ==================================================================================================
        # script for calculating X_flow_sep from dX_flow_sep_by_dt and dt
        # ==================================================================================================
        # X_flow_sep_file_path = "X_flow_sep.txt"
        # with open(X_flow_sep_file_path, "r") as file:
        #     X_flow_sep = float(file.read().strip()) + dX_flow_sep_by_dt * dt_time_step
        # # Write the new value back to the file without formatting to preserve accuracy
        # with open(X_flow_sep_file_path, "w") as file:
        #     file.write(str(X_flow_sep))
        # print("Current value:", X_flow_sep)
        # ==================================================================================================

        # ==================================================================================================
        # script for locating the home directory
        # ==================================================================================================
        # # Data to write
        # text_to_write = "Hello, world!"

        # # Open the file in write mode
        # with open('newfile.txt', 'w') as file:
        #     file.write(text_to_write)
        # ==================================================================================================

        # ==================================================================================================
        # script for writing X_flow_sep for checking
        # ==================================================================================================
        # whether_write_separate_files_flag = 0
        # if whether_write_separate_files_flag == 1:
        #     # Open the file in write mode
        #     with open('Vt.txt', 'w') as file:
        #         file.write(str(Vt))
        #     with open('X_flow_sep.txt', 'w') as file:
        #         file.write(str(X_flow_sep))
        #     with open('tau_1.txt', 'w') as file:
        #         file.write(str(tau_1))
        #     with open('tau_2.txt', 'w') as file:
        #         file.write(str(tau_2))
        #     with open('alpha_star.txt', 'w') as file:
        #         file.write(str(alpha_star))
        #     with open('global_alpha.txt', 'w') as file:
        #         file.write(str(global_alpha))
        #     with open('alpha.txt', 'w') as file:
        #         file.write(str(alpha))
        #     with open('alpha_dot.txt', 'w') as file:
        #         file.write(str(alpha_dot))
        #     with open('dX_flow_sep_by_dt.txt', 'w') as file:
        #         file.write(str(dX_flow_sep_by_dt))
        #     with open('richard_checked_terms\\term_0.txt', 'w') as file:
        #         file.write(str(term_0))
        #     with open('richard_checked_terms\\term_1.txt', 'w') as file:
        #         file.write(str(term_1))
        #     with open('richard_checked_terms\\term_2.txt', 'w') as file:
        #         file.write(str(term_2))
        #     with open('richard_checked_terms\\term_3.txt', 'w') as file:
        #         file.write(str(term_3))
        #     with open('richard_checked_terms\\term_4.txt', 'w') as file:
        #         file.write(str(term_4))
        #     with open('richard_checked_terms\\term_5.txt', 'w') as file:
        #         file.write(str(term_5))
        #     with open('richard_checked_terms\\global_time.txt', 'w') as file:
        #         file.write(str(global_time))
        #     with open('richard_checked_terms\\term_5_history.txt', 'a') as file:
        #         file.write(str(term_5) + '\n')
        #     with open('richard_checked_terms\\Vt_history.txt', 'a') as file:
        #         file.write(str(Vt) + '\n')
        #     with open('richard_checked_terms\\alpha_history.txt', 'a') as file:
        #         file.write(str(alpha) + '\n')
        #     with open('richard_checked_terms\\alpha_dot_history.txt', 'a') as file:
        #         file.write(str(alpha_dot) + '\n')
        #     with open('richard_checked_terms\\X_flow_sep_history.txt', 'a') as file:
        #         file.write(str(X_flow_sep) + '\n')
            
        # # write combined    
        # # with open('richard_checked_terms\\combined_history.txt', 'a') as file:
        # #     file.write(f"{alpha_dot}\t{X_flow_sep}\n")
        # file_path = 'richard_checked_terms\\combined_history.txt'
        # # column width in terms of character string length
        # cwidth = 25

        # # Function to format the columns
        # def format_line(global_time, global_dt, dX_flow_sep_by_dt, term_0, X_flow_sep, term_1, term_2, term_3, term_4, term_5, alpha, alpha_dot, a_1, alpha_star, tau_1, tau_2, Vt, cwidth):
        #     # Format global_time to 5 significant figures
        #     formatted_global_time = f"{global_time:.5g}"
        #     return f"{formatted_global_time:<{cwidth}}{global_dt:<{cwidth}}{dX_flow_sep_by_dt:<{cwidth}}{term_0:<{cwidth}}{X_flow_sep:<{cwidth}}{term_1:<{cwidth}}{term_2:<{cwidth}}{term_3:<{cwidth}}{term_4:<{cwidth}}{term_5:<{cwidth}}{alpha:<{cwidth}}{alpha_dot:<{cwidth}}{a_1:<{cwidth}}{alpha_star:<{cwidth}}{tau_1:<{cwidth}}{tau_2:<{cwidth}}{Vt:<{cwidth}}\n"

        # try:
        #     # Check if the file is empty
        #     with open(file_path, 'r') as file:
        #         is_empty = file.read().strip() == ""
        # except FileNotFoundError:
        #     # If the file doesn't exist, it will be considered empty
        #     is_empty = True

        # # Write the column titles if the file is empty
        # with open(file_path, 'a') as file:
        #     if is_empty:
        #         file.write(f"{'global_time[s]':<{cwidth}}{'global_dt[s]':<{cwidth}}{'dX_flow_sep_by_dt[-/s]':<{cwidth}}{'term_0[-/s]':<{cwidth}}{'X_flow_sep[-]':<{cwidth}}{'term_1[-/s]':<{cwidth}}{'term_2[-/s]':<{cwidth}}{'term_3[-]':<{cwidth}}{'term_4[-]':<{cwidth}}{'term_5[-]':<{cwidth}}{'alpha[rad]':<{cwidth}}{'alpha_dot[rad/s]':<{cwidth}}{'a_1[-]':<{cwidth}}{'alpha_star[rad]':<{cwidth}}{'tau_1[s]':<{cwidth}}{'tau_2[s]':<{cwidth}}{'Vt[m/s]':<{cwidth}}\n")
        #     file.write(format_line(global_time, global_dt, dX_flow_sep_by_dt, term_0, X_flow_sep, term_1, term_2, term_3, term_4, term_5, alpha, alpha_dot, a_1, alpha_star, tau_1, tau_2, Vt, cwidth))
        # ===================================================================
        # Richard Additional Code Block 1 (End)
        # ===================================================================

        # Lift/Drag Force Coefficients

        # rho/rho_SSL ratio
        sig = rho / 1.2256

        #--------------------------------------------------------------------
        # Richard's extra state variables (TM18)
        global global_n_prop_rot_speed # n_prop_rot_speed from the last time step
        P_prop_power = self.FlightData_Propeller_P_max * (1.1324*sig - 0.1324) * U[0] # W
        n_prop_rot_speed = 0.0004842579 * P_prop_power + 11.6666666667 # rev/s
        # time rate of change of n_prop_rot_speed
        n_prop_rot_speed_dot = (n_prop_rot_speed - global_n_prop_rot_speed) / (global_dt if global_time > 3 else 0.01)
        global_n_prop_rot_speed = n_prop_rot_speed

        global global_q_pitch_rate
        q_pitch_rate = (X[4] if global_time > 3 else 0.00)
        try:
            # Try to print the value
            print(q_pitch_rate)
        except NameError:
            # If variable has not been defined, catch the NameError exception
            # Assign the value to variable since it was not previously defined
            q_pitch_rate = 0.00
        q_pitch_rate_dot = (q_pitch_rate - global_q_pitch_rate) / (global_dt if global_time > 3 else 0.01)
        global_q_pitch_rate = q_pitch_rate

        global global_r_yaw_rate
        r_yaw_rate = (X[5] if global_time > 3 else 0.00)
        try:
            # Try to print the value
            print(r_yaw_rate)
        except NameError:
            # If variable has not been defined, catch the NameError exception
            # Assign the value to variable since it was not previously defined
            r_yaw_rate = 0.00
        r_yaw_rate_dot = (r_yaw_rate - global_r_yaw_rate) / (global_dt if global_time > 3 else 0.01)
        global_r_yaw_rate = r_yaw_rate
        #--------------------------------------------------------------------

        #--------------------------------------------------------------------
        # Richard's propulsion experimentation (TM18)
        whether_to_add_richard_propulsion_model = 1
        #--------------------------------------------------------------------
        if whether_to_add_richard_propulsion_model == 1:
            # P_prop_power = self.FlightData_Propeller_P_max * (1.1324*sig - 0.1324) * U[0] # W
            # n_prop_rot_speed = 0.0004842579 * P_prop_power + 11.6666666667 # rev/s
            # From the Jabiru AFM, the prop diameter is stated as 1.52 m
            D_prop_diameter = 1.524 # m
            J_advance_ratio = (Vt if global_time > 3 else 52.65) / ((n_prop_rot_speed if global_time > 3 else 35) * D_prop_diameter) # if the simulation is not loaded up yet (3 sec wait), make RPM = 2100 and therefore n = 35 rev/s
            # J_advance_ratio = J_advance_ratio * 2 # Just a trick to stretch the plot (because the polynomial below is incorrect)
            # polynomial_coefficients = [-2.67438907e+08,  1.24592095e+09, -2.61422244e+09,  3.26408269e+09,
            #                            -2.69939288e+09,  1.55737948e+09, -6.43083806e+08,  1.91878929e+08,
            #                            -4.12570203e+07,  6.30205417e+06, -6.65243868e+05,  4.63408791e+04,
            #                            -1.97163351e+03,  4.07222763e+01,  2.85202480e+00,  2.45521078e-05]
            # polynomial_expression = np.poly1d(polynomial_coefficients)
            # eta_p_prop_eff = (14.2191 * J_advance_ratio**4 - 22.6185 * J_advance_ratio**3 - 0.5659 * J_advance_ratio**2 + 5.0815 * J_advance_ratio - 0.1058) * 15 # this polynomial need improvement (use graph tool)
            # eta_p_prop_eff = polynomial_expression(J_advance_ratio)
            # eta_p_prop_eff = (-2.6743890653121743e+08 * J_advance_ratio**15 +
            #                    1.2459209493015354e+09 * J_advance_ratio**14 -
            #                    2.6142224426932373e+09 * J_advance_ratio**13 +
            #                    3.2640826865733147e+09 * J_advance_ratio**12 -
            #                    2.6993928791204376e+09 * J_advance_ratio**11 +
            #                    1.5573794758376265e+09 * J_advance_ratio**10 -
            #                    6.4308380623244667e+08 * J_advance_ratio**9 +
            #                    1.9187892870937729e+08 * J_advance_ratio**8 -
            #                    4.1257020279518954e+07 * J_advance_ratio**7 +
            #                    6.3020541739022173e+06 * J_advance_ratio**6 -
            #                    6.6524386816298915e+05 * J_advance_ratio**5 +
            #                    4.6340879054777324e+04 * J_advance_ratio**4 -
            #                    1.9716335106041690e+03 * J_advance_ratio**3 +
            #                    4.0722276272755224e+01 * J_advance_ratio**2 +
            #                    2.8520248036777662 * J_advance_ratio +
            #                    2.4552107788622379e-05)
            # if (J_advance_ratio >= 0) and (J_advance_ratio <= 0.65789474):
            #     eta_p_prop_eff = (-2.5162505287317768e+02 * J_advance_ratio**5 +
            #                         3.5739712390881306e+02 * J_advance_ratio**4 +
            #                         -1.7659295608329072e+02 * J_advance_ratio**3 +
            #                         3.2577623277903101e+01 * J_advance_ratio**2 +
            #                         5.3932820682034999e-01 * J_advance_ratio +
            #                         3.7218972929306204e-02) # an approximate eta curve
            # else:
            #     eta_p_prop_eff = 0
            if (J_advance_ratio >= 0.000) and (J_advance_ratio < 0.54368204):
                eta_p_prop_eff = (-1.0339808355672748e+00 * J_advance_ratio**5 +
                                1.1513141630282613e+00 * J_advance_ratio**4 +
                                8.0980324477669519e-01 * J_advance_ratio**3 +
                                -2.9358586169893126e+00 * J_advance_ratio**2 +
                                2.7568991029961047e+00 * J_advance_ratio +
                                3.0600838902929569e-06)
            elif (J_advance_ratio >= 0.54368204) and (J_advance_ratio < 0.82489689):
                eta_p_prop_eff = (-414.56709062760984 * J_advance_ratio**5 +
                                1341.4163943124554 * J_advance_ratio**4 +
                                -1733.7299437622642 * J_advance_ratio**3 +
                                1117.049389105175 * J_advance_ratio**2 +
                                -357.9291035961872 * J_advance_ratio +
                                46.33549936554715)
            elif (J_advance_ratio >= 0.82489689) and (J_advance_ratio <= 0.8814218387820809):
                eta_p_prop_eff = (-94738.66870330359 * J_advance_ratio**4 +
                                314836.3898499948 * J_advance_ratio**3 +
                                -392340.3248280777 * J_advance_ratio**2 +
                                217290.48582230246 * J_advance_ratio +
                                -45125.05373275733)
            else:
                eta_p_prop_eff = 0
            # eta_p_prop_eff = 5.0

            # Thrust Coefficient
            C_t = (self.FlightData_Propeller_P_max * (1.1324*sig - 0.1324) * eta_p_prop_eff * U[0] / Vt) / qbar / self.FlightData_Geometric_S # now prop efficiency is not a constant
        else:
            # Assume just some typical values for the data recorder module (and these variables will not be used down the line in the simulation)
            J_advance_ratio = 0.8
            eta_p_prop_eff = 0.8
            # Thrust Coefficient
            C_t = (self.FlightData_Propeller_P_max * (1.1324*sig - 0.1324) * self.FlightData_Propeller_eta * U[0] / Vt) / qbar / self.FlightData_Geometric_S
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Turn the stall model on (1) or off (0)
        whether_activate_stall_model = 1
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if whether_activate_stall_model == 1:
            alpha_dot = M06_calculate_alpha_dot(global_dt, global_alpha, alpha)
            # 'M01_calculate_dX_flow_sep_by_dt' updates global_alpha automatically
            dX_flow_sep_by_dt, global_alpha= M01_xp_calculate_dX_flow_sep_by_dt(global_time, global_dt, X_flow_sep, global_alpha, alpha, FlightData_Geometric_c, Vt)
            # X_flow_sep = M02_update_X_flow_sep(X_flow_sep, dX_flow_sep_by_dt, dt)

            alpha_crit_1 = M04_calculate_alpha_crit(1, global_time, dX_flow_sep_by_dt, alpha_dot, FlightData_Geometric_c, Vt)
            # Choose alpba_crit interpretation version
            alpha_crit = alpha_crit_1
            # Combine three modules M07, M08, M09 into one M10, to ease investigation of parameters.
            CL, poststall_on, Cm, CD = M10_xp_calculate_CL_Cm_CD(Vt, FlightData_Geometric_c, alpha, alpha_dot, alpha_crit, X_flow_sep, U)
            # global_time = M03_update_global_time(global_time, dt)
            # Change variable according to convention
            Cd = CD
        else:
            CL = self.FlightData_Aero_CLo + self.FlightData_Aero_CLa * alpha + self.FlightData_Aero_CLad*da_a + self.FlightData_Aero_CLq*q_a + self.FlightData_Aero_CLde*U[1] + self.FlightData_Aero_CLdf*U[4]
            Cd = self.FlightData_Aero_Cdo + self.FlightData_Aero_k * CL**2
            Cm = self.FlightData_Aero_Cmo + self.FlightData_Aero_Cma*alpha + self.FlightData_Aero_Cmad*da_a + self.FlightData_Aero_Cmq*q_a + self.FlightData_Aero_Cmde*U[1] + self.FlightData_Aero_Cmdf*U[4]
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

        #--------------------------------------------------------------------
        # Richard's propeller wash experimentation (TM18)
        whether_to_add_richard_propeller_wash = 1
        #--------------------------------------------------------------------
        if whether_to_add_richard_propeller_wash == 1:
            CL = CL * (1 + 0.1154 * (1 + 4.396 * C_t))
            CD = CD * (1 + 0.1154 * (1 + 4.396 * C_t))
        #--------------------------------------------------------------------

        # Force Coefficients
        Sa = np.sin(alpha)
        Ca = np.cos(alpha)
        Cx = -Cd*Ca + CL*Sa + C_t
        Cy = self.FlightData_Aero_Cyb*beta + self.FlightData_Aero_Cybd*db_a + self.FlightData_Aero_Cyp*p_a + self.FlightData_Aero_Cyr*r_a + self.FlightData_Aero_Cyda*U[2] + self.FlightData_Aero_Cydr*U[3] * (-1)
        Cz = -CL*Ca - Cd*Sa
        #--------------------------------------------------------------------
        # Richard's buffet experimentation
        whether_to_add_richard_buffet_old_version = 0
        #--------------------------------------------------------------------
        if whether_to_add_richard_buffet_old_version == 1:
            # Experiment buffet model based on Delft Paper (p 13 - 14 of 28)
            if X_flow_sep < 0.89:
            # if X_flow_sep <= 1:
                # az = np.random.uniform(-0.2, 0.2)
                az_noise_mean = 0.0
                az_noise_std = 0.5
                az_noise = np.random.normal(az_noise_mean, az_noise_std)
                Cz = Cz + ((2 * self.FlightData_Inertial_m * az_noise) / (rho * Vt**2 * self.FlightData_Geometric_S))
            else:
                # Else set acceleration noise variables to zero
                az_noise = 0
                ay_noise = 0
        else:
            # Keep variables just for flight recorder
            az_noise = 0
            ay_noise = 0
        #--------------------------------------------------------------------
        #--------------------------------------------------------------------
        # Richard's buffet experimentation
        whether_to_add_richard_buffet_intermediate_version = 0
        #--------------------------------------------------------------------
        if whether_to_add_richard_buffet_intermediate_version == 1:
            # Experiment buffet model based on Delft Paper (p 13 - 14 of 28)
            if (X_flow_sep < 0.89) and (X_flow_sep >= 0):
            # if X_flow_sep <= 1:
                # # az = np.random.uniform(-0.2, 0.2)
                # try:
                #     # Try to print the value of the global variable 'global_dt'
                #     print(X_flow_sep)
                #     X_flow_sep_parameter = X_flow_sep
                # except NameError:
                #     # If 'X_flow_sep' has not been defined, catch the NameError exception
                #     X_flow_sep_parameter = 1
                az_noise_mean = 0.0
                az_noise_std = 1.0
                az_noise = np.random.normal(az_noise_mean, az_noise_std) 
                az_noise = az_noise * (1 - X_flow_sep)

                ay_noise_mean = 0.0
                ay_noise_std = 0.5
                ay_noise = np.random.normal(ay_noise_mean, ay_noise_std) 
                ay_noise = ay_noise * (1 - X_flow_sep)

                Cz = Cz + ((self.FlightData_Inertial_m * az_noise) / (qbar * self.FlightData_Geometric_S))
                Cy = Cy + ((self.FlightData_Inertial_m * ay_noise) / (qbar * self.FlightData_Geometric_S))
            else:
                # Else set acceleration noise variables to zero
                az_noise = 0
                ay_noise = 0
        else:
            # Keep variables just for flight recorder
            az_noise = 0
            ay_noise = 0
        #--------------------------------------------------------------------
        #--------------------------------------------------------------------
        # Richard's buffet experimentation
        whether_to_add_richard_buffet = 0
        #--------------------------------------------------------------------
        if whether_to_add_richard_buffet == 1:
            # Experiment buffet model based on Delft Paper (p 13 - 14 of 28)
            if X_flow_sep < 0.89:

                whether_to_use_full_formulation = 0 # Warning: full formulation is more accurate but much slower for the sim to run (Place 1 of 2)

                if whether_to_use_full_formulation == 1:
                    az_noise, ay_noise = buffet_noise_generator.M11_calculate_az_ay_from_buffet(X_flow_sep, dt)
                else:
                    # az_noise, ay_noise = M11_calculate_az_ay_from_buffet_simplified(X_flow_sep) # this option need some trouble shoot (extremely slow)
                    az_noise_mean = 0.0
                    az_noise_std = 1.0
                    az_noise = np.random.normal(az_noise_mean, az_noise_std) * (1-X_flow_sep)

                    ay_noise_mean = 0.0
                    ay_noise_std = 0.5
                    ay_noise = np.random.normal(ay_noise_mean, ay_noise_std) * (1-X_flow_sep)

                # Calculate az_noise and ay_noise using the BuffetNoiseGenerator
                Cz = Cz + ((self.FlightData_Inertial_m * az_noise) / (qbar * self.FlightData_Geometric_S))
                Cy = Cy + ((self.FlightData_Inertial_m * ay_noise) / (qbar * self.FlightData_Geometric_S))
            else:
                # Else set acceleration noise variables to zero
                az_noise = 0
                ay_noise = 0
        else:
            # Keep variables just for flight recorder
            az_noise = 0
            ay_noise = 0
        #--------------------------------------------------------------------

        ##############################################################################################################################
        # Richard Taped Stall Buffet Module (place 1 of 2)
        whether_to_add_richard_taped_stall_buffet_module = 1
        #------------------------------------------------------------------------------------------------------
        if whether_to_add_richard_taped_stall_buffet_module == 1:

            global Ay_Az_noises_data_matrix_taped

            global t_taped_for_Ay_Az_noises
            global interpolated_Ay_Az_noises_taped

            global data_recording_current_duration 

            # Experiment buffet model based on Delft Paper (p 13 - 14 of 28)
            if (X_flow_sep < 0.89) and (X_flow_sep >= 0):

                # Taped Ay_Az_noises vs t formulation by interpolation!
                Ay_Az_noises_tape_run_duration = data_recording_current_duration

                while Ay_Az_noises_tape_run_duration >= t_taped_for_Ay_Az_noises[-1]:
                    # Subtract the last element until the duration falls within our desired range
                    Ay_Az_noises_tape_run_duration = Ay_Az_noises_tape_run_duration - t_taped_for_Ay_Az_noises[-1]*(4/5)

                if t_taped_for_Ay_Az_noises[1] < Ay_Az_noises_tape_run_duration < t_taped_for_Ay_Az_noises[-1]:

                    Ay_Az_noises = interpolated_Ay_Az_noises_taped(Ay_Az_noises_tape_run_duration)

                az_noise = Ay_Az_noises[1]
                ay_noise = Ay_Az_noises[0]

                az_noise = az_noise * (1 - X_flow_sep)
                ay_noise = ay_noise * (1 - X_flow_sep)

                Cz = Cz + ((self.FlightData_Inertial_m * az_noise) / (qbar * self.FlightData_Geometric_S))
                Cy = Cy + ((self.FlightData_Inertial_m * ay_noise) / (qbar * self.FlightData_Geometric_S))
            else:
                # Else set acceleration noise variables to zero
                az_noise = 0
                ay_noise = 0
        else:
            # Keep variables just for flight recorder
            az_noise = 0
            ay_noise = 0

        ##############################################################################################################################

        ForceCoeff = np.array([Cx, Cy, Cz])

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Richard's aileron_total_deflection_ratio = (13+24)/(24+24)
        # This factor is to account for the asymmetry of aileron deflection
        # From Manual: Down 13 deg, Up 24 deg
        # The following relation should be considered:
        # delta_a:                  -24 deg ~ 24 deg
        # port side aileron:        -13 deg ~ 24 deg
        # starboard side aileron:    24 deg ~ -13 deg
        aileron_total_deflection_ratio = 0.770833
        # This should be used in Cl, Cn equation
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # Moment Coefficients
        Cl = self.FlightData_Aero_Clb*beta + self.FlightData_Aero_Clbd*db_a + self.FlightData_Aero_Clp*p_a + self.FlightData_Aero_Clr*r_a + self.FlightData_Aero_Clda * aileron_total_deflection_ratio * U[2] + self.FlightData_Aero_Cldr*U[3] * (-1)
        # Cm = self.FlightData_Aero_Cmo + self.FlightData_Aero_Cma*alpha + self.FlightData_Aero_Cmad*da_a + self.FlightData_Aero_Cmq*q_a + self.FlightData_Aero_Cmde*U[1] + self.FlightData_Aero_Cmdf*U[4]
        Cn = self.FlightData_Aero_Cnb*beta + self.FlightData_Aero_Cnbd*db_a + self.FlightData_Aero_Cnp*p_a + self.FlightData_Aero_Cnr*r_a + self.FlightData_Aero_Cnda * aileron_total_deflection_ratio * U[2] + self.FlightData_Aero_Cndr*U[3] * (-1)
        #--------------------------------------------------------------------

        #--------------------------------------------------------------------
        # Richard's torque effect experimentation (TM18)
        whether_to_add_richard_torque_effect = 1
        #--------------------------------------------------------------------
        if whether_to_add_richard_torque_effect == 1:

            # utilizes n_prop_rot_speed_dot
            I_prop = 0.3 # axial mass moment of inertia of propeller (kg.m^2)
            Cl = Cl - ((I_prop * (2*np.pi) * n_prop_rot_speed_dot)/(qbar * self.FlightData_Geometric_S * self.FlightData_Geometric_b)) # clockwise propeller generates negative Cl via torque effect
        #--------------------------------------------------------------------

        #--------------------------------------------------------------------
        # Richard's spiral slipstream effect experimentation (TM18)
        whether_to_add_richard_spiral_slipstream_effect = 1
        #--------------------------------------------------------------------
        if whether_to_add_richard_spiral_slipstream_effect == 1:

            c_bar_vtail = 0.7939 # mean chord of vertical tail (m)
            l_vtail = 3.788 # moment arm of the center of pressure of the vertail tail (m)
            Cn = Cn - 0.1 * C_t * (c_bar_vtail/((Vt * np.sqrt(1 + 4.396 * C_t))/n_prop_rot_speed)) * (l_vtail / self.FlightData_Geometric_b)
        #--------------------------------------------------------------------

        #--------------------------------------------------------------------
        # Richard's precession (gyroscopic) effect experimentation (TM18)
        whether_to_add_richard_precession_effect = 1
        #--------------------------------------------------------------------
        if whether_to_add_richard_precession_effect == 1:

            m_prop = 3.63 # mass of propeller (on supplier website) (kg)
            l_prop = 1.9848 # moment arm of propeller around y axis (which in turn goes through center of gravity of J400) (m)
            # q_pitch_rate_dot = (Xdot[4] if global_time > 3 else 0.00) # time rate of change of pitch rate (pitch angular acceleration) (rad/(s^2))
            # r_pitch_rate_dot = (Xdot[5] if global_time > 3 else 0.00) # time rate of change of yaw rate (yaw angular acceleration) (rad/(s^2))
            # q_pitch_rate_dot = 0.00 # time rate of change of pitch rate (pitch angular acceleration) (rad/(s^2))
            # r_yaw_rate_dot = 0.00 # time rate of change of yaw rate (yaw angular acceleration) (rad/(s^2))
            Cm = Cm - ((m_prop * (l_prop**2) * r_yaw_rate_dot)/((qbar if global_time > 3 else 980) * self.FlightData_Geometric_S * self.FlightData_Geometric_c))
            Cn = Cn + ((m_prop * (l_prop**2) * q_pitch_rate_dot)/((qbar if global_time > 3 else 980) * self.FlightData_Geometric_S * self.FlightData_Geometric_b))
            # Cm = Cm - ((m_prop * (l_prop**2) * r_yaw_rate_dot)/(980 * self.FlightData_Geometric_S * self.FlightData_Geometric_c))
            # Cn = Cn + ((m_prop * (l_prop**2) * q_pitch_rate_dot)/(980 * self.FlightData_Geometric_S * self.FlightData_Geometric_b))
        #--------------------------------------------------------------------

        #--------------------------------------------------------------------
        # Richard's p-factor (asymmetric thrust) effect experimentation (TM18)
        whether_to_add_richard_p_factor_effect = 1
        #--------------------------------------------------------------------
        if whether_to_add_richard_p_factor_effect == 1:

            l_p_factor = 0.323403 # moment arm of asymmetric thrust vector 
            k_asymmetric_thrust_factor = 0.2864788976 # slope of asymmetric thrust factor \xi vs alpha (in rad)
            Cm = Cm - C_t * k_asymmetric_thrust_factor * beta * (l_p_factor/self.FlightData_Geometric_c)
            Cn = Cn - C_t * k_asymmetric_thrust_factor * alpha * (l_p_factor/self.FlightData_Geometric_b)
        #--------------------------------------------------------------------

        #--------------------------------------------------------------------
        # Richard's auto-rotation experimentation
        whether_to_add_richard_autorotation = 0
        #--------------------------------------------------------------------
        if whether_to_add_richard_autorotation == 1:
            # Very complicated and related to moments of lift on each wing
            # And also need CL - alpha curve considered
            # Need alpha values on left and right wings
            # Need to consider auto-roll and auto-yaw
            # Question 1: Can we obtain the alpha on each wing?
            if (CL > 2 or alpha > 20 * (np.pi/180)):
                if X[6] > 0:
                    Cl = 0.01 # In reality, this value is relatively constant & need some calculation
                    if (np.abs(X[6]) < 20 * (np.pi/180)):
                        Cn = 0.0001
                elif X[6] < 0:
                    Cl = -0.01
                    if (np.abs(X[6]) < 20 * (np.pi/180)):
                        Cn = -0.0001
        #--------------------------------------------------------------------


        #--------------------------------------------------------------------
        # Richard's auto-rotation experimentation v2
        whether_to_add_richard_autorotation_v2 = 0
        #--------------------------------------------------------------------
        if whether_to_add_richard_autorotation_v2 == 1:

            global global_CL
            global global_CD
            
            CL_dot = (CL - global_CL) / global_dt
            CD_dot = (CD - global_CD) / global_dt    

            CL_derivative_wrt_alpha = CL_dot / alpha_dot
            CD_derivative_wrt_alpha = CD_dot / alpha_dot      

            # alpha_average_of_left_wing = np.arctan2((X[2]-Xg[2]+(X[3]-Xg[3])*(-(self.FlightData_Geometric_b/4))), (X[0]-Xg[0]))
            # alpha_average_of_right_wing = np.arctan2((X[2]-Xg[2]+(X[3]-Xg[3])*(self.FlightData_Geometric_b/4)), (X[0]-Xg[0]))

            # alpha_average_of_left_wing = np.arctan2((X[2]-Xg[2]+(X[3])*(-(self.FlightData_Geometric_b/4))), (X[0]-Xg[0]))
            # alpha_average_of_right_wing = np.arctan2((X[2]-Xg[2]+(X[3])*(self.FlightData_Geometric_b/4)), (X[0]-Xg[0]))

            alpha_average_of_left_wing = np.arctan2((X[2]+(X[3])*(-(self.FlightData_Geometric_b/4))), (X[0]))
            alpha_average_of_right_wing = np.arctan2((X[2]+(X[3])*(self.FlightData_Geometric_b/4)), (X[0]))

            # By interpreting the results of R05 and T02, we assume 86.4 % of J400's lift come from its main wing (this becomes just a rough approximation if we consider large AoA.)
            CL_wing = CL * (1.6 / 1.85)
            # At slow speeds just after take off and in the initial climb, it is of maximum importance and may produce as much as 70% of total drag
            CD_wing = CD * 0.7

            CL_wing_derivative_wrt_alpha = CL_derivative_wrt_alpha * (1.6 / 1.85)
            CD_wing_derivative_wrt_alpha = CD_derivative_wrt_alpha * 0.7

            alpha_average_of_left_wing_minus_alpha = alpha_average_of_left_wing - alpha
            alpha_average_of_right_wing_minus_alpha = alpha_average_of_right_wing - alpha

            CL_average_of_left_wing = (CL_wing / 2) +  (CL_wing_derivative_wrt_alpha / 2) * alpha_average_of_left_wing_minus_alpha
            CL_average_of_right_wing = (CL_wing / 2) + (CL_wing_derivative_wrt_alpha / 2) * alpha_average_of_right_wing_minus_alpha

            CD_average_of_left_wing = (CD_wing / 2) +  (CD_wing_derivative_wrt_alpha / 2) * alpha_average_of_left_wing_minus_alpha
            CD_average_of_right_wing = (CD_wing / 2) + (CD_wing_derivative_wrt_alpha / 2) * alpha_average_of_right_wing_minus_alpha

            # Some simplified assumptions
            moment_arm_length_of_left_wing_lift_around_cg = self.FlightData_Geometric_b / 4
            moment_arm_length_of_right_wing_lift_around_cg = self.FlightData_Geometric_b / 4

            moment_arm_length_of_left_wing_drag_around_cg = self.FlightData_Geometric_b / 4
            moment_arm_length_of_right_wing_drag_around_cg = self.FlightData_Geometric_b / 4

            # The resulting Cl, Cn increment (or decrement)

            Cl_increment_due_to_autorotation = (((((CL_average_of_left_wing) * (moment_arm_length_of_left_wing_lift_around_cg)) - ((CL_average_of_right_wing) * (moment_arm_length_of_right_wing_lift_around_cg))) * np.cos(alpha) + (((CD_average_of_left_wing) * (moment_arm_length_of_left_wing_drag_around_cg)) - ((CD_average_of_right_wing) * (moment_arm_length_of_right_wing_drag_around_cg))) * np.sin(alpha)) / self.FlightData_Geometric_b)

            Cn_increment_due_to_autorotation = (((((CD_average_of_right_wing) * (moment_arm_length_of_right_wing_drag_around_cg)) - ((CD_average_of_left_wing) * (moment_arm_length_of_left_wing_drag_around_cg))) * np.cos(alpha) + (((CL_average_of_left_wing) * (moment_arm_length_of_left_wing_lift_around_cg)) - ((CL_average_of_right_wing) * (moment_arm_length_of_right_wing_lift_around_cg))) * np.sin(alpha)) / self.FlightData_Geometric_b)

            # The impact on Cl and Cn

            # Check if the numbers are reasonable, especially for the start of the sim. If not, don't modify moments.
            if (np.abs(Cl_increment_due_to_autorotation) < 1) and (np.abs(Cn_increment_due_to_autorotation) < 1):

                Cl = Cl + Cl_increment_due_to_autorotation

                Cn = Cn + Cn_increment_due_to_autorotation
        #--------------------------------------------------------------------


        #--------------------------------------------------------------------
        # Richard's auto-rotation experimentation v3
        whether_to_add_richard_autorotation_v3 = 1
        #--------------------------------------------------------------------
        if whether_to_add_richard_autorotation_v3 == 1:   

            # Calculate the average alpha of left wing and the average alpha of right wing 
            alpha_average_of_left_wing = np.arctan2((X[2]-Xg[2]+(X[3]-Xg[3])*(-(self.FlightData_Geometric_b/4))), (X[0]-Xg[0]))
            alpha_average_of_right_wing = np.arctan2((X[2]-Xg[2]+(X[3]-Xg[3])*(self.FlightData_Geometric_b/4)), (X[0]-Xg[0]))

            # This is an approximation
            CL_average_of_left_wing, poststall_on_left_wing_useless, Cm_left_wing_useless, CD_average_of_left_wing = M10_xp_calculate_CL_Cm_CD(Vt, FlightData_Geometric_c, alpha_average_of_left_wing, alpha_dot, alpha_crit, X_flow_sep, U)
            CL_average_of_right_wing, poststall_on_right_wing_useless, Cm_right_wing_useless, CD_average_of_right_wing = M10_xp_calculate_CL_Cm_CD(Vt, FlightData_Geometric_c, alpha_average_of_right_wing, alpha_dot, alpha_crit, X_flow_sep, U)

            # By interpreting the results of R05 and T02, we assume 86.4 % of J400's lift come from its main wing (this becomes just a rough approximation if we consider large AoA.)
            CL_average_of_left_wing = (1/2) * CL_average_of_left_wing * (1.6 / 1.85)
            CL_average_of_right_wing = (1/2) * CL_average_of_right_wing * (1.6 / 1.85)
            # At slow speeds just after take off and in the initial climb, it is of maximum importance and may produce as much as 70% of total drag
            CD_average_of_left_wing = (1/2) * CD_average_of_left_wing * 0.7
            CD_average_of_right_wing = (1/2) * CD_average_of_right_wing * 0.7

            # Some simplified assumptions
            moment_arm_length_of_left_wing_lift_around_cg = self.FlightData_Geometric_b / 4
            moment_arm_length_of_right_wing_lift_around_cg = self.FlightData_Geometric_b / 4

            moment_arm_length_of_left_wing_drag_around_cg = self.FlightData_Geometric_b / 4
            moment_arm_length_of_right_wing_drag_around_cg = self.FlightData_Geometric_b / 4

            # The resulting Cl, Cn increment (or decrement)

            Cl_increment_due_to_autorotation = (((((CL_average_of_left_wing) * (moment_arm_length_of_left_wing_lift_around_cg)) - ((CL_average_of_right_wing) * (moment_arm_length_of_right_wing_lift_around_cg))) * np.cos(alpha) + (((CD_average_of_left_wing) * (moment_arm_length_of_left_wing_drag_around_cg)) - ((CD_average_of_right_wing) * (moment_arm_length_of_right_wing_drag_around_cg))) * np.sin(alpha)) / self.FlightData_Geometric_b)

            Cn_increment_due_to_autorotation = (((((CD_average_of_right_wing) * (moment_arm_length_of_right_wing_drag_around_cg)) - ((CD_average_of_left_wing) * (moment_arm_length_of_left_wing_drag_around_cg))) * np.cos(alpha) + (((CL_average_of_left_wing) * (moment_arm_length_of_left_wing_lift_around_cg)) - ((CL_average_of_right_wing) * (moment_arm_length_of_right_wing_lift_around_cg))) * np.sin(alpha)) / self.FlightData_Geometric_b)

            # The impact on Cl and Cn

            # Check if the numbers reasonable, especially for the start of the sim. If not, don't modify moments.
            if (np.abs(Cl_increment_due_to_autorotation) < 5) and (np.abs(Cn_increment_due_to_autorotation) < 5):

                Cl = Cl + Cl_increment_due_to_autorotation

                Cn = Cn + Cn_increment_due_to_autorotation
        #--------------------------------------------------------------------


        MomentCoeff = np.array([Cl, Cm, Cn])

        #######################################################################################################
        # Richard Flight Data Recorder Module (place 1 of 3)
        whether_to_add_richard_flight_data_recorder = 1
        #------------------------------------------------------------------------------------------------------
        if whether_to_add_richard_flight_data_recorder == 1:
            # User input field:
            whether_to_add_unused_legacy_variables = 1

            if whether_to_add_unused_legacy_variables == 1:
                # Add additional unused legacy variables
                alpha_crit_2 = M04_calculate_alpha_crit(2, global_time, dX_flow_sep_by_dt, alpha_dot, FlightData_Geometric_c, Vt)
                alpha_crit_3 = M04_calculate_alpha_crit(3, global_time, dX_flow_sep_by_dt, alpha_dot, FlightData_Geometric_c, Vt)
            # Fill in sub data vectors
            global stage2_combined_history
            # stage2_combined_history = np.zeros(14)
            stage2_combined_history[0] = global_time
            stage2_combined_history[1] = global_dt
            stage2_combined_history[2] = dX_flow_sep_by_dt
            stage2_combined_history[3] = X_flow_sep
            stage2_combined_history[4] = alpha
            stage2_combined_history[5] = alpha_crit_1
            stage2_combined_history[6] = alpha_crit_2
            stage2_combined_history[7] = alpha_crit_3
            stage2_combined_history[8] = alpha_dot
            stage2_combined_history[9] = Vt # This is used in precision timing
            stage2_combined_history[10] = CL
            stage2_combined_history[11] = poststall_on
            stage2_combined_history[12] = Cm
            stage2_combined_history[13] = CD

            global stage3_combined_history
            # stage3_combined_history = np.zeros(17)
            stage3_combined_history[0] = dt
            stage3_combined_history[1] = rho
            stage3_combined_history[2] = a_sound_speed
            stage3_combined_history[3] = beta
            stage3_combined_history[4] = Cx
            stage3_combined_history[5] = Cy
            stage3_combined_history[6] = Cz
            stage3_combined_history[7] = Cl
            stage3_combined_history[8] = Cn
            stage3_combined_history[9] = P_prop_power
            stage3_combined_history[10] = n_prop_rot_speed
            stage3_combined_history[11] = global_n_prop_rot_speed
            stage3_combined_history[12] = n_prop_rot_speed_dot
            stage3_combined_history[13] = J_advance_ratio
            stage3_combined_history[14] = eta_p_prop_eff
            stage3_combined_history[15] = ay_noise
            stage3_combined_history[16] = az_noise

            global stage4_combined_history
            # stage4_combined_history = np.zeros(8)
            stage4_combined_history[0] = alpha_average_of_left_wing
            stage4_combined_history[1] = alpha_average_of_right_wing
            stage4_combined_history[2] = CL_average_of_left_wing
            stage4_combined_history[3] = CL_average_of_right_wing
            stage4_combined_history[4] = CD_average_of_left_wing
            stage4_combined_history[5] = CD_average_of_right_wing
            stage4_combined_history[6] = Cl_increment_due_to_autorotation
            stage4_combined_history[7] = Cn_increment_due_to_autorotation

        #######################################################################################################

        return ForceCoeff, MomentCoeff

    """
    ======================================================================================================================
     \brief   Calculates the State Rate Vector of the Aircraft at the given State, with the given force and moments 
              (obtained from aero4560_aero)
              NOTE Usage: Xdot = aero4560_motion(X, ForceCoeff, MomentCoeff, FlightData)
    ======================================================================================================================
    """
    def system_rates(self, X: np.ndarray, ForceCoeff : np.ndarray, MomentCoeff : np.ndarray):
        # Input:
        # X             = state vector
        # ForceCoeff    = force coefficients [Cx; Cy; Cz]
        # MomentCoeff   = moment coefficients [Cl; Cm; Cn]
        # FlightData    = dictionary containing flight data
        # Output:
        # Xdot          = state derivative vector


        ## Shoudn't this be Vt = np.sqrt((X[0]+Xg[0])**2 + (X[1]+Xg[1])**2 + (X[2]+Xg[2])**2), used to be Vt = np.sqrt(X[0]**2 + X[1]**2 + X[2]**2) (outdated line)
        # Richard Note: u=X[0], v=X[1], w=X[2] are aircraft speeds relative to still air in body axes (still air is considered still relative to ground beneath)
        # Richard Note: ug=Xg[0], vg=Xg[1], wg=Xg[2] are wind gust speeds relative to hypothetical "still air" in body axes (still air is considered still relative to ground beneath)
        # Richard Note: pg=Xg[3], qg=Xg[4], rg=Xg[5] are rotational wind gust speeds relative to hypothetical "still air" in body axes (think of them as vortices in turbulent air)
        # Richard Note: For example, ug=Xg[0]= 10 kn roughly represents a 10 kn tail wind
        # Richard Note: According to this formulation, all expressions like X[0]+Xg[0] need sign flips, becoming like X[0]-Xg[0]. Currently no gust model, so not necessary.
        # Richard Note: From a system dynamics perspective, Xdot = A_aero * (X - Xg) + A_kin * X + B * u (if we think in terms of linearly simplified dynamics), this would be A_kin.
        #               Only Vt and qbar are still part of A_aero.
        Vt = np.sqrt((X[0]-Xg[0])**2 + (X[1]-Xg[1])**2 + (X[2]-Xg[2])**2)

        #--------------------------------------------------------
        # Custom atmosphere models (place 2 of 2)
        # User input 
        # 0: original Earth atmosphere model
        # 1: NASA Earth atmosphere model
        # 2: NASA Mars atmosphere model
        which_atmosphere_model_to_use = 1
        #--------------------------------------------------------
        if which_atmosphere_model_to_use == 1: # option 1
            rho, a = M12_nasa_earth_atmosphere_model(-X[11])
        elif which_atmosphere_model_to_use == 2: # option 2
            rho, a = M12_nasa_mars_atmosphere_model(-X[11])
        else: # option 0
            rho, a = self.atomosphere_model(-X[11])
        #--------------------------------------------------------

        qbar = 0.5 * rho * Vt**2
        
        # Calculate Forces and Moments
        F_x = qbar * self.FlightData_Geometric_S * ForceCoeff[0]
        F_y = qbar * self.FlightData_Geometric_S * ForceCoeff[1]
        F_z = qbar * self.FlightData_Geometric_S * ForceCoeff[2]
        M_x = qbar * self.FlightData_Geometric_S * self.FlightData_Geometric_b * MomentCoeff[0]
        M_y = qbar * self.FlightData_Geometric_S * self.FlightData_Geometric_c * MomentCoeff[1]
        M_z = qbar * self.FlightData_Geometric_S * self.FlightData_Geometric_b * MomentCoeff[2]
        
        # Calculate Linear Acceleration
        Xdot = np.zeros(12)
        Xdot[0] = -self.FlightData_Inertial_g * np.sin(X[7]) + X[5] * X[1] - X[4] * X[2] + F_x / self.FlightData_Inertial_m
        Xdot[1] = self.FlightData_Inertial_g * np.sin(X[6]) * np.cos(X[7]) - X[5] * X[0] + X[3] * X[2] + F_y / self.FlightData_Inertial_m
        Xdot[2] = self.FlightData_Inertial_g * np.cos(X[6]) * np.cos(X[7]) + X[4] * X[0] - X[3] * X[1] + F_z / self.FlightData_Inertial_m
        
        # # Calculate Inertial Coefficients
        c0 = self.FlightData_Inertial_Ixx * self.FlightData_Inertial_Izz - self.FlightData_Inertial_Ixz**2
        c1 = self.FlightData_Inertial_Izz / c0
        c2 = self.FlightData_Inertial_Ixz / c0
        c3 = c2 * (self.FlightData_Inertial_Ixx - self.FlightData_Inertial_Iyy + self.FlightData_Inertial_Izz)
        c4 = c1 * (self.FlightData_Inertial_Iyy - self.FlightData_Inertial_Izz) - c2 * self.FlightData_Inertial_Ixz
        c5 = 1 / self.FlightData_Inertial_Iyy
        c6 = c5 * self.FlightData_Inertial_Ixz
        c7 = c5 * (self.FlightData_Inertial_Izz - self.FlightData_Inertial_Ixx)
        c8 = self.FlightData_Inertial_Ixx / c0
        c9 = c8 * (self.FlightData_Inertial_Ixx - self.FlightData_Inertial_Iyy) + c2 * self.FlightData_Inertial_Ixz
        
        # # Calculate Rotation Acceleration
        Xdot[3] = c3 * X[3] * X[4] + c4 * X[4] * X[5] + c1 * M_x + c2 * M_z
        Xdot[4] = c7 * X[3] * X[5] - c6 * (X[3]**2 - X[5]**2) + c5 * M_y
        Xdot[5] = c9 * X[3] * X[4] - c3 * X[4] * X[5] + c2 * M_x + c8 * M_z
        
        # Calculate Euler Rates
        Xdot[6] = X[3] + X[4] * np.sin(X[6]) * np.tan(X[7]) + X[5] * np.cos(X[6]) * np.tan(X[7])
        Xdot[7] = X[4] * np.cos(X[6]) - X[5] * np.sin(X[6])
        Xdot[8] = X[4] * np.sin(X[6]) / np.cos(X[7]) + X[5] * np.cos(X[6]) / np.cos(X[7])
        
        # Calculate Lbe (now C_eb, the correct naming, because the convention of the lower index is from right to left, which is the conventional direction in which we multiply vector by matrix)
        X = X.flatten()

        # Rotation matrix from body axes to earth axes, C_eb = (C_x(phi) Cy(theta) Cz(psi))^T (transpose operator), the result by hand calculation is the following.
        C_eb = np.array([[np.cos(X[8]) * np.cos(X[7]), np.cos(X[8]) * np.sin(X[7]) * np.sin(X[6]) - np.sin(X[8]) * np.cos(X[6]), np.cos(X[8]) * np.sin(X[7]) * np.cos(X[6]) + np.sin(X[8]) * np.sin(X[6])],
                        [np.sin(X[8]) * np.cos(X[7]), np.sin(X[8]) * np.sin(X[7]) * np.sin(X[6]) + np.cos(X[8]) * np.cos(X[6]), np.sin(X[8]) * np.sin(X[7]) * np.cos(X[6]) - np.cos(X[8]) * np.sin(X[6])],
                        [-np.sin(X[7]), np.cos(X[7]) * np.sin(X[6]), np.cos(X[7]) * np.cos(X[6])]])

        # # Calculate Position Rates
        C_eb = C_eb.reshape((3, 3))
        Vel = np.array([X[0], X[1], X[2]])
        Xdot[9:12] = C_eb @ Vel
        

        return Xdot

    """
    ======================================================================================================================
     \brief   Calculates the ISA density and speed of sound at the given altitude
              NOTE Usage: rho, a = aero4560_atmos(alt)
    ======================================================================================================================
    """
    def atomosphere_model(self, alt):
        # Input:
        # alt   = altitude in meters
        # Output:
        # rho   = air density in kg/m^3
        # a     = speed of sound in m/s

        lap_rate = 0.0065  # (C)deg/meter
        t0 = 15            # (C)deg
        rho0 = 1.225       # Kg/m^3
        p0 = 101310        # Pa

        if alt >= 11000:
            temp = -56.4
            pres = p0 * (0.2189 * (np.exp(9.81 * (11000 - alt) / (287.1 * (temp + 273)))))
            rho = rho0 * (0.2972 * (np.exp(9.81 * (11000 - alt) / (287.1 * (temp + 273)))))
        else:
            temp = (t0 - (alt * lap_rate))
            pres = p0 * (((temp + 273) / (t0 + 273)) ** (9.81 / (lap_rate * 287.1)))
            rho = rho0 * (((temp + 273) / (t0 + 273)) ** ((9.81 - (lap_rate * 287.1)) / (lap_rate * 287.1)))
        
        a = (1.4 * 287.1 * (temp + 273)) ** 0.5

        return rho, a

    """
    ======================================================================================================================
     \brief   Function to perform simple Euler integration of state equations for generic flight dynamic model
              NOTE Usage: X = aero4560_euler(dt, X, Xg, U, FlightData)
    ======================================================================================================================
    """
    def euler_integrate(self, dt: float, X : np.ndarray, Xg : np.ndarray, U : np.ndarray):
        ''''
        Euler integrates the flight model 

        Parameters:
        ----------
        dt : float
            Time step for integration
        X : np.ndarray (12,)
            State vector of the flight model
        Xg : np.ndarray (12,)
            State vector of the Gust 
        U : np.ndarray (5,)
            Control vector of the flight model
        
        Returns:
        ----------
        X_out : np.ndarray (12,)
            State vector of the flight model after integration
        
        Example:
        ----------
        >>> fm = flightModel()

        >>> X_out = fm.euler_integrate(dt, X, Xg, U)
        '''

        X = X.flatten()
        Xg = Xg.flatten()
        U = U.flatten()

        U = np.clip(U, self.FlightData_ControlLimits_Lower, self.FlightData_ControlLimits_Upper)
        xp.log("Euler Integration Step"); xp.log(f"Control Vector Before Tweak: {U}"); U[0] = (U[0] - 0.4375) * 8; xp.log(f"Control Vector After Tweak: {U}"); # Richard generating custom log in XPPython3Log; # Richard temporary fixing U control vector bug
        if not self.checkIntegratable(X, Xg, U):
            # If the flight model is not integratable, return the state after 
            # increasing the velocity by a small amount
            ##### HACK: THIS IS A FLIGHT MODEL HACK #####
            # NOTE: This is a hacky solution to avoid the flight model from
            #       crashing when the velocity is zero. It works but might not be 
            #       the best solution .
            ##### HACK: THIS IS A FLIGHT MODEL HACK #####
            X[0:3] += 1e-1
            return X

        Xdot = self.iter_solve(X, Xg, U)
        
        # Integrate states
        X += Xdot * dt

        # manually enforce t = 0.01 s (not guarateed to be the correct approach)
        # dt = 0.01
        # ===================================================================
        # Richard Additional Code Block 2 (Start)
        # ===================================================================
        # Update X_flow_sep as well
        global X_flow_sep  # Declare that we will use the global variable
        global dX_flow_sep_by_dt
        global global_time
        X_flow_sep += dX_flow_sep_by_dt * dt  # Modify the global variable
        global_time += dt # Let the global clock tick

        # # check the value of dt for euler integration
        # with open('richard_checked_terms\\dt_euler_integration.txt', 'a') as file:
        #     file.write(str(dt) + '\n')
        # with open('richard_checked_terms\\dt.txt', 'a') as file:
        #     file.write(str(dt) + '\n')

        # use global_dt to record dt in euler integration
        global global_dt
        global_dt = dt
        # ===================================================================
        # Richard Additional Code Block 2 (End)
        # ===================================================================

        return X
    
    """
    ======================================================================================================================
     \brief   Function to perform midpoint integration of state equations for generic flight dynamic model
              NOTE Usage: X = midpoint_integration(dt, X, Xg, U, FlightData)
    ======================================================================================================================
    """
    def midpoint_integration(self, dt: float, X : np.ndarray, Xg : np.ndarray, U : np.ndarray):
        ''''
        RK2 integrate the flight model 

        Parameters:
        ----------
        dt : float
            Time step for integration
        X : np.ndarray (12,)
            State vector of the flight model
        Xg : np.ndarray (12,)
            State vector of the Gust 
        U : np.ndarray (5,)
            Control vector of the flight model
        
        Returns:
        ----------
        X_out : np.ndarray (12,)
            State vector of the flight model after integration
        
        Example:
        ----------
        >>> fm = flightModel()

        >>> X_out = fm.midpoint_integrate(dt, X, Xg, U)
        '''

        X = X.flatten()
        Xg = Xg.flatten()
        U = U.flatten()

        U = np.clip(U, self.FlightData_ControlLimits_Lower, self.FlightData_ControlLimits_Upper)
        xp.log("Midpoint Integration Step:"); xp.log(f"Control Vector Before Tweak: {U}"); U[0] = (U[0] - 0.4375) * 8; xp.log(f"Control Vector After Tweak: {U}"); # Richard generating custom log in XPPython3Log; # Richard temporary fixing U control vector bug

        #------------------------------------------------------------------------------------------------------------------------------------------------------------
        # Declare this important global variable first (for Richard Taped Control Input Module, Richard Taped Gust Module, Richard Flight Data Recorder Module)
        global data_recording_current_duration
        #------------------------------------------------------------------------------------------------------------------------------------------------------------



        # ##############################################################################################################################
        # # # Richard Flap Experimentation (place 1 of 2) (should be useless now)

        # # U from last time step
        # global global_U 
        # global global_flap_is_going_down
        # global global_flap_is_going_up

        # # if flap deflection is -0 deg last time step and -20 this time step (means flap down is pressed)
        # if (global_U[4] == np.deg2rad(-0)) and (U[4] == np.deg2rad(-20)):
        #     U[4] = 0
        #     global_flap_is_going_down = True
        # elif (global_U[4] == np.deg2rad(-19)) and (U[4] == np.deg2rad(-20)):
        #     global_flap_is_going_down = False

        # # if flap deflection is -20 deg last time step and -40 this time step (means flap down is pressed)
        # elif (global_U[4] == np.deg2rad(-20)) and (U[4] == np.deg2rad(-40)):
        #     U[4] = -20
        #     global_flap_is_going_down = True
        # elif (global_U[4] == np.deg2rad(-39)) and (U[4] == np.deg2rad(-40)):
        #     global_flap_is_going_down = False

        # # if flap deflection is -40 deg last time step and -20 this time step (means flap up is pressed)
        # elif (global_U[4] == np.deg2rad(-40)) and (U[4] == np.deg2rad(-20)):
        #     U[4] = -40
        #     global_flap_is_going_up = True
        # elif (global_U[4] == np.deg2rad(-21)) and (U[4] == np.deg2rad(-20)):
        #     global_flap_is_going_up = False

        # # if flap deflection is -20 deg last time step and -0 this time step (means flap up is pressed)
        # elif (global_U[4] == np.deg2rad(-20)) and (U[4] == np.deg2rad(-0)):
        #     U[4] = -20
        #     global_flap_is_going_up = True
        # elif (global_U[4] == np.deg2rad(-1)) and (U[4] == np.deg2rad(-0)):
        #     global_flap_is_going_up = False

        # # Check whether flap is going up, going down, or idle (not mentioned)
        # if global_flap_is_going_down == True:
        #     U[4] = U[4] - np.deg2rad(1)

        # elif global_flap_is_going_up == True:
        #     U[4] = U[4] + np.deg2rad(1)

        # global_U = U
        # ##############################################################################################################################

        ##############################################################################################################################
        # Richard Taped Control Input Module (place 1 of 3) # Relies on the recorded time by Richard Flight Data Recorder Module
        whether_to_add_richard_taped_control_inputs = 0
        #------------------------------------------------------------------------------------------------------
        # Custom input U (for remote testing) (can be pre-recorded with R32_extract...) (place 1 of 3)
        if whether_to_add_richard_taped_control_inputs == 1:
            global U_data_matrix_taped
            global U_data_matrix_taped_row_index
            global t_taped_for_U
            global interpolated_U_taped
            # global data_recording_current_duration 

            # User input
            whether_to_use_timewise_accurate_taped_control_inputs = 1

            if whether_to_use_timewise_accurate_taped_control_inputs == 1:
                # Taped U vs t formulation by interpolation!
                if t_taped_for_U[1] < data_recording_current_duration < t_taped_for_U[-1]:
                    U = interpolated_U_taped(data_recording_current_duration)
            else:
                # U_data_matrix = np.load('Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\U_data_matrix.npy')
                # Check if the row index is within the valid range
                if 0 <= U_data_matrix_taped_row_index < U_data_matrix_taped.shape[0]:
                    U = U_data_matrix_taped[U_data_matrix_taped_row_index, :]
        ##############################################################################################################################


        ##############################################################################################################################
        # Richard Taped Gust Module (place 1 of 2)
        whether_to_add_richard_taped_gust_module = 1
        #------------------------------------------------------------------------------------------------------
        if whether_to_add_richard_taped_gust_module == 1:
            global Xg_data_matrix_taped
            # global Xg_data_matrix_taped_row_index
            global t_taped_for_Xg
            global interpolated_Xg_taped
            # if whether_to_add_richard_taped_control_inputs != 1:
            #     global data_recording_current_duration 

            # Taped Xg vs t formulation by interpolation!
            # Run the track for the 1st time
            if t_taped_for_Xg[1] < data_recording_current_duration < t_taped_for_Xg[-1]:
                Xg = interpolated_Xg_taped(data_recording_current_duration)

            # User input
            # If repeats are not necessary, mute the following lines.
            # Repeat the track the 1st time
            elif t_taped_for_Xg[1] < (data_recording_current_duration - t_taped_for_Xg[-1]) < t_taped_for_Xg[-1]:
                Xg = interpolated_Xg_taped(data_recording_current_duration - t_taped_for_Xg[-1])

            # Repeat the track the 2nd time
            elif t_taped_for_Xg[1] < (data_recording_current_duration - 2 * t_taped_for_Xg[-1]) < t_taped_for_Xg[-1]:
                Xg = interpolated_Xg_taped(data_recording_current_duration - 2 * t_taped_for_Xg[-1])

            # Repeat the track the 3rd time
            elif t_taped_for_Xg[1] < (data_recording_current_duration - 3 * t_taped_for_Xg[-1]) < t_taped_for_Xg[-1]:
                Xg = interpolated_Xg_taped(data_recording_current_duration - 3 * t_taped_for_Xg[-1])

            # Repeat the track the 4th time
            elif t_taped_for_Xg[1] < (data_recording_current_duration - 4 * t_taped_for_Xg[-1]) < t_taped_for_Xg[-1]:
                Xg = interpolated_Xg_taped(data_recording_current_duration - 4 * t_taped_for_Xg[-1])

            # Repeat the track the 5th time
            elif t_taped_for_Xg[1] < (data_recording_current_duration - 5 * t_taped_for_Xg[-1]) < t_taped_for_Xg[-1]:
                Xg = interpolated_Xg_taped(data_recording_current_duration - 5 * t_taped_for_Xg[-1])
                
        ##############################################################################################################################


        ##############################################################################################################################
        # Richard Real Time Gust Module (place 1 of 2)
        whether_to_add_richard_real_time_gust_module = 0
        #------------------------------------------------------------------------------------------------------
        if whether_to_add_richard_real_time_gust_module == 1:

            global interpolate_stan_dev_u
            global interpolate_stan_dev_v
            global interpolate_stan_dev_w
            global interpolate_stan_dev_p
            global interpolate_stan_dev_q
            global interpolate_stan_dev_r

            # User input
            whether_to_use_altitude_dependent_gust = 0

            # Acquire current aircraft altitude
            current_aircraft_altitude = -X[11]

            try:
                # Try to print the value of the global variable
                print(current_aircraft_altitude)
            except NameError:
                # If variable has not been defined, catch the NameError exception
                # Assign the value 0.01 to 'global_dt' since it was not previously defined
                current_aircraft_altitude = 640.350647 # m (for 10 nm approach)

            if whether_to_use_altitude_dependent_gust == 1 and (-10 < current_aircraft_altitude < 45720):

                # Gust Xg to be calculated with Dryden spectra (alt dependent)
                u_g_mean = 0.0
                u_g_std = interpolate_stan_dev_u(current_aircraft_altitude)
                u_g = np.random.normal(u_g_mean, u_g_std)

                v_g_mean = 0.0
                v_g_std = interpolate_stan_dev_v(current_aircraft_altitude)
                v_g = np.random.normal(v_g_mean, v_g_std)

                w_g_mean = 0.0
                w_g_std = interpolate_stan_dev_w(current_aircraft_altitude)
                w_g = np.random.normal(w_g_mean, w_g_std)

                p_g_mean = 0.0
                p_g_std = interpolate_stan_dev_p(current_aircraft_altitude)
                p_g = np.random.normal(p_g_mean, p_g_std)

                q_g_mean = 0.0
                q_g_std = interpolate_stan_dev_q(current_aircraft_altitude)
                q_g = np.random.normal(q_g_mean, q_g_std)

                r_g_mean = 0.0
                r_g_std = interpolate_stan_dev_r(current_aircraft_altitude)
                r_g = np.random.normal(r_g_mean, r_g_std)

            else:

                # Gust Xg to be calculated with Dryden spectra (Currently alt independent)
                u_g_mean = 0.0
                u_g_std = 6.7271452489e-01
                u_g = np.random.normal(u_g_mean, u_g_std)

                v_g_mean = 0.0
                v_g_std = 6.7272174542e-01
                v_g = np.random.normal(v_g_mean, v_g_std)

                w_g_mean = 0.0
                w_g_std = 8.7071167904e-01
                w_g = np.random.normal(w_g_mean, w_g_std)

                p_g_mean = 0.0
                p_g_std = 2.3900903098e-02
                p_g = np.random.normal(p_g_mean, p_g_std)

                q_g_mean = 0.0
                q_g_std = 1.2983696230e-02
                q_g = np.random.normal(q_g_mean, q_g_std)

                r_g_mean = 0.0
                r_g_std = 1.6944854495e-02
                r_g = np.random.normal(r_g_mean, r_g_std)

            Xg[0] = u_g
            Xg[1] = v_g
            Xg[2] = w_g
            Xg[3] = p_g
            Xg[4] = q_g
            Xg[5] = r_g
        ##############################################################################################################################


        ##############################################################################################################################
        # Richard Directional Wind Module (place 1 of 1)
        whether_to_add_richard_directional_wind_module = 0
        #------------------------------------------------------------------------------------------------------
        if whether_to_add_richard_directional_wind_module == 1:

            # Rotation matrix from body axes to earth axes (custom physics engine uses NED), C_eb = (C_x(phi) Cy(theta) Cz(psi))^T (transpose operator), the result by hand calculation is the following.
            C_eb = np.array([[np.cos(X[8]) * np.cos(X[7]), np.cos(X[8]) * np.sin(X[7]) * np.sin(X[6]) - np.sin(X[8]) * np.cos(X[6]), np.cos(X[8]) * np.sin(X[7]) * np.cos(X[6]) + np.sin(X[8]) * np.sin(X[6])],
                            [np.sin(X[8]) * np.cos(X[7]), np.sin(X[8]) * np.sin(X[7]) * np.sin(X[6]) + np.cos(X[8]) * np.cos(X[6]), np.sin(X[8]) * np.sin(X[7]) * np.cos(X[6]) - np.cos(X[8]) * np.sin(X[6])],
                            [-np.sin(X[7]), np.cos(X[7]) * np.sin(X[6]), np.cos(X[7]) * np.cos(X[6])]])

            # Add a reshape just to be safe (and redundant)
            C_eb = C_eb.reshape((3, 3))
            # Transpose C_eb to obtain C_be (the inverse of a direct cosine martix (or any rotation) is its transpose)
            C_be = C_eb.T
            # # Rotation matrix from NED (conventional North-East-Down) axes to SWD (South-West-Down) axes, simply flip the X and Y coordinates will do.
            # C_SWD_NED = np.array([[-1, 0, 0],
            #                       [0, -1, 0],
            #                       [0, 0, 1]])
            # # Add a reshape just to be safe (and redundant)
            # C_SWD_NED = C_SWD_NED.reshape((3, 3))

            k_karman_constant = 0.4 # Karman constant
            H_0_surface_roughness_height = 3 # m
            V_a0_friction_velocity = 0.4 # m/s

            # Acquire aircraft altitude
            aircraft_altitude = -X[11]
            local_average_terrain_altitude = 0 # Assume local average terrain altitude is 0 m for the moment (in the future, could be obtained from some kind of DataRef in XPlane).
            aircraft_altitude_relative_to_local_average_terrain = aircraft_altitude - local_average_terrain_altitude

            # Calculate V_a
            # Below ~ 3 m, the average wind speed is effectively 0 due to surface roughness blockage.
            # Above 40000 m, Earth air density rho ~ 0, aero effect is too weak to be accounted for.
            if (aircraft_altitude_relative_to_local_average_terrain > H_0_surface_roughness_height) and (aircraft_altitude_relative_to_local_average_terrain < 40000): 
                V_a_average_wind_speed = (V_a0_friction_velocity / k_karman_constant) * np.log(aircraft_altitude_relative_to_local_average_terrain / H_0_surface_roughness_height) # m/s
            else: 
                V_a_average_wind_speed = 0 # m/s

            # User input
            # Set wind direction:
            # psi_wind = 0 is blowing towards north (which is likely magnetic north, slightly different from true north)
            # psi_wind = 90 is blowing towards east
            # psi_wind = 180 is blowing towards south
            # psi_wind = 270 is blowing towards west

            # X-Plane Datarefs:

            # sim/flightmodel/position/psi
            # "The true heading of the aircraft in degrees from the Z axis - OpenGL coordinates" (Richard correction: "from the -Z axis")

            # sim/flightmodel/position/true_psi
            # "The heading of the aircraft relative to the earth precisely below the aircraft - true degrees north, always"

            # We try a wind blowing towards the magnetic north (aka "South Wind")
            psi_wind = np.deg2rad(0)

            # Directional wind component in North-East-Down axis system (local "earth axis")
            northward_wind_component_in_mps = V_a_average_wind_speed * np.cos(psi_wind) # m/s # negative is southward
            eastward_wind_component_in_mps = V_a_average_wind_speed * np.sin(psi_wind) # m/s # negative is westward
            downward_wind_component_in_mps = 0 # m/s # negative is upward

            directional_wind_in_NED_axes = np.array([northward_wind_component_in_mps, eastward_wind_component_in_mps, downward_wind_component_in_mps])

            # # Convert from NED (conventional North-East-Down) axes to SWD (South-West-Down) axes, using a rotation matrix.
            # directional_wind_in_SWD_axes = C_SWD_NED @ directional_wind_in_NED_axes

            Xg_directional_wind_component = C_be @ directional_wind_in_NED_axes

            # Add the direction wind to u_g, v_g, w_g
            Xg[0:3] = Xg[0:3] + Xg_directional_wind_component

        ##############################################################################################################################


        if not self.checkIntegratable(X, Xg, U):
            # If the flight model is not integratable, return the state after 
            # increasing the velocity by a small amount
            ##### HACK: THIS IS A FLIGHT MODEL HACK #####
            # NOTE: This is a hacky solution to avoid the flight model from
            #       crashing when the velocity is zero. It works but might not be 
            #       the best solution .
            ##### HACK: THIS IS A FLIGHT MODEL HACK #####

            X[0:3] += 1e-1
            return X

        Xdot = np.zeros(12)
        for _ in range(100):
            ## THIS IS AN ALTERNATIVE IMPLEMENTATION OF THE MIDPOINT INTEGRATION
            ForceCoeff, MomentCoeff = self.calculate_aero_forces(X, Xg, Xdot, U)
            ydash = self.system_rates(X, ForceCoeff, MomentCoeff)
            ForceCoeff, MomentCoeff = self.calculate_aero_forces(X + dt/2 * ydash, Xg, Xdot, U)
            Xdot_new = self.system_rates(X + dt/2 * ydash, ForceCoeff, MomentCoeff)
            ## END OF ALTERNATIVE IMPLEMENTATION

            ## THIS IS THE ORIGINAL IMPLEMENTATION OF THE MIDPOINT INTEGRATION
            # ydash = self.iter_solve(X, Xg, U)
            # Xdot_new = self.iter_solve(X + dt/2 * ydash, Xg, U)
            ## END OF ORIGINAL IMPLEMENTATION

            err = Xdot_new - Xdot
            err = err.dot(err)**0.5
            if err < 1e-2:
                break
            Xdot = Xdot_new
        
        # manually enforce t = 0.01 s (not guarateed to be the correct approach)
        # dt = 0.01
        # ===================================================================
        # Richard Additional Code Block 3 (Start)
        # ===================================================================
        # Update X_flow_sep as well
        global X_flow_sep  # Declare that we will use the global variable
        global dX_flow_sep_by_dt
        global global_time
        # use global_dt to record dt in midpoint integration
        global global_dt


        ##############################################################################################################################
        # Richard Precision Timing Module (place 1 of 1)
        Vt = stage2_combined_history[9] #Vt # This is used in precision timing # Vt need to be speed relative to ground (check)
        c_speed_of_light = 299792458 # m/s
        # SR frequency shift

        # GR frequency shift

        # Sagnac shift

        # This module is developed in M13
        # But decided to only include in data analysis section R32 (since it won't affect the aircraft behavior in the simulation)

        ### This module is ultimately not run in real time but is incorporated in R32 post processing of flight data ***
        ##############################################################################################################################
 

        #######################################################################################################
        # Richard Flight Data Recorder Module (place 2 of 3)
        whether_to_add_richard_flight_data_recorder = 1
        #------------------------------------------------------------------------------------------------------
        if whether_to_add_richard_flight_data_recorder == 1:
            ## Initialize an empty list to store data vectors
            # flight_data_matrix_as_list = []
            global flight_data_matrix_as_list
            global whether_flight_data_recorder_is_on

            # if (whether_to_add_richard_taped_control_inputs != 1) and (whether_to_add_richard_taped_gust_module != 1):
            #     global data_recording_current_duration 
            # global data_recording_current_duration

            # global global_time
            # global global_dt

            # At which moment in time to stop recording data? (Record data for how long since pressing '9'?)
            global_time_to_stop_recording_data = 50 # sec

            if len(flight_data_matrix_as_list) > 2:
                data_recording_start_time = flight_data_matrix_as_list[2][0]
                data_recording_current_duration = global_time - data_recording_start_time
            else:
                data_recording_current_duration = global_time

            # declare it here just to be on the safe side (should be declared inside calculate_aero_forces() function calls)
            # global stage2_combined_history

            # Create a new flight_data_row each iteration to avoid overwriting previous data
            flight_data_row = np.zeros(81) # 1 + 12 + 12 + 5 + 14 + 17 + 12 + 8
            flight_data_row[0] = data_recording_current_duration
            flight_data_row[1:13] = X[0:12]  # Adjusted to fill the entire row
            flight_data_row[13:25] = Xdot[0:12]
            flight_data_row[25:30] = U[0:5]
            flight_data_row[30:44] = stage2_combined_history[0:14]
            flight_data_row[44:61] = stage3_combined_history[0:17]
            flight_data_row[61:73] = Xg[0:12]
            flight_data_row[73:81] = stage4_combined_history[0:8]
            xp.log(f"flight_data_row: {flight_data_row}")
            flight_data_matrix_as_list.append(flight_data_row)
            xp.log(f"flight_data_matrix_as_list Length: {len(flight_data_matrix_as_list)}")

            # if len(flight_data_matrix_as_list) > 1000 and len(flight_data_matrix_as_list) < 1100 and whether_flight_data_recorder_is_on == 1:
            if data_recording_current_duration > global_time_to_stop_recording_data and data_recording_current_duration < (global_time_to_stop_recording_data+1) and whether_flight_data_recorder_is_on == 1:
                # Convert the list of vectors to a NumPy array (matrix)
                flight_data_matrix_as_list[2][0] = 0 # fix the starting time of recording
                flight_data_matrix = np.array(flight_data_matrix_as_list)
                flight_data_matrix_without_first_two_rows = flight_data_matrix[2:,:]
                # np.save('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\flight_data_matrix.npy', flight_data_matrix_without_first_two_rows)
                # np.save('C:\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\flight_data_matrix.npy', flight_data_matrix_without_first_two_rows) # Motion Chair 2
                np.save('Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\flight_data_matrix.npy', flight_data_matrix_without_first_two_rows) # Relative path approach that works on all machines
                whether_flight_data_recorder_is_on = 0
        #######################################################################################################

        ##############################################################################################################################
        # Richard Taped Control Input Module (place 2 of 3)
        # No need for user input
        if whether_to_add_richard_taped_control_inputs == 1:
            # Custom input U (for remote testing) (can be pre-recorded with R32_extract...) (place 2 of 3)
            if whether_to_use_timewise_accurate_taped_control_inputs != 1:
                # global U_data_matrix
                # global U_data_matrix_taped_row_index
                # U_data_matrix = np.load('Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\U_data_matrix.npy')
                # U = U_data_matrix[U_data_matrix_row_index, :]
                U_data_matrix_taped_row_index = U_data_matrix_taped_row_index + 1 # proceed to the next U input vector
        ##############################################################################################################################

        X_flow_sep += dX_flow_sep_by_dt * dt  # Modify the global variable
        global_time += dt # Let the global clock tick

        # # check the value of dt for midpoint integration
        # with open('richard_checked_terms\\dt_midpoint_integration.txt', 'a') as file:
        #     file.write(str(dt) + '\n')
        # with open('richard_checked_terms\\dt.txt', 'a') as file:
        #     file.write(str(dt) + '\n')
        global_dt = dt
        # global global_flight_recorder
        # global_flight_recorder[0:11,] = 
        # ===================================================================
        # Richard Additional Code Block 3 (End)
        # ===================================================================

        X += Xdot*dt

        ##############################################################################################################################
        # Richard Crash Prevention Module (place 1 of 1)
        whether_to_add_richard_crash_prevention_module = 1
        #------------------------------------------------------------------------------------------------------
        # If we turn on this module, don't press '9' when altitude is below 100 m.
        if whether_to_add_richard_crash_prevention_module == 1:

            if (-X[11]) < 100:

                # Rotation matrix from body axes to earth axes, C_eb = (C_x(phi) Cy(theta) Cz(psi))^T (transpose operator), the result by hand calculation is the following.
                C_eb = np.array([[np.cos(X[8]) * np.cos(X[7]), np.cos(X[8]) * np.sin(X[7]) * np.sin(X[6]) - np.sin(X[8]) * np.cos(X[6]), np.cos(X[8]) * np.sin(X[7]) * np.cos(X[6]) + np.sin(X[8]) * np.sin(X[6])],
                                [np.sin(X[8]) * np.cos(X[7]), np.sin(X[8]) * np.sin(X[7]) * np.sin(X[6]) + np.cos(X[8]) * np.cos(X[6]), np.sin(X[8]) * np.sin(X[7]) * np.cos(X[6]) - np.cos(X[8]) * np.sin(X[6])],
                                [-np.sin(X[7]), np.cos(X[7]) * np.sin(X[6]), np.cos(X[7]) * np.cos(X[6])]])

                # # Calculate Position Rates
                C_eb = C_eb.reshape((3, 3))

                # Vel before the bounce
                Vel = np.array([X[0], X[1], X[2]])
                Xdot[9:12] = C_eb @ Vel

                # # Bounce the Z component of velocity vector in NED axes
                # if Xdot[11] > 0: # if the aircraft is descending
                #     Xdot[11] = - Xdot[11]
                #     if Xdot[11] < 3: # if the descend velocity component is less than 3 m/s
                #         Xdot[11] = 0

                # Bounce the Z component of velocity vector in NED axes
                Xdot[11] = - Xdot[11]

                # Vel after the bounce
                C_be = C_eb.T
                Vel = C_be @ Xdot[9:12]

                # Calculate Energy dissipation
                Kinetic_Energy_after_bounce_to_Kinetic_Energy_before_bounc_ratio = 1.2 # Create an unphysical bounce to avoid crashing # (previous version) assume 50% kinetic energy loss in each bounce
                Vel_magnitude_after_bounce_to_Vel_magnitude_before_bounce = np.sqrt(Kinetic_Energy_after_bounce_to_Kinetic_Energy_before_bounc_ratio)

                # Reduce the magnitude of the velocity after the bounce
                Vel_magnitude = np.sqrt(Vel[0]**2 + Vel[1]**2 + Vel[2]**2)
                if Vel_magnitude < 61.73: # if |Vel| < 120 kn
                    Vel = Vel * Vel_magnitude_after_bounce_to_Vel_magnitude_before_bounce
                
                # This will bounce the aircraft
                X[0:3] = Vel

        ##############################################################################################################################

        return X
    
    def iter_solve(self, X: np.ndarray, Xg: np.ndarray, U : np.ndarray):
        '''
        Iteratively solve the flight model
        '''
        Xdot = np.zeros(12)
        for _ in range(100):
            ForceCoeff, MomentCoeff = self.calculate_aero_forces(X, Xg, Xdot, U)
            Xdot_new = self.system_rates(X, ForceCoeff, MomentCoeff)
            err = Xdot_new - Xdot
            err = err.dot(err)**0.5
            if err < 1e-2:
                break
            Xdot = Xdot_new
        return Xdot
    
    def checkIntegratable (self, X : np.ndarray, Xg : np.ndarray, U : np.ndarray) -> bool:
        ''''
        Check if the flight model is integratable 

        Parameters:
        ----------
        X : np.ndarray (12,)
            State vector of the flight model
        Xg : np.ndarray (12,)
            State vector of the Gust 
        U : np.ndarray (5,)
            Control vector of the flight model
        
        Returns:
        ----------
        integratable : bool
            True if integratable, False otherwise
        
        Example:
        ----------
        >>> fm = flightModel()

        >>> integratable = fm.checkIntegratable(X, Xg, U)
        '''

        Vt = np.sqrt((X[0]-Xg[0])**2 + (X[1]-Xg[1])**2 + (X[2]-Xg[2])**2)

        if Vt < 1e-1:
            return False
        
        return True

# ===================================================================
# Richard Additional Code Block 4 (Start)
# ===================================================================
# Define and initialize global variables 
# (if any global variable is not initiated, Python will not run)
X_flow_sep = 1
dX_flow_sep_by_dt = 0
global_alpha = np.arctan2((1+1), (1+1))
global_time = 0

global_n_prop_rot_speed = 35 # rev/s (2100 RPM)

global_q_pitch_rate = 0 # rad/s
global_r_yaw_rate = 0 # rad/s

whether_to_use_full_formulation = 0 # Warning: full formulation is more accurate but much slower for the sim to run (Place 2 of 2)
if whether_to_use_full_formulation == 1:
    # Initialize the BuffetNoiseGenerator
    buffet_noise_generator = BuffetNoiseGenerator(fs=100, series_duration=0.5)

#---------------------------------
# For autorotation
global_CL = 0.3
global_CD = 0.05
#---------------------------------

#######################################################################################################
# Richard Flight Data Recorder Module (place 3 of 3)
# Initialize an empty list to store data vectors
whether_to_add_richard_flight_data_recorder = 1
#------------------------------------------------------------------------------------------------------
if whether_to_add_richard_flight_data_recorder == 1:
    flight_data_matrix_as_list = []
    whether_flight_data_recorder_is_on = 1 # initialize with the flight data recorder on
    stage2_combined_history = np.zeros(14)
    data_recording_current_duration = 0 # Shared with "Custom input U"

    stage3_combined_history = np.zeros(17)
    stage4_combined_history = np.zeros(8)
#######################################################################################################

##############################################################################################################################
# Richard Taped Control Input Module (place 3 of 3)
whether_to_add_richard_taped_control_inputs = 0
#------------------------------------------------------------------------------------------------------
if whether_to_add_richard_taped_control_inputs == 1:
    # Custom input U (for remote testing) (can be pre-recorded with R32_extract...) (place 3 of 3)
    U_data_matrix_taped = np.load('Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\U_data_matrix_taped.npy')
    t_taped_for_U = np.load('Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\t_taped_for_U.npy')
    U_data_matrix_taped_row_index = 0
    # Create an interpolation function for each column in U_data_matrix
    interpolated_U_taped = interp1d(t_taped_for_U, U_data_matrix_taped, axis=0, kind='linear', fill_value='extrapolate')
    data_recording_current_duration = 0 # Borrow a variable from "Richard Flight Data Recorder Module"
##############################################################################################################################


##############################################################################################################################
# Richard Taped Gust Module (place 2 of 2)
whether_to_add_richard_taped_gust_module = 1
#------------------------------------------------------------------------------------------------------
if whether_to_add_richard_taped_gust_module == 1:
    # Custom Xg (pre-generated)
    Xg_data_matrix_taped = np.load('Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\Xg_data_matrix_taped.npy')
    t_taped_for_Xg = np.load('Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\t_taped_for_Xg.npy')
    # Xg_data_matrix_taped_row_index = 0
    # Create an interpolation function for each column in Xg_data_matrix
    interpolated_Xg_taped = interp1d(t_taped_for_Xg, Xg_data_matrix_taped, axis=0, kind='linear', fill_value='extrapolate')
    data_recording_current_duration = 0 # Borrow a variable from "Richard Flight Data Recorder Module"
##############################################################################################################################



##############################################################################################################################
# Richard Taped Stall Buffet Module (place 2 of 2)
whether_to_add_richard_taped_stall_buffet_module = 1
#------------------------------------------------------------------------------------------------------
if whether_to_add_richard_taped_stall_buffet_module == 1:
    # Custom Ay_Az_noises (pre-generated)
    Ay_Az_noises_data_matrix_taped = np.load('Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\Ay_Az_noises_data_matrix_taped.npy')
    t_taped_for_Ay_Az_noises = np.load('Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\t_taped_for_Ay_Az_noises.npy')
    # Ay_Az_noises_data_matrix_taped_row_index = 0
    # Create an interpolation function for each column in Xg_data_matrix
    interpolated_Ay_Az_noises_taped = interp1d(t_taped_for_Ay_Az_noises, Ay_Az_noises_data_matrix_taped, axis=0, kind='linear', fill_value='extrapolate')
    data_recording_current_duration = 0 # Borrow a variable from "Richard Flight Data Recorder Module"
##############################################################################################################################



##############################################################################################################################
# Richard Real Time Gust Module (place 2 of 2)
whether_to_add_richard_real_time_gust_module = 0
#------------------------------------------------------------------------------------------------------
if whether_to_add_richard_real_time_gust_module == 1:
    altitude_points = np.array([-1.00000000e+02,  0.00000000e+00,  1.00000000e+00,  5.00000000e+01,
                                1.00000000e+02,  1.50000000e+02,  2.00000000e+02,  2.50000000e+02,
                                3.04800000e+02,  6.09600000e+02,  1.17199466e+03,  2.34398932e+03,
                                4.55775701e+03,  1.05805073e+04,  1.38034927e+04,  1.68311455e+04,
                                1.97936876e+04,  4.57200000e+04])
    stan_dev_u_points = np.array([0, 0, 1.7300465747703184, 1.3874462887500087, 
                                    1.201581198127368, 1.081200335621849, 0.9946432929792891, 0.9283228269822059, 
                                    0.87072649, 1.62608392, 1.8188234167332495, 1.7304083895309388, 
                                    1.376748280721696, 0.8652041947654694, 0.7199509357902445, 0.4547058541833124, 
                                    0.0, 0.0])
    stan_dev_v_points = np.array([0, 0, 1.7300320302943264, 1.387456241528783, 
                                    1.2015926104134622, 1.0812116077686194, 0.9946540750850655, 0.9283330606659893, 
                                    0.87073615, 1.62611592, 1.8188592053593668, 1.730442438432175, 
                                    1.3767753707234096, 0.8652212192160875, 0.7199651021214158, 0.4547148013398417, 
                                    0.0, 0.0])
    stan_dev_w_points = np.array([0, 0, 0.8705543858535112, 0.8707274776263864, 
                                0.870727643480018, 0.87072662117927, 0.8707253018381317, 0.8707238636791084, 
                                0.87072222, 1.62611592, 1.8188592053593668, 1.730442438432175, 
                                1.3767753707234096, 0.8652212192160875, 0.7199651021214158, 0.4547148013398417, 
                                0.0, 0.0])
    stan_dev_p_points = np.array([0, 0, 0.20600934850341626, 0.05591954047609127, 
                                0.0443833686886628, 0.038772443845604534, 0.03522710307313788, 0.03270194565296064, 
                                0.03061131, 0.04744, 0.053063057949548145, 0.050483603743667324, 
                                0.040165786920144074, 0.025241801871833662, 0.02100412710502947, 0.013265764487387036, 
                                0.0, 0.0])
    stan_dev_q_points = np.array([0, 0, 0.08235526487557372, 0.04152173303832524, 
                                0.03111987019699544, 0.025943871161635775, 0.022709870650808452, 0.020445225731463526, 
                                0.01860422, 0.02651186, 0.029654304407694003, 0.028212775721208873, 
                                0.022446660975268377, 0.014106387860604436, 0.011738162161378874, 0.007413576101923501, 
                                0.0, 0.0])
    stan_dev_r_points = np.array([0, 0, 0.17632270766234787, 0.041899603356915, 
                                0.03201676805602458, 0.027602966162907725, 0.024939226085019205, 0.023095573241348973, 
                                0.02159947, 0.03071001, 0.03435006670116053, 0.03268027179207634, 
                                0.02600109215573957, 0.01634013589603817, 0.013596901402542707, 0.008587516675290133, 
                                0.0, 0.0])

    # Create interpolation functions, these should be global variables in the simulation
    interpolate_stan_dev_u = interp1d(altitude_points, stan_dev_u_points, kind='linear', fill_value="extrapolate")
    interpolate_stan_dev_v = interp1d(altitude_points, stan_dev_v_points, kind='linear', fill_value="extrapolate")
    interpolate_stan_dev_w = interp1d(altitude_points, stan_dev_w_points, kind='linear', fill_value="extrapolate")
    interpolate_stan_dev_p = interp1d(altitude_points, stan_dev_p_points, kind='linear', fill_value="extrapolate")
    interpolate_stan_dev_q = interp1d(altitude_points, stan_dev_q_points, kind='linear', fill_value="extrapolate")
    interpolate_stan_dev_r = interp1d(altitude_points, stan_dev_r_points, kind='linear', fill_value="extrapolate")

##############################################################################################################################


# ##############################################################################################################################
# # Richard Flap Experimentation (place 2 of 2) (should be useless now)
# global_U = np.ones(5)
# global_flap_is_going_down = False
# global_flap_is_going_up = False
# ##############################################################################################################################

# ===================================================================
# Richard Additional Code Block 4 (End)
# ===================================================================

# Force the compilation of the code 
# Functions in this class need to be called continously otherwise they will go cold
# and require recompilation
test = flightModel()
dt = 0.01 # Ric changed it from 0.1 to 0.01
X = np.ones(12)
Xg = np.ones(12)
U = np.ones(5)
test.euler_integrate(dt, X, Xg, U)
test.midpoint_integration(dt, X, Xg, U)

X = np.ones((12, 1))
Xg = np.ones((12, 1))
U = np.ones((5, 1))
test.euler_integrate(dt, X, Xg, U)