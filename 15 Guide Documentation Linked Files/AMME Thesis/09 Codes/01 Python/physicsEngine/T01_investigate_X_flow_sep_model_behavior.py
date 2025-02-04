# test_script.py
import numpy as np
import matplotlib.pyplot as plt
from M01_calculate_dX_flow_sep_by_dt import M01_calculate_dX_flow_sep_by_dt
from M02_update_X_flow_sep import M02_update_X_flow_sep
from M03_update_global_time import M03_update_global_time
from M04_calculate_alpha_crit import M04_calculate_alpha_crit
from M06_calculate_alpha_dot import M06_calculate_alpha_dot
from M05_write_variables_3_interp import M05_write_variables_3_interp
from M07_calculate_CL import M07_calculate_CL
from M08_calculate_Cm import M08_calculate_Cm
from M09_calculate_CD import M09_calculate_CD
from M10_calculate_CL_Cm_CD import M10_calculate_CL_Cm_CD
from T02_automation import dt, FlightData_Geometric_c, Vt, alpha_series

# # Test cases (This is now imported from T02_automation for convenience)
# #++++++++++++++++++++++++++++++++++++++++++++
# selected_test_case = 2
# #++++++++++++++++++++++++++++++++++++++++++++

# if selected_test_case == 1:
#     # Parameters
#     time_steps = 100
#     dt = 0.01
#     FlightData_Geometric_c = 0.99
#     Vt = 40.0
#     alpha_constant_value = 10 * (np.pi/180)
#     alpha_series = np.full(time_steps, alpha_constant_value)  # Constant alpha of 0.5

# if selected_test_case == 2:
#     # Parameters
#     time_steps = 100
#     dt = 0.01
#     FlightData_Geometric_c = 0.99
#     Vt = 40.0
#     alpha_series = np.linspace(0, 1, time_steps)  # Linear increase from 0 to 1

# initialization of other variables
global_time = 0
global_alpha = 0
X_flow_sep = 1

#+++++++++++++++++++++++++++++++++++++++++++++
# Change the simulation time step size
#+++++++++++++++++++++++++++++++++++++++++++++
global_dt = 0.1 # default is 0.01
#+++++++++++++++++++++++++++++++++++++++++++++

# # Lists to store results for plotting
# global_time_series = []
# X_flow_sep_series = []

# Simulation loop
for alpha in alpha_series:
    alpha_dot = M06_calculate_alpha_dot(global_dt, global_alpha, alpha)
    # 'M01_calculate_dX_flow_sep_by_dt' updates global_alpha automatically
    dX_flow_sep_by_dt, global_alpha= M01_calculate_dX_flow_sep_by_dt(global_time, global_dt, X_flow_sep, global_alpha, alpha, FlightData_Geometric_c, Vt)
    X_flow_sep = M02_update_X_flow_sep(X_flow_sep, dX_flow_sep_by_dt, dt)

    alpha_crit_1 = M04_calculate_alpha_crit(1, global_time, dX_flow_sep_by_dt, alpha_dot, FlightData_Geometric_c, Vt)
    alpha_crit_2 = M04_calculate_alpha_crit(2, global_time, dX_flow_sep_by_dt, alpha_dot, FlightData_Geometric_c, Vt)
    alpha_crit_3 = M04_calculate_alpha_crit(3, global_time, dX_flow_sep_by_dt, alpha_dot, FlightData_Geometric_c, Vt)

    #+++++++++++++++++++++++++++++++++++++++++++++
    # Choose alpba_crit interpretation version
    #+++++++++++++++++++++++++++++++++++++++++++++
    alpha_crit = alpha_crit_1
    #+++++++++++++++++++++++++++++++++++++++++++++

    # calculate CL, Cm, CD

    # CL, poststall_on = M07_calculate_CL(Vt, FlightData_Geometric_c, alpha, alpha_dot, alpha_crit, X_flow_sep)
    # Cm = M08_calculate_Cm(Vt, FlightData_Geometric_c, alpha, alpha_dot, alpha_crit, X_flow_sep)
    # CD = M09_calculate_CD(Vt, FlightData_Geometric_c, alpha, alpha_dot, alpha_crit, X_flow_sep)

    # Combine three modules M07, M08, M09 into one M10, to ease investigation of parameters.
    CL, poststall_on, Cm, CD = M10_calculate_CL_Cm_CD(Vt, FlightData_Geometric_c, alpha, alpha_dot, alpha_crit, X_flow_sep)
    
    # Write parameters into .txt file stage2_combined_history
    M05_write_variables_3_interp(global_time, global_dt, dX_flow_sep_by_dt, X_flow_sep, alpha, alpha_crit_1, alpha_crit_2, alpha_crit_3, alpha_dot, Vt, CL, poststall_on, Cm, CD)
    # update global_time at the very end of a loop
    global_time = M03_update_global_time(global_time, dt)

    # print(global_time)
    
    # # Store results
    # global_time_series.append(global_time)
    # X_flow_sep_series.append(X_flow_sep)

# # Print the results
# print("Final global_time:", global_time)
# print("Final X_flow_sep:", X_flow_sep)













