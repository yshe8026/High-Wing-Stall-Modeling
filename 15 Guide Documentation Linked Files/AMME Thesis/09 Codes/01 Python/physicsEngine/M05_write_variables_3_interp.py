#==============================================
# This version of M05 is currently in use
#==============================================

import numpy as np

# global_time = 0
# global_alpha = 0
# global_dt = 0.01
# X_flow_sep = 0
# dt = 0.01

def M05_write_variables_3_interp(global_time, global_dt, dX_flow_sep_by_dt, X_flow_sep, alpha, alpha_crit_1, alpha_crit_2, alpha_crit_3, alpha_dot, Vt, CL, poststall_on, Cm, CD):
    """
    Write given parameters into file.

    Parameters:
    - global_time: Current global time (float or int)
    - global_dt: Time step (float or int)
    - dX_flow_sep_by_dt: rate of change of X_flow_sep (float or int)
    - X_flow_sep: flow separation paramter (float or int)
    - alpha: Current alpha angle (float or int)
    - alpha_crit: Critical AoA (float or int)
    - alpha_dot: rate of change of AoA (float or int)
    - Vt: Velocity parameter (float or int)

    Returns:
    - Nothing
    """
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Whether to write angles in deg (1: Yes) (0: no, write angles in rad)
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    flag_write_angles_in_deg = 1
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # write combined    
    # with open('richard_checked_terms\\combined_history.txt', 'a') as file:
    #     file.write(f"{alpha_dot}\t{X_flow_sep}\n")
    file_path = 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\stage2_combined_history.txt'
    # column width in terms of character string length
    cwidth = 25

    # Function to format the columns
    def format_line(global_time, global_dt, dX_flow_sep_by_dt, X_flow_sep, alpha, alpha_crit_1, alpha_crit_2, alpha_crit_3, alpha_dot, Vt, CL, poststall_on, Cm, CD, cwidth):
        # Format global_time to 5 significant figures
        formatted_global_time = f"{global_time:.5g}"
        return f"{formatted_global_time:<{cwidth}}{global_dt:<{cwidth}}{dX_flow_sep_by_dt:<{cwidth}}{X_flow_sep:<{cwidth}}{alpha:<{cwidth}}{alpha_crit_1:<{cwidth}}{alpha_crit_2:<{cwidth}}{alpha_crit_3:<{cwidth}}{alpha_dot:<{cwidth}}{Vt:<{cwidth}}{CL:<{cwidth}}{poststall_on:<{cwidth}}{Cm:<{cwidth}}{CD:<{cwidth}}\n"

    try:
        # Check if the file is empty
        with open(file_path, 'r') as file:
            is_empty = file.read().strip() == ""
    except FileNotFoundError:
        # If the file doesn't exist, it will be considered empty
        is_empty = True

    #+++++++++++++++++++++++++++++++++++++++++++++
    # Write angles in rad
    #+++++++++++++++++++++++++++++++++++++++++++++
    if flag_write_angles_in_deg == 1:
        alpha_in_deg = alpha * (180/np.pi)

        alpha_crit_1_in_deg = alpha_crit_1 * (180/np.pi)
        alpha_crit_2_in_deg = alpha_crit_2 * (180/np.pi)
        alpha_crit_3_in_deg = alpha_crit_3 * (180/np.pi)

        alpha_dot_in_deg_per_sec = alpha_dot * (180/np.pi)
        # Write the column titles if the file is empty
        with open(file_path, 'a') as file:
            if is_empty:
                file.write(f"{'global_time[s]':<{cwidth}}{'global_dt[s]':<{cwidth}}{'dX_flow_sep_by_dt[-/s]':<{cwidth}}{'X_flow_sep[-]':<{cwidth}}{'alpha[deg]':<{cwidth}}{'alpha_crit_1[deg]':<{cwidth}}{'alpha_crit_2[deg]':<{cwidth}}{'alpha_crit_3[deg]':<{cwidth}}{'alpha_dot[deg/s]':<{cwidth}}{'Vt[m/s]':<{cwidth}}{'CL[-]':<{cwidth}}{'poststall_on[-]':<{cwidth}}{'Cm[-]':<{cwidth}}{'CD[-]':<{cwidth}}\n")
            file.write(format_line(global_time, global_dt, dX_flow_sep_by_dt, X_flow_sep, alpha_in_deg, alpha_crit_1_in_deg, alpha_crit_2_in_deg, alpha_crit_3_in_deg, alpha_dot_in_deg_per_sec, Vt, CL, poststall_on, Cm, CD, cwidth))

    elif flag_write_angles_in_deg == 0:
        # Write the column titles if the file is empty
        with open(file_path, 'a') as file:
            if is_empty:
                file.write(f"{'global_time[s]':<{cwidth}}{'global_dt[s]':<{cwidth}}{'dX_flow_sep_by_dt[-/s]':<{cwidth}}{'X_flow_sep[-]':<{cwidth}}{'alpha[rad]':<{cwidth}}{'alpha_crit_1[rad]':<{cwidth}}{'alpha_crit_2[rad]':<{cwidth}}{'alpha_crit_3[rad]':<{cwidth}}{'alpha_dot[rad/s]':<{cwidth}}{'Vt[m/s]':<{cwidth}}{'CL[-]':<{cwidth}}{'poststall_on[-]':<{cwidth}}{'Cm[-]':<{cwidth}}{'CD[-]':<{cwidth}}\n")
            file.write(format_line(global_time, global_dt, dX_flow_sep_by_dt, X_flow_sep, alpha, alpha_crit_1, alpha_crit_2, alpha_crit_3, alpha_dot, Vt, CL, poststall_on, Cm, CD, cwidth))