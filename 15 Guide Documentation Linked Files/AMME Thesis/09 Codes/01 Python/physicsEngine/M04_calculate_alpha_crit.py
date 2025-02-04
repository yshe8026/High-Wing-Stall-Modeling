import numpy as np

# global_time = 0
# global_alpha = 0
# global_dt = 0.01
# X_flow_sep = 0
# dt = 0.01

def M04_calculate_alpha_crit(flag_which_alpha_crit_interpretation, global_time, dX_flow_sep_by_dt, alpha_dot, FlightData_Geometric_c, Vt):
    """
    Calculate alpha_crit based on the given parameters.

    Parameters:
    - global_time: Current global time (float or int)
    - dt: Time step (float or int)
    - alpha: Current alpha angle (float or int)
    - global_alpha: Previous global alpha angle (float or int)
    - FlightData_Geometric_c: Geometric parameter (float or int)
    - Vt: Velocity parameter (float or int)

    Returns:
    - Calculated value of alpha_crit (float)
    """

    # #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # # Select which interpretation of alpha_crit to run (See Overleaf) (1, 2, 3)
    # #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # flag_which_alpha_crit_interpretation = 1
    # #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # note: initialization lasts around 2.3 seconds, or 230 iterations, during which Vt and term_5 can be nan
    init_wait_time = 0
    # # # Declare the variable as global to ensure it is recognized at the global scope (time step size for the LAST timestep)
    # # global global_dt
    # try:
    #     # Try to see if global variable 'global_dt' is defined (Modified)
    #     global_dt
    # except NameError:
    #     # If 'global_dt' has not been defined, catch the NameError exception
    #     # Assign the value 0.01 to 'global_dt' since it was not previously defined
    #     global_dt = 0.01

    # dt_time_step = global_dt
    # alpha_dot = (alpha - global_alpha) / dt_time_step

    X_crit = 0.95

    a_1 = (1.00) * 22.5
    alpha_star = 20 * (np.pi / 180)

    if flag_which_alpha_crit_interpretation == 1:

        tau_2 = (1.00) * 6.66 * (FlightData_Geometric_c / (Vt if global_time > init_wait_time else 52.65)) # does not impact how X behaves if alpha is constant
        # interpretation 1 formula for alpha_crit
        alpha_crit = 0.8 * ((1/a_1) * np.arctanh(1 - 2 * X_crit) + alpha_star + tau_2 * alpha_dot)

    elif flag_which_alpha_crit_interpretation == 2:

        tau_1 = (1.00) * 11.93 * (FlightData_Geometric_c / (Vt if global_time > init_wait_time else 52.65)) # determines how fast X converges if alpha is constant
        tau_2 = (1.00) * 6.66 * (FlightData_Geometric_c / (Vt if global_time > init_wait_time else 52.65)) # does not impact how X behaves if alpha is constant
        # interpretation 2 formula for alpha_crit
        alpha_crit = 0.8 * ((1/a_1) * np.arctanh(1 - 2 * X_crit - 2 * tau_1 * dX_flow_sep_by_dt) + alpha_star + tau_2 * alpha_dot)

    elif flag_which_alpha_crit_interpretation == 3:

        # interpretation 3 formula for alpha_crit
        alpha_crit = 0.8 * ((1/a_1) * np.arctanh(1 - 2 * X_crit) + alpha_star)

    # return the value of "alpha_crit" at the very end of this function
    return alpha_crit