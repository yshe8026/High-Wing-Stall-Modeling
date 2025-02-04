import numpy as np

# global_time = 0
# global_alpha = 0
# global_dt = 0.01
# X_flow_sep = 0
# dt = 0.01

def M06_calculate_alpha_dot(global_dt, global_alpha, alpha):
    """
    Calculate alpha_dot based on the given parameters.

    Parameters:
    - global_dt: Time step (float or int)
    - alpha: Current alpha angle (float or int)
    - global_alpha: Previous global alpha angle (float or int)

    Returns:
    - Calculated value of alpha_dot (float)
    """

    # # Declare the variable as global to ensure it is recognized at the global scope (time step size for the LAST timestep)
    # global global_dt
    try:
        # Try to see if global variable 'global_dt' is defined (Modified)
        global_dt
    except NameError:
        # If 'global_dt' has not been defined, catch the NameError exception
        # Assign the value 0.01 to 'global_dt' since it was not previously defined
        global_dt = 0.01

    dt_time_step = global_dt
    alpha_dot = (alpha - global_alpha) / dt_time_step

    # return the value of "alpha_dot" at the very end of this function
    return alpha_dot