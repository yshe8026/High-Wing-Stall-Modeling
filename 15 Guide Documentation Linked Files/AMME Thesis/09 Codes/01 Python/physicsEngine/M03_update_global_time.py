def M03_update_global_time(global_time, dt):
    """
    Update global_time based on the formula:
    global_time += dt

    Parameters:
    - global_time: Current value of global_time (float or int)
    - dt: Time step (float or int)

    Returns:
    - Updated value of global_time (float or int)
    """

    global_time += dt

    return global_time
