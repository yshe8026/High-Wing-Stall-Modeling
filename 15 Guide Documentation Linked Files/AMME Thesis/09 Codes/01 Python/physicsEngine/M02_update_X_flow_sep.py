def M02_update_X_flow_sep(X_flow_sep, dX_flow_sep_by_dt, dt):
    """
    Update X_flow_sep using the formula:
    X_flow_sep += dX_flow_sep_by_dt * dt

    Parameters:
    - X_flow_sep: Current value of X_flow_sep (float or int)
    - dX_flow_sep_by_dt: Rate of change of X_flow_sep with respect to time (float or int)
    - dt: Time step (float or int)

    Returns:
    - Updated value of X_flow_sep (float or int)
    """

    X_flow_sep += dX_flow_sep_by_dt * dt

    return X_flow_sep