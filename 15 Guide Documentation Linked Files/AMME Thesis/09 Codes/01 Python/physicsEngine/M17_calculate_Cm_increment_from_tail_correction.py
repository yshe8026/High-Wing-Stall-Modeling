import numpy as np

'''This function is called by M10 function modules for Cm corrections from tail stall and tail drag'''


def M17_calculate_Cm_increment_from_tail_correction(alpha, CL):
# def M17_calculate_Cm_increment_from_tail_correction(global_time, global_dt, X_flow_sep, global_alpha, alpha, FlightData_Geometric_c, Vt, CL):

    # i_horizontal_tail_setting_angle = np.deg2rad(-2) # rad # J400's tail setting angle (maybe need a more precise number if possible to acquire, but the value is close enough for now.)
    i_horizontal_tail_setting_angle = np.deg2rad(-3) # rad # J400's tail setting angle (maybe need a more precise number if possible to acquire, but the value is close enough for now.)
    AR_wing = 8.181818
    # Philips Mechanics of Flight 1st Ed, page 391 to page 397.
    # The following values are for J400 by reading the charts, Figure 4.5.2 to Figure 4.5.6.
    # kappa_v = 0.88
    # kappa_p = 0.48
    # kappa_s = 1.02
    # kappa_b = 0.887
    # proportional_factor_downwash_angle_to_wing_CL = (kappa_v * kappa_p * kappa_s)/(kappa_b * AR_wing)
    proportional_factor_downwash_angle_to_wing_CL = 0.05936776
    CL_wing = CL * (1.6 / 1.85) # By interpreting the results of R05 and T02, we assume 86.4 % of J400's lift come from its main wing (this becomes just a rough approximation if we consider large AoA.)
    epsilon_d_down_wash_angle = proportional_factor_downwash_angle_to_wing_CL * CL_wing # rad

    alpha_horizontal_tail = alpha - epsilon_d_down_wash_angle + i_horizontal_tail_setting_angle
    # alpha_horizontal_tail = alpha

    # convert alpha of horizontal tail to deg
    alpha_horizontal_tail_in_deg = np.rad2deg(alpha_horizontal_tail)

    # Default value of moment correction
    Cm_increment_from_tail_correction = 0

    # # Define the alpha range for each segment of lookup table, which is in turn calculated using Philips and R42

    # alpha_range_segment1 = alpha[(alpha >= -12) & (alpha <= 12)]
    if (alpha_horizontal_tail_in_deg >= -12) and (alpha_horizontal_tail_in_deg <= 12):

        # Segment 1: Linear fit (-12 to 12 degrees)
        Cm_increment_from_tail_correction = -0.0008343000243785809 * alpha_horizontal_tail_in_deg - 2.514225649552972e-18

    # alpha_range_segment2 = alpha[(alpha > 12) & (alpha <= 15)]
    elif (alpha_horizontal_tail_in_deg > 12) and (alpha_horizontal_tail_in_deg <= 15):

        # Segment 2: 2nd Order Polynomial (12 to 15 degrees)
        Cm_increment_from_tail_correction = (0.031980941153777045 * alpha_horizontal_tail_in_deg**2 
                                           - 0.7710063139336805 * alpha_horizontal_tail_in_deg 
                                           + 4.652377288950034)

    # alpha_range_segment3 = alpha[(alpha > 15) & (alpha <= 20)]
    elif (alpha_horizontal_tail_in_deg > 15) and (alpha_horizontal_tail_in_deg <= 20):

        # Segment 3: Linear fit (15 to 20 degrees)
        Cm_increment_from_tail_correction = 0.19313442481940202 * alpha_horizontal_tail_in_deg - 2.6390583347540058

    # alpha_range_segment4 = alpha[(alpha > 20) & (alpha <= 90)] correct to 77.1 deg
    elif (alpha_horizontal_tail_in_deg > 20) and (alpha_horizontal_tail_in_deg <= 77.1):    

        # Segment 4: 3rd Order Polynomial (20 to 90 degrees)
        Cm_increment_from_tail_correction = (
                                            -5.326176479063682e-06 * alpha_horizontal_tail_in_deg**3
                                            - 0.0005525779113814487 * alpha_horizontal_tail_in_deg**2
                                            + 0.07448956312117563 * alpha_horizontal_tail_in_deg
                                            + 0.0005870900914903243)
    else:
        # Assume otherwise no Cm correction is needed
        Cm_increment_from_tail_correction = 0


    #--------------------------------------------------------------------------------------------------------------------------------
    whether_add_richard_Cm_increment_tune_down_factor = 1
    #--------------------------------------------------------------------------------------------------------------------------------
    if whether_add_richard_Cm_increment_tune_down_factor == 1:

        # User input
        richard_Cm_increment_tune_down_factor = 0.5

        Cm_increment_from_tail_correction = Cm_increment_from_tail_correction * richard_Cm_increment_tune_down_factor
    #--------------------------------------------------------------------------------------------------------------------------------

    return alpha_horizontal_tail, Cm_increment_from_tail_correction