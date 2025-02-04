#===================================================
# This version of M10 is for calculating CL, Cm, CD
#===================================================

import numpy as np
from M17_calculate_Cm_increment_from_tail_correction import M17_calculate_Cm_increment_from_tail_correction

# Define the nested function (once incorporated into flightModel.py) M10_calculate_CL_Cm_CD
def M10_calculate_CL_Cm_CD(Vt, FlightData_Geometric_c, alpha, alpha_dot, alpha_crit, X_flow_sep):
    """
    Calculate CL, Cm, CD based on the given parameters.

    Parameters:
    - alpha: Angle of attack (float or int)
    - da_a: (alpha_dot * c_bar) / (2 * V0) (float or int)
    - q_a: (q * c_bar) / (2 * V0) (float or int)
    - U: Control input array containing [U1, U2, U3, U4, U5]
    - CLo: lift coefficent (CL) at zero alpha (o) (float)
    - CLa: Lift coefficient (CL) wrt alpha (a) (float)
    - CLad: Lift coefficient (CL) wrt alpha dot (ad) (float)
    - CLq: Lift coefficient (CL) wrt pitch rate (q) (float)
    - CLde: Lift coefficent (CL) wrt elevator deflection (de) (float)
    - CLdf: Lift coefficent (CL) wrt flap deflection (df) (float)
    - Cmo: pitching moment coefficent (Cm) at zero alpha (o) (float)
    - Cma: pitching moment coefficent (Cm) wrt alpha (a) (float)
    - Cmad: pitching moment coefficent (Cm) wrt alpha dot (ad) (float)
    - Cmq: pitching moment coefficent (Cm) wrt pitch rate (q) (float)
    - Cmde: pitching moment coefficent (Cm) wrt elevator deflection (de) (float)
    - Cmdf: pitching moment coefficent (Cm) wrt flap deflection (df) (float)
    - CD0: zero lift (o) drag coefficient (float)
    - k: induced drag coefficient (float)

    Returns:
    - Calculated value of CL, Cm, CD (float)
    """
    # Add placeholder values for parameters that should be feed to the function later on

    # For CL
    CLo = 0.36486 #* (0.8254/0.36486) # 0.8254
    CLa = 4.9635
    CLad = -1.987
    CLq = 9.8197 * (8/9.8197) # 8
    CLde = 0.0007976 * 1000 # 0.7976
    CLdf = -1.51 #-1.75 #0.0229 * 10 * (0) # 0

    # For Cm
    Cmo = 0.07575
    Cma = -2.6001
    Cmad = -2.210
    Cmq = -23.6813
    Cmde = -0.036 * 100 # -3.6
    Cmdf = -0.002

    # For CD
    CDo = 0.0061 * 8 * (10/9.81) # 0.049745 
    k = 0.0483 * (0.09/0.0483)* (10/9.81) # 0.091743

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    whether_investigate_effect_of_parameters = 0
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if whether_investigate_effect_of_parameters == 1:
        # CLa = 6
        k = 0.05
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # add an empty U vector as a placeholder, feed U to the function later on
    U = [0, 0, 0, 0, 0]

    # assume alpha_dot = q at this stage (for completeness, check out the document named "alpha dot q relation" on iPad)
    q = alpha_dot # assume roughly level flight

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    whether_insert_delft_paper_numbers = 0
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if whether_insert_delft_paper_numbers == 1:
        # a_1 = 40
        # alpha_star = 10 * (np.pi / 180)
        # tau_1 = 0.3 # assume c_bar / Vt is constant for now
        # tau_2 = 0.4 # assume c_bar / Vt is constant for now
        CLo = 0.1
        CLa = 1.6 * np.pi # 5.0265
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    da_a = (alpha_dot * FlightData_Geometric_c) / (2 * Vt)
    q_a = (q * FlightData_Geometric_c) / (2 * Vt)


    #--------------------------------------------------------------------
    whether_to_add_chris_high_aoa_correction_chris_version = 0 # Chris Version (Latest Version) Necessary Toggle Location: 1 of 2
    #--------------------------------------------------------------------
    if whether_to_add_chris_high_aoa_correction_chris_version == 1:
        alpha_crit = -10 * (np.pi/180) # Please enter in deg
    #--------------------------------------------------------------------
    

    # next incorporate the (prestall/poststall) model switch criteria
    if alpha < alpha_crit:
        # For CL
        CL = (CLo + 
              CLa * alpha + 
              CLad * da_a + 
              CLq * q_a + 
              CLde * U[1] + 
              CLdf * U[4])
        
        # For Cm
        Cm = (Cmo + 
              Cma * alpha + 
              Cmad * da_a + 
              Cmq * q_a + 
              Cmde * U[1] + 
              Cmdf * U[4])        
        
        # For CD
        CD = (CDo + 
              k * CL**2)
        
        # indicate that the poststall model is off (0)
        poststall_on = 0

    else:

        # For CL, Cm
        # conditions at alpha_crit, assume stable flight
        da_a_cr = 0
        q_a_cr = 0
        U_1_trim_cr = 0 # assume 0 for now, modify later on
        U_4_trim_cr = 0 # assume 0 for now, modify later on

        # For CD
        CDo_cr = 0.0061 * 8 * (10/9.81) # 0.049745 
        k_cr = 0.0483 * (0.09/0.0483)* (10/9.81) # 0.091743

        # For CL
        CL_cr = (CLo + 
                 CLa * alpha_crit + 
                 CLad * da_a_cr + 
                 CLq * q_a_cr + 
                 CLde * U_1_trim_cr + 
                 CLdf * U_4_trim_cr)

        # For Cm
        Cm_cr = (Cmo + 
                 Cma * alpha_crit + 
                 Cmad * da_a_cr + 
                 Cmq * q_a_cr + 
                 Cmde * U_1_trim_cr + 
                 Cmdf * U_4_trim_cr)
        
        # For CD
        CD_cr = (CDo_cr + 
                 k_cr * CL_cr**2)      
        
        # For CL
        CLawb_cr = CLa
        CLad_cr = CLad
        CLq_cr = CLq

        # For Cm
        Cmawb_cr = Cma
        Cmad_cr = Cmad
        Cmq_cr = Cmq       

        # For Cm
        # Page 18/25 in Wang Paper
        # 1st, 2nd, 3rd order pitching moment coefficient correction terms
        if X_flow_sep >= 0.6:

            CmX1 = 0.34
            CmX2 = -1.5
            CmX3 = 2.2
        
        elif X_flow_sep < 0.6 and X_flow_sep >= 0.43:

            CmX1 = 1.32 - 1.64 * X_flow_sep
            CmX2 = - 9.84 + 13.91 * X_flow_sep
            CmX3 = 9.21 - 11.69 * X_flow_sep

        elif X_flow_sep < 0.43 and X_flow_sep >= 0.16:

            CmX1 = - 4.02 * (X_flow_sep ** 2) + 3.1 * X_flow_sep + 0.03
            CmX2 = 140.95 * (X_flow_sep ** 3) - 144.23 * (X_flow_sep ** 2) + 39.56 * X_flow_sep -5.36
            CmX3 = - 208.32 * (X_flow_sep ** 3) + 225.26 * (X_flow_sep ** 2) - 65.38 * X_flow_sep + 7.21

        else:

            CmX1 = 0.42
            CmX2 = -2.2
            CmX3 = 1.66

        # For CD
        # Page 17/25 in Wang Paper
        # drag coefficient correction term
        CDX = 0.18 * (X_flow_sep**2) - 0.16 * (X_flow_sep) + 0.213

        # assume linear mapping of stick position to control surface deflection angle (as in the default flight simulation physics logic)

        # For CL
        CL = (CL_cr + 
              CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * (alpha - alpha_crit) + 
              CLad_cr * da_a + 
              CLq_cr * q_a + 
              CLde * (U[1] - U_1_trim_cr) + 
              CLdf * (U[4] - U_4_trim_cr))
        
        # For Cm
        Cm = (Cm_cr + 
              Cmawb_cr * (alpha - alpha_crit) + 
              Cmad_cr * da_a + 
              Cmq_cr * q_a + 
              Cmde * (U[1] - U_1_trim_cr) + 
              Cmdf * (U[4] - U_4_trim_cr) + 
              CmX1 * (1 - X_flow_sep) +
              CmX2 * ((1 - X_flow_sep) ** 2) +
              CmX3 * ((1 - X_flow_sep) ** 3))

        # For CD
        CD = (CD_cr +
              k * (CL - CL_cr)**2 +
              CDX * (1 - X_flow_sep)) 
        
        # indicate that the poststall model is on (1)
        poststall_on = 1


        #--------------------------------------------------------------------
        # Richard's heavily modified version of:
        # Eq.(49) model 2 in "High Angle of Attack Parameter Estimation of 
        # Cascaded Fins Using Neural Network" by Sinha et al.
        whether_to_add_richard_high_aoa_correction = 0
        #--------------------------------------------------------------------
        if whether_to_add_richard_high_aoa_correction == 1:
            # For CL
            CL = (CL_cr + 
                  CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * ((alpha - alpha_crit)/(alpha if alpha != 0 else 0.000000001)) * np.sin(alpha) * np.cos(alpha) * (np.cos(alpha) + (2/np.pi) * np.abs(np.sin(alpha))) + 
                  CLad_cr * da_a + 
                  CLq_cr * q_a + 
                  CLde * (U[1] - U_1_trim_cr) * (1 / (1 + np.exp(20 * (X_flow_sep - 0.5)))) + # Adding effect of X_flow_sep on elevator control effectiveness
                  CLdf * (U[4] - U_4_trim_cr)) * (1 / (1 + np.exp(20 * (alpha - 0.5))))
        #--------------------------------------------------------------------


        #--------------------------------------------------------------------
        # Eq.(48) model 1 in "High Angle of Attack Parameter Estimation of 
        # Cascaded Fins Using Neural Network" by Sinha et al.
        whether_to_add_sinha_model_1_high_aoa_correction = 0
        #--------------------------------------------------------------------
        if whether_to_add_sinha_model_1_high_aoa_correction == 1:
            # # For CL
            # CL = (CLo + 
            #       CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * np.sin(alpha) * np.cos(alpha) + # Need careful scrutiny
            #       CLad_cr * da_a + 
            #       CLq_cr * q_a + 
            #       CLde * U[1] + 
            #       CLdf * U[4])
            
            # For CL
            CL = (CL_cr + 
                  CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * ((alpha - alpha_crit)/(alpha if alpha != 0 else 0.000000001)) * np.sin(alpha) * np.cos(alpha) + 
                  CLad_cr * da_a + 
                  CLq_cr * q_a + 
                  CLde * (U[1] - U_1_trim_cr) + 
                  CLdf * (U[4] - U_4_trim_cr))
        #--------------------------------------------------------------------


        #--------------------------------------------------------------------
        # Eq.(49) model 2 in "High Angle of Attack Parameter Estimation of 
        # Cascaded Fins Using Neural Network" by Sinha et al.
        whether_to_add_sinha_model_2_high_aoa_correction = 0
        #--------------------------------------------------------------------
        if whether_to_add_sinha_model_2_high_aoa_correction == 1:
            # # For CL
            # CL = (CLo + 
            #       CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * np.sin(alpha) * np.cos(alpha) * (np.cos(alpha) + (2/np.pi) * np.abs(np.sin(alpha))) + # Need careful scrutiny
            #       CLad_cr * da_a + 
            #       CLq_cr * q_a + 
            #       CLde * U[1] + 
            #       CLdf * U[4])
            
            # # For CL
            # CL = (CL_cr + 
            #       CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * ((alpha - alpha_crit)/(alpha)) * np.sin(alpha) * np.cos(alpha) * (np.cos(alpha) + (2/np.pi) * np.abs(np.sin(alpha))) + 
            #       CLad_cr * da_a + 
            #       CLq_cr * q_a + 
            #       CLde * (U[1] - U_1_trim_cr) + 
            #       CLdf * (U[4] - U_4_trim_cr))
            
            # For CL
            CL = (CL_cr + 
                  CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * ((alpha - alpha_crit)/(alpha if alpha != 0 else 0.000000001)) * np.sin(alpha) * np.cos(alpha) * (np.cos(alpha) + (2/np.pi) * np.abs(np.sin(alpha))) + 
                  CLad_cr * da_a + 
                  CLq_cr * q_a + 
                  CLde * (U[1] - U_1_trim_cr) + 
                  CLdf * (U[4] - U_4_trim_cr)) * (1 / (1 + np.exp(20 * (alpha - 0.5))))
        #--------------------------------------------------------------------


        #--------------------------------------------------------------------
        whether_to_add_chris_high_aoa_correction_richard_version = 0 # Richard Version (Early Version)
        #--------------------------------------------------------------------
        if whether_to_add_chris_high_aoa_correction_richard_version == 1:
            #--------------------------------------------------------------------
            # choose the angle to switch to high aoa aerodynamic model
            switch_angle_for_high_aoa = 25 * (np.pi/180)
            #--------------------------------------------------------------------
            # Add this high AoA correction
            if alpha >= switch_angle_for_high_aoa:
                # For CL
                CL_at_switch = (CL_cr + 
                                CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * (switch_angle_for_high_aoa - alpha_crit) + 
                                CLad_cr * da_a + 
                                CLq_cr * q_a + 
                                CLde * (U[1] - U_1_trim_cr) + 
                                CLdf * (U[4] - U_4_trim_cr))
                
                alpha_in_deg = alpha * (180/np.pi)
                switch_angle_for_high_aoa_in_deg = switch_angle_for_high_aoa * (180/np.pi)
    
                CL_unscaled_after_switch = (-1.19711067e-13 * switch_angle_for_high_aoa_in_deg**8 
                                           + 5.20082162e-11 * switch_angle_for_high_aoa_in_deg**7 
                                           - 9.56411265e-09 * switch_angle_for_high_aoa_in_deg**6 
                                           + 9.67495221e-07 * switch_angle_for_high_aoa_in_deg**5 
                                           - 5.84818549e-05 * switch_angle_for_high_aoa_in_deg**4 
                                           + 2.14165289e-03 * switch_angle_for_high_aoa_in_deg**3 
                                           - 4.57673865e-02 * switch_angle_for_high_aoa_in_deg**2 
                                           + 4.99191880e-01 * switch_angle_for_high_aoa_in_deg 
                                           - 9.02802242e-01)
                
                CL_scale_factor_after_switch = CL_at_switch / CL_unscaled_after_switch
                
                CL = (-1.19711067e-13 * alpha_in_deg**8 
                     + 5.20082162e-11 * alpha_in_deg**7 
                     - 9.56411265e-09 * alpha_in_deg**6 
                     + 9.67495221e-07 * alpha_in_deg**5 
                     - 5.84818549e-05 * alpha_in_deg**4 
                     + 2.14165289e-03 * alpha_in_deg**3 
                     - 4.57673865e-02 * alpha_in_deg**2 
                     + 4.99191880e-01 * alpha_in_deg 
                     - 9.02802242e-01) * CL_scale_factor_after_switch
                
                CD = (-1.10924209e-12 * alpha_in_deg**7 
                      + 4.08803237e-10 * alpha_in_deg**6 
                      - 6.17031482e-08 * alpha_in_deg**5 
                      + 4.93635653e-06 * alpha_in_deg**4 
                      - 2.32257424e-04 * alpha_in_deg**3 
                      + 6.77591188e-03 * alpha_in_deg**2 
                      - 8.78088130e-02 * alpha_in_deg 
                      + 5.59383006e-01)
        #--------------------------------------------------------------------


        #--------------------------------------------------------------------
        whether_to_add_chris_high_aoa_correction_chris_version = 1 # Chris Version (Latest Version) Necessary Toggle Location: 2 of 2
        #--------------------------------------------------------------------
        if whether_to_add_chris_high_aoa_correction_chris_version == 1:
            #--------------------------------------------------------------------
            # choose the angle to switch to high aoa aerodynamic model
            switch_angle_for_high_aoa = 25 * (np.pi/180)
            #--------------------------------------------------------------------
            # Add this high AoA correction
            if alpha >= switch_angle_for_high_aoa:
                # For CL
                CL_at_switch = (CL_cr + 
                                CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * (switch_angle_for_high_aoa - alpha_crit) + 
                                CLad_cr * da_a + 
                                CLq_cr * q_a + 
                                CLde * (U[1] - U_1_trim_cr) + 
                                CLdf * (U[4] - U_4_trim_cr))
                
                alpha_in_deg = alpha * (180/np.pi)
                switch_angle_for_high_aoa_in_deg = switch_angle_for_high_aoa * (180/np.pi)

                if switch_angle_for_high_aoa_in_deg >= 0 and switch_angle_for_high_aoa_in_deg < 10:
    
                    CL_unscaled_after_switch = (0.084 * switch_angle_for_high_aoa_in_deg 
                                              + 0.43)
                    
                elif switch_angle_for_high_aoa_in_deg >= 10 and switch_angle_for_high_aoa_in_deg < 16:

                    CL_unscaled_after_switch = (-0.0117 * switch_angle_for_high_aoa_in_deg**2 
                                               + 0.315 * switch_angle_for_high_aoa_in_deg
                                               - 0.7133)
                    
                if switch_angle_for_high_aoa_in_deg >= 16 and switch_angle_for_high_aoa_in_deg < 21:
    
                    CL_unscaled_after_switch = (-0.052 * switch_angle_for_high_aoa_in_deg 
                                               + 2.17)
                    
                if switch_angle_for_high_aoa_in_deg >= 21 and switch_angle_for_high_aoa_in_deg < 30:
    
                    CL_unscaled_after_switch = (-0.017 * switch_angle_for_high_aoa_in_deg 
                                               + 1.44)                 

                elif switch_angle_for_high_aoa_in_deg >= 30 and switch_angle_for_high_aoa_in_deg <= 90:

                    CL_unscaled_after_switch = (-0.0000801 * switch_angle_for_high_aoa_in_deg**2 
                                               - 0.00599 * switch_angle_for_high_aoa_in_deg
                                               + 1.1876)
                    
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                whether_enable_CL_scale_factor = 1
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                if whether_enable_CL_scale_factor == 1:
                    CL_scale_factor_after_switch = CL_at_switch / CL_unscaled_after_switch
                else:
                    CL_scale_factor_after_switch = 1 # effectively no scale
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                if alpha_in_deg >= 0 and alpha_in_deg < 10:
    
                    CL = (0.084 * alpha_in_deg 
                                              + 0.43)
                    
                elif alpha_in_deg >= 10 and alpha_in_deg < 16:

                    CL = (-0.0117 * alpha_in_deg**2 
                                               + 0.315 * alpha_in_deg
                                               - 0.7133)
                    
                if alpha_in_deg >= 16 and alpha_in_deg < 21:
    
                    CL = (-0.052 * alpha_in_deg 
                                               + 2.17)
                    
                if alpha_in_deg >= 21 and alpha_in_deg < 30:
    
                    CL = (-0.017 * alpha_in_deg 
                                               + 1.44)                 

                elif alpha_in_deg >= 30 and alpha_in_deg <= 90:

                    CL = (-0.0000801 * alpha_in_deg**2 
                                               - 0.00599 * alpha_in_deg
                                               + 1.1876)
                                    
                CL = CL * CL_scale_factor_after_switch


                # For CD
                CD_at_switch = (CD_cr +
                      k * (CL_at_switch - CL_cr)**2 +
                      CDX * (1 - X_flow_sep)) 
                
                # Drag model up to 180 deg AoA
                # CD90 = 2 # Hoerner
                CD90 = 4.9322919 # J300: CD90 = (np.pi * AR) / (1 + np.sqrt((AR / 2)**2 + 1)), AR = b**2 / S = 8.181818...

                # CD_unscaled_after_switch = CD90 * (np.sin(switch_angle_for_high_aoa)**2)
                CD_unscaled_after_switch = CDo * (np.cos(switch_angle_for_high_aoa)**2) + CD90 * (np.sin(switch_angle_for_high_aoa)**2)
                
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                whether_enable_CD_scale_factor = 1
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                if whether_enable_CD_scale_factor == 1:
                    CD_scale_factor_after_switch = CD_at_switch / CD_unscaled_after_switch
                    #----------------------------------------------------------------------
                    if alpha < (np.pi/3):
                        # #..................................................................
                        # # The following line applies a linear scale to the scale factor from switch_angle_for_high_aoa (scale = scale_factor) to np.pi/3 (scale = 1).
                        # CD_scale_factor_after_switch = ((1 - CD_scale_factor_after_switch)/((np.pi/3) - switch_angle_for_high_aoa)) * (alpha - (np.pi/3)) + 1
                        # #..................................................................
                        #..................................................................
                        # The following line applies a quadratic scale to the scale factor from switch_angle_for_high_aoa (scale = scale_factor) to np.pi/3 (scale = 1).
                        alternative_coordinate_system_x = ((1 - 0)/(switch_angle_for_high_aoa - (np.pi/3))) * (alpha - switch_angle_for_high_aoa) + 1
                        alternative_coordinate_system_y = alternative_coordinate_system_x ** 2 # quadraic, to have smoother CD
                        CD_scale_factor_after_switch = ((CD_scale_factor_after_switch - 1)/(1 - 0)) * (alternative_coordinate_system_y - 1) + CD_scale_factor_after_switch
                        #..................................................................
                    else:
                        CD_scale_factor_after_switch = 1 # no scale from np.pi/4 to np.pi/2
                    #----------------------------------------------------------------------
                else:
                    CD_scale_factor_after_switch = 1 # effectively no scale
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # consider only scaling the part between 25 deg and 60 deg using some kind of method (keeping CD90 the same) (This is already done.)
                
                # Drag model up to 180 deg AoA
                # CD90 = 2 # Hoerner
                CD90 = 4.9322919 # J300: CD90 = (np.pi * AR) / (1 + np.sqrt((AR / 2)**2 + 1)), AR = b**2 / S = 8.181818...

                # CD = CD90 * (np.sin(alpha)**2)
                CD = CDo * (np.cos(alpha)**2) + CD90 * (np.sin(alpha)**2)

                CD = CD * CD_scale_factor_after_switch    

        #--------------------------------------------------------------------            

    # # For Cm
    #....................................................................
    whether_to_add_richard_Cm_increment_from_tail_correction = 1 # Richard Cm increment from tail correction, derived based on Philips (Mechanics of Flight 1st Ed)
    #....................................................................
    if whether_to_add_richard_Cm_increment_from_tail_correction == 1:

        alpha_horizontal_tail, Cm_increment_from_tail_correction = M17_calculate_Cm_increment_from_tail_correction(alpha, CL)

        Cm = Cm + Cm_increment_from_tail_correction
    #....................................................................


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    whether_insert_delft_paper_static_CL = 0
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if whether_insert_delft_paper_static_CL == 1:
        # a_1 = 40
        # alpha_star = 10 * (np.pi / 180)
        # tau_1 = 0.3 # assume c_bar / Vt is constant for now
        # tau_2 = 0.4 # assume c_bar / Vt is constant for now
        # CLo = 0.1
        # CLa = 1.6 * np.pi # 5.0265
        CL = (CLo +
              CLa * (((1 + np.sqrt(X_flow_sep))/2)**2) * alpha) # add bracket to avoid the need of a breakpoint in Python
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    return CL, poststall_on, Cm, CD