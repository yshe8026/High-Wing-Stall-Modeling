#============================================
# This version of M08 is for calculating Cm
#============================================

import numpy as np

# Define the nested function (once incorporated into flightModel.py) M07_calculate_CL
def M08_calculate_Cm(Vt, FlightData_Geometric_c, alpha, alpha_dot, alpha_crit, X_flow_sep):
    """
    Calculate Cm based on the given parameters.

    Parameters:
    - alpha: Angle of attack (float or int)
    - da_a: (alpha_dot * c_bar) / (2 * V0) (float or int)
    - q_a: (q * c_bar) / (2 * V0) (float or int)
    - U: Control input array containing [U1, U2, U3, U4, U5]
    - Cmo: pitching moment coefficent (Cm) at zero alpha (o) (float)
    - Cma: pitching moment coefficent (Cm) wrt alpha (a) (float)
    - Cmad: pitching moment coefficent (Cm) wrt alpha dot (ad) (float)
    - Cmq: pitching moment coefficent (Cm) wrt pitch rate (q) (float)
    - Cmde: pitching moment coefficent (Cm) wrt elevator deflection (de) (float)
    - Cmdf: pitching moment coefficent (Cm) wrt flap deflection (df) (float)

    Returns:
    - Calculated value of Cm (float)
    """

    # Cm = self.FlightData_Aero_Cmo + 
    #      self.FlightData_Aero_Cma*alpha + 
    #      self.FlightData_Aero_Cmad*da_a + 
    #      self.FlightData_Aero_Cmq*q_a + 
    #      self.FlightData_Aero_Cmde*U[1] + 
    #      self.FlightData_Aero_Cmdf*U[4]


    # #  ,CLo, CLa, CLad, CLq, CLde, CLdf
    # # add placeholder values for parameters that should be feed to the function later on
    # CLo = 0.36486 * (0.8254/0.36486) # 0.8254
    # CLa = 4.9635
    # CLad = -1.987
    # CLq = 9.8197 * (8/9.8197) # 8
    # CLde = 0.0007976 * 1000 # 0.7976
    # CLdf = 0.0229 * 10 * (0) # 0

    # 'Cmo' : 0.07575,                  # pitching moment coefficent (Cm) at zero alpha (o)
    # 'Cma' : -2.6001,                  # pitching moment coefficent (Cm) wrt alpha (a)
    # 'Cmq' : -23.6813,                 # pitching moment coefficent (Cm) wrt pitch rate (p)
    # 'Cmad' : -2.210,                  # pitching moment coefficent (Cm) wrt alpha dot(ad)
    # 'Cmde' : -0.036 * 100,            # pitching moment coefficent (Cm) wrt elevator angle delta(de)
    # 'Cmdf' : -0.002,                  # pitching moment coefficent (Cm) wrt flap deflection (df)

    Cmo = 0.07575
    Cma = -2.6001
    Cmad = -2.210
    Cmq = -23.6813
    Cmde = -0.036 * 100 # -3.6
    Cmdf = -0.002

    # add an empty U vector as a placeholder, feed U to the function later on
    U = [0, 0, 0, 0, 0]

    # assume alpha_dot = q at this stage (for completeness, check out the document named "alpha dot q relation" on iPad)
    q = alpha_dot # assume roughly level flight

    # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # whether_insert_delft_paper_numbers = 0
    # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # if whether_insert_delft_paper_numbers == 1:
    #     # a_1 = 40
    #     # alpha_star = 10 * (np.pi / 180)
    #     # tau_1 = 0.3 # assume c_bar / Vt is constant for now
    #     # tau_2 = 0.4 # assume c_bar / Vt is constant for now
    #     CLo = 0.1
    #     CLa = 1.6 * np.pi # 5.0265
    # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    da_a = (alpha_dot * FlightData_Geometric_c) / (2 * Vt)
    q_a = (q * FlightData_Geometric_c) / (2 * Vt)

    # next incorporate the (prestall/poststall) model switch criteria
    if alpha < alpha_crit:
        # CL = (CLo + 
        #       CLa * alpha + 
        #       CLad * da_a + 
        #       CLq * q_a + 
        #       CLde * U[1] + 
        #       CLdf * U[4])
        
        Cm = (Cmo + 
              Cma * alpha + 
              Cmad * da_a + 
              Cmq * q_a + 
              Cmde * U[1] + 
              Cmdf * U[4])        
        
        # # indicate that the poststall model is off (0)
        # poststall_on = 0

    else:

        # conditions at alpha_crit, assume stable flight
        da_a_cr = 0
        q_a_cr = 0
        U_1_trim_cr = 0 # assume 0 for now, modify later on
        U_4_trim_cr = 0 # assume 0 for now, modify later on

        # CL_cr = (CLo + 
        #          CLa * alpha_crit + 
        #          CLad * da_a_cr + 
        #          CLq * q_a_cr + 
        #          CLde * U_1_trim_cr + 
        #          CLdf * U_4_trim_cr)
        
        Cm_cr = (Cmo + 
                 Cma * alpha_crit + 
                 Cmad * da_a_cr + 
                 Cmq * q_a_cr + 
                 Cmde * U_1_trim_cr + 
                 Cmdf * U_4_trim_cr)
        
        # CLawb_cr = CLa
        # CLad_cr = CLad
        # CLq_cr = CLq

        Cmawb_cr = Cma
        Cmad_cr = Cmad
        Cmq_cr = Cmq

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

        # assume linear mapping of stick position to control surface deflection angle (as in the default flight simulation physics logic)

        # CL = (CL_cr + 
        #       CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * (alpha - alpha_crit) + 
        #       CLad_cr * da_a + 
        #       CLq_cr * q_a + 
        #       CLde * (U[1] - U_1_trim_cr) + 
        #       CLdf * (U[4] - U_4_trim_cr))

        Cm = (Cm_cr + 
              Cmawb_cr * (alpha - alpha_crit) + 
              Cmad_cr * da_a + 
              Cmq_cr * q_a + 
              Cmde * (U[1] - U_1_trim_cr) + 
              Cmdf * (U[4] - U_4_trim_cr) + 
              CmX1 * (1 - X_flow_sep) +
              CmX2 * ((1 - X_flow_sep) ** 2) +
              CmX3 * ((1 - X_flow_sep) ** 3))
        
        # # indicate that the poststall model is on (1)
        # poststall_on = 1

    # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # whether_insert_delft_paper_static_CL = 0
    # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # if whether_insert_delft_paper_static_CL == 1:
    #     # a_1 = 40
    #     # alpha_star = 10 * (np.pi / 180)
    #     # tau_1 = 0.3 # assume c_bar / Vt is constant for now
    #     # tau_2 = 0.4 # assume c_bar / Vt is constant for now
    #     # CLo = 0.1
    #     # CLa = 1.6 * np.pi # 5.0265
    #     CL = (CLo +
    #           CLa * (((1 + np.sqrt(X_flow_sep))/2)**2) * alpha) # add bracket to avoid the need of a breakpoint in Python
    # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    return Cm
