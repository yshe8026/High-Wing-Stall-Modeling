#============================================
# This version of M07 is for calculating CL 
#============================================

import numpy as np

# Define the nested function (once incorporated into flightModel.py) M07_calculate_CL
def M09_calculate_CD(Vt, FlightData_Geometric_c, alpha, alpha_dot, alpha_crit, X_flow_sep):
    """
    Calculate CL based on the given parameters.

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
    - CD0: zero lift (o) drag coefficient (float)
    - k: induced drag coefficient (float)

    Returns:
    - Calculated value of CD (float)
    """
    #  ,CLo, CLa, CLad, CLq, CLde, CLdf
    # add placeholder values for parameters that should be feed to the function later on
    CLo = 0.36486 * (0.8254/0.36486) # 0.8254
    CLa = 4.9635
    CLad = -1.987
    CLq = 9.8197 * (8/9.8197) # 8
    CLde = 0.0007976 * 1000 # 0.7976
    CLdf = 0.0229 * 10 * (0) # 0

    # 'Cdo' : 0.0061 * 8 * (10/9.81),               # zero lift (o) drag coefficient (Cd)
    # 'k' : 0.0483 * (0.09/0.0483)* (10/9.81),     # induced drag coefficient (k)
    CDo = 0.0061 * 8 * (10/9.81) # 0.049745 
    k = 0.0483 * (0.09/0.0483)* (10/9.81) # 0.091743

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

    # next incorporate the (prestall/poststall) model switch criteria
    if alpha < alpha_crit:

        CL = (CLo + 
              CLa * alpha + 
              CLad * da_a + 
              CLq * q_a + 
              CLde * U[1] + 
              CLdf * U[4])
        
        CD = (CDo + 
              k * CL**2)
        
        # indicate that the poststall model is off (0)
        poststall_on = 0

    else:

        # conditions at alpha_crit, assume stable flight
        da_a_cr = 0
        q_a_cr = 0
        U_1_trim_cr = 0 # assume 0 for now, modify later on
        U_4_trim_cr = 0 # assume 0 for now, modify later on

        CDo_cr = 0.0061 * 8 * (10/9.81) # 0.049745 
        k_cr = 0.0483 * (0.09/0.0483)* (10/9.81) # 0.091743

        CL_cr = (CLo + 
                 CLa * alpha_crit + 
                 CLad * da_a_cr + 
                 CLq * q_a_cr + 
                 CLde * U_1_trim_cr + 
                 CLdf * U_4_trim_cr)
        
        CD_cr = (CDo_cr + 
                 k_cr * CL_cr**2)      
        
        CLawb_cr = CLa
        CLad_cr = CLad
        CLq_cr = CLq

        # Page 17/25 in Wang Paper
        # drag coefficient correction term
        CDX = 0.18 * (X_flow_sep**2) - 0.16 * (X_flow_sep) + 0.213

        # assume linear mapping of stick position to control surface deflection angle (as in the default flight simulation physics logic)

        CL = (CL_cr + 
              CLawb_cr * (((1 + np.sqrt(X_flow_sep))/2)**2) * (alpha - alpha_crit) + 
              CLad_cr * da_a + 
              CLq_cr * q_a + 
              CLde * (U[1] - U_1_trim_cr) + 
              CLdf * (U[4] - U_4_trim_cr))
        
        CD = (CD_cr +
              k * (CL - CL_cr)**2 +
              CDX * (1 - X_flow_sep)) 
        
        # indicate that the poststall model is on (1)
        poststall_on = 1

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

    return CD
