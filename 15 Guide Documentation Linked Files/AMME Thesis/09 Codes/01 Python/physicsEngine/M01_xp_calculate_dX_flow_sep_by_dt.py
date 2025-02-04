import numpy as np

# global_time = 0
# global_alpha = 0
# global_dt = 0.01
# X_flow_sep = 0
# dt = 0.01

def M01_xp_calculate_dX_flow_sep_by_dt(global_time, global_dt, X_flow_sep, global_alpha, alpha, FlightData_Geometric_c, Vt):
    """
    Calculate dX_flow_sep_by_dt based on the given parameters.

    Parameters:
    - global_time: Current global time (float or int)
    - global_dt: Time step (float or int)
    - alpha: Current alpha angle (float or int)
    - global_alpha: Previous global alpha angle (float or int)
    - FlightData_Geometric_c: Geometric parameter (float or int)
    - Vt: Velocity parameter (float or int)

    Returns:
    - Calculated value of dX_flow_sep_by_dt (float)
    """
    # ===================================================================
    # Richard Additional Code Block 1 (Start)
    # ===================================================================
    # # Declare that we will use the global clock
    # global global_time
    # note: initialization lasts around 2.3 seconds, or 230 iterations, during which Vt and term_5 can be nan
    init_wait_time = 5

    # # Calculate the time rate of change of X_flow_sep
    # global global_alpha
    # # use the global dt determined by xplane
    # global dt # declare dt as a global variable so that I could print it in calculate_aero_forces() (Maybe useless now)

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
    global_alpha = alpha # this will in the output along with dX_flow_sep_by_dt

    a_1 = (1.00) * 22.5
    alpha_star = 20 * (np.pi / 180)
    # tau_1 = 11.93 * (self.FlightData_Geometric_c / Vt)
    # tau_2 = 6.66 * (self.FlightData_Geometric_c / Vt)
    # tau_1 = 11.93 * (self.FlightData_Geometric_c / 52.65)
    # tau_2 = 6.66 * (self.FlightData_Geometric_c / 52.65)
    # tau_1 = 11.93 * (self.FlightData_Geometric_c / (Vt if Vt == Vt else 52.65)) #uses the trick that nan is not equal to nan
    # tau_2 = 6.66 * (self.FlightData_Geometric_c / (Vt if Vt == Vt else 52.65))
    # tau_1 = (1.00) * 11.93 * (self.FlightData_Geometric_c / (Vt if global_time > init_wait_time else 52.65))
    # tau_2 = (1.00) * 6.66 * (self.FlightData_Geometric_c / (Vt if global_time > init_wait_time else 52.65))
    tau_1 = (1.00) * 11.93 * (FlightData_Geometric_c / (Vt if global_time > init_wait_time else 52.65)) # determines how fast X converges if alpha is constant
    tau_2 = (1.00) * 6.66 * (FlightData_Geometric_c / (Vt if global_time > init_wait_time else 52.65)) # does not impact how X behaves if alpha is constant
    # global X_flow_sep  # Declare that we will use the global variable
    global dX_flow_sep_by_dt
    # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep
    # dX_flow_sep_by_dt = -(1/tau_1) * X_flow_sep + (1/(2*tau_1)) * (1 - np.tanh(a_1 * (alpha - alpha_star))) * (alpha - tau_2 * alpha_dot)
    
    # global dt_time_step
    # dt_time_step = 0.01


    # # old term naming convention
    # term_0 = 1/(2*tau_1)
    # term_1 = a_1 * (alpha - alpha_star)
    # # term_2 = 1 - np.tanh(a_1 * (alpha - alpha_star))
    # term_2 = 1 - np.tanh(a_1 * ((alpha if global_time > init_wait_time else 0.0873) - alpha_star))
    # # term_3 = alpha - tau_2 * alpha_dot
    # term_3 = (alpha if global_time > init_wait_time else 0.0873) - tau_2 * (alpha_dot if global_time > init_wait_time else 0)
    # term_4 = -(1 / tau_1)
    # term_5 = term_0 * term_2 * term_3

    # Wrong Version:
    # # new term naming convention
    # term_0 = -(1 / tau_1)

    # term_5 = a_1 * (alpha - alpha_star)

    # term_2 = 1/(2*tau_1)
    # # term_3 = 1 - np.tanh(a_1 * (alpha - alpha_star))
    # term_3 = 1 - np.tanh(a_1 * ((alpha if global_time > init_wait_time else 0.0873) - alpha_star))
    # # term_4 = alpha - tau_2 * alpha_dot
    # term_4 = (alpha if global_time > init_wait_time else 0.0873) - tau_2 * (alpha_dot if global_time > init_wait_time else 0)

    # term_1 = term_2 * term_3 * term_4

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    whether_insert_delft_paper_numbers = 0
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if whether_insert_delft_paper_numbers == 1:
        a_1 = 40
        alpha_star = 10 * (np.pi / 180)
        tau_1 = 0.3 # assume c_bar / Vt is constant for now
        tau_2 = 0.4 # assume c_bar / Vt is constant for now
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    whether_investigate_effect_of_parameters = 0
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if whether_investigate_effect_of_parameters == 1:
        # a_1 = 5
        # alpha_star = 50 * (np.pi / 180)
        # tau_1 = 0.3 # assume c_bar / Vt is constant for now
        # tau_2 = 0.4 # assume c_bar / Vt is constant for now
        # tau_1 = 0.1
        # tau_2 = 0.0

        a_1 = 40
        alpha_star = 20 * (np.pi / 180)
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Corrected Version:
    # new term naming convention
    term_0 = -(1 / tau_1)

    term_5 = a_1 * (alpha - alpha_star)

    term_2 = 1/(2*tau_1)

    # term_4 = alpha - tau_2 * alpha_dot
    term_4 = (alpha if global_time > init_wait_time else 0.0873) - tau_2 * (alpha_dot if global_time > init_wait_time else 0)

    # term_3 = 1 - np.tanh(a_1 * (alpha - alpha_star))
    term_3 = 1 - np.tanh(a_1 * ((term_4 if global_time > init_wait_time else 0.0873) - alpha_star))

    term_1 = term_2 * term_3

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    whether_add_ric_multiplier_flag = 0 # just a relic now
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if whether_add_ric_multiplier_flag == 1:
        # record trimmed value of alpha (through flight test) (rad)
        alpha_trimmed = 0.029276918705747886
        # alpha_trimmed = 6.4 * (np.pi / 180)
        ric_multiplier = 1/alpha_trimmed
        # Enforce modification on differential equation to achieve dX_flow_sep_by_dt = 0 when X = 1 (not guarateed to be physical)
        term_3 = ric_multiplier * term_3
        # revise term_1 accordingly
        term_1 = term_2 * term_3 * term_4

    # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + term_0 # ok
    # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + (-1) # ok
    # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + np.random.uniform(-0.2, 2.8) % ok
    # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + (term_5 if global_time > 3 else np.random.uniform(-0.2, 2.8))

    # # with old term naming convention
    # dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + (term_5 if global_time > init_wait_time else 2.5)

    # with new term naming convention
    dX_flow_sep_by_dt = -(1 / tau_1) * X_flow_sep + (term_1 if global_time > init_wait_time else 2.5)

    # #---------------------------------------
    # dX_flow_sep_by_dt = - dX_flow_sep_by_dt
    # #---------------------------------------
    # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ # manual enforcements

    # flag for whether manually enforce X_flow_sep stay between 0 and 1. no (0), yes (1).
    manual_enforcement_flag = 0 # just a relic now

    if manual_enforcement_flag == 1:

        # prevent unphysically large X_flow_sep due to large global_dt
        if global_time <= init_wait_time or global_dt >= 0.25:
            X_flow_sep = 1
            dX_flow_sep_by_dt = 0

        # make sure X_flow_sep stay between 0 and 1
        if X_flow_sep > 1:
            X_flow_sep = 1
            dX_flow_sep_by_dt = 0

        if X_flow_sep < 0:
            X_flow_sep = 0
            dX_flow_sep_by_dt = 0

    # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # dX_flow_sep_by_dt = term_4 + term_0 * term_2 * term_3 # not ok

    # X_flow_sep = X_flow_sep + dX_flow_sep_by_dt * dt_time_step  # Modify the global variable

    # ==================================================================================================
    # script for calculating X_flow_sep from dX_flow_sep_by_dt and dt
    # ==================================================================================================
    # X_flow_sep_file_path = "X_flow_sep.txt"
    # with open(X_flow_sep_file_path, "r") as file:
    #     X_flow_sep = float(file.read().strip()) + dX_flow_sep_by_dt * dt_time_step
    # # Write the new value back to the file without formatting to preserve accuracy
    # with open(X_flow_sep_file_path, "w") as file:
    #     file.write(str(X_flow_sep))
    # print("Current value:", X_flow_sep)
    # ==================================================================================================

    # ==================================================================================================
    # script for locating the home directory
    # ==================================================================================================
    # # Data to write
    # text_to_write = "Hello, world!"

    # # Open the file in write mode
    # with open('newfile.txt', 'w') as file:
    #     file.write(text_to_write)
    # ==================================================================================================

    # ==================================================================================================
    # script for writing X_flow_sep for checking
    # ==================================================================================================
    whether_write_separate_files_flag = 0
    if whether_write_separate_files_flag == 1:
        # Open the file in write mode
        with open('Vt.txt', 'w') as file:
            file.write(str(Vt))
        with open('X_flow_sep.txt', 'w') as file:
            file.write(str(X_flow_sep))
        with open('tau_1.txt', 'w') as file:
            file.write(str(tau_1))
        with open('tau_2.txt', 'w') as file:
            file.write(str(tau_2))
        with open('alpha_star.txt', 'w') as file:
            file.write(str(alpha_star))
        with open('global_alpha.txt', 'w') as file:
            file.write(str(global_alpha))
        with open('alpha.txt', 'w') as file:
            file.write(str(alpha))
        with open('alpha_dot.txt', 'w') as file:
            file.write(str(alpha_dot))
        with open('dX_flow_sep_by_dt.txt', 'w') as file:
            file.write(str(dX_flow_sep_by_dt))
        with open('richard_checked_terms\\term_0.txt', 'w') as file:
            file.write(str(term_0))
        with open('richard_checked_terms\\term_1.txt', 'w') as file:
            file.write(str(term_1))
        with open('richard_checked_terms\\term_2.txt', 'w') as file:
            file.write(str(term_2))
        with open('richard_checked_terms\\term_3.txt', 'w') as file:
            file.write(str(term_3))
        with open('richard_checked_terms\\term_4.txt', 'w') as file:
            file.write(str(term_4))
        with open('richard_checked_terms\\term_5.txt', 'w') as file:
            file.write(str(term_5))
        with open('richard_checked_terms\\global_time.txt', 'w') as file:
            file.write(str(global_time))
        with open('richard_checked_terms\\term_5_history.txt', 'a') as file:
            file.write(str(term_5) + '\n')
        with open('richard_checked_terms\\Vt_history.txt', 'a') as file:
            file.write(str(Vt) + '\n')
        with open('richard_checked_terms\\alpha_history.txt', 'a') as file:
            file.write(str(alpha) + '\n')
        with open('richard_checked_terms\\alpha_dot_history.txt', 'a') as file:
            file.write(str(alpha_dot) + '\n')
        with open('richard_checked_terms\\X_flow_sep_history.txt', 'a') as file:
            file.write(str(X_flow_sep) + '\n')
        
    # # write combined    
    # # with open('richard_checked_terms\\combined_history.txt', 'a') as file:
    # #     file.write(f"{alpha_dot}\t{X_flow_sep}\n")
    # file_path = 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\combined_history.txt'
    # # column width in terms of character string length
    # cwidth = 25

    # # Function to format the columns
    # def format_line(global_time, global_dt, dX_flow_sep_by_dt, term_0, X_flow_sep, term_1, term_2, term_3, term_4, term_5, alpha, alpha_dot, a_1, alpha_star, tau_1, tau_2, Vt, cwidth):
    #     # Format global_time to 5 significant figures
    #     formatted_global_time = f"{global_time:.5g}"
    #     return f"{formatted_global_time:<{cwidth}}{global_dt:<{cwidth}}{dX_flow_sep_by_dt:<{cwidth}}{term_0:<{cwidth}}{X_flow_sep:<{cwidth}}{term_1:<{cwidth}}{term_2:<{cwidth}}{term_3:<{cwidth}}{term_4:<{cwidth}}{term_5:<{cwidth}}{alpha:<{cwidth}}{alpha_dot:<{cwidth}}{a_1:<{cwidth}}{alpha_star:<{cwidth}}{tau_1:<{cwidth}}{tau_2:<{cwidth}}{Vt:<{cwidth}}\n"

    # try:
    #     # Check if the file is empty
    #     with open(file_path, 'r') as file:
    #         is_empty = file.read().strip() == ""
    # except FileNotFoundError:
    #     # If the file doesn't exist, it will be considered empty
    #     is_empty = True

    # # Write the column titles if the file is empty
    # with open(file_path, 'a') as file:
    #     if is_empty:
    #         file.write(f"{'global_time[s]':<{cwidth}}{'global_dt[s]':<{cwidth}}{'dX_flow_sep_by_dt[-/s]':<{cwidth}}{'term_0[-/s]':<{cwidth}}{'X_flow_sep[-]':<{cwidth}}{'term_1[-/s]':<{cwidth}}{'term_2[-/s]':<{cwidth}}{'term_3[-]':<{cwidth}}{'term_4[-]':<{cwidth}}{'term_5[-]':<{cwidth}}{'alpha[rad]':<{cwidth}}{'alpha_dot[rad/s]':<{cwidth}}{'a_1[-]':<{cwidth}}{'alpha_star[rad]':<{cwidth}}{'tau_1[s]':<{cwidth}}{'tau_2[s]':<{cwidth}}{'Vt[m/s]':<{cwidth}}\n")
    #     file.write(format_line(global_time, global_dt, dX_flow_sep_by_dt, term_0, X_flow_sep, term_1, term_2, term_3, term_4, term_5, alpha, alpha_dot, a_1, alpha_star, tau_1, tau_2, Vt, cwidth))
    # ===================================================================
    # Richard Additional Code Block 1 (End)
    # ===================================================================

    # print(global_time)

    # return the value of "dX_flow_sep_by_dt" at the very end of this function
    return dX_flow_sep_by_dt, global_alpha

