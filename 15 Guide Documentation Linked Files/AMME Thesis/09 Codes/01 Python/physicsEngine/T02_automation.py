import subprocess
import numpy as np

def run_script(script_name):
    try:
        # Runs the script
        result = subprocess.run(['python', script_name], capture_output=True, text=True)
        print(f"Output of {script_name}:\n{result.stdout}")
        if result.stderr:
            print(f"Errors from {script_name}:\n{result.stderr}")
    except Exception as e:
        print(f"An error occurred while running {script_name}: {e}")

def main():
    scripts_set_1 = [
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\T01_investigate_X_flow_sep_model_behavior.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R10_manually_format_combined_history.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R11_plot_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R12_count_number_of_lines_in_txt.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R13_plot_alpha_vs_global_time_and_find_average_alpha.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R14_plot_alpha_and_X_flow_sep_vs_global_time.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R15_plot_alpha_alpha_crit_and_X_flow_sep_vs_global_time.py',
    ] # List scripts to be automated

    scripts_set_2 = [
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\T01_investigate_X_flow_sep_model_behavior.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R10_manually_format_combined_history.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R11_plot_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R12_count_number_of_lines_in_txt.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R13_plot_alpha_vs_global_time_and_find_average_alpha.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R14_plot_alpha_and_X_flow_sep_vs_global_time.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R15_plot_alpha_alpha_crit_and_X_flow_sep_vs_global_time.py',
    ] # List scripts to be automated

    scripts_set_3 = [
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\T01_investigate_X_flow_sep_model_behavior.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R10_manually_format_combined_history.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R11_plot_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R12_count_number_of_lines_in_txt.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R13_plot_alpha_vs_global_time_and_find_average_alpha.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R14_plot_alpha_and_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R15_plot_alpha_alpha_crit_and_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R17_plot_CL_vs_global_time.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R18_plot_CL_vs_global_time_with_poststall_model_indicator.py',
    ] # List scripts to be automated

    scripts_set_4 = [
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\T01_investigate_X_flow_sep_model_behavior.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R10_manually_format_combined_history.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R11_plot_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R12_count_number_of_lines_in_txt.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R13_plot_alpha_vs_global_time_and_find_average_alpha.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R14_plot_alpha_and_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R15_plot_alpha_alpha_crit_and_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R17_plot_CL_vs_global_time.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R19_plot_CL_vs_alpha.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R18_plot_CL_vs_global_time_with_poststall_model_indicator.py',
    ] # List scripts to be automated

    scripts_set_5 = [
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\T01_investigate_X_flow_sep_model_behavior.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R10_manually_format_combined_history.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R11_plot_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R12_count_number_of_lines_in_txt.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R13_plot_alpha_vs_global_time_and_find_average_alpha.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R14_plot_alpha_and_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R15_plot_alpha_alpha_crit_and_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R17_plot_CL_vs_global_time.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R20_plot_Cm_vs_alpha.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R19_plot_CL_vs_alpha.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R18_plot_CL_vs_global_time_with_poststall_model_indicator.py',
    ] # List scripts to be automated

    scripts_set_6 = [
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\T01_investigate_X_flow_sep_model_behavior.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R10_manually_format_combined_history.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R11_plot_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R12_count_number_of_lines_in_txt.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R13_plot_alpha_vs_global_time_and_find_average_alpha.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R14_plot_alpha_and_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R15_plot_alpha_alpha_crit_and_X_flow_sep_vs_global_time.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R17_plot_CL_vs_global_time.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R21_plot_CD_vs_alpha.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R20_plot_Cm_vs_alpha.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R19_plot_CL_vs_alpha.py',
    # 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R18_plot_CL_vs_global_time_with_poststall_model_indicator.py',
    'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\R22_plot_Cm_CD_CL_vs_global_time.py',
    ] # List scripts to be automated

    #+++++++++++++++++++++++++++++++++++++++++++++
    # Select which set of scripts to run
    #+++++++++++++++++++++++++++++++++++++++++++++
    scripts = scripts_set_6
    #+++++++++++++++++++++++++++++++++++++++++++++

    for script in scripts:
        run_script(script)


# Enter test case number below
#++++++++++++++++++++++++++++++++++++++++++++
selected_test_case = 8
#++++++++++++++++++++++++++++++++++++++++++++

# Constant alpha with time
#=====================================================================
# note: since alpha_star is the value of stable alpha when X = 0.5,
# if we set the value of alpha to alpha_star here, then X must be 0.5
# otherwise: 1. constants are wrong
#            2. code is wrong
#            3. model is wrong
#=====================================================================
if selected_test_case == 1:
    # Parameters
    time_steps = 1000
    dt = 0.01
    FlightData_Geometric_c = 0.99
    Vt = 40.0
    alpha_trimmed = 0.029276918705747886 # record trimmed value of alpha (through flight test) (rad)
    alpha_trimmed_in_deg = 0.029276918705747886 * (180/np.pi) # record trimmed value of alpha (through flight test) (deg)
    print(f"Alpha trimmed: {alpha_trimmed_in_deg} deg")
    # Enter the constant value of alpha below
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_constant_value_in_deg = 18
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_constant_value = alpha_constant_value_in_deg * (np.pi/180)
    alpha_series = np.full(time_steps, alpha_constant_value)  # Constant alpha of 0.5

# Linear alpha with time
if selected_test_case == 2:
    # Parameters
    time_steps = 1000
    dt = 0.01
    FlightData_Geometric_c = 0.99
    Vt = 40.0
    # Enter the maximum value of alpha below
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value_in_deg = 25
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value = alpha_maximum_value_in_deg * (np.pi/180)
    alpha_series = np.linspace(0, alpha_maximum_value, time_steps)  # Linear increase from 0 to 1

# Sinusoidal alpha with time
if selected_test_case == 3:
    # Parameters
    time_steps = 1000
    dt = 0.01
    FlightData_Geometric_c = 0.99
    Vt = 40.0
    # Enter the oscillation magnitude value of alpha below
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_oscillation_magnitude_value_in_deg = 25
    frequency_multiplier = 5
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_oscillation_magnitude_value = alpha_oscillation_magnitude_value_in_deg * (np.pi/180)
    alpha_series = 0.029276918705747886 + alpha_oscillation_magnitude_value * np.sin(np.linspace(0, 2 * np.pi * frequency_multiplier, time_steps))  # Oscillates between 0 and 1


# Random flucturating alpha with time
if selected_test_case == 4:
    # Parameters
    time_steps = 1000
    dt = 0.01
    FlightData_Geometric_c = 0.99
    Vt = 40.0
    np.random.seed(42)  # Seed for reproducibility
    mean = 0.029276918705747886
    fluctuation_size = 0.015  # Smaller fluctuation size
    alpha_series = np.random.uniform(mean - fluctuation_size, mean + fluctuation_size, time_steps)


# Step increase of alpha with time
if selected_test_case == 5:
    # Parameters
    time_steps = 1000
    dt = 0.01
    FlightData_Geometric_c = 0.99
    Vt = 40.0
    # Enter the initial constant value of alpha below
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_initial_constant_value_in_deg = 20
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_initial_constant_value = alpha_initial_constant_value_in_deg * (np.pi/180)
    initial_value = alpha_initial_constant_value     # Value before the step change
    step_value = initial_value * (50/20)            # Value after the step change
    step_time = 500 # at which time step that it changes
    alpha_series = np.full(time_steps, initial_value)  # Start with all values at initial_value
    alpha_series[step_time:] = step_value              # Change to step_value at step_time


# Test case 6: Sawtooth pattern
if selected_test_case == 6:
    # Parameters
    time_steps = 1000
    dt = 0.01
    FlightData_Geometric_c = 0.99
    Vt = 40.0
    # Enter the maximum value of alpha below
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value_in_deg = 25
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value = alpha_maximum_value_in_deg * (np.pi/180)
    
    # Sawtooth pattern
    repeats = 5
    single_sawtooth = np.linspace(0, alpha_maximum_value, time_steps // repeats)
    alpha_series = np.tile(single_sawtooth, repeats)
    alpha_series = np.append(alpha_series, 0)  # Ensure it returns to initial value at the end


# Test case 7: Sawtooth pattern with gradual return
if selected_test_case == 7:
    # Parameters
    time_steps = 1000
    dt = 0.01
    FlightData_Geometric_c = 0.99
    Vt = 40.0
    # Enter the maximum value of alpha below
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value_in_deg = 25
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value = alpha_maximum_value_in_deg * (np.pi/180)
    
    # Sawtooth pattern with gradual return
    repeats = 5
    single_sawtooth_up = np.linspace(0, alpha_maximum_value, time_steps // (2 * repeats))
    single_sawtooth_down = np.linspace(alpha_maximum_value, 0, time_steps // (2 * repeats))
    single_sawtooth = np.concatenate((single_sawtooth_up, single_sawtooth_down))
    alpha_series = np.tile(single_sawtooth, repeats)
    
    # Ensure it returns to initial value at the end
    if len(alpha_series) < time_steps:
        alpha_series = np.append(alpha_series, np.linspace(alpha_series[-1], 0, time_steps - len(alpha_series)))

# Test case 8: Single Sawtooth pattern (quasi-static)
if selected_test_case == 8:
    # Parameters
    time_steps = 3000
    dt = 0.1
    FlightData_Geometric_c = 0.99
    Vt = 40.0
    # Enter the maximum value of alpha below
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value_in_deg = 90 # default: 60
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value = alpha_maximum_value_in_deg * (np.pi/180)
    
    # Sawtooth pattern with gradual return
    repeats = 1
    single_sawtooth_up = np.linspace(0, alpha_maximum_value, time_steps // (2 * repeats))
    single_sawtooth_down = np.linspace(alpha_maximum_value, 0, time_steps // (2 * repeats))
    single_sawtooth = np.concatenate((single_sawtooth_up, single_sawtooth_down))
    alpha_series = np.tile(single_sawtooth, repeats)
    
    # Ensure it returns to initial value at the end
    if len(alpha_series) < time_steps:
        alpha_series = np.append(alpha_series, np.linspace(alpha_series[-1], 0, time_steps - len(alpha_series)))

# Test case 9: Single Sawtooth pattern (fast cycle)
if selected_test_case == 9:
    # Parameters
    time_steps = 300 # 10 times faster
    dt = 0.1
    FlightData_Geometric_c = 0.99
    Vt = 40.0
    # Enter the maximum value of alpha below
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value_in_deg = 90 # default: 60
    #++++++++++++++++++++++++++++++++++++++++++++
    alpha_maximum_value = alpha_maximum_value_in_deg * (np.pi/180)
    
    # Sawtooth pattern with gradual return
    repeats = 1
    single_sawtooth_up = np.linspace(0, alpha_maximum_value, time_steps // (2 * repeats))
    single_sawtooth_down = np.linspace(alpha_maximum_value, 0, time_steps // (2 * repeats))
    single_sawtooth = np.concatenate((single_sawtooth_up, single_sawtooth_down))
    alpha_series = np.tile(single_sawtooth, repeats)
    
    # Ensure it returns to initial value at the end
    if len(alpha_series) < time_steps:
        alpha_series = np.append(alpha_series, np.linspace(alpha_series[-1], 0, time_steps - len(alpha_series)))

# use Ctrl + c if you want to stop the running Python script


# Run the automation
if __name__ == "__main__":
    main()
