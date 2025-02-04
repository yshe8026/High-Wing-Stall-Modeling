import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter1d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.interpolate import interp1d
import os

'''Each Row in the Flight Data Matrix is a data entry at a specific moment in time'''


# Please read the following options and choose what output is wanted
#################################################################################################
### User preferences:

## Search for "User input" using "Shift + F" will show more input fields throughout the script  

## Print options (0: no, 1: yes)

# input field
whether_to_enable_extra_wide_and_full_print_setting = 1
# input field
whether_to_print_matrix = 0
whether_to_print_full_matrix = 0



## Plot options (0: no, 1: yes)

# input field
whether_to_plot_X_X_dot = 0
# input field
whether_to_plot_3D_flight_path = 0
# input field
whether_to_create_3D_flight_path_animation = 0
# input field
whether_to_apply_smoothing_for_Xdot_data = 0
# input field
whether_to_plot_U = 0
# input field
whether_to_plot_stage2_combined_history = 1
#################################################################################################



# Load the flight data matrix

# flight_data_matrix = np.load(xplane11_path + 'Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\flight_data_matrix.npy')
# flight_data_matrix = np.load('C:\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\flight_data_matrix.npy')

# xplane11_path_list = [xplane11_path + '', 'C:\\X-Plane 11\\', 'C:\\X-Plane 11\\']
# xplane11_path = xplane11_path_list[2]

current_work_directory_path = os.getcwd()
known_part_of_path_to_remove = 'Resources\plugins\PythonPlugins'
xplane11_path = current_work_directory_path.replace(known_part_of_path_to_remove, "")
# Print Python work directory
print(current_work_directory_path)
# Print xplane 11 directory on this machine
print(xplane11_path) 
# Path to the directory where the Python file resides
script_dir = os.path.dirname(os.path.abspath(__file__))
print(script_dir)
# Full path to the current Python file
file_path = os.path.abspath(__file__)
print(file_path)

flight_data_matrix = np.load(xplane11_path + 'Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\flight_data_matrix.npy')

# Wide print setting
if whether_to_enable_extra_wide_and_full_print_setting == 1:
    np.set_printoptions(threshold=np.inf, linewidth=500)


# Print the full array (optional)
if whether_to_print_matrix == 1 and whether_to_print_full_matrix == 1:
    np.set_printoptions(threshold=np.inf, linewidth=300)
    print('Flight Data Matrix:')
    print(flight_data_matrix)
    np.set_printoptions(threshold=1000, linewidth=100)

# Extract the time (t) and u, v, w columns
t = flight_data_matrix[:, 0]  # 0th column (time)
u = flight_data_matrix[:, 1]  # 1st column (u)
v = flight_data_matrix[:, 2]  # 2nd column (v)
w = flight_data_matrix[:, 3]  # 3rd column (w)
Vt = np.sqrt(u**2 + v**2 + w**2)
p = flight_data_matrix[:, 4]  # 4th column (p)
q = flight_data_matrix[:, 5]  # 5th column (q)
r = flight_data_matrix[:, 6]  # 6th column (r)
phi = flight_data_matrix[:, 7]  # 7th column (phi)
theta = flight_data_matrix[:, 8]  # 8th column (theta)
psi = flight_data_matrix[:, 9]  # 9th column (psi)
X_e = flight_data_matrix[:, 10] # 10th column (X_e)
Y_e = flight_data_matrix[:, 11] # 11th column (Y_e)
Z_e = flight_data_matrix[:, 12] # 12th column (Z_e)
Z_e_raw = Z_e # record the raw data of Z_e as Z_e_raw
Z_e_raw_initial = Z_e_raw[0]
Z_e_raw_final = Z_e_raw[-1]

u_dot = flight_data_matrix[:, 13]  # 13th column (u_dot)
v_dot = flight_data_matrix[:, 14]  # 14th column (v_dot)
w_dot = flight_data_matrix[:, 15]  # 15th column (w_dot)
Vt_dot = np.sqrt(u_dot**2 + v_dot**2 + w_dot**2)
p_dot = flight_data_matrix[:, 16]  # 16th column (p_dot)
q_dot = flight_data_matrix[:, 17]  # 17th column (q_dot)
r_dot = flight_data_matrix[:, 18]  # 18th column (r_dot)
phi_dot = flight_data_matrix[:, 19]  # 19th column (phi_dot)
theta_dot = flight_data_matrix[:, 20]  # 20th column (theta_dot)
psi_dot = flight_data_matrix[:, 21]  # 21th column (psi_dot)
X_e_dot = flight_data_matrix[:, 22] # 22th column (X_e_dot)
Y_e_dot = flight_data_matrix[:, 23] # 23th column (Y_e_dot)
Z_e_dot = flight_data_matrix[:, 24] # 24th column (Z_e_dot)

U_data_matrix = flight_data_matrix[:, 25:30]
# print(U_data_matrix[-10:-1, :]) # print the last 10 rows
delta_t = flight_data_matrix[:, 25] # 25th column (delta_t) (Throttle input) [-]
delta_e = flight_data_matrix[:, 26] # 25th column (delta_t) (elevator input) [rad]
delta_a = flight_data_matrix[:, 27] # 25th column (delta_t) (aileron input) [rad]
delta_r = flight_data_matrix[:, 28] # 25th column (delta_t) (rudder input) [rad]
delta_f = flight_data_matrix[:, 29] # 25th column (delta_t) (flap input) [rad]

stage2_combined_history_data_matrix = flight_data_matrix[:, 30:44]
print(stage2_combined_history_data_matrix[-10:-1, :]) # print the last 10 rows
# This list of variables are featured back in R21, R22
global_time = flight_data_matrix[:, 30]
global_dt = flight_data_matrix[:, 31]
dX_flow_sep_by_dt = flight_data_matrix[:, 32]
X_flow_sep = flight_data_matrix[:, 33]
alpha = flight_data_matrix[:, 34]
alpha_crit_1 = flight_data_matrix[:, 35]
alpha_crit_2 = flight_data_matrix[:, 36]
alpha_crit_3 = flight_data_matrix[:, 37]
alpha_dot = flight_data_matrix[:, 38]
Vt = flight_data_matrix[:, 39]
CL = flight_data_matrix[:, 40]
poststall_on = flight_data_matrix[:, 41]
Cm = flight_data_matrix[:, 42]
CD = flight_data_matrix[:, 43]



#####################################################################################
# Extract Input
# User input
U_data_matrix_taped = U_data_matrix
t_taped_for_U = t

np.save(xplane11_path + 'Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\U_data_matrix_taped', U_data_matrix_taped) # Relative path approach that works on all machines
np.save(xplane11_path + 'Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\t_taped_for_U', t_taped_for_U) # Relative path approach that works on all machines

#####################################################################################

# Variables used in stage2_combined_history plot 1
# X_flow_sep = flight_data_matrix[:, 33]
# alpha = flight_data_matrix[:, 34]
# alpha_crit_1 = flight_data_matrix[:, 35]
# alpha_crit_2 = flight_data_matrix[:, 36]
# alpha_crit_3 = flight_data_matrix[:, 37]
# CL = flight_data_matrix[:, 40]
# poststall_on = flight_data_matrix[:, 41]
# Cm = flight_data_matrix[:, 42]
# CD = flight_data_matrix[:, 43]
# Vt = flight_data_matrix[:, 39]

# # Mid-point integration method
# def midpoint_integration(y, t):
#     dt = np.diff(t)
#     midpoints = (y[:-1] + y[1:]) / 2
#     integral = np.zeros_like(y)
#     integral[1:] = np.cumsum(midpoints * dt)
#     return integral

# # Integrate using mid-point method
# X_e = midpoint_integration(X_e_dot, t)
# Y_e = midpoint_integration(Y_e_dot, t)
# Z_e = midpoint_integration(Z_e_dot, t)

# Integrate using the cumulative trapezoidal rule
X_e = cumulative_trapezoid(X_e_dot, t, initial=0)  # Integrate X_e_dot to get X_e
Y_e = cumulative_trapezoid(Y_e_dot, t, initial=0)  # Integrate Y_e_dot to get Y_e
Z_e = cumulative_trapezoid(Z_e_dot, t, initial=0)  # Integrate Z_e_dot to get Z_e

# # Basic rectangular integration method (left Riemann sum)
# def rectangular_integration(y, t):
#     dt = np.diff(t)
#     integral = np.zeros_like(y)
#     integral[1:] = np.cumsum(y[:-1] * dt)
#     return integral

# # Integrate using rectangular method
# X_e = rectangular_integration(X_e_dot, t)
# Y_e = rectangular_integration(Y_e_dot, t)
# Z_e = rectangular_integration(Z_e_dot, t)

Z_e = Z_e + Z_e_raw_initial
Z_e_final = Z_e[-1]
Z_e = Z_e + (t/t[-1]) * (Z_e_raw_final - Z_e_final) # Need a mysterious correction, might be the deviation between altitude of default and custom physics since part of the journey is on default physics (t< approx 10s)


# Invert Z_e to obtain a more intuitive alitude expression
Z_e_raw = - Z_e_raw
Z_e = -Z_e


# Convert u from m/s to knots (1 m/s = 1.943844 knots)
u = 1.943844 * u # kn
v = 1.943844 * v # kn
w = 1.943844 * w # kn
Vt = 1.943844 * Vt # kn
p = (180/np.pi) * p
q = (180/np.pi) * q
r = (180/np.pi) * r
phi = (180/np.pi) * phi
theta = (180/np.pi) * theta
psi = (180/np.pi) * psi

delta_e = (180/np.pi) * delta_e
delta_a = (180/np.pi) * delta_a
delta_r = (180/np.pi) * delta_r
delta_f = (180/np.pi) * delta_f

alpha = (180/np.pi) * alpha
alpha_crit_1 = (180/np.pi) * alpha_crit_1
alpha_crit_2 = (180/np.pi) * alpha_crit_2
alpha_crit_3 = (180/np.pi) * alpha_crit_3


# Custom colors for the plot
color_u = [61/255, 142/255, 190/255]  # Customized Blue for u
color_v = [255/255, 36/255, 10/255]   # Customized Red for v
color_w = [21/255, 110/255, 72/255]   # Customized Green for w
color_p = [255/255, 159/255, 64/255]  # Customized Orange for p
color_q = [75/255, 192/255, 192/255]  # Customized Cyan for q
color_r = [153/255, 102/255, 255/255]  # Customized Purple for r
color_phi = [255/255, 127/255, 80/255]    # Coral for phi (roll)
color_theta = [100/255, 149/255, 237/255] # Cornflower blue for theta (pitch)
color_psi = [50/255, 205/255, 50/255]     # Lime green for psi (yaw)
color_X_e = [255/255, 105/255, 180/255]  # Hot Pink for X_e
color_Y_e = [30/255, 144/255, 255/255]   # Dodger Blue for Y_e
color_Z_e = [34/255, 139/255, 34/255]    # Forest Green for Z_e

color_u_dot = [255/255, 69/255, 0/255]   # Red-Orange for u_dot
color_v_dot = [34/255, 139/255, 34/255]  # Forest Green for v_dot
color_w_dot = [65/255, 105/255, 225/255] # Royal Blue for w_dot
color_p_dot = [255/255, 127/255, 80/255]   # Coral for p_dot
color_q_dot = [100/255, 149/255, 237/255]  # Cornflower Blue for q_dot
color_r_dot = [75/255, 0/255, 130/255]     # Indigo for r_dot
color_phi_dot = [255/255, 165/255, 0/255]   # Orange for phi_dot
color_theta_dot = [0/255, 191/255, 255/255] # Deep Sky Blue for theta_dot
color_psi_dot = [220/255, 20/255, 60/255]   # Crimson for psi_dot
color_X_e_dot = [255/255, 105/255, 180/255]  # Hot Pink for X_e_dot
color_Y_e_dot = [30/255, 144/255, 255/255]   # Dodger Blue for Y_e_dot
color_Z_e_dot = [34/255, 139/255, 34/255]    # Forest Green for Z_e_dot

# Colors for plotting stage2_combined_history
color1 = [255/255, 36/255, 10/255]  # Customized Red
color2 = [61/255, 142/255, 190/255]  # Customized Blue
color3 = [21/255, 110/255, 72/255]   # Customized Green
color4 = [217/255, 11/255, 125/255]  # Customized Pink
color5 = [228/255, 100/255, 0/255]   # Customized Orange
color6 = [120/255, 81/255, 169/255]  # Customized Purple
color7 = [239/255, 68/255, 68/255]   # Contrast Red
color8 = [0/255, 159/255, 117/255]   # Contrast Green
color9 = [136/255, 198/255, 237/255] # Contrast Cyan


if whether_to_apply_smoothing_for_Xdot_data == 1:
    # Apply Gaussian smoothing (adjust sigma for more or less smoothing)
    sigma_value = 10 # 3 is the statistical recommendation for averaging over 3-7 time steps

    u_dot = gaussian_filter1d(u_dot, sigma=sigma_value)
    v_dot = gaussian_filter1d(v_dot, sigma=sigma_value)
    w_dot = gaussian_filter1d(w_dot, sigma=sigma_value)
    Vt_dot = gaussian_filter1d(Vt_dot, sigma=sigma_value)
    p_dot = gaussian_filter1d(p_dot, sigma=sigma_value)
    q_dot = gaussian_filter1d(q_dot, sigma=sigma_value)
    r_dot = gaussian_filter1d(r_dot, sigma=sigma_value)
    phi_dot = gaussian_filter1d(phi_dot, sigma=sigma_value)
    theta_dot = gaussian_filter1d(theta_dot, sigma=sigma_value)
    psi_dot = gaussian_filter1d(psi_dot, sigma=sigma_value)
    X_e_dot = gaussian_filter1d(X_e_dot, sigma=sigma_value)
    Y_e_dot = gaussian_filter1d(Y_e_dot, sigma=sigma_value)
    Z_e_dot = gaussian_filter1d(Z_e_dot, sigma=sigma_value)


# ################################################################################################################################################################################
# ################################################################################################################################################################################
# # import matplotlib.pyplot as plt

# def plot_variables_vs_variable_1_y_axis(t, variables, labels, colors, linestyles, xplane11_path, xlim = [0, 50], ylim = None, title=None, xlabel='Time [s]', ylabel=None, legend_loc_relative_to_anchor='upper right', save_path=None, legend_anchor_loc = (0.9, 0.9)):
#     """
#     Plots multiple variables against time.
    
#     Parameters:
#     t (array-like): Time variable (x-axis).
#     variables (list of array-like): List of variables to plot against t (y-axes).
#     labels (list of str): List of labels for each variable.
#     colors (list of str): List of colors for each variable.
#     linestyles (list of str): List of linestyles for each variable.
#     title (str, optional): Title of the plot. Default is None.
#     xlabel (str, optional): Label for the x-axis. Default is 'Time [s]'.
#     ylabel (str, optional): Label for the y-axis. Default is None.
#     legend_loc (str, optional): Location of the legend. Default is 'upper right'.
#     save_path (str, optional): Path to save the plot image. Default is None.
#     """
    
#     # Create the plot
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot each variable
#     for var, label, color, linestyle in zip(variables, labels, colors, linestyles):
#         ax1.plot(t, var, label=label, color=color, linestyle=linestyle, linewidth=1.5)
    
#     # Set axis labels
#     ax1.set_xlabel(xlabel)
#     ax1.set_xlim(xlim)
#     if ylabel:
#         ax1.set_ylabel(ylabel, color='k')
#         ax1.set_ylim(ylim)
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)
    
#     # Add title if provided
#     if title:
#         plt.title(title)
    
#     # Add legend
#     fig.legend(loc=legend_loc_relative_to_anchor, bbox_to_anchor=legend_anchor_loc)
    
#     # Adjust layout
#     fig.tight_layout()
    
#     # Save the figure if a save path is provided
#     if save_path:
#         plt.savefig(xplane11_path + 'richard_checked_terms\\' + save_path)
    
#     # Show the plot
#     plt.show()

# # User input
# whether_show_plot_variables_vs_variable_1_y_axis_example_A = 1

# if whether_show_plot_variables_vs_variable_1_y_axis_example_A == 1:

#     # Example tests:

#     variables = [u, v, w, Vt]
#     labels = [r'$u$ vs. $t$', r'$v$ vs. $t$', r'$w$ vs. $t$', 'TAS vs. $t$']
#     colors = ['b', 'r', 'g', 'k']
#     linestyles = ['-', '-', '-', ':']

#     plot_variables_vs_variable_1_y_axis(t, variables, labels, colors, linestyles, xplane11_path, xlim = [0, 50], ylim = None, 
#                         title=r'$u$, $v$, $w$ and TAS vs. $t$', xlabel='Time [s]', ylabel='Airspeeds [kn]', save_path='custom_plot_function_example.png')

#     # variables = [1/global_dt, poststall_on * 5]
#     # labels = [r'Physics Frame Rate vs. $t$', r'poststall$\_$on vs. $t$']
#     # colors = ['b', 'r', 'g', 'k']
#     # linestyles = ['-', '-', '-', ':']
#     # plot_variables_vs_variable_1_y_axis(t, variables, labels, colors, linestyles, xplane11_path, xlim = [0, 50], ylim = [0, 105], 
#     #                     title=r'Custom Physics Engine Frame Rate (FPS)', xlabel='Time [s]', ylabel='Physics Frame Rate [Hz or fps]', save_path='Physics Frame Rate vs t.png')
    
#     # Test motion chair 2 frame rate
#     variables = [1/global_dt, poststall_on * 20]
#     labels = [r'Physics Frame Rate vs. $t$', r'poststall$\_$on vs. $t$']
#     colors = ['b', 'r', 'g', 'k']
#     linestyles = ['-', '-', '-', ':']
#     plot_variables_vs_variable_1_y_axis(t, variables, labels, colors, linestyles, xplane11_path, xlim = [0, 50], ylim = [0, 250], 
#                         title=r'Custom Physics Engine Frame Rate (FPS)', xlabel='Time [s]', ylabel='Physics Frame Rate [Hz or fps]', save_path='Physics Frame Rate vs t.png')
    
# # User input
# whether_show_plot_variables_vs_variable_1_y_axis_example_B = 0

# if whether_show_plot_variables_vs_variable_1_y_axis_example_B == 1:
    
#     variables = [u_dot]
#     labels = [r'$\dot{u}$ vs. $u$', 'xxx']
#     colors = [color_u, 'r', 'g', 'k']
#     linestyles = ['-', '-', '-', ':']
#     plot_variables_vs_variable_1_y_axis(u, variables, labels, colors, linestyles, xplane11_path, xlim = None, ylim = None, 
#                         title=r'$\dot{u}$ vs. $u$', xlabel=r'$u$ [m/s]', ylabel=r'$\dot{u}$ [m/s$^{2}$]', save_path='u_dot vs u.png')
    
#     variables = [v_dot]
#     labels = [r'$\dot{v}$ vs. $v$', 'xxx']
#     colors = [color_v, 'r', 'g', 'k']
#     linestyles = ['-', '-', '-', ':']
#     plot_variables_vs_variable_1_y_axis(v, variables, labels, colors, linestyles, xplane11_path, xlim = None, ylim = None, 
#                         title=r'$\dot{v}$ vs. $v$', xlabel=r'$v$ [m/s]', ylabel=r'$\dot{v}$ [m/s$^{2}$]', save_path='v_dot vs v.png')
    
#     variables = [w_dot]
#     labels = [r'$\dot{w}$ vs. $w$', 'xxx']
#     colors = [color_w, 'r', 'g', 'k']
#     linestyles = ['-', '-', '-', ':']
#     plot_variables_vs_variable_1_y_axis(w, variables, labels, colors, linestyles, xplane11_path, xlim = None, ylim = None, 
#                         title=r'$\dot{w}$ vs. $w$', xlabel=r'$w$ [m/s]', ylabel=r'$\dot{w}$ [m/s$^{2}$]', save_path='w_dot vs w.png')

# # Now we can plot any variables against each other
# ################################################################################################################################################################################
# ################################################################################################################################################################################
# # import matplotlib.pyplot as plt

# import matplotlib.pyplot as plt

# def plot_variables_vs_variable_2_y_axis(t, left_variables, right_variables, left_labels, right_labels, left_colors, right_colors, left_linestyles, right_linestyles, 
#                                         xplane11_path, xlim=[0, 50], left_ylim=None, right_ylim=None, title=None, xlabel='Time [s]', 
#                                         left_ylabel=None, right_ylabel=None, legend_loc_relative_to_anchor='upper right', 
#                                         save_path=None, legend_anchor_loc=(0.9, 0.9), left_yaxis_color='k', right_yaxis_color='k'):
#     """
#     Plots multiple variables against time on two y-axes with optional colored y-axes.
    
#     Parameters:
#     t (array-like): Time variable (x-axis).
#     left_variables (list of array-like): List of variables to plot against t on the left y-axis.
#     right_variables (list of array-like): List of variables to plot against t on the right y-axis.
#     left_labels (list of str): List of labels for each variable on the left y-axis.
#     right_labels (list of str): List of labels for each variable on the right y-axis.
#     left_colors (list of str): List of colors for each variable on the left y-axis.
#     right_colors (list of str): List of colors for each variable on the right y-axis.
#     left_linestyles (list of str): List of linestyles for each variable on the left y-axis.
#     right_linestyles (list of str): List of linestyles for each variable on the right y-axis.
#     xlim (tuple, optional): x-limits for the x-axis. Default is [0, 50].
#     left_ylim (tuple, optional): y-limits for the left y-axis. Default is None.
#     right_ylim (tuple, optional): y-limits for the right y-axis. Default is None.
#     title (str, optional): Title of the plot. Default is None.
#     xlabel (str, optional): Label for the x-axis. Default is 'Time [s]'.
#     left_ylabel (str, optional): Label for the left y-axis. Default is None.
#     right_ylabel (str, optional): Label for the right y-axis. Default is None.
#     legend_loc_relative_to_anchor (str, optional): Location of the legend. Default is 'upper right'.
#     legend_anchor_loc (tuple, optional): Anchor location for the legend. Default is (0.9, 0.9).
#     save_path (str, optional): Path to save the plot image. Default is None.
#     left_yaxis_color (str, optional): Color for the left y-axis labels and ticks. Default is 'k'.
#     right_yaxis_color (str, optional): Color for the right y-axis labels and ticks. Default is 'k'.
#     """
    
#     # Create the plot
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot each variable on the left y-axis
#     for var, label, color, linestyle in zip(left_variables, left_labels, left_colors, left_linestyles):
#         ax1.plot(t, var, label=label, color=color, linestyle=linestyle, linewidth=1.5)
    
#     # Set axis labels and limits for the left y-axis
#     ax1.set_xlabel(xlabel)
#     ax1.set_xlim(xlim)
#     if left_ylabel:
#         ax1.set_ylabel(left_ylabel, color=left_yaxis_color)
#         ax1.set_ylim(left_ylim)
#     ax1.tick_params(axis='y', labelcolor=left_yaxis_color)
#     ax1.grid(True)
    
#     # Create a second y-axis for the right side
#     ax2 = ax1.twinx()
    
#     # Plot each variable on the right y-axis
#     for var, label, color, linestyle in zip(right_variables, right_labels, right_colors, right_linestyles):
#         ax2.plot(t, var, label=label, color=color, linestyle=linestyle, linewidth=1.5)
    
#     # Set axis labels and limits for the right y-axis
#     if right_ylabel:
#         ax2.set_ylabel(right_ylabel, color=right_yaxis_color)
#         ax2.set_ylim(right_ylim)
#     ax2.tick_params(axis='y', labelcolor=right_yaxis_color)
    
#     # Add title if provided
#     if title:
#         plt.title(title)
    
#     # Add legend
#     fig.legend(loc=legend_loc_relative_to_anchor, bbox_to_anchor=legend_anchor_loc)
    
#     # Adjust layout
#     fig.tight_layout()
    
#     # Save the figure if a save path is provided
#     if save_path:
#         plt.savefig(xplane11_path + 'richard_checked_terms\\' + save_path)
    
#     # Show the plot
#     plt.show()

# # User input
# whether_show_plot_variables_vs_variable_2_y_axis_example = 0

# if whether_show_plot_variables_vs_variable_2_y_axis_example == 1:

#     # Example test:
#     Framerate = 1/global_dt

#     left_variables = [u, v, w, Vt]
#     right_variables = [Framerate]
#     left_labels = [r'$u$ vs. $t$', r'$v$ vs. $t$', r'$w$ vs. $t$', 'TAS vs. $t$']
#     right_labels = ['Frame Rate vs. $t$']
#     left_colors = ['b', 'r', 'g', 'k']
#     right_colors = ['orange']
#     left_linestyles = ['-', '-', '-', ':']
#     right_linestyles = ['-']

#     plot_variables_vs_variable_2_y_axis(t, left_variables, right_variables, left_labels, right_labels, left_colors, right_colors, left_linestyles, right_linestyles,
#                             xplane11_path, xlim=[0, 25], left_ylim=[-20, 50], right_ylim=[0, 50],
#                             title=r'$u$, $v$, $w$ and TAS vs. $t$', xlabel='Time [s]', left_ylabel='Airspeeds [kn]', right_ylabel='Frame Rate [fps]',
#                             legend_loc_relative_to_anchor='upper right', legend_anchor_loc=(0.3, 0.9),
#                             left_yaxis_color=left_colors[0], right_yaxis_color=right_colors[0],
#                             save_path='custom_plot_with_two_y_axes.png')

# # Now we can plot any variables against each other with two available y-scales
# ################################################################################################################################################################################
# ################################################################################################################################################################################


# if whether_to_plot_X_X_dot == 1:
#     # Plotting u, v, w vs t
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot u on the left y-axis
#     ax1.plot(t, u, label=r'$u$ vs. $t$', color=color_u, linestyle='-', linewidth=1.5)
#     ax1.set_xlabel(r'Time $t$ [s]')
#     ax1.set_ylabel(r'Airspeeds [kn]', color='k')
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)

#     # Plot v on the left y-axis
#     ax1.plot(t, v, label=r'$v$ vs. $t$', color=color_v, linestyle='-', linewidth=1.5)

#     # Plot w on the left y-axis
#     ax1.plot(t, w, label=r'$w$ vs. $t$', color=color_w, linestyle='-', linewidth=1.5)

#     # Plot Vt on the left y-axis
#     ax1.plot(t, Vt, label=r'TAS vs. $t$', color='k', linestyle=':', linewidth=1.5)

#     # Title and legends
#     plt.title(r'$u$, $v$, $w$ and TAS vs. $t$')
#     fig.tight_layout()
#     fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\u_v_w_Vt_vs_t_plot.png')

#     # Show the plot
#     plt.show()


#     # Plotting p, q, r vs t
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot p on the left y-axis
#     ax1.plot(t, p, label=r'$p$ vs. $t$', color=color_p, linestyle='-', linewidth=1.5)
#     ax1.set_xlabel(r'Time $t$ [s]')
#     ax1.set_ylabel(r'Angular Rates [$^{\circ}$/s]', color='k')
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)

#     # Plot q on the left y-axis
#     ax1.plot(t, q, label=r'$q$ vs. $t$', color=color_q, linestyle='-', linewidth=1.5)

#     # Plot r on the left y-axis
#     ax1.plot(t, r, label=r'$r$ vs. $t$', color=color_r, linestyle='-', linewidth=1.5)

#     # Title and legends
#     plt.title(r'$p$, $q$, $r$ vs. $t$')
#     fig.tight_layout()
#     fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\p_q_r_vs_t_plot.png')

#     # Show the plot
#     plt.show()


#     # Plotting phi, theta, psi vs t
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot phi on the left y-axis
#     ax1.plot(t, phi, label=r'$\phi$ vs. $t$ (Roll)', color=color_phi, linestyle='-', linewidth=1.5)
#     ax1.set_xlabel(r'Time $t$ [s]')
#     ax1.set_ylabel(r'Angles [$^{\circ}$]', color='k')
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)

#     # Plot theta on the left y-axis
#     ax1.plot(t, theta, label=r'$\theta$ vs. $t$ (Pitch)', color=color_theta, linestyle='-', linewidth=1.5)

#     # Plot psi on the left y-axis
#     ax1.plot(t, psi, label=r'$\psi$ vs. $t$ (Yaw)', color=color_psi, linestyle='-', linewidth=1.5)

#     # Title and legends
#     plt.title(r'$\phi$, $\theta$, $\psi$ vs. $t$')
#     fig.tight_layout()
#     fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\phi_theta_psi_vs_t_plot.png')

#     # Show the plot
#     plt.show()


#     # Plotting X_e, Y_e, Z_e vs t
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot X_e on the left y-axis
#     ax1.plot(t, X_e, label=r'$X_e$ vs. $t$', color=color_X_e, linestyle='-', linewidth=1.5)
#     ax1.set_xlabel(r'Time $t$ [s]')
#     ax1.set_ylabel(r'Position [m]', color='k')
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)

#     # Plot Y_e on the left y-axis
#     ax1.plot(t, Y_e, label=r'$Y_e$ vs. $t$', color=color_Y_e, linestyle='-', linewidth=1.5)

#     # Plot Z_e on the left y-axis
#     ax1.plot(t, Z_e, label=r'$Z_e$ vs. $t$', color=color_Z_e, linestyle='-', linewidth=1.5)

#     # Title and legends
#     plt.title(r'$X_e$, $Y_e$, $Z_e$ vs. $t$')
#     fig.tight_layout()
#     fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\X_e_Y_e_Z_e_vs_t_plot.png')

#     # Show the plot
#     plt.show()


#     # Plotting u_dot, v_dot, w_dot vs t
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot u_dot on the left y-axis
#     ax1.plot(t, u_dot, label=r'$\dot{u}$ vs. $t$', color=color_u_dot, linestyle='-', linewidth=1.5)
#     ax1.set_xlabel(r'Time $t$ [s]')
#     ax1.set_ylabel(r'Accelerations [m/s²]', color='k')
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)

#     # Plot v_dot on the left y-axis
#     ax1.plot(t, v_dot, label=r'$\dot{v}$ vs. $t$', color=color_v_dot, linestyle='-', linewidth=1.5)

#     # Plot w_dot on the left y-axis
#     ax1.plot(t, w_dot, label=r'$\dot{w}$ vs. $t$', color=color_w_dot, linestyle='-', linewidth=1.5)

#     # Plot Vt on the left y-axis
#     ax1.plot(t, Vt_dot, label=r'$\dot{Vt}$ vs. $t$', color='k', linestyle=':', linewidth=1.5)
#     # plt.xlim([40,50])
#     # Title and legends
#     plt.title(r'$\dot{u}$, $\dot{v}$, $\dot{w}$, and $\dot{Vt}$ vs. $t$')
#     fig.tight_layout()
#     fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\u_dot_v_dot_w_dot_Vt_dot_vs_t_plot.png')

#     # Show the plot
#     plt.show()


#     # Plotting smoothed p_dot, q_dot, r_dot vs t
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot p_dot on the left y-axis
#     ax1.plot(t, p_dot, label=r'$\dot{p}$ vs. $t$', color=color_p_dot, linestyle='-', linewidth=1.5)
#     ax1.set_xlabel(r'Time $t$ [s]')
#     ax1.set_ylabel(r'Angular Accelerations [rad/s²]', color='k')
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)

#     # Plot q_dot on the left y-axis
#     ax1.plot(t, q_dot, label=r'$\dot{q}$ vs. $t$', color=color_q_dot, linestyle='-', linewidth=1.5)

#     # Plot r_dot on the left y-axis
#     ax1.plot(t, r_dot, label=r'$\dot{r}$ vs. $t$', color=color_r_dot, linestyle='-', linewidth=1.5)

#     # Title and legends
#     plt.title(r'$\dot{p}$, $\dot{q}$, $\dot{r}$ vs. $t$')
#     fig.tight_layout()
#     fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\smoothed_p_dot_q_dot_r_dot_vs_t_plot.png')

#     # Show the plot
#     plt.show()


#     # Plotting smoothed phi_dot, theta_dot, psi_dot vs t
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot phi_dot on the left y-axis
#     ax1.plot(t, phi_dot, label=r'$\dot{\phi}$ vs. $t$', color=color_phi_dot, linestyle='-', linewidth=1.5)
#     ax1.set_xlabel(r'Time $t$ [s]')
#     ax1.set_ylabel(r'Angular Rates [rad/s]', color='k')
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)

#     # Plot theta_dot on the left y-axis
#     ax1.plot(t, theta_dot, label=r'$\dot{\theta}$ vs. $t$', color=color_theta_dot, linestyle='-', linewidth=1.5)

#     # Plot psi_dot on the left y-axis
#     ax1.plot(t, psi_dot, label=r'$\dot{\psi}$ vs. $t$', color=color_psi_dot, linestyle='-', linewidth=1.5)

#     # Title and legends
#     plt.title(r'$\dot{\phi}$, $\dot{\theta}$, $\dot{\psi}$ vs. $t$')
#     fig.tight_layout()
#     fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\smoothed_phi_dot_theta_dot_psi_dot_vs_t_plot.png')

#     # Show the plot
#     plt.show()


#     # Plotting smoothed X_e_dot, Y_e_dot, Z_e_dot vs t
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot X_e_dot on the left y-axis
#     ax1.plot(t, X_e_dot, label=r'$\dot{X}_e$ vs. $t$', color=color_X_e_dot, linestyle='-', linewidth=1.5)
#     ax1.set_xlabel(r'Time $t$ [s]')
#     ax1.set_ylabel(r'Position Rates [m/s]', color='k')
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)

#     # Plot Y_e_dot on the left y-axis
#     ax1.plot(t, Y_e_dot, label=r'$\dot{Y}_e$ vs. $t$', color=color_Y_e_dot, linestyle='-', linewidth=1.5)

#     # Plot Z_e_dot on the left y-axis
#     ax1.plot(t, Z_e_dot, label=r'$\dot{Z}_e$ vs. $t$', color=color_Z_e_dot, linestyle='-', linewidth=1.5)

#     # Title and legends
#     plt.title(r'$\dot{X}_e$, $\dot{Y}_e$, $\dot{Z}_e$ vs. $t$')
#     fig.tight_layout()
#     fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\smoothed_X_e_dot_Y_e_dot_Z_e_dot_vs_t_plot.png')

#     # Show the plot
#     plt.show()


#     # Plotting smoothed X_e_dot, Y_e_dot, Z_e_dot vs t
#     fig, ax1 = plt.subplots(figsize=(10, 6))

#     # Plot X_e_dot on the left y-axis
#     ax1.plot(t, Z_e, label=r'$Z_e$ vs. $t$', color=color_X_e_dot, linestyle='-', linewidth=1.5)
#     ax1.set_xlabel(r'Time $t$ [s]')
#     ax1.set_ylabel(r'Altitude [m/s]', color='k')
#     ax1.tick_params(axis='y', labelcolor='k')
#     ax1.grid(True)

#     # Plot Y_e_dot on the left y-axis
#     ax1.plot(t, Z_e_raw, label=r'$Z_e$ Raw Data vs. $t$', color='k', linestyle='--', linewidth=1.5)

#     # Title and legends
#     plt.title(r'$Z_e$ vs. $t$')
#     fig.tight_layout()
#     fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\Ze_Ze_raw_vs_t.png')

#     # Show the plot
#     plt.show()


# ################################################################################################################################################################################
# ################################################################################################################################################################################
# if whether_to_plot_3D_flight_path == 1:

#     # Plotting the 3D flight path
#     fig = plt.figure(figsize=(10, 7))
#     ax = fig.add_subplot(111, projection='3d')

#     # Plot the flight path
#     ax.plot(X_e, Y_e, Z_e, label='Flight Path', color='b', linewidth=2)

#     # Set labels and title
#     ax.set_xlabel(r'$X_e$ [m]')
#     ax.set_ylabel(r'$Y_e$ [m]')
#     ax.set_zlabel(r'$Z_e$ [m]')
#     ax.set_title('3D Flight Path')

#     whether_to_plot_3D_equal_x_y_z_scale = 0

#     if whether_to_plot_3D_equal_x_y_z_scale == 1:

#         # Ensure the aspect ratio is equal by setting the limits for each axis
#         max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0

#         mid_x = (X_e.max() + X_e.min()) * 0.5
#         mid_y = (Y_e.max() + Y_e.min()) * 0.5
#         mid_z = (Z_e.max() + Z_e.min()) * 0.5

#         ax.set_xlim(mid_x - max_range, mid_x + max_range)
#         ax.set_ylim(mid_y - max_range, mid_y + max_range)
#         ax.set_zlim(mid_z - max_range, mid_z + max_range)

#     # # Invert the Z-axis for better visualization (if necessary)
#     # ax.invert_zaxis()

#     # Show grid
#     ax.grid(True)

#     # Show legend
#     ax.legend()

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\3D_flight_path.png')

#     # Show the plot
#     plt.show()


# ################################################################################################################################################################################
# ################################################################################################################################################################################
# if whether_to_create_3D_flight_path_animation == 1:

#     # # Animation v5
#     # Interpolate the data for smoother animation
#     frames = 1000  # Increase the number of frames for smoother animation
#     interp_func_x = interp1d(t, X_e, kind='linear')
#     interp_func_y = interp1d(t, Y_e, kind='linear')
#     interp_func_z = interp1d(t, Z_e, kind='linear')
#     t_interp = np.linspace(t[0], t[-1], frames)
#     X_e_interp = interp_func_x(t_interp)
#     Y_e_interp = interp_func_y(t_interp)
#     Z_e_interp = interp_func_z(t_interp)

#     # Calculate the interval between frames to match real-time duration
#     total_duration = t[-1] - t[0]  # Total flight time in seconds
#     interval = total_duration / frames * 1000  # Interval in milliseconds

#     # Create a figure for the animation
#     fig = plt.figure(figsize=(10, 7))
#     ax = fig.add_subplot(111, projection='3d')

#     # Set up the plot limits to be equal
#     max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
#     mid_x = (X_e.max() + X_e.min()) * 0.5
#     mid_y = (Y_e.max() + Y_e.min()) * 0.5
#     mid_z = (Z_e.max() + Z_e.min()) * 0.5

#     ax.set_xlim(mid_x - max_range, mid_x + max_range)
#     ax.set_ylim(mid_y - max_range, mid_y + max_range)
#     ax.set_zlim(mid_z - max_range, mid_z + max_range)

#     # # Invert the Z-axis for better visualization (if necessary)
#     # ax.invert_zaxis()

#     # Initialize a line object and a point object for the aircraft position
#     line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)
#     point, = ax.plot([], [], [], 'ro', label='Aircraft Position')  # Red dot for the aircraft

#     # Function to initialize the background of the animation
#     def init():
#         line.set_data([], [])
#         line.set_3d_properties([])
#         point.set_data([], [])
#         point.set_3d_properties([])
#         return line, point

#     # Function to update the animation at each frame
#     def update(num, X_e_interp, Y_e_interp, Z_e_interp, line, point):
#         line.set_data(X_e_interp[:num], Y_e_interp[:num])
#         line.set_3d_properties(Z_e_interp[:num])
#         point.set_data(X_e_interp[num], Y_e_interp[num])
#         point.set_3d_properties(Z_e_interp[num])
#         return line, point

#     # Create the animation with the correct interval
#     ani = animation.FuncAnimation(fig, update, frames=frames, fargs=(X_e_interp, Y_e_interp, Z_e_interp, line, point),
#                                 init_func=init, blit=True, interval=interval)

#     # Add labels and title
#     ax.set_xlabel(r'$X_e$ [m]')
#     ax.set_ylabel(r'$Y_e$ [m]')
#     ax.set_zlabel(r'$Z_e$ [m]')
#     ax.set_title('3D Flight Path Animation with Smooth Motion')
#     ax.legend()

#     # Save the animation as an mp4 file
#     ani.save(xplane11_path + 'richard_checked_terms\\3D_flight_path_animation_smooth.mp4', writer='ffmpeg')

#     # Display the animation
#     plt.show()




#     # # Animation v1
#     # # Create a figure for the animation
#     # fig = plt.figure(figsize=(10, 7))
#     # ax = fig.add_subplot(111, projection='3d')

#     # # Set up the plot limits to be equal
#     # max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
#     # mid_x = (X_e.max() + X_e.min()) * 0.5
#     # mid_y = (Y_e.max() + Y_e.min()) * 0.5
#     # mid_z = (Z_e.max() + Z_e.min()) * 0.5

#     # ax.set_xlim(mid_x - max_range, mid_x + max_range)
#     # ax.set_ylim(mid_y - max_range, mid_y + max_range)
#     # ax.set_zlim(mid_z - max_range, mid_z + max_range)

#     # # # Invert the Z-axis for better visualization (if necessary)
#     # # ax.invert_zaxis()

#     # # Initialize a line object which will be updated during the animation
#     # line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)

#     # # Function to initialize the background of the animation
#     # def init():
#     #     line.set_data([], [])
#     #     line.set_3d_properties([])
#     #     return line,

#     # # Function to update the animation at each frame
#     # def update(num, X_e, Y_e, Z_e, line):
#     #     line.set_data(X_e[:num], Y_e[:num])
#     #     line.set_3d_properties(Z_e[:num])
#     #     return line,

#     # # Create the animation
#     # ani = animation.FuncAnimation(fig, update, frames=len(t), fargs=(X_e, Y_e, Z_e, line),
#     #                               init_func=init, blit=True, interval=30)

#     # # Add labels and title
#     # ax.set_xlabel(r'$X_e$ [m]')
#     # ax.set_ylabel(r'$Y_e$ [m]')
#     # ax.set_zlabel(r'$Z_e$ [m]')
#     # ax.set_title('3D Flight Path Animation')
#     # ax.legend()

#     # # Save the animation as an mp4 file
#     # ani.save(xplane11_path + 'richard_checked_terms\\3D_flight_path_animation.mp4', writer='ffmpeg')

#     # # Display the animation
#     # plt.show()







#     # # Animation v2
#     # # Calculate the time intervals between frames in milliseconds
#     # time_intervals = np.diff(t) * 1000  # Convert from seconds to milliseconds

#     # # Create a figure for the animation
#     # fig = plt.figure(figsize=(10, 7))
#     # ax = fig.add_subplot(111, projection='3d')

#     # # Set up the plot limits to be equal
#     # max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
#     # mid_x = (X_e.max() + X_e.min()) * 0.5
#     # mid_y = (Y_e.max() + Y_e.min()) * 0.5
#     # mid_z = (Z_e.max() + Z_e.min()) * 0.5

#     # ax.set_xlim(mid_x - max_range, mid_x + max_range)
#     # ax.set_ylim(mid_y - max_range, mid_y + max_range)
#     # ax.set_zlim(mid_z - max_range, mid_z + max_range)

#     # # # Invert the Z-axis for better visualization (if necessary)
#     # # ax.invert_zaxis()

#     # # Initialize a line object which will be updated during the animation
#     # line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)

#     # # Function to initialize the background of the animation
#     # def init():
#     #     line.set_data([], [])
#     #     line.set_3d_properties([])
#     #     return line,

#     # # Function to update the animation at each frame
#     # def update(num, X_e, Y_e, Z_e, line):
#     #     line.set_data(X_e[:num], Y_e[:num])
#     #     line.set_3d_properties(Z_e[:num])
#     #     return line,

#     # # Create the animation with real-time progression
#     # ani = animation.FuncAnimation(fig, update, frames=len(t), fargs=(X_e, Y_e, Z_e, line),
#     #                               init_func=init, blit=True, interval=time_intervals[0])

#     # # Add labels and title
#     # ax.set_xlabel(r'$X_e$ [m]')
#     # ax.set_ylabel(r'$Y_e$ [m]')
#     # ax.set_zlabel(r'$Z_e$ [m]')
#     # ax.set_title('3D Flight Path Animation')
#     # ax.legend()

#     # # Save the animation as an mp4 file
#     # ani.save(xplane11_path + 'richard_checked_terms\\3D_flight_path_animation_real_time.mp4', writer='ffmpeg')

#     # # Display the animation
#     # plt.show()






#     # # Animation v3
#     # # Calculate the total duration and number of frames
#     # total_duration = t[-1] - t[0]  # Total flight time in seconds (50 seconds in this case)
#     # num_frames = len(t)  # Total number of frames

#     # # Calculate the interval between frames to match real-time duration
#     # interval = total_duration / num_frames * 1000  # Interval in milliseconds

#     # # Create a figure for the animation
#     # fig = plt.figure(figsize=(10, 7))
#     # ax = fig.add_subplot(111, projection='3d')

#     # # Set up the plot limits to be equal
#     # max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
#     # mid_x = (X_e.max() + X_e.min()) * 0.5
#     # mid_y = (Y_e.max() + Y_e.min()) * 0.5
#     # mid_z = (Z_e.max() + Z_e.min()) * 0.5

#     # ax.set_xlim(mid_x - max_range, mid_x + max_range)
#     # ax.set_ylim(mid_y - max_range, mid_y + max_range)
#     # ax.set_zlim(mid_z - max_range, mid_z + max_range)

#     # # # Invert the Z-axis for better visualization (if necessary)
#     # # ax.invert_zaxis()

#     # # Initialize a line object which will be updated during the animation
#     # line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)

#     # # Function to initialize the background of the animation
#     # def init():
#     #     line.set_data([], [])
#     #     line.set_3d_properties([])
#     #     return line,

#     # # Function to update the animation at each frame
#     # def update(num, X_e, Y_e, Z_e, line):
#     #     line.set_data(X_e[:num], Y_e[:num])
#     #     line.set_3d_properties(Z_e[:num])
#     #     return line,

#     # # Create the animation with the correct interval
#     # ani = animation.FuncAnimation(fig, update, frames=num_frames, fargs=(X_e, Y_e, Z_e, line),
#     #                               init_func=init, blit=True, interval=interval)

#     # # Add labels and title
#     # ax.set_xlabel(r'$X_e$ [m]')
#     # ax.set_ylabel(r'$Y_e$ [m]')
#     # ax.set_zlabel(r'$Z_e$ [m]')
#     # ax.set_title('3D Flight Path Animation')
#     # ax.legend()

#     # # Save the animation as an mp4 file
#     # ani.save(xplane11_path + 'richard_checked_terms\\3D_flight_path_animation_real_time_fixed.mp4', writer='ffmpeg')

#     # # Display the animation
#     # plt.show()






#     # # Animation v4
#     # # Calculate the total duration and number of frames
#     # total_duration = t[-1] - t[0]  # Total flight time in seconds (50 seconds in this case)
#     # num_frames = len(t)  # Total number of frames

#     # # Calculate the interval between frames to match real-time duration
#     # interval = total_duration / num_frames * 1000  # Interval in milliseconds

#     # # Create a figure for the animation
#     # fig = plt.figure(figsize=(10, 7))
#     # ax = fig.add_subplot(111, projection='3d')

#     # # Set up the plot limits to be equal
#     # max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
#     # mid_x = (X_e.max() + X_e.min()) * 0.5
#     # mid_y = (Y_e.max() + Y_e.min()) * 0.5
#     # mid_z = (Z_e.max() + Z_e.min()) * 0.5

#     # ax.set_xlim(mid_x - max_range, mid_x + max_range)
#     # ax.set_ylim(mid_y - max_range, mid_y + max_range)
#     # ax.set_zlim(mid_z - max_range, mid_z + max_range)

#     # # Invert the Z-axis for better visualization (if necessary)
#     # ax.invert_zaxis()

#     # # Initialize a line object and a point object for the aircraft position
#     # line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)
#     # point, = ax.plot([], [], [], 'ro', label='Aircraft Position')  # Red dot for the aircraft

#     # # Function to initialize the background of the animation
#     # def init():
#     #     line.set_data([], [])
#     #     line.set_3d_properties([])
#     #     point.set_data([], [])
#     #     point.set_3d_properties([])
#     #     return line, point

#     # # Function to update the animation at each frame
#     # def update(num, X_e, Y_e, Z_e, line, point):
#     #     line.set_data(X_e[:num], Y_e[:num])
#     #     line.set_3d_properties(Z_e[:num])
#     #     point.set_data(X_e[num], Y_e[num])
#     #     point.set_3d_properties(Z_e[num])
#     #     return line, point

#     # # Create the animation with the correct interval
#     # ani = animation.FuncAnimation(fig, update, frames=num_frames, fargs=(X_e, Y_e, Z_e, line, point),
#     #                               init_func=init, blit=True, interval=interval)

#     # # Add labels and title
#     # ax.set_xlabel(r'$X_e$ [m]')
#     # ax.set_ylabel(r'$Y_e$ [m]')
#     # ax.set_zlabel(r'$Z_e$ [m]')
#     # ax.set_title('3D Flight Path Animation')
#     # ax.legend()

#     # # Save the animation as an mp4 file
#     # ani.save(xplane11_path + 'richard_checked_terms\\3D_flight_path_animation_real_time_with_dot.mp4', writer='ffmpeg')

#     # # Display the animation
#     # plt.show()


# if whether_to_plot_U == 1:

#     # Create a figure and axis objects with 5 subplots stacked vertically
#     fig, axs = plt.subplots(5, 1, figsize=(10, 10), sharex=True)

#     U_plots_axes_0_upper_buffer = (0.007/1.00) * (1-0)
#     U_plots_axes_0_lower_buffer = (0.001/1.00) * (1-0)

#     U_plots_axes_1_upper_buffer = (0.007/1.00) * (9-(-18.4))
#     U_plots_axes_1_lower_buffer = (0.001/1.00) * (9-(-18.4))

#     U_plots_axes_2_upper_buffer = (0.007/1.00) * (24-(-24))
#     U_plots_axes_2_lower_buffer = (0.001/1.00) * (24-(-24))

#     U_plots_axes_3_upper_buffer = (0.007/1.00) * (23.21-(-23.21))
#     U_plots_axes_3_lower_buffer = (0.001/1.00) * (23.21-(-23.21))

#     U_plots_axes_4_upper_buffer = (0.0065/1.00) * (-0-(-40))
#     U_plots_axes_4_lower_buffer = (0.001/1.00) * (-0-(-40))



#     # Plot Throttle Input
#     axs[0].plot(t, delta_t, color='r')
#     axs[0].set_ylabel(r'Throttle [-]')
#     axs[0].set_ylim([0 - U_plots_axes_0_lower_buffer, 1 + U_plots_axes_0_upper_buffer])
#     axs[0].grid(True)
#     axs[0].set_title('Control Inputs vs. Time')

#     # Plot Elevator Input
#     axs[1].plot(t, delta_e, color='g')
#     axs[1].set_ylabel(r'Elevator [$^{\circ}$]')
#     axs[1].set_ylim([-18.4 - U_plots_axes_1_lower_buffer, 9 + U_plots_axes_1_upper_buffer])
#     axs[1].grid(True)

#     # Plot Aileron Input
#     axs[2].plot(t, delta_a, color='b')
#     axs[2].set_ylabel(r'Aileron [$^{\circ}$]')
#     axs[2].set_ylim([-24 - U_plots_axes_2_lower_buffer, 24 + U_plots_axes_2_upper_buffer])
#     axs[2].grid(True)

#     # Plot Rudder Input
#     axs[3].plot(t, delta_r, color='m')
#     axs[3].set_ylabel(r'Rudder [$^{\circ}$]')
#     axs[3].set_ylim([-23.21 - U_plots_axes_3_lower_buffer, 23.21 + U_plots_axes_3_upper_buffer])
#     axs[3].grid(True)

#     # Plot Flap Input
#     axs[4].plot(t, delta_f, color='c')
#     axs[4].set_ylabel(r'Flap [$^{\circ}$]')
#     axs[4].set_ylim([-40 - U_plots_axes_4_lower_buffer, -0 + U_plots_axes_4_upper_buffer])
#     axs[4].set_xlabel('Time [s]')
#     axs[4].grid(True)

#     # Adjust layout to prevent overlap
#     plt.tight_layout()

#     # Save the figure
#     plt.savefig(xplane11_path + 'richard_checked_terms\\control_inputs_vs_time.png')

#     # Show the plot
#     plt.show()

#     # print(min(delta_e))








# # # Plotting the data
# # fig, ax1 = plt.subplots(figsize=(10, 6))

# # # Plot X_flow_sep on the left y-axis
# # ax1.plot(t, X_flow_sep, label=r'$X_{flow\_sep}$ vs. Time', color=color2)
# # ax1.set_xlabel('Time [s]')
# # # ax1.set_xlim([15, 30])
# # ax1.set_ylabel(r'$X_{flow\_sep}$ [-]', color=color2)
# # # if np.all((X_flow_sep >= 0) & (X_flow_sep <= 1)):
# # #     ax1.set_ylim([-0.1, 1.1])
# # ax1.set_ylim([-0.1, 1.1])
# # ax1.tick_params(axis='y', labelcolor=color2)
# # ax1.grid(True)

# # # Create a second y-axis for alpha, CL, poststall_on, Cm, and CD
# # ax2 = ax1.twinx()
# # ax2.plot(t, alpha, label=r'$\alpha$ vs. Time', color=color1, linestyle='--')
# # ax2.plot(t, alpha_crit_1, label=r'$\alpha_{crit\_1}$ vs. Time', color=color3, linestyle='-.')  
# # ax2.plot(t, alpha_crit_2, label=r'$\alpha_{crit\_2}$ vs. Time', color=color4, linestyle='-.')  
# # ax2.plot(t, alpha_crit_3, label=r'$\alpha_{crit\_3}$ vs. Time', color=color5, linestyle='-.')  
# # ax2.plot(t, 2 * CL, label=r'2 $C_L$ vs. Time', color=color6, linestyle='-')  # CL is multiplied by 2 and displayed on ax2 scale
# # ax2.plot(t, poststall_on, label='Post-Stall On vs. Time', color=color7, linestyle='-')  # Pre-stall/Post-stall indicator
# # ax2.plot(t, 2 * Cm, label=r'2 $C_m$ vs. Time', color=color8, linestyle='-')  # Cm
# # ax2.plot(t, 5 * CD, label=r'5 $C_D$ vs. Time', color=color9, linestyle='-')  # CD
# # ax2.set_ylabel(r'$\alpha$ [$^{\circ}$]', color=color1)
# # ax2.tick_params(axis='y', labelcolor=color1)

# # # Title and legends
# # plt.title(r'$X_{flow\_sep}$ and $C_L$, $C_m$, $C_D$ vs. Time')
# # fig.tight_layout()
# # fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))
# # plt.savefig('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\X_flow_sep_CL_Cm_CD_vs_time.png')  # Save as PNG
# # plt.show()

# # print(X_flow_sep[0])




# # # Plotting the data
# # fig, ax1 = plt.subplots(figsize=(10, 6))

# # # Plot X_flow_sep on the left y-axis
# # ax1.plot(t, X_flow_sep, label=r'$X_{flow\_sep}$', color=color2)
# # ax1.plot(t, Vt/1.943844/50, label=r'TAS/$V_{stall}$', color='k', linestyle=':', linewidth=1.5)
# # ax1.set_xlabel('Time [s]')
# # ax1.set_xlim([15, 26])
# # ax1.set_ylabel(r'$X_{flow\_sep}$ [-], TAS/$V_{stall}$ [-]', color=color2)
# # ax1.set_ylim([-0.1, 1.1])
# # ax1.tick_params(axis='y', labelcolor=color2)
# # ax1.grid(True)

# # # Create a second y-axis for alpha, alpha_crit, and poststall_on
# # ax2 = ax1.twinx()
# # ax2.plot(t, alpha, label=r'$\alpha$', color=color1, linestyle='--')
# # ax2.plot(t, alpha_crit_1, label=r'$\alpha_{crit\_1}$', color=color3, linestyle='-.')
# # ax2.plot(t, alpha_crit_2, label=r'$\alpha_{crit\_2}$', color=color4, linestyle='-.')
# # ax2.plot(t, alpha_crit_3, label=r'$\alpha_{crit\_3}$', color=color5, linestyle='-.')
# # ax2.plot(t, poststall_on, label='Post-Stall On', color=color7, linestyle='-')
# # ax2.set_ylabel(r'$\alpha$ [$^{\circ}$]', color=color1)
# # ax2.set_ylim([-15, 45])
# # ax2.tick_params(axis='y', labelcolor=color1)

# # # Create a third y-axis for CL, Cm, CD
# # ax3 = ax1.twinx()
# # ax3.spines["right"].set_position(("outward", 60))  # Move the third y-axis outward
# # ax3.plot(t, CL, label=r'$C_L$', color=color6, linestyle='-')  # CL is multiplied by 2
# # ax3.plot(t, Cm, label=r'$C_m$', color=color8, linestyle='-')  # Cm is multiplied by 2
# # ax3.plot(t, CD, label=r'$C_D$', color=color9, linestyle='-')  # CD is multiplied by 5
# # ax3.set_ylabel(r'$C_L$, $C_m$, $C_D$', color=color6)
# # ax3.set_ylim([-15/10/1.5, 45/10/1.5])
# # ax3.tick_params(axis='y', labelcolor=color6)

# # # Adding legends
# # ax1.legend(loc='upper left')
# # ax2.legend(loc='upper center')
# # ax3.legend(loc='upper right')

# # # Title and layout
# # plt.title(r'$X_{flow\_sep}$, TAS/$V_{stall}$, $\alpha$, $C_L$, $C_m$, $C_D$ vs. Time')
# # fig.tight_layout()

# # # Save and display
# # plt.savefig('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\X_flow_sep_alpha_CL_Cm_CD_vs_time.png')  # Save as PNG
# # plt.show()