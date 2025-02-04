import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter1d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.interpolate import interp1d
import os

from M13_precision_timing_v2 import calculate_and_print_effects_v2


'''Each Row in the Flight Data Matrix is a data entry at a specific moment in time'''


# Please read the following options and choose what output is wanted
#################################################################################################
### User preferences:

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
whether_to_plot_U = 1
# input field
whether_to_plot_stage2_combined_history = 1
# input field
whether_to_plot_stage3_combined_history = 0
# input field
whether_to_plot_precision_timing = 0
# input field
whether_to_plot_Xg = 0
# input field
whether_to_plot_stage4_combined_history = 0

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

X = flight_data_matrix[:, 1:13]
# print(X[-10:-1, :]) # print the last 10 rows
# print(X[0:100, :]) # print the first 100 rows

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
# print(U_data_matrix[0:100, :]) # print the first 100 rows

delta_t = flight_data_matrix[:, 25] # 25th column (delta_t) (Throttle input) [-]
delta_e = flight_data_matrix[:, 26] # 25th column (delta_t) (elevator input) [rad]
delta_a = flight_data_matrix[:, 27] # 25th column (delta_t) (aileron input) [rad]
delta_r = flight_data_matrix[:, 28] # 25th column (delta_t) (rudder input) [rad]
delta_f = flight_data_matrix[:, 29] # 25th column (delta_t) (flap input) [rad]

stage2_combined_history_data_matrix = flight_data_matrix[:, 30:44]
# print(stage2_combined_history_data_matrix[-10:-1, :]) # print the last 10 rows

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

# Variables used in stage2_combined_history plot 2
# X_flow_sep = flight_data_matrix[:, 33]
# CL = flight_data_matrix[:, 40]
# poststall_on = flight_data_matrix[:, 41]
# Cm = flight_data_matrix[:, 42]
# CD = flight_data_matrix[:, 43]

stage3_combined_history_data_matrix = flight_data_matrix[:, 44:61]

dt = flight_data_matrix[:, 44]
rho = flight_data_matrix[:, 45]
a_sound_speed = flight_data_matrix[:, 46]
beta = flight_data_matrix[:, 47]
Cx = flight_data_matrix[:, 48]
Cy = flight_data_matrix[:, 49]
Cz = flight_data_matrix[:, 50]
Cl = flight_data_matrix[:, 51]
Cn = flight_data_matrix[:, 52]
P_prop_power = flight_data_matrix[:, 53]
n_prop_rot_speed = flight_data_matrix[:, 54]
global_n_prop_rot_speed = flight_data_matrix[:, 55]
n_prop_rot_speed_dot = flight_data_matrix[:, 56]
J_advance_ratio = flight_data_matrix[:, 57]
eta_p_prop_eff = flight_data_matrix[:, 58]
ay_noise = flight_data_matrix[:, 59]
az_noise = flight_data_matrix[:, 60]

Xg = flight_data_matrix[:, 61:73]

u_g = flight_data_matrix[:, 61]
v_g = flight_data_matrix[:, 62]
w_g = flight_data_matrix[:, 63]
p_g = flight_data_matrix[:, 64]
q_g = flight_data_matrix[:, 65]
r_g = flight_data_matrix[:, 66]
phi_g = flight_data_matrix[:, 67]  # 67th column (phi_g)
theta_g = flight_data_matrix[:, 68]  # 68th column (theta_g)
psi_g = flight_data_matrix[:, 69]  # 69th column (psi_g)
X_e_g = flight_data_matrix[:, 70] # 70th column (X_e_g)
Y_e_g = flight_data_matrix[:, 71] # 71th column (Y_e_g)
Z_e_g = flight_data_matrix[:, 72] # 72th column (Z_e_g)

stage4_combined_history = flight_data_matrix[:, 73:81]

alpha_average_of_left_wing = flight_data_matrix[:, 73]
alpha_average_of_right_wing = flight_data_matrix[:, 74]
CL_average_of_left_wing = flight_data_matrix[:, 75]
CL_average_of_right_wing = flight_data_matrix[:, 76]
CD_average_of_left_wing = flight_data_matrix[:, 77]
CD_average_of_right_wing = flight_data_matrix[:, 78]
Cl_increment_due_to_autorotation = flight_data_matrix[:, 79]
Cn_increment_due_to_autorotation = flight_data_matrix[:, 80]

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

#------------------------------------------------------------------------------------------------------
# Record in NED axes as well for precision timing module M13
X_e_in_NED = X_e
Y_e_in_NED = Y_e
Z_e_in_NED = Z_e_raw # use the Z_e from raw data instead of from integration of Z_e_dot
altitude_in_NED = -Z_e_in_NED # Simply the altitude as a positive number
#------------------------------------------------------------------------------------------------------

Z_e = Z_e + Z_e_raw_initial
Z_e_final = Z_e[-1]
Z_e = Z_e + (t/t[-1]) * (Z_e_raw_final - Z_e_final) # Need a mysterious correction, might be the deviation between altitude of default and custom physics since part of the journey is on default physics (t< approx 10s)

# Currently, X_e, Y_e, Z_e represents North-East-Down (NED) axes

# Invert Z_e to obtain a more intuitive alitude expression
Z_e_raw = - Z_e_raw
Z_e = -Z_e

# Invert Y_e accordingly to perserve right-hand coordinate system
Y_e = -Y_e

# Now X_e, Y_e, Z_e represents North-West-Up (NWU) axes (Just like how you read a map)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# from M13
theta_of_earth_spherical_coordinate_system = np.deg2rad(90)  # (deg) At equator (0 deg is North pole, 180 deg is South pole)
nanosecond_per_day_drift_for_atomic_clock = np.full(len(t), 10)  # Time series for atomic clock drift over time

(delta_f_to_f_clock, delta_f_to_f_doppler, delta_f_to_f_special_relativity,
delta_f_to_f_gravity, delta_f_to_f_acceleration, delta_T_to_T_frame_dragging,
delta_T_to_T_light_speed_reduction, delta_T_to_T_ionospheric, delta_T_to_T_tropospheric, delta_T_to_T_sagnac) = calculate_and_print_effects_v2(altitude_in_NED, Vt, theta, psi, X_e_in_NED, Y_e_in_NED, Z_e_in_NED, rho, u_dot, theta_of_earth_spherical_coordinate_system, nanosecond_per_day_drift_for_atomic_clock)

# Summing up all apparent effects  
total_apparent_delta_f_to_f = (delta_f_to_f_clock
                             + delta_f_to_f_doppler)

# Summing up all real effects
total_real_delta_f_to_f = (delta_f_to_f_special_relativity
                         + delta_f_to_f_gravity
                         + delta_f_to_f_acceleration)

# Summing up all delay effects
total_delta_T_to_T = (delta_T_to_T_frame_dragging
                    + delta_T_to_T_light_speed_reduction
                    + delta_T_to_T_ionospheric
                    + delta_T_to_T_tropospheric
                    + delta_T_to_T_sagnac)

#----------------------------------------------------------------------------------------------
# User input
# The range of time that you want data plotted
user_defined_t_initial = 0 # sec
user_defined_t_final = 5 # sec
# user_defined_t_initial = 0 # sec
# user_defined_t_final = 50 # sec

# Find the indices where values in t range from 15 to 20
t_index_range = np.where((t >= user_defined_t_initial) & (t <= user_defined_t_final))[0]
# If there are valid indices, find the first and last
if t_index_range.size > 0:
    first_t_index = t_index_range[0]
    last_t_index = t_index_range[-1]
    print(f"[{first_t_index}:{last_t_index}]")
#----------------------------------------------------------------------------------------------
# Change the initial global dt to something sensible (start of the sim always come with large global_dt)
global_dt[first_t_index:last_t_index] = 0.03

total_apparent_time_drift = cumulative_trapezoid(total_apparent_delta_f_to_f * global_dt, t, initial=0)  # Integrate total_apparent_delta_f_to_f to get total_apparent_time_drift
total_real_time_drift = cumulative_trapezoid(total_real_delta_f_to_f * global_dt, t, initial=0)  # Integrate total_real_delta_f_to_f * global_dt to get total_real_time_drift
total_signal_path_induced_time_drift = total_delta_T_to_T * (np.sqrt(X_e**2 + Y_e**2 + Z_e**2)/299792458)

combined_time_drift = total_apparent_time_drift + total_real_time_drift + total_signal_path_induced_time_drift
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Unit changes for plotting purpose
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

beta = (180/np.pi) * beta
n_prop_rot_speed = n_prop_rot_speed * 60
global_n_prop_rot_speed = global_n_prop_rot_speed * 60
n_prop_rot_speed_dot = n_prop_rot_speed_dot * 60

alpha_average_of_left_wing = (180/np.pi) * alpha_average_of_left_wing
alpha_average_of_right_wing = (180/np.pi) * alpha_average_of_right_wing

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

# Colors for plotting stage3_combined_history
color_dt = [0/255, 128/255, 255/255]  # Customized Blue for dt
color_rho = [0/255, 100/255, 0/255]  # Customized Green for rho
color_a = [0/255, 0/255, 255/255]    # Customized Blue for a_sound_speed
color_alpha = [255/255, 36/255, 10/255]  # Customized Red for alpha # 2nd appearance
color_beta = [0/255, 0/255, 255/255]     # Customized Blue for beta
color_Cx = [255/255, 36/255, 10/255]  # Customized Red for Cx
color_Cy = [61/255, 142/255, 190/255]  # Customized Blue for Cy
color_Cz = [21/255, 110/255, 72/255]   # Customized Green for Cz
color_Cl = [139/255, 68/255, 68/255]   # Contrast Red for Cl
color_Cm = [0/255, 209/255, 167/255]   # Contrast Green for Cm # 2nd appearance
color_Cn = [136/255, 148/255, 187/255] # Contrast Cyan for Cn

# Define custom colors for plotting
color_u_g = [61/255, 142/255, 190/255]  # Customized Blue for u_g
color_v_g = [255/255, 36/255, 10/255]   # Customized Red for v_g
color_w_g = [21/255, 110/255, 72/255]   # Customized Green for w_g
color_p_g = [255/255, 159/255, 64/255]  # Customized Orange for p_g
color_q_g = [75/255, 192/255, 192/255]  # Customized Cyan for q_g
color_r_g = [153/255, 102/255, 255/255]  # Customized Purple for r_g


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


################################################################################################################################################################################
################################################################################################################################################################################
if whether_to_plot_X_X_dot == 1:
    # Plotting u, v, w vs t
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot u on the left y-axis
    ax1.plot(t, u, label=r'$u$ vs. $t$', color=color_u, linestyle='-', linewidth=1.5)
    ax1.set_xlabel(r'Time $t$ [s]')
    ax1.set_ylabel(r'Airspeeds [kn]', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True)

    # Plot v on the left y-axis
    ax1.plot(t, v, label=r'$v$ vs. $t$', color=color_v, linestyle='-', linewidth=1.5)

    # Plot w on the left y-axis
    ax1.plot(t, w, label=r'$w$ vs. $t$', color=color_w, linestyle='-', linewidth=1.5)

    # Plot Vt on the left y-axis
    ax1.plot(t, Vt, label=r'TAS vs. $t$', color='k', linestyle=':', linewidth=1.5)

    # Title and legends
    plt.title(r'$u$, $v$, $w$ and TAS vs. $t$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A01_u_v_w_Vt_vs_t_plot.png')

    # Show the plot
    plt.show()


    # Plotting p, q, r vs t
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot p on the left y-axis
    ax1.plot(t, p, label=r'$p$ vs. $t$', color=color_p, linestyle='-', linewidth=1.5)
    ax1.set_xlabel(r'Time $t$ [s]')
    ax1.set_ylabel(r'Angular Rates [$^{\circ}$/s]', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True)

    # Plot q on the left y-axis
    ax1.plot(t, q, label=r'$q$ vs. $t$', color=color_q, linestyle='-', linewidth=1.5)

    # Plot r on the left y-axis
    ax1.plot(t, r, label=r'$r$ vs. $t$', color=color_r, linestyle='-', linewidth=1.5)

    # Title and legends
    plt.title(r'$p$, $q$, $r$ vs. $t$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A02_p_q_r_vs_t_plot.png')

    # Show the plot
    plt.show()


    # Plotting phi, theta, psi vs t
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot phi on the left y-axis
    ax1.plot(t, phi, label=r'$\phi$ vs. $t$ (Roll)', color=color_phi, linestyle='-', linewidth=1.5)
    ax1.set_xlabel(r'Time $t$ [s]')
    ax1.set_ylabel(r'Angles [$^{\circ}$]', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True)

    # Plot theta on the left y-axis
    ax1.plot(t, theta, label=r'$\theta$ vs. $t$ (Pitch)', color=color_theta, linestyle='-', linewidth=1.5)

    # Plot psi on the left y-axis
    ax1.plot(t, psi, label=r'$\psi$ vs. $t$ (Yaw)', color=color_psi, linestyle='-', linewidth=1.5)

    # Title and legends
    plt.title(r'$\phi$, $\theta$, $\psi$ vs. $t$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A03_phi_theta_psi_vs_t_plot.png')

    # Show the plot
    plt.show()


    # Plotting X_e, Y_e, Z_e vs t
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot X_e on the left y-axis
    ax1.plot(t, X_e, label=r'$X_e$ vs. $t$', color=color_X_e, linestyle='-', linewidth=1.5)
    ax1.set_xlabel(r'Time $t$ [s]')
    ax1.set_ylabel(r'Position [m]', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True)

    # Plot Y_e on the left y-axis
    ax1.plot(t, Y_e, label=r'$Y_e$ vs. $t$', color=color_Y_e, linestyle='-', linewidth=1.5)

    # Plot Z_e on the left y-axis
    ax1.plot(t, Z_e, label=r'$Z_e$ vs. $t$', color=color_Z_e, linestyle='-', linewidth=1.5)

    # Title and legends
    plt.title(r'$X_e$, $Y_e$, $Z_e$ vs. $t$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A04_X_e_Y_e_Z_e_vs_t_plot.png')

    # Show the plot
    plt.show()


    # Plotting u_dot, v_dot, w_dot vs t
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot u_dot on the left y-axis
    ax1.plot(t, u_dot, label=r'$\dot{u}$ vs. $t$', color=color_u_dot, linestyle='-', linewidth=1.5)
    ax1.set_xlabel(r'Time $t$ [s]')
    ax1.set_ylabel(r'Accelerations [m/s²]', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True)

    # Plot v_dot on the left y-axis
    ax1.plot(t, v_dot, label=r'$\dot{v}$ vs. $t$', color=color_v_dot, linestyle='-', linewidth=1.5)

    # Plot w_dot on the left y-axis
    ax1.plot(t, w_dot, label=r'$\dot{w}$ vs. $t$', color=color_w_dot, linestyle='-', linewidth=1.5)

    # Plot Vt on the left y-axis
    ax1.plot(t, Vt_dot, label=r'$\dot{Vt}$ vs. $t$', color='k', linestyle=':', linewidth=1.5)
    # plt.xlim([40,50])
    # Title and legends
    plt.title(r'$\dot{u}$, $\dot{v}$, $\dot{w}$, and $\dot{Vt}$ vs. $t$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A05_u_dot_v_dot_w_dot_Vt_dot_vs_t_plot.png')

    # Show the plot
    plt.show()


    # Plotting smoothed p_dot, q_dot, r_dot vs t
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot p_dot on the left y-axis
    ax1.plot(t, p_dot, label=r'$\dot{p}$ vs. $t$', color=color_p_dot, linestyle='-', linewidth=1.5)
    ax1.set_xlabel(r'Time $t$ [s]')
    ax1.set_ylabel(r'Angular Accelerations [rad/s²]', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True)

    # Plot q_dot on the left y-axis
    ax1.plot(t, q_dot, label=r'$\dot{q}$ vs. $t$', color=color_q_dot, linestyle='-', linewidth=1.5)

    # Plot r_dot on the left y-axis
    ax1.plot(t, r_dot, label=r'$\dot{r}$ vs. $t$', color=color_r_dot, linestyle='-', linewidth=1.5)

    # Title and legends
    plt.title(r'$\dot{p}$, $\dot{q}$, $\dot{r}$ vs. $t$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A06_smoothed_p_dot_q_dot_r_dot_vs_t_plot.png')

    # Show the plot
    plt.show()


    # Plotting smoothed phi_dot, theta_dot, psi_dot vs t
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot phi_dot on the left y-axis
    ax1.plot(t, phi_dot, label=r'$\dot{\phi}$ vs. $t$', color=color_phi_dot, linestyle='-', linewidth=1.5)
    ax1.set_xlabel(r'Time $t$ [s]')
    ax1.set_ylabel(r'Angular Rates [rad/s]', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True)

    # Plot theta_dot on the left y-axis
    ax1.plot(t, theta_dot, label=r'$\dot{\theta}$ vs. $t$', color=color_theta_dot, linestyle='-', linewidth=1.5)

    # Plot psi_dot on the left y-axis
    ax1.plot(t, psi_dot, label=r'$\dot{\psi}$ vs. $t$', color=color_psi_dot, linestyle='-', linewidth=1.5)

    # Title and legends
    plt.title(r'$\dot{\phi}$, $\dot{\theta}$, $\dot{\psi}$ vs. $t$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A07_smoothed_phi_dot_theta_dot_psi_dot_vs_t_plot.png')

    # Show the plot
    plt.show()


    # Plotting smoothed X_e_dot, Y_e_dot, Z_e_dot vs t
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot X_e_dot on the left y-axis
    ax1.plot(t, X_e_dot, label=r'$\dot{X}_e$ vs. $t$', color=color_X_e_dot, linestyle='-', linewidth=1.5)
    ax1.set_xlabel(r'Time $t$ [s]')
    ax1.set_ylabel(r'Position Rates [m/s]', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True)

    # Plot Y_e_dot on the left y-axis
    ax1.plot(t, Y_e_dot, label=r'$\dot{Y}_e$ vs. $t$', color=color_Y_e_dot, linestyle='-', linewidth=1.5)

    # Plot Z_e_dot on the left y-axis
    ax1.plot(t, Z_e_dot, label=r'$\dot{Z}_e$ vs. $t$', color=color_Z_e_dot, linestyle='-', linewidth=1.5)

    # Title and legends
    plt.title(r'$\dot{X}_e$, $\dot{Y}_e$, $\dot{Z}_e$ vs. $t$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A08_smoothed_X_e_dot_Y_e_dot_Z_e_dot_vs_t_plot.png')

    # Show the plot
    plt.show()


    # Plotting smoothed X_e_dot, Y_e_dot, Z_e_dot vs t
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot X_e_dot on the left y-axis
    ax1.plot(t, Z_e, label=r'$Z_e$ vs. $t$', color=color_X_e_dot, linestyle='-', linewidth=1.5)
    ax1.set_xlabel(r'Time $t$ [s]')
    ax1.set_ylabel(r'Altitude [m]', color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True)

    # Plot Y_e_dot on the left y-axis
    ax1.plot(t, Z_e_raw, label=r'$Z_e$ Raw Data vs. $t$', color='k', linestyle='--', linewidth=1.5)

    # Title and legends
    plt.title(r'$Z_e$ vs. $t$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.5))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A09_Ze_Ze_raw_vs_t.png')

    # Show the plot
    plt.show()


################################################################################################################################################################################
################################################################################################################################################################################
if whether_to_plot_3D_flight_path == 1:

    # Plotting the 3D flight path in North-East-Down (NED) axes
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the flight path
    ax.plot(X_e_in_NED, Y_e_in_NED, Z_e_in_NED, label='Flight Path', color='b', linewidth=2)

    # Set labels and title
    ax.set_xlabel(r'$X_e$ (+ for North) [m]')
    ax.set_ylabel(r'$Y_e$ (+ for East) [m]')
    ax.set_zlabel(r'$Z_e$ (+ for Down) [m]')
    ax.set_title('3D Flight Path (in NED axes)')

    # User input
    whether_to_plot_3D_equal_x_y_z_scale = 1

    if whether_to_plot_3D_equal_x_y_z_scale == 1:

        # Ensure the aspect ratio is equal by setting the limits for each axis
        max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0

        mid_x = (X_e_in_NED.max() + X_e_in_NED.min()) * 0.5
        mid_y = (Y_e_in_NED.max() + Y_e_in_NED.min()) * 0.5
        mid_z = (Z_e_in_NED.max() + Z_e_in_NED.min()) * 0.5

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # # Invert the Z-axis for better visualization (if necessary)
    # ax.invert_zaxis()

    # Show grid
    ax.grid(True)

    # Show legend
    ax.legend()

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A10_3D_flight_path_in_NED.png')

    # Show the plot
    plt.show()


    # Plotting the 3D flight path in North-West-Up (NWU) axes
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the flight path
    ax.plot(X_e, Y_e, Z_e, label='Flight Path', color='b', linewidth=2)

    # Set labels and title
    ax.set_xlabel(r'$X_e$ (+ for North) [m]')
    ax.set_ylabel(r'$Y_e$ (+ for West) [m]')
    ax.set_zlabel(r'$Z_e$ (+ for Up) [m]')
    ax.set_title('3D Flight Path (in NWU axes)')

    if whether_to_plot_3D_equal_x_y_z_scale == 1:

        # Ensure the aspect ratio is equal by setting the limits for each axis
        max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0

        mid_x = (X_e.max() + X_e.min()) * 0.5
        mid_y = (Y_e.max() + Y_e.min()) * 0.5
        mid_z = (Z_e.max() + Z_e.min()) * 0.5

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # # Invert the Z-axis for better visualization (if necessary)
    # ax.invert_zaxis()

    # Show grid
    ax.grid(True)

    # Show legend
    ax.legend()

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A10_3D_flight_path_in_NWU.png')

    # Show the plot
    plt.show()


################################################################################################################################################################################
################################################################################################################################################################################
if whether_to_create_3D_flight_path_animation == 1:

    # # Animation v5
    # Interpolate the data for smoother animation
    frames = 1000  # Increase the number of frames for smoother animation
    interp_func_x = interp1d(t, X_e, kind='linear')
    interp_func_y = interp1d(t, Y_e, kind='linear')
    interp_func_z = interp1d(t, Z_e, kind='linear')
    t_interp = np.linspace(t[0], t[-1], frames)
    X_e_interp = interp_func_x(t_interp)
    Y_e_interp = interp_func_y(t_interp)
    Z_e_interp = interp_func_z(t_interp)

    # Calculate the interval between frames to match real-time duration
    total_duration = t[-1] - t[0]  # Total flight time in seconds
    interval = total_duration / frames * 1000  # Interval in milliseconds

    # Create a figure for the animation
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Set up the plot limits to be equal
    max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
    mid_x = (X_e.max() + X_e.min()) * 0.5
    mid_y = (Y_e.max() + Y_e.min()) * 0.5
    mid_z = (Z_e.max() + Z_e.min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # # Invert the Z-axis for better visualization (if necessary)
    # ax.invert_zaxis()

    # Initialize a line object and a point object for the aircraft position
    line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)
    point, = ax.plot([], [], [], 'ro', label='Aircraft Position')  # Red dot for the aircraft

    # Function to initialize the background of the animation
    def init():
        line.set_data([], [])
        line.set_3d_properties([])
        point.set_data([], [])
        point.set_3d_properties([])
        return line, point

    # Function to update the animation at each frame
    def update(num, X_e_interp, Y_e_interp, Z_e_interp, line, point):
        line.set_data(X_e_interp[:num], Y_e_interp[:num])
        line.set_3d_properties(Z_e_interp[:num])
        point.set_data(X_e_interp[num], Y_e_interp[num])
        point.set_3d_properties(Z_e_interp[num])
        return line, point

    # Create the animation with the correct interval
    ani = animation.FuncAnimation(fig, update, frames=frames, fargs=(X_e_interp, Y_e_interp, Z_e_interp, line, point),
                                init_func=init, blit=True, interval=interval)

    # Add labels and title
    ax.set_xlabel(r'$X_e$ [m]')
    ax.set_ylabel(r'$Y_e$ [m]')
    ax.set_zlabel(r'$Z_e$ [m]')
    ax.set_title('3D Flight Path Animation with Smooth Motion')
    ax.legend()

    # Save the animation as an mp4 file
    ani.save(xplane11_path + 'richard_checked_terms\\A11_3D_flight_path_animation_smooth.mp4', writer='ffmpeg')

    # Display the animation
    plt.show()




    # # Animation v1
    # # Create a figure for the animation
    # fig = plt.figure(figsize=(10, 7))
    # ax = fig.add_subplot(111, projection='3d')

    # # Set up the plot limits to be equal
    # max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
    # mid_x = (X_e.max() + X_e.min()) * 0.5
    # mid_y = (Y_e.max() + Y_e.min()) * 0.5
    # mid_z = (Z_e.max() + Z_e.min()) * 0.5

    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_y - max_range, mid_y + max_range)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # # # Invert the Z-axis for better visualization (if necessary)
    # # ax.invert_zaxis()

    # # Initialize a line object which will be updated during the animation
    # line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)

    # # Function to initialize the background of the animation
    # def init():
    #     line.set_data([], [])
    #     line.set_3d_properties([])
    #     return line,

    # # Function to update the animation at each frame
    # def update(num, X_e, Y_e, Z_e, line):
    #     line.set_data(X_e[:num], Y_e[:num])
    #     line.set_3d_properties(Z_e[:num])
    #     return line,

    # # Create the animation
    # ani = animation.FuncAnimation(fig, update, frames=len(t), fargs=(X_e, Y_e, Z_e, line),
    #                               init_func=init, blit=True, interval=30)

    # # Add labels and title
    # ax.set_xlabel(r'$X_e$ [m]')
    # ax.set_ylabel(r'$Y_e$ [m]')
    # ax.set_zlabel(r'$Z_e$ [m]')
    # ax.set_title('3D Flight Path Animation')
    # ax.legend()

    # # Save the animation as an mp4 file
    # ani.save(xplane11_path + 'richard_checked_terms\\3D_flight_path_animation.mp4', writer='ffmpeg')

    # # Display the animation
    # plt.show()







    # # Animation v2
    # # Calculate the time intervals between frames in milliseconds
    # time_intervals = np.diff(t) * 1000  # Convert from seconds to milliseconds

    # # Create a figure for the animation
    # fig = plt.figure(figsize=(10, 7))
    # ax = fig.add_subplot(111, projection='3d')

    # # Set up the plot limits to be equal
    # max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
    # mid_x = (X_e.max() + X_e.min()) * 0.5
    # mid_y = (Y_e.max() + Y_e.min()) * 0.5
    # mid_z = (Z_e.max() + Z_e.min()) * 0.5

    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_y - max_range, mid_y + max_range)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # # # Invert the Z-axis for better visualization (if necessary)
    # # ax.invert_zaxis()

    # # Initialize a line object which will be updated during the animation
    # line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)

    # # Function to initialize the background of the animation
    # def init():
    #     line.set_data([], [])
    #     line.set_3d_properties([])
    #     return line,

    # # Function to update the animation at each frame
    # def update(num, X_e, Y_e, Z_e, line):
    #     line.set_data(X_e[:num], Y_e[:num])
    #     line.set_3d_properties(Z_e[:num])
    #     return line,

    # # Create the animation with real-time progression
    # ani = animation.FuncAnimation(fig, update, frames=len(t), fargs=(X_e, Y_e, Z_e, line),
    #                               init_func=init, blit=True, interval=time_intervals[0])

    # # Add labels and title
    # ax.set_xlabel(r'$X_e$ [m]')
    # ax.set_ylabel(r'$Y_e$ [m]')
    # ax.set_zlabel(r'$Z_e$ [m]')
    # ax.set_title('3D Flight Path Animation')
    # ax.legend()

    # # Save the animation as an mp4 file
    # ani.save(xplane11_path + 'richard_checked_terms\\3D_flight_path_animation_real_time.mp4', writer='ffmpeg')

    # # Display the animation
    # plt.show()






    # # Animation v3
    # # Calculate the total duration and number of frames
    # total_duration = t[-1] - t[0]  # Total flight time in seconds (50 seconds in this case)
    # num_frames = len(t)  # Total number of frames

    # # Calculate the interval between frames to match real-time duration
    # interval = total_duration / num_frames * 1000  # Interval in milliseconds

    # # Create a figure for the animation
    # fig = plt.figure(figsize=(10, 7))
    # ax = fig.add_subplot(111, projection='3d')

    # # Set up the plot limits to be equal
    # max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
    # mid_x = (X_e.max() + X_e.min()) * 0.5
    # mid_y = (Y_e.max() + Y_e.min()) * 0.5
    # mid_z = (Z_e.max() + Z_e.min()) * 0.5

    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_y - max_range, mid_y + max_range)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # # # Invert the Z-axis for better visualization (if necessary)
    # # ax.invert_zaxis()

    # # Initialize a line object which will be updated during the animation
    # line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)

    # # Function to initialize the background of the animation
    # def init():
    #     line.set_data([], [])
    #     line.set_3d_properties([])
    #     return line,

    # # Function to update the animation at each frame
    # def update(num, X_e, Y_e, Z_e, line):
    #     line.set_data(X_e[:num], Y_e[:num])
    #     line.set_3d_properties(Z_e[:num])
    #     return line,

    # # Create the animation with the correct interval
    # ani = animation.FuncAnimation(fig, update, frames=num_frames, fargs=(X_e, Y_e, Z_e, line),
    #                               init_func=init, blit=True, interval=interval)

    # # Add labels and title
    # ax.set_xlabel(r'$X_e$ [m]')
    # ax.set_ylabel(r'$Y_e$ [m]')
    # ax.set_zlabel(r'$Z_e$ [m]')
    # ax.set_title('3D Flight Path Animation')
    # ax.legend()

    # # Save the animation as an mp4 file
    # ani.save(xplane11_path + 'richard_checked_terms\\3D_flight_path_animation_real_time_fixed.mp4', writer='ffmpeg')

    # # Display the animation
    # plt.show()






    # # Animation v4
    # # Calculate the total duration and number of frames
    # total_duration = t[-1] - t[0]  # Total flight time in seconds (50 seconds in this case)
    # num_frames = len(t)  # Total number of frames

    # # Calculate the interval between frames to match real-time duration
    # interval = total_duration / num_frames * 1000  # Interval in milliseconds

    # # Create a figure for the animation
    # fig = plt.figure(figsize=(10, 7))
    # ax = fig.add_subplot(111, projection='3d')

    # # Set up the plot limits to be equal
    # max_range = np.array([X_e.max()-X_e.min(), Y_e.max()-Y_e.min(), Z_e.max()-Z_e.min()]).max() / 2.0
    # mid_x = (X_e.max() + X_e.min()) * 0.5
    # mid_y = (Y_e.max() + Y_e.min()) * 0.5
    # mid_z = (Z_e.max() + Z_e.min()) * 0.5

    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_y - max_range, mid_y + max_range)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # # Invert the Z-axis for better visualization (if necessary)
    # ax.invert_zaxis()

    # # Initialize a line object and a point object for the aircraft position
    # line, = ax.plot([], [], [], label='Flight Path', color='b', linewidth=2)
    # point, = ax.plot([], [], [], 'ro', label='Aircraft Position')  # Red dot for the aircraft

    # # Function to initialize the background of the animation
    # def init():
    #     line.set_data([], [])
    #     line.set_3d_properties([])
    #     point.set_data([], [])
    #     point.set_3d_properties([])
    #     return line, point

    # # Function to update the animation at each frame
    # def update(num, X_e, Y_e, Z_e, line, point):
    #     line.set_data(X_e[:num], Y_e[:num])
    #     line.set_3d_properties(Z_e[:num])
    #     point.set_data(X_e[num], Y_e[num])
    #     point.set_3d_properties(Z_e[num])
    #     return line, point

    # # Create the animation with the correct interval
    # ani = animation.FuncAnimation(fig, update, frames=num_frames, fargs=(X_e, Y_e, Z_e, line, point),
    #                               init_func=init, blit=True, interval=interval)

    # # Add labels and title
    # ax.set_xlabel(r'$X_e$ [m]')
    # ax.set_ylabel(r'$Y_e$ [m]')
    # ax.set_zlabel(r'$Z_e$ [m]')
    # ax.set_title('3D Flight Path Animation')
    # ax.legend()

    # # Save the animation as an mp4 file
    # ani.save(xplane11_path + 'richard_checked_terms\\3D_flight_path_animation_real_time_with_dot.mp4', writer='ffmpeg')

    # # Display the animation
    # plt.show()


################################################################################################################################################################################
################################################################################################################################################################################
if whether_to_plot_U == 1:

    # Create a figure and axis objects with 5 subplots stacked vertically
    fig, axs = plt.subplots(5, 1, figsize=(10, 10), sharex=True)

    U_plots_axes_0_upper_buffer = (0.007/1.00) * (1-0)
    U_plots_axes_0_lower_buffer = (0.001/1.00) * (1-0)

    U_plots_axes_1_upper_buffer = (0.007/1.00) * (9-(-18.4))
    U_plots_axes_1_lower_buffer = (0.001/1.00) * (9-(-18.4))

    U_plots_axes_2_upper_buffer = (0.007/1.00) * (24-(-24))
    U_plots_axes_2_lower_buffer = (0.001/1.00) * (24-(-24))

    U_plots_axes_3_upper_buffer = (0.007/1.00) * (23.21-(-23.21))
    U_plots_axes_3_lower_buffer = (0.001/1.00) * (23.21-(-23.21))

    U_plots_axes_4_upper_buffer = (0.0065/1.00) * (-0-(-40))
    U_plots_axes_4_lower_buffer = (0.001/1.00) * (-0-(-40))



    # Plot Throttle Input
    axs[0].plot(t, delta_t, color='r')
    axs[0].set_ylabel(r'Throttle [-]')
    axs[0].set_ylim([0 - U_plots_axes_0_lower_buffer, 1 + U_plots_axes_0_upper_buffer])
    axs[0].grid(True)
    axs[0].set_title('Control Inputs vs. Time')

    # Plot Elevator Input
    axs[1].plot(t, delta_e, color='g')
    axs[1].set_ylabel(r'Elevator [$^{\circ}$]')
    axs[1].set_ylim([-18.4 - U_plots_axes_1_lower_buffer, 9 + U_plots_axes_1_upper_buffer])
    axs[1].grid(True)

    # Plot Aileron Input
    axs[2].plot(t, delta_a, color='b')
    axs[2].set_ylabel(r'Aileron [$^{\circ}$]')
    axs[2].set_ylim([-24 - U_plots_axes_2_lower_buffer, 24 + U_plots_axes_2_upper_buffer])
    axs[2].grid(True)

    # Plot Rudder Input
    axs[3].plot(t, delta_r, color='m')
    axs[3].set_ylabel(r'Rudder [$^{\circ}$]')
    axs[3].set_ylim([-23.21 - U_plots_axes_3_lower_buffer, 23.21 + U_plots_axes_3_upper_buffer])
    axs[3].grid(True)

    # Plot Flap Input
    axs[4].plot(t, delta_f, color='c')
    axs[4].set_ylabel(r'Flap [$^{\circ}$]')
    axs[4].set_ylim([-40 - U_plots_axes_4_lower_buffer, -0 + U_plots_axes_4_upper_buffer])
    axs[4].set_xlabel('Time [s]')
    axs[4].grid(True)

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A12_control_inputs_vs_time.png')

    # Show the plot
    plt.show()

    # print(min(delta_e))



################################################################################################################################################################################
################################################################################################################################################################################
if whether_to_plot_stage2_combined_history == 1:

    ## stage2_combined_history plot 1

    # Plotting the data
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot X_flow_sep on the left y-axis
    ax1.plot(t, X_flow_sep, label=r'$X_{flow\_sep}$', color=color2)
    ax1.plot(t, Vt/1.943844/50, label=r'TAS/$V_{stall}$', color='k', linestyle=':', linewidth=1.5)
    ax1.set_xlabel('Time [s]')

    # User input
    ax1.set_xlim([0, 50])
    # ax1.set_xlim([15, 26])

    ax1.set_ylabel(r'$X_{flow\_sep}$ [-], TAS/$V_{stall}$ [-]', color=color2)
    ax1.set_ylim([-0.1, 1.1])
    ax1.tick_params(axis='y', labelcolor=color2)
    ax1.grid(True)

    # Create a second y-axis for alpha, alpha_crit, and poststall_on
    ax2 = ax1.twinx()
    ax2.plot(t, alpha, label=r'$\alpha$', color=color1, linestyle='--')
    ax2.plot(t, alpha_crit_1, label=r'$\alpha_{crit\_1}$', color=color3, linestyle='-.')
    ax2.plot(t, alpha_crit_2, label=r'$\alpha_{crit\_2}$', color=color4, linestyle='-.')
    ax2.plot(t, alpha_crit_3, label=r'$\alpha_{crit\_3}$', color=color5, linestyle='-.')
    ax2.plot(t, poststall_on, label=r'poststall$\_$on', color=color7, linestyle='-')
    ax2.set_ylabel(r'$\alpha$ [$^{\circ}$]', color=color1)
    ax2.set_ylim([-15, 45])
    ax2.tick_params(axis='y', labelcolor=color1)

    # Create a third y-axis for CL, Cm, CD
    ax3 = ax1.twinx()
    ax3.spines["right"].set_position(("outward", 60))  # Move the third y-axis outward
    ax3.plot(t, CL, label=r'$C_L$', color=color6, linestyle='-')  # CL is multiplied by 2
    ax3.plot(t, Cm, label=r'$C_m$', color=color8, linestyle='-')  # Cm is multiplied by 2
    ax3.plot(t, CD, label=r'$C_D$', color=color9, linestyle='-')  # CD is multiplied by 5
    ax3.set_ylabel(r'$C_L$, $C_m$, $C_D$', color=color6)
    ax3.set_ylim([-15/10/1.5, 45/10/1.5])
    ax3.tick_params(axis='y', labelcolor=color6)

    # Adding legends
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper center')
    ax3.legend(loc='upper right')

    # Title and layout
    plt.title(r'$X_{flow\_sep}$, TAS/$V_{stall}$, $\alpha$, $C_L$, $C_m$, $C_D$ vs. Time')
    fig.tight_layout()

    # Save and display
    plt.savefig(xplane11_path + 'richard_checked_terms\\A13_X_flow_sep_alpha_CL_Cm_CD_vs_time.png')  # Save as PNG
    plt.show()



    ## stage2_combined_history plot 2

    # User input
    # The range of time that you want data plotted
    # user_defined_t_initial = 15 # sec
    # user_defined_t_final = 20 # sec
    user_defined_t_initial = 0 # sec
    user_defined_t_final = 50 # sec

    # Find the indices where values in t range from 15 to 20
    t_index_range = np.where((t >= user_defined_t_initial) & (t <= user_defined_t_final))[0]
    # If there are valid indices, find the first and last
    if t_index_range.size > 0:
        first_t_index = t_index_range[0]
        last_t_index = t_index_range[-1]
        print(f"[{first_t_index}:{last_t_index}]")

    # # Output the index range
    # print(t_index_range)

    # User input
    # Determine whether to plot in degrees or radians
    plot_alpha_unit_flag = 0

    if plot_alpha_unit_flag == 0:
        alpha_label = r'$\alpha$ [$^{\circ}$]'
        alpha_values = alpha[first_t_index:last_t_index]
    else:
        alpha_label = r'$\alpha$ [rad]'
        alpha_values = alpha * (np.pi/180)

    X_flow_sep_values = X_flow_sep[first_t_index:last_t_index]
    CL_values = CL[first_t_index:last_t_index]
    poststall_on_values = poststall_on[first_t_index:last_t_index]
    Cm_values = Cm[first_t_index:last_t_index]
    CD_values = CD[first_t_index:last_t_index]

    # Plotting the data
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # User input
    ax1.set_xlim([-5, 95])

    # Check if all values are within the range [0, 1]
    if np.all((X_flow_sep_values >= 0) & (X_flow_sep_values <= 1)):
        ax1.set_ylim([-0.1, 1.1])
    else:
        # Set values outside the range of [0, 1] to 1
        X_flow_sep_values[0] = 1
        X_flow_sep_values = np.clip(X_flow_sep_values, 0, 1)
        ax1.set_ylim([-0.1, 1.1])

    # Plot X_flow_sep on the left y-axis
    ax1.plot(alpha_values, X_flow_sep_values, label=r'$X_{flow\_sep}$ vs. $\alpha$', color=color2)
    ax1.set_xlabel(alpha_label)
    ax1.set_ylabel(r'$X_{flow\_sep}$ [-]', color=color2)
    ax1.tick_params(axis='y', labelcolor=color2)
    ax1.grid(True)

    # Create a second y-axis for CL, Cm, CD, and poststall_on
    ax2 = ax1.twinx()

    # User input
    # ax2.set_ylim([-4.5, 2.9]) # For static comparison
    ax2.set_ylim([-4.6, 3.9]) # For dynamic comparison

    ax2.plot(alpha_values, CL_values, label=r'$C_L$ vs. $\alpha$', color=color6, linestyle='-')
    ax2.plot(alpha_values, poststall_on_values, label=r'poststall$\_$on vs. $\alpha$', color=color7, linestyle='-')
    ax2.plot(alpha_values, Cm_values, label=r'$C_m$ vs. $\alpha$', color=color8, linestyle='-')
    ax2.plot(alpha_values, CD_values, label=r'$C_D$ vs. $\alpha$', color=color9, linestyle='-')
    ax2.set_ylabel(r'$C_L$, $C_m$, $C_D$ [-]', color=color1)
    ax2.tick_params(axis='y', labelcolor=color1)

    # Title and legends
    plt.title(r'$X_{flow\_sep}$, $C_L$, $C_m$, $C_D$ vs. $\alpha$')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.35))

    # Save and display
    plt.savefig(xplane11_path + 'richard_checked_terms\\A14_X_flow_sep_CL_Cm_CD_vs_alpha.png')  # Save as PNG
    plt.show()


################################################################################################################################################################################
################################################################################################################################################################################
if whether_to_plot_stage3_combined_history == 1:

    #----------------------------------------------------------------------------------------------
    # User input
    # The range of time that you want data plotted
    user_defined_t_initial = 15 # sec
    user_defined_t_final = 20 # sec
    # user_defined_t_initial = 0 # sec
    # user_defined_t_final = 50 # sec

    # Find the indices where values in t range from 15 to 20
    t_index_range = np.where((t >= user_defined_t_initial) & (t <= user_defined_t_final))[0]
    # If there are valid indices, find the first and last
    if t_index_range.size > 0:
        first_t_index = t_index_range[0]
        last_t_index = t_index_range[-1]
        print(f"[{first_t_index}:{last_t_index}]")
    #----------------------------------------------------------------------------------------------

    # whether showing a now useless variable dt that is not actually used
    # User input
    whether_to_plot_dt = 0

    if whether_to_plot_dt == 1:
        # Plotting dt over time
        plt.figure(figsize=(10, 6))
        plt.plot(t, dt, label=r'$\Delta t$ vs. Time', color=color_dt, linestyle='-', linewidth=1.5)
        plt.xlabel('Time [s]')
        plt.ylabel(r'$\Delta t$ [s]')
        plt.title(r'$\Delta t$ vs. Time')
        plt.grid(True)

        # Add legend
        plt.legend()

        # Save the figure
        plt.savefig(xplane11_path + 'richard_checked_terms\\A15_dt_vs_time.png')

        # Show the plot
        plt.show()

    # Create the figure and axis objects
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot rho on the left y-axis
    ax1.plot(t, rho, label=r'$\rho$ vs. Time', color=color_rho, linestyle='-', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'$\rho$ [kg/m$^3$]', color=color_rho)
    ax1.tick_params(axis='y', labelcolor=color_rho)
    ax1.grid(True)

    # Create a second y-axis for a_sound_speed
    ax2 = ax1.twinx()
    ax2.plot(t, a_sound_speed, label=r'$a$ (Speed of Sound) vs. Time', color=color_a, linestyle='--', linewidth=1.5)
    ax2.set_ylabel(r'$a$ [m/s]', color=color_a)
    ax2.tick_params(axis='y', labelcolor=color_a)

    # Title and legends
    plt.title(r'$\rho$ and $a$ (Speed of Sound) vs. Time')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A16_rho_vs_a_sound_speed_vs_time.png')

    # Show the plot
    plt.show()


    # Create the figure and axis objects
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot alpha on the left y-axis
    ax1.plot(t, alpha, label=r'$\alpha$ vs. Time', color=color_alpha, linestyle='-', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'$\alpha$ [$^{\circ}$]', color=color_alpha)  # Assuming alpha is in degrees
    ax1.tick_params(axis='y', labelcolor=color_alpha)
    ax1.grid(True)

    # Create a second y-axis for beta
    ax2 = ax1.twinx()
    ax2.plot(t, beta, label=r'$\beta$ vs. Time', color=color_beta, linestyle='--', linewidth=1.5)
    ax2.set_ylabel(r'$\beta$ [$^{\circ}$]', color=color_beta)  # Assuming beta is in degrees
    ax2.tick_params(axis='y', labelcolor=color_beta)

    # Title and legends
    plt.title(r'$\alpha$ and $\beta$ vs. Time')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A17_alpha_vs_beta_vs_time.png')

    # Show the plot
    plt.show()


    # Create the figure and axis objects
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot Cx, Cy, Cz on the left y-axis
    ax1.plot(t, Cx, label=r'$C_x$ vs. Time', color=color_Cx, linestyle='-', linewidth=1.5)
    ax1.plot(t, Cy, label=r'$C_y$ vs. Time', color=color_Cy, linestyle='-', linewidth=1.5)
    ax1.plot(t, Cz, label=r'$C_z$ vs. Time', color=color_Cz, linestyle='-', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'Force Coefficients ($C_x$, $C_y$, $C_z$)', color=color_Cx)
    ax1.tick_params(axis='y', labelcolor=color_Cx)
    ax1.grid(True)

    # Create a second y-axis for Cl, Cm, Cn
    ax2 = ax1.twinx()
    ax2.plot(t, Cl, label=r'$C_l$ vs. Time', color=color_Cl, linestyle='--', linewidth=1.5)
    ax2.plot(t, Cm, label=r'$C_m$ vs. Time', color=color_Cm, linestyle=':', linewidth=1.5)
    ax2.plot(t, Cn, label=r'$C_n$ vs. Time', color=color_Cn, linestyle='--', linewidth=1.5)
    ax2.set_ylabel(r'Moment Coefficients ($C_l$, $C_m$, $C_n$)', color=color_Cl)
    ax2.tick_params(axis='y', labelcolor=color_Cl)

    # Title and legends
    plt.title(r'Force and Moment Coefficients vs. Time')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A18_force_moment_coefficients_vs_time.png')

    # Show the plot
    plt.show()


    # Create the figure and axis objects
    plt.figure(figsize=(10, 6))

    # Plot Cl vs beta
    plt.plot(beta[first_t_index:last_t_index], Cl[first_t_index:last_t_index], label=r'$C_l$ vs. $\beta$', color=color_Cl, linestyle='-', linewidth=1.5)

    # Plot Cn vs beta
    plt.plot(beta[first_t_index:last_t_index], Cn[first_t_index:last_t_index], label=r'$C_n$ vs. $\beta$', color=color_Cn, linestyle='--', linewidth=1.5)

    # Add labels and title
    plt.xlabel(r'$\beta$ [deg]')  # Assuming beta is in degrees
    plt.ylabel(r'$C_l$, $C_n$ [-]')
    plt.title(r'$C_l$ and $C_n$ vs. $\beta$')

    # Add grid and legend
    plt.grid(True)
    plt.legend()

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A19_Cl_Cn_vs_beta.png')

    # Show the plot
    plt.show()

    
    # Plot 1: P_prop_power vs. time and n_prop_rot_speed vs. time
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(t, P_prop_power, label=r'$P_{prop\_power}$', color='r', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'$P_{prop\_power}$ [W]', color='r')
    ax1.tick_params(axis='y', labelcolor='r')
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(t, n_prop_rot_speed, label=r'$n_{prop\_rot\_speed}$', color='b', linestyle='--', linewidth=1.5)
    ax2.set_ylabel(r'$n_{prop\_rot\_speed}$ [rpm]', color='b')
    ax2.tick_params(axis='y', labelcolor='b')

    plt.title('Plot 1: Propeller Power and Rotational Speed vs. Time')
    fig.tight_layout()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A20_plot1.png')
    plt.show()

    # Plot 2: P_prop_power vs. n_prop_rot_speed
    plt.figure(figsize=(10, 6))
    plt.plot(n_prop_rot_speed[first_t_index:last_t_index], P_prop_power[first_t_index:last_t_index], color='r', linewidth=1.5, label=r'$P_{prop\_power}$ vs. $n_{prop\_rot\_speed}$')
    plt.xlabel(r'$n_{prop\_rot\_speed}$ [rpm]')
    plt.ylabel(r'$P_{prop\_power}$ [W]')
    plt.grid(True)
    plt.title('Plot 2: P_prop_power vs. n_prop_rot_speed')
    plt.legend()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A21_plot2.png')
    plt.show()

    # Plot 3: n_prop_rot_speed and global_n_prop_rot_speed on left axis, n_prop_rot_speed_dot on right axis
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(t, n_prop_rot_speed, label=r'$n_{prop\_rot\_speed}$', color='b', linewidth=1.5)
    ax1.plot(t, global_n_prop_rot_speed, label=r'$global\_n_{prop\_rot\_speed}$', color='g', linestyle='--', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'$n_{prop\_rot\_speed}$ [rpm]', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(t, n_prop_rot_speed_dot, label=r'$n_{prop\_rot\_speed\_dot}$', color='r', linestyle='-', linewidth=1.5)
    ax2.set_ylabel(r'$n_{prop\_rot\_speed\_dot}$ [rpm/s]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    plt.title('Plot 3: Propeller Speeds and Rotational Acceleration vs. Time')
    fig.tight_layout()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A22_plot3.png')
    plt.show()

    # Plot 4: J_advance_ratio on left axis, eta_p_prop_eff on right axis
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(t, J_advance_ratio, label=r'$J_{advance\_ratio}$', color='b', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'$J_{advance\_ratio}$ [-]', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(t, eta_p_prop_eff, label=r'$\eta_{p\_prop\_eff}$', color='r', linestyle='-', linewidth=1.5)
    ax2.set_ylabel(r'$\eta_{p\_prop\_eff}$ [-]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    plt.title('Plot 4: Advance Ratio and Propeller Efficiency vs. Time')
    fig.tight_layout()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A23_plot4.png')
    plt.show()

    # Plot 5: eta_p_prop_eff, Vt vs. J_advance_ratio
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(J_advance_ratio[first_t_index:last_t_index], eta_p_prop_eff[first_t_index:last_t_index], label=r'$\eta_{p\_prop\_eff}$ vs. $J_{advance\_ratio}$', color='b', linewidth=1.5)
    ax1.set_xlabel(r'$J_{advance\_ratio}$ [-]')
    ax1.set_ylabel(r'$\eta_{p\_prop\_eff}$ [-]', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(J_advance_ratio[first_t_index:last_t_index], Vt[first_t_index:last_t_index], label=r'$V_t$ vs. $J_{advance\_ratio}$', color='r', linestyle='-', linewidth=1.5)
    ax2.set_ylabel(r'$V_t$ [m/s]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    plt.title('Plot 5: Propeller Efficiency and Airspeed vs. Advance Ratio')
    fig.tight_layout()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A24_plot5.png')
    plt.show()

    # Plot 6: ay_noise and az_noise vs. time
    plt.figure(figsize=(10, 6))
    plt.plot(t[first_t_index:last_t_index], ay_noise[first_t_index:last_t_index], label=r'$a_y\_noise$', color='b', linewidth=1.5)
    plt.plot(t[first_t_index:last_t_index], az_noise[first_t_index:last_t_index], label=r'$a_z\_noise$', color='r', linestyle='--', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.ylabel('Acceleration Noise [m/s²]')
    plt.grid(True)
    plt.title('Plot 6: Lateral and Vertical Acceleration Noise vs. Time')
    plt.legend()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A25_plot6.png')
    plt.show()


################################################################################################################################################################################
################################################################################################################################################################################
if whether_to_plot_precision_timing == 1:

    #----------------------------------------------------------------------------------------------
    # User input
    # The range of time that you want data plotted
    user_defined_t_initial = 5 # sec
    user_defined_t_final = 50 # sec
    # user_defined_t_initial = 0 # sec
    # user_defined_t_final = 50 # sec

    # Find the indices where values in t range from 15 to 20
    t_index_range = np.where((t >= user_defined_t_initial) & (t <= user_defined_t_final))[0]
    # If there are valid indices, find the first and last
    if t_index_range.size > 0:
        first_t_index = t_index_range[0]
        last_t_index = t_index_range[-1]
        print(f"[{first_t_index}:{last_t_index}]")
    #----------------------------------------------------------------------------------------------

    # User input:
    whether_to_plot_effects_separately = 0

    if whether_to_plot_effects_separately == 1:

        # Plot B1: delta_f_to_f_clock on the left axis and delta_f_to_f_doppler on the right axis
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(t, delta_f_to_f_clock, label=r'$\Delta f/f_{clock}$', color='r', linewidth=1.5)
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel(r'$\Delta f/f_{clock}$', color='r')
        ax1.tick_params(axis='y', labelcolor='r')
        ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.plot(t, delta_f_to_f_doppler, label=r'$\Delta f/f_{doppler}$', color='b', linestyle='--', linewidth=1.5)
        ax2.set_ylabel(r'$\Delta f/f_{doppler}$', color='b')
        ax2.tick_params(axis='y', labelcolor='b')

        plt.title('Plot B1: Atomic Clock Drift and Doppler Shift vs. Time')
        fig.tight_layout()
        plt.savefig(xplane11_path + 'richard_checked_terms\\A26_B1.png')
        plt.show()

        # Plot B2: delta_f_to_f_special_relativity on the left axis and delta_f_to_f_gravity, delta_f_to_f_acceleration on the right axis
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(t, delta_f_to_f_special_relativity, label=r'$\Delta f/f_{special\_relativity}$', color='r', linewidth=1.5)
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel(r'$\Delta f/f_{special\_relativity}$', color='r')
        ax1.tick_params(axis='y', labelcolor='r')
        ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.plot(t, delta_f_to_f_gravity, label=r'$\Delta f/f_{gravity}$', color='b', linestyle='--', linewidth=1.5)
        ax2.plot(t, delta_f_to_f_acceleration, label=r'$\Delta f/f_{acceleration}$', color='g', linestyle='--', linewidth=1.5)
        ax2.set_ylabel(r'$\Delta f/f_{gravity}$, $\Delta f/f_{acceleration}$', color='b')
        ax2.tick_params(axis='y', labelcolor='b')

        plt.title('Plot B2: Special Relativity, Gravity, and Acceleration Effects vs. Time')
        fig.tight_layout()
        plt.savefig(xplane11_path + 'richard_checked_terms\\A27_B2.png')
        plt.show()

        # Plot B3: delta_T_to_T_frame_dragging vs. time
        plt.figure(figsize=(10, 6))
        plt.plot(t, delta_T_to_T_frame_dragging, label=r'$\Delta T/T_{frame\_dragging}$', color='r', linewidth=1.5)
        plt.xlabel('Time [s]')
        plt.ylabel(r'$\Delta T/T_{frame\_dragging}$')
        plt.grid(True)
        plt.title('Plot B3: Frame Dragging Effect vs. Time')
        plt.legend()
        plt.savefig(xplane11_path + 'richard_checked_terms\\A28_B3.png')
        plt.show()

        # Plot B4: delta_T_to_T_light_speed_reduction and delta_T_to_T_tropospheric vs. time
        plt.figure(figsize=(10, 6))
        plt.plot(t[first_t_index:last_t_index], delta_T_to_T_light_speed_reduction[first_t_index:last_t_index], label=r'$\Delta T/T_{light\_speed\_reduction}$', color='r', linewidth=1.5)
        plt.plot(t[first_t_index:last_t_index], delta_T_to_T_tropospheric[first_t_index:last_t_index], label=r'$\Delta T/T_{tropospheric}$', color='b', linestyle='--', linewidth=1.5)
        plt.xlabel('Time [s]')
        plt.ylabel(r'$\Delta T/T$')
        plt.grid(True)
        plt.title('Plot B4: Light Speed Reduction and Tropospheric Delay vs. Time')
        plt.legend()
        plt.savefig(xplane11_path + 'richard_checked_terms\\A29_B4.png')
        plt.show()

        # Plot B5: delta_T_to_T_sagnac on the left axis, delta_T_to_T_ionospheric on the right axis
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(t, delta_T_to_T_sagnac, label=r'$\Delta T/T_{sagnac}$', color='r', linewidth=1.5)
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel(r'$\Delta T/T_{sagnac}$', color='r')
        ax1.tick_params(axis='y', labelcolor='r')
        ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.plot(t, delta_T_to_T_ionospheric, label=r'$\Delta T/T_{ionospheric}$', color='b', linestyle='--', linewidth=1.5)
        ax2.set_ylabel(r'$\Delta T/T_{ionospheric}$', color='b')
        ax2.tick_params(axis='y', labelcolor='b')

        plt.title('Plot B5: Sagnac Effect and Ionospheric Delay vs. Time')
        fig.tight_layout()
        plt.savefig(xplane11_path + 'richard_checked_terms\\A30_B5.png')
        plt.show()

    # User input
    whether_to_plot_instantaneous_effects_at_each_moment = 0

    if whether_to_plot_instantaneous_effects_at_each_moment == 1:

        # Plotting the data
        fig, ax1 = plt.subplots(figsize=(10, 6))

        # Plot total_real_delta_f_to_f on the left y-axis
        ax1.plot(t, total_real_delta_f_to_f, label=r'Total Real $\Delta f/f$', color='b', linewidth=1.5)
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel(r'Total Real $\Delta f/f$', color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        ax1.grid(True)

        # Create the first right y-axis for total_apparent_delta_f_to_f
        ax2 = ax1.twinx()
        ax2.plot(t, total_apparent_delta_f_to_f, label=r'Total Apparent $\Delta f/f$', color='r', linestyle='--', linewidth=1.5)
        ax2.set_ylabel(r'Total Apparent $\Delta f/f$', color='r')
        ax2.tick_params(axis='y', labelcolor='r')

        # Create a second right y-axis for total_delta_T_to_T
        ax3 = ax1.twinx()
        ax3.spines["right"].set_position(("outward", 60))  # Move the second right y-axis outward
        ax3.plot(t, total_delta_T_to_T, label=r'Total $\Delta T/T$', color='g', linestyle='-.', linewidth=1.5)
        ax3.set_ylabel(r'Total $\Delta T/T$', color='g')
        ax3.tick_params(axis='y', labelcolor='g')

        # Title and legends
        plt.title('Total Real and Apparent Frequency Shift, and Total Delay vs. Time')
        fig.tight_layout()

        # Add legends
        ax1.legend(loc='upper left')
        ax2.legend(loc='upper right')
        ax3.legend(loc='lower right')

        # Save the figure
        plt.savefig(xplane11_path + 'richard_checked_terms\\A31_total_effects_plot.png')

        # Show the plot
        plt.show()


    # Unit: second
    # Plotting the data
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot total_real_time_drift on the left y-axis
    ax1.plot(t, total_real_time_drift, label=r'Total Real Time Drift', color='b', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'Total Real Time Drift [s]', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True)

    # # Adjust the position of the left y-axis label
    # ax1.yaxis.set_label_coords(-0.1, 0.5)

    # Create the first right y-axis for total_apparent_time_drift
    ax2 = ax1.twinx()
    ax2.plot(t, total_apparent_time_drift, label=r'Total Apparent Time Drift', color='r', linestyle='--', linewidth=1.5)
    ax2.plot(t, combined_time_drift, label=r'Combined Time Drift', color='k', linestyle=':', linewidth=1.5)
    ax2.set_ylabel(r'Total Apparent Time Drift [s]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    # # Adjust the position of the first right y-axis label
    # ax2.yaxis.set_label_coords(1.05, 0.5)

    # Create a second right y-axis for total_signal_path_induced_time_drift
    ax3 = ax1.twinx()
    ax3.spines["right"].set_position(("outward", 60))  # Move the second right y-axis outward
    ax3.plot(t, total_signal_path_induced_time_drift, label=r'Signal Path Induced Time Drift', color='g', linestyle='-.', linewidth=1.5)
    ax3.set_ylabel(r'Signal Path Induced Time Drift [s]', color='g')
    ax3.tick_params(axis='y', labelcolor='g')

    # # Adjust the position of the second right y-axis label
    # ax3.yaxis.set_label_coords(1.15, 0.5)

    # Adjust scientific notation offset text positions
    # ax1.get_yaxis().get_offset_text().set_x(-0.15)  # Move the offset text of the left y-axis
    ax2.get_yaxis().get_offset_text().set_x(1.05)   # Move the offset text of the first right y-axis
    ax3.get_yaxis().get_offset_text().set_x(1.15)   # Move the offset text of the second right y-axis

    # Title and legends
    plt.title('Time Drift: Real, Apparent, and Signal Path Induced vs. Time')
    fig.tight_layout()

    # Add legends
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax3.legend(loc='lower right')

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A32_time_drift_plot_in_seconds.png')

    # Show the plot
    plt.show()


    # Unit: nanosecond
    # Plotting the data
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot total_real_time_drift on the left y-axis
    ax1.plot(t, (10**9) * total_real_time_drift, label=r'Total Real Time Drift', color='b', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'Total Real Time Drift [ns]', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True)

    # # Adjust the position of the left y-axis label
    # ax1.yaxis.set_label_coords(-0.1, 0.5)

    # Create the first right y-axis for total_apparent_time_drift
    ax2 = ax1.twinx()
    ax2.plot(t, (10**9) * total_apparent_time_drift, label=r'Total Apparent Time Drift', color='r', linestyle='--', linewidth=1.5)
    ax2.plot(t, (10**9) * combined_time_drift, label=r'Combined Time Drift', color='k', linestyle=':', linewidth=1.5)
    ax2.set_ylabel(r'Total Apparent Time Drift [ns]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    # # Adjust the position of the first right y-axis label
    # ax2.yaxis.set_label_coords(1.05, 0.5)

    # Create a second right y-axis for total_signal_path_induced_time_drift
    ax3 = ax1.twinx()
    ax3.spines["right"].set_position(("outward", 60))  # Move the second right y-axis outward
    ax3.plot(t, (10**9) * total_signal_path_induced_time_drift, label=r'Signal Path Induced Time Drift', color='g', linestyle='-.', linewidth=1.5)
    ax3.set_ylabel(r'Signal Path Induced Time Drift [ns]', color='g')
    ax3.tick_params(axis='y', labelcolor='g')

    # # Adjust the position of the second right y-axis label
    # ax3.yaxis.set_label_coords(1.15, 0.5)

    # Adjust scientific notation offset text positions
    # ax1.get_yaxis().get_offset_text().set_x(-0.15)  # Move the offset text of the left y-axis
    ax2.get_yaxis().get_offset_text().set_x(1.05)   # Move the offset text of the first right y-axis
    ax3.get_yaxis().get_offset_text().set_x(1.15)   # Move the offset text of the second right y-axis

    # Title and legends
    plt.title('Time Drift: Real, Apparent, and Signal Path Induced vs. Time')
    fig.tight_layout()

    # Add legends
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax3.legend(loc='lower right')

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A33_time_drift_plot_in_nanosecond.png')

    # Show the plot
    plt.show()

    # User input
    time_duration_of_the_recorded_flight_test = 50 # s



    #----------------------------------------------------------------------------------------------
    # User input
    # The range of time that you want data plotted
    user_defined_t_initial = 1 # sec
    user_defined_t_final = 50 # sec
    # user_defined_t_initial = 0 # sec
    # user_defined_t_final = 50 # sec

    # Find the indices where values in t range from 15 to 20
    t_index_range = np.where((t >= user_defined_t_initial) & (t <= user_defined_t_final))[0]
    # If there are valid indices, find the first and last
    if t_index_range.size > 0:
        first_t_index = t_index_range[0]
        last_t_index = t_index_range[-1]
        print(f"[{first_t_index}:{last_t_index}]")
    #----------------------------------------------------------------------------------------------

    # Unit: nanosecond / day (meaning at this potential time drift rate)
    # Plotting the data
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot total_real_time_drift on the left y-axis
    ax1.plot(t[first_t_index:last_t_index], (86400/t[first_t_index:last_t_index]) * (10**9) * total_real_time_drift[first_t_index:last_t_index], label=r'Total Real Time Drift', color='b', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'Total Real Time Drift [ns/day]', color='b')
    ax1.set_ylim([0, 0.20])
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True)

    # # Adjust the position of the left y-axis label
    # ax1.yaxis.set_label_coords(-0.1, 0.5)

    # Create the first right y-axis for total_apparent_time_drift
    ax2 = ax1.twinx()
    ax2.plot(t[first_t_index:last_t_index], (86400/t[first_t_index:last_t_index]) * (10**6) * total_apparent_time_drift[first_t_index:last_t_index], label=r'Total Apparent Time Drift', color='r', linestyle='--', linewidth=1.5)
    ax2.plot(t[first_t_index:last_t_index], (86400/t[first_t_index:last_t_index]) * (10**6) * combined_time_drift[first_t_index:last_t_index], label=r'Combined Time Drift', color='k', linestyle=':', linewidth=1.5)
    ax2.set_ylabel(r'Total Apparent Time Drift [$\mu$s/day]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    # # Adjust the position of the first right y-axis label
    # ax2.yaxis.set_label_coords(1.05, 0.5)

    # Create a second right y-axis for total_signal_path_induced_time_drift
    ax3 = ax1.twinx()
    ax3.spines["right"].set_position(("outward", 60))  # Move the second right y-axis outward
    ax3.plot(t[first_t_index:last_t_index], (86400/t[first_t_index:last_t_index]) * (10**6) * total_signal_path_induced_time_drift[first_t_index:last_t_index], label=r'Signal Path Induced Time Drift', color='g', linestyle='-.', linewidth=1.5)
    ax3.set_ylabel(r'Signal Path Induced Time Drift [$\mu$s/day]', color='g')
    ax3.tick_params(axis='y', labelcolor='g')

    # # Adjust the position of the second right y-axis label
    # ax3.yaxis.set_label_coords(1.15, 0.5)

    # Adjust scientific notation offset text positions
    # ax1.get_yaxis().get_offset_text().set_x(-0.15)  # Move the offset text of the left y-axis
    ax2.get_yaxis().get_offset_text().set_x(1.05)   # Move the offset text of the first right y-axis
    ax3.get_yaxis().get_offset_text().set_x(1.15)   # Move the offset text of the second right y-axis

    # Title and legends
    plt.title('Time Drift: Real, Apparent, and Signal Path Induced vs. Time')
    fig.tight_layout()

    # Add legends
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax3.legend(loc='lower right')

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A34_time_drift_plot_in_ns_ms_per_day.png')

    # Show the plot
    plt.show()


    # Unit: m / day (time equivalent distance) (meaning at this potential time or time-equiv distance drift rate)
    # Plotting the data
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot total_real_time_drift on the left y-axis
    ax1.plot(t[first_t_index:last_t_index], (3 * (10**8)) * (86400/t[first_t_index:last_t_index]) * total_real_time_drift[first_t_index:last_t_index], label=r'Total Real Distance Drift', color='b', linewidth=1.5)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'Total Real Distance Drift [m/day]', color='b')
    ax1.set_ylim([0, 0.06])
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True)

    # # Adjust the position of the left y-axis label
    # ax1.yaxis.set_label_coords(-0.1, 0.5)

    # Create the first right y-axis for total_apparent_time_drift
    ax2 = ax1.twinx()
    ax2.plot(t[first_t_index:last_t_index], (10**(-3)) * (3 * (10**8)) * (86400/t[first_t_index:last_t_index]) * total_apparent_time_drift[first_t_index:last_t_index], label=r'Total Apparent Distance Drift', color='r', linestyle='--', linewidth=1.5)
    ax2.plot(t[first_t_index:last_t_index], (10**(-3)) * (3 * (10**8)) * (86400/t[first_t_index:last_t_index]) * combined_time_drift[first_t_index:last_t_index], label=r'Combined Distance Drift', color='k', linestyle=':', linewidth=1.5)
    ax2.set_ylabel(r'Total Apparent Distance Drift [km/day]', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    # # Adjust the position of the first right y-axis label
    # ax2.yaxis.set_label_coords(1.05, 0.5)

    # Create a second right y-axis for total_signal_path_induced_time_drift
    ax3 = ax1.twinx()
    ax3.spines["right"].set_position(("outward", 60))  # Move the second right y-axis outward
    ax3.plot(t[first_t_index:last_t_index], (10**(-3)) * (3 * (10**8)) * (86400/t[first_t_index:last_t_index]) * total_signal_path_induced_time_drift[first_t_index:last_t_index], label=r'Signal Path Induced Distance Drift', color='g', linestyle='-.', linewidth=1.5)
    ax3.set_ylabel(r'Signal Path Induced Distance Drift [km/day]', color='g')
    ax3.tick_params(axis='y', labelcolor='g')

    # # Adjust the position of the second right y-axis label
    # ax3.yaxis.set_label_coords(1.15, 0.5)

    # Adjust scientific notation offset text positions
    # ax1.get_yaxis().get_offset_text().set_x(-0.15)  # Move the offset text of the left y-axis
    ax2.get_yaxis().get_offset_text().set_x(1.05)   # Move the offset text of the first right y-axis
    ax3.get_yaxis().get_offset_text().set_x(1.15)   # Move the offset text of the second right y-axis

    # Title and legends
    plt.title('Estimated Distance Drift: Real, Apparent, and Signal Path Induced vs. Time')
    fig.tight_layout()

    # Add legends
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax3.legend(loc='lower right')

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A35_time_drift_plot_in_m_km_per_day.png')

    # Show the plot
    plt.show()


    # # Plotting the data
    # fig, ax1 = plt.subplots(figsize=(10, 6))

    # # Plot total_real_time_drift on the left y-axis
    # ax1.plot(t, total_real_time_drift, label=r'Total Real Time Drift', color='b', linewidth=1.5)
    # ax1.set_xlabel('Time [s]')
    # ax1.set_ylabel(r'Total Real Time Drift [s]', color='b')
    # ax1.tick_params(axis='y', labelcolor='b')
    # ax1.grid(True)

    # # Create the first right y-axis for total_apparent_time_drift
    # ax2 = ax1.twinx()
    # ax2.plot(t, total_apparent_time_drift, label=r'Total Apparent Time Drift', color='r', linestyle='--', linewidth=1.5)
    # ax2.set_ylabel(r'Total Apparent Time Drift [s]', color='r')
    # ax2.tick_params(axis='y', labelcolor='r')

    # # Create a second right y-axis for total_signal_path_induced_time_drift
    # ax3 = ax1.twinx()
    # ax3.spines["right"].set_position(("outward", 60))  # Move the second right y-axis outward
    # ax3.plot(t, total_signal_path_induced_time_drift, label=r'Signal Path Induced Time Drift', color='g', linestyle='-.', linewidth=1.5)
    # ax3.set_ylabel(r'Signal Path Induced Time Drift [s]', color='g')
    # ax3.tick_params(axis='y', labelcolor='g')

    # # Title and legends
    # plt.title('Time Drift: Real, Apparent, and Signal Path Induced vs. Time')
    # fig.tight_layout()

    # # Add legends
    # ax1.legend(loc='upper left')
    # ax2.legend(loc='upper right')
    # ax3.legend(loc='lower right')

    # # Save the figure
    # plt.savefig(xplane11_path + 'richard_checked_terms\\time_drift_plot.png')

    # # Show the plot
    # plt.show()

################################################################################################################################################################################
################################################################################################################################################################################
if whether_to_plot_Xg == 1:

    # Create the figure and plot u_g, v_g, and w_g vs time
    plt.figure(figsize=(10, 6))

    plt.plot(t, u_g, label=r'$u_g$ (Longitudinal Gust)', color=color_u_g, linewidth=1.5)
    plt.plot(t, v_g, label=r'$v_g$ (Lateral Gust)', color=color_v_g, linestyle='--', linewidth=1.5)
    plt.plot(t, w_g, label=r'$w_g$ (Vertical Gust)', color=color_w_g, linestyle='-.', linewidth=1.5)

    # Add labels and title
    plt.xlabel('Time [s]')
    plt.ylabel('Gust Components [m/s]')
    plt.title('Plot of Longitudinal, Lateral, and Vertical Gusts vs. Time')

    # Add grid and legend
    plt.grid(True)
    plt.legend()

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A36_gust_components_plot.png')

    # Show the plot
    plt.show()

    # Create the figure and plot p_g, q_g, and r_g vs time
    plt.figure(figsize=(10, 6))

    plt.plot(t, p_g, label=r'$p_g$ (Roll Rate Gust)', color=color_p_g, linewidth=1.5)
    plt.plot(t, q_g, label=r'$q_g$ (Pitch Rate Gust)', color=color_q_g, linestyle='--', linewidth=1.5)
    plt.plot(t, r_g, label=r'$r_g$ (Yaw Rate Gust)', color=color_r_g, linestyle='-.', linewidth=1.5)

    # Add labels and title
    plt.xlabel('Time [s]')
    plt.ylabel('Rotational Gust Components [rad/s]')
    plt.title('Plot of Roll, Pitch, and Yaw Rate Gusts vs. Time')

    # Add grid and legend
    plt.grid(True)
    plt.legend()

    # Save the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A37_rotational_gust_components_plot.png')

    # Show the plot
    plt.show()


    # Fourier transform
    def fourier_transform(gust_component, t):
        gust_ft = np.fft.fft(gust_component)
        freqs = np.fft.fftfreq(len(t), d=t[1] - t[0])
        return freqs, gust_ft

    # Perform Fourier transforms
    freq_u, u_g_ft = fourier_transform(u_g, t)
    freq_v, v_g_ft = fourier_transform(v_g, t)
    freq_w, w_g_ft = fourier_transform(w_g, t)
    freq_p, p_g_ft = fourier_transform(p_g, t)
    freq_q, q_g_ft = fourier_transform(q_g, t)
    freq_r, r_g_ft = fourier_transform(r_g, t)

    # # Frequency range (logarithmic scale for omega)
    # omega = np.logspace(-5, 5, 10000)

    # Plotting Fourier transforms of u_g, v_g, and w_g
    plt.figure(figsize=(10, 6))

    plt.plot(np.abs(freq_u), np.abs(u_g_ft), label=r'$|u_g(\omega)|$', linewidth=1.5)
    plt.plot(np.abs(freq_v), np.abs(v_g_ft), label=r'$|v_g(\omega)|$', linestyle='--', linewidth=1.5)
    plt.plot(np.abs(freq_w), np.abs(w_g_ft), label=r'$|w_g(\omega)|$', linestyle='-.', linewidth=1.5)

    # Add labels and title
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Frequency $\omega$ [rad/s]')
    plt.ylabel('Magnitude')
    plt.title('Fourier Transforms of Longitudinal, Lateral, and Vertical Gusts vs. Frequency')

    # Add grid and legend
    plt.grid(True)
    plt.legend()

    # Save and show the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A38_gust_components_fourier_plot.png')
    plt.show()

    # Plotting Fourier transforms of p_g, q_g, and r_g
    plt.figure(figsize=(10, 6))

    plt.plot(np.abs(freq_p), np.abs(p_g_ft), label=r'$|p_g(\omega)|$', linewidth=1.5)
    plt.plot(np.abs(freq_q), np.abs(q_g_ft), label=r'$|q_g(\omega)|$', linestyle='--', linewidth=1.5)
    plt.plot(np.abs(freq_r), np.abs(r_g_ft), label=r'$|r_g(\omega)|$', linestyle='-.', linewidth=1.5)

    # Add labels and title
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Frequency $\omega$ [rad/s]')
    plt.ylabel('Magnitude')
    plt.title('Fourier Transforms of Roll, Pitch, and Yaw Rate Gusts vs. Frequency')

    # Add grid and legend
    plt.grid(True)
    plt.legend()

    # Save and show the figure
    plt.savefig(xplane11_path + 'richard_checked_terms\\A39_rotational_gust_components_fourier_plot.png')
    plt.show()



################################################################################################################################################################################
################################################################################################################################################################################
if whether_to_plot_stage4_combined_history == 1:

    # Plot 1: Alpha averages of left and right wings
    plt.figure(figsize=(10, 6))
    plt.plot(t, alpha_average_of_left_wing, label=r'$\alpha_{left\_wing}$', color='b', linewidth=1.5)
    plt.plot(t, alpha_average_of_right_wing, label=r'$\alpha_{right\_wing}$', color='r', linestyle='--', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\alpha$ [deg]')
    plt.title('Plot 1: Angle of Attack Averages of Left and Right Wings vs. Time')
    plt.grid(True)
    plt.legend()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A41_alpha_of_each_wing_plot.png')
    plt.show()

    # Plot 2: CL averages of left and right wings
    plt.figure(figsize=(10, 6))
    plt.plot(t, CL_average_of_left_wing, label=r'$C_{L_{left\_wing}}$', color='b', linewidth=1.5)
    plt.plot(t, CL_average_of_right_wing, label=r'$C_{L_{right\_wing}}$', color='r', linestyle='--', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$C_L$')
    plt.title('Plot 2: Lift Coefficient Averages of Left and Right Wings vs. Time')
    plt.grid(True)
    plt.legend()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A42_CL_of_each_wing_plot.png')
    plt.show()

    # Plot 3: CD averages of left and right wings
    plt.figure(figsize=(10, 6))
    plt.plot(t, CD_average_of_left_wing, label=r'$C_{D_{left\_wing}}$', color='b', linewidth=1.5)
    plt.plot(t, CD_average_of_right_wing, label=r'$C_{D_{right\_wing}}$', color='r', linestyle='--', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$C_D$')
    plt.title('Plot 3: Drag Coefficient Averages of Left and Right Wings vs. Time')
    plt.grid(True)
    plt.legend()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A43_CD_of_each_wing_plot.png')
    plt.show()

    # Plot 4: Cl and Cn increments due to autorotation
    plt.figure(figsize=(10, 6))
    plt.plot(t, Cl_increment_due_to_autorotation, label=r'$\Delta C_l$ (Autorotation)', color='b', linewidth=1.5)
    plt.plot(t, Cn_increment_due_to_autorotation, label=r'$\Delta C_n$ (Autorotation)', color='r', linestyle='--', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.ylabel(r'Increments of $C_l$, $C_n$')
    plt.title('Plot 4: Cl and Cn Increments Due to Autorotation vs. Time')
    plt.grid(True)
    plt.legend()
    plt.savefig(xplane11_path + 'richard_checked_terms\\A44_autorotation_Cl_Cn_increments_plot.png')
    plt.show()

################################################################################################################################################################################
################################################################################################################################################################################