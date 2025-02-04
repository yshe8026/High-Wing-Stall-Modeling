import numpy as np
from scipy.signal import lfilter, freqz, bilinear
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import os

'''Essentially similar as R34, but removed X_flow_sep term from the generated buffet tape. This X_flow_sep term will be incorporated when we run the XPlane simulation instead.'''

#==============================================================================================
# For extracting taped gust
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
#==============================================================================================

# St_J400 = 0.24
# St_cessna_citation_II = 0.25
# V_stall_J400 = 30.8667 # m/s
# V_stall_cessna_citation_II = 56.5889 # m/s
# c_bar_J400 = 0.99 # m/s
# c_bar_cessna_citation_II = 2.06 # m/s
# f_J400_to_f_cessna_citation_II_ratio = (St_J400 / St_cessna_citation_II) * (V_stall_J400 / V_stall_cessna_citation_II) * (c_bar_cessna_citation_II / c_bar_J400)

f_J400_to_f_cessna_citation_II_ratio = 1.08958774 # f_J400_to_f_cessna_citation_II_ratio = (St_J400 / St_cessna_citation_II) * (V_stall_J400 / V_stall_cessna_citation_II) * (c_bar_cessna_citation_II / c_bar_J400)

# Given parameters for lateral acceleration
H0_lateral = 0.02
w0_lateral = 36.43 * f_J400_to_f_cessna_citation_II_ratio # rad/s
Q0_lateral = 4.19

H1_lateral = 0.01
w1_lateral = 64.71 * f_J400_to_f_cessna_citation_II_ratio # rad/s
Q1_lateral = 11.99

# Given parameters for vertical acceleration (az)
H0_vertical = 0.05
w0_vertical = 75.92 * f_J400_to_f_cessna_citation_II_ratio # rad/s
Q0_vertical = 8.28

# Sampling frequency
fs = 100  # Hz

# Define the transfer function
def transfer_function(H, w, Q, fs):
    # Continuous transfer function coefficients
    b = [H * w**2]
    a = [1, w / Q, w**2]
    
    # Discretize using bilinear transform
    b_discrete, a_discrete = bilinear(b, a, fs=fs)
    return b_discrete, a_discrete

# Discretize the transfer functions for lateral acceleration
b0_lat, a0_lat = transfer_function(H0_lateral, w0_lateral, Q0_lateral, fs)
b1_lat, a1_lat = transfer_function(H1_lateral, w1_lateral, Q1_lateral, fs)

# Discretize the transfer function for vertical acceleration
b_vert, a_vert = transfer_function(H0_vertical, w0_vertical, Q0_vertical, fs)

# Generate the frequency response of the filters
w, h0_lat = freqz(b0_lat, a0_lat, fs=fs)
_, h1_lat = freqz(b1_lat, a1_lat, fs=fs)
_, h_vert = freqz(b_vert, a_vert, fs=fs)

# Apply lateral gain factor to match real test flight data
lateral_gain_factor = np.sqrt(0.0231/0.0150)  # Delft, adjust as needed
h0_lat *= lateral_gain_factor
h1_lat *= lateral_gain_factor

# Apply lateral gain factor to match real test flight data
vertical_gain_factor = np.sqrt(0.205/0.173)  # Delft, adjust as needed
h_vert *= vertical_gain_factor

# Plot the magnitude response with log scale for frequency axis
plt.figure(figsize=(12, 6))

plt.subplot(2, 1, 1)
plt.semilogx(w * (2*np.pi), 20 * np.log10(abs(h0_lat + h1_lat)), label='Lateral', color='r')
plt.semilogx(w * (2*np.pi), 20 * np.log10(abs(h_vert)), label='Vertical', color='b')
plt.title(r'Magnitude Response of the Filters')
plt.xlabel(r'$\omega$ [rad/s]')
plt.xlim([10, 200])
plt.ylabel(r'$H(j \omega)$ Magnitude (dB)')
plt.ylim([-60, 0])
plt.legend()
plt.grid(True, which="both", ls="--")

# Plot the phase response with log scale for frequency axis
plt.subplot(2, 1, 2)
plt.semilogx(w * (2*np.pi), np.angle(h0_lat + h1_lat)*(180/np.pi), label='Lateral', color='r')
plt.semilogx(w * (2*np.pi), np.angle(h_vert)*(180/np.pi), label='Vertical', color='b')
plt.title(r'Phase Response of the Filters')
plt.xlabel(r'$\omega$ [rad/s]')
plt.xlim([10, 200])
plt.ylabel(r'$H(j \omega)$ Phase ($^{\circ}$)')
plt.ylim([-180, 0])
plt.legend()
plt.grid(True, which="both", ls="--")

plt.tight_layout()
plt.show()

# Plot Power Spectral Density (PSD)
plt.figure(figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.plot(w * (2*np.pi), abs(h_vert)**2, label='Vertical (Az) PSD', color='b')
plt.title(r'$a_{z}$ Power Spectral Density Function (PSD)')
plt.xlabel(r'$\omega$ [rad/s]')
plt.xlim([0, 157])
plt.ylabel(r'$S_{a_{z} a_{z}}$ [(m$^2$/$s^4$)/(rad/s)]')
plt.ylim([0, 0.25])
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(w * (2*np.pi), abs(h0_lat + h1_lat)**2, label='Lateral (Ay) PSD', color='r')
plt.title(r'$a_{y}$ Power Spectral Density Function (PSD)')
plt.xlabel(r'$\omega$ [rad/s]')
plt.xlim([0, 157])
plt.ylabel(r'$S_{a_{y} a_{y}}$ [(m$^2$/$s^4$)/(rad/s)]')
plt.ylim([0, 0.025])
plt.grid()

plt.tight_layout()
plt.show()



# Sampling frequency and time vector
fs = 100  # Sampling frequency (Hz)
# t = np.arange(0, 10, 1/fs)  # 10 seconds of data
t = np.arange(0, 100, 1/fs)  # 100 seconds of data
# Generate white noise time series
u_white_noise = np.random.normal(0, 1, len(t))

# Discretize the transfer functions for lateral acceleration
b0_lat, a0_lat = transfer_function(H0_lateral, w0_lateral, Q0_lateral, fs)
b1_lat, a1_lat = transfer_function(H1_lateral, w1_lateral, Q1_lateral, fs)

# Discretize the transfer function for vertical acceleration
b_vert, a_vert = transfer_function(H0_vertical, w0_vertical, Q0_vertical, fs)

# Apply the transfer functions to the white noise signal
y0_lat = lfilter(b0_lat, a0_lat, u_white_noise)
y1_lat = lfilter(b1_lat, a1_lat, u_white_noise)
y_vert = lfilter(b_vert, a_vert, u_white_noise)

# Combine the responses for lateral acceleration
y_total_lat = y0_lat + y1_lat

# Define X and calculate Ay and Az
X = 0.7  # Example value, adjust as needed

# Initialize X with zeros
X = np.zeros_like(t)

# Define the segments
X[(t >= 0) & (t < 3.5)] = np.linspace(0.9, 0.86, num=len(t[(t >= 0) & (t < 3.5)]))
X[(t >= 3.5) & (t < 4.5)] = np.linspace(0.86, 0.75, num=len(t[(t >= 3.5) & (t < 4.5)]))
X[(t >= 4.5) & (t < 7.5)] = np.linspace(0.75, 0.96, num=len(t[(t >= 4.5) & (t < 7.5)]))
X[(t >= 7.5)] = 0.96
# Smoothing the X time series using Savitzky-Golay filter
window_length = 101  # Must be an odd number
polyorder = 3  # Polynomial order
X = savgol_filter(X, window_length, polyorder)

data_correction_factor = 10 # use this to correlate with Delft paper

# Create the output vector based on the condition
# whether_buffet_is_on = np.where(X < 0.89, 1, 0)
whether_buffet_is_on = 1
# Ay_noise = data_correction_factor * lateral_gain_factor * 10 * (1 - X) * y_total_lat * whether_buffet_is_on
# Az_noise = data_correction_factor * vertical_gain_factor * 10 * (1 - X) * y_vert * whether_buffet_is_on

Ay_noise = data_correction_factor * lateral_gain_factor * 10 * y_total_lat * whether_buffet_is_on #* (1 - X)
Az_noise = data_correction_factor * vertical_gain_factor * 10 * y_vert * whether_buffet_is_on #* (1 - X)

# Construct Xg matrix: first 6 columns with the time series, last 6 columns with zeros
Ay_Az_noises = np.zeros((len(t), 2))
Ay_Az_noises[:, 0] = Ay_noise
Ay_Az_noises[:, 1] = Az_noise

# Print or save Xg as needed
np.set_printoptions(threshold=np.inf, linewidth=500)
print("Ay_Az_noises:")
print(Ay_Az_noises)

#####################################################################################
# Extract Input
# User input
Ay_Az_noises_data_matrix_taped = Ay_Az_noises
t_taped_for_Ay_Az_noises = t

np.save(xplane11_path + 'Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\Ay_Az_noises_data_matrix_taped', Ay_Az_noises_data_matrix_taped) # Relative path approach that works on all machines
np.save(xplane11_path + 'Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\t_taped_for_Ay_Az_noises', t_taped_for_Ay_Az_noises) # Relative path approach that works on all machines

#####################################################################################

# Plot the Ay_noise time series (Lateral Acceleration)
plt.figure(figsize=(12, 9))

plt.subplot(3, 1, 1)
plt.plot(t, Az_noise, label='Az Noise (Vertical)', color='b')
plt.title(r'Vertical Acceleration Time Series ($A_z$)')
plt.xlabel(r't [s]')
plt.xlim([0, 100])
plt.ylabel(r'$A_z$ [m/s$^{2}$]')
plt.ylim([-10, 10])
plt.legend()
plt.grid(True)

# Plot the Az_noise time series (Vertical Acceleration)
plt.subplot(3, 1, 2)
plt.plot(t, Ay_noise, label='Ay Noise (Lateral)', color='r')
plt.title(r'Lateral Acceleration Time Series ($A_y$)')
plt.xlabel(r't [s]')
plt.xlim([0, 100])
plt.ylabel(r'$A_y$ [m/s$^{2}$]')
plt.ylim([-5, 5])
plt.legend()
plt.grid(True)

# Plot X
plt.subplot(3, 1, 3)
plt.plot(t, X)
# Add a horizontal line at X = 0.89
plt.axhline(y=0.89, color='k', linestyle='--', label='X = 0.89')
plt.xlabel(r't [s]')
plt.xlim([0, 100])
plt.ylabel(r'$X$ [-]')
plt.ylim([0.7, 1.0])
plt.title(r'Flow Separation Parameter Times Series $X$')
plt.grid(True)

plt.tight_layout()
plt.show()






