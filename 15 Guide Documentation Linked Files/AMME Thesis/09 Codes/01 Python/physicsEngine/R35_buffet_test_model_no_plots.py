import numpy as np
from scipy.signal import lfilter, freqz, bilinear
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# Given parameters for lateral acceleration
H0_lateral = 0.02
w0_lateral = 36.43  # rad/s
Q0_lateral = 4.19

H1_lateral = 0.01
w1_lateral = 64.71  # rad/s
Q1_lateral = 11.99

# Given parameters for vertical acceleration (az)
H0_vertical = 0.05
w0_vertical = 75.92  # rad/s
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


# Sampling frequency and time vector
fs = 100  # Sampling frequency (Hz)
t = np.arange(0, 10, 1/fs)  # 10 seconds of data

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

data_correction_factor = 10 # use this to correlate with Delft paper

X_flow_sep = 0.7 # placeholder

Ay_noise = data_correction_factor * lateral_gain_factor * 10 * (1 - X_flow_sep) * y_total_lat
Az_noise = data_correction_factor * vertical_gain_factor * 10 * (1 - X_flow_sep) * y_vert

