import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


# from M11_calculate_az_ay_from_buffet import M11_calculate_az_ay_from_buffet


from scipy.signal import lfilter, bilinear
# Define the function as previously discussed
def M11_calculate_az_ay_from_buffet(X_flow_sep):
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
    fs = 50  # Hz
    duration = 0.5  # Duration of the signal in seconds

    # Define the transfer function
    def transfer_function(H, w, Q, fs):
        b = [H * w**2]
        a = [1, w / Q, w**2]
        b_discrete, a_discrete = bilinear(b, a, fs=fs)
        return b_discrete, a_discrete

    # Discretize the transfer functions for lateral acceleration
    b0_lat, a0_lat = transfer_function(H0_lateral, w0_lateral, Q0_lateral, fs)
    b1_lat, a1_lat = transfer_function(H1_lateral, w1_lateral, Q1_lateral, fs)

    # Discretize the transfer function for vertical acceleration
    b_vert, a_vert = transfer_function(H0_vertical, w0_vertical, Q0_vertical, fs)

    # Generate white noise time series for the duration specified
    t = np.arange(0, duration, 1/fs)
    u_white_noise = np.random.normal(0, 1, len(t))

    # Apply the transfer functions to the white noise signal
    y0_lat = lfilter(b0_lat, a0_lat, u_white_noise)
    y1_lat = lfilter(b1_lat, a1_lat, u_white_noise)
    y_vert = lfilter(b_vert, a_vert, u_white_noise)

    # Combine the responses for lateral acceleration
    y_total_lat = y0_lat + y1_lat

    # Apply gain factors based on Delft paper
    lateral_gain_factor = np.sqrt(0.0231/0.0150)  # Delft, adjust as needed
    vertical_gain_factor = np.sqrt(0.205/0.173)   # Delft, adjust as needed

    data_correction_factor = 10  # use this to correlate with Delft paper

    # Calculate Ay_noise and Az_noise for the entire time series
    Ay_noise_series = data_correction_factor * lateral_gain_factor * 10 * (1 - X_flow_sep) * y_total_lat
    Az_noise_series = data_correction_factor * vertical_gain_factor * 10 * (1 - X_flow_sep) * y_vert

    # # Apply Savitzky-Golay filter to smooth the noise series
    # Ay_noise_series = savgol_filter(Ay_noise_series, window_length=11, polyorder=2)
    # Az_noise_series = savgol_filter(Az_noise_series, window_length=11, polyorder=2)

    # Return the most recent value of Ay_noise and Az_noise
    ay_noise = Ay_noise_series[-1]
    az_noise = Az_noise_series[-1]

    return az_noise, ay_noise, 


# Sampling frequency and time vector
fs = 100  # Sampling frequency (Hz)
t = np.arange(0, 10, 1/fs)  # 10 seconds of data

# Initialize X with zeros
X = np.zeros_like(t)

# Define the segments for X
X[(t >= 0) & (t < 3.5)] = np.linspace(0.9, 0.86, num=len(t[(t >= 0) & (t < 3.5)]))
X[(t >= 3.5) & (t < 4.5)] = np.linspace(0.86, 0.75, num=len(t[(t >= 3.5) & (t < 4.5)]))
X[(t >= 4.5) & (t < 7.5)] = np.linspace(0.75, 0.96, num=len(t[(t >= 4.5) & (t < 7.5)]))
X[(t >= 7.5)] = 0.96
# Smoothing the X time series using Savitzky-Golay filter
window_length = 101  # Must be an odd number
polyorder = 3  # Polynomial order
X = savgol_filter(X, window_length, polyorder)

# Initialize arrays to store the noise values
Az_noise = np.zeros_like(t)
Ay_noise = np.zeros_like(t)

# Run the "flight simulation" and generate noise values at each time step
for i in range(len(t)):
    az_noise, ay_noise = M11_calculate_az_ay_from_buffet(X[i])
    Az_noise[i] = az_noise
    Ay_noise[i] = ay_noise

# Plot the Az_noise time series (Vertical Acceleration)
plt.figure(figsize=(12, 9))

plt.subplot(3, 1, 1)
plt.plot(t, Az_noise, label='Az Noise (Vertical)', color='b')
plt.title(r'Vertical Acceleration Time Series ($A_z$)')
plt.xlabel(r't [s]')
plt.xlim([0, 10])
plt.ylabel(r'$A_z$ [m/s$^{2}$]')
plt.ylim([-10, 10])
plt.legend()
plt.grid(True)

# Plot the Ay_noise time series (Lateral Acceleration)
plt.subplot(3, 1, 2)
plt.plot(t, Ay_noise, label='Ay Noise (Lateral)', color='r')
plt.title(r'Lateral Acceleration Time Series ($A_y$)')
plt.xlabel(r't [s]')
plt.xlim([0, 10])
plt.ylabel(r'$A_y$ [m/s$^{2}$]')
plt.ylim([-5, 5])
plt.legend()
plt.grid(True)

# Plot X
plt.subplot(3, 1, 3)
plt.plot(t, X, label='X')
# Add a horizontal line at X = 0.89
plt.axhline(y=0.89, color='k', linestyle='--', label='X = 0.89')
plt.xlabel(r't [s]')
plt.xlim([0, 10])
plt.ylabel(r'$X$ [-]')
plt.ylim([0.7, 1.0])
plt.title(r'Flow Separation Parameter Times Series $X$')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

