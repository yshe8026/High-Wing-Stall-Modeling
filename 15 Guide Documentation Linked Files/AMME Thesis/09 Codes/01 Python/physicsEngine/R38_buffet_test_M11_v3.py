import numpy as np
from scipy.signal import lfilter, bilinear, savgol_filter
import time

class BuffetNoiseGenerator:
    def __init__(self, fs=100, series_duration=0.5):
        self.fs = fs  # Sampling frequency in Hz
        self.series_duration = series_duration  # Duration of the noise series in seconds
        self.series_length = int(self.series_duration * self.fs)  # Number of samples in each series
        self.current_index = 0  # Initialize index to start of series
        self.Az_noise_series = np.zeros(self.series_length)
        self.Ay_noise_series = np.zeros(self.series_length)
        self.seed_counter = int(time.time())  # Initialize seed counter based on current time

    def M11_calculate_az_ay_from_buffet(self, X_flow_sep, dt):
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

        # Define the transfer function
        def transfer_function(H, w, Q, fs):
            b = [H * w**2]
            a = [1, w / Q, w**2]
            b_discrete, a_discrete = bilinear(b, a, fs=fs)
            return b_discrete, a_discrete

        # Generate a new noise series if we've reached the end of the current one
        if self.current_index >= self.series_length:
            # Update the seed for randomness
            np.random.seed(self.seed_counter)
            self.seed_counter += 1  # Increment seed counter to ensure different seed every time

            # Generate white noise time series
            t = np.arange(0, self.series_duration, 1/self.fs)
            u_white_noise = np.random.normal(0, 1, len(t))

            # Discretize the transfer functions for lateral and vertical acceleration
            b0_lat, a0_lat = transfer_function(H0_lateral, w0_lateral, Q0_lateral, self.fs)
            b1_lat, a1_lat = transfer_function(H1_lateral, w1_lateral, Q1_lateral, self.fs)
            b_vert, a_vert = transfer_function(H0_vertical, w0_vertical, Q0_vertical, self.fs)

            # Apply the transfer functions to the white noise signal
            y0_lat = lfilter(b0_lat, a0_lat, u_white_noise)
            y1_lat = lfilter(b1_lat, a1_lat, u_white_noise)
            y_vert = lfilter(b_vert, a_vert, u_white_noise)

            # Combine the responses for lateral acceleration
            y_total_lat = y0_lat + y1_lat

            # Apply gain factors based on Delft paper
            lateral_gain_factor = np.sqrt(0.0231/0.0150)  # Delft, adjust as needed
            vertical_gain_factor = np.sqrt(0.205/0.173)   # Delft, adjust as needed

            data_correction_factor = 10  # use this to correlate with Delft paper (this is revised to compensate for short generation duration of 0.5 sec)

            # Create the output vector based on the condition (this line might be redundant after all)
            whether_buffet_is_on = np.where(X_flow_sep < 0.89, 1, 0)

            # Calculate Ay_noise and Az_noise for the entire time series
            self.Ay_noise_series = data_correction_factor * lateral_gain_factor * 10 * (1 - X_flow_sep) * y_total_lat * whether_buffet_is_on
            self.Az_noise_series = data_correction_factor * vertical_gain_factor * 10 * (1 - X_flow_sep) * y_vert * whether_buffet_is_on

            # # Apply smoothing if needed (optional)
            # self.Ay_noise_series = savgol_filter(self.Ay_noise_series, window_length=11, polyorder=2)
            # self.Az_noise_series = savgol_filter(self.Az_noise_series, window_length=11, polyorder=2)

            # Reset index
            self.current_index = 0

        # Return the current values and increment the index
        az_noise = self.Az_noise_series[self.current_index]
        ay_noise = self.Ay_noise_series[self.current_index]
        self.current_index += int(dt * self.fs)  # Increment index based on time step

        return az_noise, ay_noise
    
import matplotlib.pyplot as plt

# Initialize the noise generator
buffet_noise_generator = BuffetNoiseGenerator(fs=100, series_duration=0.5)

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
    dt = 1 / fs  # Assume constant time step for simplicity
    az_noise, ay_noise = buffet_noise_generator.M11_calculate_az_ay_from_buffet(X[i], dt)
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

