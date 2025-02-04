import numpy as np
from scipy.signal import lfilter, bilinear, savgol_filter
import time

# Developed through R33 - R38
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

            # data_correction_factor = 10  # use this to correlate with Delft paper (this is revised to compensate for short generation duration of 0.5 sec) # Maybe this should be just "1"
            data_correction_factor = 1  # use this to correlate with Delft paper (this is revised to compensate for short generation duration of 0.5 sec) # Maybe this should be just "1"

            # Create the output vector based on the condition (this line is redundant, and only reserved for using R38 to test this module)
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