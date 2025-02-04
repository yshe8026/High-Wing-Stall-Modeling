import numpy as np
from scipy.signal import lfilter, bilinear

# Stall buffet effect developed through R33, R34, R35 # The output is more like white noise, but our simulation time step is 0.03 s and it won't make a big difference
def M11_calculate_az_ay_from_buffet_simplified(X_flow_sep):
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
    duration = 1  # Duration of the signal in seconds

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

    # Apply the transfer functions to the white noise signal % This is likely what makes the sim very slow since it generates time series every time step
    y0_lat = lfilter(b0_lat, a0_lat, u_white_noise)
    y1_lat = lfilter(b1_lat, a1_lat, u_white_noise)
    y_vert = lfilter(b_vert, a_vert, u_white_noise)

    # Combine the responses for lateral acceleration
    y_total_lat = y0_lat + y1_lat

    # Apply gain factors based on Delft paper
    lateral_gain_factor = np.sqrt(0.0231/0.0150)  # Delft, adjust as needed
    vertical_gain_factor = np.sqrt(0.205/0.173)   # Delft, adjust as needed

    # data_correction_factor = 10  # use this to correlate with Delft paper # Maybe this should be just "1"
    data_correction_factor = 1  # use this to correlate with Delft paper # Maybe this should be just "1"

    # Calculate Ay_noise and Az_noise for the entire time series
    Ay_noise_series = data_correction_factor * lateral_gain_factor * 10 * (1 - X_flow_sep) * y_total_lat
    Az_noise_series = data_correction_factor * vertical_gain_factor * 10 * (1 - X_flow_sep) * y_vert

    # Return the most recent value of Ay_noise and Az_noise
    ay_noise = Ay_noise_series[-1]
    az_noise = Az_noise_series[-1]

    return az_noise, ay_noise

