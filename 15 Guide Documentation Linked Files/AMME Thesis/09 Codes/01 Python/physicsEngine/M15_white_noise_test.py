import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Generate white noise (with standard deviation of 0.1)
np.random.seed(42)
white_noise = np.random.normal(0, 0.1, 1000)

# Design a low-pass filter using a Butterworth filter
sampling_rate = 100  # Sampling rate in Hz (adjust as needed)
cutoff_freq = 10     # Cutoff frequency in Hz for the low-pass filter
order = 4            # Filter order

# Design the filter
b, a = signal.butter(order, cutoff_freq / (0.5 * sampling_rate), btype='low')

# Apply the filter to create band-limited noise
band_limited_noise = signal.filtfilt(b, a, white_noise)

# Plot the white noise and band-limited noise
plt.figure(figsize=(10, 6))
plt.plot(white_noise, label="White Noise")
plt.plot(band_limited_noise, label="Band-Limited Noise")
plt.title("Comparison of White Noise and Band-Limited Noise")
plt.legend()
plt.grid()
plt.show()
