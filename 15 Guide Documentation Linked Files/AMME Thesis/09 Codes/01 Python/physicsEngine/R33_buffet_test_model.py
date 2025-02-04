import numpy as np
from scipy.signal import lfilter, freqz, bilinear
import matplotlib.pyplot as plt

# Given parameters
H0 = 0.02
w0 = 36.43  # rad/s
Q0 = 4.19

H1 = 0.01
w1 = 64.71  # rad/s
Q1 = 11.99

# Sampling frequency
fs = 100  # Hz # smaller this number is, the smoother the transfer function looks, but the worse the time response looks (Uncertainty Principle)

# Define the transfer function for G0(jw) and G1(jw)
def transfer_function(H, w, Q, fs):
    # Continuous transfer function coefficients
    b = [H * w**2]
    a = [1, w / Q, w**2]
    
    # Discretize using bilinear transform
    b_discrete, a_discrete = bilinear(b, a, fs=fs)
    return b_discrete, a_discrete

# Discretize the transfer functions
b0, a0 = transfer_function(H0, w0, Q0, fs)
b1, a1 = transfer_function(H1, w1, Q1, fs)

# Generate the frequency response of the filters
w, h0 = freqz(b0, a0, fs=fs)
_, h1 = freqz(b1, a1, fs=fs)


lateral_gain_factor = np.sqrt(0.0231/0.0150) # To correlate with real test flight data from Delft paper
h0 = h0 * lateral_gain_factor
h1 = h1 * lateral_gain_factor

# Vertical PSD
plt.figure(figsize=(12, 6))

plt.subplot(2, 1, 1)
plt.semilogx(w * (2*np.pi), 20 * np.log10(abs(h0+h1)), label='Lateral', color = 'r')
plt.title(r'Magnitude Response of the Filters')
plt.xlabel(r'$\omega$ [rad/s]')
plt.xlim([10, 200])
plt.ylabel(r'$H(j \omega)$ Magnitude (dB)')
plt.ylim([-60, 0])
plt.legend()
plt.grid(True, which="both", ls="--")

# Lateral PSD
plt.subplot(2, 1, 2)
plt.semilogx(w * (2*np.pi), np.angle(h0+h1)*(180/np.pi), label='Lateral', color = 'r')
plt.title(r'Phase Response of the Filters')
plt.xlabel(r'$\omega$ [rad/s]')
plt.xlim([10, 200])
plt.ylabel(r'$H(j \omega)$ Phase ($^{\circ}$)')
plt.ylim([-180, 0])
plt.legend()
plt.grid(True, which="both", ls="--")

plt.tight_layout()
plt.show()


# Plot magnitude response
plt.figure(figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.plot(w * (2*np.pi), abs(h0+h1)**2, label='G(jw) Magnitude Response')
plt.title(r'$a_{z}$ Power Spectral Density Function (PSD)')
plt.xlabel(r'$\omega$ [rad/s]')
plt.xlim([0, 157])
plt.ylabel(r'$S_{a_{z} a_{z}}$ [(m$^2$/$s^4$)/(rad/s)]')
plt.ylim([0, 0.025])
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(w * (2*np.pi), abs(h0+h1)**2, label='G(jw) Magnitude Response', color = 'r')
plt.title(r'$a_{y}$ Power Spectral Density Function (PSD)')
plt.xlabel(r'$\omega$ [rad/s]')
plt.xlim([0, 157])
plt.ylabel(r'$S_{a_{y} a_{y}}$ [(m$^2$/$s^4$)/(rad/s)]')
plt.ylim([0, 0.025])
plt.grid()
plt.tight_layout()
plt.show()


# Time vector
fs = 100  # Sampling frequency (Hz)
t = np.arange(0, 10, 1/fs)  # 10 seconds of data

# Generate white noise time series
u_white_noise = np.random.normal(0, 1, len(t))

# Define the transfer function for G0(jw) and G1(jw)
def transfer_function(H, w, Q, fs):
    # Continuous transfer function coefficients
    b = [H * w**2]
    a = [1, w / Q, w**2]
    
    # Discretize using bilinear transform
    b_discrete, a_discrete = bilinear(b, a, fs=fs)
    return b_discrete, a_discrete

# Discretize the transfer functions
from scipy.signal import bilinear

b0, a0 = transfer_function(H0, w0, Q0, fs)
b1, a1 = transfer_function(H1, w1, Q1, fs)

# Apply the transfer functions to the white noise signal
y0 = lfilter(b0, a0, u_white_noise)
y1 = lfilter(b1, a1, u_white_noise)

# Combine the responses
y_total = y0 + y1

# Define X and calculate Ay
X = 0.5  # Example value, this should be provided or calculated elsewhere in your code
Ay_noise = lateral_gain_factor * 10 * (1 - X) * y_total

# Output the Ay_noise time series
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(t, Ay_noise, label='Ay Noise')
plt.title('Lateral Acceleration Time Series (Ay)')
plt.xlabel('Time (s)')
plt.ylabel('Ay (Lateral Acceleration)')
plt.legend()
plt.show()
