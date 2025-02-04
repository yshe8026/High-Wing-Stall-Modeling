import numpy as np
import matplotlib.pyplot as plt

# This is a concept that is not adopted due to its computational expensive nature

# Parameters
fs = 1000  # Sampling frequency (Hz)
T = 100    # Total time (seconds)
t = np.linspace(0, T, T * fs, endpoint=False)

# Define the frequency band
lowcut = 0.01  # Lower bound of the frequency band (Hz)
highcut = 0.1  # Upper bound of the frequency band (Hz)

# Number of sinusoids to approximate the noise
num_sinusoids = 1000

# Generate random frequencies within the band
frequencies = np.random.uniform(lowcut, highcut, num_sinusoids)

# Generate random phases and amplitudes
phases = np.random.uniform(0, 2 * np.pi, num_sinusoids)
amplitudes = np.random.normal(0, 1, num_sinusoids)

# Sum of sinusoids
u_g_bandlimited = np.zeros_like(t)
for freq, phase, amplitude in zip(frequencies, phases, amplitudes):
    u_g_bandlimited += amplitude * np.sin(2 * np.pi * freq * t + phase)

# Plot the generated time-domain signal
plt.figure(figsize=(10, 6))
plt.plot(t, u_g_bandlimited, label=f'Band-Limited Noise ({lowcut} Hz - {highcut} Hz)', color='orange')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title(f'Band-Limited Noise in Time Domain ({lowcut} Hz - {highcut} Hz)')
plt.grid(True)
plt.show()

# Fourier Transform to check frequency content
def plot_fourier_transform(signal, t, title):
    signal_ft = np.fft.fft(signal)
    freqs = np.fft.fftfreq(len(t), d=t[1] - t[0])

    plt.figure(figsize=(10, 6))
    plt.plot(np.abs(freqs), np.abs(signal_ft), label=title, linewidth=1.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Frequency $\omega$ [Hz]')
    plt.ylabel('Magnitude')
    plt.title(title)
    plt.grid(True)
    plt.show()

# Plot Fourier Transforms to verify the band-limited nature of the noise
plot_fourier_transform(u_g_bandlimited, t, f'Fourier Transform of Band-Limited Noise ({lowcut} Hz - {highcut} Hz)')
