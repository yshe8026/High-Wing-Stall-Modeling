import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Parameters for the transfer function
sigma_u = 1.0  # Example turbulence intensity, adjust as needed
Lu = 100.0     # Turbulence scale length (longitudinal), adjust as needed
V = 50.0       # Velocity in m/s, adjust as needed
b = 8.10       # m
# Parameters for vertical component
sigma_v = 1.0  # Example turbulence intensity, adjust as needed
Lv = 80.0      # Turbulence scale length (vertical), adjust as needed
# Parameters for vertical component
sigma_w = 1.0  # Example turbulence intensity, adjust as needed
Lw = 80.0      # Turbulence scale length (vertical), adjust as needed


# Frequency range for plotting the Bode plot
w = np.logspace(-2, 2, 100)


# Define the transfer function for longitudinal H_u(s)
numerator_u = [sigma_u * np.sqrt(2 * Lu / (np.pi * V))]
denominator_u = [1, Lu / V]

# Create the transfer function using scipy
Hu = signal.TransferFunction(numerator_u, denominator_u)


# Define the transfer function for vertical H_w(s)
numerator_v = [sigma_v * np.sqrt(Lv / (np.pi * V)), (sigma_v * np.sqrt(Lv / (np.pi * V))) * (np.sqrt(3) * Lv / V)]
denominator_v = [1, 2 * Lv / V, (Lv / V)**2]

Hv = signal.TransferFunction(numerator_v, denominator_v)


# Define the transfer function for vertical H_w(s)
numerator_w = [sigma_w * np.sqrt(Lw / (np.pi * V)), (sigma_w * np.sqrt(Lw / (np.pi * V))) * (np.sqrt(3) * Lw / V)]
denominator_w = [1, 2 * Lw / V, (Lw / V)**2]

Hw = signal.TransferFunction(numerator_w, denominator_w)

# Plot the Bode plot for H_u(s)
w, mag, phase = signal.bode(Hu, w)

plt.figure()
plt.subplot(2, 1, 1)
plt.semilogx(w, mag)
plt.title(r'Bode plot of $H_u(s)$ (Longitudinal)')
plt.ylabel('Magnitude [dB]')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.semilogx(w, phase)
plt.ylabel(r'Phase [$^{\circ}$]')
plt.xlabel('Frequency [rad/s]')
plt.grid(True)

plt.tight_layout()
# plt.show()

# Plot the Bode plot for H_w(s)
w, mag, phase = signal.bode(Hv, w)

plt.figure()
plt.subplot(2, 1, 1)
plt.semilogx(w, mag)
plt.title(r'Bode plot of $H_v(s)$ (Lateral)')
plt.ylabel('Magnitude [dB]')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.semilogx(w, phase)
plt.ylabel(r'Phase [$^{\circ}$]')
plt.xlabel('Frequency [rad/s]')
plt.grid(True)

plt.tight_layout()
# plt.show()

# Plot the Bode plot for H_w(s)
w, mag, phase = signal.bode(Hw, w)

plt.figure()
plt.subplot(2, 1, 1)
plt.semilogx(w, mag)
plt.title(r'Bode plot of $H_w(s)$ (Vertical)')
plt.ylabel('Magnitude [dB]')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.semilogx(w, phase)
plt.ylabel(r'Phase [$^{\circ}$]')
plt.xlabel('Frequency [rad/s]')
plt.grid(True)

plt.tight_layout()
# plt.show()



#----------------------------------------------------------------------------------------------------------------------------------

# Define H_p(s) transfer function
numerator_hp = [(sigma_w * np.sqrt(0.8 / V) * ((np.pi/(4*b))**(1/6))) / (Lw**(1/3))]
denominator_hp = [1, (4 * b) / (np.pi * V)]

Hp = signal.TransferFunction(numerator_hp, denominator_hp)

# Define H_q(s) transfer function (H_q(s) uses H_w(s))
numerator_hq = np.polymul([0, 1/V], numerator_w) # Represents +s/V
denominator_hq = np.polymul([1, (4 * b) / (np.pi * V)], denominator_w)

Hq = signal.TransferFunction(numerator_hq, denominator_hq)

# Define H_r(s) transfer function (H_r(s) uses H_v(s))
numerator_hr = np.polymul([0, 1/V], numerator_v) # Represents +s/V
denominator_hr = np.polymul([1, (3 * b) / (np.pi * V)], denominator_v)

Hr = signal.TransferFunction(numerator_hr, denominator_hr)

# Frequency range for plotting the Bode plot
w = np.logspace(-2, 2, 100)

# Plot Bode plot for H_p(s)
w, mag, phase = signal.bode(Hp, w)
plt.figure()
plt.subplot(2, 1, 1)
plt.semilogx(w, mag)
plt.title(r'Bode plot of $H_p(s)$ (Longitudinal)')
plt.ylabel('Magnitude [dB]')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.semilogx(w, phase)
plt.ylabel(r'Phase [$^{\circ}$]')
plt.xlabel('Frequency [rad/s]')
plt.grid(True)

plt.tight_layout()

# Plot Bode plot for H_q(s)
w, mag, phase = signal.bode(Hq, w)
plt.figure()
plt.subplot(2, 1, 1)
plt.semilogx(w, mag)
plt.title(r'Bode plot of $H_q(s)$ (Vertical)')
plt.ylabel('Magnitude [dB]')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.semilogx(w, phase)
plt.ylabel(r'Phase [$^{\circ}$]')
plt.xlabel('Frequency [rad/s]')
plt.grid(True)

plt.tight_layout()

# Plot Bode plot for H_r(s)
w, mag, phase = signal.bode(Hr, w)
plt.figure()
plt.subplot(2, 1, 1)
plt.semilogx(w, mag)
plt.title(r'Bode plot of $H_r(s)$ (Lateral)')
plt.ylabel('Magnitude [dB]')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.semilogx(w, phase)
plt.ylabel(r'Phase [$^{\circ}$]')
plt.xlabel('Frequency [rad/s]')
plt.grid(True)

plt.tight_layout()

plt.show()