import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Function: Generate the Dryden filters
def dryden_filters(omega, V, b, L_u, L_v, L_w, sigma_u, sigma_v, sigma_w, plot_filters=0):

    # Define the transfer function for longitudinal H_u(s)
    numerator_u = [sigma_u * np.sqrt(2 * L_u / (np.pi * V))]
    denominator_u = [1, L_u / V]

    # Create the transfer function using scipy
    Hu = signal.TransferFunction(numerator_u, denominator_u)


    # Define the transfer function for vertical H_w(s)
    numerator_v = [sigma_v * np.sqrt(L_v / (np.pi * V)), (sigma_v * np.sqrt(L_v / (np.pi * V))) * (np.sqrt(3) * L_v / V)]
    denominator_v = [1, 2 * L_v / V, (L_v / V)**2]

    Hv = signal.TransferFunction(numerator_v, denominator_v)


    # Define the transfer function for vertical H_w(s)
    numerator_w = [sigma_w * np.sqrt(L_w / (np.pi * V)), (sigma_w * np.sqrt(L_w / (np.pi * V))) * (np.sqrt(3) * L_w / V)]
    denominator_w = [1, 2 * L_w / V, (L_w / V)**2]

    Hw = signal.TransferFunction(numerator_w, denominator_w)


    # Define H_p(s) transfer function
    numerator_hp = [(sigma_w * np.sqrt(0.8 / V) * ((np.pi/(4*b))**(1/6))) / (L_w**(1/3))]
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

    if plot_filters:
        #----------------------------------------------------------------------------------------------------------------------------------
        # Plot the Bode plot for H_u(s)
        omega, mag, phase = signal.bode(Hu, omega)

        plt.figure()
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, mag)
        plt.title(r'Bode plot of $H_u(s)$')
        plt.ylabel('Magnitude [dB]')
        plt.grid(True)

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, phase)
        plt.ylabel(r'Phase [$^{\circ}$]')
        plt.xlabel('Frequency [rad/s]')
        plt.grid(True)

        plt.tight_layout()
        # plt.show()

        # Plot the Bode plot for H_w(s)
        omega, mag, phase = signal.bode(Hv, omega)

        plt.figure()
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, mag)
        plt.title(r'Bode plot of $H_v(s)$')
        plt.ylabel('Magnitude [dB]')
        plt.grid(True)

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, phase)
        plt.ylabel(r'Phase [$^{\circ}$]')
        plt.xlabel('Frequency [rad/s]')
        plt.grid(True)

        plt.tight_layout()
        # plt.show()

        # Plot the Bode plot for H_w(s)
        omega, mag, phase = signal.bode(Hw, omega)

        plt.figure()
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, mag)
        plt.title(r'Bode plot of $H_w(s)$')
        plt.ylabel('Magnitude [dB]')
        plt.grid(True)

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, phase)
        plt.ylabel(r'Phase [$^{\circ}$]')
        plt.xlabel('Frequency [rad/s]')
        plt.grid(True)

        plt.tight_layout()
        # plt.show()

        #----------------------------------------------------------------------------------------------------------------------------------

        # Plot Bode plot for H_p(s)
        omega, mag, phase = signal.bode(Hp, omega)
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, mag)
        plt.title(r'Bode plot of $H_p(s)$')
        plt.ylabel('Magnitude [dB]')
        plt.grid(True)

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, phase)
        plt.ylabel(r'Phase [$^{\circ}$]')
        plt.xlabel('Frequency [rad/s]')
        plt.grid(True)

        plt.tight_layout()

        # Plot Bode plot for H_q(s)
        omega, mag, phase = signal.bode(Hq, omega)
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, mag)
        plt.title(r'Bode plot of $H_q(s)$')
        plt.ylabel('Magnitude [dB]')
        plt.grid(True)

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, phase)
        plt.ylabel(r'Phase [$^{\circ}$]')
        plt.xlabel('Frequency [rad/s]')
        plt.grid(True)

        plt.tight_layout()

        # Plot Bode plot for H_r(s)
        omega, mag, phase = signal.bode(Hr, omega)
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, mag)
        plt.title(r'Bode plot of $H_r(s)$')
        plt.ylabel('Magnitude [dB]')
        plt.grid(True)

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, phase)
        plt.ylabel(r'Phase [$^{\circ}$]')
        plt.xlabel('Frequency [rad/s]')
        plt.grid(True)

        plt.tight_layout()

        plt.show()
        #----------------------------------------------------------------------------------------------------------------------------------

    return Hu, Hv, Hw, Hp, Hq, Hr

# Example use of the function
if __name__ == "__main__":
    # Parameters for the transfer function
    sigma_u = 1.0  # Example turbulence intensity, adjust as needed
    L_u = 100.0     # Turbulence scale length (longitudinal), adjust as needed
    V = 50.0       # Velocity in m/s, adjust as needed
    b = 8.10       # m
    # Parameters for vertical component
    sigma_v = 1.0  # Example turbulence intensity, adjust as needed
    L_v = 80.0      # Turbulence scale length (vertical), adjust as needed
    # Parameters for vertical component
    sigma_w = 1.0  # Example turbulence intensity, adjust as needed
    L_w = 80.0      # Turbulence scale length (vertical), adjust as needed

    # Frequency range for plotting the Bode plot
    omega = np.logspace(-2, 2, 100)

    Hu, Hv, Hw, Hp, Hq, Hr = dryden_filters(omega, V, b, L_u, L_v, L_w, sigma_u, sigma_v, sigma_w, plot_filters=1)













