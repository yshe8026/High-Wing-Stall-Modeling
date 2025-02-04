import numpy as np
import matplotlib.pyplot as plt

# A2 Q1(c) Determine the Probability Density Functions for each turbulence component
# Estimate the maximum gust magnitudes for each gust component

# Function: Get trim condition
def trim_condition():
    # Replace this with the appropriate way to load your trim conditions
    # For example, you could load them from a file or define them directly
    X_trim = np.array([4.28509979e+01, -3.34313946e-01, 1.16436917e+00, -1.69348740e-03, -2.52844822e-02, 2.04244591e-02, -3.38758335e-02, -1.34899263e-01, 2.92678634e+00, 0.00000000e+00, 0.00000000e+00, -6.40350647e+02])
    # X_trim[11] = -60.9
    X_trim[11] = -200
    U_trim = np.array([0.5, 0, 0, 0, 0])
    return X_trim, U_trim

# Function: Get Turbulent Scales & Intensities
def turb_scale_inten(X0):
    u_20_in_kn = 30  # Reference mean wind speed in knots (moderate intensity)
    u_20 = kn2mps(u_20_in_kn)  # Convert to m/s

    h = -X0[11]  # Trimmed flight altitude
    L_u = h / (0.177 + 0.0027 * h)**1.2
    L_v = h / (0.177 + 0.0027 * h)**1.2
    L_w = h

    sigma_u = (0.1 * u_20) * (1 / (0.177 + 0.0027 * h)**0.4)
    sigma_v = (0.1 * u_20) * (1 / (0.177 + 0.0027 * h)**0.4)
    sigma_w = 0.1 * u_20

    return L_u, L_v, L_w, sigma_u, sigma_v, sigma_w

# Function: Get aircraft wing span and trim velocity (b, V)
def geometry_and_trim():
    # FlightData = aero4560_LoadFlightData()  # Implement your aircraft data loading
    b = 8.10 # m
    # b = 10.2

    X0 = np.array([4.28509979e+01, -3.34313946e-01, 1.16436917e+00, -1.69348740e-03, -2.52844822e-02, 2.04244591e-02, -3.38758335e-02, -1.34899263e-01, 2.92678634e+00, 0.00000000e+00, 0.00000000e+00, -6.40350647e+02])
    V = np.sqrt(X0[0]**2 + X0[1]**2 + X0[2]**2)
    print(V)
    print(-X0[11])
    return b, V

# Function: Generate the Dryden Spectra PSDs
def dryden_psd(omega, V, b, L_u, L_v, L_w, sigma_u, sigma_v, sigma_w, plot_psd=0):
    Phi_u_g = sigma_u**2 * (2 * L_u / (np.pi * V)) * (1 / (1 + (L_u / V * omega)**2))
    Phi_v_g = sigma_v**2 * (L_v / (np.pi * V)) * (1 + 3 * (L_v / V * omega)**2) / ((1 + (L_v / V * omega)**2)**2)
    Phi_w_g = sigma_w**2 * (L_w / (np.pi * V)) * (1 + 3 * (L_w / V * omega)**2) / ((1 + (L_w / V * omega)**2)**2)
    Phi_p_g = 0.8 * sigma_w**2 / (L_w * V) * ((np.pi * L_w / (4 * b))**(1/3)) / (1 + ((4 * b / (np.pi * V)) * omega)**2)
    Phi_q_g = ((omega / V)**2 / (1 + (4 * b / (np.pi * V) * omega)**2)) * Phi_w_g
    Phi_r_g = ((omega / V)**2 / (1 + (3 * b / (np.pi * V) * omega)**2)) * Phi_v_g

    if plot_psd:
        plt.figure()
        plt.loglog(omega, Phi_u_g, label=r'$\Phi_{u_g}$')
        plt.loglog(omega, Phi_v_g, label=r'$\Phi_{v_g}$')
        plt.loglog(omega, Phi_w_g, label=r'$\Phi_{w_g}$')
        plt.loglog(omega, Phi_p_g, label=r'$\Phi_{p_g}$')
        plt.loglog(omega, Phi_q_g, label=r'$\Phi_{q_g}$')
        plt.loglog(omega, Phi_r_g, label=r'$\Phi_{r_g}$')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\Phi$')
        plt.grid(True)
        plt.legend()
        plt.show()

    return Phi_u_g, Phi_v_g, Phi_w_g, Phi_p_g, Phi_q_g, Phi_r_g

# Function: Calculate PDFs
def pdf(x, stan_dev_u, stan_dev_v, stan_dev_w, stan_dev_p, stan_dev_q, stan_dev_r, plot_pdf=0):
    pdf_u_g = (1 / (stan_dev_u * np.sqrt(2 * np.pi))) * np.exp(-(1/2) * (x / stan_dev_u)**2)
    pdf_v_g = (1 / (stan_dev_v * np.sqrt(2 * np.pi))) * np.exp(-(1/2) * (x / stan_dev_v)**2)
    pdf_w_g = (1 / (stan_dev_w * np.sqrt(2 * np.pi))) * np.exp(-(1/2) * (x / stan_dev_w)**2)
    pdf_p_g = (1 / (stan_dev_p * np.sqrt(2 * np.pi))) * np.exp(-(1/2) * (x / stan_dev_p)**2)
    pdf_q_g = (1 / (stan_dev_q * np.sqrt(2 * np.pi))) * np.exp(-(1/2) * (x / stan_dev_q)**2)
    pdf_r_g = (1 / (stan_dev_r * np.sqrt(2 * np.pi))) * np.exp(-(1/2) * (x / stan_dev_r)**2)

    if plot_pdf:
        plt.figure()
        plt.plot(x, pdf_u_g, label='$f_{u_g}(u_g)$')
        plt.plot(x, pdf_v_g, '--', label='$f_{v_g}(v_g)$')
        plt.plot(x, pdf_w_g, label='$f_{w_g}(w_g)$')
        plt.xlabel('$x$ [m/s]')
        plt.xlim([-5, 5])
        plt.ylabel('$f(x)$ [1/(m/s)]')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.figure()
        plt.plot(x, pdf_p_g, label='$f_{p_g}(p_g)$')
        plt.plot(x, pdf_q_g, label='$f_{q_g}(q_g)$')
        plt.plot(x, pdf_r_g, label='$f_{r_g}(r_g)$')
        plt.xlabel('$x$ [rad/s]')
        plt.xlim([-0.1, 0.1])
        plt.ylabel('$f(x)$ [1/(rad/s)]')
        plt.legend()
        plt.grid(True)
        plt.show()

    return pdf_u_g, pdf_v_g, pdf_w_g, pdf_p_g, pdf_q_g, pdf_r_g

# Function: Convert [kn] to [m/s]
def kn2mps(kn):
    return kn / 1.9438452

# Example workflow
X0, U0 = trim_condition()
L_u, L_v, L_w, sigma_u, sigma_v, sigma_w = turb_scale_inten(X0)
b, V = geometry_and_trim()
omega = np.logspace(-5, 5, 10000)
plot_psd = 1
Phi_u_g, Phi_v_g, Phi_w_g, Phi_p_g, Phi_q_g, Phi_r_g = dryden_psd(omega, V, b, L_u, L_v, L_w, sigma_u, sigma_v, sigma_w, plot_psd)

# Integrate Dryden Spectra PSDs from (-inf,inf) to get PDFs
stan_dev_u = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_u_g, omega))
stan_dev_v = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_v_g, omega))
stan_dev_w = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_w_g, omega))
stan_dev_p = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_p_g, omega))
stan_dev_q = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_q_g, omega))
stan_dev_r = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_r_g, omega))

# Maximum gust magnitudes with 3 sigma confidence
max_ug = 3 * stan_dev_u
max_vg = 3 * stan_dev_v
max_wg = 3 * stan_dev_w
max_pg = 3 * stan_dev_p
max_qg = 3 * stan_dev_q
max_rg = 3 * stan_dev_r

# Generate x vector for PDF f(x)
x = np.linspace(-10, 10, 10000)
plot_pdf = 1
pdf_u_g, pdf_v_g, pdf_w_g, pdf_p_g, pdf_q_g, pdf_r_g = pdf(x, stan_dev_u, stan_dev_v, stan_dev_w, stan_dev_p, stan_dev_q, stan_dev_r, plot_pdf)

print(f"------------------------------------------------------------------------------------------------------------------------------------------")
print(f"Standard Deviation of u_g: {stan_dev_u:.10e}")
print(f"Standard Deviation of v_g: {stan_dev_v:.10e}")
print(f"Standard Deviation of w_g: {stan_dev_w:.10e}")
print(f"Standard Deviation of p_g: {stan_dev_p:.10e}")
print(f"Standard Deviation of q_g: {stan_dev_q:.10e}")
print(f"Standard Deviation of r_g: {stan_dev_r:.10e}")
print(f"------------------------------------------------------------------------------------------------------------------------------------------")
