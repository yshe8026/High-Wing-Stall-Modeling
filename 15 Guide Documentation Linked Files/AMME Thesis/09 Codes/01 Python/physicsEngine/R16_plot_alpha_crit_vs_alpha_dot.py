# Plot alpha_crit vs alpha_dot as Dr. Lawson have suggested
import numpy as np
import matplotlib.pyplot as plt

# Define constants
a1 = (1.00) * 22.5  # Example value, adjust as needed
alpha_star = 20 * (np.pi / 180)  # Example value, adjust as needed
FlightData_Geometric_c = 0.99 # m
Vt = 40 # m/s
tau2 = (1.00) * 6.66 * (FlightData_Geometric_c / Vt)  # Example value, adjust as needed

# Define the relationship as a function
def alpha_cr(dot_alpha):
    return 0.8 * (np.arctanh(-0.9) / a1 + alpha_star + tau2 * dot_alpha)

# Generate data for plotting
dot_alpha_values = np.linspace(-4, 4, 400)  # Adjust range and number of points as needed
alpha_cr_values = alpha_cr(dot_alpha_values)

dot_alpha_values_in_deg = dot_alpha_values * (180/np.pi)
alpha_cr_values_in_deg = alpha_cr_values * (180/np.pi)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(dot_alpha_values_in_deg, alpha_cr_values_in_deg, label=r'$\alpha_{cr}$')
plt.xlabel(r'$\dot{\alpha}$', fontsize=14)
plt.ylabel(r'$\alpha_{cr}$', fontsize=14)
plt.title(r'$\alpha_{cr}$ vs $\dot{\alpha}$', fontsize=16)
plt.legend()
plt.grid(True)
plt.savefig('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\alpha_crit_vs_alpha_dot.png')  # Save as PNG
plt.show()
