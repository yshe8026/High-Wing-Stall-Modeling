import numpy as np
import matplotlib.pyplot as plt

# Data points for CL
alpha_cl = [21.073825503355703, 29.06040268456376, 37.852348993288594, 46.97986577181208, 
            56.77852348993289, 66.04026845637584, 75.57046979865771, 83.28859060402685, 90]
CL_values = [1.075780089153046, 0.9450222882615156, 0.8469539375928677, 0.7340267459138187, 
             0.5913818722139673, 0.44576523031203563, 0.27637444279346207, 0.1337295690936107, 0]

# Fit a polynomial to the CL data
CL_fit_poly = np.polyfit(alpha_cl, CL_values, 8)

# Define the CL function
def CL(alpha):
    return (-1.19711067e-13 * alpha**8 
            + 5.20082162e-11 * alpha**7 
            - 9.56411265e-09 * alpha**6 
            + 9.67495221e-07 * alpha**5 
            - 5.84818549e-05 * alpha**4 
            + 2.14165289e-03 * alpha**3 
            - 4.57673865e-02 * alpha**2 
            + 4.99191880e-01 * alpha 
            - 9.02802242e-01)

# Data points for CD
alpha_cd = [19.59731543624161, 30.536912751677853, 40.80536912751678, 51.006711409395976, 
            61.61073825503356, 72.01342281879195, 81.74496644295301, 89.93288590604027]
CD_values = [0.26448736998514116, 0.5408618127786032, 0.8618127786032689, 1.200594353640416, 
             1.5245170876671619, 1.7682020802377414, 1.9108469539375927, 1.949479940564636]

# Fit a polynomial to the CD data
CD_fit_poly = np.polyfit(alpha_cd, CD_values, 7)

# Define the CD function
def CD(alpha):
    return (-1.10924209e-12 * alpha**7 
            + 4.08803237e-10 * alpha**6 
            - 6.17031482e-08 * alpha**5 
            + 4.93635653e-06 * alpha**4 
            - 2.32257424e-04 * alpha**3 
            + 6.77591188e-03 * alpha**2 
            - 8.78088130e-02 * alpha 
            + 5.59383006e-01)

# Generate a smooth curve for plotting
alpha_smooth = np.linspace(0, 90, 500)
CL_smooth = CL(alpha_smooth)
CD_smooth = CD(alpha_smooth)

# Plot the results
plt.figure(figsize=(12, 6))

# Plot CL
plt.subplot(1, 2, 1)
plt.plot(alpha_cl, CL_values, 'ro', label='Provided $C_L$ points')
plt.plot(alpha_smooth, CL_smooth, 'r-', label='$C_L$ (Fitted Curve)')
plt.xlabel(r'$\alpha$ (°)')
plt.ylabel('Lift Coefficient $C_L$')
plt.legend()
plt.grid(True)

# Plot CD
plt.subplot(1, 2, 2)
plt.plot(alpha_cd, CD_values, 'bo', label='Provided $C_D$ points')
plt.plot(alpha_smooth, CD_smooth, 'b-', label='$C_D$ (Fitted Curve)')
plt.xlabel(r'$\alpha$ (°)')
plt.ylabel('Drag Coefficient $C_D$')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
