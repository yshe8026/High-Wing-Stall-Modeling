import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# flag for whether you are running xplane: yes (1), no (0)
flag_whether_running_xplane = 0

# Specify the correct file path
file_path = 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\formatted_combined_history.txt'

# Load the data
df = pd.read_csv(file_path, sep=r'\s+', skiprows=1)

# decide whether to trim the data according due to xplane scene initialization
if flag_whether_running_xplane == 1:
    # Extract the necessary columns
    global_time = df.iloc[150:, 0]  # assuming the first column is global_time
    alpha = df.iloc[150:, 10]  # assuming the 11th column is alpha
else:
    global_time = df.iloc[:, 0]  # assuming the first column is global_time
    alpha = df.iloc[:, 10]  # assuming the 11th column is alpha


# alpha ploting unit flag for choosing either deg (0) or rad (1)
plot_alpha_unit_flag = 0

if plot_alpha_unit_flag == 0:
    # unit conversion: rad to deg
    alpha_in_deg = alpha * (180/np.pi)
    # Calculate the average of alpha
    alpha_avg = alpha_in_deg.mean()
    print(f"Average of alpha: {alpha_avg} deg")
    # Plotting the data
    plt.figure(figsize=(10, 6))
    plt.plot(global_time, alpha_in_deg, label='alpha vs. global_time')
    plt.xlabel('Global Time [s]')
    plt.ylabel('Alpha [deg]')
    plt.title('Alpha vs. Global Time')
    plt.legend()
    plt.grid(True)
    plt.savefig('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\alpha.png')  # Save as PNG
    plt.show()


if plot_alpha_unit_flag == 1:
    # Calculate the average of alpha
    alpha_avg = alpha.mean()
    print(f"Average of alpha: {alpha_avg} rad")
    # Plotting the data
    plt.figure(figsize=(10, 6))
    plt.plot(global_time, alpha, label='alpha vs. global_time')
    plt.xlabel('Global Time [s]')
    plt.ylabel('Alpha [rad]')
    plt.title('Alpha vs. Global Time')
    plt.legend()
    plt.grid(True)
    plt.savefig('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\alpha.png')  # Save as PNG
    plt.show()
