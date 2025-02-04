import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# flag for whether you are running xplane: yes (1), no (0)
flag_whether_running_xplane = 0

# Define custom colors for plotting

# Customized internship poster blue
color00 = [35/255, 113/255, 148/255]
# Customized Red
color1 = [255/255, 36/255, 10/255]
# Customized Blue
color2 = [61/255, 142/255, 190/255]
# Customized Green
color3 = [21/255, 110/255, 72/255]
# Customized Pink
color4 = [217/255, 11/255, 125/255]
# Customized Orange
color5 = [228/255, 100/255, 0/255]
# Customized Purple
color6 = [120/255, 81/255, 169/255]
# Contrast Red
color7 = [239/255, 68/255, 68/255]
# Contrast Green
color8 = [0/255, 159/255, 117/255]
# Contrast Cyan
color9 = [136/255, 198/255, 237/255]
# Contrast Orange
color10 = [250/255, 163/255, 27/255]
# Contrast Purple
color11 = [57/255, 75/255, 160/255]
# Contrast Green
color12 = [130/255, 195/255, 65/255]
# Contrast Pink
color13 = [213/255, 71/255, 153/255]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# flag for whether loading more colors (1:yes) (0:no)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
flag_whether_load_more_colors = 0
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

if flag_whether_load_more_colors == 1:
    # Customized Red Shades
    color1_100 = [255/255, 36/255, 10/255]
    color1_80 = [255/255, 36/255, 10/255 * 0.8]
    color1_60 = [255/255, 36/255, 10/255 * 0.6]
    color1_40 = [255/255, 36/255, 10/255 * 0.4]
    color1_20 = [255/255, 36/255, 10/255 * 0.2]

    # Customized Pink Shades
    color4_100 = [217/255, 11/255, 125/255]
    color4_80 = [217/255, 11/255, 125/255 * 0.8]
    color4_60 = [217/255, 11/255, 125/255 * 0.6]
    color4_40 = [217/255, 11/255, 125/255 * 0.4]
    color4_20 = [217/255, 11/255, 125/255 * 0.2]

    # Customized Purple Shades
    color6_7 = [120/255, 81/255, 169/255 * 1.0]
    color6_6 = [120/255, 81/255, 169/255 * 0.9]
    color6_5 = [120/255, 81/255, 169/255 * 0.8]
    color6_4 = [120/255, 81/255, 169/255 * 0.7]
    color6_3 = [120/255, 81/255, 169/255 * 0.6]
    color6_2 = [120/255, 81/255, 169/255 * 0.5]
    color6_1 = [120/255, 81/255, 169/255 * 0.4]

# Specify the correct file path
# file_path = 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\formatted_combined_history.txt'
# ** Temporary solution ** need to change to 'formatted_stage2_combined_history'
file_path = 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\stage2_combined_history.txt'

# Load the data
df = pd.read_csv(file_path, sep=r'\s+', skiprows=1)


# decide whether to trim the data according due to xplane scene initialization
if flag_whether_running_xplane == 1:
    # Extract the necessary columns
    global_time = df.iloc[150:, 0]  # assuming the first column is global_time
    alpha = df.iloc[150:, 4]  # assuming the 5th column is alpha
    x_flow_sep = df.iloc[150:, 3]  # assuming the 4th column is X_flow_sep
    alpha_crit_1 = df.iloc[150:, 5]  # assuming the 6th column is alpha_crit
    alpha_crit_2 = df.iloc[150:, 6]  # assuming the 7th column is alpha_crit
    alpha_crit_3 = df.iloc[150:, 7]  # assuming the 8th column is alpha_crit
    CL = df.iloc[150:, 10]  # assuming the 9th column is CL
else:
    global_time = df.iloc[:, 0]  # assuming the first column is global_time
    alpha = df.iloc[:, 4]  # assuming the 5th column is alpha
    x_flow_sep = df.iloc[:, 3]  # assuming the 4th column is X_flow_sep
    alpha_crit_1 = df.iloc[:, 5]  # assuming the 6th column is alpha_crit
    alpha_crit_2 = df.iloc[:, 6]  # assuming the 7th column is alpha_crit
    alpha_crit_3 = df.iloc[:, 7]  # assuming the 8th column is alpha_crit
    CL = df.iloc[:, 10]  # assuming the 9th column is CL

# Calculate the average of alpha
alpha_avg = alpha.mean()
print(f"Average of alpha: {alpha_avg}")

# Plotting the data
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot X_flow_sep on the left y-axis
ax1.plot(global_time, x_flow_sep, label='X_flow_sep vs. global_time', color=color2)
ax1.set_xlabel('Global Time [s]')
ax1.set_ylabel('X_flow_sep [-]', color=color2)
if np.all((x_flow_sep >= 0) & (x_flow_sep <= 1)):
    ax1.set_ylim([-0.1, 1.1])
ax1.tick_params(axis='y', labelcolor=color2)
ax1.grid(True)

# Create a second y-axis for alpha
# alpha ploting unit flag for choosing either deg (0) or rad (1)
plot_alpha_unit_flag = 0

if plot_alpha_unit_flag == 0:
    # unit conversion: deg to deg (no conversion because in data sheet alpha is already in deg)
    alpha_in_deg = alpha 
    alpha_crit_1_in_deg = alpha_crit_1
    alpha_crit_2_in_deg = alpha_crit_2
    alpha_crit_3_in_deg = alpha_crit_3
    ax2 = ax1.twinx()
    ax2.plot(global_time, alpha_in_deg, label='alpha vs. global_time', color=color1, linestyle='--')
    ax2.plot(global_time, alpha_crit_1_in_deg, label='alpha_crit (interp 1) vs. global_time', color=color3, linestyle='-.')  
    ax2.plot(global_time, alpha_crit_2_in_deg, label='alpha_crit (interp 2) vs. global_time', color=color4, linestyle='-.')  
    ax2.plot(global_time, alpha_crit_3_in_deg, label='alpha_crit (interp 3) vs. global_time', color=color5, linestyle='-.')  
    ax2.plot(global_time, 2 * CL, label='CL vs. global_time', color=color6, linestyle='-') # CL is multiplied by 2 and displayed on ax2 scale
    ax2.set_ylabel('Alpha [deg]', color=color1)
    ax2.tick_params(axis='y', labelcolor=color1)

if plot_alpha_unit_flag == 1:
    # unit conversion: deg to rad
    alpha_in_rad = alpha * (np.pi/180)
    alpha_crit_1_in_rad = alpha_crit_1 * (np.pi/180)
    alpha_crit_2_in_rad = alpha_crit_2 * (np.pi/180)
    alpha_crit_3_in_rad = alpha_crit_3 * (np.pi/180)
    ax2 = ax1.twinx()
    ax2.plot(global_time, alpha_in_rad, label='alpha vs. global_time', color=color1, linestyle='--')
    ax2.plot(global_time, alpha_crit_1_in_rad, label='alpha_crit (interp 1) vs. global_time', color=color3, linestyle='-.')
    ax2.plot(global_time, alpha_crit_2_in_rad, label='alpha_crit (interp 2) vs. global_time', color=color4, linestyle='-.')  
    ax2.plot(global_time, alpha_crit_3_in_rad, label='alpha_crit (interp 3) vs. global_time', color=color5, linestyle='-.') 
    ax2.plot(global_time, 2 * CL, label='CL vs. global_time', color=color6, linestyle='-') # CL is multiplied by 2 and displayed on ax2 scale
    ax2.set_ylabel('Alpha [rad]', color=color1)
    ax2.tick_params(axis='y', labelcolor=color1)

# Title and legends
plt.title('X_flow_sep and Alpha vs. Global Time')
fig.tight_layout()
fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))
plt.savefig('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\alpha_crit_X_flow_sep_and_alpha.png')  # Save as PNG
plt.show()

