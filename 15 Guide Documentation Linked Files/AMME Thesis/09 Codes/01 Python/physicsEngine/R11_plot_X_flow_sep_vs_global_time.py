# An improvement from R09
# Plot X_flow_sep vs. global_time from formatted_combined_history

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Specify the correct file path
file_path = 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\formatted_combined_history.txt'

# Load the data
df = pd.read_csv(file_path, sep=r'\s+', skiprows=0)

# Print the column names to identify the exact names
print("Column names:", df.columns)

# Extract the necessary columns (adjust the column names based on the actual output)
global_time = df.iloc[:, 0]  # assuming the first column is global_time
x_flow_sep = df.iloc[:, 4]  # assuming the fifth column is X_flow_sep

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(global_time, x_flow_sep, label='X_flow_sep vs. global_time')
plt.xlabel('Global Time [s]')
plt.ylabel('X Flow Sep [-]')
if np.all((x_flow_sep >= 0) & (x_flow_sep <= 1)):
    plt.ylim([-0.1, 1.1])
plt.title('X Flow Sep vs. Global Time')
plt.legend()
plt.grid(True)
plt.savefig('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\X_flow_sep.png')  # Save as PNG
plt.show()
