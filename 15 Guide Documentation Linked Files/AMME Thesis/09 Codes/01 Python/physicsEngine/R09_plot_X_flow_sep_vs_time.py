import numpy as np
import matplotlib.pyplot as plt

# Load the data from the file
data = np.loadtxt('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\X_flow_sep_history.txt')

# Define the time interval and create a time array
time_interval = 0.01
time = np.arange(0, len(data) * time_interval, time_interval)

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time, data)
plt.xlabel('Time (s)')
plt.ylabel('X')
plt.title('Time Series of X')
plt.ylim(-1, 3)  # Restrict the vertical axis
plt.grid(True)
plt.show()
