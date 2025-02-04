import numpy as np
import matplotlib.pyplot as plt

# This is used as the development concept of the frozen field for Richard Directional Wind Module

k_karman_constant = 0.4 # Karman constant
H_0_surface_roughness_height = 3 # m
V_a0_friction_velocity = 0.4 # m/s

# Generate aircraft altitude from 0 m to 10000 m
aircraft_altitude = np.linspace(0.1, 10000, 1000)  # Avoid log(0) by starting from 0.1 m

# Calculate V_a
# if aircraft_altitude > 3:
V_a = (V_a0_friction_velocity / k_karman_constant) * np.log(aircraft_altitude / H_0_surface_roughness_height)
# else: 
#     V_a = 0

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(aircraft_altitude, V_a, label=r'$V_a$ vs Aircraft Altitude')
plt.xlabel('Aircraft Altitude [m]')
plt.ylabel(r'$V_a$ [m/s]')
plt.title('Average Wind Velocity vs Aircraft Altitude')
plt.legend()
plt.grid(True)
plt.show()

# "In a flight simulator, you can usually do not consider the average wind changes with time, that so-called "frozen field" hypothesis.[3]" - Modeling and Simulation of Atmospheric Synthetic Wind Field in Flight Simulation (Weiting Cui and Xiaoyong Lei, BUAA (Beijing University of Aeronautics and Astronautics)).
# Ref [3] w. L. Wang. Atmospheric Wind Field Modeling and its Application [D].National University of Defense Technology, 2009.
# In the simulation, we take the Average Wind Velocity as our "frozen field."