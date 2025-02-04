import numpy as np
import matplotlib.pyplot as plt

# NASA Mars Atmosphere Model function
def M12_nasa_mars_atmosphere_model(alt):
    if alt < 7000:  # Lower Martian Atmosphere: h < 7000 m
        temp = -31 - 0.000998 * alt
        pres = 699 * np.exp(-0.00009 * alt)
        rho = pres / (192.1 * (temp + 273.15))
    else:  # Upper Martian Atmosphere: h > 7000 m
        temp = -23.4 - 0.00222 * alt
        pres = 699 * np.exp(-0.00009 * alt)
        rho = pres / (192.1 * (temp + 273.15))

    a = (1.2941 * 188.92 * (temp + 273.15)) ** 0.5  # Mars gamma and Rgas

    return rho, a, pres, temp

# NASA Earth Atmosphere Model function
def M12_nasa_earth_atmosphere_model(alt):
    lap_rate = 0.00649  # (C)deg/meter
    t0 = 15.04          # (C)deg
    rho0 = 1.225        # Kg/m^3
    p0 = 101290         # Pa

    if alt < 11000:  # Troposphere: h < 11000 m
        temp = (t0 - (alt * lap_rate))
        pres = p0 * (((temp + 273.15) / (t0 + 273.15)) ** (9.81 / (lap_rate * 287.1)))
        rho = rho0 * (((temp + 273.15) / (t0 + 273.15)) ** ((9.81 - (lap_rate * 287.1)) / (lap_rate * 287.1)))

    elif alt >= 11000 and alt < 25000:  # Lower Stratosphere: 11000 m < h < 25000 m
        temp = -56.46
        pres = 22650 * np.exp(1.73 - 0.000157 * alt)
        rho = pres / (286.9 * (temp + 273.15))

    else:  # Upper Stratosphere: h > 25000 m
        temp = -131.21 + 0.00299 * alt
        pres = 2488 * (((temp + 273.15) / (216.6)) ** (-11.388))
        rho = pres / (286.9 * (temp + 273.15))

    a = (1.4 * 287.1 * (temp + 273)) ** 0.5

    return rho, a, pres, temp

# Original Earth Atmosphere model included with the custom physics engine
def M12_atmosphere_model(alt):
    lap_rate = 0.0065  # (C)deg/meter
    t0 = 15            # (C)deg
    rho0 = 1.225       # Kg/m^3
    p0 = 101310        # Pa

    if alt >= 11000:
        temp = -56.4
        pres = p0 * (0.2189 * (np.exp(9.81 * (11000 - alt) / (287.1 * (temp + 273)))))
        rho = rho0 * (0.2972 * (np.exp(9.81 * (11000 - alt) / (287.1 * (temp + 273)))))
    else:
        temp = (t0 - (alt * lap_rate))
        pres = p0 * (((temp + 273) / (t0 + 273)) ** (9.81 / (lap_rate * 287.1)))
        rho = rho0 * (((temp + 273) / (t0 + 273)) ** ((9.81 - (lap_rate * 287.1)) / (lap_rate * 287.1)))
    
    a = (1.4 * 287.1 * (temp + 273)) ** 0.5

    return rho, a, pres, temp

#===========================================================================================================================

# Generate altitude data
altitudes = np.linspace(0, 40000, 2000)  # Altitude range from 0 to 40000 meters

# Calculate the atmospheric properties at each altitude for Mars, NASA Earth, and Original Earth models
rho_values_mars, a_values_mars, pres_values_mars, temp_values_mars = [], [], [], []
rho_values_earth, a_values_earth, pres_values_earth, temp_values_earth = [], [], [], []
rho_values, a_values, pres_values, temp_values = [], [], [], []

#===========================================================================================================================

for alt in altitudes:
    # Mars atmosphere
    rho_mars, a_mars, pres_mars, temp_mars = M12_nasa_mars_atmosphere_model(alt)
    rho_values_mars.append(rho_mars)
    a_values_mars.append(a_mars)
    pres_values_mars.append(pres_mars)
    temp_values_mars.append(temp_mars)

    # NASA Earth atmosphere
    rho_earth, a_earth, pres_earth, temp_earth = M12_nasa_earth_atmosphere_model(alt)
    rho_values_earth.append(rho_earth)
    a_values_earth.append(a_earth)
    pres_values_earth.append(pres_earth)
    temp_values_earth.append(temp_earth)

    # Original Earth atmosphere
    rho, a, pres, temp = M12_atmosphere_model(alt)
    rho_values.append(rho)
    a_values.append(a)
    pres_values.append(pres)
    temp_values.append(temp)

############################################################################################################################
############################################################################################################################

# Plot for Mars Atmosphere Model
fig, ax1 = plt.subplots(figsize=(8, 10))

ax1.plot(rho_values_mars, altitudes, 'b-', label='Density (rho) - Mars')
ax1.grid(True)
ax1.set_xlabel(r'Density [kg/m$^3$]', color='b')
ax1.set_ylabel(r'Altitude [m]')
ax1.set_ylim([0, 40000])
ax1.tick_params('x', colors='b')

ax2 = ax1.twiny()
ax2.plot(a_values_mars, altitudes, 'r-', label='Speed of Sound (a) - Mars')
ax2.set_xlabel(r'Speed of Sound [m/s]', color='r')
ax2.tick_params('x', colors='r')

ax3 = ax1.twiny()
ax3.spines['top'].set_position(('outward', 40))
ax3.plot(pres_values_mars, altitudes, 'g-', label='Pressure (pres) - Mars')
ax3.set_xlabel(r'Pressure [Pa]', color='g')
ax3.tick_params('x', colors='g')

ax4 = ax1.twiny()
ax4.spines['top'].set_position(('outward', 80))
ax4.plot(temp_values_mars, altitudes, 'orange', label='Temperature (temp) - Mars')
ax4.set_xlabel(r'Temperature [°C]', color='orange')
ax4.tick_params('x', colors='orange')

# Adding a legend to the plot
fig.legend(loc="upper left", bbox_to_anchor=(0.45, 0.85), bbox_transform=ax1.transAxes)

plt.subplots_adjust(left=0.15, right=0.9, top=0.8, bottom=0.1, hspace=0.4, wspace=0.4)
# fig.tight_layout()
plt.title('Martian Atmospheric Properties vs Altitude')
plt.show()

############################################################################################################################
############################################################################################################################

# Plot for NASA Earth Atmosphere Model
fig, ax1 = plt.subplots(figsize=(8, 10))

ax1.plot(rho_values_earth, altitudes, 'b-', label='Density (rho) - NASA Earth')
ax1.grid(True)
ax1.set_xlabel(r'Density [kg/m$^3$]', color='b')
ax1.set_ylabel(r'Altitude [m]')
ax1.set_ylim([0, 40000])
ax1.tick_params('x', colors='b')

ax2 = ax1.twiny()
ax2.plot(a_values_earth, altitudes, 'r-', label='Speed of Sound (a) - NASA Earth')
ax2.set_xlabel(r'Speed of Sound [m/s]', color='r')
ax2.tick_params('x', colors='r')

ax3 = ax1.twiny()
ax3.spines['top'].set_position(('outward', 40))
ax3.plot(pres_values_earth, altitudes, 'g-', label='Pressure (pres) - NASA Earth')
ax3.set_xlabel(r'Pressure [Pa]', color='g')
ax3.tick_params('x', colors='g')

ax4 = ax1.twiny()
ax4.spines['top'].set_position(('outward', 80))
ax4.plot(temp_values_earth, altitudes, 'orange', label='Temperature (temp) - NASA Earth')
ax4.set_xlabel(r'Temperature [°C]', color='orange')
ax4.tick_params('x', colors='orange')

# Adding a legend to the plot
fig.legend(loc="upper left", bbox_to_anchor=(0.45, 0.85), bbox_transform=ax1.transAxes)

plt.subplots_adjust(left=0.15, right=0.9, top=0.8, bottom=0.1, hspace=0.4, wspace=0.4)
# fig.tight_layout()
plt.title('NASA Earth Atmospheric Properties vs Altitude')
plt.show()

############################################################################################################################
############################################################################################################################

# Plot for Original Earth Atmosphere Model
fig, ax1 = plt.subplots(figsize=(8, 10))

ax1.plot(rho_values, altitudes, 'b-', label='Density (rho) - Original Earth')
ax1.grid(True)
ax1.set_xlabel(r'Density [kg/m$^3$]', color='b')
ax1.set_ylabel(r'Altitude [m]')
ax1.set_ylim([0, 40000])
ax1.tick_params('x', colors='b')

ax2 = ax1.twiny()
ax2.plot(a_values, altitudes, 'r-', label='Speed of Sound (a) - Original Earth')
ax2.set_xlabel(r'Speed of Sound [m/s]', color='r')
ax2.tick_params('x', colors='r')

ax3 = ax1.twiny()
ax3.spines['top'].set_position(('outward', 40))
ax3.plot(pres_values, altitudes, 'g-', label='Pressure (pres) - Original Earth')
ax3.set_xlabel(r'Pressure [Pa]', color='g')
ax3.tick_params('x', colors='g')

ax4 = ax1.twiny()
ax4.spines['top'].set_position(('outward', 80))
ax4.plot(temp_values, altitudes, 'orange', label='Temperature (temp) - Original Earth')
ax4.set_xlabel(r'Temperature [°C]', color='orange')
ax4.tick_params('x', colors='orange')

# Adding a legend to the plot
fig.legend(loc="upper left", bbox_to_anchor=(0.45, 0.85), bbox_transform=ax1.transAxes)

plt.subplots_adjust(left=0.15, right=0.9, top=0.8, bottom=0.1, hspace=0.4, wspace=0.4)
# fig.tight_layout()
plt.title('Original Earth Atmospheric Properties vs Altitude')
plt.show()

############################################################################################################################
############################################################################################################################

# Plot for NASA Mars and NASA Earth Atmosphere Models Comparison
fig, ax1 = plt.subplots(figsize=(8, 10))

# Plotting rho (Density)
ax1.plot(rho_values_mars, altitudes, 'b--', label='Density (rho) - Mars')
ax1.plot(rho_values_earth, altitudes, 'b-', label='Density (rho) - Earth')
ax1.grid(True)
ax1.set_xlabel(r'Density [kg/m$^3$]', color='b')
ax1.set_ylabel(r'Altitude [m]')
ax1.set_ylim([0, 40000])
ax1.tick_params('x', colors='b')

ax2 = ax1.twiny()
# Plotting a (Speed of Sound)
ax2.plot(a_values_mars, altitudes, 'r--', label='Speed of Sound (a) - Mars')
ax2.plot(a_values_earth, altitudes, 'r-', label='Speed of Sound (a) - Earth')
ax2.set_xlabel(r'Speed of Sound [m/s]', color='r')
ax2.tick_params('x', colors='r')

ax3 = ax1.twiny()
ax3.spines['top'].set_position(('outward', 40))
# Plotting pres (Pressure)
ax3.plot(pres_values_mars, altitudes, 'g--', label='Pressure (pres) - Mars')
ax3.plot(pres_values_earth, altitudes, 'g-', label='Pressure (pres) - Earth')
ax3.set_xlabel(r'Pressure [Pa]', color='g')
ax3.tick_params('x', colors='g')

ax4 = ax1.twiny()
ax4.spines['top'].set_position(('outward', 80))
# Plotting temp (Temperature)
ax4.plot(temp_values_mars, altitudes, 'orange', linestyle='--', label='Temperature (temp) - Mars')
ax4.plot(temp_values_earth, altitudes, 'orange', linestyle='-', label='Temperature (temp) - Earth')
ax4.set_xlabel(r'Temperature [°C]', color='orange')
ax4.tick_params('x', colors='orange')

# Adding a legend to the plot
fig.legend(loc="upper left", bbox_to_anchor=(0.22, 1.00), bbox_transform=ax1.transAxes)

plt.subplots_adjust(left=0.15, right=0.9, top=0.8, bottom=0.1, hspace=0.4, wspace=0.4)
plt.title('NASA Mars vs NASA Earth Atmospheric Properties Comparison')
plt.show()

############################################################################################################################
############################################################################################################################


# import numpy as np
# import matplotlib.pyplot as plt

# # NASA Mars Atmosphere Model function
# def M12_nasa_mars_atmosphere_model(alt):
#     if alt < 7000:  # Lower Martian Atmosphere: h < 7000 m
#         temp = -31 - 0.000998 * alt
#         pres = 0.699 * np.exp(-0.00009 * alt)
#         rho = pres / (0.1921 * (temp + 273.15))
#     else:  # Upper Martian Atmosphere: h > 7000 m
#         temp = -23.4 - 0.00222 * alt
#         pres = 0.699 * np.exp(-0.00009 * alt)
#         rho = pres / (0.1921 * (temp + 273.15))

#     a = (1.2941 * 188.92 * (temp + 273.15)) ** 0.5  # Mars gamma and Rgas

#     return rho, a, pres, temp

# # NASA Earth Atmosphere Model function
# def M12_nasa_earth_atmosphere_model(alt):
#     lap_rate = 0.00649  # (C)deg/meter
#     t0 = 15.04          # (C)deg
#     rho0 = 1.225        # Kg/m^3
#     p0 = 101290         # Pa

#     if alt < 11000:  # Troposphere: h < 11000 m
#         temp = (t0 - (alt * lap_rate))
#         pres = p0 * (((temp + 273.15) / (t0 + 273.15)) ** (9.81 / (lap_rate * 287.1)))
#         rho = rho0 * (((temp + 273.15) / (t0 + 273.15)) ** ((9.81 - (lap_rate * 287.1)) / (lap_rate * 287.1)))

#     elif alt >= 11000 and alt < 25000:  # Lower Stratosphere: 11000 m < h < 25000 m
#         temp = -56.46
#         pres = 22.65 * np.exp(1.73 - 0.000157 * alt)
#         rho = pres / (0.2869 * (temp + 273.15))

#     else:  # Upper Stratosphere: h > 25000 m
#         temp = -131.21 + 0.00299 * alt
#         pres = 2.488 * (((temp + 273.15) / (216.6)) ** (-11.388))
#         rho = pres / (0.2869 * (temp + 273.15))

#     a = (1.4 * 287.1 * (temp + 273)) ** 0.5

#     return rho, a, pres, temp

# # Original Earth Atmosphere model included with the custom physics engine
# def M12_atmosphere_model(alt):
#     lap_rate = 0.0065  # (C)deg/meter
#     t0 = 15            # (C)deg
#     rho0 = 1.225       # Kg/m^3
#     p0 = 101310        # Pa

#     if alt >= 11000:
#         temp = -56.4
#         pres = p0 * (0.2189 * (np.exp(9.81 * (11000 - alt) / (287.1 * (temp + 273)))))
#         rho = rho0 * (0.2972 * (np.exp(9.81 * (11000 - alt) / (287.1 * (temp + 273)))))
#     else:
#         temp = (t0 - (alt * lap_rate))
#         pres = p0 * (((temp + 273) / (t0 + 273)) ** (9.81 / (lap_rate * 287.1)))
#         rho = rho0 * (((temp + 273) / (t0 + 273)) ** ((9.81 - (lap_rate * 287.1)) / (lap_rate * 287.1)))
    
#     a = (1.4 * 287.1 * (temp + 273)) ** 0.5

#     return rho, a, pres, temp

# # Generate altitude data for Earth
# altitudes = np.linspace(0, 40000, 2000)  # Altitude range from 0 to 12000 meters

# # Calculate the atmospheric properties at each altitude for Mars and Earth
# rho_values_mars, a_values_mars, pres_values_mars, temp_values_mars = [], [], [], []
# rho_values_earth, a_values_earth, pres_values_earth, temp_values_earth = [], [], [], []

# rho_values = []
# a_values = []
# pres_values = []
# temp_values = []


# for alt in altitudes:
#     # Mars atmosphere
#     rho_mars, a_mars, pres_mars, temp_mars = M12_nasa_mars_atmosphere_model(alt)
#     rho_values_mars.append(rho_mars)
#     a_values_mars.append(a_mars)
#     pres_values_mars.append(pres_mars)
#     temp_values_mars.append(temp_mars)

#     # Earth atmosphere
#     rho_earth, a_earth, pres_earth, temp_earth = M12_nasa_earth_atmosphere_model(alt)
#     rho_values_earth.append(rho_earth)
#     a_values_earth.append(a_earth)
#     pres_values_earth.append(pres_earth)
#     temp_values_earth.append(temp_earth)

#     rho, a, pres, temp = M12_atmosphere_model(alt)
#     rho_values.append(rho)
#     a_values.append(a)
#     pres_values.append(pres)
#     temp_values.append(temp)

# # Create the plot for Mars
# fig, ax1 = plt.subplots(figsize=(8, 10))

# # Plotting rho on the first x-axis for Mars
# ax1.plot(rho_values_mars, altitudes, 'b-', label='Density (rho) - Mars')
# ax1.set_xlabel(r'Density [kg/m$^3$]', color='b')
# ax1.set_ylabel(r'Altitude [m]')
# ax1.tick_params('x', colors='b')

# # Create the second x-axis and plot a (speed of sound) for Mars
# ax2 = ax1.twiny()
# ax2.plot(a_values_mars, altitudes, 'r-', label='Speed of Sound (a) - Mars')
# ax2.set_xlabel(r'Speed of Sound [m/s]', color='r')
# ax2.tick_params('x', colors='r')

# # Create the third x-axis and plot pres (pressure) for Mars
# ax3 = ax1.twiny()
# ax3.spines['top'].set_position(('outward', 40))  # Move the axis upwards
# ax3.plot(pres_values_mars, altitudes, 'g-', label='Pressure (pres) - Mars')
# ax3.set_xlabel(r'Pressure [kPa]', color='g')
# ax3.tick_params('x', colors='g')

# # Create the fourth x-axis and plot temp (temperature) for Mars
# ax4 = ax1.twiny()
# ax4.spines['top'].set_position(('outward', 80))  # Move the axis further upwards
# ax4.plot(temp_values_mars, altitudes, 'orange', label='Temperature (temp) - Mars')
# ax4.set_xlabel(r'Temperature [°C]', color='orange')
# ax4.tick_params('x', colors='orange')

# fig.tight_layout()

# # Show the plot for Mars
# # plt.title('Martian Atmospheric Properties vs Altitude')
# plt.show()

# # Create the plot for Earth
# fig, ax1 = plt.subplots(figsize=(8, 10))

# # Plotting rho on the first x-axis for Earth
# ax1.plot(rho_values_earth, altitudes, 'b-', label='Density (rho) - Earth')
# ax1.set_xlabel(r'Density [kg/m$^3$]', color='b')
# ax1.set_ylabel(r'Altitude [m]')
# ax1.tick_params('x', colors='b')

# # Create the second x-axis and plot a (speed of sound) for Earth
# ax2 = ax1.twiny()
# ax2.plot(a_values_earth, altitudes, 'r-', label='Speed of Sound (a) - Earth')
# ax2.set_xlabel(r'Speed of Sound [m/s]', color='r')
# ax2.tick_params('x', colors='r')

# # Create the third x-axis and plot pres (pressure) for Earth
# ax3 = ax1.twiny()
# ax3.spines['top'].set_position(('outward', 40))  # Move the axis upwards
# ax3.plot(pres_values_earth, altitudes, 'g-', label='Pressure (pres) - Earth')
# ax3.set_xlabel(r'Pressure [kPa]', color='g')
# ax3.tick_params('x', colors='g')

# # Create the fourth x-axis and plot temp (temperature) for Earth
# ax4 = ax1.twiny()
# ax4.spines['top'].set_position(('outward', 80))  # Move the axis further upwards
# ax4.plot(temp_values_earth, altitudes, 'orange', label='Temperature (temp) - Earth')
# ax4.set_xlabel(r'Temperature [°C]', color='orange')
# ax4.tick_params('x', colors='orange')

# fig.tight_layout()

# # Show the plot for Earth
# # plt.title('Earth Atmospheric Properties vs Altitude')
# plt.show()




# # Create the plot
# fig, ax1 = plt.subplots(figsize=(8, 10))

# # Plotting rho on the first x-axis
# ax1.plot(rho_values, altitudes, 'b-', label='Density (rho)')
# ax1.set_xlabel(r'Density [kg/m$^3$]', color='b')
# ax1.set_ylabel(r'Altitude [m]')
# ax1.tick_params('x', colors='b')

# # Create the second x-axis and plot a (speed of sound)
# ax2 = ax1.twiny()
# ax2.plot(a_values, altitudes, 'r-', label='Speed of Sound (a)')
# ax2.set_xlabel(r'Speed of Sound [m/s]', color='r')
# ax2.tick_params('x', colors='r')

# # Create the third x-axis and plot pres (pressure)
# ax3 = ax1.twiny()
# ax3.spines['top'].set_position(('outward', 40))  # Move the axis upwards
# ax3.plot(pres_values, altitudes, 'g-', label='Pressure (pres)')
# ax3.set_xlabel('Pressure (Pa)', color='g')
# ax3.tick_params('x', colors='g')

# # Create the fourth x-axis and plot temp (temperature)
# ax4 = ax1.twiny()
# ax4.spines['top'].set_position(('outward', 80))  # Move the axis further upwards
# ax4.plot(temp_values, altitudes, 'orange', label='Temperature (temp)')
# ax4.set_xlabel(r'Temperature [°C]', color='orange')
# ax4.tick_params('x', colors='orange')

# fig.tight_layout()

# # Show the plot
# # plt.title('Original Earth Atmospheric Properties vs Altitude')
# plt.show()



