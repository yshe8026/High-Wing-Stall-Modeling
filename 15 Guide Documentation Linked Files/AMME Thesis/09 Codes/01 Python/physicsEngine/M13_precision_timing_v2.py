import numpy as np

# Warning: convert X_e, Y_e, Z_e (default in SWD axes) to in NED axes before feeding them to some of the functions below.

# Inherent aircraft atomic clock drift (time series input)
def atomic_clock_drift(nanosecond_per_day_drift_for_atomic_clock):
    delta_f_to_f = (nanosecond_per_day_drift_for_atomic_clock * (10**(-9)))/(24*3600) 
    return delta_f_to_f

# Doppler shift (Special Relativity) for time series
def doppler_shift(Vt, theta, psi, X_e, Y_e, Z_e):
    c = 299792458  # Speed of light in m/s
    ground_station_to_aircraft_distance = np.sqrt((X_e**2)+(Y_e**2)+(Z_e**2))
    relative_location_direction_NED = np.array([-X_e, -Y_e, -Z_e]).T
    aircraft_velocity_direction_NED = np.array([np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), -np.sin(theta)]).T
    delta_f_to_f = (1/ground_station_to_aircraft_distance) * np.einsum('ij,ij->i', relative_location_direction_NED, aircraft_velocity_direction_NED) * (Vt / c)
    return delta_f_to_f

# Time dilation calculations (Special Relativity)
def special_relativity(Vt):
    c = 299792458  # Speed of light in m/s
    return - (Vt**2) / (2 * c**2)

# Gravitational time dilation (General Relativity)
def gravitational_time_dilation(alt):
    c = 299792458  # Speed of light in m/s
    return (9.81 * alt) / (c**2)

def aircraft_acceleration_time_dilation(theta, psi, u_dot, X_e, Y_e, Z_e):
    c = 299792458  # Speed of light in m/s
    # Reshape u_dot to (100, 1) to broadcast correctly with (100, 3)
    u_dot = u_dot[:, np.newaxis]
    # Create the NED acceleration direction vector (shape (100, 3))
    a_aircraft_acceleration_in_NED = u_dot * np.array([np.cos(psi)*np.cos(theta), 
                                                       np.sin(psi)*np.cos(theta), 
                                                       -np.sin(theta)]).T  # shape (100, 3)
    # Create the aircraft position vector in NED (shape (100, 3))
    R_aircraft_position_in_NED = np.array([X_e, Y_e, Z_e]).T  # shape (100, 3)
    # Dot product over the two vectors (for each row), result is shape (100,)
    delta_f_to_f = np.einsum('ij,ij->i', R_aircraft_position_in_NED, a_aircraft_acceleration_in_NED) / (c**2)
    return delta_f_to_f  # Returns a time series result (shape (100,))

# Frame-dragging (Lense-Thirring Effect)
def frame_dragging(alt, theta_of_earth_spherical_coordinate_system):
    c = 299792458  # Speed of light in m/s
    G = 6.67430 * 10**-11  # Gravitational constant in m^3/kg/s^2
    M_earth = 5.972 * 10**24  # Mass of Earth in kg
    R_earth = 6.37 * 10**6  # Radius of Earth in meters
    J_earth = 7.08 * 10**33  # Earth's angular momentum in kg m^2/s
    
    r = R_earth + alt
    r_s = (2 * G * M_earth)/(c**2)
    alpha = J_earth/(M_earth * c)
    rho_squared_at_r = (r**2) + ((alpha**2) * (np.cos(theta_of_earth_spherical_coordinate_system)**2))
    Omega_angular_speed_of_inertial_frame_at_r = (r_s * alpha * r * c)/(rho_squared_at_r * ((r**2) + (alpha**2)) + (r_s * (alpha**2) * r * (np.sin(theta_of_earth_spherical_coordinate_system)**2)))
    v_space_flow_at_r = Omega_angular_speed_of_inertial_frame_at_r * r
    rho_squared_at_R_earth = (R_earth**2) + ((alpha**2) * (np.cos(theta_of_earth_spherical_coordinate_system)**2))
    Omega_angular_speed_of_inertial_frame_at_R_earth = (r_s * alpha * R_earth * c)/(rho_squared_at_R_earth * ((R_earth**2) + (alpha**2)) + (r_s * (alpha**2) * R_earth * (np.sin(theta_of_earth_spherical_coordinate_system)**2)))
    v_space_flow_at_R_earth = Omega_angular_speed_of_inertial_frame_at_R_earth * R_earth
    delta_T_to_T = (1/4)*((v_space_flow_at_R_earth - v_space_flow_at_r)/c)
    return delta_T_to_T

# Effective light speed reduction in air delay
def effective_light_speed_reduction_in_air_delay(rho):
    c = 299792458  # Speed of light in m/s
    rho_at_ground_station = 1.293 # kg/m^3
    effective_c_at_ground_station = c * (1 - 0.00029 * (rho_at_ground_station / 1.293))
    effective_c_at_aircraft = c * (1 - 0.00029 * (rho / 1.293))
    average_effective_c_during_signal_travel = (effective_c_at_aircraft + effective_c_at_ground_station)/2
    delta_T_to_T = (c - average_effective_c_during_signal_travel)/c
    return delta_T_to_T

# Ionospheric delay (not required for ground station communication)
def ionospheric_delay(alt):
    TEC = 10  # Total electron content
    f_L1 = 1.0 * 10**9  # GPS L1 signal frequency in Hz
    ionospheric_delay_meters = (40.3 * TEC * (10**16)) / (f_L1**2)
    ionospheric_total_meters = (35786000 - np.abs(alt))
    return ionospheric_delay_meters / ionospheric_total_meters

# Tropospheric delay (based on Saastamoinen model)
def tropospheric_delay(X_e, Y_e, Z_e):
    c = 299792458  # Speed of light in m/s
    zenith_angle = np.arctan2(np.sqrt(X_e**2 + Y_e**2), np.abs(Z_e))
    zenith_delay_meters = 2.3 * (1 - np.exp(-0.000116 * np.abs(Z_e)))
    delay_meters = zenith_delay_meters / np.cos(zenith_angle)
    total_meters = np.sqrt(X_e**2 + Y_e**2 + Z_e**2)
    return delay_meters / total_meters

# Sagnac effect
def sagnac_effect(X_e, Y_e, Z_e, theta_of_earth_spherical_coordinate_system):
    c = 299792458  # Speed of light in m/s
    R_earth = 6.37 * 10**6  # Radius of Earth in meters
    Omega_earth = 7.2921159 * 10**-5  # Earth's angular velocity in rad/s
    # Create R_vector as a time-series array (shape (100, 3))
    R_vector = np.array([-X_e, -Y_e, -Z_e]).T  # (100, 3)
    # Earth’s rotation velocity vector (constant for all time steps) (shape (3,))
    v_earth_rotation = Omega_earth * R_earth * np.sin(theta_of_earth_spherical_coordinate_system) * np.array([0, 1, 0])
    # Broadcast v_earth_rotation to have shape (100, 3) so it matches R_vector
    v_earth_rotation = np.tile(v_earth_rotation, (R_vector.shape[0], 1))  # shape (100, 3)
    # Compute the dot product element-wise over time
    delta_T_to_T = np.einsum('ij,ij->i', R_vector, v_earth_rotation) / c
    return delta_T_to_T

# Function to calculate and print total effects for time series
def calculate_and_print_effects_v2(alt, Vt, theta, psi, X_e, Y_e, Z_e, rho, u_dot, theta_of_earth_spherical_coordinate_system, nanosecond_per_day_drift_for_atomic_clock):
    # theta_of_earth_spherical_coordinate_system = np.deg2rad(90)  # (deg) At equator (0 deg is North pole, 180 deg is South pole)
    # nanosecond_per_day_drift_for_atomic_clock = np.full(100, 10)  # Time series for atomic clock drift over time

    # # Example input time series (you can replace these with real-time series data)
    # time_series_length = 100  # Define the length of your time series
    # time = np.linspace(0, 100, time_series_length)  # Example time vector (100 seconds)

    # # Define each parameter as a function of time
    # theta = np.deg2rad(np.full(time_series_length, 0))  # Level flight over time
    # psi = np.deg2rad(np.full(time_series_length, 90))  # Fly towards east over time
    # X_e = np.zeros(time_series_length)  # Example time series for X_e
    # Y_e = -1000 + 10 * time  # Example time series for Y_e (changing position over time)
    # Z_e = np.full(time_series_length, 1000)  # Example time series for Z_e
    # rho = np.full(time_series_length, 1.225)  # Air density as a constant for this example

    # # Accelerations as functions of time
    # u_dot = 9.81 * np.full(time_series_length, 0.1)  # Forward acceleration (constant in this example)
    # v_dot = 9.81 * np.zeros(time_series_length)  # No lateral acceleration in this example
    # w_dot = 9.81 * 0.1 * np.deg2rad(3) * np.ones(time_series_length)  # Example AoA changes

    delta_f_to_f_clock = atomic_clock_drift(nanosecond_per_day_drift_for_atomic_clock)
    delta_f_to_f_doppler = doppler_shift(Vt, theta, psi, X_e, Y_e, Z_e)
    delta_f_to_f_special_relativity = special_relativity(Vt)
    delta_f_to_f_gravity = gravitational_time_dilation(alt)
    delta_f_to_f_acceleration = aircraft_acceleration_time_dilation(theta, psi, u_dot, X_e, Y_e, Z_e)
    delta_T_to_T_frame_dragging = frame_dragging(alt, theta_of_earth_spherical_coordinate_system)
    delta_T_to_T_light_speed_reduction = effective_light_speed_reduction_in_air_delay(rho)
    delta_T_to_T_ionospheric = ionospheric_delay(alt)
    delta_T_to_T_tropospheric = tropospheric_delay(X_e, Y_e, Z_e)
    delta_T_to_T_sagnac = sagnac_effect(X_e, Y_e, Z_e, theta_of_earth_spherical_coordinate_system)

    print("Total calculations over time series are complete.")
    return (delta_f_to_f_clock, delta_f_to_f_doppler, delta_f_to_f_special_relativity,
            delta_f_to_f_gravity, delta_f_to_f_acceleration, delta_T_to_T_frame_dragging,
            delta_T_to_T_light_speed_reduction, delta_T_to_T_ionospheric, delta_T_to_T_tropospheric, delta_T_to_T_sagnac)

# Example usage
theta_of_earth_spherical_coordinate_system = np.deg2rad(90)  # (deg) At equator (0 deg is North pole, 180 deg is South pole)
nanosecond_per_day_drift_for_atomic_clock = np.full(100, 10)  # Time series for atomic clock drift over time

alt = np.linspace(1000, 5000, 100)
Vt = np.linspace(69.44, 80, 100)

# Example input time series (you can replace these with real-time series data)
time_series_length = 100  # Define the length of your time series
time = np.linspace(0, 100, time_series_length)  # Example time vector (100 seconds)

# Define each parameter as a function of time
theta = np.deg2rad(np.full(time_series_length, 0))  # Level flight over time
psi = np.deg2rad(np.full(time_series_length, 90))  # Fly towards east over time
X_e = np.zeros(time_series_length)  # Example time series for X_e
Y_e = -1000 + 10 * time  # Example time series for Y_e (changing position over time)
Z_e = np.full(time_series_length, 1000)  # Example time series for Z_e
rho = np.full(time_series_length, 1.225)  # Air density as a constant for this example

# Accelerations as functions of time
u_dot = 9.81 * np.full(time_series_length, 0.1)  # Forward acceleration (constant in this example)
v_dot = 9.81 * np.zeros(time_series_length)  # No lateral acceleration in this example
w_dot = 9.81 * 0.1 * np.deg2rad(3) * np.ones(time_series_length)  # Example AoA changes

calculate_and_print_effects_v2(alt, Vt, theta, psi, X_e, Y_e, Z_e, rho, u_dot, theta_of_earth_spherical_coordinate_system, nanosecond_per_day_drift_for_atomic_clock)

