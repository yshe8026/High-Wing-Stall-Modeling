# Precision Time Keeping Module that include functions that account for deviation of measured aircraft time from ground station time (ground station is at (X_e = 0, Y_e = 0, Z_e = 0))
'''Since this module currently does not have a direct effect in the behavior of the flight simulation, it will be called at data analysis stage R32 for the purpose of additional insight.'''

import numpy as np

# # Constants
# c = 299792458  # Speed of light in m/s
# G = 6.67430 * 10**-11  # Gravitational constant in m^3/kg/s^2
# M_earth = 5.972 * 10**24  # Mass of Earth in kg
# R_earth = 6.37 * 10**6  # Radius of Earth in meters
# Omega_earth = 7.2921159 * 10**-5  # Earth's angular velocity in rad/s
# J_earth = 7.08 * 10**33  # Earth's angular momentum in kg m^2/s

# TEC = 10 #* 10**16  /m^2 # Total electron content (in TECU)
# f_L1 = 1.0 * 10**9  # GPS L1 signal frequency in Hz

theta_of_earth_spherical_coordinate_system = np.deg2rad(90) # (deg) At equator (0 deg is North pole, 180 deg is South pole)
nanosecond_per_day_drift_for_atomic_clock = 10 # Positive: clock is faster than real time; Negative: clock is slower than real time

theta = np.deg2rad(0) # level flight
psi = np.deg2rad(90) # fly towards east
X_e = 0
Y_e = -1000 # distance to the west of the ground station (blueshift) (approach)
Z_e = -1000
rho = 1.225 # air density of location of the aircraft

u_dot = 9.81 * 0.1 # 0.1 g forward acceleration
v_dot = 9.81 * 0.0
w_dot = 9.81 * 0.1 * np.deg2rad(3) # roughly 3 deg of AoA at level flight

# Inherent aircraft atomic clock drift
def atomic_clock_drift(nanosecond_per_day_drift_for_atomic_clock):
    delta_f_to_f = (nanosecond_per_day_drift_for_atomic_clock * (10**(-9)))/(24*3600) # Assume a common state-of-the-art 10 ns drift per day (clock faster than actual time: positive drift) 
    return delta_f_to_f

# Function to calculate Doppler shift (Special Relativity)
def doppler_shift(Vt, theta, psi, X_e,Y_e,Z_e):
    c = 299792458  # Speed of light in m/s
    # calculate the maximum doppler frequency drift (The real values are very different from situation to situation
    # Distance
    ground_station_to_aircraft_distance = np.sqrt((X_e**2)+(Y_e**2)+(Z_e**2))
    # Approximation
    # Ground start point relative to aircraft
    relative_location_direction_NED = np.array([-X_e, -Y_e, -Z_e])
    # Aircraft velocity (approximate) relative to ground start point
    aircraft_velocity_direction_NED = np.array([np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), -np.sin(theta)])
    c = 299792458  # Speed of light in m/s
    delta_f_to_f = (1/ground_station_to_aircraft_distance) * np.dot(relative_location_direction_NED, aircraft_velocity_direction_NED) * (Vt / c)
    # delta_f_to_f = (Vt / c)
    return delta_f_to_f

# Function to calculate special relativistic time dilation
def special_relativity(Vt):
    c = 299792458  # Speed of light in m/s
    delta_f_to_f = - (Vt**2) / (2 * c**2)
    return delta_f_to_f

# Function to calculate gravitational time dilation (General Relativity)
def gravitational_time_dilation(alt):
    c = 299792458  # Speed of light in m/s
    delta_h = alt  # altitude in meters
    delta_f_to_f = (9.81 * delta_h) / (c**2)
    return delta_f_to_f

def aircraft_acceleration_time_dilation(theta, psi, u_dot):
    c = 299792458  # Speed of light in m/s
    # Assume only u_dot has major impact on time dilation over the course of a flight trip (Here we approximate the acceleration/deceleration direction as the u direction)
    a_aircraft_acceleration_in_NED = u_dot * np.array([np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), -np.sin(theta)])
    R_aircraft_position_in_NED = np.array([X_e, Y_e, Z_e])
    # print(a_aircraft_acceleration_in_NED)
    # print(R_aircraft_position_in_NED)
    delta_f_to_f = np.dot(R_aircraft_position_in_NED, a_aircraft_acceleration_in_NED) / (c**2)
    return delta_f_to_f # delta_f_to_f


# Function to calculate frame-dragging (Lense-Thirring Effect)
def frame_dragging(alt, theta_of_earth_spherical_coordinate_system):
    # r = R_earth + alt
    # Omega_LT = (2 * G * J_earth) / (c**2 * r**3)

    c = 299792458  # Speed of light in m/s
    G = 6.67430 * 10**-11  # Gravitational constant in m^3/kg/s^2
    M_earth = 5.972 * 10**24  # Mass of Earth in kg
    R_earth = 6.37 * 10**6  # Radius of Earth in meters
    J_earth = 7.08 * 10**33  # Earth's angular momentum in kg m^2/s
    
    # Aircraft distance to the center of Earth
    r = R_earth + alt
    # Schwarzschild radius of the Earth
    r_s = (2 * G * M_earth)/(c**2)
    alpha = J_earth/(M_earth * c)
    rho_squared_at_r = (r**2) + ((alpha**2) * (np.cos(theta_of_earth_spherical_coordinate_system)**2))
    Omega_angular_speed_of_inertial_frame_at_r = (r_s * alpha * r * c)/(rho_squared_at_r * ((r**2) + (alpha**2)) + (r_s * (alpha**2) * r * (np.sin(theta_of_earth_spherical_coordinate_system)**2))) # rad/s
    v_space_flow_at_r = Omega_angular_speed_of_inertial_frame_at_r * r
    rho_squared_at_R_earth = (R_earth**2) + ((alpha**2) * (np.cos(theta_of_earth_spherical_coordinate_system)**2))
    Omega_angular_speed_of_inertial_frame_at_R_earth = (r_s * alpha * R_earth * c)/(rho_squared_at_R_earth * ((R_earth**2) + (alpha**2)) + (r_s * (alpha**2) * R_earth * (np.sin(theta_of_earth_spherical_coordinate_system)**2))) # rad/s
    v_space_flow_at_R_earth = Omega_angular_speed_of_inertial_frame_at_R_earth * R_earth
    # A rough estimation
    delta_T_to_T = (1/4)*((v_space_flow_at_R_earth - v_space_flow_at_r)/c)

    return delta_T_to_T #((2 *np.pi * R_earth) / (Omega_angular_speed_of_inertial_frame_at_r * r))/(365*24*3600)  #(Omega_angular_speed_of_inertial_frame * r) / c

def effective_light_speed_reduction_in_air_delay(rho):
    c = 299792458  # Speed of light in m/s
    rho_at_ground_station = 1.293 # kg/m^3
    effective_c_at_ground_station = c * (1 - 0.00029 * (rho_at_ground_station / 1.293))
    effective_c_at_aircraft = c * (1 - 0.00029 * (rho / 1.293))
    average_effective_c_during_signal_travel = (effective_c_at_aircraft + effective_c_at_ground_station)/2
    delta_T_to_T = (c - average_effective_c_during_signal_travel)/c
    return delta_T_to_T
    
# Function to calculate ionospheric delay
def ionospheric_delay(alt):
    TEC = 10 #* 10**16  /m^2 # Total electron content (in TECU)
    f_L1 = 1.0 * 10**9  # GPS L1 signal frequency in Hz
    ionospheric_delay_meters = (40.3 * TEC * (10**16)) / (f_L1**2)
    ionospheric_total_maters = (35786000 - np.abs(alt))
    delta_T_to_T = ionospheric_delay_meters / ionospheric_total_maters
    return delta_T_to_T # Not needed for communication with ground station

# Function to calculate tropospheric delay (based on Saastamoinen model)
def tropospheric_delay(X_e, Y_e, Z_e):
    # Assume the weather is dry (not wet or rainy)
    c = 299792458  # Speed of light in m/s
    # Zenith in radians
    zenith_angle_of_aircraft_relative_to_ground_station = np.arctan2( np.sqrt(X_e**2 + Y_e**2), np.abs(Z_e))
    # Zenith delay meters (delay of a signal in zenith direction (going vertically down))
    zenith_delay_meters = 2.3 * (1 - np.exp(-0.000116 * np.abs(Z_e)))
    # Z = np.deg2rad(zenith_angle_of_aircraft_relative_to_ground_station)
    delay_meters = zenith_delay_meters / np.cos(zenith_angle_of_aircraft_relative_to_ground_station)  # Approximate tropospheric delay in meters
    total_meters = np.sqrt((X_e**2)+(Y_e**2)+(Z_e**2))
    delta_T_to_T = delay_meters/total_meters
    return delta_T_to_T

# Function to calculate Sagnac effect
def sagnac_effect(X_e,Y_e,Z_e, theta_of_earth_spherical_coordinate_system):
    c = 299792458  # Speed of light in m/s
    R_earth = 6.37 * 10**6  # Radius of Earth in meters
    Omega_earth = 7.2921159 * 10**-5  # Earth's angular velocity in rad/s
    R_vector_from_aircraft_to_ground_station_in_NED = np.array([-X_e, -Y_e, -Z_e])
    v_vector_of_the_earth_rotation_velocity_at_ground_station_in_NED = Omega_earth * R_earth * np.sin(theta_of_earth_spherical_coordinate_system) * np.array([0, 1, 0]) # Earth rotation from west to east
    R_vector_from_aircraft_to_ground_station_in_NED_magnitude = np.linalg.norm(R_vector_from_aircraft_to_ground_station_in_NED)  # Use np.linalg.norm to calculate magnitude
    R_vector_from_aircraft_to_ground_station_in_NED_direction = R_vector_from_aircraft_to_ground_station_in_NED / R_vector_from_aircraft_to_ground_station_in_NED_magnitude
    delta_T_to_T = (np.dot(R_vector_from_aircraft_to_ground_station_in_NED_direction, v_vector_of_the_earth_rotation_velocity_at_ground_station_in_NED))/(c)

    # r = R_earth + alt
    return  delta_T_to_T # ((2 * Omega_earth * r * Vt) / c**2) / (R_earth/c)

# Main function to calculate total frequency shift and print individual effects
def calculate_and_print_effects(alt, Vt):

    # Inherent aircraft atomic clock drift
    delta_f_to_f_clock = atomic_clock_drift(nanosecond_per_day_drift_for_atomic_clock)

    # Doppler shift
    delta_f_to_f_doppler = doppler_shift(Vt, theta, psi, X_e, Y_e, Z_e)
    
    #----------------------------------------------------------------------------------------------------

    # Special relativity time dilation
    delta_f_to_f_special_relativity = special_relativity(Vt)
    
    # Gravitational time dilation
    delta_f_to_f_gravity = gravitational_time_dilation(alt)

    # Aircraft acceleration time dilation
    delta_f_to_f_acceleration = aircraft_acceleration_time_dilation(theta, psi, u_dot)

    #----------------------------------------------------------------------------------------------------
    
    # Frame-dragging effect
    delta_T_to_T_frame_dragging = frame_dragging(alt, theta_of_earth_spherical_coordinate_system)

    # Light speed reduction (effective)
    delta_T_to_T_light_speed_reduction = effective_light_speed_reduction_in_air_delay(rho)
    
    # Ionospheric delay
    delta_T_to_T_ionospheric = ionospheric_delay(alt)
    
    # Tropospheric delay
    delta_T_to_T_tropospheric = tropospheric_delay(X_e, Y_e, Z_e)
    
    # Sagnac effect
    delta_T_to_T_sagnac = sagnac_effect(X_e,Y_e,Z_e, theta_of_earth_spherical_coordinate_system)

    # Print individual effects
    print(f"------------------------------------------------------------------------------------------------------------------------------------------")
    print(f"Inherent aircraft atomic clock drift (artificial) (Δf / f): {delta_f_to_f_clock:.10e}")
    print(f"Doppler shift (optical illusion) (Δf / f): {delta_f_to_f_doppler:.10e}")
    print(f"------------------------------------------------------------------------------------------------------------------------------------------")
    print(f"Special relativity time dilation (real but not permanent) (Δf / f): {delta_f_to_f_special_relativity:.10e}")
    print(f"Gravitational time dilation (real and permanent) (Δf / f): {delta_f_to_f_gravity:.10e}")
    print(f"Aircraft acceleration time dilation (real and permanent) (Δf / f): {delta_f_to_f_acceleration:.10e}")
    print(f"------------------------------------------------------------------------------------------------------------------------------------------")
    print(f"Frame-dragging (Lense-Thirring effect) (signal delay) (ΔT / T): {delta_T_to_T_frame_dragging:.10e}")
    print(f"Effective light speed reduction delay contribution (signal delay) (ΔT / T): {delta_T_to_T_light_speed_reduction:.10e}")
    print(f"Ionospheric delay contribution (no need to include for ground communication) (signal delay) (ΔT / T): {delta_T_to_T_ionospheric:.10e}")
    print(f"Tropospheric delay contribution (signal delay) (ΔT / T): {delta_T_to_T_tropospheric:.10e}")
    print(f"Sagnac effect (signal delay) (ΔT / T): {delta_T_to_T_sagnac:.10e}")
    print(f"==========================================================================================================================================")
    
    # Summing up all apparent effects  
    total_apparent_delta_f_to_f = (delta_f_to_f_clock
                                 + delta_f_to_f_doppler)

    # Summing up all real effects
    total_real_delta_f_to_f = (delta_f_to_f_special_relativity
                             + delta_f_to_f_gravity
                             + delta_f_to_f_acceleration)
    
    # Summing up all delay effects
    total_delta_T_to_T = (delta_T_to_T_frame_dragging
                        + delta_T_to_T_light_speed_reduction
                        + delta_T_to_T_ionospheric
                        + delta_T_to_T_tropospheric
                        + delta_T_to_T_sagnac)
    
    print(f"Total frequency apparent shift (Δf / f): {total_apparent_delta_f_to_f:.10e}")
    print(f"Total frequency real shift (Δf / f): {total_real_delta_f_to_f:.10e}")
    print(f"Total time delay ratio (ΔT / T): {total_delta_T_to_T:.10e}")
    print(f"==========================================================================================================================================")
    return total_apparent_delta_f_to_f, total_real_delta_f_to_f, total_delta_T_to_T

# Example usage:
if __name__ == "__main__":
    altitude = 1000  # Example altitude in meters (1000 meters)
    aircraft_speed = 60  # Example speed in m/s (about 117 kn)
    
    total_apparent_delta_f_to_f, total_real_delta_f_to_f, total_delta_T_to_T = calculate_and_print_effects(altitude, aircraft_speed)

