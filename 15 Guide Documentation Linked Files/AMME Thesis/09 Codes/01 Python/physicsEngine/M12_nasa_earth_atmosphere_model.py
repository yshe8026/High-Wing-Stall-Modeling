# NASA Earth Atmosphere Model
# https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html

import numpy as np

"""
======================================================================================================================
    \brief   Calculates the ISA density and speed of sound at the given altitude
            NOTE Usage: rho, a = aero4560_atmos(alt)
======================================================================================================================
"""

def M12_nasa_earth_atmosphere_model(alt):
    # Input:
    # alt   = altitude in meters
    # Output:
    # rho   = air density in kg/m^3
    # a     = speed of sound in m/s

    lap_rate = 0.00649  # (C)deg/meter
    t0 = 15.04            # (C)deg
    rho0 = 1.225       # Kg/m^3
    p0 = 101290        # Pa

    if alt < 11000: # Troposphere: h < 11000 m
        temp = (t0 - (alt * lap_rate))
        pres = p0 * (((temp + 273.15) / (t0 + 273.15)) ** (9.81 / (lap_rate * 287.1)))
        rho = rho0 * (((temp + 273.15) / (t0 + 273.15)) ** ((9.81 - (lap_rate * 287.1)) / (lap_rate * 287.1)))
        
        # # Much Simpler Formulation, but less fexible on potential adjustments on lapse rate:
        # temp = 15.04 - 0.00649 * alt
        # pres = 101.29 * (((temp + 273.15)/(288.08)) ** 5.256)
        # rho = pres / (0.2869 * (temp + 273.15))

    elif alt >= 11000 and alt < 25000: # Lower Stratosphere: 11000 m < h < 25000 m
        temp = -56.46
        pres = 22650 * np.exp(1.73 - 0.000157 * alt)
        rho = pres / (286.9 * (temp + 273.15))

    else: # Upper Stratosphere: h > 25000 m
        temp = -131.21 + 0.00299 * alt
        p = 2488 * (((temp + 273.15)/(216.6)) ** (-11.388))
        rho = pres / (286.9 * (temp + 273.15))
    
    a = (1.4 * 287.1 * (temp + 273)) ** 0.5

    return rho, a
