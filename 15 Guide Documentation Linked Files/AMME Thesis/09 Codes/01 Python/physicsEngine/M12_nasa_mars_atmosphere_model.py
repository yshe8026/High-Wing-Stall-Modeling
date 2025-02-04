# NASA Earth Atmosphere Model
# https://www.grc.nasa.gov/www/k-12/airplane/atmosmrm.html

import numpy as np

"""
======================================================================================================================
    \brief   Calculates the ISA density and speed of sound at the given altitude
            NOTE Usage: rho, a = aero4560_atmos(alt)
======================================================================================================================
"""

def M12_nasa_mars_atmosphere_model(alt):
    # Input:
    # alt   = altitude in meters
    # Output:
    # rho   = air density in kg/m^3
    # a     = speed of sound in m/s

    if alt < 7000: # Lower Martian Atmosphere: h < 7000 m
        
        temp = -31 - 0.000998 * alt
        pres = 699 * np.exp(-0.00009 * alt)
        rho = pres / (192.1 * (temp + 273.15))

    else: # Upper Martian Atmosphere: h > 7000 m

        temp = -23.4 - 0.00222 * alt
        pres = 699 * np.exp(-0.00009 * alt)
        rho = pres / (192.1 * (temp + 273.15))
    
    
    a = (1.2941 * 188.92 * (temp + 273.15)) ** 0.5 # Mars gamma and Rgas, please refer to M12_gas_constant_of_planets

    return rho, a