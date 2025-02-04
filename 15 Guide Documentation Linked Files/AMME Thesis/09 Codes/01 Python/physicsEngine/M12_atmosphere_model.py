# Original Earth Atmosphere model included with the original custom physics engine

import numpy as np

"""
======================================================================================================================
    \brief   Calculates the ISA density and speed of sound at the given altitude
            NOTE Usage: rho, a = aero4560_atmos(alt)
======================================================================================================================
"""
def M12_atmosphere_model(alt):
    # Input:
    # alt   = altitude in meters
    # Output:
    # rho   = air density in kg/m^3
    # a     = speed of sound in m/s

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

    return rho, a