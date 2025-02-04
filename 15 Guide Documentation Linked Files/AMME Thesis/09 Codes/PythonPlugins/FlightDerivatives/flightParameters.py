import numpy as np

"""
Description:
    - This python script is used to define lookup tables for custom aircraft.
      + Each dictionary contains all the flight data required to simulate a specific aircraft.

Use:
    - To create a new aircraft, specify the following:

        'Name': {
            'FlightData' : {
                'Inertial' : {
                    'g'   : x ,
                    'm'   : x ,
                    'Ixx' : x ,
                    'Iyy' : x ,
                    'Izz' : x ,
                    'Ixz' : x ,
                },
                'Geometric' : {
                    'S' : x ,
                    'c' : x ,
                    'b' : x 
                },
                'Propeller' : {
                    'P_max' : x ,
                    'eta'   : x 
                },
                'ControlLimits' : {
                    'Lower' : np.array([0, np.deg2rad(-x), np.deg2rad(-x), np.deg2rad(-x), np.deg2rad(-x)]),
                    'Upper' : np.array([1, np.deg2rad(x) , np.deg2rad(x) , np.deg2rad(x), np.deg2rad(x) ])
                },
                'Aero' : {
                    'alpha_o' : x ,
                    'Cdo'  : x ,
                    'k'    : x ,
                    'CLa'  : x ,
                    'CLq'  : x ,
                    'CLad' : x ,
                    'CLde' : x ,
                    'CLdf' : x ,
                    'CLo'  :  x ,  # defined as: CLo = -CLa * alpha_o 
                    'Cyb'  : x ,
                    'Cybd' : x ,
                    'Cyp'  : x ,
                    'Cyr'  : x ,
                    'Cyda' : x ,
                    'Cydr' : x ,
                    'Cmo'  : x ,
                    'Cma'  : x ,
                    'Cmq'  : x ,
                    'Cmad' : x ,
                    'Cmde' : x ,
                    'Cmdf' : x ,
                    'Cnb'  : x ,
                    'Cnbd' : x ,
                    'Cnp'  : x ,
                    'Cnr'  : x ,
                    'Cnda' : x ,
                    'Cndr' : x ,
                    'Clb'  : x ,
                    'Clbd' : x ,
                    'Clp'  : x ,
                    'Clr'  : x ,
                    'Clda' : x ,
                    'Cldr' : x 
                }
            }
            },

    - You must use the config file to select that aircraft in the plug-in
        + In config.ini:
          [Aircraft_Settings]
          aircraft = Name
"""

aircraft_info = {
    'J400': {
        'FlightData' : {
            'Inertial' : {
                'g'   : 9.81 * (9.81/9.81),
                'm'   : 700 ,
                'Ixx' : 734.8230 ,
                'Iyy' : 1022.700 ,
                'Izz' : 1067.100,
                'Ixz' : 31.7037 * 0
            },
            'Geometric' : {
                'S' : 8.019,
                'c' : 0.99,
                'b' : 8.10
            },
            'Propeller' : {
                'P_max' : 89484,
                'eta'   : 0.8
            },
            'ControlLimits' : {                 # -18.4 (but this need a calibration correction: 18.4 * 2 - 9.0 = 27.8)
                'Lower' : np.array([0, np.deg2rad(-27.8), np.deg2rad(-24), np.deg2rad(-23.21), np.deg2rad(-40)]), # dT, da, de, dr, df (Likely wrong label*) # Richard think it should actually be [dT, de, da, dr, df] instead *
                'Upper' : np.array([1, np.deg2rad(9.0) , np.deg2rad(24) , np.deg2rad(23.21), np.deg2rad(-0)])  # dT, da, de, dr, df (Likely wrong label*) # Richard think it should actually be [dT, de, da, dr, df] instead *
            },
            'Aero' : {
                'alpha_o' : -0.0698,              # zero lift (o) alpha
                'Cdo' : 0.0061 * 8 * (10/9.81),               # zero lift (o) drag coefficient (Cd)     # 0.049745
                'k' : 0.0483 * (0.09/0.0483)* (10/9.81),     # induced drag coefficient (k)             # 0.0917
                'CLa' : 4.9635,                   # Lift coefficient (CL) wrt alpha (a)
                'CLq' : 9.8197 * (8/9.8197),                   # Lift coefficient (CL) wrt pitch rate (q)
                'CLad' : -1.987,                  # Lift coefficient (CL) wrt alpha dot (ad)
                'CLde' : 0.0007976 * 1000,        # Lift coefficent (CL) wrt elevator deflection (de)
                'CLdf' : -1.75, # 0.0229 * 10,             # Lift coefficent (CL) wrt flap deflection (df)
                'CLo' : 0.36486, #* (0.8254/0.36486),                  # lift coefficent (CL) at zero alpha (o)
                'Cyb' : -0.23,                    # side force coefficent (Cy) wrt beta(b)
                'Cybd' : -0.0032,                 # side force coefficent (Cy) wrt beta dot(bd)
                'Cyp' : 0.1624,                   # side force coefficent (Cy) wrt roll rate(p)
                'Cyr' : 0.1889,                   # side force coefficent (Cy) wrt yaw rate(r)
                'Cyda' : 0.00033,                 # side force coefficent (Cy) wrt aileron deflection(da)
                'Cydr' : -0.0031 * (-10),         # side force coefficent (Cy) wrt rudder deflection(dr)
                'Cmo' : 0.07575,                  # pitching moment coefficent (Cm) at zero alpha (o)
                'Cma' : -2.6001,                  # pitching moment coefficent (Cm) wrt alpha (a)
                'Cmq' : -23.6813,                 # pitching moment coefficent (Cm) wrt pitch rate (p)
                'Cmad' : -2.210,                  # pitching moment coefficent (Cm) wrt alpha dot(ad)
                'Cmde' : -0.036 * 100,            # pitching moment coefficent (Cm) wrt elevator angle delta(de)
                'Cmdf' : -0.002,                  # pitching moment coefficent (Cm) wrt flap deflection (df)
                'Cnb' : 0.11,                     # yawing moment coefficent (Cn) wrt beta (b)
                'Cnbd' : 0.0015,                  # yawing moment coefficent (Cn) wrt beta dot(bd)
                'Cnp' : -0.0753,                  # yawing moment coefficent (Cn) wrt pitch rate (p)
                'Cnr' : -0.1210,                  # yawing moment coefficent (Cn) wrt yaw rate (r)
                'Cnda' : -0.000042 * 100,         # yawing moment coefficent (Cn) wrt aileron deflection (da)
                'Cndr' : 0.0014 * 100,            # yawing moment coefficent (Cn) wrt rudder deflection (dr)
                'Clb' : -0.0838,                  # rolling moment coefficent (Cl) wrt beta(b)
                'Clbd' : -0.0004,                 # rolling moment coefficent (Cl) wrt beta dot(bd)
                'Clp' : -0.5050,                  # rolling moment coefficent (Cl) wrt roll rate (p)
                'Clr' : 0.2608 * (1/3),           # rolling moment coefficent (Cl) wrt yaw rate (r)
                'Clda' : -0.003 * 100,            # rolling moment coefficent (Cl) wrt aileron deflection (da)
                'Cldr' : -0.0002                  # rolling moment coefficent (Cl) wrt rudder deflection (dr)
            }
        }
        },
    'J400_Mars': {
        'FlightData' : {
            'Inertial' : {
                'g'   : 9.81 * (3.71/9.81),
                'm'   : 700 ,
                'Ixx' : 734.8230 ,
                'Iyy' : 1022.700 ,
                'Izz' : 1067.100,
                'Ixz' : 31.7037 * 0
            },
            'Geometric' : {
                'S' : 8.019,
                'c' : 0.99,
                'b' : 8.10
            },
            'Propeller' : {
                'P_max' : 89484 * (0.13/100),
                'eta'   : 0.8
            },
            'ControlLimits' : {                 # -18.4 (but this need a calibration correction: 18.4 * 2 - 9.0 = 27.8)
                'Lower' : np.array([0, np.deg2rad(-27.8), np.deg2rad(-24), np.deg2rad(-23.21), np.deg2rad(-40)]), # dT, da, de, dr, df (Likely wrong label*) # Richard think it should actually be [dT, de, da, dr, df] instead *
                'Upper' : np.array([1, np.deg2rad(9.0) , np.deg2rad(24) , np.deg2rad(23.21), np.deg2rad(-0)])  # dT, da, de, dr, df (Likely wrong label*) # Richard think it should actually be [dT, de, da, dr, df] instead *
            },
            'Aero' : {
                'alpha_o' : -0.0698,              # zero lift (o) alpha
                'Cdo' : 0.0061 * 8 * (10/9.81),               # zero lift (o) drag coefficient (Cd)     # 0.049745
                'k' : 0.0483 * (0.09/0.0483)* (10/9.81),     # induced drag coefficient (k)             # 0.0917
                'CLa' : 4.9635,                   # Lift coefficient (CL) wrt alpha (a)
                'CLq' : 9.8197 * (8/9.8197),                   # Lift coefficient (CL) wrt pitch rate (q)
                'CLad' : -1.987,                  # Lift coefficient (CL) wrt alpha dot (ad)
                'CLde' : 0.0007976 * 1000,        # Lift coefficent (CL) wrt elevator deflection (de)
                'CLdf' : -1.75, # 0.0229 * 10,             # Lift coefficent (CL) wrt flap deflection (df)
                'CLo' : 0.36486, #* (0.8254/0.36486),                  # lift coefficent (CL) at zero alpha (o)
                'Cyb' : -0.23,                    # side force coefficent (Cy) wrt beta(b)
                'Cybd' : -0.0032,                 # side force coefficent (Cy) wrt beta dot(bd)
                'Cyp' : 0.1624,                   # side force coefficent (Cy) wrt roll rate(p)
                'Cyr' : 0.1889,                   # side force coefficent (Cy) wrt yaw rate(r)
                'Cyda' : 0.00033,                 # side force coefficent (Cy) wrt aileron deflection(da)
                'Cydr' : -0.0031 * (-10),         # side force coefficent (Cy) wrt rudder deflection(dr)
                'Cmo' : 0.07575,                  # pitching moment coefficent (Cm) at zero alpha (o)
                'Cma' : -2.6001,                  # pitching moment coefficent (Cm) wrt alpha (a)
                'Cmq' : -23.6813,                 # pitching moment coefficent (Cm) wrt pitch rate (p)
                'Cmad' : -2.210,                  # pitching moment coefficent (Cm) wrt alpha dot(ad)
                'Cmde' : -0.036 * 100,            # pitching moment coefficent (Cm) wrt elevator angle delta(de)
                'Cmdf' : -0.002,                  # pitching moment coefficent (Cm) wrt flap deflection (df)
                'Cnb' : 0.11,                     # yawing moment coefficent (Cn) wrt beta (b)
                'Cnbd' : 0.0015,                  # yawing moment coefficent (Cn) wrt beta dot(bd)
                'Cnp' : -0.0753,                  # yawing moment coefficent (Cn) wrt pitch rate (p)
                'Cnr' : -0.1210,                  # yawing moment coefficent (Cn) wrt yaw rate (r)
                'Cnda' : -0.000042 * 100,         # yawing moment coefficent (Cn) wrt aileron deflection (da)
                'Cndr' : 0.0014 * 100,            # yawing moment coefficent (Cn) wrt rudder deflection (dr)
                'Clb' : -0.0838,                  # rolling moment coefficent (Cl) wrt beta(b)
                'Clbd' : -0.0004,                 # rolling moment coefficent (Cl) wrt beta dot(bd)
                'Clp' : -0.5050,                  # rolling moment coefficent (Cl) wrt roll rate (p)
                'Clr' : 0.2608 * (1/3),           # rolling moment coefficent (Cl) wrt yaw rate (r)
                'Clda' : -0.003 * 100,            # rolling moment coefficent (Cl) wrt aileron deflection (da)
                'Cldr' : -0.0002                  # rolling moment coefficent (Cl) wrt rudder deflection (dr)
            }
        }
        },    
    'PC-9': {
        'FlightData' : {
            'Inertial' : {
                'g'   : 9.81 ,
                'm'   : 2087 ,
                'Ixx' : 5066 ,
                'Iyy' : 6578 ,
                'Izz' : 10975,
                'Ixz' : 203
            },
            'Geometric' : {
                'S' : 16.29,
                'c' : 1.652,
                'b' : 10.12
            },
            'Propeller' : {
                'P_max' : 950000,
                'eta'   : 0.8
            },
            'ControlLimits' : {
                'Lower' : np.array([0, np.deg2rad(-15), np.deg2rad(-20), np.deg2rad(-30), np.deg2rad(-45)]), # dT, da, de, dr, df
                'Upper' : np.array([1, np.deg2rad(20) , np.deg2rad(25) , np.deg2rad(30),  np.deg2rad(-0) ])  # dT, da, de, dr, df
            },
            'Aero' : {
                'alpha_o' : -3.0 / 57.3,
                'Cdo' : 0.020,
                'k' : 0.050,
                'CLa' : 5.827,
                'CLq' : 7.960,
                'CLad' : -1.987,
                'CLde' : 0.532,
                'CLdf' : 0,
                'CLo' : -5.827 * (-3.0 / 57.3),  # defined as: CLo = -CLa * alpha_o 
                'Cyb' : -0.507,
                'Cybd' : -0.0032,
                'Cyp' : -0.128,
                'Cyr' : 0.336,
                'Cyda' : 0.000,
                'Cydr' : 0.050,
                'Cmo' : 0.06,
                'Cma' : -0.802,
                'Cmq' : -17.72,
                'Cmad' : -2.210,
                'Cmde' : -1.822,
                'Cmdf' : 0,
                'Cnb' : 0.107,
                'Cnbd' : 0.0015,
                'Cnp' : -0.0226,
                'Cnr' : -0.160,
                'Cnda' : 0.0048,
                'Cndr' : -0.115,
                'Clb' : -0.0852,
                'Clbd' : -0.0004,
                'Clp' : -0.328,
                'Clr' : 0.0776,
                'Clda' : -0.164,
                'Cldr' : 0.0302
            }
        }
        },
    'Aircraft1': {
        'FlightData' : {
            'Inertial' : {
                'g'   : 9.81 ,
                'm'   : 2087 ,
                'Ixx' : 5066 ,
                'Iyy' : 6578 ,
                'Izz' : 10975,
                'Ixz' : 203
            },
            'Geometric' : {
                'S' : 16.29,
                'c' : 1.652,
                'b' : 10.119
            },
            'Propeller' : {
                'P_max' : 950000,
                'eta'   : 0.8
            },
            'ControlLimits' : {
                'Lower' : np.array([0, np.deg2rad(-25), np.deg2rad(-25), np.deg2rad(-25),np.deg2rad(-45)]), # dT, da, de, dr, df
                'Upper' : np.array([1, np.deg2rad(25) , np.deg2rad(25) , np.deg2rad(25), np.deg2rad(-0) ])  # dT, da, de, dr, df
            },
            'Aero' : {
                'alpha_o' : -3.0 / 57.3,
                'Cdo' : 0.020,
                'k' : 0.050,
                'CLa' : 5.827,
                'CLq' : 7.960,
                'CLad' : -1.987,
                'CLde' : 0.532,
                'CLdf' : 0.5,
                'CLo' : -5.827 * (-3.0 / 57.3),  # defined as: CLo = -CLa * alpha_o 
                'Cyb' : -0.507,
                'Cybd' : -0.0032,
                'Cyp' : -0.128,
                'Cyr' : 0.336,
                'Cyda' : 0.000,
                'Cydr' : 0.050,
                'Cmo' : 0.06,
                'Cma' : -0.802,
                'Cmq' : -17.72,
                'Cmad' : -2.210,
                'Cmde' : -1.822,
                'Cmdf' : -0.25,
                'Cnb' : 0.107,
                'Cnbd' : 0.0015,
                'Cnp' : -0.0226,
                'Cnr' : -0.160,
                'Cnda' : 0.0048,
                'Cndr' : -0.115,
                'Clb' : -0.0852,
                'Clbd' : -0.0004,
                'Clp' : -0.328,
                'Clr' : 0.0776,
                'Clda' : -0.164,
                'Cldr' : 0.0302
            }
        }
        },
    'Aircraft2': {
        'FlightData' : {
            'Inertial' : {
                'g'   : 9.81 ,
                'm'   : 1700 ,
                'Ixx' : 4000 ,
                'Iyy' : 5000 ,
                'Izz' : 9000 ,
                'Ixz' : 180 ,
            },
            'Geometric' : {
                'S' : 16.29 ,
                'c' : 1.3 ,
                'b' : 8.12 
            },
            'Propeller' : {
                'P_max' : 950000 ,
                'eta'   : 0.8 
            },
            'ControlLimits' : {
                'Lower' : np.array([0, np.deg2rad(-25), np.deg2rad(-25), np.deg2rad(-25),np.deg2rad(-45)]), # dT, da, de, dr, df
                'Upper' : np.array([1, np.deg2rad(25) , np.deg2rad(25) , np.deg2rad(25), np.deg2rad(-0) ])  # dT, da, de, dr, df
            },
            'Aero' : {
                'alpha_o' : -0.052356020942408 ,
                'Cdo'  :  0.019 ,
                'k'    :  0.05 ,
                'CLa'  :  5.927 ,
                'CLq'  :  8.012 ,
                'CLad' : -1.987 ,
                'CLde' :  0.632 ,
                'CLdf' :  0.6 ,
                'CLo'  :  0.310314136125654 , 
                'Cyb'  : -0.107 ,
                'Cybd' : -0.0032 ,
                'Cyp'  : -0.178 ,
                'Cyr'  :  0.436 ,
                'Cyda' :  0 ,
                'Cydr' :  0.05 ,
                'Cmo'  :  0.06 ,
                'Cma'  : -0.802 ,
                'Cmq'  : -17.72 ,
                'Cmad' : -2.21 ,
                'Cmde' : -1.822 ,
                'Cmdf' : -0.25 ,
                'Cnb'  :  0.127 ,
                'Cnbd' :  0.0015 ,
                'Cnp'  : -0.0226 ,
                'Cnr'  : -0.16 ,
                'Cnda' :  0.0048 ,
                'Cndr' : -0.115 ,
                'Clb'  : -0.0852 ,
                'Clbd' : -0.0004 ,
                'Clp'  : -0.428 ,
                'Clr'  :  0.0776 ,
                'Clda' : -0.164 ,
                'Cldr' :  0.0302 
            }
        }
        },
    'Aircraft3': {
        'FlightData' : {
            'Inertial' : {
                'g'   : 9.81 ,
                'm'   : 968 ,
                'Ixx' : 792.45 ,
                'Iyy' : 1527.46 ,
                'Izz' : 2070.81 ,
                'Ixz' : 37.75 ,
            },
            'Geometric' : {
                'S' : 9.5 ,
                'c' : 1.231 ,
                'b' : 8 
            },
            'Propeller' : {
                'P_max' : 260960 ,
                'eta'   : 0.82 
            },
            'ControlLimits' : {
                'Lower' : np.array([0, np.deg2rad(-25), np.deg2rad(-25), np.deg2rad(-25),np.deg2rad(-45)]), # dT, da, de, dr, df
                'Upper' : np.array([1, np.deg2rad(25) , np.deg2rad(25) , np.deg2rad(25), np.deg2rad(-0) ])  # dT, da, de, dr, df
            },
            'Aero' : {
                'alpha_o' : -0.04612 ,
                'Cdo'  :  0.0213 ,
                'k'    :  0.0888 ,
                'CLa'  :  5.6488 ,
                'CLq'  :  9.0662 ,
                'CLad' :  0 ,
                'CLde' :  0.753 ,
                'CLdf' :  0.5 ,
                'CLo'  :  0.26053 ,
                'Cyb'  : -0.3947 ,
                'Cybd' :  0 ,
                'Cyp'  :  0.0166 ,
                'Cyr'  :  0.1858 ,
                'Cyda' : -0.009053 ,
                'Cydr' :  0.2347 ,
                'Cmo'  :  0.04418 ,
                'Cma'  : -0.2937 ,
                'Cmq'  : -19.6161 ,
                'Cmad' :  0 ,
                'Cmde' : -2.2203 ,
                'Cmdf' : -0.25 ,
                'Cnb'  :  0.07082 ,
                'Cnbd' :  0 ,
                'Cnp'  : -0.02447 ,
                'Cnr'  : -0.1384 ,
                'Cnda' : -0.0106 ,
                'Cndr' : -0.1097 ,
                'Clb'  : -0.01179 ,
                'Clbd' :  0 ,
                'Clp'  : -0.4549 ,
                'Clr'  :  0.06467 ,
                'Clda' : -0.3545 ,
                'Cldr' :  0.009856 
            }
        }
        },
    'Aircraft4': {
        'FlightData' : {
            'Inertial' : {
                'g'   : 9.81 ,
                'm'   : 750 ,
                'Ixx' : 605.29 ,
                'Iyy' : 1766.52 ,
                'Izz' : 2168.83 ,
                'Ixz' : 84.71 ,
            },
            'Geometric' : {
                'S' : 11.9 ,
                'c' : 0.85 ,
                'b' : 7 
            },
            'Propeller' : {
                'P_max' : 253504 ,
                'eta'   : 0.85 
            },
            'ControlLimits' : {
                'Lower' : np.array([0, np.deg2rad(-25), np.deg2rad(-25), np.deg2rad(-25),np.deg2rad(-45)]), # dT, da, de, dr, df
                'Upper' : np.array([1, np.deg2rad(25) , np.deg2rad(25) , np.deg2rad(25), np.deg2rad(-0) ])  # dT, da, de, dr, df
            },
            'Aero' : {
                'alpha_o' : -0.013772916726263 ,
                'Cdo'  :  0.020342 ,
                'k'    :  0.063776 ,
                'CLa'  :  4.7860 ,
                'CLq'  :  9.2690 ,
                'CLad' : -1.1439 ,
                'CLde' :  0.61398157326219 ,
                'CLdf' :  0.5 ,
                'CLo'  :  0.065917179451894 , 
                'Cyb'  : -0.379162 ,
                'Cybd' :  0.117 ,
                'Cyp'  :  0.016304 ,
                'Cyr'  :  0.1222 ,
                'Cyda' :  0.036153636872755 ,
                'Cydr' :  0.202941651035338 ,
                'Cmo'  :  0.03337 ,
                'Cma'  : -0.295926 ,
                'Cmq'  : -25.094866 ,
                'Cmad' :  0.9531 ,
                'Cmde' : -2.394161442733658 ,
                'Cmdf' : -0.25 ,
                'Cnb'  :  0.027705 ,
                'Cnbd' : -0.0051 ,
                'Cnp'  : -0.000421 ,
                'Cnr'  : -0.140618 ,
                'Cnda' : -0.015126085791454 ,
                'Cndr' : -0.112299727845641 ,
                'Clb'  : -0.024446 ,
                'Clbd' :  0.0068 ,
                'Clp'  : -0.456441 ,
                'Clr'  :  0.0601 ,
                'Clda' : -0.357353776823094 ,
                'Cldr' :  0.000394 
            }
        }
        },
    'Aircraft5': {
        'FlightData' : {
            'Inertial' : {
                'g'   : 9.81 ,
                'm'   : 2000 ,
                'Ixx' : 500.45 ,
                'Iyy' : 600.46 ,
                'Izz' : 800.81 ,
                'Ixz' : 25.75 ,
            },
            'Geometric' : {
                'S' : 9.5 ,
                'c' : 1.231 ,
                'b' : 8 
            },
            'Propeller' : {
                'P_max' : 260960 ,
                'eta'   : 0.82 
            },
            'ControlLimits' : {
                'Lower' : np.array([0, np.deg2rad(-25), np.deg2rad(-25), np.deg2rad(-25),np.deg2rad(-45)]), # dT, da, de, dr, df
                'Upper' : np.array([1, np.deg2rad(25) , np.deg2rad(25) , np.deg2rad(25), np.deg2rad(-0) ])  # dT, da, de, dr, df
            },
            'Aero' : {
                'alpha_o' : -0.04612 ,
                'Cdo'  :  0.0193 ,
                'k'    :  0.0888 ,
                'CLa'  :  5.9488 ,
                'CLq'  :  9.0662 ,
                'CLad' :  0 ,
                'CLde' :  0.753 ,
                'CLdf' :  0.5 ,
                'CLo'  :  0.26053 , 
                'Cyb'  : -0.6947 ,
                'Cybd' :  0 ,
                'Cyp'  :  0.0376 ,
                'Cyr'  :  0.3558 ,
                'Cyda' : -0.012053 ,
                'Cydr' :  0.2947 ,
                'Cmo'  :  0.04418 ,
                'Cma'  : -0.2937,
                'Cmq'  : -19.6161 ,
                'Cmad' :  0 ,
                'Cmde' : -2.2203 ,
                'Cmdf' : -0.25 ,
                'Cnb'  :  0.1 ,
                'Cnbd' :  0 ,
                'Cnp'  : -0.02447 ,
                'Cnr'  : -0.06 ,
                'Cnda' : -0.0426 ,
                'Cndr' : -0.0697 ,
                'Clb'  : -0.01179 ,
                'Clbd' :  0 ,
                'Clp'  : -0.8549 ,
                'Clr'  :  0.05467 ,
                'Clda' : -0.3545 ,
                'Cldr' :  0.009856 
            }
        }
        },
    'Aircraft6': {
        'FlightData' : {
            'Inertial' : {
                'g'   : 9.81 ,
                'm'   : 1200 ,
                'Ixx' : 500.29 ,
                'Iyy' : 1250.52 ,
                'Izz' : 1750.83 ,
                'Ixz' : 75.71 ,
            },
            'Geometric' : {
                'S' : 10.9 ,
                'c' : 0.85 ,
                'b' : 6 
            },
            'Propeller' : {
                'P_max' : 253504 ,
                'eta'   : 0.85 
            },
            'ControlLimits' : {
                'Lower' : np.array([0, np.deg2rad(-25), np.deg2rad(-25), np.deg2rad(-25),np.deg2rad(-45)]), # dT, da, de, dr, df
                'Upper' : np.array([1, np.deg2rad(25) , np.deg2rad(25) , np.deg2rad(25), np.deg2rad(-0) ])  # dT, da, de, dr, df
            },
            'Aero' : {
                'alpha_o' : -0.013772916726263 ,
                'Cdo'  :  0.020342 ,
                'k'    :  0.063776 ,
                'CLa'  :  4.786 ,
                'CLq'  :  9.269 ,
                'CLad' : -1.1439 ,
                'CLde' :  0.613981573262190 ,
                'CLdf' :  0.5 ,
                'CLo'  :  0.065917179451894 ,
                'Cyb'  : -0.179162 ,
                'Cybd' :  0.07117 ,
                'Cyp'  :  0.016304 ,
                'Cyr'  :  0.1222 ,
                'Cyda' :  0.036153636872755 ,
                'Cydr' :  0.202941651035338 ,
                'Cmo'  :  0.03337 ,
                'Cma'  : -0.295926 ,
                'Cmq'  : -25.094866 ,
                'Cmad' :  0.9531 ,
                'Cmde' : -2.394161442733658 ,
                'Cmdf' : -0.25 ,
                'Cnb'  :  0.0207705 ,
                'Cnbd' : -0.0551 ,
                'Cnp'  : -0.000421 ,
                'Cnr'  : -0.0650618 ,
                'Cnda' : -0.015126085791454 ,
                'Cndr' : -0.112299727845641 ,
                'Clb'  : -0.024446 ,
                'Clbd' :  0.0068 ,
                'Clp'  : -0.456441 ,
                'Clr'  :  0.0601 ,
                'Clda' : -0.357353776823094 ,
                'Cldr' :  0.000394 
            }
        }
        },     
}
