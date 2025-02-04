import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# This should be in the initialization phase of the simulation to save run time work

altitude_points = np.array([-1.00000000e+02,  0.00000000e+00,  1.00000000e+00,  5.00000000e+01,
                            1.00000000e+02,  1.50000000e+02,  2.00000000e+02,  2.50000000e+02,
                            3.04800000e+02,  6.09600000e+02,  1.17199466e+03,  2.34398932e+03,
                            4.55775701e+03,  1.05805073e+04,  1.38034927e+04,  1.68311455e+04,
                            1.97936876e+04,  4.57200000e+04])
stan_dev_u_points = np.array([0, 0, 1.7300465747703184, 1.3874462887500087, 
                                1.201581198127368, 1.081200335621849, 0.9946432929792891, 0.9283228269822059, 
                                0.87072649, 1.62608392, 1.8188234167332495, 1.7304083895309388, 
                                1.376748280721696, 0.8652041947654694, 0.7199509357902445, 0.4547058541833124, 
                                0.0, 0.0])
stan_dev_v_points = np.array([0, 0, 1.7300320302943264, 1.387456241528783, 
                                1.2015926104134622, 1.0812116077686194, 0.9946540750850655, 0.9283330606659893, 
                                0.87073615, 1.62611592, 1.8188592053593668, 1.730442438432175, 
                                1.3767753707234096, 0.8652212192160875, 0.7199651021214158, 0.4547148013398417, 
                                0.0, 0.0])
stan_dev_w_points = np.array([0, 0, 0.8705543858535112, 0.8707274776263864, 
                            0.870727643480018, 0.87072662117927, 0.8707253018381317, 0.8707238636791084, 
                            0.87072222, 1.62611592, 1.8188592053593668, 1.730442438432175, 
                            1.3767753707234096, 0.8652212192160875, 0.7199651021214158, 0.4547148013398417, 
                            0.0, 0.0])
stan_dev_p_points = np.array([0, 0, 0.20600934850341626, 0.05591954047609127, 
                            0.0443833686886628, 0.038772443845604534, 0.03522710307313788, 0.03270194565296064, 
                            0.03061131, 0.04744, 0.053063057949548145, 0.050483603743667324, 
                            0.040165786920144074, 0.025241801871833662, 0.02100412710502947, 0.013265764487387036, 
                            0.0, 0.0])
stan_dev_q_points = np.array([0, 0, 0.08235526487557372, 0.04152173303832524, 
                            0.03111987019699544, 0.025943871161635775, 0.022709870650808452, 0.020445225731463526, 
                            0.01860422, 0.02651186, 0.029654304407694003, 0.028212775721208873, 
                            0.022446660975268377, 0.014106387860604436, 0.011738162161378874, 0.007413576101923501, 
                            0.0, 0.0])
stan_dev_r_points = np.array([0, 0, 0.17632270766234787, 0.041899603356915, 
                            0.03201676805602458, 0.027602966162907725, 0.024939226085019205, 0.023095573241348973, 
                            0.02159947, 0.03071001, 0.03435006670116053, 0.03268027179207634, 
                            0.02600109215573957, 0.01634013589603817, 0.013596901402542707, 0.008587516675290133, 
                            0.0, 0.0])

# Create interpolation functions, these should be global variables in the simulation
interpolate_stan_dev_u = interp1d(altitude_points, stan_dev_u_points, kind='linear', fill_value="extrapolate")
interpolate_stan_dev_v = interp1d(altitude_points, stan_dev_v_points, kind='linear', fill_value="extrapolate")
interpolate_stan_dev_w = interp1d(altitude_points, stan_dev_w_points, kind='linear', fill_value="extrapolate")
interpolate_stan_dev_p = interp1d(altitude_points, stan_dev_p_points, kind='linear', fill_value="extrapolate")
interpolate_stan_dev_q = interp1d(altitude_points, stan_dev_q_points, kind='linear', fill_value="extrapolate")
interpolate_stan_dev_r = interp1d(altitude_points, stan_dev_r_points, kind='linear', fill_value="extrapolate")


def calculate_std_vs_altitude(altitudes):
    std_u = []
    std_v = []
    std_w = []
    std_p = []
    std_q = []
    std_r = []

    for altitude in altitudes:

        stan_dev_u = interpolate_stan_dev_u(altitude)
        stan_dev_v = interpolate_stan_dev_v(altitude)
        stan_dev_w = interpolate_stan_dev_w(altitude)
        stan_dev_p = interpolate_stan_dev_p(altitude)
        stan_dev_q = interpolate_stan_dev_q(altitude)
        stan_dev_r = interpolate_stan_dev_r(altitude)


        std_u.append(stan_dev_u)
        std_v.append(stan_dev_v)
        std_w.append(stan_dev_w)
        std_p.append(stan_dev_p)
        std_q.append(stan_dev_q)
        std_r.append(stan_dev_r)

    return std_u, std_v, std_w, std_p, std_q, std_r


altitudes = np.linspace(-100, 20000, 200)  # Example altitudes in meters

# Calculate standard deviations for each altitude
std_u, std_v, std_w, std_p, std_q, std_r = calculate_std_vs_altitude(altitudes)


# User input
whether_to_plot_altitude_dependence = 1

if whether_to_plot_altitude_dependence == 1:

    # Plotting
    plt.figure()
    plt.plot(altitudes, std_u, label='std_u')
    plt.plot(altitudes, std_v, '--', label='std_v')
    plt.plot(altitudes, std_w, ':', label='std_w')
    # plt.plot(altitudes, std_p, label='std_p')
    # plt.plot(altitudes, std_q, label='std_q')
    # plt.plot(altitudes, std_r, label='std_r')
    plt.xlabel('Altitude [m]')
    plt.xlim(-100, 20000)    
    plt.ylabel('Standard Deviation')
    plt.ylim(0, 2.0)
    plt.title(r'Standard Deviations of $u_g$, $v_g$, $w_g$ vs Altitude')
    plt.legend()
    plt.grid(True)
    plt.show()

    # # Plotting
    # plt.figure()
    # # plt.plot(altitudes, std_u, label='std_u')
    # # plt.plot(altitudes, std_v, '--', label='std_v')
    # plt.plot(altitudes, std_w, label='std_w')
    # # plt.plot(altitudes, std_p, label='std_p')
    # # plt.plot(altitudes, std_q, label='std_q')
    # # plt.plot(altitudes, std_r, label='std_r')
    # plt.xlabel('Altitude [m]')
    # plt.ylabel('Standard Deviation')
    # plt.title(r'Standard Deviations of $w_g$ vs Altitude')
    # plt.legend()
    # plt.grid(True)
    # plt.show()

    # Plotting
    plt.figure()
    # plt.plot(altitudes, std_u, label='std_u')
    # plt.plot(altitudes, std_v, '--', label='std_v')
    # plt.plot(altitudes, std_w, label='std_w')
    plt.plot(altitudes, std_p, label='std_p')
    plt.plot(altitudes, std_q, label='std_q')
    plt.plot(altitudes, std_r, label='std_r')
    plt.xlim(-100, 20000)    
    plt.xlabel('Altitude [m]')
    plt.ylabel('Standard Deviation')
    plt.title(r'Standard Deviations of $p_g$, $q_g$, $r_g$ vs Altitude')
    plt.legend()
    plt.grid(True)
    plt.show()