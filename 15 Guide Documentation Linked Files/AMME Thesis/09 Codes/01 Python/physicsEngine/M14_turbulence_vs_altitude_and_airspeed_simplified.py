import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Function: Convert [kn] to [m/s]
def kn2mps(kn):
    return kn / 1.9438452

# Function: Get aircraft wing span and trim velocity (b, V)
def geometry_and_trim():
    # FlightData = aero4560_LoadFlightData()  # Implement your aircraft data loading
    b = 8.10 # m
    # b = 10.2

    X0 = np.array([4.28509979e+01, -3.34313946e-01, 1.16436917e+00, -1.69348740e-03, -2.52844822e-02, 2.04244591e-02, -3.38758335e-02, -1.34899263e-01, 2.92678634e+00, 0.00000000e+00, 0.00000000e+00, -6.40350647e+02])
    V = np.sqrt(X0[0]**2 + X0[1]**2 + X0[2]**2)
    # print(V)
    # print(-X0[11])
    return b, V

# Function: Get Turbulent Scales & Intensities
def turb_scale_inten(X0):

    # Acquire altitude
    alt = -X0[11]

    if alt <= 304.8: # m = (1000 ft)

        # Typically, at 20 feet (6 meters) the wind speed is 15 knots in light turbulence, 
        #                                                    30 knots in moderate turbulence, 
        #                                                and 45 knots for severe turbulence.
        # u_20 is the wind speed at 20 feet (6 meters).

        u_20_in_kn = 30  # Reference mean wind speed in knots (moderate intensity)
        u_20 = kn2mps(u_20_in_kn)  # Convert to m/s

        h = -X0[11]  # Trimmed flight altitude
        L_u = h / (0.177 + 0.0027 * h)**1.2
        L_v = h / (0.177 + 0.0027 * h)**1.2
        L_w = h

        sigma_u = (0.1 * u_20) * (1 / (0.177 + 0.0027 * h)**0.4)
        sigma_v = (0.1 * u_20) * (1 / (0.177 + 0.0027 * h)**0.4)
        sigma_w = 0.1 * u_20

    elif alt >= 609.6: # m = (2000 ft)

        L_u = 0.3048 * 1750 # m
        L_v = (0.3048 * 1750) # m
        L_w = (0.3048 * 1750) # m

        # User input
        which_turbulence_category_to_use = 3 # Category: 2 is light, 3 is moderate, 5 is severe

        if which_turbulence_category_to_use == 1:

            altitude_1 = [0.5340453938584779, 1.9225634178905207, 4.913217623497998, 7.369826435246996, 150.00]
            sigma_1 = [3.1951731374606505, 2.093389296956978, 1.0650577124868836, 0.00, 0.00]
            altitude_table = np.array(altitude_1) * (1000 * 0.3048)
            sigma_table = np.array(sigma_1) * (0.3048)

        elif which_turbulence_category_to_use == 2:

            altitude_light_2 = [0.6408544726301736, 3.7383177570093458, 7.369826435246996, 14.739652870493993, 35.24699599465955, 43.89853137516689, 150.00]
            sigma_light_2 = [6.610703043022036, 7.381951731374607, 6.720881427072403, 4.664218258132214, 0.40398740818468, 0.00, 0.00]
            altitude_table = np.array(altitude_light_2) * (1000 * 0.3048)
            sigma_table = np.array(sigma_light_2) * (0.3048)

        elif which_turbulence_category_to_use == 3:

            altitude_moderate_3 = [0.6408544726301736, 3.8451268357810413, 7.690253671562083, 14.953271028037383, 34.71295060080107, 45.287049399198935, 55.22029372496662, 64.93991989319092, 150.00]
            sigma_moderate_3 = [8.630640083945435, 10.577124868835257, 10.06295907660021, 8.006295907660022, 5.031479538300105, 4.186778593913956, 2.6442812172088144, 0.00, 0.00]
            altitude_table = np.array(altitude_moderate_3) * (1000 * 0.3048)
            sigma_table = np.array(sigma_moderate_3) * (0.3048)

        elif which_turbulence_category_to_use == 4:

            altitude_moderate_4 = [0.5340453938584779, 1.8157543391188251, 3.7383177570093458, 7.4766355140186915, 14.846461949265688, 34.92656875834446, 55.11348464619493, 65.04672897196262, 80.00, 90.00, 150.00]
            sigma_moderate_4 = [11.752360965372509, 12.964323189926548, 16.012591815320043, 15.094438614900316, 11.605456453305353, 8.0797481636936, 7.932843651626443, 4.921301154249738, 2.1301154249737673, 0.00, 0.00]
            altitude_table = np.array(altitude_moderate_4) * (1000 * 0.3048)
            sigma_table = np.array(sigma_moderate_4) * (0.3048)

        elif which_turbulence_category_to_use == 5:

            altitude_severe_5 = [0.6408544726301736, 3.7383177570093458, 7.2630173564753004, 24.99332443257677, 35.14018691588785, 44.64619492656876, 65.26034712950602, 80.00, 100.00, 150.00]
            sigma_severe_5 = [15.645330535152151, 22.95383001049318, 23.614900314795385, 20.05246589716684, 16.012591815320043, 15.204616998950682, 7.859391395592865, 5.141657922350472, 0.00, 0.00]
            altitude_table = np.array(altitude_severe_5) * (1000 * 0.3048)
            sigma_table = np.array(sigma_severe_5) * (0.3048)

        elif which_turbulence_category_to_use == 6:

            altitude_6 = [0.5340453938584779, 3.6315086782376502, 7.369826435246996, 24.886515353805073, 35.033377837116156, 45.07343124165554, 65.04672897196262, 79.89319092122831, 110.00, 150.00]
            sigma_6 = [18.803777544596013, 28.389296956977965, 30.22560335781742, 31.070304302203567, 25.267576075550892, 23.100734522560337, 10.724029380902413, 7.271773347324239, 0.00, 0.00]
            altitude_table = np.array(altitude_6) * (1000 * 0.3048)
            sigma_table = np.array(sigma_6) * (0.3048)

        # Create an interpolation function
        interpolate_sigma = interp1d(altitude_table, sigma_table, kind='linear', fill_value="extrapolate")
        sigma_values = interpolate_sigma(alt)

        sigma_u = sigma_values
        sigma_v = sigma_values
        sigma_w = sigma_values

    else: # between 1000 ft and 2000 ft

        # The values below are not used since u_g, v_g, w_g, p_g, q_g, r_g will be interpolated instead

        L_u = 0.0001
        L_v = 0.0001
        L_w = 0.0001

        sigma_u = 0.0001
        sigma_v = 0.0001
        sigma_w = 0.0001     
    
    return L_u, L_v, L_w, sigma_u, sigma_v, sigma_w

# Function: Generate the Dryden Spectra PSDs
def dryden_psd(omega, V, b, L_u, L_v, L_w, sigma_u, sigma_v, sigma_w, plot_psd=0):

    Phi_u_g = sigma_u**2 * (2 * L_u / (np.pi * V)) * (1 / (1 + (L_u / V * omega)**2))
    Phi_v_g = sigma_v**2 * (L_v / (np.pi * V)) * (1 + 3 * (L_v / V * omega)**2) / ((1 + (L_v / V * omega)**2)**2)
    Phi_w_g = sigma_w**2 * (L_w / (np.pi * V)) * (1 + 3 * (L_w / V * omega)**2) / ((1 + (L_w / V * omega)**2)**2)
    Phi_p_g = 0.8 * sigma_w**2 / (L_w * V) * ((np.pi * L_w / (4 * b))**(1/3)) / (1 + ((4 * b / (np.pi * V)) * omega)**2)
    Phi_q_g = ((omega / V)**2 / (1 + (4 * b / (np.pi * V) * omega)**2)) * Phi_w_g
    Phi_r_g = ((omega / V)**2 / (1 + (3 * b / (np.pi * V) * omega)**2)) * Phi_v_g

    if plot_psd:
        plt.figure()
        plt.loglog(omega, Phi_u_g, label=r'$\Phi_{u_g}$')
        plt.loglog(omega, Phi_v_g, label=r'$\Phi_{v_g}$')
        plt.loglog(omega, Phi_w_g, label=r'$\Phi_{w_g}$')
        plt.loglog(omega, Phi_p_g, label=r'$\Phi_{p_g}$')
        plt.loglog(omega, Phi_q_g, label=r'$\Phi_{q_g}$')
        plt.loglog(omega, Phi_r_g, label=r'$\Phi_{r_g}$')
        plt.xlabel(r'$\omega$ [rad/s]')
        plt.ylabel(r'$\Phi$')
        plt.grid(True)
        plt.legend()
        plt.show()

    return Phi_u_g, Phi_v_g, Phi_w_g, Phi_p_g, Phi_q_g, Phi_r_g

##################################################################################################################################################################
##################################################################################################################################################################

# Function to calculate standard deviations as a function of altitude
def calculate_std_vs_altitude(altitudes):
    std_u = []
    std_v = []
    std_w = []
    std_p = []
    std_q = []
    std_r = []
    
    for altitude in altitudes:
        # Simulate trim condition for each altitude
        X0 = np.array([4.28509979e+01, -3.34313946e-01, 1.16436917e+00, 
                       -1.69348740e-03, -2.52844822e-02, 2.04244591e-02, 
                       -3.38758335e-02, -1.34899263e-01, 2.92678634e+00, 
                       0.00000000e+00, 0.00000000e+00, -altitude])


        if (0 < altitude < 304.8) or (altitude > 609.6):

            L_u, L_v, L_w, sigma_u, sigma_v, sigma_w = turb_scale_inten(X0)
            b, V = geometry_and_trim()

            omega = np.logspace(-5, 5, 10000)
            Phi_u_g, Phi_v_g, Phi_w_g, Phi_p_g, Phi_q_g, Phi_r_g = dryden_psd(omega, V, b, L_u, L_v, L_w, sigma_u, sigma_v, sigma_w, plot_psd=0)

            stan_dev_u = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_u_g, omega))
            stan_dev_v = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_v_g, omega))
            stan_dev_w = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_w_g, omega))
            stan_dev_p = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_p_g, omega))
            stan_dev_q = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_q_g, omega))
            stan_dev_r = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_r_g, omega))

        elif altitude <= 0:

            stan_dev_u = 0
            stan_dev_v = 0
            stan_dev_w = 0
            stan_dev_p = 0
            stan_dev_q = 0
            stan_dev_r = 0

        else:

            X0_at_1000ft = X0
            X0_at_1000ft[11] = -304.8

            L_u, L_v, L_w, sigma_u, sigma_v, sigma_w = turb_scale_inten(X0_at_1000ft)
            b, V = geometry_and_trim()

            omega = np.logspace(-5, 5, 10000)
            Phi_u_g, Phi_v_g, Phi_w_g, Phi_p_g, Phi_q_g, Phi_r_g = dryden_psd(omega, V, b, L_u, L_v, L_w, sigma_u, sigma_v, sigma_w, plot_psd=0)

            stan_dev_u_at_1000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_u_g, omega))
            stan_dev_v_at_1000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_v_g, omega))
            stan_dev_w_at_1000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_w_g, omega))
            stan_dev_p_at_1000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_p_g, omega))
            stan_dev_q_at_1000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_q_g, omega))
            stan_dev_r_at_1000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_r_g, omega))


            X0_at_2000ft = X0
            X0_at_2000ft[11] = -609.6

            L_u, L_v, L_w, sigma_u, sigma_v, sigma_w = turb_scale_inten(X0_at_2000ft)
            b, V = geometry_and_trim()

            omega = np.logspace(-5, 5, 10000)
            Phi_u_g, Phi_v_g, Phi_w_g, Phi_p_g, Phi_q_g, Phi_r_g = dryden_psd(omega, V, b, L_u, L_v, L_w, sigma_u, sigma_v, sigma_w, plot_psd=0)

            stan_dev_u_at_2000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_u_g, omega))
            stan_dev_v_at_2000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_v_g, omega))
            stan_dev_w_at_2000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_w_g, omega))
            stan_dev_p_at_2000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_p_g, omega))
            stan_dev_q_at_2000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_q_g, omega))
            stan_dev_r_at_2000ft = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_r_g, omega))

            altitude_points = np.array([304.8, 609.6])
            stan_dev_u_points = np.array([stan_dev_u_at_1000ft, stan_dev_u_at_2000ft])
            stan_dev_v_points = np.array([stan_dev_v_at_1000ft, stan_dev_v_at_2000ft])
            stan_dev_w_points = np.array([stan_dev_w_at_1000ft, stan_dev_w_at_2000ft])
            stan_dev_p_points = np.array([stan_dev_p_at_1000ft, stan_dev_p_at_2000ft])
            stan_dev_q_points = np.array([stan_dev_q_at_1000ft, stan_dev_q_at_2000ft])
            stan_dev_r_points = np.array([stan_dev_r_at_1000ft, stan_dev_r_at_2000ft])

            # Create an interpolation function
            interpolate_stan_dev_u = interp1d(altitude_points, stan_dev_u_points, kind='linear', fill_value="extrapolate")
            interpolate_stan_dev_v = interp1d(altitude_points, stan_dev_v_points, kind='linear', fill_value="extrapolate")
            interpolate_stan_dev_w = interp1d(altitude_points, stan_dev_w_points, kind='linear', fill_value="extrapolate")
            interpolate_stan_dev_p = interp1d(altitude_points, stan_dev_p_points, kind='linear', fill_value="extrapolate")
            interpolate_stan_dev_q = interp1d(altitude_points, stan_dev_q_points, kind='linear', fill_value="extrapolate")
            interpolate_stan_dev_r = interp1d(altitude_points, stan_dev_r_points, kind='linear', fill_value="extrapolate")
            # sigma_values = interpolate_sigma(alt)


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

# Altitude range to test
# altitudes = np.linspace(-100, 3000, 100)  # Example altitudes in meters
altitudes = np.array([-100, 0, 1, 50, 100, 150, 200, 250, 304.8, 609.6, 3.8451268357810413 * (1000 * 0.3048), 7.690253671562083 * (1000 * 0.3048), 14.953271028037383 * (1000 * 0.3048), 34.71295060080107 * (1000 * 0.3048), 45.287049399198935 * (1000 * 0.3048), 55.22029372496662 * (1000 * 0.3048), 64.93991989319092 * (1000 * 0.3048), 150.00 * (1000 * 0.3048)])  # Example altitudes in meters

# Calculate standard deviations for each altitude
std_u, std_v, std_w, std_p, std_q, std_r = calculate_std_vs_altitude(altitudes)


# Print individual variable tables to be used in flightModel.py 
print(f"------------------------------------------------------------------------------------------------------------------------------------------")
print(f"altitudes: {altitudes}")
# print(f"std_u: {std_u}")
# print(f"std_v: {std_v}")
# print(f"std_w: {std_w}")
# print(f"std_p: {std_p}")
# print(f"std_q: {std_q}")
# print(f"std_r: {std_r}")
altitudes = np.array([-1.00000000e+02,  0.00000000e+00,  1.00000000e+00,  5.00000000e+01,
                       1.00000000e+02,  1.50000000e+02,  2.00000000e+02,  2.50000000e+02,
                       3.04800000e+02,  6.09600000e+02,  1.17199466e+03,  2.34398932e+03,
                       4.55775701e+03,  1.05805073e+04,  1.38034927e+04,  1.68311455e+04,
                       1.97936876e+04,  4.57200000e+04])
std_u = np.array([0, 0, 1.7300465747703184, 1.3874462887500087, 
                 1.201581198127368, 1.081200335621849, 0.9946432929792891, 0.9283228269822059, 
                 0.87072649, 1.62608392, 1.8188234167332495, 1.7304083895309388, 
                 1.376748280721696, 0.8652041947654694, 0.7199509357902445, 0.4547058541833124, 
                 0.0, 0.0])
std_v = np.array([0, 0, 1.7300320302943264, 1.387456241528783, 
                  1.2015926104134622, 1.0812116077686194, 0.9946540750850655, 0.9283330606659893, 
                  0.87073615, 1.62611592, 1.8188592053593668, 1.730442438432175, 
                  1.3767753707234096, 0.8652212192160875, 0.7199651021214158, 0.4547148013398417, 
                  0.0, 0.0])
std_w = np.array([0, 0, 0.8705543858535112, 0.8707274776263864, 
                  0.870727643480018, 0.87072662117927, 0.8707253018381317, 0.8707238636791084, 
                  0.87072222, 1.62611592, 1.8188592053593668, 1.730442438432175, 
                  1.3767753707234096, 0.8652212192160875, 0.7199651021214158, 0.4547148013398417, 
                  0.0, 0.0])
std_p = np.array([0, 0, 0.20600934850341626, 0.05591954047609127, 
                  0.0443833686886628, 0.038772443845604534, 0.03522710307313788, 0.03270194565296064, 
                  0.03061131, 0.04744, 0.053063057949548145, 0.050483603743667324, 
                  0.040165786920144074, 0.025241801871833662, 0.02100412710502947, 0.013265764487387036, 
                  0.0, 0.0])
std_q = np.array([0, 0, 0.08235526487557372, 0.04152173303832524, 
                  0.03111987019699544, 0.025943871161635775, 0.022709870650808452, 0.020445225731463526, 
                  0.01860422, 0.02651186, 0.029654304407694003, 0.028212775721208873, 
                  0.022446660975268377, 0.014106387860604436, 0.011738162161378874, 0.007413576101923501, 
                  0.0, 0.0])
std_r = np.array([0, 0, 0.17632270766234787, 0.041899603356915, 
                  0.03201676805602458, 0.027602966162907725, 0.024939226085019205, 0.023095573241348973, 
                  0.02159947, 0.03071001, 0.03435006670116053, 0.03268027179207634, 
                  0.02600109215573957, 0.01634013589603817, 0.013596901402542707, 0.008587516675290133, 
                  0.0, 0.0])
print(f"==========================================================================================================================================")

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
    plt.xlim(-100, 10000)    
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
    plt.xlim(-100, 10000)    
    plt.xlabel('Altitude [m]')
    plt.ylabel('Standard Deviation')
    plt.title(r'Standard Deviations of $p_g$, $q_g$, $r_g$ vs Altitude')
    plt.legend()
    plt.grid(True)
    plt.show()

##################################################################################################################################################################
##################################################################################################################################################################

# Function to calculate standard deviations as a function of airspeed
def calculate_std_vs_airspeed(airspeeds):
    std_u = []
    std_v = []
    std_w = []
    std_p = []
    std_q = []
    std_r = []
    
    for V in airspeeds:
        # Simulate trim condition for each airspeed
        # X0 = np.array([V, -3.34313946e-01, 1.16436917e+00, 
        #                -1.69348740e-03, -2.52844822e-02, 2.04244591e-02, 
        #                -3.38758335e-02, -1.34899263e-01, 2.92678634e+00, 
        #                0.00000000e+00, 0.00000000e+00, -6.40350647e+02])
        X0 = np.array([V, -3.34313946e-01, 1.16436917e+00, 
                    -1.69348740e-03, -2.52844822e-02, 2.04244591e-02, 
                    -3.38758335e-02, -1.34899263e-01, 2.92678634e+00, 
                    0.00000000e+00, 0.00000000e+00, -3.0e+02])
        L_u, L_v, L_w, sigma_u, sigma_v, sigma_w = turb_scale_inten(X0)
        b, _ = geometry_and_trim()  # Airspeed V is already defined

        omega = np.logspace(-5, 5, 10000)
        Phi_u_g, Phi_v_g, Phi_w_g, Phi_p_g, Phi_q_g, Phi_r_g = dryden_psd(omega, V, b, L_u, L_v, L_w, sigma_u, sigma_v, sigma_w, plot_psd=0)

        stan_dev_u = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_u_g, omega))
        stan_dev_v = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_v_g, omega))
        stan_dev_w = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_w_g, omega))
        stan_dev_p = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_p_g, omega))
        stan_dev_q = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_q_g, omega))
        stan_dev_r = np.sqrt(2 * (1 / (2 * np.pi)) * np.trapz(Phi_r_g, omega))

        std_u.append(stan_dev_u)
        std_v.append(stan_dev_v)
        std_w.append(stan_dev_w)
        std_p.append(stan_dev_p)
        std_q.append(stan_dev_q)
        std_r.append(stan_dev_r)

    return std_u, std_v, std_w, std_p, std_q, std_r

# Airspeed range to test
airspeeds = np.linspace(1, 100, 50)  # Example airspeeds in m/s

# Calculate standard deviations for each airspeed
std_u, std_v, std_w, std_p, std_q, std_r = calculate_std_vs_airspeed(airspeeds)

# User input
whether_to_plot_velocity_dependence = 0

if whether_to_plot_velocity_dependence == 1:

    # Plotting
    plt.figure()
    plt.plot(airspeeds, std_u, label='std_u')
    plt.plot(airspeeds, std_v, label='std_v')
    # plt.plot(airspeeds, std_w, label='std_w')
    # plt.plot(airspeeds, std_p, label='std_p')
    # plt.plot(airspeeds, std_q, label='std_q')
    # plt.plot(airspeeds, std_r, label='std_r')
    plt.xlabel('Airspeed [m/s]')
    plt.ylabel('Standard Deviation')
    plt.title(r'Standard Deviations of $u_g$, $v_g$ vs Airspeed')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plotting
    plt.figure()
    # plt.plot(airspeeds, std_u, label='std_u')
    # plt.plot(airspeeds, std_v, label='std_v')
    plt.plot(airspeeds, std_w, label='std_w')
    # plt.plot(airspeeds, std_p, label='std_p')
    # plt.plot(airspeeds, std_q, label='std_q')
    # plt.plot(airspeeds, std_r, label='std_r')
    plt.xlabel('Airspeed [m/s]')
    plt.ylabel('Standard Deviation')
    plt.title(r'Standard Deviations of $w_g$ vs Airspeed')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plotting
    plt.figure()
    # plt.plot(airspeeds, std_u, label='std_u')
    # plt.plot(airspeeds, std_v, label='std_v')
    # plt.plot(airspeeds, std_w, label='std_w')
    plt.plot(airspeeds, std_p, label='std_p')
    plt.plot(airspeeds, std_q, label='std_q')
    plt.plot(airspeeds, std_r, label='std_r')
    plt.xlabel('Airspeed [m/s]')
    plt.ylabel('Standard Deviation')
    plt.title(r'Standard Deviations of $p_g$, $q_g$, $r_g$ vs Airspeed')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Airspeed dependence is very small, negligible.