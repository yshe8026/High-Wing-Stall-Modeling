import numpy as np
import matplotlib.pyplot as plt
#============================================================================================================

# Piece 1 + Piece 2 + Piece 3: propeller diameter is 1.52 m, pitch is 1.346 m

J = [0.000,
    0.01874766, 0.03749531, 0.05624297, 0.07499063, 0.09373828, 0.11248594,
    0.1312336,  0.14998125, 0.16872891, 0.18747657, 0.20622422, 0.22497188,
    0.24371954, 0.26246719, 0.28121485, 0.2999625,  0.31871016, 0.33745782,
    0.35620547, 0.37495313, 0.39370079, 0.41244844, 0.4311961,  0.44994376,
    0.46869141, 0.48743907, 0.50618673, 0.52493438, 0.54368204, 0.5624297,
    0.58117735, 0.59992501, 0.61867267, 0.63742032, 0.65616798, 0.67491564,
    0.69366329, 0.71241095, 0.73115861, 0.74990626, 0.76865392, 0.78740157,
    0.80614923, 0.82489689, 0.84364454, 0.8623922,  0.88113986]

eff = [0.000,
        0.05066336, 0.09929375, 0.14592944, 0.19061103, 0.23338125, 0.27428488,
        0.31336851, 0.35068013, 0.38626937, 0.42018657, 0.45248358, 0.48321171,
        0.51242301, 0.54017017, 0.56650397, 0.59147658, 0.61513809, 0.63753735,
        0.65872229, 0.67873928, 0.69763355, 0.71544374, 0.732212,   0.74797312,
        0.76276091, 0.776604,   0.78952665, 0.80154821, 0.81268147, 0.82293242,
        0.83229563, 0.84075536, 0.84828156, 0.85482226, 0.86029943, 0.86459484,
        0.86753653, 0.86886971, 0.86821151, 0.86495917, 0.85814345, 0.84611525,
        0.82579236, 0.79077725, 0.72533662, 0.57750443, 0.01443093]

Jmax = max(J)

# Piece 1
J_p1 = [0.000,
            0.01874766, 0.03749531, 0.05624297, 0.07499063, 0.09373828, 0.11248594,
            0.1312336,  0.14998125, 0.16872891, 0.18747657, 0.20622422, 0.22497188,
            0.24371954, 0.26246719, 0.28121485, 0.2999625,  0.31871016, 0.33745782,
            0.35620547, 0.37495313, 0.39370079, 0.41244844, 0.4311961,  0.44994376,
            0.46869141, 0.48743907, 0.50618673, 0.52493438, 0.54368204]

eff_p1 = [0.000,
                0.05066336, 0.09929375, 0.14592944, 0.19061103, 0.23338125, 0.27428488,
                0.31336851, 0.35068013, 0.38626937, 0.42018657, 0.45248358, 0.48321171,
                0.51242301, 0.54017017, 0.56650397, 0.59147658, 0.61513809, 0.63753735,
                0.65872229, 0.67873928, 0.69763355, 0.71544374, 0.732212,   0.74797312,
                0.76276091, 0.776604,   0.78952665, 0.80154821, 0.81268147]

Jmax_p1 = max(J_p1)

# Piece 2
J_p2 = [0.54368204, 0.5624297,
        0.58117735, 0.59992501, 0.61867267, 0.63742032, 0.65616798, 0.67491564,
        0.69366329, 0.71241095, 0.73115861, 0.74990626, 0.76865392, 0.78740157,
        0.80614923, 0.82489689]

eff_p2 = [0.81268147, 0.82293242,
        0.83229563, 0.84075536, 0.84828156, 0.85482226, 0.86029943, 0.86459484,
        0.86753653, 0.86886971, 0.86821151, 0.86495917, 0.85814345, 0.84611525,
        0.82579236, 0.79077725]

Jmax_p2 = max(J_p2)

# Piece 3
J_p3 = [0.80614923, 0.82489689, 0.84364454, 0.8623922,  0.88113986]

eff_p3 = [0.82579236, 0.79077725, 0.72533662, 0.57750443, 0.01443093]

Jmax_p3 = max(J_p3)


#=============================================================================================================

# Piece 1 + Piece 2 + Piece 3
# Fit a 5th order polynomial
polynomial_coefficients = np.polyfit(J, eff, 5)
polynomial_expression = np.poly1d(polynomial_coefficients)

# # Set the print precision to 16 significant digits
np.set_printoptions(precision=16)
print(polynomial_coefficients)


polynomial_coefficients = [-6.2488040647582864e+01,  1.2029855562514737e+02, -8.0155951548989478e+01,
                            1.9742625027929112e+01,  4.4711017865478048e-01,  4.6476147646879323e-02]

# print(polynomial_expression)


# Piece 1
# Fit a 5th order polynomial
polynomial_coefficients_p1 = np.polyfit(J_p1, eff_p1, 5)
polynomial_expression_p1 = np.poly1d(polynomial_coefficients_p1)

# # Set the print precision to 16 significant digits
np.set_printoptions(precision=16)
print(polynomial_coefficients_p1)


polynomial_coefficients_p1 = [-1.0339808355672748e+00,  1.1513141630282613e+00,  8.0980324477669519e-01,
                              -2.9358586169893126e+00,  2.7568991029961047e+00,  3.0600838902929569e-06]


# Piece 2
# Fit a 5th order polynomial
polynomial_coefficients_p2 = np.polyfit(J_p2, eff_p2, 5)
polynomial_expression_p2 = np.poly1d(polynomial_coefficients_p2)

# # Set the print precision to 16 significant digits
np.set_printoptions(precision=16)
print(polynomial_coefficients_p2)


polynomial_coefficients_p2 = [-414.56709062760984,  1341.4163943124554,  -1733.7299437622642,
                                1117.049389105175,    -357.9291035961872,     46.33549936554715]


# Piece 3
# Fit a 4th order polynomial
polynomial_coefficients_p3 = np.polyfit(J_p3, eff_p3, 4)
polynomial_expression_p3 = np.poly1d(polynomial_coefficients_p3)

# # Set the print precision to 16 significant digits
np.set_printoptions(precision=16)
print(polynomial_coefficients_p3)


polynomial_coefficients_p3 = [-94738.66870330359, 314836.3898499948, -392340.3248280777,
                              217290.48582230246, -45125.05373275733]

#=============================================================================================================

# Piece 1 + Piece 2 + Piece 3
# Generate values for plotting the 15th order polynomial fit
# J_fit = np.linspace(min(J), Jmax, 500)
J_fit = np.linspace(-0.5, Jmax+0.5, 500)
eff_fit = polynomial_expression(J_fit)


eff_fit = (-6.2488040647582864e+01 * J_fit**5 +
       1.2029855562514737e+02 * J_fit**4 +
       -8.0155951548989478e+01 * J_fit**3 +
       1.9742625027929112e+01 * J_fit**2 +
       4.4711017865478048e-01 * J_fit +
       4.6476147646879323e-02)

# Piece 1
# Generate values for plotting the 15th order polynomial fit
# J_fit = np.linspace(min(J), Jmax, 500)
J_fit_p1 = np.linspace(-0.5, Jmax_p1+0.5, 500)
eff_fit_p1 = polynomial_expression(J_fit_p1)


eff_fit_p1 = (-1.0339808355672748e+00 * J_fit_p1**5 +
       1.1513141630282613e+00 * J_fit_p1**4 +
       8.0980324477669519e-01 * J_fit_p1**3 +
       -2.9358586169893126e+00 * J_fit_p1**2 +
       2.7568991029961047e+00 * J_fit_p1 +
       3.0600838902929569e-06)

# Piece 2
# Generate values for plotting the 15th order polynomial fit
# J_fit = np.linspace(min(J), Jmax, 500)
J_fit_p2 = np.linspace(-0.5, Jmax_p2+0.5, 500)
eff_fit_p2 = polynomial_expression(J_fit_p2)


eff_fit_p2 = (-414.56709062760984 * J_fit_p2**5 +
       1341.4163943124554 * J_fit_p2**4 +
       -1733.7299437622642 * J_fit_p2**3 +
       1117.049389105175 * J_fit_p2**2 +
       -357.9291035961872 * J_fit_p2 +
       46.33549936554715)

# Piece 3
# Generate values for plotting the 15th order polynomial fit
# J_fit = np.linspace(min(J), Jmax, 500)
J_fit_p3 = np.linspace(-0.5, Jmax_p3+0.5, 500)
eff_fit_p3 = polynomial_expression(J_fit_p3)


eff_fit_p3 = (-94738.66870330359 * J_fit_p3**4 +
              314836.3898499948 * J_fit_p3**3 +
               -392340.3248280777 * J_fit_p3**2 +
                217290.48582230246 * J_fit_p3 +
                -45125.05373275733)

#=============================================================================================================
# Plot the original data
plt.figure(figsize=(10, 6))
plt.plot(J, eff, 'o', label='Original Data')

# Piece 1 + Piece 2 + Piece 3
# Plot the 5th order polynomial fit
plt.plot(J_fit, eff_fit, '-', label='5th Order Polynomial Fit')

# Piece 1
# Plot the 5th order polynomial fit
plt.plot(J_fit_p1, eff_fit_p1, '-', label='Piece1: 5th Order Polynomial Fit')

# Piece 2
# Plot the 5th order polynomial fit
plt.plot(J_fit_p2, eff_fit_p2, '-', label='Piece2: 5th Order Polynomial Fit')

# Piece 3
# Plot the 3rd order polynomial fit
plt.plot(J_fit_p3, eff_fit_p3, '-', label='Piece3: 4th Order Polynomial Fit')


plt.title('Propeller Efficiency')
plt.xlabel('Advance Ratio ($J$)')
plt.ylabel('Efficiency')
plt.axis([-0.5, Jmax+0.5, -1, 1])
plt.legend()
plt.grid(True)
plt.show()

