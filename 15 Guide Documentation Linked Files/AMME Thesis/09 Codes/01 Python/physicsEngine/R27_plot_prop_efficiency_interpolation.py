import numpy as np
import matplotlib.pyplot as plt

# propeller diameter is 1.6 m

J = [
    0.01785714, 0.03571429, 0.05357143, 0.07142857, 0.08928571, 0.10714286,
    0.125, 0.14285714, 0.16071429, 0.17857143, 0.19642857, 0.21428571,
    0.23214286, 0.25, 0.26785714, 0.28571429, 0.30357143, 0.32142857,
    0.33928571, 0.35714286, 0.375, 0.39285714, 0.41071429, 0.42857143,
    0.44642857, 0.46428571, 0.48214286, 0.5, 0.51785714, 0.53571429,
    0.55357143, 0.57142857, 0.58928571, 0.60714286, 0.625, 0.64285714,
    0.66071429, 0.67857143, 0.69642857, 0.71428571, 0.73214286, 0.75,
    0.76785714, 0.78571429, 0.80357143, 0.82142857, 0.83928571, 0.85714286,
    0.875, 0.89285714, 0.91071429, 0.92857143, 0.94642857, 0.96428571,
    0.98214286, 1.0, 1.01785714, 1.03571429, 1.05357143, 1.07142857
]



eff = [
    0.05915038, 0.11558659, 0.16936015, 0.22052606, 0.26914243, 0.31527017,
    0.35897298, 0.40031568, 0.43936447, 0.47618648, 0.51084819, 0.54341613,
    0.57395401, 0.60252332, 0.62918193, 0.65398157, 0.6769681, 0.69817755,
    0.71763545, 0.73535201, 0.75131859, 0.76549979, 0.77782654, 0.78817943,
    0.7963704, 0.80210603, 0.80493168, 0.80413188, 0.79854485, 0.78617448,
    0.76336591, 0.7226753, 0.64633477, 0.47922336, -0.07676112, 11.58912973,
    1.91083444, 1.49040356, 1.34621077, 1.27548983, 1.23476453, 1.20910426,
    1.19203273, 1.18027759, 1.17201194, 1.16613953, 1.16197046, 1.15902912,
    1.15700645, 1.15567842, 1.15488305, 1.15450078, 1.15444187, 1.1546383,
    1.155037, 1.15559769, 1.15628827, 1.15708889, 1.15796668, 1.15891201
]

# propeller diameter is 1.52 m

J = [
    0.01879699, 0.03759398, 0.05639098, 0.07518797, 0.09398496, 0.11278195,
    0.13157895, 0.15037594, 0.16917293, 0.18796992, 0.20676692, 0.22556391,
    0.2443609, 0.26315789, 0.28195489, 0.30075188, 0.31954887, 0.33834586,
    0.35714286, 0.37593985, 0.39473684, 0.41353383, 0.43233083, 0.45112782,
    0.46992481, 0.4887218, 0.5075188, 0.52631579, 0.54511278, 0.56390977,
    0.58270677, 0.60150376, 0.62030075, 0.63909774, 0.65789474, 0.67669173,
    0.69548872, 0.71428571, 0.73308271, 0.7518797, 0.77067669, 0.78947368,
    0.80827068, 0.82706767, 0.84586466, 0.86466165, 0.88345865, 0.90225564,
    0.92105263, 0.93984962, 0.95864662, 0.97744361, 0.9962406, 1.01503759,
    1.03383459, 1.05263158, 1.07142857, 1.09022556, 1.10902256, 1.12781955
]


eff = [
    0.05955886, 0.11637408, 0.17049907, 0.22199082, 0.27090961, 0.31731859,
    0.36128352, 0.40287216, 0.44215286, 0.47919485, 0.51406841, 0.54684146,
    0.57758148, 0.60635237, 0.63321508, 0.6582252, 0.68143103, 0.70287318,
    0.72258152, 0.74057134, 0.75683925, 0.77135824, 0.78406627, 0.79485631,
    0.80355403, 0.80988511, 0.81342125, 0.81348382, 0.80895862, 0.7979328,
    0.7768597, 0.73844679, 0.66500263, 0.50009088, -0.0854262, 6.27273046,
    1.78626917, 1.43922115, 1.31528891, 1.25373688, 1.21819651, 1.1958613,
    1.18110697, 1.17106066, 1.16411316, 1.15929417, 1.15598434, 1.15376533,
    1.1523596, 1.15156717, 1.15124332, 1.15128125, 1.15160089, 1.15214147,
    1.15285602, 1.15370873, 1.15467214, 1.15572593, 1.15684349, 1.15801581
]


# Improving the realism of the BEM efficiency result (propeller diameter is 1.6 m)

J = [0,
    0.01785714, 0.03571429, 0.05357143, 0.07142857, 0.08928571, 0.10714286,
    0.125, 0.14285714, 0.16071429, 0.17857143, 0.19642857, 0.21428571,
    0.23214286, 0.25, 0.26785714, 0.28571429, 0.30357143, 0.32142857,
    0.33928571, 0.35714286, 0.375, 0.39285714, 0.41071429, 0.42857143,
    0.44642857, 0.46428571, 0.48214286, 0.5, 0.51785714, 0.53571429,
    0.55357143, 0.57142857, 0.58928571, 0.60714286, 0.625
]



eff = [0,
    0.05915038, 0.11558659, 0.16936015, 0.22052606, 0.26914243, 0.31527017,
    0.35897298, 0.40031568, 0.43936447, 0.47618648, 0.51084819, 0.54341613,
    0.57395401, 0.60252332, 0.62918193, 0.65398157, 0.6769681, 0.69817755,
    0.71763545, 0.73535201, 0.75131859, 0.76549979, 0.77782654, 0.78817943,
    0.7963704, 0.80210603, 0.80493168, 0.80413188, 0.79854485, 0.78617448,
    0.76336591, 0.7226753, 0.64633477, 0.47922336, 0.000
]


# Improving the realism of the BEM efficiency result (propeller diameter is 1.52 m)

J = [0.000,
    0.01879699, 0.03759398, 0.05639098, 0.07518797, 0.09398496, 0.11278195,
    0.13157895, 0.15037594, 0.16917293, 0.18796992, 0.20676692, 0.22556391,
    0.2443609, 0.26315789, 0.28195489, 0.30075188, 0.31954887, 0.33834586,
    0.35714286, 0.37593985, 0.39473684, 0.41353383, 0.43233083, 0.45112782,
    0.46992481, 0.4887218, 0.5075188, 0.52631579, 0.54511278, 0.56390977,
    0.58270677, 0.60150376, 0.62030075, 0.63909774, 0.65789474
]


eff = [0.000,
    0.05955886, 0.11637408, 0.17049907, 0.22199082, 0.27090961, 0.31731859,
    0.36128352, 0.40287216, 0.44215286, 0.47919485, 0.51406841, 0.54684146,
    0.57758148, 0.60635237, 0.63321508, 0.6582252, 0.68143103, 0.70287318,
    0.72258152, 0.74057134, 0.75683925, 0.77135824, 0.78406627, 0.79485631,
    0.80355403, 0.80988511, 0.81342125, 0.81348382, 0.80895862, 0.7979328,
    0.7768597, 0.73844679, 0.66500263, 0.50009088, 0.000
]

Jmax = max(J)

# Fit a 15th order polynomial
polynomial_coefficients = np.polyfit(J, eff, 15)
polynomial_expression = np.poly1d(polynomial_coefficients)

# Set the print precision to 16 significant digits
np.set_printoptions(precision=16)
print(polynomial_coefficients)

polynomial_coefficients = [-2.67438907e+08,  1.24592095e+09, -2.61422244e+09,  3.26408269e+09,
                           -2.69939288e+09,  1.55737948e+09, -6.43083806e+08,  1.91878929e+08,
                           -4.12570203e+07,  6.30205417e+06, -6.65243868e+05,  4.63408791e+04,
                           -1.97163351e+03,  4.07222763e+01,  2.85202480e+00,  2.45521078e-05]
# improving significant figures
polynomial_coefficients = [-2.6743890653121743e+08,  1.2459209493015354e+09, -2.6142224426932373e+09,
                            3.2640826865733147e+09, -2.6993928791204376e+09,  1.5573794758376265e+09,
                           -6.4308380623244667e+08,  1.9187892870937729e+08, -4.1257020279518954e+07,
                            6.3020541739022173e+06, -6.6524386816298915e+05,  4.6340879054777324e+04,
                           -1.9716335106041690e+03,  4.0722276272755224e+01,  2.8520248036777662e+00,
                            2.4552107788622379e-05]
polynomial_expression = np.poly1d(polynomial_coefficients)



# print(polynomial_expression)

# Generate values for plotting the 15th order polynomial fit
J_fit = np.linspace(min(J), Jmax, 500)
eff_fit = polynomial_expression(J_fit)

eff_fit = (-2.6743890653121743e+08 * J_fit**15 +
       1.2459209493015354e+09 * J_fit**14 -
       2.6142224426932373e+09 * J_fit**13 +
       3.2640826865733147e+09 * J_fit**12 -
       2.6993928791204376e+09 * J_fit**11 +
       1.5573794758376265e+09 * J_fit**10 -
       6.4308380623244667e+08 * J_fit**9 +
       1.9187892870937729e+08 * J_fit**8 -
       4.1257020279518954e+07 * J_fit**7 +
       6.3020541739022173e+06 * J_fit**6 -
       6.6524386816298915e+05 * J_fit**5 +
       4.6340879054777324e+04 * J_fit**4 -
       1.9716335106041690e+03 * J_fit**3 +
       4.0722276272755224e+01 * J_fit**2 +
       2.8520248036777662 * J_fit +
       2.4552107788622379e-05)


# eff = -2.674389e+08*J**15 + 1.245921e+09*J**14 - 2.614223e+09*J**13 + 3.264083e+09*J**12 - 2.699393e+09*J**11 
#        + 1.557380e+09*J**10 - 6.430839e+08*J**9 + 1.918789e+08*J**8 - 4.125702e+07*J**7 + 6.302055e+06*J**6 
#        - 6.652439e+05*J**5 + 4.634088e+04*J**4 - 1.971634e+03*J**3 + 4.072228e+01*J**2 + 2.852025e+00*J + 2.455536e-05

# eff_fit_15th = (
#     -2.674389e+08*J**15 + 1.245921e+09*J**14 - 2.614223e+09*J**13 + 3.264083e+09*J**12 - 2.699393e+09*J**11 +
#     1.557380e+09*J**10 - 6.430839e+08*J**9 + 1.918789e+08*J**8 - 4.125702e+07*J**7 + 6.302055e+06*J**6 -
#     6.652439e+05*J**5 + 4.634088e+04*J**4 - 1.971634e+03*J**3 + 4.072228e+01*J**2 + 2.852025e+00*J + 2.455536e-05
# )

# eff_fit_15th = (
#     -2.67438907e+08*J**15 + 1.24592095e+09*J**14 - 2.61422244e+09*J**13 + 3.26408269e+09*J**12 - 2.69939288e+09*J**11 +
#     1.55737948e+09*J**10 - 6.43083806e+08*J**9 + 1.91878929e+08*J**8 - 4.12570203e+07*J**7 + 6.30205417e+06*J**6 -
#     6.65243868e+05*J**5 + 4.63408791e+04*J**4 - 1.97163351e+03*J**3 + 4.07222763e+01*J**2 + 2.85202480e+00*J + 2.45521078e-05
# )

# eff_fit_15th = (
#     -2.67438907e+08*J_fit**15 + 1.24592095e+09*J_fit**14 - 2.61422244e+09*J_fit**13 + 3.26408269e+09*J_fit**12 - 2.69939288e+09*J_fit**11 +
#     1.55737948e+09*J_fit**10 - 6.43083806e+08*J_fit**9 + 1.91878929e+08*J_fit**8 - 4.12570203e+07*J_fit**7 + 6.30205417e+06*J_fit**6 -
#     6.65243868e+05*J_fit**5 + 4.63408791e+04*J_fit**4 - 1.97163351e+03*J_fit**3 + 4.07222763e+01*J_fit**2 + 2.85202480e+00*J_fit + 2.45521078e-05
# )

# Plot the original data
plt.figure(figsize=(10, 6))
plt.plot(J, eff, 'o', label='Original Data')

# Plot the 15th order polynomial fit
plt.plot(J_fit, eff_fit, '-', label='15th Order Polynomial Fit')

plt.title('Propeller Efficiency')
plt.xlabel('Advance Ratio ($J$)')
plt.ylabel('Efficiency')
plt.axis([0, Jmax, 0, 1])
plt.legend()
plt.grid(True)
plt.show()

# plt.figure()
# plt.plot(J, eff)
# plt.title('Propeller Efficiency')
# plt.xlabel('Advance Ratio ($J$)')
# plt.ylabel('Efficiency')
# plt.axis([0, Jmax, 0, 1])
# plt.show()

# # Plot results
# fig, ax1 = plt.subplots()

# ax1.set_xlabel('Advance Ratio ($J$)')
# ax1.set_ylabel('$C_t$', color='tab:blue')
# line1, = ax1.plot(J, t, color='tab:blue', label='$C_t$')
# ax1.tick_params(axis='y', labelcolor='tab:blue')

# ax2 = ax1.twinx()
# ax2.set_ylabel('$C_q$', color='tab:red')
# line2, = ax2.plot(J, q, color='tab:red', label='$C_q$')
# ax2.tick_params(axis='y', labelcolor='tab:red')

# plt.title('Thrust and Torque Coefficients')

# # Create a single legend for both lines from different axes
# lines = [line1, line2]
# labels = [line.get_label() for line in lines]
# ax1.legend(lines, labels, loc='upper right')

# fig.tight_layout()
# plt.show()

