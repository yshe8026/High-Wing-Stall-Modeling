import numpy as np
import matplotlib.pyplot as plt

# Constants
cd0 = 0.008 # zero-lift drag coefficient of the 2D aerofoil section used in the propeller blade
k1 = 0.003  # coefficient used in the drag polar equation to account for the linear relationship between the lift coefficient (Cl) and the drag coefficient (Cd)
k2 = 0.01   # another coefficient used in the drag polar equation, but it accounts for the quadratic relationship between the lift coefficient (Cl) and the drag coefficient (Cd)
chord = 0.10 # chord length of the propeller blade, which is the distance from the leading edge to the trailing edge of the blade cross-section
pitch = 1.346 # pitch of the propeller, which represents the distance the propeller would move forward in one complete rotation in an ideal scenario with no slip # J400 Manual 1.346 m
dia = 1.524 # diameter of the propeller # default: 1.6 m, Jabiru J400: 1.52 m # J400 Manual 1.524 m
R = dia / 2.0 # radius of the propeller
RPM = 2100.0 # rotational speed of the engine driving the propeller, measured in revolutions per minute (RPM)
tonc = 0.12 * chord # thickness-to-chord ratio of the propeller blade section, which represents the relative thickness of the blade cross-section compared to its chord length
rho = 1.225 # air density at standard sea-level conditions (kg/m³)
n = RPM / 60.0 # rotational speed of the propeller in revolutions per second (RPS), converted from the given RPM
omega = n * 2.0 * np.pi # angular velocity of the propeller in radians per second, derived from the RPS
B = 2 # number of blades on the propeller
xs = 0.1 * R # starting radial position of the blade elements used in the calculation, set at 10% of the propeller's radius from the hub (defines the inner boundary of the blade element analysis, ignoring the hub section where the blade shape may be more complex)
xt = R # tip radius of the blade elements used in the calculation, corresponding to the outer boundary of the propeller
rstep = (xt - xs) / 10 # step size for the radial positions of the blade elements (blade is divided into 10 segments for analysis, and this step size determines the spacing between the analyzed points along the radius)

xt_plus_rstep = xt + rstep - 0.00001 # spinning pigeon head: make it ever slightly smaller than "xt + rstep" to fix a issue with np.arange()
r1 = np.arange(xs, xt_plus_rstep, rstep)

# Initialize arrays to store results
t = np.zeros(60)
q = np.zeros(60)
J = np.zeros(60)
eff = np.zeros(60)
icheck = np.zeros(60)

# Loop over a range of velocities from 1 to 60 m/s
for V in range(1, 61):
    thrust = 0.0
    torque = 0.0

    for rad in r1:
        theta = np.arctan(pitch / (2.0 * np.pi * rad))
        sigma = 2.0 * chord / (2.0 * np.pi * rad)
        a = 0.1
        b = 0.01
        finished = False
        sum_iter = 1
        itercheck = 0

        while not finished:
            V0 = V * (1 + a)
            V2 = omega * rad * (1 - b)
            phi = np.arctan2(V0, V2)
            alpha = theta - phi
            cl = 6.2 * alpha
            cd = cd0 - k1 * cl + k2 * cl**2
            Vlocal = np.sqrt(V0**2 + V2**2)
            DtDr = 0.5 * rho * Vlocal**2 * B * chord * (cl * np.cos(phi) - cd * np.sin(phi))
            DqDr = 0.5 * rho * Vlocal**2 * B * chord * rad * (cd * np.cos(phi) + cl * np.sin(phi))
            tem1 = DtDr / (4.0 * np.pi * rad * rho * V**2 * (1 + a))
            tem2 = DqDr / (4.0 * np.pi * rad**3 * rho * V * (1 + a) * omega)
            anew = 0.5 * (a + tem1)
            bnew = 0.5 * (b + tem2)

            if abs(anew - a) < 1.0e-5 and abs(bnew - b) < 1.0e-5:
                finished = True

            a = anew
            b = bnew
            sum_iter += 1

            if sum_iter > 500:
                finished = True
                itercheck = 1

        thrust += DtDr * rstep
        torque += DqDr * rstep

    t[V-1] = thrust / (rho * n**2 * dia**4)
    q[V-1] = torque / (rho * n**2 * dia**5)
    J[V-1] = V / (n * dia)
    eff[V-1] = J[V-1] / 2.0 / np.pi * t[V-1] / q[V-1]
    icheck[V-1] = itercheck

Jmax = np.max(J)
Tmax = np.max(t)

# print(r1)
# print(rstep)
# print(thrust)

# # Plot results
# fig, ax1 = plt.subplots()

# ax1.set_xlabel('Advance Ratio (J)')
# ax1.set_ylabel('Ct', color='tab:blue')
# ax1.plot(J, t, color='tab:blue', label='Ct')
# ax1.tick_params(axis='y', labelcolor='tab:blue')

# ax2 = ax1.twinx()
# ax2.set_ylabel('Cq', color='tab:red')
# ax2.plot(J, q, color='tab:red', label='Cq')
# ax2.tick_params(axis='y', labelcolor='tab:red')

# plt.title('Thrust and Torque Coefficients')
# fig.tight_layout()
# plt.legend(['Ct', 'Cq'], loc='upper left')
# plt.show()

# plt.figure()
# plt.plot(J, eff)
# plt.title('Propeller Efficiency')
# plt.xlabel('Advance Ratio (J)')
# plt.ylabel('Efficiency')
# plt.axis([0, Jmax, 0, 1])
# plt.show()

# Plot results
fig, ax1 = plt.subplots()

ax1.set_xlabel('Advance Ratio ($J$)')
ax1.set_ylabel('$C_t$', color='tab:blue')
line1, = ax1.plot(J, t, color='tab:blue', label='$C_t$')
ax1.tick_params(axis='y', labelcolor='tab:blue')

ax2 = ax1.twinx()
ax2.set_ylabel('$C_q$', color='tab:red')
line2, = ax2.plot(J, q, color='tab:red', label='$C_q$')
ax2.tick_params(axis='y', labelcolor='tab:red')

plt.title('Thrust and Torque Coefficients')

# Create a single legend for both lines from different axes
lines = [line1, line2]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc='upper right')

fig.tight_layout()
plt.show()

plt.figure()
plt.plot(J, eff)
plt.title('Propeller Efficiency')
plt.xlabel('Advance Ratio ($J$)')
plt.ylabel('Efficiency')
plt.axis([0, Jmax, 0, 1])
plt.show()


# print(eff)