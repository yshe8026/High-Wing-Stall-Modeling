import matplotlib.pyplot as plt
import numpy as np

# Example aerodynamic coefficient lookup tables for NACA 4412 (wing) and NACA 0012 (tail)
# Note: The keys are angles of attack in degrees and the values are aerodynamic coefficients
# You should replace the example values with actual aerodynamic data

wing_aero_data = {
    # NACA 4412 wing aerodynamic data
    "CL": {
        -10: -0.5, 0: 0.0, 10: 1.0, 15: 1.3, 16: 1.4, 17: 1.2, 20: 0.9,
    },
    "CD": {
        -10: 0.02, 0: 0.006, 10: 0.02, 15: 0.03, 16: 0.031, 17: 0.033, 20: 0.045,
    },
    "CM": {
        -10: -0.1, 0: 0.0, 10: -0.02, 15: -0.05, 16: -0.055, 17: -0.06, 20: -0.1,
    }
}

tail_aero_data = {
    # NACA 0012 tail aerodynamic data
    "CL": {
        -10: -0.6, 0: 0.0, 10: 0.6, 15: 0.8, 20: 0.7,
    },
    "CD": {
        -10: 0.025, 0: 0.007, 10: 0.025, 15: 0.04, 20: 0.05,
    },
    "CM": {
        -10: -0.02, 0: 0.0, 10: -0.01, 15: -0.03, 20: -0.05,
    }
}

# Function to perform linear interpolation
def linear_interpolate(x, x0, y0, x1, y1):
    return y0 + (y1 - y0) * ((x - x0) / (x1 - x0))

# Updated function to get aerodynamic coefficients with interpolation
def get_aero_coefficients(aero_data, alpha):
    """
    Retrieve the aerodynamic coefficients for a given angle of attack (alpha) from the lookup table.
    If the specific alpha is not in the table, performs linear interpolation between the two closest points.
    """
    # Directly return the value if the exact match is found
    if alpha in aero_data["CL"]:
        cl = aero_data["CL"][alpha]
        cd = aero_data["CD"][alpha]
        cm = aero_data["CM"][alpha]
    else:
        # Find the closest angles larger and smaller than alpha, for interpolation
        all_angles = sorted(aero_data["CL"].keys())
        larger_angles = [angle for angle in all_angles if angle > alpha]
        smaller_angles = [angle for angle in all_angles if angle < alpha]

        if larger_angles and smaller_angles:
            lower_alpha = max(smaller_angles)
            upper_alpha = min(larger_angles)

            # Perform linear interpolation for each coefficient
            cl = linear_interpolate(alpha, lower_alpha, aero_data["CL"][lower_alpha], upper_alpha, aero_data["CL"][upper_alpha])
            cd = linear_interpolate(alpha, lower_alpha, aero_data["CD"][lower_alpha], upper_alpha, aero_data["CD"][upper_alpha])
            cm = linear_interpolate(alpha, lower_alpha, aero_data["CM"][lower_alpha], upper_alpha, aero_data["CM"][upper_alpha])
        else:
            # If alpha is outside the range of known angles, return None or raise an error
            return None, None, None # Or handle this scenario as needed

    return cl, cd, cm
# Example usage
alpha_wing = 10.1  # Example angle of attack for the wing
alpha_tail = -5.1  # Example angle of attack for the tail, assuming you'd interpolate or adjust this example

# Retrieve coefficients for the wing and tail at specified angles of attack
wing_cl, wing_cd, wing_cm = get_aero_coefficients(wing_aero_data, alpha_wing)
tail_cl, tail_cd, tail_cm = get_aero_coefficients(tail_aero_data, alpha_tail)

print(f"Wing coefficients at {alpha_wing} degrees: CL={wing_cl}, CD={wing_cd}, CM={wing_cm}")
print(f"Tail coefficients at {alpha_tail} degrees: CL={tail_cl}, CD={tail_cd}, CM={tail_cm}")





# Assuming wing_aero_data is defined as before and contains CL data for NACA 4412

# Generate interpolated CL values over a range of angles
def generate_interpolated_values(aero_data, angle_range):
    interpolated_values = []
    for alpha in angle_range:
        cl, _, _ = get_aero_coefficients(aero_data, alpha)  # We're focusing on CL here
        interpolated_values.append(cl)
    return interpolated_values

# Define the range of angles of attack to consider
angle_range = np.linspace(-10, 20, 300)  # From -10 to 20 degrees, 300 points

# Generate interpolated CL values
interpolated_cl_values = generate_interpolated_values(wing_aero_data, angle_range)

# Plotting
plt.figure(figsize=(10, 6))

# Plot actual data points
actual_angles = sorted(wing_aero_data['CL'].keys())
actual_cl = [wing_aero_data['CL'][angle] for angle in actual_angles]
plt.scatter(actual_angles, actual_cl, color='red', label='Actual Data Points', zorder=5)

# Plot interpolated CL curve
plt.plot(angle_range, interpolated_cl_values, label='Interpolated CL Curve', linestyle='--', color='blue')

plt.title('Coefficient of Lift (CL) vs. Angle of Attack for NACA 4412')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Coefficient of Lift (CL)')
plt.legend()
plt.grid(True)
plt.show()