import matplotlib.pyplot as plt
import numpy as np

def read_xplane_airfoil_data(filepath):
    """
    Reads an X-Plane 11 airfoil file, extracting aerodynamic coefficients.
    Assumes the coefficient section starts with "alpha cl cd cm:".
    """
    # Flags to track whether we're in the aerodynamic data section
    in_aero_section = False
    
    # Storage for the aerodynamic coefficients
    aero_data = {
        "alpha": [],
        "cl": [],
        "cd": [],
        "cm": []
    }

    with open(filepath, 'r') as file:
        for line in file:
            # Trim whitespace and ignore comments
            line = line.strip()
            if line.startswith('#') or not line:
                continue

            # Check if we're entering the aerodynamic data section
            if line.lower() == "alpha cl cd cm:":
                in_aero_section = True
                continue

            # Process data if we're in the correct section
            if in_aero_section:
                # Attempt to parse the aerodynamic data
                try:
                    parts = line.split()
                    # Check if we're transitioning out of the aerodynamic section
                    if len(parts) < 4:
                        break  # Exit loop if we don't have 4 parts, indicating end of aero data
                    alpha, cl, cd, cm = map(float, parts[:4])
                    aero_data["alpha"].append(alpha)
                    aero_data["cl"].append(cl)
                    aero_data["cd"].append(cd)
                    aero_data["cm"].append(cm)
                except ValueError:
                    # Handle the case where conversion to float fails
                    print(f"Warning: Could not process line '{line}'. Skipping.")
                    continue

    return aero_data

def check_lengths(aero_data):
    """
    Check if the lengths of alpha, cl, cd, and cm lists in the aerodynamic data are the same.
    """
    lengths = [len(aero_data[key]) for key in ['alpha', 'cl', 'cd', 'cm']]
    if len(set(lengths)) == 1:
        print("All lists are of the same length.")
    else:
        print("Warning: Lists are of differing lengths.")

# Use a relative path since the file is in the same directory as the script
filename = "physicsEngine\\NACA_4412.txt"
aero_data = read_xplane_airfoil_data(filename)
# 00 print(aero_data)

# Check if alpha, cl, cd, cm are the same length
check_lengths(aero_data)






def transform_aero_data(aero_data):
    """
    Transforms the aerodynamic data into a dictionary format where CL, CD, and CM are dictionaries
    with the angle of attack as keys and coefficients as values.
    """
    # Initialize the new structure
    transformed_data = {
        "CL": {},
        "CD": {},
        "CM": {}
    }

    # Use zip to iterate over all lists simultaneously
    for alpha, cl, cd, cm in zip(aero_data["alpha"], aero_data["cl"], aero_data["cd"], aero_data["cm"]):
        transformed_data["CL"][alpha] = cl
        transformed_data["CD"][alpha] = cd
        transformed_data["CM"][alpha] = cm

    return transformed_data

# Transform the aerodynamic data to the desired format
transformed_aero_data = transform_aero_data(aero_data)
# 01 print(transformed_aero_data)






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

wing_aero_data = transformed_aero_data

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
alpha_wing = 180.1  # Example angle of attack for the wing
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
angle_range = np.linspace(-180, 180, 1000)  # From -10 to 20 degrees, 300 points
# 03 print(angle_range)

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


