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
print(aero_data)

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
print(transformed_aero_data)