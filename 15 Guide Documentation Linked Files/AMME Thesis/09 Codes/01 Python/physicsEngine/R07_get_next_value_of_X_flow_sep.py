#==================================================================================================
# Function for calculating X_flow_sep from dX_flow_sep_by_dt and dt
# Input:
# - dX_flow_sep_by_dt
# - dt
# Output:
# - X_flow_sep
#==================================================================================================
def R07_get_next_value_of_X_flow_sep(dX_flow_sep_by_dt, dt):
    X_flow_sep_file_path = "physicsEngine\\X_flow_sep.txt"

    try:
        # Try to open the file and read the current value
        with open(X_flow_sep_file_path, "r") as file:
            X_flow_sep = float(file.read().strip()) + dX_flow_sep_by_dt * dt
    except (FileNotFoundError, ValueError):
        # If the file doesn't exist or the value is not a float, start from the increment value
        X_flow_sep = dX_flow_sep_by_dt * dt
    
    # Write the new value back to the file without formatting to preserve accuracy
    with open(X_flow_sep_file_path, "w") as file:
        file.write(str(X_flow_sep))
    
    return X_flow_sep

# # Example usage: Incrementing by 0.001
# variable = R07_get_next_value_of_X_flow_sep(0.1,0.01)
# print("Current value:", variable)