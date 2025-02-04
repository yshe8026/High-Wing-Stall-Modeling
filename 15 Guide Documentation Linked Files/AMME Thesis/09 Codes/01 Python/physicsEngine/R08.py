# calculate the time rate of change of X_flow_sep
dX_flow_sep_by_dt = 1
dt = 0.01

#==================================================================================================
# script for calculating X_flow_sep from dX_flow_sep_by_dt and dt
#==================================================================================================
X_flow_sep_file_path = "physicsEngine\\X_flow_sep.txt"
with open(X_flow_sep_file_path, "r") as file:
    X_flow_sep = float(file.read().strip()) + dX_flow_sep_by_dt * dt
# Write the new value back to the file without formatting to preserve accuracy
with open(X_flow_sep_file_path, "w") as file:
    file.write(str(X_flow_sep))
# print("Current value:", X_flow_sep)
#==================================================================================================