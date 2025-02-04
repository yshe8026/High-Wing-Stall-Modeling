import numpy as np

'''Each Row in the Flight Data Matrix is a data entry at a specific momemt in time'''

#######################################################################################################
# Flight Data Recorder Module
#------------------------------------------------------------------------------------------------------

# Initialize an empty list to store data vectors
flight_data_matrix_as_list = []

# global global_time
# global global_dt

# At when moment to stop?
global_time_to_stop_recording_data = 0.1

global_time = 0
global_dt = 0.01

# Simulate the process of recording X at each timestep
while global_time < global_time_to_stop_recording_data:  # Replace with actual stopping condition
    X = np.random.rand(12)  # Example X vector, replace with actual computation
    
    # Create a new flight_data_row each iteration to avoid overwriting previous data
    flight_data_row = np.zeros(13)
    flight_data_row[0] = global_time
    flight_data_row[1:13] = X[0:12]  # Adjusted to fill the entire row

    flight_data_matrix_as_list.append(flight_data_row)
    global_time += global_dt

# Convert the list of vectors to a NumPy array (matrix)
flight_data_matrix = np.array(flight_data_matrix_as_list)

# Assuming flight_data_matrix is your NumPy array
np.save('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\flight_data_matrix.npy', flight_data_matrix)

#######################################################################################################

#==============================================
whether_to_print_full_matrix = 1
#----------------------------------------------
if whether_to_print_full_matrix == 1:
    # Set print options to avoid truncation
    np.set_printoptions(threshold=np.inf)
    # Set print options to avoid wrapping
    np.set_printoptions(linewidth=200)

# Print the full array
print('Flight Data Matrix:')
print(flight_data_matrix)

if whether_to_print_full_matrix == 1:
    # Reset print options to default (truncated output)
    np.set_printoptions(threshold=1000)  # 1000 is the default threshold
    # Set print options to avoid wrapping
    np.set_printoptions(linewidth=100)
#==============================================

loaded_flight_data_matrix = np.load('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\flight_data_matrix.npy')

#==============================================
whether_to_print_full_matrix = 1
#----------------------------------------------
if whether_to_print_full_matrix == 1:
    # Set print options to avoid truncation
    np.set_printoptions(threshold=np.inf)
    # Set print options to avoid wrapping
    np.set_printoptions(linewidth=200)

# Print the full array
print('Loaded Flight Data Matrix:')
print(loaded_flight_data_matrix)

if whether_to_print_full_matrix == 1:
    # Reset print options to default (truncated output)
    np.set_printoptions(threshold=1000)  # 1000 is the default threshold
    # Set print options to avoid wrapping
    np.set_printoptions(linewidth=100)
#==============================================

print(flight_data_matrix_as_list[1][0]) # print first element in the 2nd np.ndarrary in the list