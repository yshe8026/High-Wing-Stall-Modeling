import numpy as np

# Initialize necessary global variables
global flight_data_matrix_as_list
global whether_flight_data_recorder_is_on
global global_time
global global_dt

flight_data_matrix_as_list = []
whether_flight_data_recorder_is_on = 1  # Start with the recorder on
global_time = 0  # Start time
global_dt = 0.01  # Time step increment

# Assume that X is generated at each time step
def generate_X_vector():
    return np.random.rand(12)

# Richard Flight Data Recorder Module (place 1 of 2)
whether_to_add_richard_flight_data_recorder = 1

# Simulate the process
while global_time < 30:  # Simulate for 30 seconds
    # Generate the X vector
    X = generate_X_vector()

    if whether_to_add_richard_flight_data_recorder == 1:
        # Check whether the flight data recorder is still on
        if whether_flight_data_recorder_is_on == 1:

            # At which moment in time to stop recording data? (Record data for how long since pressing '9'?)
            global_time_to_stop_recording_data = 0.1  # sec

            # Simulate the process of recording X at each timestep
            if global_time < global_time_to_stop_recording_data:  # Stopping condition

                # Create a new flight_data_row each iteration to avoid overwriting previous data
                flight_data_row = np.zeros(13)
                flight_data_row[0] = global_time
                flight_data_row[1:13] = X[0:12]  # Adjusted to fill the entire row

                flight_data_matrix_as_list.append(flight_data_row)
            else:
                # Convert the list of vectors to a NumPy array (matrix)
                flight_data_matrix = np.array(flight_data_matrix_as_list)

                if whether_to_add_richard_flight_data_recorder == 1:
                    # Assuming flight_data_matrix is your NumPy array
                    np.save('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\flight_data_matrix.npy', flight_data_matrix)

                    # Turn the flight data recorder off via a global switch
                    whether_flight_data_recorder_is_on = 0

    # Increment the global time
    global_time += global_dt

# Check the results
print("Flight data recorder test completed.")
print(f"Number of recorded time steps: {len(flight_data_matrix_as_list)}")

# Optional: Load the saved .npy file to verify its contents
loaded_flight_data_matrix = np.load('C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\flight_data_matrix.npy')
print("Loaded matrix shape:", loaded_flight_data_matrix.shape)
# print("First row of the loaded matrix:", loaded_flight_data_matrix[0])

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
