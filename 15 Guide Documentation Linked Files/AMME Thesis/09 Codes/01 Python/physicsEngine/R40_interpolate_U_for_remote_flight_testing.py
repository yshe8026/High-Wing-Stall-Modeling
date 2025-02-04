import numpy as np
from scipy.interpolate import interp1d
import os

# This is used as the development concept of "Richard Taped Control Input Module" in "flightModel.py"

current_work_directory_path = os.getcwd()
known_part_of_path_to_remove = 'Resources\plugins\PythonPlugins'
xplane11_path = current_work_directory_path.replace(known_part_of_path_to_remove, "")
# Print Python work directory
print(current_work_directory_path)
# Print xplane 11 directory on this machine
print(xplane11_path) 
# Path to the directory where the Python file resides
script_dir = os.path.dirname(os.path.abspath(__file__))
print(script_dir)
# Full path to the current Python file
file_path = os.path.abspath(__file__)
print(file_path)

U_data_matrix_taped = np.load(xplane11_path + 'Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\U_data_matrix_taped.npy')
t_taped = np.load(xplane11_path + 'Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\t_taped.npy')

t = 15  # Example time for interpolation

# Create an interpolation function for each column in U_data_matrix
interpolated_U = interp1d(t_taped, U_data_matrix_taped, axis=0, kind='linear', fill_value='extrapolate')

# Get the interpolated value of U at time t
U_at_t = interpolated_U(t)

print("Interpolated U at t =", t, "is", U_at_t)
