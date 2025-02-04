# import matplotlib.pyplot as plt

# # # Try to enable LaTeX rendering
# # try:
# #     plt.rc('text', usetex=True)
# #     print("LaTeX is available for Matplotlib.")
# # except Exception as e:
# #     print(f"An error occurred: {e}")

# import numpy as np
# import matplotlib.pyplot as plt
# import os


# current_work_directory_path = os.getcwd()
# known_part_of_path_to_remove = 'Resources\plugins\PythonPlugins'
# xplane11_path = current_work_directory_path.replace(known_part_of_path_to_remove, "")
# # Print Python work directory
# print(current_work_directory_path)
# # Print xplane 11 directory on this machine
# print(xplane11_path) 
# # Path to the directory where the Python file resides
# script_dir = os.path.dirname(os.path.abspath(__file__))
# print(script_dir)
# # Full path to the current Python file
# file_path = os.path.abspath(__file__)
# print(file_path)

# # flight_data_matrix = np.load(xplane11_path + 'Resources\\plugins\\PythonPlugins\\physicsEngine\\F01_flight_data_recorder_1\\flight_data_matrix.npy')


# # #Direct input 
# # # plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
# # plt.rcParams['text.usetex'] = True
# # # #Options
# # # params = {'text.usetex' : True,
# # #           'font.size' : 11,
# # #           'font.family' : 'lmodern',
# # #           'text.latex.unicode': True,
# # #           }
# # # plt.rcParams.update(params) 

# fig = plt.figure()

# #You must select the correct size of the plot in advance
# fig.set_size_inches(3.54,3.54) 

# plt.plot([1,2,3,4])
# plt.xlabel("Excitation-Energy")
# plt.ylabel("Intensität")
# plt.savefig(xplane11_path + "richard_quality_plots\\graph.pdf", 
#             #This is simple recomendation for publication plots
#             dpi=1000, 
#             # Plot will be occupy a maximum of available space
#             bbox_inches='tight', 
#             )

# Failed
############################################################################################################################

# from cycler import cycler

# import matplotlib.pyplot as plt
# import numpy as np

# import matplotlib as mpl

# mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.linestyle'] = '--'
# data = np.random.randn(50)
# mpl.rcParams['axes.prop_cycle'] = cycler(color=['r', 'g', 'b', 'y'])
# plt.plot(data)  # first color is red
# plt.show()

# Failed
############################################################################################################################
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

# plt.style.use(['science','ieee'])
# plt.rcParams.update({'figure.dpi': '100'})

# plt.style.use(['science','nature'])

# with plt.style.context('science'):
#     plt.figure()
#     plt.plot(x, y)
#     plt.show()

# plt.style.use(['science','no-latex'])
# plt.style.use(['science','grid'])
# plt.style.use(['science','scatter'])
# plt.style.use(['science','notebook'])

plt.plot([1,2,3,4])
plt.xlabel("Excitation-Energy")
plt.ylabel("Intensität")

plt.show()

# Successful
############################################################################################################################