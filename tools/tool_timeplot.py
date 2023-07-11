# Imports ====================================================================
import numpy as np
# import pyqtgraph as pg

from time import time
from datetime import datetime
from copy import deepcopy

from GaussmeterAnalysis import (
    # TimeplotPyQt,
    EulerRotation,
    read_data,
    vector_normal_distribution,
    data_savgol_filter,
    zipBA,
    data_rotate,
    data_add_vector,
    time_plot,
    time_plot2,
)



# EMF vector =================================================================
# R_G2C = Rz(-23, deg=True)
# R_S2C = Rx(-180, deg=True)

R_G2C = EulerRotation().rz(-23, deg=True)
R_S2C = EulerRotation().rx(-180, deg=True)

# EMF vector at AE building (WMM-2020 model)
B_EMF_WMM2020 = np.array([0.7216, 19.1629, -45.4592])
# EMF vector at AE building (IGRF2020 model)
B_EMF_IGRF2020 = np.array([0.7313, 19.1870, -45.4557])

# Take average of the two as baseline B_EMF:
B_EMF_G = 0.5 * (B_EMF_WMM2020 + B_EMF_IGRF2020)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
B_EMF = np.dot(R_G2C, B_EMF_G)


# Data processing ============================================================
filename = "WeekdayTest_2023-07-03_10.41.03_SPARSE50.dat"

data_r = read_data(filename, header=True)

data_r = data_savgol_filter(data_r, 51, 4)


# data = []
# for i in range(len(data_r[0])):
#     data.append([
#         data_r[0][i],
#         np.dot(R_S2C, np.array([data_r[1][i],
#                                 data_r[2][i],
#                                 data_r[3][i]]))
#     ])

# print(sum(data[1]), sum(data[2]), sum(data[3]))

data = data_rotate(data_r, R_S2C)

# print(sum(data[1]), sum(data[2]), sum(data[3]))



# test_data = [
#     [1, 2, 3, 4, 5, 6, 7],
#     [0, 1, 1, 0, 0, 0, 0],
#     [0, 0, 0, 1, 1, 0, 0],
#     [0, 0, 0, 0, 0, 1, 1],
# ]
# test_vector = np.array([10, 10, 10])
# # test_data2 = data_add_vector(test_data, -test_vector)
# test_data3 = data_rotate(test_data, EulerRotation().rz(-90, deg=True))
#
# for line in test_data3:
#     print(line)

# Put all the vectors in a list
# data_vectors = []
# for d in data:
#     data_vectors.append(d[1])

# mean_data, sd_data = vector_normal_distribution(data_vectors)
# print("C-frame data: mean, sd:", mean_data.round(3), sd_data.round(3))

# Normalize the data (in C frame) with respect to the   EMF (in C frame):
# data_n = deepcopy(data)
#
# for i in range(len(data_n)):
#     data_n[i][1] = data_n[i][1] - B_EMF
#
# data_n = zipBA(data_n)

data_n = data_add_vector(data, -B_EMF)

plot_title = """
<h3><b> Magnetic field components of a weekday test performed from {} to {}.
</b> <br> Data fetched from file:  <samp> {} </samp> </h3>
""".format(datetime.utcfromtimestamp(data_n[0][0]).strftime("%d/%m/%y %H:%M:%S"),
           datetime.utcfromtimestamp(data_n[0][-1]).strftime("%d/%m/%y %H:%M:%S"),
           filename)

side_label = "<h3>Magnetic field component B (uT)<h3>"


# plot_object = TimeplotPyQt()
# plot_object.set_plot_title(plot_title)
# plot_object.set_side_label(side_label)

# plot_object.timeplot_pyqtgraph(data_n)


persistent = time_plot2(data_n)
