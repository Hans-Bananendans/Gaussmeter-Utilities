# Imports ====================================================================
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from datetime import datetime
from time import time

from local_emf import local_emf
from TimeplotPyQt import TimeplotPyQt
from GaussmeterAnalysis import (
    EulerRotation,
    read_data,
    data_normal_distribution,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
    time_plot,
    VectorPlot
)

time0 = time()
# Rotation matrix from geographic frame to cage frame
R_G2C = EulerRotation().rz(-23, deg=True)
# Rotation matrix from sensor frame to cage frame
R_S2C = EulerRotation().rx(-180, deg=True)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
local_emf = np.dot(R_G2C, local_emf())

# Data processing ============================================================

# Filenames in chronological order
filenames = (
    "SpatialX+50Y+0_2023-07-12_06.49.34_SPARSE25.dat",
    "SpatialX+50Y-50_2023-07-12_07.38.19_SPARSE25.dat",
    "SpatialX0Y-50_2023-07-12_08.37.24_SPARSE25.dat",
    "SpatialX-50Y-50_2023-07-12_09.46.26_SPARSE25.dat",
    "SpatialX-50Y0_2023-07-12_10.35.41_SPARSE25.dat",
    "SpatialX-50Y+50_2023-07-12_11.29.49_SPARSE25.dat",
    "SpatialX0Y+50_2023-07-12_12.21.00_SPARSE25.dat",
    "SpatialX+50Y+50_2023-07-12_13.18.39_SPARSE25.dat",
    "SpatialX0Y0_2023-07-12_14.15.31_SPARSE25.dat"  # Control test
)

# Coordinates in [m] of test locations relative to center of HHC
origins = (
    ( 0.5,    0, 0),
    ( 0.5, -0.5, 0),
    (   0, -0.5, 0),
    (-0.5, -0.5, 0),
    (-0.5,    0, 0),
    (-0.5,  0.5, 0),
    (   0,  0.5, 0),
    ( 0.5,  0.5, 0),
    (   0,    0, 0),
)

# Pre-allocate storage lists
mean_data_n = [0.]*len(filenames)
sd_data_n = [0.]*len(filenames)

for i, filename in enumerate(filenames):
    # Readout files
    data = read_data(filename, header=True)

    # Apply Savitzky-Golay smoothing
    data = data_savgol_filter(data, 51, 4)

    # Rotate data from the sensor frame S into the cage frame C
    data = data_rotate(data, R_S2C)

    # Normalize the data (in C frame) with respect to the EMF (in C frame):
    data = data_add_vector(data, -local_emf)

    # Distribution properties of normalized disturbance datapoints in C frame
    mean_data_n[i], sd_data_n[i] = data_normal_distribution(data)

# Easy Copy-Paste readout
easy_copy_paste_readout = False
if easy_copy_paste_readout:
    for j in (0, 1, 2):
        for i in range(len(filenames)):
            print(round(mean_data_n[i][j], 3))

        print(" ")

    for j in (0, 1, 2):
        for i in range(len(filenames)):
            print(round(sd_data_n[i][j], 3))

        print(" ")

    for i in range(len(filenames)):
        print(round(np.linalg.norm(mean_data_n[i]), 3))

    print(" ")

# Vector plot
plot_title = """
    Average disturbance vectors at different locations in the HHC. \n
    Measured on 12-07-2023, measurements each 45+ minutes"""

vectorplot = VectorPlot()
vectorplot.set_plot_title(plot_title)
vectorplot.autoplot(tripod=True, coils=True, walls=True,
                    table=True, emf_vector=True)
for i in range(len(mean_data_n)):
    vectorplot.plot_vector(mean_data_n[i], origin=origins[i], color="black",
                           linewidth=3, alr=0.2)
vectorplot.show()

print("Elapsed time:", time()-time0)
