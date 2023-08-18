"""
Instructions:
1. Run the file splitter tool first to generate test segments
2. If needed (and possible), run desample tool to reduce processing time
3. Set the correct breakpoint so the filenames get chopped up correcly (probably 13)
4. Fill the number of generated segments into generate_filenames() in line 52
5. Run this script
"""

# Imports ====================================================================
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from datetime import datetime
from time import time

from tools.local_emf import local_emf
from tools.plotting import TimeplotPyQt, VectorPlot
from tools.GaussmeterLib import (
    EulerRotation,
    read_data,
    data_normal_distribution,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
    multipart_filenames
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
filename = "./data/Weekday2500Hz_LightsOff_2023-07-06_19.12.02_SPARSE25.dat"

def generate_filenames(filename_base, n_files, breakpoint: int = 4):
    filenames = []
    for i in range(1, n_files+1):
        filenames.append(
            filename_base[:-breakpoint]
            + "_s{}of{}".format(i, n_files)
            + filename_base[-breakpoint:]
        )
    return filenames


filenames = multipart_filenames(filename_base, 25, breakpoint=9+4)

print(filenames)

# Pre-allocate storage lists
mean_data_n = [0.]*len(filenames)
sd_data_n = [0.]*len(filenames)

print("Starting data processing...")
for i, filename in enumerate(filenames):

    print("Starting segment {} of {} (time: {} s)"
          .format(i+1, len(filenames), time()-time0))
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
easy_copy_paste_readout = True
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
print("Starting vector plot")

plot_title = """
    Average disturbance vectors at different locations in the HHC. \n
    Measured on 12-07-2023, measurements each 45+ minutes"""

vectorplot = VectorPlot()
vectorplot.set_plot_title(plot_title)
vectorplot.autoplot(tripod=True, coils=True, walls=True,
                    table=True)
for i in range(len(mean_data_n)):
    vectorplot.plot_vector(mean_data_n[i], color="black",
                           linewidth=1.5, alr=0.1)
vectorplot.show()

print("Elapsed time:", time()-time0)
