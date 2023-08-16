"""
Instructions:
1. Run Weekday4_preprocessing.py first
2. Run this script
"""

# Imports ====================================================================
import numpy as np
from datetime import datetime
from time import time

from tools.local_emf import local_emf
from tools.plotting import TimeplotPyQt
from tools.GaussmeterLib import (
    EulerRotation,
    read_data,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
)

time0 = time()
# Rotation matrix from geographic frame to cage frame
R_G2C = EulerRotation().rz(-23, deg=True)
# Rotation matrix from sensor frame to cage frame
R_S2C = EulerRotation().rx(-180, deg=True)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
local_emf = np.dot(R_G2C, local_emf())

# Data processing ============================================================

# Filename
filename = "./data/Weekday4_2023-07-13_08.29.31_SPARSE25.dat"

# Import data
data = read_data(filename, header=True)

# Apply Savitzky-Golay smoothing
data = data_savgol_filter(data, 51, 4)

# Rotate data from the sensor frame S into the cage frame C
data = data_rotate(data, R_S2C)

# Normalize the data (in C frame) with respect to the EMF (in C frame):
data = data_add_vector(data, -local_emf)

# Set up plot ================================================================
plot_title = """
<h3><b> Magnetic field components of a weekday test performed from {} to {}.
</b> <br> Data fetched from file:  <samp> {} </samp> </h3>
""".format(
    datetime.utcfromtimestamp(data[0][0]).strftime("%d/%m/%y %H:%M:%S"),
    datetime.utcfromtimestamp(data[0][-1]).strftime("%d/%m/%y %H:%M:%S"),
    filename
)
side_label = "<h3>Magnetic field component B (uT)<h3>"

plot_filename = "./figures/Weekday4_timeplot.png"

plot_object = TimeplotPyQt()
plot_object.set_plot_title(plot_title)
plot_object.set_side_label(side_label)
plot_object.set_window_size(1920, 1080)

plot_object.save_plot(plot_filename)

# Generate plot
plot_object.timeplot_pyqtgraph(data)

# 3axis plot (separate to ensure that the class object is fresh)
plot_object2 = TimeplotPyQt()
plot_object2.set_plot_title(plot_title)
plot_object2.set_side_label(side_label)
plot_object2.set_window_size(1920, 1080)

plot_object2.save_plot(plot_filename[:-4]+"3axis"+plot_filename[-4:])

# Generate plot
plot_object2.timeplot_3axis_pyqtgraph(data)
