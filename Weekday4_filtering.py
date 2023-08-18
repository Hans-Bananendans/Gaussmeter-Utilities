"""
Instructions:
1. Run Weekday4_preprocessing.py first if these files do not exist yet
2. Then run Weekday4_preprocessing2.py
3. Run this script
"""

# Imports ====================================================================
import numpy as np
from time import time

from tools.local_emf import local_emf
from tools.plotting import TimeplotPyQt, SpectralPlotPyQt
from tools.GaussmeterLib import (
    EulerRotation,
    read_data,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
    data_notch_filter,
    data_rfft,
    data_normal_distribution
)

# Select plotset
plotset = {
    "unfiltered": {
        "skip": (True, True),
        "plot_title_insert": "UNFILTERED",
        "filename_insert": "_unfiltered"
    },
    "sgf": {        
        "skip": (False, True),
        "plot_title_insert": "SGF",
        "filename_insert": "_sgf"
    },
    # Savitsky-Golay filter and 
    "sgf_notch": {
        "skip": (False, False),
        "plot_title_insert": "SGF & notch",
        "filename_insert": "_sgf_notch"
    }
}

_selection = "sgf"

time0 = time()
# Rotation matrix from geographic frame to cage frame
R_G2C = EulerRotation().rz(-23, deg=True)
# Rotation matrix from sensor frame to cage frame
R_S2C = EulerRotation().rx(-180, deg=True)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
local_emf = np.dot(R_G2C, local_emf())


# Data processing ============================================================

# Filename
filename_d = "./data/Weekday4_2023-07-13_08.29.31_sample_day.dat"
filename_n = "./data/Weekday4_2023-07-13_08.29.31_sample_night.dat"

# Import data
data_d = read_data(filename_d, header=True)
data_n = read_data(filename_n, header=True)

if not plotset[_selection]["skip"][0]:
    # Apply Savitzky-Golay smoothing
    data_d = data_savgol_filter(data_d, 51, 4)
    data_n = data_savgol_filter(data_n, 51, 4)

if not plotset[_selection]["skip"][1]:
    # Apply notch filter at 50 Hz
    nw = 1  # Approximate width of peak in Hz based on spectral plot
    fs = 2500
    fn = (50, 100, 150, 200, 250, 500, 750, 1000)
    for fn_i in fn:
        data_d = data_notch_filter(data_d, fn_i, fs, nw=nw)
        data_n = data_notch_filter(data_n, fn_i, fs, nw=nw)

# Rotate data from the sensor frame S into the cage frame C
data_d = data_rotate(data_d, R_S2C)
data_n = data_rotate(data_n, R_S2C)

# Normalize the data (in C frame) with respect to the EMF (in C frame):
data_d = data_add_vector(data_d, -local_emf)
data_n = data_add_vector(data_n, -local_emf)

# Since notch filter performance is poor at edges, cut away 10% of data at
# start and end of dataset when analysing the distribution of the filtered data
data_d_truncated = []
data_n_truncated = []

for i in (0, 1, 2, 3):
    data_d_truncated.append(data_d[i][round(0.1*len(data_d[0])):round(0.1005*len(data_d[0]))])
    data_n_truncated.append(data_n[i][round(0.1*len(data_n[0])):round(0.1005*len(data_n[0]))])

# Calculate standard deviation of data
mean_data_d, sd_data_d = data_normal_distribution(data_d_truncated)
mean_data_n, sd_data_n = data_normal_distribution(data_n_truncated)

print(
    "sd_day{}:   {} {} {}".format(
        plotset[_selection]["filename_insert"],
        round(sd_data_d[0], 4),
        round(sd_data_d[1], 4),
        round(sd_data_d[2], 4),
    )
)
print(
    "sd_night{}: {} {} {}".format(
        plotset[_selection]["filename_insert"],
        round(sd_data_n[0], 4),
        round(sd_data_n[1], 4),
        round(sd_data_n[2], 4),
    )
)

# Spectral analysis
[f, x_t, y_t, z_t] = data_rfft(data_d)
data_f_d = data_rfft(data_d)
data_f_n = data_rfft(data_n)


# Timeplot ===================================================================
timeplot_title_d = """
<h3><b> Short daytime sample of weekday test - {}
</b> <br> Data fetched from file:  <samp> {} </samp> </h3>
""".format(plotset[_selection]["plot_title_insert"], filename_d)

timeplot_title_n = """
<h3><b> Short nighttime sample of weekday test - {}
</b> <br> Data fetched from file:  <samp> {} </samp> </h3>
""".format(plotset[_selection]["plot_title_insert"], filename_n)

side_label = "<h3>Magnetic field component B (uT)<h3>"

timeplot_filename_d = "./figures/Weekday4_day_timeplot{}.png".format(plotset[_selection]["filename_insert"])
timeplot_filename_n = "./figures/Weekday4_night_timeplot{}.png".format(plotset[_selection]["filename_insert"])

# 3axis plot (separate to ensure that the class object is fresh)
# Generate daytime plot
plot_object = TimeplotPyQt()
plot_object.set_plot_title(timeplot_title_d)
plot_object.set_side_label(side_label)
plot_object.set_window_size(1920, 1080)
plot_object.save_plot(timeplot_filename_d)
plot_object.timeplot_3axis_pyqtgraph(data_d)

# Generate nighttime plot
plot_object = TimeplotPyQt()
plot_object.set_plot_title(timeplot_title_n)
plot_object.set_side_label(side_label)
plot_object.set_window_size(1920, 1080)
plot_object.save_plot(timeplot_filename_n)
plot_object.timeplot_3axis_pyqtgraph(data_n)

# Spectral Plot ==============================================================
spectralplot_title_d = """
<h3><b> Frequency spectrum of a short daytime sample of weekday test - {}
</b> <br> Data fetched from file:  <samp> {} </samp> </h3>
""".format(plotset[_selection]["plot_title_insert"], filename_d)

spectralplot_title_n = """
<h3><b> Frequency spectrum of a short nighttime sample of weekday test - {}
</b> <br> Data fetched from file:  <samp> {} </samp> </h3>
""".format(plotset[_selection]["plot_title_insert"], filename_n)

spectralplot_filename_d = "./figures/Weekday4_day_spectralplot{}.png".format(plotset[_selection]["filename_insert"])
spectralplot_filename_n = "./figures/Weekday4_night_spectralplot{}.png".format(plotset[_selection]["filename_insert"])

# Generate daytime plot
spectralplot = SpectralPlotPyQt()
spectralplot.set_window_size(1200, 800)
spectralplot.set_pen_alpha(1)
spectralplot.set_plot_title(spectralplot_title_d)
spectralplot.set_xlabel("Frequency [Hz]")
spectralplot.set_ylabel("Field magnitude [\u03BCT / \u221AHz]")
spectralplot.save_plot(spectralplot_filename_d)
spectralplot.spectralplot_pyqtgraph(data_f_d)

# Generate nighttime plot
spectralplot = SpectralPlotPyQt()
spectralplot.set_window_size(1200, 800)
spectralplot.set_pen_alpha(1)
spectralplot.set_plot_title(spectralplot_title_n)
spectralplot.set_xlabel("Frequency [Hz]")
spectralplot.set_ylabel("Field magnitude [\u03BCT / \u221AHz]")
spectralplot.save_plot(spectralplot_filename_n)
spectralplot.spectralplot_pyqtgraph(data_f_n)