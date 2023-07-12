# Imports ====================================================================
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
# import mpl_toolkits.mplot3d as mp3d
from scipy.fft import rfft, fftshift

from cp import (
    Vertex,
    Face,
    plot_face,
    plot_arrow,
    plot_global_tripod,
)

from local_emf import local_emf
from TimeplotPyQt import TimeplotPyQt
from SpectralplotPyQt import SpectralplotPyQt
from GaussmeterAnalysis import (
    EulerRotation,
    data_rfft,
    read_data,
    vector_normal_distribution,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
    time_plot,
    create_hhc_elements,
    create_lab_walls,
)

# Rotation matrix from geographic frame to cage frame
R_G2C = EulerRotation().rz(-23, deg=True)
# Rotation matrix from sensor frame to cage frame
R_S2C = EulerRotation().rx(-180, deg=True)

# Rotate EMF vector from geographic frame (G) into cage frame (C)
local_emf = np.dot(R_G2C, local_emf())

# Data processing ============================================================
filename = "2500HzTestLightOn_2023-07-05_11.04.28.dat"

# Import data
data = read_data(filename, header=True)

# Apply Savitzky-Golay smoothing
data = data_savgol_filter(data, 51, 4)

# Rotate data from the sensor frame S into the cage frame C
data = data_rotate(data, R_S2C)

# Normalize the data (in C frame) with respect to the EMF (in C frame):
data = data_add_vector(data, -local_emf)

# Obtain one-sided FFT data
# [[x_t1, x_f1], [y_t1, y_f1], [z_t1, z_f1]] = data_rfft(data)
[f, x_t, y_t, z_t] = data_rfft(data)
data_f = data_rfft(data)


spectralplot = SpectralplotPyQt()
spectralplot.set_window_size(1920, 1080)
spectralplot.set_pen_alpha(1)
spectralplot.spectralplot_pyqtgraph(data_f)

# # Plot Spectrogram LightOn XYZ
# fig1 = plt.figure(figsize=(8, 8))
# gs = GridSpec(1,1)
# ax = fig1.add_subplot(gs[0,0])
# # ax.set_ylim(ylim[0], ylim[1])
#
#
# filter_dc: int = 0
#
# ax.plot(f[filter_dc:], np.abs(x_t[filter_dc:]), color='r', label="X", alpha=0.5)
# ax.plot(f[filter_dc:], np.abs(y_t[filter_dc:]), color='g', label="Y", alpha=0.4)
# ax.plot(f[filter_dc:], np.abs(z_t[filter_dc:]), color='b', label="Z", alpha=0.3)
#
# ax.set_title("Spectrogram of magnetic measurements")
# ax.set_xlabel("f [Hz]")
# ax.set_ylabel("Linear magnitude")
# ax.legend()
# ax.grid(True)
# # ax.set_xscale('log')
# ax.set_yscale('log')
# plt.show()


# # Set up plot ================================================================
# plot_title = """
# <h3><b> Magnetic field components of a weekday test performed from {} to {}.
# </b> <br> Data fetched from file:  <samp> {} </samp> </h3>
# """.format(
#     datetime.utcfromtimestamp(data[0][0]).strftime("%d/%m/%y %H:%M:%S"),
#     datetime.utcfromtimestamp(data[0][-1]).strftime("%d/%m/%y %H:%M:%S"),
#     filename
# )
# side_label = "<h3>Magnetic field component B (uT)<h3>"
#
# plot_filename = "WeekdayTest1_timeplot.png"
#
# plot_object = TimeplotPyQt()
# plot_object.set_plot_title(plot_title)
# plot_object.set_side_label(side_label)
# plot_object.set_window_size(1920, 1080)
#
# plot_object.save_plot(plot_filename)
#
# # Generate plot
# plot_object.timeplot_pyqtgraph(data)
