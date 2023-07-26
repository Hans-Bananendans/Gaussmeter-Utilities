# Place to test VectorPlot class

# Imports ====================================================================
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from datetime import datetime

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

coil_sides = [1.85, 1.95, 2.05]                 # [m]
coil_spacings = [1.0073, 1.0618, 1.1162]        # [m]

vectorplot = VectorPlot()
vectorplot.plot_global_tripod()
vectorplot.plot_hhc_coils(coil_sides, coil_spacings)
vectorplot.plot_emf()
vectorplot.plot_x_wall(3, 2, -1.5, 0, -0.1)
vectorplot.plot_y_wall(3, 2,  0, 1.5, -0.1)
vectorplot.plot_table()

vectorplot.show()
