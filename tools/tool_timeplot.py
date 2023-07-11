# Imports ====================================================================
import numpy as np
import pyqtgraph as pg

from datetime import datetime


from GaussmeterAnalysis import (
    TimeplotPyQt,
    EulerRotation,
    read_data,
    vector_normal_distribution,
    data_savgol_filter,
    zipAB,
)


# Helper functions =========================================================


# class TimeplotPyQt:
#     def __init__(self):
#         self.app = pg.mkQApp("Timeplot")
#         self.view = pg.GraphicsView()
#         self.gl = pg.GraphicsLayout()
#
#         self.view.setCentralItem(self.gl)
#
#         # Default settings
#         self.gl.setBorder(100, 100, 100)
#         self.view.setWindowTitle("Timeplot")
#         self.view.resize(800, 600)
#
#         self.gl.layout.setColumnPreferredWidth(0, 25)
#         self.gl.layout.setColumnPreferredWidth(1, 1200)
#
#         pg.setConfigOptions(antialias=True)  # Antialiasing for nicer plots
#
#         self.set_defaults()
#
#
#     def set_defaults(self):
#         self.light_theme = False
#
#         self.padding = [1.1, 1.0]  # Padding factor
#         self.force_identical_scale = True
#         self.pen_alpha = 0.8
#
#         pen_colors = ((255,   0,   0),
#                       (  0, 255,   0),
#                       (  0,  75, 255))
#         self.generate_pens(pen_colors)
#
#         # Colours of the draggable lines, depending on theme
#         if self.light_theme:
#             self.label_fill = (255, 255, 255, 180)  # White label backgrounds
#             self.label_text_color = (100, 100, 100)  # Darker label text
#         else:
#             self.label_fill = (0, 0, 0, 180)
#             self.label_text_color = (200, 200, 100)
#
#         # Background colour:
#         if self.light_theme:
#             self.view.setBackground("white")  # White background
#         else:
#             self.view.setBackground("black")
#
#         self.plot_title = "Plot Title"
#         self.side_label = "Side label"
#
#         self.grids = True
#
#
#     def set_light_theme(self):
#         self.light_theme = True
#
#     def set_dark_theme(self):
#         self.light_theme = False
#
#     def generate_pens(self, pc):
#         self.pen_rgb = pc
#
#         self.pen_rgba = (
#             (pc[0][0], pc[0][1], pc[0][2], self.pen_alpha*100),
#             (pc[1][0], pc[1][1], pc[1][2], self.pen_alpha*100),
#             (pc[2][0], pc[2][1], pc[2][2], self.pen_alpha*100),
#         )
#
#     def set_window_size(self, w: int, h: int):
#         self.view.resize(w, h)
#
#     def set_window_title(self, window_title: str):
#         self.view.setWindowTitle(window_title)
#
#     def set_padding(self, padding_h: float, padding_v: float):
#         self.padding = [padding_h, padding_v]
#
#     def set_force_identical_scale(self, force_identical_scale: bool):
#         self.force_identical_scale = force_identical_scale
#
#     def set_pen_alpha(self, alpha: float):
#         """
#         Set alpha of the pen.
#         Important: set in a range between 0.0 - 1.0 !
#         """
#         self.pen_alpha = alpha
#         generate_pens(self, self.pen_rgb)
#
#     def set_pen_colors(self, pen_colors):
#         """
#         Set colors of the pen.
#         Must be supplied as a list of listed rgb values, e.g.:
#         pen_colors = ((255,   0,   0),
#                       (  0, 255,   0),
#                       (  0,   0, 255))
#         """
#         self.generate_pens(pen_colors)
#
#     def set_label_fill(self, rgba):
#         """
#         Set a label fill RGBA, for example: (255, 255, 255, 180)
#         """
#         self.label_fill = rgba
#
#     def set_label_text_color(self, rgb):
#         """
#         Set a label text RGB, for example: (255, 255, 255)
#         """
#         self.label_text_color = rgb
#
#     def set_background(self, rgb):
#         self.view.setBackground(rgb)
#
#     def set_plot_title(self, plot_title: str):
#         self.plot_title = plot_title
#
#     def set_side_label(self, side_label: str):
#         self.side_label = side_label
#
#     def show_grids(self):
#         self.grids = True
#
#     def hide_grids(self):
#         self.grids = False
#
#     def timeplot_pyqtgraph(self, data: list):
#         # TODO: Plotting stuff
#
#         # # [DEBUG]
#         # for data_list in data:
#         #     print(type(data_list))
#         #     print(type(data_list[0]))
#
#         # Ensure data is of the same length, and of the right type(s)
#         for data_list in (data[1], data[2], data[3]):
#             assert len(data[0]) == len(data_list)
#             assert type(data_list) in (list, np.ndarray)
#             assert type(data_list[0]) in (float, int, np.float64)
#
#         # Unpacking data for plotting
#         # self.data_t = np.array(data[0])
#         # self.data_t_date = []
#         # for i in range(len(data[0])):
#         #     self.data_t_date.append(self.data[0][i])
#         # self.data_x = np.array(data[1])
#         # self.data_y = np.array(data[2])
#         # self.data_z = np.array(data[3])
#
#         self.trange = [data[0][0], data[0][-1]]
#         self.xrange = [min(data[1]), max(data[1])]
#         self.yrange = [min(data[2]), max(data[2])]
#         self.zrange = [min(data[3]), max(data[3])]
#
#         # Force identical scale if set to True:
#         if self.force_identical_scale is True:
#             xmid = self.xrange[0] + (self.xrange[1] - self.xrange[0]) / 2
#             ymid = self.yrange[0] + (self.yrange[1] - self.yrange[0]) / 2
#             zmid = self.zrange[0] + (self.zrange[1] - self.zrange[0]) / 2
#
#             xscale = self.xrange[1] - self.xrange[0]
#             yscale = self.yrange[1] - self.yrange[0]
#             zscale = self.zrange[1] - self.zrange[0]
#
#             scale = max([xscale, yscale, zscale])
#
#             self.xrange = [xmid - scale / 2, xmid + scale / 2]
#             self.yrange = [ymid - scale / 2, ymid + scale / 2]
#             self.zrange = [zmid - scale / 2, zmid + scale / 2]
#
#         # Start constructing plot ============================================
#
#         # Title and side label
#         self.gl.addLabel(self.plot_title, col=1, colspan=2)
#         self.gl.nextRow()
#         self.gl.addLabel(self.side_label, angle=-90, rowspan=3)
#
#         # Data plots
#
#         plotlabels = ("X", "Y", "Z")
#         p_range = (self.xrange, self.yrange, self.zrange)
#         for i, data_array in enumerate((data[1], data[2], data[3])):
#
#             # Main data plot
#             p = self.gl.addPlot(
#                 title=plotlabels[i],
#                 x=data[0],
#                 y=np.array(data_array),
#                 pen=pg.mkPen(color=self.pen_rgba[i]),
#                 xmin=self.trange[0],
#                 xmax=self.trange[1],
#                 axisItems={'bottom': pg.DateAxisItem()}
#             )
#             p.setYRange(p_range[i][0], p_range[i][1])
#             p.showGrid(x=self.grids, y=self.grids)
#
#             # Horizontal nfinite line object (min)
#             infline_min = pg.InfiniteLine(
#                 angle=0,
#                 label="Min: {value:.2f} uT",
#                 pen=pg.mkPen(color=self.pen_rgb[i]),
#                 pos=p_range[i][0],
#                 movable=False,
#                 labelOpts={"position": 0.05,
#                            "color": self.pen_rgb[i],
#                            "fill": self.label_fill}
#             )
#             p.addItem(infline_min)
#
#             # Horizontal nfinite line object (max)
#             infline_max = pg.InfiniteLine(
#                 angle=0,
#                 label="Max: {value:.2f} uT",
#                 pen=pg.mkPen(color=self.pen_rgb[i]),
#                 pos=p_range[i][1],
#                 movable=False,
#                 labelOpts={"position": 0.05,
#                            "color": self.pen_rgb[i],
#                            "fill": self.label_fill}
#             )
#             p.addItem(infline_max)
#
#             # # [DEBUG]
#             # print(self.label_text_color)
#             # print(type(self.label_text_color))
#
#             # Horizontal scanline
#             infline_h = pg.InfiniteLine(
#                 angle=90,
#                 label="{value:.0f}",
#                 pos=self.trange[1],
#                 movable=True,
#                 bounds=[self.trange[0], self.trange[1]],
#                 labelOpts={"position": 0.12,
#                            "color": self.label_text_color,
#                            "fill": self.label_fill}
#             )
#             p.addItem(infline_h)
#
#             # Vertical scanline
#             infline_v = pg.InfiniteLine(
#                 angle=0,
#                 label="{value:.02f}",
#                 pos=p_range[i][1],
#                 movable=True,
#                 bounds=[p_range[i][0], p_range[i][1]],
#                 labelOpts={"position": 0.88,
#                            "color": self.label_text_color,
#                            "fill": self.label_fill}
#             )
#             p.addItem(infline_v)
#
#             self.gl.nextRow()
#
#         self.view.show()
#
#         # Executing PyQtGraph
#         pg.exec()



# EMF vector =================================================================
R_G2C = Rz(-23, deg=True)
R_S2C = Rx(-180, deg=True)

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

data = []
for i in range(len(data_r[0])):
    data.append([
        data_r[0][i],
        np.dot(R_S2C, np.array([data_r[1][i],
                                data_r[2][i],
                                data_r[3][i]]))
    ])

# Put all the vectors in a list
data_vectors = []
for d in data:
    data_vectors.append(d[1])

# mean_data, sd_data = vector_normal_distribution(data_vectors)
# print("C-frame data: mean, sd:", mean_data.round(3), sd_data.round(3))

# Normalize the data (in C frame) with respect to the   EMF (in C frame):
data_n = deepcopy(data)

for i in range(len(data_n)):
    data_n[i][1] = data_n[i][1] - B_EMF

data_n = zipBA(data_n)



plot_title = """
<h3><b> Magnetic field components of a weekday test performed from {} to {}.
</b> <br> Data fetched from file:  <samp> {} </samp> </h3>
""".format(datetime.utcfromtimestamp(data_n[0][0]).strftime("%d/%m/%y %H:%M:%S"),
           datetime.utcfromtimestamp(data_n[0][-1]).strftime("%d/%m/%y %H:%M:%S"),
           filename)

side_label = "<h3>Magnetic field component B (uT)<h3>"


plot_object = TimeplotPyQt()
plot_object.set_plot_title(plot_title)
plot_object.set_side_label(side_label)

# plot_object.timeplot_pyqtgraph(data_n)
