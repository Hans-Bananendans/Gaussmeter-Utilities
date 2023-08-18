"""
This is a script that hardcodes the results of the Weekday4 data and plots the
standard deviation of each measured axis as a function of time. The time
resolution is in hours.
"""

# Imports ====================================================================
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from datetime import datetime

# Fill up an array with 24 unix timestamps, spaced hourly, with the first one
# being at the start of the Weekday4 test, namely 2023-13-07 09:00:00.000
timestamps = []
# start = datetime(2023, 7, 13, 9, 0, 0, 0, tzinfo=timezone.utc)
start = datetime(2023, 7, 13, 9, 0, 0, 0)  # Look, I know the timezones will be
                                           # fucked, but it's to unfuck the
                                           # horizontal axis of pyqtgraph.
timestamps.append(start.timestamp())
i = 1
while i < 24:
    timestamps.append(timestamps[-1]+3600)
    i += 1
timestamps = np.array(timestamps)

# Standard deviation data (imported from "Test Results.ods")
sd_x = np.array([0.127, 0.108, 0.119, 0.119, 0.102, 0.085, 0.094, 0.083, 0.033,
                 0.040, 0.043, 0.020, 0.020, 0.014, 0.014, 0.014, 0.014, 0.014,
                 0.014, 0.015, 0.079, 0.116, 0.127, 0.109])
sd_y = np.array([0.019, 0.019, 0.020, 0.019, 0.019, 0.019, 0.020, 0.018, 0.017,
                 0.016, 0.016, 0.015, 0.016, 0.015, 0.015, 0.015, 0.015, 0.016,
                 0.016, 0.015, 0.017, 0.019, 0.019, 0.019])
sd_z = np.array([0.105, 0.097, 0.110, 0.112, 0.100, 0.082, 0.097, 0.084, 0.038,
                 0.055, 0.054, 0.032, 0.031, 0.028, 0.024, 0.024, 0.024, 0.024,
                 0.027, 0.032, 0.085, 0.105, 0.110, 0.104])

data = [timestamps, sd_x, sd_y, sd_z]

# Plotting:
class SDPlotPyQt:
    def __init__(self):
        self.app = pg.mkQApp("SD Plot")
        self.view = pg.GraphicsView()
        self.gl = pg.GraphicsLayout()

        self.view.setCentralItem(self.gl)

        # Default settings
        self.gl.setBorder(100, 100, 100)
        self.view.setWindowTitle("SD Plot")

        self.view.resize(1280, 720)

        self.gl.layout.setColumnPreferredWidth(0, 25)
        self.gl.layout.setColumnPreferredWidth(1, 1200)

        pg.setConfigOptions(antialias=True)  # Antialiasing for nicer plots

        self.padding = [1.1, 1.0]  # Padding factor
        self.force_identical_scale = True
        self.pen_alpha = 1

        self.pen_width = 4
        pen_colors = ((255, 0, 0),
                      (0, 255, 0),
                      (0, 75, 255))
        self.generate_pens(pen_colors)

        # Colours of the draggable lines, depending on theme
        self.label_fill = (0, 0, 0, 180)
        self.label_text_color = (200, 200, 100)

        # Background colour:
        self.view.setBackground("black")

        self.plot_title = "Plot Title"
        self.side_label = "Side label"

        self.grids = True

        self.save_plot_toggle = False

    def generate_pens(self, pc):
        self.pen_rgb = pc

        self.pen_rgba = (
            (pc[0][0], pc[0][1], pc[0][2], self.pen_alpha * 100),
            (pc[1][0], pc[1][1], pc[1][2], self.pen_alpha * 100),
            (pc[2][0], pc[2][1], pc[2][2], self.pen_alpha * 100),
        )

    def set_plot_title(self, plot_title: str):
        self.plot_title = plot_title

    def set_side_label(self, side_label: str):
        self.side_label = side_label

    def set_window_size(self, w: int, h: int):
        self.view.resize(w, h)

    def save_plot(self, plot_filename: str, width: int = None):
        if width is None:
            width = self.view.size().width()
        self.save_plot_toggle = True
        self.save_plot_filename = plot_filename
        self.save_plot_width = width

    def sdplot(self, data: list, drange=None, verbose=0):
        """
        drange must be a list[float, float]
        """

        # Ensure data is of the same length, and of the right type(s)
        for data_list in (data[1], data[2], data[3]):
            assert len(data[0]) == len(data_list)
            assert type(data_list) in (list, np.ndarray)
            assert type(data_list[0]) in (float, int, np.float64)

        self.trange = [data[0][0], data[0][-1]]
        if drange is None:
            self.drange = [min([min(data[1]), min(data[2]), min(data[3])]),
                           max([max(data[1]), max(data[2]), max(data[3])])]

        # Start constructing plot ============================================

        # Title and side label
        self.gl.addLabel(self.plot_title, col=1, colspan=2)
        self.gl.nextRow()
        self.gl.addLabel(self.side_label, angle=-90, rowspan=3)

        # Data plots
        plotlabels = ("X", "Y", "Z")

        # Main data plot
        p = self.gl.addPlot(xmin=self.trange[0], xmax=self.trange[1],
                            axisItems={'bottom': pg.DateAxisItem()})
        # p.setTitle(self.plot_title)
        p.setYRange(self.drange[0], self.drange[1])
        p.showGrid(x=self.grids, y=self.grids)

        for i, data_array in enumerate((data[1], data[2], data[3])):
            p.plot(
                title=plotlabels[i],
                x=data[0],
                y=np.array(data_array),
                pen=pg.mkPen(color=self.pen_rgba[i], width=self.pen_width),
                symbolBrush=pg.mkBrush(color=self.pen_rgba[i]),
                xmin=self.trange[0],
                xmax=self.trange[1],
            )

        self.view.show()

        # Save the graph layout, if necessary
        if self.save_plot_toggle:
            exporter = pg.exporters.ImageExporter(self.gl.scene())
            exporter.parameters()['width'] = self.save_plot_width
            exporter.export(self.save_plot_filename)
            if verbose > 0:
                print("Exported image as: {}".format(self.save_plot_filename))

        # Executing PyQtGraph
        pg.exec()


plot_title = """
<h3><b> Standard deviation of hourly batches of data for each axis.
</b> <br> Data from {} to {} </h3>
""".format(
    datetime.utcfromtimestamp(data[0][0]+7200).strftime("%d/%m/%y %H:%M:%S"),
    datetime.utcfromtimestamp(data[0][-1]+7200).strftime("%d/%m/%y %H:%M:%S")
    # For the record, I am not proud of the +7200, but it works as a one-off.
)
side_label = "<h3>Standard deviation (uT)<h3>"

plot_filename = "./figures/diurnal_standard_deviation.png"

plot_object = SDPlotPyQt()
plot_object.set_plot_title(plot_title)
plot_object.set_side_label(side_label)
plot_object.set_window_size(1280, 720)

plot_object.save_plot(plot_filename)

# Generate plot
plot_object.sdplot(data)