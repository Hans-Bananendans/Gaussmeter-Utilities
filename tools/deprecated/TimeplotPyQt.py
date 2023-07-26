import pyqtgraph as pg
import pyqtgraph.exporters
import numpy as np

class TimeplotPyQt:
    def __init__(self):
        self.app = pg.mkQApp("Timeplot")
        self.view = pg.GraphicsView()
        self.gl = pg.GraphicsLayout()

        self.view.setCentralItem(self.gl)

        # Default settings
        self.gl.setBorder(100, 100, 100)
        self.view.setWindowTitle("Timeplot")

        self.view.resize(1280, 720)

        self.gl.layout.setColumnPreferredWidth(0, 25)
        self.gl.layout.setColumnPreferredWidth(1, 1200)

        pg.setConfigOptions(antialias=True)  # Antialiasing for nicer plots

        self.set_defaults()

    def set_defaults(self):
        self.light_theme = False

        self.padding = [1.1, 1.0]  # Padding factor
        self.force_identical_scale = True
        self.pen_alpha = 0.8

        pen_colors = ((255,   0,   0),
                      (  0, 255,   0),
                      (  0,  75, 255))
        self.generate_pens(pen_colors)

        # Colours of the draggable lines, depending on theme
        if self.light_theme:
            self.label_fill = (255, 255, 255, 180)  # White label backgrounds
            self.label_text_color = (100, 100, 100)  # Darker label text
        else:
            self.label_fill = (0, 0, 0, 180)
            self.label_text_color = (200, 200, 100)

        # Background colour:
        if self.light_theme:
            self.view.setBackground("white")  # White background
        else:
            self.view.setBackground("black")

        self.plot_title = "Plot Title"
        self.side_label = "Side label"

        self.grids = True

        self.save_plot_toggle = False

    def set_light_theme(self):
        self.light_theme = True

    def set_dark_theme(self):
        self.light_theme = False

    def generate_pens(self, pc):
        self.pen_rgb = pc

        self.pen_rgba = (
            (pc[0][0], pc[0][1], pc[0][2], self.pen_alpha*100),
            (pc[1][0], pc[1][1], pc[1][2], self.pen_alpha*100),
            (pc[2][0], pc[2][1], pc[2][2], self.pen_alpha*100),
        )

    def set_window_size(self, w: int, h: int):
        self.view.resize(w, h)

    def set_window_title(self, window_title: str):
        self.view.setWindowTitle(window_title)

    def set_padding(self, padding_h: float, padding_v: float):
        self.padding = [padding_h, padding_v]

    def set_force_identical_scale(self, force_identical_scale: bool):
        self.force_identical_scale = force_identical_scale

    def set_pen_alpha(self, alpha: float):
        """
        Set alpha of the pen.
        Important: set in a range between 0.0 - 1.0 !
        """
        self.pen_alpha = alpha
        self.generate_pens(self.pen_rgb)

    def set_pen_colors(self, pen_colors):
        """
        Set colors of the pen.
        Must be supplied as a list of listed rgb values, e.g.:
        pen_colors = ((255,   0,   0),
                      (  0, 255,   0),
                      (  0,   0, 255))
        """
        self.generate_pens(pen_colors)

    def set_label_fill(self, rgba):
        """
        Set a label fill RGBA, for example: (255, 255, 255, 180)
        """
        self.label_fill = rgba

    def set_label_text_color(self, rgb):
        """
        Set a label text RGB, for example: (255, 255, 255)
        """
        self.label_text_color = rgb

    def set_background(self, rgb):
        self.view.setBackground(rgb)

    def set_plot_title(self, plot_title: str):
        self.plot_title = plot_title

    def set_side_label(self, side_label: str):
        self.side_label = side_label

    def show_grids(self):
        self.grids = True

    def hide_grids(self):
        self.grids = False

    def save_plot(self, plot_filename: str, width: int = None):
        if width is None:
            width = self.view.size().width()
        self.save_plot_toggle = True
        self.save_plot_filename = plot_filename
        self.save_plot_width = width

    def timeplot_pyqtgraph(self, data: list, drange = None):  # noqa
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
                pen=pg.mkPen(color=self.pen_rgba[i]),
                xmin=self.trange[0],
                xmax=self.trange[1],
                # axisItems={'bottom': pg.DateAxisItem()}
            )

        self.view.show()

        # Save the graph layout, if necessary
        if self.save_plot_toggle:
            exporter = pg.exporters.ImageExporter(self.gl.scene())
            exporter.parameters()['width'] = self.save_plot_width
            exporter.export(self.save_plot_filename)

        # Executing PyQtGraph
        pg.exec()

    def timeplot_3axis_pyqtgraph(self, data: list):

        # Ensure data is of the same length, and of the right type(s)
        for data_list in (data[1], data[2], data[3]):
            assert len(data[0]) == len(data_list)
            assert type(data_list) in (list, np.ndarray)
            assert type(data_list[0]) in (float, int, np.float64)

        self.trange = [data[0][0], data[0][-1]]
        self.xrange = [min(data[1]), max(data[1])]
        self.yrange = [min(data[2]), max(data[2])]
        self.zrange = [min(data[3]), max(data[3])]

        # Force identical scale if set to True:
        if self.force_identical_scale is True:
            xmid = self.xrange[0] + (self.xrange[1] - self.xrange[0]) / 2
            ymid = self.yrange[0] + (self.yrange[1] - self.yrange[0]) / 2
            zmid = self.zrange[0] + (self.zrange[1] - self.zrange[0]) / 2

            xscale = self.xrange[1] - self.xrange[0]
            yscale = self.yrange[1] - self.yrange[0]
            zscale = self.zrange[1] - self.zrange[0]

            scale = max([xscale, yscale, zscale])

            self.xrange = [xmid - scale / 2, xmid + scale / 2]
            self.yrange = [ymid - scale / 2, ymid + scale / 2]
            self.zrange = [zmid - scale / 2, zmid + scale / 2]

        # Start constructing plot ============================================

        # Title and side label
        self.gl.addLabel(self.plot_title, col=1, colspan=2)
        self.gl.nextRow()
        self.gl.addLabel(self.side_label, angle=-90, rowspan=3)

        # Data plots
        plotlabels = ("X", "Y", "Z")
        p_range = (self.xrange, self.yrange, self.zrange)
        for i, data_array in enumerate((data[1], data[2], data[3])):

            # Main data plot
            p = self.gl.addPlot(
                title=plotlabels[i],
                x=data[0],
                y=np.array(data_array),
                pen=pg.mkPen(color=self.pen_rgba[i]),
                xmin=self.trange[0],
                xmax=self.trange[1],
                axisItems={'bottom': pg.DateAxisItem()}
            )
            p.setYRange(p_range[i][0], p_range[i][1])
            p.showGrid(x=self.grids, y=self.grids)

            # Horizontal nfinite line object (min)
            infline_min = pg.InfiniteLine(
                angle=0,
                label="Min: {value:.2f} uT",
                pen=pg.mkPen(color=self.pen_rgb[i]),
                pos=min(data[i+1]),
                movable=False,
                labelOpts={"position": 0.05,
                           "color": self.pen_rgb[i],
                           "fill": self.label_fill}
            )
            p.addItem(infline_min)

            # Horizontal nfinite line object (max)
            infline_max = pg.InfiniteLine(
                angle=0,
                label="Max: {value:.2f} uT",
                pen=pg.mkPen(color=self.pen_rgb[i]),
                pos=max(data[i+1]),
                movable=False,
                labelOpts={"position": 0.05,
                           "color": self.pen_rgb[i],
                           "fill": self.label_fill}
            )
            p.addItem(infline_max)

            # Horizontal scanline
            infline_h = pg.InfiniteLine(
                angle=90,
                label="{value:.0f}",
                pos=self.trange[1],
                movable=True,
                bounds=[self.trange[0], self.trange[1]],
                labelOpts={"position": 0.12,
                           "color": self.label_text_color,
                           "fill": self.label_fill}
            )
            p.addItem(infline_h)

            # Vertical scanline
            infline_v = pg.InfiniteLine(
                angle=0,
                label="{value:.02f}",
                pos=p_range[i][1],
                movable=True,
                bounds=[p_range[i][0], p_range[i][1]],
                labelOpts={"position": 0.88,
                           "color": self.label_text_color,
                           "fill": self.label_fill}
            )
            p.addItem(infline_v)

            self.gl.nextRow()

        self.view.show()

        # Save the graph layout, if necessary
        if self.save_plot_toggle:
            exporter = pg.exporters.ImageExporter(self.gl.scene())
            exporter.parameters()['width'] = self.save_plot_width
            exporter.export(self.save_plot_filename)

        # Executing PyQtGraph
        pg.exec()


