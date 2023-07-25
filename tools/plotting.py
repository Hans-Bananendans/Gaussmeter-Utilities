import pyqtgraph as pg
import pyqtgraph.exporters
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d

from .cp import (
    Vertex,
    Face,
    Geometry,
    Frame,
    plot_face,
    plot_arrow,
    plot_global_tripod,
    plot_frame
)

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


class SpectralplotPyQt:
    def __init__(self):
        self.app = pg.mkQApp("Spectralplot")
        self.view = pg.GraphicsView()
        self.gl = pg.GraphicsLayout()

        self.view.setCentralItem(self.gl)

        # Default settings
        self.gl.setBorder(100, 100, 100)
        self.view.setWindowTitle("Spectralplot")

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

        self.plot_title = "PLOT TITLE"
        self.xlabel = "XLABEL"
        self.ylabel = "YLABEL"

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

    def set_xlabel(self, xlabel: str):
        self.xlabel = xlabel

    def set_ylabel(self, ylabel: str):
        self.ylabel = ylabel

    # def set_bottom_label(self, bottom_label: str):
    #     self.bottom_label = bottom_label

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

    def spectralplot_pyqtgraph(self, data: list):

        # Ensure data is of the same length, and of the right type(s)
        # for data_list in (data[1], data[2], data[3]):
        #     assert len(data[0]) == len(data_list)
        #     assert type(data_list) in (list, np.ndarray)
        #     assert type(data_list[0]) in (float, int, np.float64)

        self.trange = [data[0][0], data[0][-1]]
        self.drange = [min([min(data[1]), min(data[2]), min(data[3])]),
                       max([max(data[1]), max(data[2]), max(data[3])])]

        # self.xrange = [min(data[1]), max(data[1])]
        # self.yrange = [min(data[2]), max(data[2])]
        # self.zrange = [min(data[3]), max(data[3])]
        #
        # # Force identical scale if set to True:
        # if self.force_identical_scale is True:
        #     xmid = self.xrange[0] + (self.xrange[1] - self.xrange[0]) / 2
        #     ymid = self.yrange[0] + (self.yrange[1] - self.yrange[0]) / 2
        #     zmid = self.zrange[0] + (self.zrange[1] - self.zrange[0]) / 2
        #
        #     xscale = self.xrange[1] - self.xrange[0]
        #     yscale = self.yrange[1] - self.yrange[0]
        #     zscale = self.zrange[1] - self.zrange[0]
        #
        #     scale = max([xscale, yscale, zscale])
        #
        #     self.xrange = [xmid - scale / 2, xmid + scale / 2]
        #     self.yrange = [ymid - scale / 2, ymid + scale / 2]
        #     self.zrange = [zmid - scale / 2, zmid + scale / 2]

        # Start constructing plot ============================================

        # p2 = self.view.addPlot(title="Multiple curves")
        # p2.plot(np.random.normal(size=100), pen=(255, 0, 0), name="Red curve")
        # p2.plot(np.random.normal(size=110) + 5, pen=(0, 255, 0), name="Green curve")
        # p2.plot(np.random.normal(size=120) + 10, pen=(0, 0, 255), name="Blue curve")
        #
        # Title and side label
        # self.gl.addLabel(self.plot_title, col=1, colspan=2)
        # self.gl.nextRow()
        # self.gl.addLabel(self.side_label, angle=-90, rowspan=1)

        # Data plots
        plotlabels = ("X", "Y", "Z")
        # p_range = (self.xrange, self.yrange, self.zrange)

        # Main data plot
        p = self.gl.addPlot(xmin=self.trange[0], xmax=self.trange[1], colspan=2)
        p.setTitle(self.plot_title)
        p.setLabel('bottom', self.xlabel)
        p.setLabel('left', self.ylabel)
        p.setYRange(self.drange[0], self.drange[1])
        p.showGrid(x=self.grids, y=self.grids)
        p.setLogMode(x=True, y=True)

        # p.plot(
        #     title=plotlabels[0],
        #     x=data[0],
        #     y=data[1],
        #     pen=pg.mkPen(color=self.pen_rgba[0]),
        #     xmin=self.trange[0],
        #     xmax=self.trange[1],
        #     # axisItems={'bottom': pg.DateAxisItem()}
        # )

        for i, data_array in enumerate((data[1], data[2], data[3])):
            p.plot(
                title=plotlabels[i],
                x=data[0],
                y=np.array(data_array),
                pen=pg.mkPen(color=self.pen_rgba[i]),
                xmin=self.trange[0],
                xmax=self.trange[1],
            )

        # self.gl.nextRow()
        # self.gl.addLabel(self.bottom_label, col=1, colspan=2)

        # print(self.pen_rgba)

        self.view.show()

        # Save the graph layout, if necessary
        if self.save_plot_toggle:
            exporter = pg.exporters.ImageExporter(self.gl.scene())
            exporter.parameters()['width'] = self.save_plot_width
            exporter.export(self.save_plot_filename)
            print("Exported file {}.".format(self.save_plot_filename))

        # Executing PyQtGraph
        pg.exec()


    def spectralplot_pyqtgraph_separateaxes(self, data: list):

        # Ensure data is of the same length, and of the right type(s)
        # for data_list in (data[1], data[2], data[3]):
        #     assert len(data[0]) == len(data_list)
        #     assert type(data_list) in (list, np.ndarray)
        #     assert type(data_list[0]) in (float, int, np.float64)

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
        self.gl.addLabel(self.ylabel, angle=-90, rowspan=3)

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
                # axisItems={'bottom': pg.DateAxisItem()}
            )
            p.setYRange(p_range[i][0], p_range[i][1])
            p.showGrid(x=self.grids, y=self.grids)
            p.setLogMode(x=False, y=True)

            self.gl.nextRow()
            self.gl.addLabel(self.xlabel, col=1, colspan=2)

        self.view.show()

        # Save the graph layout, if necessary
        if self.save_plot_toggle:
            exporter = pg.exporters.ImageExporter(self.gl.scene())
            exporter.parameters()['width'] = self.save_plot_width
            exporter.export(self.save_plot_filename)

        # Executing PyQtGraph
        pg.exec()


def create_hhc_elements(coil_sides, coil_spacings):
    # Unpacking coil sides into X/Y/Z
    [coilX_side, coilY_side, coilZ_side] = coil_sides  # [m]
    # Unpacking coil spacings into X/Y/Z
    [coilX_spacing, coilY_spacing, coilZ_spacing] = coil_spacings

    coilXp = Face(Vertex([coilX_spacing / 2, coilX_side / 2, coilX_side / 2]),
                  Vertex([coilX_spacing / 2, -coilX_side / 2, coilX_side / 2]),
                  Vertex([coilX_spacing / 2, -coilX_side / 2, -coilX_side / 2]),
                  Vertex([coilX_spacing / 2, coilX_side / 2, -coilX_side / 2]))
    coilXn = Face(Vertex([-coilX_spacing / 2, coilX_side / 2, coilX_side / 2]),
                  Vertex([-coilX_spacing / 2, -coilX_side / 2, coilX_side / 2]),
                  Vertex([-coilX_spacing / 2, -coilX_side / 2, -coilX_side / 2]),
                  Vertex([-coilX_spacing / 2, coilX_side / 2, -coilX_side / 2]))

    coilYp = Face(Vertex([coilY_side / 2, coilY_spacing / 2, coilY_side / 2]),
                  Vertex([coilY_side / 2, coilY_spacing / 2, -coilY_side / 2]),
                  Vertex([-coilY_side / 2, coilY_spacing / 2, -coilY_side / 2]),
                  Vertex([-coilY_side / 2, coilY_spacing / 2, coilY_side / 2]))
    coilYn = Face(Vertex([coilY_side / 2, -coilY_spacing / 2, coilY_side / 2]),
                  Vertex([coilY_side / 2, -coilY_spacing / 2, -coilY_side / 2]),
                  Vertex([-coilY_side / 2, -coilY_spacing / 2, -coilY_side / 2]),
                  Vertex([-coilY_side / 2, -coilY_spacing / 2, coilY_side / 2]))

    coilZp = Face(Vertex([coilZ_side / 2, coilZ_side / 2, coilZ_spacing / 2]),
                  Vertex([-coilZ_side / 2, coilZ_side / 2, coilZ_spacing / 2]),
                  Vertex([-coilZ_side / 2, -coilZ_side / 2, coilZ_spacing / 2]),
                  Vertex([coilZ_side / 2, -coilZ_side / 2, coilZ_spacing / 2]))
    coilZn = Face(Vertex([coilZ_side / 2, coilZ_side / 2, -coilZ_spacing / 2]),
                  Vertex([-coilZ_side / 2, coilZ_side / 2, -coilZ_spacing / 2]),
                  Vertex([-coilZ_side / 2, -coilZ_side / 2, -coilZ_spacing / 2]),
                  Vertex([coilZ_side / 2, -coilZ_side / 2, -coilZ_spacing / 2]))

    return [[coilXn, coilXp], [coilYn, coilYp], [coilZn, coilZp]]

def create_cuboid(lx, ly, lz):
    # Define vertices
    p1 = Vertex([ lx/2,  ly/2,  lz/2])      # .   7___ 6
    p2 = Vertex([-lx/2,  ly/2,  lz/2])      # .   |\8__\5
    p3 = Vertex([-lx/2, -ly/2,  lz/2])      # .   | |  |        Y
    p4 = Vertex([ lx/2, -ly/2,  lz/2])      # .  3| | 2| -------->
    p5 = Vertex([ lx/2,  ly/2, -lz/2])      # .    \|__|
    p6 = Vertex([-lx/2,  ly/2, -lz/2])      # .    4  \ 1
    p7 = Vertex([-lx/2, -ly/2, -lz/2])      # .        \
    p8 = Vertex([ lx/2, -ly/2, -lz/2])      # .         v X

    # Define faces
    fA = Face(p1, p4, p3, p2)  # .  E ___
    fB = Face(p3, p4, p8, p7)  # .   |\ F_\
    fC = Face(p4, p1, p5, p8)  # . B | |  | D        Y
    fD = Face(p1, p2, p6, p5)  # .   | | C|  -------->
    fE = Face(p2, p3, p7, p6)  # .    \|__|
    fF = Face(p5, p6, p7, p8)  # .      A

    # Assembling geometry
    cuboid = Geometry([fA, fB, fC, fD, fE, fF])

    return cuboid

def create_table(leg_height=1.23, leg_width=0.04, leg_spacing=0.52,
                 top_side=0.61, top_height=0.02, ox=0, oy=0, oz=-1.25):

    frame_legs = Frame()
    frame_top = Frame()

    for i in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
        leg = create_cuboid(leg_width, leg_width, leg_height)
        leg.translate(i[0]*leg_spacing/2, i[1]*leg_spacing/2, leg_height/2)
        frame_legs.add_geometry(leg)

    top = create_cuboid(top_side, top_side, top_height)
    top.translate(0, 0, leg_height+top_height/2)
    frame_top.add_geometry(top)

    for frame in (frame_legs, frame_top):
        frame.translate(ox, oy, oz)

    return frame_legs, frame_top

def create_x_wall(ly=2., lz=2., ox=0., oy=0., oz=0.):
    wall_x = Face(Vertex([ox,  ly / 2 + oy,  lz / 2 + oz]),
                  Vertex([ox, -ly / 2 + oy,  lz / 2 + oz]),
                  Vertex([ox, -ly / 2 + oy, -lz / 2 + oz]),
                  Vertex([ox,  ly / 2 + oy, -lz / 2 + oz]))
    return wall_x

def create_y_wall(lx=2., lz=2., ox=0., oy=0., oz=0.):
    wall_y = Face(Vertex([-lx / 2 + ox, oy,  lz / 2 + oz]),
                  Vertex([ lx / 2 + ox, oy,  lz / 2 + oz]),
                  Vertex([ lx / 2 + ox, oy, -lz / 2 + oz]),
                  Vertex([-lx / 2 + ox, oy, -lz / 2 + oz]))
    return wall_y

# def create_floor(ly=2, lz=2, ox=0, oy=0, oz=0):
#     pass  # TODO: Implement if needed


class VectorPlot:
    def __init__(self):
        # Define matplotlib plot structure
        self.fig = plt.figure(figsize=(15, 10.5))
        self.ax = mp3d.Axes3D(self.fig, auto_add_to_figure=False)
        self.fig.add_axes(self.ax)
        self.ax.clear()

        self.default_plotting_settings()

        self.plotscale = 1.1
        self.ax.set_xlim(-4 / 3 * self.plotscale, 4 / 3 * self.plotscale)
        self.ax.set_ylim(-4 / 3 * self.plotscale, 4 / 3 * self.plotscale)
        self.ax.set_zlim(-self.plotscale, self.plotscale)

        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_zlabel('z')

        # Default camera view
        self.ax.view_init(elev=20, azim=-50)

    def show(self):
        plt.show()

    def set_plot_title(self, plot_title: str):
        self.ax.set_title(plot_title)

    def default_plotting_settings(self):
        self.plot_settings = {
            # Tripod properties
            "show_tripod": True,  # If False, does not plot the tripod
            "tripod_scale": 1,  # Sets the scale of the tripod
            "plot_perpendiculars": True,  # Plot perpendiculars

            # Vertex plotting properties:
            "vertexfill": True,  # If False, vertex will not be plotted
            "vertexcolour": "#000",  # Specifies the vertex colour
            "vertexsize": 10,  # Size of the plotted vertex
            "vertexalpha": 1,  # Opacity of the plotted vertex

            # Face plotting properties:
            "linefill": True,  # If False, does not plot face lines
            "linecolour": "#000",  # Colour of face lines
            "linewidth": 2,  # Thickness of face lines
            "linealpha": 1,  # Opacity of face lines
            "facefill": True,  # If False, does not shade the face area
            "facecolour": "#555",  # Colour of the face area shading
            "facealpha": 1,  # Opacity of the face area shading

            # Face perpendicular arrow plotting properties:
            "perpfill": False,  # If True, plot perpendiculars
            "perpcolour": "#888",  # Specifies the perp. arrow colour
            "perpscale": 1,  # Size of the plotted perp. arrow
            "perpalpha": 0.5,  # Opacity of the plotted perp. arrow

            # Illumination:
            "illumination": False,  # If True, plots illumination intensity
            "ill_value": 0,  # Used to plot illumination intensity
            "ill_plane": None,  # If illumination is used, a plane

            # Vector plotting properties:
            "vectorfill": True,  # If False, does not plot vector arrow
            "vectorcolour": "#000",  # Colour of vector arrow
            "vectoralpha": 1,  # Opacity of vector arrow
            "vectorscale": 1,  # Scale the whole vector by a constant
            "vectorratio": 0.15  # Vector arrow length ratio
        }

    def plot_vector(self, tip, origin=(0, 0, 0), linewidth=1.5,
                    alpha=1, scaling=0.02, alr=0.1, color="purple"):
        plot_arrow(self.ax, origin, tip, linewidth=linewidth,
                   alpha=alpha, scaling=scaling, alr=alr, color=color)

    def plot_global_tripod(self, scaling: float = None):
        if scaling is None:
            scaling = self.plotscale / 2
        plot_global_tripod(self.ax, scaling=scaling)

    def plot_hhc_coils(self, coil_sides: list, coil_spacings: list,
                       coil_thickness: int = 5, coil_alpha: float = 0.5,
                       coil_colors=("#F00", "#0F0", "#00F")):

        self.coil_sides = coil_sides
        self.coil_spacings = coil_spacings

        # First make the cage elements
        cage_objects = create_hhc_elements(coil_sides, coil_spacings)

        for coilX in cage_objects[0]:
            plot_face(self.ax, coilX, linecolour=coil_colors[0],
                      linewidth=coil_thickness, linealpha=coil_alpha,
                      facefill=False, vertexfill=False)
        for coilY in cage_objects[1]:
            plot_face(self.ax, coilY, linecolour=coil_colors[1],
                      linewidth=coil_thickness, linealpha=coil_alpha,
                      facefill=False, vertexfill=False)
        for coilZ in cage_objects[2]:
            plot_face(self.ax, coilZ, linecolour=coil_colors[2],
                      linewidth=coil_thickness, linealpha=coil_alpha,
                      facefill=False, vertexfill=False)

    # Plot the lab walls
    def plot_x_wall(self, ly=2., lz=2., ox=0., oy=0., oz=0.,
                    linecolour="#333", linewidth=3, linealpha=0.4,
                    facefill=True, facealpha=0.25):
        face = create_x_wall(ly=ly, lz=lz, ox=ox, oy=oy, oz=oz)
        plot_face(self.ax, face,
                  linecolour=linecolour, linewidth=linewidth,
                  linealpha=linealpha, facefill=facefill, facealpha=facealpha,
                  vertexfill=False)

    def plot_y_wall(self, lx=2., lz=2., ox=0., oy=0., oz=0.,
                    linecolour="#333", linewidth=3, linealpha=0.4,
                    facefill=True, facealpha=0.25):
        face = create_y_wall(lx=lx, lz=lz, ox=ox, oy=oy, oz=oz)
        plot_face(self.ax, face,
                  linecolour=linecolour, linewidth=linewidth,
                  linealpha=linealpha, facefill=facefill, facealpha=facealpha,
                  vertexfill=False)

    def plot_table(self):
        frame_legs, frame_top = create_table()
        plot_frame(self.ax, frame_legs,
                   show_tripod=False,
                   vertexfill=False,
                   vertexalpha=0,
                   linecolour="#CCC",
                   linealpha=0.3,
                   facecolour="#CCC",
                   facealpha=0.2,
                   )
        plot_frame(self.ax, frame_top,
                   show_tripod=False,
                   vertexfill=False,
                   vertexalpha=0,
                   linecolour="#C95",
                   linealpha=0.3,
                   facecolour="#C95",
                   facealpha=0.2,
                   )

    def autoplot(self,
                 tripod=True,
                 coils=False,
                 walls=False,
                 table=False,
                 delfipq=False,
                 cubesat3u=False,
                 cubesat12u=False):
        if tripod:
            self.plot_global_tripod()
        if coils:
            coil_sides = [1.85, 1.95, 2.05]
            coil_spacings = [1.0073, 1.0618, 1.1162]
            self.plot_hhc_coils(coil_sides, coil_spacings)
        if walls:
            self.plot_x_wall(3, 2, -1.5, 0, -0.1)
            self.plot_y_wall(3, 2, 0, 1.5, -0.1)
        if table:
            self.plot_table()
        if delfipq:
            pass  # TODO: Implement if needed
            # geometry = create_cuboid(0.05, 0.05, 0.15)
        if cubesat3u:
            pass  # TODO: Implement if needed
            # geometry = create_cuboid(0.1, 0.1, 0.3)
        if cubesat12u:
            pass  # TODO: Implement if needed
            # geometry = create_cuboid(0.2, 0.2, 0.3)