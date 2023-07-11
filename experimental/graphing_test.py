#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filename.py

@author: Johan Monster
"""


# From Martin Fitzpatrick's book:
# https://www.pythonguis.com/tutorials/plotting-pyqtgraph/

import sys
from PySide6 import QtWidgets
import pyqtgraph as pg  # import PyQtGraph after Qt


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.graphWidget = pg.PlotWidget()
        self.setCentralWidget(self.graphWidget)
        hour = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        temperature_1 = [30, 32, 34, 32, 33, 31, 29, 32, 35, 45]
        temperature_2 = [50, 35, 44, 22, 38, 32, 27, 38, 32, 44]
        # plot data: x, y values
        pen_1 = pg.mkPen(color=(0, 0, 255))
        pen_2 = pg.mkPen(color=(255, 0, 0))
        line_1 = self.graphWidget.plot(hour, temperature_1, pen=pen_1, name="Dataset1")
        line_2 = self.graphWidget.plot(hour, temperature_2, pen=pen_2, name="Dataset2")

        # self.graphWidget.setBackground("w")  # White background
        self.graphWidget.setTitle("Title", size="24pt")
        self.graphWidget.setLabel("left", "Temperature (°C)", size="18pt")
        self.graphWidget.setLabel("bottom", "xlabel (s)", size="18pt")
        self.graphWidget.addLegend()
        self.graphWidget.showGrid(x=True, y=True)
        # self.graphWidget.setXRange(5, 20, padding=0)
        # self.graphWidget.setYRange(30, 40, padding=0)

        # # To refresh plot:
        # self.graphWidget.clear()

        # # Updating exiting plots by updating data alone:
        # temperature_1_new = [30, 32, 34, 32, 33, 31, 29, 32, 35, 45]
        # line_1.setData(hour, temperature_1_data)


app = QtWidgets.QApplication(sys.argv)
main = MainWindow()
main.show()
app.exec()

