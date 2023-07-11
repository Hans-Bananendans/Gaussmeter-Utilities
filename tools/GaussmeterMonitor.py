# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 13:29:48 2023

@author: Johan Monster

Convenient monitoring Qt6 interface. Proof-of-concept for HHC setup GUI.

"""



#%% IMPORT DEPENDENCIES

# import os
# import sys
# import nidaqmx
# import numpy as np

# from tqdm import tqdm
# from array import array
# from time import time, sleep
# from datetime import datetime
# from argparse import ArgumentParser
# from colorama import init, Fore, Style; init(autoreset=True)

# import psutil

#%% Interface framework


# from PyQt6.QtWidgets import QApplication, QLabel, QWidget, QGridLayout
# from PyQt6.QtGui import QIcon
# from PyQt6.QtCore import Qt


from PyQt6.QtWidgets import (
    QApplication,
    QDialogButtonBox,
    QFormLayout,
    QGridLayout, 
    QMainWindow,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget, 
)
from PyQt6.QtGui import QIcon
from PyQt6.QtCore import Qt
import sys
 


# def __init__(self):

#     super().__init__(parent=None)

#     self.setWindowTitle("QMainWindow")
# class Window(QMainWindow):
#     def __init__(self):
#         super().__init__()
#         self.resize(300, 300)
#         self.setWindowTitle("CodersLegacy")
#         self.setWindowIcon(QIcon("icon.jpg"))
 
#         layout = QVBoxLayout()
#         self.setLayout(layout)
 
#         label = QLabel("Hello World")
#         label.setAlignment(Qt.AlignmentFlag.AlignCenter)
#         layout.addWidget(label)
 
class Window(QMainWindow):
    def __init__(self):
        super().__init__(parent=None)
        self.setWindowTitle("QMainWindow")
        self.setGeometry(100, 100, 800, 600)
        self.centralWidget = QWidget(self)
        self.setCentralWidget(self.centralWidget)
        self.outerLayout = QVBoxLayout()
        
        self._createMonitors()
        self._createExitButton()

    def h1(self, text):
        return "<h1>"+str(text)+"</h1>"
    
    def _createMonitors(self):

        B0 = [100.0, -200.0, 300.0]
        
        monitor_layout = QGridLayout()
        
        readout_stylesheet = "font-size: 46px;"
        
        label_Bx = QLabel(self.h1(B0[0]))
        label_By = QLabel(self.h1(B0[1]))
        label_Bz = QLabel(self.h1(B0[2]))

        
        for i, label in enumerate((label_Bx, label_By, label_Bz)):
            label.setStyleSheet(readout_stylesheet)
            monitor_layout.addWidget(label, 0, i)
        self.outerLayout.addLayout(monitor_layout)
        

    # def _createRateInput(self):
    #     formLayout = QFormLayout()
    #     formLayout.addRow("Name:", QLineEdit())
        
    #     tools = QToolBar()
    #     tools.addAction("Exit", self.close)
    #     self.addToolBar(tools)
    


    def _createExitButton(self):
        grid_layout = QGridLayout()
        button_exit = QPushButton("Exit")
        grid_layout.addWidget(button_exit, 0, 2)
        self.outerLayout.addLayout(grid_layout)

def main():
    """Main function."""
    app = QApplication([])
    app.setStyle('Fusion')
    window = Window()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()

# app = QApplication(sys.argv)
# window = Window()
# window.show()
# sys.exit(app.exec())

# app = QApplication(sys.argv)



# window = QWidget()
# window.setWindowTitle("PyQt App")
# window.setGeometry(100, 100, 800, 600)





# def h1(text):
#     return "<h1>"+str(text)+"</h1>"

# B = [100.0, -200.0, 300.0]

# layout = QGridLayout()

# readout_stylesheet = "font-size: 46px;"

# label_Bx = QLabel(h1(B[0]), parent=window)
# label_By = QLabel(h1(B[1]), parent=window)
# label_Bz = QLabel(h1(B[2]), parent=window)
# # label_Bz = layout.addWidget(QLabel(h1(B[2]), parent=window), 0, 2)

# for i, label in enumerate((label_Bx, label_By, label_Bz)):
#     label.setStyleSheet(readout_stylesheet)
#     layout.addWidget(label, 0, i)


# window.setLayout(layout)

# window.show()


# sys.exit(app.exec())

# if __name__ == "__main__": 
    # pass