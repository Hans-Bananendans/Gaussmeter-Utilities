"""
Notes on processing:
 - No spectral analysis will be done on the data, so the data will be
    downsampled from 2500 S/s to 100 S/s.
"""

# Imports ====================================================================
from tools.GaussmeterLib import data_downsample

filenames = ("./data/TableElimination_2023-07-14_13.12.26.dat",
             "./data/TableEliminationControl_2023-07-14_12.19.22.dat"
             )

downsampling_factor = 25

data_downsample(filenames, downsampling_factor, verbose=2)
