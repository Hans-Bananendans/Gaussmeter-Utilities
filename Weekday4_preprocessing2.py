"""
Notes on processing:
 - The main dataset is 9.6 GB, which will be assumed too large to be handled
    by Python code on most machines. So the dataset will have to be cut up.
 - First split up the data into hour-long chunks
 - Downsample the chunks from 2500 S/s to 100 S/s
 - Merge the downsampled files into one file.
 - All of these files will be needed:
    - The raw 9.6 GB dataset is not used, but should be kept for reference
    - The hour-long chunks of 2500 S/s will be used for the spectral plots
    - The hour-long chunks downsampled to 100 S/s are used for readout
    - The merged downsampled file at 100 S/s will be used for the timeplot
"""

# Imports ====================================================================
from time import time
from tools.GaussmeterLib import (
    write_data,
    read_data
)

time0 = time()

# Base filename
filename_d = "./data/Weekday4_2023-07-13_08.29.31_s6of27.dat"
filename_n = "./data/Weekday4_2023-07-13_08.29.31_s17of27.dat"

new_filename_d = "./data/Weekday4_2023-07-13_08.29.31_sample_day.dat"
new_filename_n = "./data/Weekday4_2023-07-13_08.29.31_sample_night.dat"

data_d = read_data(filename_d, header=True)
data_n = read_data(filename_n, header=True)

len_section_d = 25000  # Take a section of ~10 seconds
len_section_n = 25000  # Take a section of ~10 seconds
start_section_d = round(25/60*len(data_d[0]))  # Start at xx:25
start_section_n = round(25/60*len(data_n[0]))  # Start at xx:25

for i in (0, 1, 2, 3):
    data_d[i] = data_d[i][start_section_d:(start_section_d + len_section_d)]
    data_n[i] = data_n[i][start_section_n:(start_section_n + len_section_n)]

write_data(new_filename_d, data_d, verbose=2)
write_data(new_filename_n, data_n, verbose=2)

print("Total elapsed time:", round(time()-time0, 3), 's')
