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
    data_merge,
    data_downsample,
    multipart_filenames,
)
from tools.large_file_processing import DataProcessor

verbose = 2
time0 = time()

# Base filename
filename_base = "./data/Weekday4_2023-07-13_08.29.31.dat"

# Split file
time1 = time()
dp = DataProcessor(filename_base, 2500)
n_segments = dp.segment_split(verbose=verbose)
print("Elapsed time:", round(1000*(time()-time1), 3), 'ms')


# Generate filenames for further processing
print("Number of segments:", n_segments)
filenames = multipart_filenames(filename_base, n_segments)


print(filenames)  # DEBUG
if input("Checkpoint reached. Press any key to continue...") is not None:
    pass

# Downsample split files
time1 = time()
data_downsample(filenames, 25, verbose=verbose)
print("Time since last checkpoint:", round(1000*(time()-time1), 3), 'ms')

# Generate filenames for further processing
for i, filename in enumerate(filenames):
    filenames[i] = filename[:-4] + "_SPARSE25" + filename[-4:]

print(filenames)  # DEBUG
if input("Checkpoint reached. Press any key to continue...") is not None:
    pass

# Merge downsampled files
time1 = time()
filename_merge = filename_base[:-4] + "_SPARSE25" + filename_base[-4:]
data_merge(filenames,
           filename_merge=filename_merge,
           len_header=10,
           verbose=verbose)
print("Elapsed time:", round(1000*(time()-time1), 3), 'ms')

print("Total elapsed time:", round(time()-time0, 3), 's')
