"""
Notes on processing:
 - The main dataset is 9 GB, which will be assumed too large to be handled by
    crappy Python code on most machines. So the dataset will have to be cut up.
 - First split up the data into hour-long chunks
 - Downsample the chunks from 2500 Hz to 100 Hz
 - Merge the downsampled files into one file.
 - The merged downsampled file should be kept, the rest can be deleted
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
filename_base = "./data/Weekday2500Hz_LightsOff_2023-07-06_19.12.02.dat"

# Split file
time1 = time()
dp = DataProcessor(filename_base, 2500)
n_segments = dp.segment_split(verbose=verbose)
print("Elapsed time:", round(1000*(time()-time1), 3), 'ms')


# Generate filenames for further processing
print("Number of files:", n_segments)
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

# Manually delete any files whose filenames include "_sXXofXX".