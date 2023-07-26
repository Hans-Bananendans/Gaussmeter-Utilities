"""
Notes on processing:
 - To generate control data, a 1 hour segment of WeekdayTest1 is used.
 - The segment chosen is the one that corresponds to the time interval
    13:00-14:00, as this is closest to the test interval of the data in
    RackElimination_2023-07-04_12.51.32.dat (albeit on a different day)
 - Sample rate is already appropriate at 100 Hz, no need to downsample.

Processing chain:
 - Perform segment split
 - Choose segment that corresponds to 13:00 in the day.
 - MANUALLY remove other files
"""

# Imports ====================================================================
from tools.large_file_processing import DataProcessor

filename = "./data/WeekdayTest_2023-07-03_10.41.03.dat"
sample_rate = 100

fp = DataProcessor(filename, sample_rate)
fp.segment_split(verbose=2)

# To do manually:
# - Rename segment 4 (corresponding to 13:00-14:00) to:
#       RackEliminationControl_2023-07-03_13.00.00.dat
# - Delete all other generated segments, as they are not needed