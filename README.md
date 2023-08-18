# Gaussmeter-Utilities
Data acquisition and processing utilities for an AlphaLab analog milligauss meter. The three analog ports are fed into a NI USB-6008 data logger and read out via USB connection over PC. To read out, the Python API for NIDAQMX is used.


## What is this?
This repository contains a number of software utilities used to collect analog magnetometer data, to process this data, and to generate useful visualizations.

This repository is primarily a place to store working files for a one-time analysis. The files are left in their original state, so that the work may be reproduced at a later date. However, the tools contained within may be of more general use.

### Examples
 `VectorPlot`
![alt text](./figures/disturbance_visualization1.png?raw=true)

 `TimePlotPyQt`
![alt text](./other/high_resolution_export.png?raw=true)

 `SpectralPlotPyQt`
![alt text](./figures/Lights%20On%20-%20Spectral%20Plot.png?raw=true)

## Notable files
```markdown
Gaussmeter-Utilities
├── data/                           Location for .dat measurement files, not on Github due to large size
├── figures/                        Generated figures in .png
├── other/                          Miscellaneous files
│   ├── Test Results.ods            Overview of test results
├── tools/                          Data acquisition and processing tools
│   ├── cp.py                       Merged Cube Projector framework file (all classes)
│   ├── GaussmeterLib.py            Various utilities
│   ├── GaussmeterRecorder.py       CLI tool for data acquisition with NI-DAQmx
│   ├── large_file_processing.py    Processing tool that can cut up >RAM datasets into smaller chunks. These chunks can be automatically synchronized to whole hours using Unix timestamps.
│   ├── local_emf.py                Convenience file for manually specifying the local Earth Magnetic Field vector
│   ├── plotting.py                 Plotting classes TimePlotPyQT, SpectralPlotPyQt, VectorPlot
├── disturbance_visualization.py    Easy-to-run 3D visualization of disturbance vector
├── *.py                            Test-specific data processing file
```

## Replication instructions
Are you here to replicate parts of this analysis, or improve the tools? Here are some instructions to help you get started:
1. (recommended) Use a Unix-based operating system. If using Windows, you may have to edit some of the hardcoded file paths to use double backslashes (`\\`), instead of single forward slashes (`/`).
2. (recommended) Open a new `venv`.
3. Clone this repository.
4. Install the Python dependencies listed in `requirements.txt`.
5. If you wish to regenerate the plots, first download the magnetic test data files from the server and place them in the cloned repository, in a folder named `/data/`. All files placed in this folder are automatically ignored by git. The total size of the data files is 43 GB. If you have trouble getting access to the measurement data, send me a message.
6. Consult `/other/Test Results.ods` to see a comprehensive list of performed tests, outcome, and what files are relevant.
7. If you are missing certain files, check to see if any of the scripts in the root containing the word `preprocessing` can generate the files for you.
8. Filenames containing the term `_SPARSE` indicate that the contents of the file have been downsampled. 
9. Filenames containing term `_s**of**` indicate that the file is part of a segmentation of a larger file.
10. All tools in `/tools/` were designed to be called with Python commands, except `GaussmeterRecorder.py`, which has a CLI and is best run from the console.

<details>
<summary>List of data files</summary>
<i>
2500HzTestLightsOff_2023-07-05_11.13.36.dat <br>
2500HzTestLightsOn_2023-07-05_11.04.28.dat <br>
CoilsClose_2023-07-04_17.31.25.dat <br>
CoilsWide_2023-07-04_16.01.53.dat <br>
FireExtinguisherAway_2023-07-06_18.03.29.dat <br>
FireExtinguisherControl_2023-07-06_18.32.02.dat <br>
LightOff_2023-07-04_10.57.37.dat <br>
LightOff_2023-07-04_11.13.19.dat <br>
LightOff_2023-07-04_11.27.14.dat <br>
LightOff_2023-07-04_11.42.11.dat <br>
LightOn_2023-07-04_10.52.22.dat <br>
LightOn_2023-07-04_11.05.10.dat <br>
LightOn_2023-07-04_11.18.38.dat <br>
LightOn_2023-07-04_11.33.19.dat <br>
NightNoRack_2023-07-04_19.04.58.dat <br>
RackElimination_2023-07-04_12.51.32.dat <br>
RackEliminationControl_2023-07-03_13.00.00.dat <br>
SpatialX0Y0_2023-07-12_14.15.31_SPARSE25.dat <br>
SpatialX0Y-50_2023-07-12_08.37.24.dat <br>
SpatialX0Y-50_2023-07-12_08.37.24_SPARSE25.dat <br>
SpatialX0Y+50_2023-07-12_12.21.00.dat <br>
SpatialX0Y+50_2023-07-12_12.21.00_SPARSE25.dat <br>
SpatialX+50Y+0_2023-07-12_06.49.34.dat <br>
SpatialX+50Y+0_2023-07-12_06.49.34_SPARSE25.dat <br>
SpatialX-50Y0_2023-07-12_10.35.41.dat <br>
SpatialX-50Y0_2023-07-12_10.35.41_SPARSE25.dat <br>
SpatialX+50Y-50_2023-07-12_07.38.19.dat <br>
SpatialX+50Y-50_2023-07-12_07.38.19_SPARSE25.dat <br>
SpatialX-50Y-50_2023-07-12_09.46.26.dat <br>
SpatialX-50Y-50_2023-07-12_09.46.26_SPARSE25.dat <br>
SpatialX-50Y+50_2023-07-12_11.29.49.dat <br>
SpatialX-50Y+50_2023-07-12_11.29.49_SPARSE25.dat <br>
SpatialX+50Y+50_2023-07-12_13.18.39.dat <br>
SpatialX+50Y+50_2023-07-12_13.18.39_SPARSE25.dat <br>
TableElimination_2023-07-14_13.12.26.dat <br>
TableElimination_2023-07-14_13.12.26_SPARSE25.dat <br>
TableEliminationControl_2023-07-14_12.19.22.dat <br>
TableEliminationControl_2023-07-14_12.19.22_SPARSE25.dat <br>
TableReturnPosition1_2023-07-12_14.15.31.dat <br>
TableReturnPosition1_2023-07-12_14.15.31_SPARSE25.dat <br>
Weekday2500Hz_LightsOff_2023-07-06_19.12.02.dat <br>
Weekday2500Hz_LightsOff_2023-07-06_19.12.02_SPARSE25.dat <br>
Weekday2500Hz_LightsOn_2023-07-05_11.46.50.dat <br>
Weekday2500Hz_LightsOn_2023-07-05_11.46.50_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31.dat <br>
Weekday4_2023-07-13_08.29.31_s10of27.dat <br>
Weekday4_2023-07-13_08.29.31_s10of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s11of27.dat <br>
Weekday4_2023-07-13_08.29.31_s11of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s12of27.dat <br>
Weekday4_2023-07-13_08.29.31_s12of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s13of27.dat <br>
Weekday4_2023-07-13_08.29.31_s13of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s14of27.dat <br>
Weekday4_2023-07-13_08.29.31_s14of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s15of27.dat <br>
Weekday4_2023-07-13_08.29.31_s15of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s16of27.dat <br>
Weekday4_2023-07-13_08.29.31_s16of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s17of27.dat <br>
Weekday4_2023-07-13_08.29.31_s17of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s18of27.dat <br>
Weekday4_2023-07-13_08.29.31_s18of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s19of27.dat <br>
Weekday4_2023-07-13_08.29.31_s19of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s1of27.dat <br>
Weekday4_2023-07-13_08.29.31_s1of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s20of27.dat <br>
Weekday4_2023-07-13_08.29.31_s20of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s21of27.dat <br>
Weekday4_2023-07-13_08.29.31_s21of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s22of27.dat <br>
Weekday4_2023-07-13_08.29.31_s22of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s23of27.dat <br>
Weekday4_2023-07-13_08.29.31_s23of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s24of27.dat <br>
Weekday4_2023-07-13_08.29.31_s24of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s25of27.dat <br>
Weekday4_2023-07-13_08.29.31_s25of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s26of27.dat <br>
Weekday4_2023-07-13_08.29.31_s26of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s27of27.dat <br>
Weekday4_2023-07-13_08.29.31_s27of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s2of27.dat <br>
Weekday4_2023-07-13_08.29.31_s2of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s3of27.dat <br>
Weekday4_2023-07-13_08.29.31_s3of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s4of27.dat <br>
Weekday4_2023-07-13_08.29.31_s4of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s5of27.dat <br>
Weekday4_2023-07-13_08.29.31_s5of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s6of27.dat <br>
Weekday4_2023-07-13_08.29.31_s6of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s7of27.dat <br>
Weekday4_2023-07-13_08.29.31_s7of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s8of27.dat <br>
Weekday4_2023-07-13_08.29.31_s8of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_s9of27.dat <br>
Weekday4_2023-07-13_08.29.31_s9of27_SPARSE25.dat <br>
Weekday4_2023-07-13_08.29.31_sample_day.dat <br>
Weekday4_2023-07-13_08.29.31_sample_night.dat <br>
Weekday4_2023-07-13_08.29.31_SPARSE25.dat <br>
WeekdayTest_2023-07-03_10.41.03.dat <br>
WeekdayTest_2023-07-03_10.41.03_SPARSE50.dat <br>
WeekendTest_2023-06-30_15.07.51.dat <br>
WeekendTest_2023-06-30_15.07.51_SPARSE50.dat
</i>
</details>


## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-sa/3.0/"><img alt="Creative Commons License" style="border-width:0" src="https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png" /></a><br />
The contents of this repository are licensed under [CC BY-NC-SA 3.0](https://creativecommons.org/licenses/by-nc-sa/3.0/), except where indicated otherwise.
