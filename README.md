# Gaussmeter-Utilities
Data acquisition and processing utilities for an AlphaLab analog milligauss meter. The three analog ports are fed into a NI USB-6008 data logger and read out via USB connection over PC. To read out, the Python API for NIDAQMX is used.


## What is this?
This repository contains a number of software utilities used to collect analog magnetometer data, to process this data, and to generate useful visualizations.

This repository is primarily a place to store working files for a one-time analysis. The files are left in their original state, so that the work may be reproduced at a later date. However, the tools contained within may be of more general use.

### Examples
 `VectorPlot`
![alt text](./figures/disturbance_visualization.png?raw=true)

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

## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-sa/3.0/"><img alt="Creative Commons License" style="border-width:0" src="https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png" /></a><br />
The contents of this repository are licensed under [CC BY-NC-SA 3.0](https://creativecommons.org/licenses/by-nc-sa/3.0/), except where indicated otherwise.
