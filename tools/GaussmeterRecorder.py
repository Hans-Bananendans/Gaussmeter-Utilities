# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:10:15 2023

@author: Main


Goals:
1. argparse CLI handler for issuing long term measurements (LTMs) on the 
    Gaussmeter. Input arguments should include verbosity, duration, 
    sample rate, file name prefixes, etc.

2. Progress bar interface: Function to calculate samples, calculate percentage
    points, and display a tqdm progress bar.
    
3. Hook in measurement functions ffrom nidaqmx_test3.py

4. (OPTIONAL) Insert benchmark from data_optimizations

5. Test and polish

6. Checks: Detect sample rate slowdowns!

7. (OPTIONAL) Auto speed-match corrector that experimentally tunes 'cf'

"""
NAME = "Gaussmeter Recorder"
VERSION = "1.0.0"
VERSION_DATE = "28-06-2023"


#%% IMPORT DEPENDENCIES

import sys
# import nidaqmx
# import numpy as np
# from tqdm import tqdm
# from array import array
# from time import time, sleep
# from datetime import datetime
from argparse import ArgumentParser
# from colorama import init, Fore, Style; init(autoreset=True)

from GaussmeterLib import (DataRecorder,
                           print_version)

# import psutil
                

#%% MAIN PROGRAM
if __name__ == "__main__":
    
    #%% ARGUMENT PARSING
    
    parser = ArgumentParser(
            prog="GaussmeterRecorder",
            description="CLI interface for long-term Gaussmeter recorder tool.")


    # SIMULATION ARGUMENTS
    parser.add_argument("-R", "-r", "--record", nargs=2, type=int, action="store",
                        help="Initiate a data recording session. Supply recording DURATION in SECONDS and SAMPLE RATE in S/s as arguments.")    
    parser.add_argument("-bs", "--buffer-size", nargs=1, type=int, action="store",
                        default=[-1],
                        help="Size of the buffer in whole samples.")
    parser.add_argument("-pd", "--predelay", nargs=1, type=float, action="store",
                        default=[0.0],
                        help="Predelay before simulation starts.")
    parser.add_argument("-n", "--name", nargs=1, type=str, action="store",
                        default=[""],
                        help="Name of the simulation.")
    parser.add_argument("-ft", "--filetype", nargs=1, type=str, action="store",
                        default=[".dat"],
                        help="Specify a file extension for the output data file. (default: .dat)")
    # parser.add_argument("-rp", "--readout-period", nargs=1, type=float, action="store",
    #                     default=[0.0],
    #                     help="Period in seconds to print readout values to the console. A frequency of 0 seconds means periodic readouts are OFF (default: 0)")
    # parser.add_argument("-f", "--folder", nargs=1, type=str, action="store",
    #                     default=[".\\recorded\\"],
    #                     help="Intended to specify a folder, but DOES NOT WORK CURRENTLY.")

    # DEBUGGING ARGUMENTS
    parser.add_argument("--echo", nargs=1, type=str, action="store",
                        help="Echos the argument.")
    parser.add_argument("--print_arguments", default=False, action="store_true",
                        help="Prints all arguments supplied to the CLI")
    parser.add_argument("--test-progress-bar", default=False, action="store_true",
                        help="Displays a test progress bar.")

    # UTILITY ARGUMENTS
    parser.add_argument("-q", "--quiet", default=False, action="store_true",
                        help="Suppresses all outputs. (default: false)")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="Increases verbosity. (default: 0)")
    parser.add_argument("--version", default=False, action="store_true",
                        help="Prints current version of program.")
    parser.add_argument("-y", "--skipconfirm", default=False, action="store_true",
                        help="Answers 'YES' to all y/n confirmations. (default: true)")
    
    args = parser.parse_args()
    
    #%% DECISION TREE FOR ARGUMENTS
    verbose = args.verbose
    if args.quiet:
        verbose = -1
 
        
    if args.echo: # DEBUG
        print(args.echo[0])
    if args.print_arguments: # DEBUG
        print(parser.parse_known_args())
    
    if args.record:
        recording_time = args.record[0]
        sr = args.record[1]
        n_samples = recording_time*sr
        
        # if verbose > 0:
        #     SimulationSummary(duration, sr, 
        #                           name=args.name[0],
        #                           folder=args.folder[0],
        #                           filetype=args.filetype[0], 
        #                           verbose=verbose)
        
        if not args.skipconfirm:
            carryon = input("Start simulation? (Y/N)")
            if carryon not in ["Y", "y", "yes", "Yes", "YES"]:
                print("Program cancelled.")
                sys.exit(0) # Not great but it works here.
        # if args.test_progress_bar:
            # TestProgressBar(n_samples, sr)
        
        # Create DataRecorder object
        recorder = DataRecorder(verbose=verbose)

        recorder.set_sampling_settings(
            sr=sr,
            n_channels=3,
            n_buffers=5,
            buffer_size=args.buffer_size[0],
            v_min=-1.25,
            v_max=1.25
            )
        
        recorder.set_data_settings(
            name=args.name[0],
            filetype=args.filetype[0],
            scale_v=100
            )
        
        recorder.set_recording_time(recording_time)
        
        recorder.record(predelay=args.predelay[0])
        
        # Problem: how to deterministically determine filename here?
        # data = ReadData(recording_file)
        
        # print(f"Mean = {round(1000*rate_mean,3)} ms, sd = {round(1000*rate_sd,3)} ms")
        # print(f"Target = {1000/sr}, delay factor = {round(100*(rate_mean*sr-1),1)} %")

    # Print version if requested:
    if args.version:
        print_version(VERSION, VERSION_DATE, NAME, verbose=verbose)
        