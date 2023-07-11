# -*- coding: utf-8 -*-

# Source:
# 

# from PyDAQmx import *
from time import time
import numpy as np
import PyDAQmx
import ctypes

# Declaration of variable passed by reference
taskHandle = PyDAQmx.TaskHandle()
read = PyDAQmx.DAQmxTypes.int32()
data = np.zeros((1000,), dtype=np.float64)

stopwatch = time()

try:
    # DAQmx Configure Code
    PyDAQmx.DAQmxCreateTask("",ctypes.byref(taskHandle))
    PyDAQmx.DAQmxCreateAIVoltageChan(
        taskHandle,
        "Dev1/ai0",
        "",
        PyDAQmx.DAQmx_Val_Cfg_Default,
        -10.0,
        10.0,
        PyDAQmx.DAQmx_Val_Volts,
        None
        )
    
    PyDAQmx.DAQmxCfgSampClkTiming(
        taskHandle,
        "",
        100.0,
        PyDAQmx.DAQmx_Val_Rising,
        PyDAQmx.DAQmx_Val_FiniteSamps,
        1000
        )

    # DAQmx Start Code
    PyDAQmx.DAQmxStartTask(taskHandle)

    # DAQmx Read Code
    PyDAQmx.DAQmxReadAnalogF64(
        taskHandle,
        1000,
        10.0,
        PyDAQmx.DAQmx_Val_GroupByChannel,
        data,
        1000,
        ctypes.byref(read),
        None
        )

    print(f"Acquired {read.value} points")
except PyDAQmx.DAQError as err:
    print("DAQmx Error:",err)
finally:
    if taskHandle:
        # DAQmx Stop Code
        PyDAQmx.DAQmxStopTask(taskHandle)
        PyDAQmx.DAQmxClearTask(taskHandle)
        
print(1000*(time()-stopwatch))