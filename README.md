# Data Processing Library
## Load DATA

This a small library developed to help HEP that work with ROOT files and want to charge this files in a python environment where is more friendly use libraries like PanDas, NumPy, etc.

- `DataProcessing()` function read the ROOT files, convert them into RootDataFrames and is capable to implement suitbale cuts previously known by the user (to reduce the memory used), transform the RDataFrames in PanDas DataFrames selecting or not the leafs desired and process this step using threads (parellising the process).