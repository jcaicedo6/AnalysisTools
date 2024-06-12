# Data Processing Library

This a small library developed to help HEP that work with ROOT files and want to charge this files in a python environment where is more friendly use libraries like PanDas, NumPy, etc.

## Load DATA
- `getData()` function read the ROOT files, convert them into RootDataFrames and is capable to implement suitbale cuts previously known by the user (to reduce the memory used), transform the RDataFrames in PanDas DataFrames selecting or not the leafs desired and process this step using threads (parellising the process), and finally create the Pandas Data Frames.

## Create ROOT histograms

- `create_hist()` is a function that creates a ROOT histograms from Pandas DF

## Implement cuts found in a optimization process

- `BestCutsPandasDF()` filters the given Pandas DataFrame based on specified cuts and return a Pandas DataFrame