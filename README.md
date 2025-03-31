# Data Processing Library

This a small library developed to help HEP analyzers that work with ROOT files and want to charge this files in a python environment where is more friendly use libraries like PanDas, NumPy, etc.

## `DataProcessing` module
In this module we can find various useful functions that help us to facilite the work with big datasets.
### Load DATA
- `getData()` function read the ROOT files, convert them into RootDataFrames and is capable to implement suitable cuts previously known by the user (to reduce the memory used), transform the RDataFrames in PanDas DataFrames selecting or not the leafs desired. This process uses threads (parellising the process), and finally create the Pandas Data Frames.

### Create ROOT histograms

- `create_hist()` is a function that creates ROOT histograms from Pandas DF

### Implement cuts found in a optimization process

- `BestCutsPandasDF()` filters the given Pandas DataFrame based on specified cuts and return a Pandas DataFrame

### Saving files in csv format
- A amicable format after apply all the necessary cuts is the *csv* format. The `save_to_csv()` is used to create those kind of files.

## `kin_variables` module
`kin_variables` is an especific module useful for particular purposes in the study of decays that involve undetectable products in decays produced in leptonic colliders. Here we can find multiple kinematical variables such as ARGUS normalized energy, $q^2$, $M^2_{min}$ and $M^2_{max}$.
- ARGUS normalized energy: This variable was developed by ARGUS Collaboration. It is a well-known variable used to study tau decays in the pseudo-rest-frame.
- $q^2$: Analisis based on $B^+ \to K^+ \nu\nu$. Formula found in next link, [page-28](https://docs.belle2.org/record/3785/files/BELLE2-TALK-DRAFT-2023-117.pdf)
- $M^2_{min}$ and $M^2_{max}$: These variables could be found in the [paper](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.102.115001)