# Specter

Specter is an algorithm for the analysis of DIA mass spectra based on linear deconvolution. It can analyze data from any instrument type and window acquisition scheme. The required user inputs are

1. A DIA data file in mzML format.
2. A spectral library in blib format.
3. A mass accuracy parameter, specified in parts-per-million. 

Specter currently requires the cluster-computing framework Apache Spark; a cloud framework and accompanying website will appear in the future. 

The general syntax for running an analysis with Specter is



