# Specter

## Description

Specter is an algorithm for the analysis of data-independent acquisition mass spectrometry experiments based on linear deconvolution. It can analyze data from any instrument type and window acquisition scheme. The required user inputs are

1. A DIA data file in mzML format.
2. A spectral library in blib format.
3. A mass accuracy parameter, specified in parts-per-million. 

The output of Specter is a .csv file describing the total ion intensities (= sum of fragment ion intensities) of each precursor in the spectral library at each retention time point of the experiment. A typical row of this output file looks like this:
```
Scan index    Retention time (s)    Precursor sequence     Precursor charge    Total ion intensity
  10032          268.3763              ETLDASLPSDYLK               2                1,569,034
```
Specter currently requires the cluster-computing framework Apache Spark; a cloud framework and accompanying website will appear in the future. 

## Running a Specter job

The general syntax for running an analysis job with Specter is (from the command line of the Spark-enabled cluster):
```bash
spark-submit --driver-memory <mem> Specter_Spark.py <mzMLname> <blibName> <index> <StartOrEnd> <numPartitions> <instrumentType> <tol>
```
where the bracketed arguments are as follows:

* mem: The amount of memory to be provisioned to the Spark driver node. 
* mzMLname: the full path to the mzML file containing the DIA data to be analyzed, without the mzML extension.
* blibName: the name of the blib file containing the spectral library, without the blib extension. 
* index: the first or last index of the subset of MS2 spectra to be analyzed 
* StartOrEnd: Should index be interpreted as the first (StartOrEnd = "start") or last (StartOrEnd = "end") index of the spectra to be analyzed? This is useful for breaking jobs into smaller pieces to respect cluster memory constraints. 
* numPartitions: The number of partitions Spark will use to parallelize the MS2 spectra. A reasonable starting choice is five times the number of cluster CPUs. 
* instrumentType: This can be one of 'orbitrap','tof', or 'other'. Use of this argument in the first two cases helps avoid certain known issues with mzMLs coming from data converted from these instrument types. 
* tol: The instrument mass tolerance, in parts-per-million.

For example, the command
```bash
spark-submit --driver-memory 30g Specter_Spark.py /rpeckner/data/20170501_PhosphoDIARun1 /rpeckner/libs/HumanPhosphoLib 100000 end 200 orbitrap 1e-5
```
would tell Specter to analyze the first 100,000 MS2 spectra in the Orbitrap DIA experiment file /rpeckner/data/20170501_PhosphoDIARun1.mzML using the spectral library /rpeckner/libs/HumanPhosphoLib.blib with a mass tolerance parameter of 10 p.p.m. and 200 partitions of the associated Spark RDD (the elements of this RDD correspond to the individual MS2 spectra).

## System requirements

Specter requires the use of Apache Spark (with the PySpark API) and Python >= 2.7.9 (Specter hasn't been tested with Python 3). The python packages pymzml and cvxopt are also required; depending on the administrative permissions for your cluster, this may necessitate use of an Anaconda environment. The Specter job commands above must be run from the directory containing the scripts Specter_Spark.py and sparse_nnls.py from this repository. 


