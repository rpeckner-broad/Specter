# Specter

Specter is an algorithm for the targeted analysis of data-independent acquisition mass spectrometry proteomics experiments. It can analyze data from any instrument type and window acquisition scheme. The required user inputs are a DIA data file in centroided mzML format, a spectral library in blib format, and a mass accuracy parameter, specified in parts-per-million. 

The raw output of Specter is a .csv file describing the total ion intensities (= sum of fragment ion intensities) of each precursor in the spectral library at each retention time point of the experiment. A snippet of the typical raw output file looks like this:
```
Scan index    Retention time (s)    Precursor sequence     Precursor charge    Total ion intensity
  10032          268.3763              ETLDASLPSDYLK               2                  1,569,034
  10032          268.3763             NPAADAGSNNASKK               2                  3,112,580
  10033          268.4273                 IVLVDDSIVR               2                    722,175
```
with the identifications and quantifications based on these scan-by-scan total ion intensities reported separately:
```
Precursor sequence    Precursor charge      Quant
   ETLDASLPSDYLK              2          148,110,338
  NPAADAGSNNASKK              2           32,234,856
      IVLVDDSIVR              2           11,768,772
```
Specter currently requires the cluster-computing framework Apache Spark; a cloud framework and accompanying website will appear in the future. See ```SpecterUserGuide.pdf``` for detailed instructions on how to set up and use Specter. To learn more about the algorithm, you can read about Specter in *Nature Methods*: dx.doi.org/10.1038/nmeth.4643.


