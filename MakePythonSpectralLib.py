import sqlite3
import struct
import zlib
import csv
import os
import sys
import pickle
import numpy as np
import pandas as pd
 

args = sys.argv
inputPath = args[1]

libPath = os.path.expanduser(inputPath)
    if os.path.exists(libPath):
        Lib = sqlite3.connect(libPath)
        LibPrecursorInfo = pd.read_sql("SELECT * FROM RefSpectra",Lib)
    
        SpectraLibrary = {}
    
        for i in range(len(LibPrecursorInfo)):
            precursorID = str(LibPrecursorInfo['id'][i])
            precursorKey = (LibPrecursorInfo['peptideModSeq'][i],LibPrecursorInfo['precursorCharge'][i]) 
            NumPeaks = pd.read_sql("SELECT numPeaks FROM RefSpectra WHERE id = "+precursorID,Lib)['numPeaks'][0]
            
            SpectrumMZ = pd.read_sql("SELECT peakMZ FROM RefSpectraPeaks WHERE RefSpectraID = " + precursorID,Lib)['peakMZ'][0]
            SpectrumIntensities = pd.read_sql("SELECT peakIntensity FROM RefSpectraPeaks WHERE RefSpectraID = "+precursorID,Lib)['peakIntensity'][0]
            
            if len(SpectrumMZ) == 8*NumPeaks and len(SpectrumIntensities) == 4*NumPeaks:
                SpectraLibrary.setdefault(precursorKey,{})
                SpectrumMZ = struct.unpack('d'*NumPeaks,SpectrumMZ)
                SpectrumIntensities = struct.unpack('f'*NumPeaks,SpectrumIntensities)
                SpectraLibrary[precursorKey]['Spectrum'] = np.array((SpectrumMZ,SpectrumIntensities)).T
                SpectraLibrary[precursorKey]['PrecursorMZ'] = LibPrecursorInfo['precursorMZ'][i]
                SpectraLibrary[precursorKey]['PrecursorRT'] = LibPrecursorInfo['retentionTime'][i]      #The library retention time is given in minutes
            elif len(SpectrumIntensities) == 4*NumPeaks:
                SpectraLibrary.setdefault(precursorKey,{})
                SpectrumMZ = struct.unpack('d'*NumPeaks,zlib.decompress(SpectrumMZ))
                SpectrumIntensities = struct.unpack('f'*NumPeaks,SpectrumIntensities)
                SpectraLibrary[precursorKey]['Spectrum'] = np.array((SpectrumMZ,SpectrumIntensities)).T
                SpectraLibrary[precursorKey]['PrecursorMZ'] = LibPrecursorInfo['precursorMZ'][i]
                SpectraLibrary[precursorKey]['PrecursorRT'] = LibPrecursorInfo['retentionTime'][i]
            elif len(SpectrumMZ) == 8*NumPeaks:
                SpectraLibrary.setdefault(precursorKey,{})
                SpectrumMZ = struct.unpack('d'*NumPeaks,SpectrumMZ)
                SpectrumIntensities = struct.unpack('f'*NumPeaks,zlib.decompress(SpectrumIntensities))
                SpectraLibrary[precursorKey]['Spectrum'] = np.array((SpectrumMZ,SpectrumIntensities)).T
                SpectraLibrary[precursorKey]['PrecursorMZ'] = LibPrecursorInfo['precursorMZ'][i]
                SpectraLibrary[precursorKey]['PrecursorRT'] = LibPrecursorInfo['retentionTime'][i]
            elif len(zlib.decompress(SpectrumMZ)) == 8*NumPeaks and len(zlib.decompress(SpectrumIntensities)) == 4*NumPeaks:
                SpectraLibrary.setdefault(precursorKey,{})
                SpectrumMZ = struct.unpack('d'*NumPeaks,zlib.decompress(SpectrumMZ))
                SpectrumIntensities = struct.unpack('f'*NumPeaks,zlib.decompress(SpectrumIntensities))
                SpectraLibrary[precursorKey]['Spectrum'] = np.array((SpectrumMZ,SpectrumIntensities)).T
                SpectraLibrary[precursorKey]['PrecursorMZ'] = LibPrecursorInfo['precursorMZ'][i]
                SpectraLibrary[precursorKey]['PrecursorRT'] = LibPrecursorInfo['retentionTime'][i]

outputPath = inputPath+"_PythonLibrary"
pickle.dump(SpectraLibrary,open(outputPath,"wb"))
