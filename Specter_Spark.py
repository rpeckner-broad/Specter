import sys 
import time
import pandas as pd
import pymzml
import numpy as np
from numpy import linalg
from scipy import stats
from scipy import integrate
from scipy import sparse
from scipy.optimize import fsolve
from scipy.stats import mvn
import os
import random
import itertools
import cPickle as pickle
from functools import partial
import sparse_nnls 
from pyspark import SparkConf,SparkContext
import sqlite3
import struct
import zlib
import csv
import re

def RegressSpectraOntoLibrary(DIASpectraIterator,Library,tol,maxWindowOffset):
        
        RefSpectraLibrary = Library.value 
        
        for DIASpectrum in DIASpectraIterator:        
        
            precMZ = float(DIASpectrum[1])     
            precRT = float(DIASpectrum[2])  #MS2 scan retention time, in minutes
            index = DIASpectrum[3]
            windowWidth = DIASpectrum[4]
            
            DIASpectrum = np.array(DIASpectrum[0])
                       
            LibraryCoeffs = []
            
            if len(DIASpectrum.shape) == 2:
                
                if windowWidth > 0:
                    CandidateRefSpectraLibrary = [spectrum['Spectrum'] for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                abs(float(spectrum['PrecursorMZ']) - precMZ) < windowWidth/2]
                    MassWindowCandidates = [key for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                abs(float(spectrum['PrecursorMZ']) - precMZ) < windowWidth/2]
                else:
                    CandidateRefSpectraLibrary = [spectrum['Spectrum'] for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                        float(spectrum['PrecursorMZ']) > precMZ - maxWindowOffset/2]   
                    MassWindowCandidates = [key for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                        float(spectrum['PrecursorMZ']) > precMZ - maxWindowOffset/2]
                

                #MERGING OF POINTS IN ACQUIRED SPECTRUM WITH NEARBY M/Z COORDINATES
                MergedDIASpecCoordIndices = np.searchsorted(DIASpectrum[:,0]+tol*DIASpectrum[:,0],DIASpectrum[:,0])
                MergedDIASpecCoords = DIASpectrum[np.unique(MergedDIASpecCoordIndices),0]
                MergedDIASpecIntensities = [np.mean(DIASpectrum[np.where(MergedDIASpecCoordIndices == i)[0],1]) for i in np.unique(MergedDIASpecCoordIndices)]
                DIASpectrum = np.array((MergedDIASpecCoords,MergedDIASpecIntensities)).transpose()
              
                #FILTER LIBRARY SPECTRA BY THE CONDITION THAT SOME NUMBER OF THEIR 10 MOST INTENSE PEAKS BELONG TO THE DIA SPECTRUM
                CentroidBreaks = np.concatenate((DIASpectrum[:,0]-tol*DIASpectrum[:,0],DIASpectrum[:,0]+tol*DIASpectrum[:,0]))               
                CentroidBreaks = np.sort(CentroidBreaks)
                
                LocateReferenceCoordsInDIA = [np.searchsorted(CentroidBreaks,M[:,0]) for M in CandidateRefSpectraLibrary]
                   
                TopTenPeaksCoordsInDIA = [np.searchsorted(CentroidBreaks,M[np.argsort(-M[:,1])[0:min(10,M.shape[0])],0]) for M in CandidateRefSpectraLibrary] 
                ReferencePeaksInDIA = [i for i in range(len(MassWindowCandidates)) if 
                                            len([a for a in TopTenPeaksCoordsInDIA[i] if a % 2 == 1]) > 5] #min(3,CandidateRefSpectraLibrary[i].shape[0])]     
                ProportionOfReferencePeaksInDIA = [len([a for a in TopTenPeaksCoordsInDIA[i] if a % 2 == 1])/CandidateRefSpectraLibrary[i].shape[0]
                                                             for i in range(len(MassWindowCandidates))]                          
                
                RefPeptideCandidatesLocations = [LocateReferenceCoordsInDIA[i] for i in ReferencePeaksInDIA]   
                RefPeptideCandidateList = [CandidateRefSpectraLibrary[i] for i in ReferencePeaksInDIA]               
                RefPeptideCandidates = [MassWindowCandidates[i] for i in ReferencePeaksInDIA]                
                NormalizedRefPeptideCandidateList = [M[:,1]/sum(M[:,1]) for M in RefPeptideCandidateList]
                           
                RefSpectraLibrarySparseRowIndices = (np.array([i for v in RefPeptideCandidatesLocations for i in v if i % 2 == 1]) + 1)/2                 
                RefSpectraLibrarySparseRowIndices = RefSpectraLibrarySparseRowIndices - 1 #Respect the 0-indexing                
                RefSpectraLibrarySparseColumnIndices = np.array([i for j in range(len(RefPeptideCandidates)) for i in [j]*len([k for k in RefPeptideCandidatesLocations[j] if k % 2 == 1])]) 
                RefSpectraLibrarySparseMatrixEntries = np.array([NormalizedRefPeptideCandidateList[k][i] for k in range(len(NormalizedRefPeptideCandidateList)) for i in range(len(NormalizedRefPeptideCandidateList[k])) 
                                                                         if RefPeptideCandidatesLocations[k][i] % 2 == 1])
                
                if (len(RefSpectraLibrarySparseRowIndices) > 0 and len(RefSpectraLibrarySparseColumnIndices) > 0 and len(RefSpectraLibrarySparseMatrixEntries) > 0):                                                             
                      
                    UniqueRowIndices = [i for i in set(RefSpectraLibrarySparseRowIndices)]
                    UniqueRowIndices.sort()
                    
                    DIASpectrumIntensities=DIASpectrum[UniqueRowIndices,1]  #Project the spectrum to those m/z bins at which at least one column of the coefficient matrix has a nonzero entry
                    DIASpectrumIntensities=np.append(DIASpectrumIntensities,[0])    #Add a zero to the end of the DIA data vector to penalize 
                                                                                    #peaks of library spectra not present in the DIA spectrum                
                    
                    
                    #AUGMENT THE LIBRARY MATRIX WITH TOTAL ION INTENSITIES OF PEAKS OF LIBRARY SPECTRA THAT DON'T CORRESPOND TO PEAKS IN DIA SPECTRUM
                    ReferencePeaksNotInDIA = np.array([k for v in RefPeptideCandidatesLocations for k in range(len(v)) if v[k] % 2 == 0])                             
                    SparseColumnIndicesForPeaksNotInDIA = np.arange(len(RefPeptideCandidates))
                    NumRowsOfLibraryMatrix = max(UniqueRowIndices)
                    SparseRowIndicesForPeaksNotInDIA = [NumRowsOfLibraryMatrix+1]*len(SparseColumnIndicesForPeaksNotInDIA)                                   
                    #Duplicate (i,j) entries are summed together, yielding total ion intensities                
                    SparseMatrixEntriesForPeaksNotInDIA = np.array([np.sum([NormalizedRefPeptideCandidateList[j][k] 
                                                                for k in range(len(NormalizedRefPeptideCandidateList[j])) 
                                                                if RefPeptideCandidatesLocations[j][k] % 2 == 0]) 
                                                                for j in range(len(NormalizedRefPeptideCandidateList))])
                    
                    SparseRowIndices=np.append(RefSpectraLibrarySparseRowIndices,SparseRowIndicesForPeaksNotInDIA)
                    SparseColumnIndices=np.append(RefSpectraLibrarySparseColumnIndices,SparseColumnIndicesForPeaksNotInDIA)
                    SparseMatrixEntries=np.append(RefSpectraLibrarySparseMatrixEntries,SparseMatrixEntriesForPeaksNotInDIA)    
                                  
                    SparseRowIndices = stats.rankdata(SparseRowIndices,method='dense').astype(int) - 1 #Renumber the row indices according to the projected spectrum,
                                                                                                    #respecting the 0-indexing                
                    LibrarySparseMatrix = sparse.coo_matrix((SparseMatrixEntries,(SparseRowIndices,SparseColumnIndices)))
                    LibraryCoeffs = sparse_nnls.lsqnonneg(LibrarySparseMatrix,DIASpectrumIntensities,{'show_progress': False})               
                    LibraryCoeffs = LibraryCoeffs['x']
                                
            NonzeroCoeffs = [c for c in LibraryCoeffs if c != 0]
            NonzeroCoeffsAboveThreshold = NonzeroCoeffs
            
            Output = [[0,index,0,0,0,0]]   
        
            if len(NonzeroCoeffs) > 0:        
                RefSpectraIDs = [RefPeptideCandidates[j] for j in range(len(RefPeptideCandidates)) if LibraryCoeffs[j] != 0]
                Output = [[NonzeroCoeffsAboveThreshold[i],index,RefSpectraIDs[i][0],RefSpectraIDs[i][1],precMZ,precRT] for i in range(len(NonzeroCoeffsAboveThreshold))]
                
            yield Output

def RegressSpectraOntoLibraryWithDecoys(DIASpectraIterator,Library,tol,maxWindowOffset):
        
        RefSpectraLibrary = Library.value 
        
        for DIASpectrum in DIASpectraIterator:        
        
            precMZ = float(DIASpectrum[1])     
            precRT = float(DIASpectrum[2])  #MS2 scan retention time, in minutes
            index = DIASpectrum[3]
            windowWidth = DIASpectrum[4]            
            
            DIASpectrum = np.array(DIASpectrum[0])
            
            LibraryCoeffs = []
            
            if len(DIASpectrum.shape) == 2:
                
                if windowWidth > 0:
                    CandidateRefSpectraLibrary = [spectrum['Spectrum'] for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                abs(float(spectrum['PrecursorMZ']) - precMZ) < windowWidth/2]
                    MassWindowCandidates = [key for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                abs(float(spectrum['PrecursorMZ']) - precMZ) < windowWidth/2]
                    CandidateDecoyLibrary = [spectrum['Spectrum'] for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                        windowWidth/2 <= abs(float(spectrum['PrecursorMZ']) - precMZ) <= windowWidth]
                    MassWindowDecoyCandidates = [("DECOY_"+key[0],key[1]) for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                        windowWidth/2 <= abs(float(spectrum['PrecursorMZ']) - precMZ) <= windowWidth]
                else:
                    CandidateRefSpectraLibrary = [spectrum['Spectrum'] for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                        float(spectrum['PrecursorMZ']) > precMZ - maxWindowOffset/2]
                    MassWindowCandidates = [key for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                        float(spectrum['PrecursorMZ']) > precMZ - maxWindowOffset/2]
                    CandidateDecoyLibrary = [spectrum['Spectrum'] for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                       precMZ - maxWindowOffset <= float(spectrum['PrecursorMZ']) <= precMZ - maxWindowOffset/2]                 
                    MassWindowDecoyCandidates = [("DECOY_"+key[0],key[1]) for key,spectrum in RefSpectraLibrary.iteritems() if 
                                                                        precMZ - maxWindowOffset <= float(spectrum['PrecursorMZ']) <= precMZ - maxWindowOffset/2]
                 
                #FILTER LIBRARY SPECTRA BY THE CONDITION THAT SOME NUMBER OF THEIR 10 MOST INTENSE PEAKS BELONG TO THE DIA SPECTRUM
                CentroidBreaks = np.concatenate((DIASpectrum[:,0]-tol*DIASpectrum[:,0],DIASpectrum[:,0]+tol*DIASpectrum[:,0]))               
                CentroidBreaks.sort()
                
                LocateReferenceCoordsInDIA = [np.searchsorted(CentroidBreaks,M[:,0]) for M in CandidateRefSpectraLibrary]
                #Hard cutoff - at least 5 of the 10 most intense peaks (or all peaks if there are fewer than 3) of reference spectrum must appear in acquired spectrum
                
                TopTenPeaksCoordsInDIA = [np.searchsorted(CentroidBreaks,M[np.argsort(-M[:,1])[0:min(10,M.shape[0])],0]) for M in CandidateRefSpectraLibrary] 
                ReferencePeaksInDIA = [i for i in range(len(MassWindowCandidates)) if 
                                            len([a for a in TopTenPeaksCoordsInDIA[i] if a % 2 == 1]) > 5] #min(3,CandidateRefSpectraLibrary[i].shape[0])]     
                

                #SHIFT ALL FRAGMENT ION PEAKS OF ALL DECOY SPECTRA BY 20 M/Z TO ENSURE DISSIMILARITY FROM REAL SPECTRA
                LocateDecoyCoordsInDIA = [np.searchsorted(CentroidBreaks,M[:,0]+20) for M in CandidateDecoyLibrary]
                TopTenPeaksCoordsInDIA = [np.searchsorted(CentroidBreaks,M[np.argsort(-M[:,1])[0:min(10,M.shape[0])],0]+20) for M in CandidateDecoyLibrary]
                DecoyPeaksInDIA = [i for i in range(len(MassWindowDecoyCandidates)) if 
                                            len([a for a in TopTenPeaksCoordsInDIA[i] if a % 2 == 1]) > 5] #min(3,CandidateRefSpectraLibrary[i].shape[0])]     
                               
                
                RefPeptideCandidatesLocations = [LocateReferenceCoordsInDIA[i] for i in ReferencePeaksInDIA]   
                RefPeptideCandidateList = [CandidateRefSpectraLibrary[i] for i in ReferencePeaksInDIA]               
                RefPeptideCandidates = [MassWindowCandidates[i] for i in ReferencePeaksInDIA]                
                NormalizedRefPeptideCandidateList = [M[:,1]/sum(M[:,1]) for M in RefPeptideCandidateList]
                
                DecoyCandidatesLocations = [LocateDecoyCoordsInDIA[i] for i in DecoyPeaksInDIA]
                DecoyCandidateList = [CandidateDecoyLibrary[i] for i in DecoyPeaksInDIA]               
                DecoyCandidates = [MassWindowDecoyCandidates[i] for i in DecoyPeaksInDIA]                
                NormalizedDecoyCandidateList = [M[:,1]/sum(M[:,1]) for M in DecoyCandidateList]
                
                RefSpectraLibrarySparseRowIndices = (np.array([i for v in RefPeptideCandidatesLocations for i in v if i % 2 == 1]) + 1)/2                 
                RefSpectraLibrarySparseRowIndices = RefSpectraLibrarySparseRowIndices - 1 #Respect the 0-indexing                
                RefSpectraLibrarySparseColumnIndices = np.array([i for j in range(len(RefPeptideCandidates)) for i in [j]*len([k for k in RefPeptideCandidatesLocations[j] if k % 2 == 1])]) 
                RefSpectraLibrarySparseMatrixEntries = np.array([NormalizedRefPeptideCandidateList[k][i] for k in range(len(NormalizedRefPeptideCandidateList)) for i in range(len(NormalizedRefPeptideCandidateList[k])) 
                                                                         if RefPeptideCandidatesLocations[k][i] % 2 == 1])
                
                if (len(RefSpectraLibrarySparseRowIndices) > 0 and len(RefSpectraLibrarySparseColumnIndices) > 0 and len(RefSpectraLibrarySparseMatrixEntries) > 0):                                                             
                    DecoyLibrarySparseRowIndices = (np.array([i for v in DecoyCandidatesLocations for i in v if i % 2 == 1]) + 1)/2                 
                    DecoyLibrarySparseRowIndices = DecoyLibrarySparseRowIndices - 1 #Respect the 0-indexing                
                    DecoyLibrarySparseColumnIndices = max(RefSpectraLibrarySparseColumnIndices) + 1 + np.array([i for j in range(len(DecoyCandidates)) for i in [j]*len([k for k in DecoyCandidatesLocations[j] if k % 2 == 1])]) 
                    DecoyLibrarySparseMatrixEntries = np.array([NormalizedDecoyCandidateList[k][i] for k in range(len(NormalizedDecoyCandidateList)) for i in range(len(DecoyCandidatesLocations[k])) 
                                                                             if DecoyCandidatesLocations[k][i] % 2 == 1])
                                                                 
                             
                    UniqueRowIndices = [i for i in set(np.concatenate((RefSpectraLibrarySparseRowIndices,DecoyLibrarySparseRowIndices)))]
                    UniqueRowIndices.sort()
                    
                    DIASpectrumIntensities=DIASpectrum[UniqueRowIndices,1]  #Project the spectrum to those m/z bins at which at least one column of the coefficient matrix has a nonzero entry

                    DIASpectrumIntensities=np.append(DIASpectrumIntensities,[0])    #Add a zero to the end of the DIA data vector to penalize 
                                                                                    #peaks of library spectra not present in the DIA spectrum                
                       
                    #AUGMENT THE LIBRARY MATRIX WITH TOTAL ION INTENSITIES OF PEAKS OF LIBRARY SPECTRA THAT DON'T CORRESPOND TO PEAKS IN DIA SPECTRUM
                    ReferencePeaksNotInDIA = np.array([k for v in RefPeptideCandidatesLocations for k in range(len(v)) if v[k] % 2 == 0])                             
                    SparseColumnIndicesForPeaksNotInDIA = np.arange(len(RefPeptideCandidates))
                    NumRowsOfLibraryMatrix = max(UniqueRowIndices)
                    SparseRowIndicesForPeaksNotInDIA = [NumRowsOfLibraryMatrix+1]*len(SparseColumnIndicesForPeaksNotInDIA)                                   
                    #Duplicate (i,j) entries are summed together, yielding total ion intensities                
                    SparseMatrixEntriesForPeaksNotInDIA = np.array([np.sum([NormalizedRefPeptideCandidateList[j][k] 
                                                                for k in range(len(NormalizedRefPeptideCandidateList[j])) 
                                                                if RefPeptideCandidatesLocations[j][k] % 2 == 0]) 
                                                                for j in range(len(NormalizedRefPeptideCandidateList))])
                    
                    RefSpectraLibrarySparseRowIndices=np.append(RefSpectraLibrarySparseRowIndices,SparseRowIndicesForPeaksNotInDIA)
                    RefSpectraLibrarySparseColumnIndices=np.append(RefSpectraLibrarySparseColumnIndices,SparseColumnIndicesForPeaksNotInDIA)
                    RefSpectraLibrarySparseMatrixEntries=np.append(RefSpectraLibrarySparseMatrixEntries,SparseMatrixEntriesForPeaksNotInDIA)    
                    
                    DecoyPeaksNotInDIA = np.array([k for v in DecoyCandidatesLocations for k in range(len(v)) if v[k] % 2 == 0])                             
                    SparseColumnIndicesForDecoyPeaksNotInDIA = np.arange(len(DecoyCandidates))
                    NumRowsOfLibraryMatrix = max(UniqueRowIndices)
                    SparseRowIndicesForDecoyPeaksNotInDIA = [NumRowsOfLibraryMatrix+1]*len(SparseColumnIndicesForDecoyPeaksNotInDIA)                                   
                    #Duplicate (i,j) entries are summed together, yielding total ion intensities                
                    SparseMatrixEntriesForDecoyPeaksNotInDIA = np.array([np.sum([NormalizedDecoyCandidateList[j][k] 
                                                                for k in range(len(NormalizedDecoyCandidateList[j])) 
                                                                if DecoyCandidatesLocations[j][k] % 2 == 0]) 
                                                                for j in range(len(NormalizedDecoyCandidateList))])
                   
                    DecoyLibrarySparseRowIndices = np.append(DecoyLibrarySparseRowIndices,SparseRowIndicesForDecoyPeaksNotInDIA) 
                    DecoyLibrarySparseColumnIndices=np.append(DecoyLibrarySparseColumnIndices,max(RefSpectraLibrarySparseColumnIndices)+SparseColumnIndicesForDecoyPeaksNotInDIA + 1)
                    DecoyLibrarySparseMatrixEntries=np.append(DecoyLibrarySparseMatrixEntries,SparseMatrixEntriesForDecoyPeaksNotInDIA) 
                    
                    SparseRowIndices = np.concatenate((RefSpectraLibrarySparseRowIndices,DecoyLibrarySparseRowIndices))
                    SparseColumnIndices = np.concatenate((RefSpectraLibrarySparseColumnIndices,DecoyLibrarySparseColumnIndices))
                    SparseMatrixEntries = np.concatenate((RefSpectraLibrarySparseMatrixEntries,DecoyLibrarySparseMatrixEntries))
                    
                    SparseRowIndices = stats.rankdata(SparseRowIndices,method='dense').astype(int) - 1 #Renumber the row indices according to the projected spectrum,
                                                                                                    #respecting the 0-indexing
                  
                    LibrarySparseMatrix = sparse.coo_matrix((SparseMatrixEntries,(SparseRowIndices,SparseColumnIndices)))
                    LibraryCoeffs = sparse_nnls.lsqnonneg(LibrarySparseMatrix,DIASpectrumIntensities,{'show_progress': False})               
                    LibraryCoeffs = LibraryCoeffs['x']
                             
            NonzeroCoeffs = [c for c in LibraryCoeffs if c != 0]
            NonzeroCoeffsAboveThreshold = NonzeroCoeffs
            
            Output = [[0,index,0,0,0,0]]   
        
            if len(NonzeroCoeffs) > 0:        
                RefSpectraIDs = [RefPeptideCandidates[j] for j in range(len(RefPeptideCandidates)) if LibraryCoeffs[j] != 0]
                DecoyIDs = [DecoyCandidates[j] for j in range(len(DecoyCandidates)) if LibraryCoeffs[max(RefSpectraLibrarySparseColumnIndices)+1+j] != 0]
                
                RefSpectraIDs = RefSpectraIDs+DecoyIDs
                Output = [[NonzeroCoeffsAboveThreshold[i],index,RefSpectraIDs[i][0],RefSpectraIDs[i][1],precMZ,precRT] for i in range(len(NonzeroCoeffsAboveThreshold))]
                
            yield Output


if __name__ == "__main__":
    
    args = sys.argv

    mzMLname = args[1]  #dirName = args[1]
    
    libName = args[2]
    
    Index = 0       
    if len(args) >= 4:
        Index = int(args[3])
    
    StartOrEnd = "start" 
    if len(args) >= 5:
        StartOrEnd = args[4]
    
    numPartitions = 200
    if len(args) >= 6:
        numPartitions = int(args[5])
    
    instrument = 'orbitrap'    
    if len(args) == 7:
        instrument=args[6]

    delta = 10
    if len(args) == 8:
    	delta=float(args[7])
        
    #Cast the spectral library as a dictionary  
        
    start = time.time()
    
    libPath = os.path.expanduser(libName+'.blib')
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

    else:
        SpectraLibrary = pickle.load(open(libName,"rb"))
        
    print "Library loaded in {} minutes".format(round((time.time()-start)/60,1))
    
    path = os.path.expanduser(mzMLname+'.mzML')  
    DIArun = pymzml.run.Reader(path)
    E = enumerate(DIArun)

    start = time.time()     

    if StartOrEnd == "end":
        if instrument == 'orbitrap':
                res = [[spectrum.peaks,spectrum['precursors'][0]['mz'],spectrum['scan time'],i] for i,spectrum in E if spectrum['ms level'] == 2.0 and i < Index]
        if instrument == 'tof':
                res = [[spectrum.peaks,spectrum['precursors'][0]['mz'],spectrum['scan time'],i] for i,spectrum in E if 'precursors' in spectrum.keys() and i < Index]
    else:
        if instrument == 'orbitrap':
                res = [[spectrum.peaks,spectrum['precursors'][0]['mz'],spectrum['scan time'],i] for i,spectrum in E if spectrum['ms level'] == 2.0 and i >= Index]
        if instrument == 'tof':
                res = [[spectrum.peaks,spectrum['precursors'][0]['mz'],spectrum['scan time'],i] for i,spectrum in E if 'precursors' in spectrum.keys() and i >= Index]


    print "Loaded {} MS2 spectra from {} in {} minutes.".format(len(res),path,round((time.time()-start)/60,1))
    
    res=[[res[i][0],float(res[i][1]),float(res[i][2]),float(res[i][3]),float(res[i+1][1]) - float(res[i][1])] for i in range(len(res)-1)]
    header=[[x[1],x[2],x[3]] for x in res]

    absolutePath = mzMLname.rsplit('/',1)[0]
    noPathName = mzMLname.rsplit('/',1)[1]
    if '/' in libName:
        libName = libName.rsplit('/',1)[1]

    outputDir = os.path.expanduser(absolutePath+'/SpecterResults')
    if not os.path.exists(outputDir):
    	os.makedirs(outputDir)

    headerPath = os.path.expanduser(absolutePath+'/SpecterResults/'+noPathName+'_'+libName+'_header.csv')     
    
    with open(headerPath, "ab") as f:
        writer = csv.writer(f)
        writer.writerows(header)      
           
    print "Output written to {}.".format(headerPath)

    MaxWindowPrecMZ = max(np.array([x[1] for x in res])) + max(np.array([x[4] for x in res]))
    MaxOffset = max(np.array([x[4] for x in res]))

    SpectraLibrary = {k:SpectraLibrary[k] for k in SpectraLibrary if SpectraLibrary[k]['PrecursorMZ'] < MaxWindowPrecMZ}

    conf = (SparkConf().set("spark.driver.maxResultSize", "25g"))
        
    sc = SparkContext(conf=conf,appName="Specter",pyFiles=['sparse_nnls.py'])
    
    #Recast the library as a broadcast variable to improve performance
    BroadcastLibrary = sc.broadcast(SpectraLibrary)  
    
    res = sc.parallelize(res, numPartitions)
    
    output = res.mapPartitions(partial(RegressSpectraOntoLibrary, Library=BroadcastLibrary, tol=delta*1e-6, maxWindowOffset = MaxOffset)).collect()  
    
    output = [[output[i][j][0],output[i][j][1],output[i][j][2],output[i][j][3],
                        output[i][j][4],output[i][j][5]] for i in range(len(output)) for j in range(len(output[i]))]
    
    outputPath = os.path.expanduser(absolutePath+'/SpecterResults/'+noPathName+'_'+libName+'_SpecterCoeffs.csv')     
    
    with open(outputPath, "ab") as f:
        writer = csv.writer(f)
        writer.writerows(output)      
           
    print "Output written to {}.".format(outputPath)

    output = res.mapPartitions(partial(RegressSpectraOntoLibraryWithDecoys, Library=BroadcastLibrary, tol=delta*1e-6, maxWindowOffset = MaxOffset)).collect()  
    
    output = [[output[i][j][0],output[i][j][1],output[i][j][2],output[i][j][3],
                        output[i][j][4],output[i][j][5]] for i in range(len(output)) for j in range(len(output[i]))]
    
    outputPath = os.path.expanduser(absolutePath+'/SpecterResults/'+noPathName+'_'+libName+'_SpecterCoeffsDecoys.csv')     
    
    with open(outputPath, "ab") as f:
        writer = csv.writer(f)
        writer.writerows(output)      
           
    print "Output written to {}.".format(outputPath)
