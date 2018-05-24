SparkSpecterCoeffs = function(rawCoeffsFilepath,ExperimentHeader=NULL, lib=NULL,MS1=FALSE)
{ library(data.table)
  
  DIA_RefSpectraCoeffs = read.csv(rawCoeffsFilepath,header=FALSE,stringsAsFactors = FALSE)
  
  names(DIA_RefSpectraCoeffs) = c("Coeff","Scan","RefPeptideSeq", "RefPrecursorCharge","PrecursorMZ","RetentionTime") 
  
  DIA_RefSpectraCoeffs = DIA_RefSpectraCoeffs[c(2,1,3,4,5,6)]
  
  DIA_RefSpectraCoeffs$Coeff = as.numeric(DIA_RefSpectraCoeffs$Coeff)
  DIA_RefSpectraCoeffs$Scan = strtoi(DIA_RefSpectraCoeffs$Scan) + 1
  DIA_RefSpectraCoeffs$RefPeptideSeq = as.character(DIA_RefSpectraCoeffs$RefPeptideSeq)
  DIA_RefSpectraCoeffs$RefPrecursorCharge = strtoi(DIA_RefSpectraCoeffs$RefPrecursorCharge)
  DIA_RefSpectraCoeffs$PrecursorMZ = as.numeric(DIA_RefSpectraCoeffs$PrecursorMZ)
  DIA_RefSpectraCoeffs$RetentionTime = as.numeric(DIA_RefSpectraCoeffs$RetentionTime)
  
  DIA_RefSpectraCoeffs = na.omit(DIA_RefSpectraCoeffs)
  
  DIA_RefSpectraCoeffs = DIA_RefSpectraCoeffs[with(DIA_RefSpectraCoeffs,order(Scan, RefPeptideSeq, RefPrecursorCharge)),]
  
  #Convert minutes to seconds if need be
  if (max(DIA_RefSpectraCoeffs$RetentionTime) < 600)
    DIA_RefSpectraCoeffs$RetentionTime = 60*DIA_RefSpectraCoeffs$RetentionTime
  
  DIA_RefSpectraCoeffs = DIA_RefSpectraCoeffs[which(DIA_RefSpectraCoeffs$Coeff != 0),]
  DIA_RefSpectraCoeffs = data.table(DIA_RefSpectraCoeffs)
  setkey(DIA_RefSpectraCoeffs,RefPeptideSeq,RefPrecursorCharge)
  
  return(DIA_RefSpectraCoeffs)
  
}



QuantifyAllFromSpecterCoeffs = function(SpecterData,header) {
  require(kza)
  require(pracma)
  require(zoo)
  require(moments)
  require(data.table)
  

  IDs = data.frame(unique(SpecterData[,c('RefPeptideSeq','RefPrecursorCharge'),with=FALSE]))
  IsThereEnoughData = unlist(Map(function(i) EnoughData(SpecterData,IDs,i,RTWindow=FALSE),
                                 seq(1,nrow(IDs))))
  
  f = function(k) {if (IsThereEnoughData[k]) {
    SpecterQuantifyPeptides(SpecterData,IDs,k,header,RTWindow=FALSE)} else {list(0,1,0,0,0,0,1,0,0,0)}}
  
  Quant = Map(f,seq(1,nrow(IDs)))
  Quants = IDs
  
  Quants$SpecterAreaQuant = as.numeric(unlist(lapply(Quant,function(l) l[[1]][1])))
  Quants$SpecterLBPval = as.numeric(unlist(lapply(Quant,function(l) l[[2]][1])))
  Quants$SpecterSNR = as.numeric(unlist(lapply(Quant,function(l) l[[3]][1])))
  Quants$MaxSpecterCoeffTime = as.numeric(unlist(lapply(Quant,function(l) l[[4]][1])))
  Quants$MaxSpecterCoeff = as.numeric(unlist(lapply(Quant,function(l) l[[5]][1])))
  Quants$SpecterPeakWidth = as.numeric(unlist(lapply(Quant,function(l) l[[6]][1])))
  Quants$SpecterLBPvalOnUniformGrid = as.numeric(unlist(lapply(Quant,function(l) l[[7]][1])))
  Quants$PeakVariance = as.numeric(unlist(lapply(Quant,function(l) l[[8]][1])))
  Quants$PeakSkewness = as.numeric(unlist(lapply(Quant,function(l) l[[9]][1])))
  Quants$PeakKurtosis = as.numeric(unlist(lapply(Quant,function(l) l[[10]][1])))
  return(Quants)
}

###################### PROFILE SMOOTHING AND PLOTTING ###########################

SmoothProfile = function(Data,Identifiers,i,header,IntensityCutoff=0,QuantileCutoff=0,RTWindow=TRUE,param=1e-2,type="rollmean",
                              Plot=TRUE,FilterWindow= 3,KZiters=3,verbose=FALSE,MS1=FALSE,Experiment=NULL,Header=NULL) {
  PeptideData = data.frame(Data[list(Identifiers$RefPeptideSeq[i],Identifiers$RefPrecursorCharge[i])])
  
  PrecMZ = unique(PeptideData$RefPrecursorMZ)
  
  PeptideDataNoCutoff = data.frame(PeptideData)
  
  if (RTWindow)
  {PeptideData = PeptideData[which(abs(PeptideData$RetentionTime - PeptideData$RefRetentionTime) < 300),]
  PeptideDataNoCutoff = PeptideDataNoCutoff[which(abs(PeptideDataNoCutoff$RetentionTime - PeptideDataNoCutoff$RefRetentionTime) < 300),]}
  
  if (verbose)
    print(PeptideData)
    
  PeptideData = PeptideData[which(PeptideData$Coeff >= IntensityCutoff),]
  PeptideData = PeptideData[which(PeptideData$Coeff >= quantile(PeptideData$Coeff,QuantileCutoff)),]
  
  PeptideDataNoCutoff = PeptideDataNoCutoff[c("RetentionTime","Coeff","PrecursorMZ")]  
  PeptideDataNoCutoff = PeptideDataNoCutoff[!duplicated(PeptideDataNoCutoff$RetentionTime),]
  
  PeptideData = PeptideData[c("Scan","RetentionTime","Coeff","PrecursorMZ")]
  PeptideData = PeptideData[!duplicated(PeptideData$RetentionTime),]
  
  SmoothFrame = data.frame(RetentionTime=numeric(),Coeff=numeric())
  KZData = NULL
  
  if (type == "SG")
  {SmoothedData = sgolayfilt(PeptideDataOnGrid[,2],p=0,n=3)
  SmoothFrame = data.frame(RetentionTime=PeptideDataOnGrid[,1],Coeff=SmoothedData)
  }
  
  if (type == "rollmean")
  {ZooData = zoo(PeptideData$Coeff, PeptideData$RetentionTime)
  SmoothedData = zoo(0,0)

  
  PeptideData$Coeff[which(PeptideData$Coeff < 1)] = 0

  
  if (nrow(PeptideData) > 3) {
   
  MinScanIndex = min(PeptideData$Scan,na.rm=TRUE)
  MaxScanIndex = max(PeptideData$Scan,na.rm=TRUE)
  
  HeaderRange = data.frame(header[J(MinScanIndex:MaxScanIndex)])  #header needs to be a data.table with key seqNum
  HeaderRange = subset(HeaderRange,msLevel == 2 & floor(precursorMZ) %in% floor(unique(PeptideData$PrecursorMZ)))
  
  DataOnUniformGrid = HeaderRange[c("seqNum","retentionTime")]
  DataOnUniformGrid$Coeff = numeric(nrow(DataOnUniformGrid))
  DataOnUniformGrid$Coeff[which(DataOnUniformGrid$seqNum %in% PeptideData$Scan)] = PeptideData$Coeff  
  
  DataOnUniformGrid = DataOnUniformGrid$Coeff

  KZData = kz(DataOnUniformGrid,m=FilterWindow,k=KZiters)}

  }
  
  BoxTestPval = 1 

  if (Plot) {
    n = nrow(PeptideDataNoCutoff)
    m = nrow(PeptideData)

    CombinedData = rbind(PeptideDataNoCutoff,PeptideData) 
    CombinedData = rbind(CombinedData,KZFrame)
    
    CombinedData$Type = paste0("Kolmogorov-Zurbenko"," (",KZiters," iterations)")
    CombinedData$Type[1:n] = "Raw"
    CombinedData$Type[(n+1):(n+m)] = "Cutoff"
    
    CombinedData$Type = factor(CombinedData$Type, levels=c("Raw","Cutoff",paste0("Kolmogorov-Zurbenko"," (",KZiters," iterations)")))
    if (MS1)
    {
      MS1Data = Header[which(Header$msLevel == 1),]
      PrecursorMS1Candidates = findInterval(PeptideDataNoCutoff$RetentionTime,MS1Data$retentionTime)
      PrecursorMS1CandidateScans = mzR::peaks(Experiment,scans=which(Header$msLevel == 1)[PrecursorMS1Candidates])
      MS1Intensity = data.frame(RetentionTime=numeric(),Coeff=numeric())
      for (k in 1:length(PrecursorMS1Candidates))  
      {if (any(abs(PrecursorMS1CandidateScans[[k]][,1] - PrecMZ) < 1e-5*PrecMZ))
      {ClosestPeak = which(abs(PrecursorMS1CandidateScans[[k]][,1] - PrecMZ) == min(abs(PrecursorMS1CandidateScans[[k]][,1] - PrecMZ)))
      ClosestPeakIntensity = PrecursorMS1CandidateScans[[k]][ClosestPeak,2]
      MS1Intensity = rbind(MS1Intensity,data.frame(RetentionTime=MS1Data$retentionTime[PrecursorMS1Candidates[k]],Coeff=ClosestPeakIntensity,
                                                   PrecursorMZ = PrecursorMS1CandidateScans[[k]][ClosestPeak,1]))
      }
      }
      

      MS1Intensity$Type = "MS1 (closest precursor m/z)"
      
      CombinedData = rbind(CombinedData,MS1Intensity)
      CombinedData$Type = factor(CombinedData$Type, levels=c("Raw","Cutoff","Moving average",paste0("Kolmogorov-Zurbenko"," (",KZiters," iterations)"),
                                                             "MS1 (closest precursor m/z)"))
      
    }
    
    p = ggplot(CombinedData,aes(x=RetentionTime,y=Coeff,group=1)) + 
      geom_line(size=1) + facet_wrap(~Type,ncol=1) + theme(legend.title=element_blank()) +
      ggtitle(paste0(Identifiers$RefPeptideSeq[i],", Charge = ",Identifiers$RefPrecursorCharge[i],
                     ", RefPrecursorMZ = ",round(PrecMZ,digits=4),
                     ", Cutoff = ",quantile(PeptideData$Coeff,QuantileCutoff),
                     ", Box test p-val = ",BoxTestPval))
    
    print(p)
  }
  
  return(KZData)
}


###################### PEAK DETECTION ###########################

FindPeaksInData = function(Data,Identifiers,i,header,Multiplexed=FALSE,IntensityCutoff=0,
                           QuantileCutoff=0,RTWindow=TRUE,smooth="rollmean",FilterWindow= 3,KZiters=3) {
  
  PeptideData = data.frame(Data[list(Identifiers$RefPeptideSeq[i],Identifiers$RefPrecursorCharge[i])])
  
   if (RTWindow)
     PeptideData = subset(PeptideData, abs(RetentionTime - RefRetentionTime) < 300)
  
  PeptideData = PeptideData[c("Scan","RetentionTime","Coeff","PrecursorMZ")]
  PeptideData = PeptideData[!duplicated(PeptideData$RetentionTime),]
  
  PeptideData$Coeff[which(PeptideData$Coeff < 1)] = 0
  
  PeptidePeaks = NULL
  if (length(which(PeptideData$Coeff > 0) > 3))  {
  
  MinScanIndex = min(PeptideData$Scan,na.rm=TRUE)
  MaxScanIndex = max(PeptideData$Scan,na.rm=TRUE)
  
  HeaderRange = data.frame(header[J(MinScanIndex:MaxScanIndex)])  #header needs to be a data.table with key seqNum
  HeaderRange = subset(HeaderRange,msLevel == 2 & floor(precursorMZ) %in% floor(unique(PeptideData$PrecursorMZ)))
    
  DataOnUniformGrid = HeaderRange[c("seqNum","retentionTime")]
  DataOnUniformGrid$Coeff = numeric(nrow(DataOnUniformGrid))
  DataOnUniformGrid$Coeff[which(DataOnUniformGrid$seqNum %in% PeptideData$Scan)] = PeptideData$Coeff  
  
  DataOnUniformGrid = DataOnUniformGrid$Coeff
  
  PeptidePeaks = findpeaks(DataOnUniformGrid, nups = 2, ndowns = 2, npeaks= 10,sortstr = TRUE)
  
  if (smooth == "rollmean")
  {
  KZData = kz(DataOnUniformGrid,m=FilterWindow,k=KZiters)
  KZData[which(DataOnUniformGrid < 1)] = 0   #Remove artifical nonzero coeffs resulting from smoothing

  PeptidePeaks = findpeaks(KZData, nups = 2, ndowns = 2, npeaks= 10,sortstr = TRUE)
  
   }
  
  }
  return(PeptidePeaks)
}

###################### QUANTITATION ###########################

SpecterQuantifyPeptides = function(Data,Identifiers,i,header,IntensityCutoff=0,QuantileCutoff=0,RTWindow=TRUE,smooth="rollmean",FilterWindow= 3,KZiters=3) {
    
    PeptideData = data.frame(Data[list(Identifiers$RefPeptideSeq[i],Identifiers$RefPrecursorCharge[i])])
    
     if (RTWindow)
       PeptideData = subset(PeptideData, abs(RetentionTime - RefRetentionTime) < 300)
    
    PeptideData = PeptideData[c("Scan","RetentionTime","Coeff","PrecursorMZ")]
    PeptideData = PeptideData[!duplicated(PeptideData$RetentionTime),]
    
    PeptideData$Coeff[which(PeptideData$Coeff < 1)] = 0  #Total ion intensities less than 1 are discarded as physically meaningless.
    
    area = 0 
        
    RawPeaks = FindPeaksInData(Data,Identifiers,i,header,IntensityCutoff=IntensityCutoff,
                               QuantileCutoff=QuantileCutoff,RTWindow=RTWindow,smooth="none",FilterWindow=FilterWindow)
    
    SmoothPeaks = FindPeaksInData(Data,Identifiers,i,header,
                                  IntensityCutoff=IntensityCutoff,QuantileCutoff=QuantileCutoff,RTWindow=RTWindow,smooth=smooth,FilterWindow=FilterWindow,KZiters=KZiters)
        
    BoxTestPval = 1
    BoxTestPvalOnGrid = 1
    MaxCoeff= 0
    TimeAtTopOfPeak = 0
    snr = 0
    PeakCentralMoments = numeric(3)
    result = list(Area=area,BoxTestPval=BoxTestPval,SNR=snr,RT=TimeAtTopOfPeak,MaxCoeff=0,PeakWidth=0,BoxTestPvalOnGrid=BoxTestPvalOnGrid,
                  Variance=0,Skewness=0,Kurtosis=0)
    
    if (length(which(PeptideData$Coeff > 1) > 3) & !is.null(SmoothPeaks)) {
      MinScanIndex = min(PeptideData$Scan,na.rm=TRUE)
      MaxScanIndex = max(PeptideData$Scan,na.rm=TRUE)
     
      HeaderRange = data.frame(header[J(MinScanIndex:MaxScanIndex)])    #header needs to be a data.table with key seqNum
      HeaderRange = subset(HeaderRange,msLevel == 2 & floor(precursorMZ) %in% floor(unique(PeptideData$PrecursorMZ)))

      DataOnUniformGrid = HeaderRange[c("seqNum","retentionTime")]
      DataOnUniformGrid$Coeff = numeric(nrow(DataOnUniformGrid))
      DataOnUniformGrid$Coeff[which(DataOnUniformGrid$seqNum %in% PeptideData$Scan)] = PeptideData$Coeff  
      
      start = SmoothPeaks[1,3]
      end = SmoothPeaks[1,4]
      
      PeakDataOnUniformGrid = DataOnUniformGrid[start:end,]
      
      area = trapz(PeakDataOnUniformGrid$retentionTime, PeakDataOnUniformGrid$Coeff) 
      
      TimeAtTopOfPeak = DataOnUniformGrid$retentionTime[SmoothPeaks[1,2]]
      MaxCoeff=DataOnUniformGrid$Coeff[SmoothPeaks[1,2]]
      
      start = SmoothPeaks[1,3]
      end = SmoothPeaks[1,4]
      
      PeakDataOnUniformGrid = DataOnUniformGrid[start:end,]
      PeakCentralMoments = all.moments(scale(PeakDataOnUniformGrid$Coeff,center=FALSE),order.max = 4,central = TRUE)[3:5]
      
      BoxTest = Box.test(PeptideData$Coeff,type="Ljung-Box")
      BoxTestPval = round(BoxTest$p.value,digits=5)
      
      BoxTestOnGrid = Box.test(DataOnUniformGrid$Coeff,type="Ljung-Box")
      BoxTestPvalOnGrid = round(BoxTestOnGrid$p.value,digits=5)
      snr=area/sd(DataOnUniformGrid$Coeff[-seq(start,end)])  #Denominator is the standard deviation of the coefficients not falling within the highest peak

      peakWidth = DataOnUniformGrid$retentionTime[end] - DataOnUniformGrid$retentionTime[start]

      result = list(Area=area,BoxTestPval=BoxTestPval,SNR=snr,RT=TimeAtTopOfPeak,MaxCoeff=MaxCoeff,PeakWidth=peakWidth,
                    BoxTestPvalOnGrid=BoxTestPvalOnGrid,Variance=PeakCentralMoments[1],Skewness=PeakCentralMoments[2],Kurtosis=PeakCentralMoments[3])
    }
    return(result)
}
  

EnoughData = function(Data,Identifiers,i,RTWindow=TRUE) {
    PeptideData = Data[list(Identifiers$RefPeptideSeq[i],Identifiers$RefPrecursorCharge[i])]
    
    if (RTWindow)
      PeptideData = subset(PeptideData, abs(RetentionTime - RefRetentionTime) < 300)
    
    if (nrow(PeptideData) > 0)
    {return(TRUE)} else {
      return(FALSE)
    }
  }


library(MASS)
library(data.table)
args=commandArgs(TRUE)
mzmlPath = paste0(args[1],".mzML")
dirName = strsplit(mzmlPath,split="/")[[1]]
mzmlName= gsub(".mzML","",dirName[length(dirName)])
libName = strsplit(args[2],split="/")[[1]]
libName = libName[length(libName)]
dirName = dirName[-length(dirName)]
dirName = paste(dirName,collapse="/")
resultsPath = file.path(dirName, "SpecterResults", paste0(mzmlName, "_", libName, "_SpecterCoeffs.csv"))
headerPath = file.path(dirName, "SpecterResults", paste0(mzmlName, "_", libName, "_header.csv"))
decoyResultsPath = paste0(gsub(".csv","",resultsPath),"Decoys.csv")

h = read.csv(headerPath,header=FALSE,stringsAsFactors=FALSE)
names(h) = c("precursorMZ","retentionTime","seqNum")
h$retentionTime = 60*h$retentionTime
h$seqNum = h$seqNum + 1
h$msLevel = 2
h = data.table(h)
setkey(h,seqNum)

results = SparkSpecterCoeffs(resultsPath)
resultsWithDecoys = SparkSpecterCoeffs(decoyResultsPath)

quants = QuantifyAllFromSpecterCoeffs(results,header = h)
quantsWithDecoys = QuantifyAllFromSpecterCoeffs(resultsWithDecoys,header=h)

quantsWithDecoys$Type = ifelse(grepl("DECOY",quantsWithDecoys$RefPeptideSeq),"Decoy","Target")
quantsWithDecoys = quantsWithDecoys[which(quantsWithDecoys$SpecterAreaQuant > 0),]

if (length(which(quantsWithDecoys$Type == "Decoy")) > 0) {                                                 
SpecterLDA = lda(Type~ MaxSpecterCoeff + PeakVariance + PeakSkewness + PeakKurtosis,
                 data = quantsWithDecoys)

plda = predict(object = SpecterLDA, newdata = quantsWithDecoys)

D = data.frame(PeptideSeq = quantsWithDecoys$RefPeptideSeq, 
               PrecursorCharge = quantsWithDecoys$RefPrecursorCharge, 
               Type=quantsWithDecoys$Type,
               Score=as.numeric(unlist(plda$x)))


FDR = function(score)  {
  return(length(which(D$Type == "Decoy" & D$Score > score))/length(which(D$Score > score)))
}

fdr = as.numeric(sapply(D$Score,FDR))

cutoff = min(D$Score[which(fdr < 0.01)],na.rm=TRUE)
D = D[which(D$Score >= cutoff),]
quants = quants[which(paste(quants$RefPeptideSeq,quants$RefPrecursorCharge) %in% paste(D$PeptideSeq,D$PrecursorCharge)),]
quants$Score = D$Score[match(paste(quants$RefPeptideSeq,quants$RefPrecursorCharge),paste(D$PeptideSeq,D$PrecursorCharge))]
} else {
  quants$Score = 1
  }
                                                 
outputPath = gsub("_SpecterCoeffs.csv","_SpecterQuants.csv",resultsPath)
write.csv(quants[c("RefPeptideSeq","RefPrecursorCharge","SpecterAreaQuant","Score")],file=outputPath,quote=FALSE,row.names=FALSE)

