#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
## ------------------------------------------------------------------------- ##
#                     READ KSAT PERMEAMETER CONSTANT HEAD DATA                #
# --------------------------------------------------------------------------- #
###############################################################################
#
#  Created by Jan De Pue
#  Ghent University, Belgium
#  May 2017


import numpy
import pylab
import os
import fnmatch
import glob
from scipy import stats
from scipy import optimize

from SoilLab_Funk import *
from HypropFunk import *
from RetentionConductivityCapacity_Funky import *


from matplotlib.backends.backend_pdf import PdfPages
figlist=[]
figsize1=[10./numpy.sqrt(2),10]
figsize2=[10,10./numpy.sqrt(2)]
#dpi=300
dpi=None

font = {'family' : 'monospace',
        'size'   : 15}
pylab.rc('font', **font)

pylab.ioff()



MainRoot = '/home/jan/git/SoilLab/TestData'
# MainRoot = 'C:\\Users\\admin\\Documents\\GitHub\\SoilLab\\TestData'
DoTempCorrection = False


#==============================================================================
# Find Identification Data
#==============================================================================

IDFile = os.path.join(MainRoot,'Identification.csv')

delimiter = '\t'
skip_header = 6
ThisAnalysisType = 'Ks_CH_Permeameter'

## Read header
fID = open(IDFile)
Project = fID.readline().split(delimiter)[1]
Date = fID.readline().split(delimiter)[1]
Collaborators = fID.readline().split(delimiter)[1:]
skipline = fID.readline()
Header1 = numpy.array(fID.readline().replace('\n','').split(delimiter))
Header2 = numpy.array(fID.readline().replace('\n','').split(delimiter))
fID.close()

iH_ID = numpy.where(Header1 == 'Identification')[0][0]
iH_Meta = numpy.where(Header1 == 'Metadata')[0][0]
iH_Ana = numpy.where(Header1 == 'Analysis')[0][0]
iH_Com = numpy.where(Header1 == 'Comments')[0][0]

iH_ThisAnalysisType= numpy.where(Header2 == ThisAnalysisType)[0][0]

## Read Identification Data
ID_AllData = numpy.genfromtxt(IDFile,delimiter = delimiter, dtype ='str', skip_header = skip_header)

FilterThisType = ID_AllData[:,iH_ThisAnalysisType] == '1'
nSamp = FilterThisType.sum()

ID_AllData = ID_AllData[FilterThisType,:]
ID_AllData = CapitalizeArray(ID_AllData)
ID_IDData = ID_AllData[:,iH_ID:iH_Meta]
ID_MetaData = ID_AllData[:,iH_Meta:iH_Ana]
ID_AnalysisData = ID_AllData[:,iH_Ana:iH_Com]
ID_CommentData = ID_AllData[:,iH_Com]

ID_IDHeader = Header2[iH_ID:iH_Meta]

CheckEmptyMetaData = ~numpy.all(ID_MetaData == '',axis=0)
ID_MetaData = ID_MetaData[:,CheckEmptyMetaData]
ID_MetaHeader = Header2[iH_Meta:iH_Ana][CheckEmptyMetaData]


ID_FullMeta = map(lambda x : '_'.join(x),ID_MetaData)


#==============================================================================
# Find SBPP data
#==============================================================================

DataRoot = os.path.join(MainRoot,'Ks_CH_Permeameter')


FileList = []
for root, dirnames, filenames in os.walk(DataRoot):
    for filename in fnmatch.filter(filenames, '*.csv'):
        FileList.append(os.path.join(root, filename))

FileList.sort()
FileList=numpy.array(FileList)
FileNames=numpy.array(map(lambda x : os.path.basename(x)[:-4].upper(),FileList)) 

## FILTER 
# FiltRef=numpy.char.find(FileList, 'Originals') > -1
# FileNames=FileNames[~FiltRef] 
# FileList=FileList[~FiltRef]

nFiles=len(FileList)
print(root+" : %s files found." %nFiles) # 


# nMetaLevel = len(FileNames[0].split('_'))
nMetaLevel = ID_MetaData.shape[1]

iF = 0 # in case this would be extended to multiple files found
filename = FileList[iF]

delimiter = '\t'
skip_header = 9


fID = open(filename)
Project = fID.readline().split(delimiter)[1]
Date = fID.readline().split(delimiter)[1]
Collaborators = fID.readline().split(delimiter)[1:]
skipline = fID.readline()
Header1 = numpy.array(fID.readline().replace('\n','').split(delimiter))
Header2 = numpy.array(fID.readline().replace('\n','').split(delimiter))
Header3 = numpy.array(fID.readline().replace('\n','').split(delimiter))
Header4 = numpy.array(fID.readline().replace('\n','').split(delimiter))
Header5 = numpy.array(fID.readline().replace('\n','').split(delimiter))
fID.close()

iH_ID_Data = numpy.where(Header4 == 'Identification')[0][0]
iH_Meta_Data = numpy.where(Header4 == 'Metadata')[0][0]
iH_Com_Data = numpy.where(Header4 == 'Comments')[0][0]
iH_Data = numpy.where(Header4 == 'SampleDimension')[0][0]

Data_All = numpy.genfromtxt(filename,delimiter = delimiter, dtype ='str', skip_header = skip_header)
Data_All = CapitalizeArray(Data_All)

Data_ID = Data_All[:,iH_ID_Data:iH_Meta_Data]
Data_Meta = Data_All[:,iH_Meta_Data:iH_Com_Data]
Data_Com = Data_All[:,iH_Com_Data:iH_Data]
Data = Data_All[:,iH_Data:].astype('float')
Data[Data_All[:,iH_Data:] == ''] = 0 # replace empty cells with zero (but not nan)
Data_Header1 = Header1[iH_Data:]
Data_Header2 = Header2[iH_Data:]
Data_Header3 = Header3[iH_Data:]
Data_Header4 = Header4[iH_Data:]
Data_Header5 = Header5[iH_Data:]

MetaCode = Data_Meta
    
MetaUnique=[]
for iM in range(nMetaLevel):
    MetaUnique.append(numpy.unique(MetaCode[:,iM]))

nMetaUnique = map(numpy.size,MetaUnique)

nDS = Data_All.shape[0]


#==============================================================================
# Check Match ID
#==============================================================================

iMatchID = numpy.zeros(nDS,dtype ='int')                # Index of match in ID database for each sample
iMatchID_R = numpy.zeros(nSamp,dtype ='int')+numpy.nan  # Index of matching sample for ID database

for iF in range(nDS):
    Meta = MetaCode[iF,:]
    FindMatch = ID_MetaData == Meta
    MatchFound = numpy.all(FindMatch, axis=1)
    if not numpy.any(MatchFound):
        raise ValueError('No match found in Identification Database for %s'%'_'.join(Meta))
    else:
        iMatchID[iF] = numpy.where(MatchFound)[0]
        iMatchID_R[MatchFound] = iF

print('%s / %s samples found'%(nDS,nSamp))


#==============================================================================
# Do Stuff
#==============================================================================

## Open, read data

Height = Data[:,Data_Header5 == 'Diameter(cm)'] * 1e-2       # m
Diameter = Data[:,Data_Header5 == 'Height(cm)'] * 1e-2   # m
HeightCorr = Data[:,Data_Header5 == 'HeightCorrection(cm)'] * 1e-2   # m

## Loop to parse measurements
Mfilt = Data_Header1 == 'measurement'
iMfilt = numpy.where(Mfilt)[0]
nM = Mfilt.sum()

Measurement = numpy.zeros(nM).astype('S256')
Date = numpy.zeros(nM).astype('S256')
Temp = numpy.zeros(nM)
Volume = numpy.zeros((nDS,nM))
Time = numpy.zeros((nDS,nM))
WaterHeightBox = numpy.zeros((nDS,nM))
WaterHeightCup = numpy.zeros((nDS,nM))

for iM in range(nM):
    Measurement[iM] = Data_Header1[iMfilt[iM]+1]
    Date[iM] = Data_Header2[iMfilt[iM]+1]
    Temp[iM] = float(Data_Header3[iMfilt[iM]+1])
    Volume[:,iM] = Data[:,iMfilt[iM]]*1e-6 # cm3 => m3
    Min = Data[:,iMfilt[iM]+1]
    Sec = Data[:,iMfilt[iM]+2]
    Time[:,iM] = Min * 60 + Sec # Sec
    WaterHeightBox[:,iM] = Data[:,iMfilt[iM]+3]*1e-2
    WaterHeightCup[:,iM] = Data[:,iMfilt[iM]+4]*1e-2


## Calculate Ks

Area = numpy.pi * (Diameter/2)**2
flux = Volume/(Time * Area)# m/s
dH = WaterHeightCup - WaterHeightBox
dZ = Height + HeightCorr

Ks = flux * dZ/dH
KsMean = numpy.mean(Ks,axis=1)
KsStd = numpy.std(Ks,axis=1)

iMeas = numpy.arange(nM)



## Temperature Correction
    

if DoTempCorrection:
    ## TEMPERATURE CORRECTION
    #http://www.viscopedia.com/viscosity-tables/substances/water/   (IAPWS 2008)
    # refT=10°C
    visc10=1.3059 # mPa/s
    # refT=20°C
    visc20=1.0016 # mPa/s
    # Ks corrected to 20°C
    Ks=Ks*visc10/visc20
    print('Ks Temperature correction to 20 °C - CORRECTION')
else:
    print('Ks - NO TEMPERATURE CORRECTION')


## Assume last metalevel is repetition => Calculate Geomean

MeanShape = nMetaUnique[:-1]
nMeta_Full = numpy.prod(MeanShape) ## total amount of classes

Ks_geomean=numpy.zeros(MeanShape)+numpy.nan
Ks_geostd=numpy.zeros(MeanShape)+numpy.nan
Ks_geomean_F=numpy.zeros(MeanShape)+numpy.nan
Ks_geostd_F=numpy.zeros(MeanShape)+numpy.nan
Ks_geomean_W=numpy.zeros(MeanShape)+numpy.nan
Ks_geostd_W=numpy.zeros(MeanShape)+numpy.nan
Ks_Min=numpy.zeros(MeanShape)+numpy.nan
Ks_Max=numpy.zeros(MeanShape)+numpy.nan
nK_BeforeFilter=numpy.zeros(MeanShape)+numpy.nan
nK_AfterFilter=numpy.zeros(MeanShape)+numpy.nan
MetaCode_2 = numpy.zeros(MeanShape+[nMetaLevel-1,]).astype('|S256')
    
factor = numpy.zeros(nMetaLevel-1)
for iL in range(nMetaLevel-1):
    factor[iL] = numpy.prod(MeanShape[:iL+1])
    
for iMFull in range(nMeta_Full):
    MetaIndex = numpy.zeros(nMetaLevel-1,dtype=int)
    for iL in range(nMetaLevel-1):
        MetaIndex[iL] = iMFull % factor[iL] - sum(MetaIndex[:iL])
    MetaIndex[1:] = MetaIndex[1:]/factor[:-1]
    # print(MetaIndex)

    MetaFilt = []
    for iL in range(nMetaLevel-1):
        MetaFilt.append(MetaUnique[iL][MetaIndex[iL]])
    # print(MetaFilt)
    MetaCode_2[tuple(MetaIndex)] = MetaFilt
    
    Filter = numpy.all(MetaCode[:,:-1] == MetaFilt,axis=1) & (~numpy.isnan(KsMean))
    nS = Filter.sum()
    print('%s : %s Samples'%(' | '.join(MetaFilt),nS))
                
    Ktemp=KsMean[Filter]
    nK_BeforeFilter[tuple(MetaIndex)] = Ktemp.size
    Ktemp=RemoveOutlier(Ktemp)
    nK_AfterFilter[tuple(MetaIndex)] = Ktemp[~numpy.isnan(Ktemp)].size
    n=Ktemp.size
    
    if n == 1: 
        Ks_geomean[tuple(MetaIndex)] = Ktemp
        Ks_geostd[tuple(MetaIndex)] = 0.0
        Ks_geomean_F[tuple(MetaIndex)] = Ktemp
        Ks_geostd_F[tuple(MetaIndex)] = 0.0
        Ks_geomean_W[tuple(MetaIndex)] = Ktemp
        Ks_geostd_W[tuple(MetaIndex)] = 0  
        Ks_Min[tuple(MetaIndex)] = Ktemp
        Ks_Max[tuple(MetaIndex)] = Ktemp
    elif n > 0: 
        Ks_geomean[tuple(MetaIndex)] = stats.gmean(Ktemp)
        Ks_geostd[tuple(MetaIndex)] = gstd(Ktemp)
        Ks_geomean_F[tuple(MetaIndex)] = gmean_Finney(Ktemp,n)
        Ks_geostd_F[tuple(MetaIndex)] = gstd_Finney(Ktemp,n)
        Ks_geomean_W[tuple(MetaIndex)] = gmean_Warrick(Ktemp,n)
        Ks_geostd_W[tuple(MetaIndex)] = gstd_Warrick(Ktemp,n)  
        Ks_Min[tuple(MetaIndex)] = numpy.nanmin(Ktemp)
        Ks_Max[tuple(MetaIndex)] = numpy.nanmax(Ktemp)


#==============================================================================
# Plot
#==============================================================================


for iF in range(nDS):
    
    fig=pylab.figure()
    ax=fig.add_subplot(111)
    ax.axhline(KsMean[iF],color = 'r',label = 'Mean')
    ax.plot(iMeas,Ks[iF,:],'o-',label = 'Measurement')
    ax.set_xticks(iMeas)
    ax.set_xticklabels(map(lambda x : '%s - %s'%(Measurement[x],Date[x]),iMeas))
    ax.set_xlabel('Measurement')
    ax.set_ylabel('Ks (m/s)')
    ax.set_ylim([0,None])
    ax.set_title('_'.join(Data_Meta[iF]))
    ax.text(0.95, 0.85,'%s'%Data_Com[iF,0],
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes)
    figlist.append(fig)

    # pylab.close()


fig=pylab.figure()
ax=fig.add_subplot(111)
ax.hist(Ks[~numpy.isnan(Ks)],bins=20)
ax.set_xlabel('Ks (m/s)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution Ks')
figlist.append(fig)

bins=20
# bins=numpy.linspace(-8,-2,20)

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.hist(numpy.log10(Ks[~numpy.isnan(Ks)]),bins=bins)
ax.set_xlabel('log10(Ks) (m/s)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution Ks')
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.hist(numpy.log10(Ks[~numpy.isnan(Ks)]*100.0*3600*24),bins=20)
ax.set_xlabel('log10(Ks) (cm/d)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution Ks')
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.hist(numpy.log10(Ks_geomean[~numpy.isnan(Ks_geomean)]),bins=bins)
ax.set_xlabel('log10(Ks) (m/s)')
ax.set_ylabel('Frequency')
ax.set_title('Geometric Mean')
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.hist(numpy.log10(Ks_geomean_F[~numpy.isnan(Ks_geomean_F)]),bins=bins)
ax.set_xlabel('log10(Ks) (m/s)')
ax.set_ylabel('Frequency')
ax.set_title('Geometric mean Ks - Finney')
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.hist(numpy.log10(Ks_geomean_W[~numpy.isnan(Ks_geomean_W)]),bins=bins)
ax.set_xlabel('log10(Ks) (m/s)')
ax.set_ylabel('Frequency')
ax.set_title('Geometric mean Ks - Warrick')
figlist.append(fig)


#==============================================================================
# Write
#==============================================================================

# WriteRoot = os.path.dirname(MainRoot)
# WriteRoot = MainRoot
WriteRoot = os.path.join(MainRoot,'Output')


FinalResults1=numpy.concatenate((ID_IDData[iMatchID,:].astype('str'),
                                 ID_MetaData[iMatchID,:],
                                 ID_CommentData[iMatchID][:,None],
                                 MetaCode,
                                 Ks, Ks*24*3600*100,
                                 KsMean[:,None], KsMean[:,None]*24*3600*100,
                                 numpy.ones(nDS)[:,None]*-1.0, numpy.log10(KsMean*24*3600*100)[:,None], numpy.ones(nDS)[:,None]
                                 ,),axis=1)

Header1L=[]
Header1L = Header1L + list(ID_IDHeader)
Header1L = Header1L + list(ID_MetaHeader)
Header1L.append((Header2[iH_Com]))
# Header1L.append('FileName')
Header1L = Header1L + list(ID_MetaHeader)
Header1L = Header1L + ['Ks (m/s)', 'Ks (cm/d)']
Header1L = Header1L + ['Ks Mean (m/s)', 'Ks Mean (cm/d)']
Header1L = Header1L + ['Hyprop pF', 'Hyprop log10(Ks) (cm/d)' ,'Hyprop Weight']

Header1 = '\t'.join(Header1L)
fmt1='%s'

resultfile1='%s_%s_1.txt'%(Project,ThisAnalysisType)
resultname1=os.path.join(WriteRoot,resultfile1)
numpy.savetxt(resultname1,FinalResults1, fmt= fmt1, delimiter='\t', newline='\n', header= Header1)



[Ks_geomean,Ks_geostd,Ks_geomean_F,Ks_geostd_F,Ks_geomean_W,Ks_geostd_W,Ks_Min,Ks_Max,nK_BeforeFilter,nK_AfterFilter] = map(lambda x : x.flatten()[:,None],[Ks_geomean,Ks_geostd,Ks_geomean_F,Ks_geostd_F,Ks_geomean_W,Ks_geostd_W,Ks_Min,Ks_Max,nK_BeforeFilter,nK_AfterFilter])
MetaCode_2 = numpy.reshape(MetaCode_2,[-1,nMetaLevel-1])
FinalResults2=numpy.concatenate((MetaCode_2,
                                 Ks_geomean,Ks_geostd,
                                 Ks_geomean_F,Ks_geostd_F,
                                 Ks_geomean_W,Ks_geostd_W,
                                 Ks_Min, Ks_Max,
                                 Ks_geomean*24*3600*100, Ks_geostd*24*3600*100,
                                 Ks_geomean_F*24*3600*100,Ks_geostd_F*24*3600*100,
                                 Ks_geomean_W*24*3600*100,Ks_geostd_W*24*3600*100,
                                 Ks_Min*24*3600*100, Ks_Max*24*3600*100,
                                 nK_BeforeFilter,nK_AfterFilter
                                 ,),axis=1)

Header2L=[]
Header2L = Header2L + list(ID_MetaHeader[:-1])
Header2L = Header2L + ['Ks Geomean (m/s)', 'Ks Geostd (m/s)']
Header2L = Header2L + ['Ks Finney (m/s)', 'Ks Finney std (m/s)']
Header2L = Header2L + ['Ks Warrick (m/s)', 'Ks Warrick std (m/s)']
Header2L = Header2L + ['Ks Min (m/s)', 'Ks Max (m/s)']
Header2L = Header2L + ['Ks Geomean (cm/d)', 'Ks Geostd (cm/d)']
Header2L = Header2L + ['Ks Finney (cm/d)', 'Ks Finney std (cm/d)']
Header2L = Header2L + ['Ks Warrick (cm/d)', 'Ks Warrick std (cm/d)']
Header2L = Header2L + ['Ks Min (cm/d)', 'Ks Max (cm/d)']
Header2L = Header2L + ['nSamples before filter', 'nSamples after filter']

Header2 = '\t'.join(Header2L)
fmt2='%s'

resultfile2='%s_%s_2.txt'%(Project,ThisAnalysisType)
resultname2=os.path.join(WriteRoot,resultfile2)
numpy.savetxt(resultname2,FinalResults2, fmt= fmt2, delimiter='\t', newline='\n', header= Header2)

## Save Plots

import sys
# basename=os.path.basename(sys.argv[0])[:-3]
basename='%s_%s'%(Project,ThisAnalysisType)


## save plots to pdf
pfdname=os.path.join(WriteRoot, basename+'.pdf')
pp = PdfPages(pfdname)
for fig in figlist:
    pp.savefig(fig)
pp.close()

# extension='.png'
# for iF in range(len(figlist)):
#     fig=figlist[iF]
#     figname = os.path.join(root, basename+'_%s'%iF + extension)
#     fig.savefig(figname, bbox_inches='tight')

# pylab.show()
