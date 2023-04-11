#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
## ------------------------------------------------------------------------- ##
#                     READ SANDBOX & PRESSUREPLATE DATA                       #
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
ThisAnalysisType = 'SBPP'

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

DataRoot = os.path.join(MainRoot,'SBPP')


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
skip_header = 6


fID = open(filename)
Project = fID.readline().split(delimiter)[1]
Date = fID.readline().split(delimiter)[1]
Collaborators = fID.readline().split(delimiter)[1:]
skipline = fID.readline()
Header1 = numpy.array(fID.readline().replace('\n','').split(delimiter))
Header2 = numpy.array(fID.readline().replace('\n','').split(delimiter))
fID.close()

iH_ID_Data = numpy.where(Header1 == 'Identification')[0][0]
iH_Meta_Data = numpy.where(Header1 == 'Metadata')[0][0]
iH_Com_Data = numpy.where(Header1 == 'Comments')[0][0]
iH_Data = numpy.where(Header1 == 'Data')[0][0]

Data_All = numpy.genfromtxt(filename,delimiter = delimiter, dtype ='str', skip_header = skip_header)
Data_All = CapitalizeArray(Data_All)

Data_ID = Data_All[:,iH_ID_Data:iH_Meta_Data]
Data_Meta = Data_All[:,iH_Meta_Data:iH_Com_Data]
Data_Com = Data_All[:,iH_Com_Data:iH_Data]
Data = Data_All[:,iH_Data:].astype('float')
Data[Data_All[:,iH_Data:] == ''] = 0 # replace empty cells with zero (but not nan)
Data_Header = Header2[iH_Data:]

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

Height = Data[:,Data_Header == 'RingHeight(m)']       # m
Diameter = Data[:,Data_Header == 'RingDiameter(m)']   # m

HeightCorr = Data[:,Data_Header == 'HeightCorrection(m)']   # m
VolumeCorr = Data[:,Data_Header == 'VolumeCorrection(m3)']   # m
MassCorr = Data[:,Data_Header == 'MassCorrection(kg)']   # m
ParrafinMassCorr = Data[:,Data_Header == 'ParrafineCorrection(kg)']   # m

TerraMass = Data[:,Data_Header == 'EmptyRingWeight(kg)'] + Data[:,Data_Header == 'NylonWeight(kg)']   # kg

filt = numpy.array(map(lambda x : 'MassSBStable' in x,Data_Header))
StableMass_SB = Data[:,filt]
Head_SB = -numpy.array(map(lambda x : float(x.replace('MassSBStable','').replace('(kg)','')),Data_Header[filt]))

MassSB_BeforeSS = Data[:,Data_Header == 'MassSBBeforeSubsample(kg)']   # kg
MassCup_SSRef = Data[:,Data_Header == 'MassSBSubsampleCup(kg)']   # kg
MassWet_SSRef = Data[:,Data_Header == 'MassSBSubsampleWet(kg)']   # kg
MassDry_SSRef = Data[:,Data_Header == 'MassSBSubsampleDry(kg)']   # kg

filt = numpy.array(map(lambda x : 'MassPPSubsample' in x,Data_Header)) & numpy.array(map(lambda x : 'Cup' in x,Data_Header))
MassCup_SS = Data[:,filt]

filt = numpy.array(map(lambda x : 'MassPPSubsample' in x,Data_Header)) & numpy.array(map(lambda x : 'Wet' in x,Data_Header))
MassWet_SS = Data[:,filt]

filt = numpy.array(map(lambda x : 'MassPPSubsample' in x,Data_Header)) & numpy.array(map(lambda x : 'Dry' in x,Data_Header))
MassDry_SS = Data[:,filt]

Head_SS = -numpy.array(map(lambda x : float(x.replace('MassPPSubsample','').replace('Dry(kg)','')),Data_Header[filt]))


## Corrections

ParrafinDensity = 912.2 # kg/m3 (W. F. Seter and K. Iuouye, 1935)
ParrafinVolCorr = ParrafinMassCorr / ParrafinDensity

RingVolume = (Height+HeightCorr) * numpy.pi*(Diameter/2)**2 - ParrafinVolCorr - VolumeCorr

TerraMass = TerraMass - MassCorr


## Calculate Water contents & bulk density

WC_g_SSRef = (MassWet_SSRef - MassDry_SSRef) / (MassDry_SSRef - MassCup_SSRef)

DryMass_SB = (MassSB_BeforeSS-TerraMass)/(1+WC_g_SSRef)
WetMass_SB = StableMass_SB - TerraMass

WC_g_SB = (WetMass_SB - DryMass_SB)/DryMass_SB
WC_g_SS = (MassWet_SS - MassDry_SS) / (MassDry_SS - MassCup_SS)

WC_g = numpy.concatenate((WC_g_SB,WC_g_SS),axis=1)
Head = numpy.concatenate((Head_SB,Head_SS))
nH = Head.size

BulkDensity = DryMass_SB/RingVolume
ParticleDensity = 2.65e3
Porosity = 1 - BulkDensity/ParticleDensity

WaterDensity = 1e3

WC_v = WC_g * BulkDensity/WaterDensity


fig = pylab.figure()
ax = fig.add_subplot(111)
ax.boxplot(WC_v)
ax.set_xticks(range(1,nH+1))
ax.set_xticklabels(map(lambda x : '%1.2e'%x,Head))
ax.set_xlabel('Matric Head (cm)')
figlist.append(fig)


fig=pylab.figure()
ax=fig.add_subplot(111)
ax.hist(BulkDensity,bins=20)
ax.set_xlabel('Bulk Density (kg/m3)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution BD')
figlist.append(fig)



## Fit

SWRCModel = Vangenuchten5
SWRCPar=numpy.zeros((nDS,5))

TensionPoints=-numpy.logspace(0,4,25)
SWCPoints=numpy.zeros((nDS,TensionPoints.size))

for iF in range(nDS):
    
    WC = WC_v[iF,:]
    
    Par0=(WC[-1],WC[0],0.01,1.5,0.5)
    bnd=((0,1),(0,1),(1e-5,None),(1.1,None),(0.05,None))
    Min=optimize.minimize(FitRetentionModel, Par0, args=(Head,WC,SWRCModel),method='L-BFGS-B',bounds=bnd)

    SWRCPar[iF,:]=Min.x
    SWCPoints[iF,:] = SWRCModel(TensionPoints,Min.x)
    
    Tensionplot = -numpy.logspace(-2,5,1000)
    
    fig=pylab.figure()
    ax=fig.add_subplot(111)
    ax.plot(-Head,WC,'.k')
    ax.plot(-Tensionplot ,SWRCModel(Tensionplot,SWRCPar[iF,:]),'-',color='r')
    ax.set_xscale('log')
    ax.set_xlabel('Matric Potential (cm)')
    ax.set_ylabel('Water Content (m3/m3)')
    ax.set_title('_'.join(Data_Meta[iF]))
    ax.text(0.95, 0.85,'%s'%Data_Com[iF,0],
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes)
    figlist.append(fig)

    # pylab.close()



## Assume last metalevel is repetition: average and grouped fit

cmap = pylab.get_cmap('Paired')

MeanShape = nMetaUnique[:-1]
nMeta_Full = numpy.prod(MeanShape) ## total amount of classes

BD_Mean = numpy.zeros(MeanShape)+numpy.nan
Por_Mean = numpy.zeros(MeanShape)+numpy.nan
WC_Mean = numpy.zeros(MeanShape + [nH,]) + numpy.nan
Par_Mean = numpy.zeros(MeanShape + [5,]) + numpy.nan
SWC_Mean = numpy.zeros(MeanShape + [TensionPoints.size,]) + numpy.nan
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
    
    Filter = numpy.all(MetaCode[:,:-1] == MetaFilt,axis=1)
    nS = Filter.sum()
    print('%s : %s Samples'%(' | '.join(MetaFilt),nS))

    if nS > 0:
        
        BD_Mean[tuple(MetaIndex)] = numpy.mean(BulkDensity[Filter])
        Por_Mean[tuple(MetaIndex)] = numpy.mean(Porosity[Filter])
        WC_Mean[tuple(MetaIndex)] = WC_v[Filter,:].mean(axis=0)


        ## FIT
        
        WC = WC_v[Filter,:]

        # WCFlat = WC.flatten()
        # HeadFlat = (numpy.ones((nS,nH))*Head[None,:]).flatten()
                               
        Par0=(WC[0,-1],WC[0,0],0.01,1.5,0.5)
        bnd=((0,1),(0,1),(1e-5,None),(1.1,None),(0.05,None))
        Min=optimize.minimize(FitRetentionModel, Par0, args=(Head[None,:],WC,SWRCModel),method='L-BFGS-B',bounds=bnd)
        MeanPar = Min.x
        
        Par_Mean[tuple(MetaIndex)] = MeanPar
        SWC_Mean[tuple(MetaIndex)] = SWRCModel(TensionPoints,MeanPar)
        
        
        Tensionplot = -numpy.logspace(-2,5,1000)

        fig=pylab.figure()
        ax=fig.add_subplot(111)
        for iS in range(nS):
            color = cmap(iS/(nS-0.99))
            ax.plot(-Head , WC[iS,:],'o',color = color,label=MetaCode[Filter,-1][iS])
        ax.plot(-Tensionplot ,SWRCModel(Tensionplot,MeanPar),'-',color='r')
        ax.set_xscale('log')
        ax.set_xlabel('Matric Potential (cm)')
        ax.set_ylabel('Water Content (m3/m3)')
        ax.set_title('_'.join(MetaFilt))
        ax.legend(loc=3,numpoints=1)
        figlist.append(fig)


# pylab.close()


#==============================================================================
# Plot
#==============================================================================



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
                                 SWRCPar,
                                 WC_v,
                                 BulkDensity,
                                 Porosity,
                                 SWCPoints,),axis=1)

Header1L=[]
Header1L = Header1L + list(ID_IDHeader)
Header1L = Header1L + list(ID_MetaHeader)
Header1L.append((Header2[iH_Com]))
# Header1L.append('FileName')
Header1L = Header1L + list(ID_MetaHeader)
Header1L = Header1L + map(lambda x : 'SWR_Par%s'%x,range(SWRCPar.shape[1]))
Header1L = Header1L + map(lambda x : 'WClab_pF%1.2f'%x,numpy.log10(-Head))
Header1L = Header1L + ['Bulk Density (kg/m3)','Porosity (m3/m3)',]
Header1L = Header1L + map(lambda x : 'SWRCPoint_pF%1.2f'%x,numpy.log10(-TensionPoints))

Header1 = '\t'.join(Header1L)
fmt1='%s'

resultfile1='%s_%s_1.txt'%(Project,ThisAnalysisType)
resultname1=os.path.join(WriteRoot,resultfile1)
numpy.savetxt(resultname1,FinalResults1, fmt= fmt1, delimiter='\t', newline='\n', header= Header1)




[BD_Mean,Por_Mean] = map(lambda x : x.flatten()[:,None],[BD_Mean,Por_Mean])
MetaCode_2 = numpy.reshape(MetaCode_2,[-1,nMetaLevel-1])
FinalResults2=numpy.concatenate((MetaCode_2,
                                 Par_Mean.reshape([-1,Par_Mean.shape[-1]]),
                                 WC_Mean.reshape([-1,WC_Mean.shape[-1]]),
                                 BD_Mean,
                                 Por_Mean,
                                 SWC_Mean.reshape([-1,SWC_Mean.shape[-1]]),),axis=1)

Header2L=[]
Header2L = Header2L + list(ID_MetaHeader[:-1])
Header2L = Header2L + map(lambda x : 'SWR_Par%s'%x,range(SWRCPar.shape[1]))
Header2L = Header2L + map(lambda x : 'WClab_pF%1.2f'%x,numpy.log10(-Head))
Header2L = Header2L + ['Bulk Density (kg/m3)','Porosity (m3/m3)']
Header2L = Header2L + map(lambda x : 'SWRCPoint_pF%1.2f'%x,numpy.log10(-TensionPoints))

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
