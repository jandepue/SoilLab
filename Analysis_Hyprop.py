#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
## ------------------------------------------------------------------------- ##
#                            READ HYPROP DATA                                 #
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
ThisAnalysisType = 'UMS_Hyprop'


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


#ID_FullMeta = list(map(lambda x : '_'.join(x),ID_MetaData))
ID_FullMeta = ['_'.join(x) for x in ID_MetaData]

#==============================================================================
# Find Hyprop data
#==============================================================================

DataRoot = os.path.join(MainRoot,'Hyprop')

FileList = []
for root, dirnames, filenames in os.walk(DataRoot):
    for filename in fnmatch.filter(filenames, '*.csv'):
        FileList.append(os.path.join(root, filename))

FileList.sort()
FileList=numpy.array(FileList)
#FileNames=numpy.array(map(lambda x : os.path.basename(x)[:-4].upper(),FileList)) 
FileNames=numpy.array([os.path.basename(x)[:-4].upper() for x in FileList])


## FILTER 
# FiltRef=numpy.char.find(FileList, 'Originals') > -1
# FileNames=FileNames[~FiltRef] 
# FileList=FileList[~FiltRef]

nFiles=len(FileList)
print(root+" : %s files found." %nFiles) # 

# nMetaLevel = len(FileNames[0].split('_'))
nMetaLevel = ID_MetaData.shape[1]


MetaCode = numpy.zeros((nFiles,nMetaLevel),dtype='S256')

for iF in range(nFiles):
    Meta = FileNames[iF].split('_')
    if len(Meta) != nMetaLevel:
        raise ValueError('file %s has a wrong name'%FileNames[iF])
    MetaCode[iF,:] = Meta

    
MetaUnique=[]
for iM in range(nMetaLevel):
    MetaUnique.append(numpy.unique(MetaCode[:,iM]))

#nMetaUnique = map(numpy.size,MetaUnique)
nMetaUnique = [x.size for x in MetaUnique]

#==============================================================================
# Check Match ID
#==============================================================================

iMatchID = numpy.zeros(nFiles,dtype ='int')
iMatchID_R = numpy.zeros(nSamp,dtype ='int')+numpy.nan

for iF in range(nFiles):
    Meta = MetaCode[iF,:]
    
    # FindMatch = ID_MetaData == Meta  ## Original
    FindMatch = ID_MetaData == numpy.tile(Meta[None,:],(nFiles,1))
    
    MatchFound = numpy.all(FindMatch, axis=1)
    if not numpy.any(MatchFound):
        raise ValueError('No match found in Identification Database for %s'%'_'.join(Meta))
    else:
        iMatchID[iF] = numpy.where(MatchFound)[0]
        iMatchID_R[MatchFound] = iF

print('%s / %s samples found'%(nFiles,nSamp))


#==============================================================================
# Open Data
#==============================================================================


BulkDens_List=numpy.zeros(nFiles)
ParHyprop_SWR_List=numpy.zeros((nFiles,7))
ParHyprop_Kh_List=numpy.zeros((nFiles,10))
AC=numpy.zeros(nFiles)
MP=numpy.zeros(nFiles)
AC2=numpy.zeros(nFiles)
MP2=numpy.zeros(nFiles)
Por=numpy.zeros(nFiles)
AC_P=numpy.zeros(nFiles)
MP_P=numpy.zeros(nFiles)
AC2_P=numpy.zeros(nFiles)
MP2_P=numpy.zeros(nFiles)
Data4_pF_All = numpy.zeros((nFiles,100))+numpy.nan
Data4_WC_All = numpy.zeros((nFiles,100))+numpy.nan
Data5_pF_All = numpy.zeros((nFiles,100))+numpy.nan
Data5_Kh_All = numpy.zeros((nFiles,100))+numpy.nan


# TensionPoints=-numpy.logspace(-1,5,25)
TensionPoints=-numpy.logspace(0,4,25)
SWCPoints=numpy.zeros((nFiles,TensionPoints.size))
KhPoints=numpy.zeros((nFiles,TensionPoints.size))
SupplementaryFitModel = Vangenuchten5
SupplementaryFitPar=numpy.zeros((nFiles,5))


for iF in range(nFiles):
    filename = FileList[iF]
    label=os.path.basename(filename).replace('.csv','')


    ## RAW DATA (for metadata, bulk density, etc)
    
    Sample,SV,DSW,Data1_Time,Data1_Time_Sec,Data1_TensionB,Data1_TensionT,Data1_Temp,Data2_Time,Data2_Time_Sec,Data2_MassNet=ReadHypropCSV_V2_RawData(filename)                       
    # Sample,SV,DSW,Data1_Time,Data1_Time_Sec,Data1_TensionB,Data1_TensionT,Data1_Temp,Data2_Time,Data2_Time_Sec,Data2_MassNet=ReadHypropCSV_V3_RawData(filename)                       
    
    ParticleDensity = 2.65e3
                        
    rhow=1.0 # g/cm3
    rhob=DSW/SV # g/cm3
    BulkDens_List[iF]=rhob*1e3
    
    Porosity=1-rhob*1e3/ParticleDensity
    Por[iF]=Porosity


    
    # ## SAME CALCULATION AS HYPROP (CAN BE SKIPPED)
    
    # AElimit=-3e-7 # hpa/s2
    # pAirEntry=8800.0 # hpa
    # nForInterpol=200
    # # nForInterpol=20
                        
    # Data1_Time_Sec_TOT,TensionMean_TOT,Data_TensionB_TOT,Data_TensionT_TOT,IStopB,IStopT,IAirEntryB,IAirEntryT = AirEntryExtrapolation(Data1_Time_Sec,Data1_TensionB,Data1_TensionT,AElimit=AElimit,pAirEntry=pAirEntry,nForInterpol=nForInterpol)
    # # [Data1_Time_Sec_TOT,TensionMean_TOT,Data_TensionB_TOT,Data_TensionT_TOT]=map(lambda x : x[Data1_Time_Sec_TOT<=Data2_Time_Sec.max()],[Data1_Time_Sec_TOT,TensionMean_TOT,Data_TensionB_TOT,Data_TensionT_TOT])


    # ## Tensiometer 
    # TensionMean_TOT_cm=TensionMean_TOT * 1e2*1e2/(1e3*9.81) # cm
    

    # ## Water Content spline/fit

    # ## Spline
    # # MassInterp=interpolate.interp1d(Data2_Time_Sec,Data2_MassNet,kind='cubic')
    # # MassInterp=interpolate.interp1d(Data2_Time_Sec,Data2_MassNet,kind='quadratic')
    # # Mass_TOT=MassInterp(Data1_Time_Sec_TOT)
                       
    # ## Polyfit
    # degree=5
    # MassInterp=numpy.polyfit(Data2_Time_Sec,Data2_MassNet,degree)
    # Mass_TOT=numpy.polyval(MassInterp,Data1_Time_Sec_TOT)
    
    # MassWater=Mass_TOT-DSW # g
    # SWC_g=MassWater/DSW # g/g
    # SWC_v_TOT=SWC_g*rhob/rhow # cm3/cm3


    
    ## HYPROPS INTERPOLATED WATER CONTENT, PF, KH
    Sample,Data4_pF,Data4_WC,Data5_pF,Data5_Kh,Data6_WC,Data6_Kh=ReadHypropCSV_V2_NiceData(filename)
    # Sample,Data4_pF,Data4_WC,Data5_pF,Data5_Kh,Data6_WC,Data6_Kh=ReadHypropCSV_V3_NiceData(filename)
    Data4_pF_All[iF,:Data4_pF.size] = Data4_pF
    Data4_WC_All[iF,:Data4_WC.size] = Data4_WC
    Data5_pF_All[iF,:Data5_pF.size] = Data5_pF
    Data5_Kh_All[iF,:Data5_Kh.size] = Data5_Kh


    
    # ## DO YOUR OWN FIT (CAN BE SKIPPED)
    # PythonModel_SWR=h2theta_VanGenuchten5
    # Par0=(SWC_v_TOT[-1]/2,SWC_v_TOT[0],0.01,1.5,0.5)
    # bnd=((0,1),(0,1),(0,None),(0,None),(0,None))
    # PythonModel_SWR=h2theta_VanGenuchten4
    # Par0=(SWC_v_TOT[-1]/2,SWC_v_TOT[0],0.01,5.0)
    # bnd=((0,1),(0,1),(0,None),(0,None))
    # PythonModel_SWR=h2theta_Durner
    # Par0_SWR=(SWC_v_TOT[-1]/2,SWC_v_TOT[0],0.5,0.01,1.5,0.01,1.5)
    # bnd=((0,1),(0,1),(0,1),(0,None),(1,None),(0,None),(1,None))

    # Min=optimize.minimize(FitRetentionModel, Par0_SWR, args=(TensionMean_TOT_cm,SWC_v_TOT,PythonModel_SWR),method='L-BFGS-B',bounds=bnd)
    # Min=optimize.minimize(FitRetentionModel, Par0_SWR, args=(TensionMean_TOT_cm,SWC_v_TOT,PythonModel_SWR),method='SLSQP',bounds=bnd)
    # ParFit_SWR=Min.x

    # PythonModel_Kh=h2K_BimodalPDIK
    # Par0_Kh=(0.0,0.5,10**numpy.mean(Data5_Kh[Data5_pF==-1]),2.0,0.5,0.5,0.01,1.5,0.01,1.5)
    # Par0_Kh=(ParFit_SWR[0],ParFit_SWR[1],10**numpy.mean(Data5_Kh[Data5_pF==-1]),2.0,0.5,0.5,ParFit_SWR[3],ParFit_SWR[4],ParFit_SWR[5],ParFit_SWR[6])
    # bnd=((0,1),(0,1),(0,None),(0,None),(0,1),(0,1),(0,None),(1,None),(0,None),(1,None))

    # Min=optimize.minimize(FitRetentionModel, Par0_Kh, args=(10**Data5_pF,10**Data5_Kh,PythonModel_Kh),method='L-BFGS-B',bounds=bnd)
    # Min=optimize.minimize(FitRetentionModel2, Par0_Kh, args=(10**Data5_pF,10**Data5_Kh,PythonModel_Kh),method='SLSQP',bounds=bnd)
    # ParFit_Kh=Min.x
    
    

    ## HYPROP FIT
    Sample,Data7_isfit,Data7_Par,Data7_ParVal,Data7_Model,NewVersion=ReadHypropCSV_V2_Fit(filename)
    # Sample,Data7_isfit,Data7_Par,Data7_ParVal,Data7_Model,NewVersion=ReadHypropCSV_V3_Fit(filename)

    Ks=Data7_ParVal[Data7_Par == 'Ks']

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
        print('Ks Temperature at 10 °C - NO CORRECTION')


    ## !!! DOUBLE CHECK THE MODEL AND PARAMETERS HERE !!! ##
    
    if NewVersion==0:
        HypropModel_SWR=Durner
        HypropModel_Kh=PetersDurnerIBimodal
        ParHyprop_SWR=(Data7_ParVal[2],Data7_ParVal[3],1.0-Data7_ParVal[6],Data7_ParVal[0],Data7_ParVal[1],Data7_ParVal[4],Data7_ParVal[5]) # thetar,thetas,w1,alpha1,n1,alpha2,n2=Par
        ParHyprop_Kh=(Ks,Data7_ParVal[7],1.0-Data7_ParVal[6],Data7_ParVal[9],Data7_ParVal[0],Data7_ParVal[1],Data7_ParVal[4],Data7_ParVal[5])

    else:
        HypropModel_SWR=h2theta_BimodalPDIRet
        HypropModel_Kh=h2K_BimodalPDIK 
        ParHyprop_SWR=(Data7_ParVal[2],Data7_ParVal[3],1.0-Data7_ParVal[6],Data7_ParVal[0],Data7_ParVal[1],Data7_ParVal[4],Data7_ParVal[5])
        ParHyprop_Kh=(Data7_ParVal[2],Data7_ParVal[3],Ks,Data7_ParVal[7],1.0-Data7_ParVal[6],Data7_ParVal[9],Data7_ParVal[0],Data7_ParVal[1],Data7_ParVal[4],Data7_ParVal[5]) # thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2=Par

    ParHyprop_SWR_List[iF,:]=ParHyprop_SWR
    ParHyprop_Kh_List[iF,:]=ParHyprop_Kh


    ## ADDITIONAL CALCULATIONS

    ## Air Content

    AC[iF]=ParHyprop_SWR[1]-HypropModel_SWR(numpy.array([-50.0,]),ParHyprop_SWR)
    MP[iF]=ParHyprop_SWR[1]-HypropModel_SWR(numpy.array([-10.0,]),ParHyprop_SWR)
    AC2[iF]=ParHyprop_SWR[1]-HypropModel_SWR(numpy.array([-100.0,]),ParHyprop_SWR)
    MP2[iF]=ParHyprop_SWR[1]-HypropModel_SWR(numpy.array([-30.0,]),ParHyprop_SWR)

    AC_P[iF]=Porosity-HypropModel_SWR(numpy.array([-50.0,]),ParHyprop_SWR)
    MP_P[iF]=Porosity-HypropModel_SWR(numpy.array([-10.0,]),ParHyprop_SWR)
    AC2_P[iF]=Porosity-HypropModel_SWR(numpy.array([-100.0,]),ParHyprop_SWR)
    MP2_P[iF]=Porosity-HypropModel_SWR(numpy.array([-30.0,]),ParHyprop_SWR)
    
    ## Points

    SWCPoints[iF,:]=HypropModel_SWR(TensionPoints,ParHyprop_SWR)
    KhPoints[iF,:]=HypropModel_Kh(TensionPoints,ParHyprop_Kh)

    ## supplementary fit against Points

    Par0=(SWCPoints[iF,-1],SWCPoints[iF,0],0.01,1.5,0.5)
    bnd=((0,1),(0,1),(1e-5,None),(1.1,None),(0.05,None))
    Minvg=optimize.minimize(FitRetentionModel, Par0, args=(TensionPoints,SWCPoints[iF,:],SupplementaryFitModel),method='L-BFGS-B',bounds=bnd)
    ParVG=Minvg.x
    SupplementaryFitPar[iF,:]=ParVG

    

    ## PLOT

    # # Plot raw data
    # fig=pylab.figure()
    # fig.suptitle(label)
    # ax=fig.add_subplot(211)
    # ax.plot(Data1_Time_Sec_TOT,Data_TensionT_TOT,'--r',label = 'python')
    # ax.plot(Data1_Time_Sec_TOT,Data_TensionB_TOT,'--b',label = 'python')
    # ax.plot(Data1_Time_Sec_TOT,TensionMean_TOT,'-r')
    # ax.plot(Data1_Time_Sec,Data1_TensionT,'k',label = 'hyprop')
    # ax.plot(Data1_Time_Sec,Data1_TensionB,'b',label = 'hyprop')
    # ax.axvline(Data1_Time_Sec[IStopT],color='k')
    # ax.axvline(Data1_Time_Sec[IStopB],color='b')
    # ax.axvline(Data1_Time_Sec[IAirEntryT],color='k')
    # ax.axvline(Data1_Time_Sec[IAirEntryB],color='b')
    # ax.set_ylim([0,88e2])
    # ax.legend()
    # xlim=ax.get_xlim()
    # ax.set_xlabel('Time (s)')
    # ax.set_ylabel('Matric Potential (hPa)')
    # # ax.set_xlabel('Tijd (s)')
    # # ax.set_ylabel('Matric Potentiaal (hPa)')
    # # fig=pylab.figure()
    # ax=fig.add_subplot(212)
    # ax.plot(Data1_Time_Sec_TOT,Mass_TOT,'--k')
    # ax.plot(Data2_Time_Sec,Data2_MassNet,'or')
    # ax.set_xlabel('Time (s)')
    # ax.set_ylabel('Mass (g)')
    # # ax.set_xlabel('Tijd (s)')
    # # ax.set_ylabel('Massa (g)')
    # ax.set_xlim(xlim)
    # figlist.append(fig)

    # fig=pylab.figure()
    # ax=fig.add_subplot(111)
    # ax.plot(TensionMean_TOT,SWC_v_TOT,'.k',label='Python')
    # ax.plot(10**Data4_pF ,Data4_WC/100.0,'.r',label='Hyprop')
    # ax.set_xscale('log')
    # ax.set_xlabel('Matric Potential (hPa)')
    # ax.set_ylabel('Volumetric Water Content (m3/m3)')
    # ax.set_xlabel('Matric Potentiaal (hPa)')
    # ax.set_ylabel('Volumetrisch Vochtgehalte (m3/m3)')
    # ax.legend(loc=1)
    # ax.set_title(label)
    # figlist.append(fig)


    Tensionplot = -numpy.logspace(-2,5,1000)
    
    fig=pylab.figure()
    fig.suptitle(label)
    ax=fig.add_subplot(211)
    ax.plot(10**Data4_pF ,Data4_WC*1e-2,'.k')
    ax.plot(-Tensionplot ,HypropModel_SWR(Tensionplot,ParHyprop_SWR),'-',color='r',label=label)
    ax.set_xscale('log')
    ax.set_xlabel('Matric Potential (cm)')
    ax.set_ylabel('Water Content (m3/m3)')
    ax=fig.add_subplot(212)
    ax.plot(10**Data5_pF ,10**Data5_Kh,'.k')
    ax.plot(-Tensionplot ,HypropModel_Kh(Tensionplot,ParHyprop_Kh),'-',color='r',label=label)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Matric Potential (cm)')
    ax.set_ylabel('Hydraulic Conductivity (cm/d)')
    figlist.append(fig)

    fig=pylab.figure()
    fig.suptitle(label)
    ax=fig.add_subplot(111)
    ax.plot(-TensionPoints,SWCPoints[iF,:],'.k')
    ax.plot(-Tensionplot ,SupplementaryFitModel(Tensionplot,SupplementaryFitPar[iF,:]),'-',color='r')
    ax.set_xscale('log')
    ax.set_xlabel('Matric Potential (cm)')
    ax.set_ylabel('Water Content (m3/m3)')
    ax.set_title('Exported SWR Points + Supplementary Fit')
    figlist.append(fig)

    pylab.close('all')




## Assume last metalevel is repetition => Make plots of the repetitions

cmap = pylab.get_cmap('Paired')

MeanShape = nMetaUnique[:-1]
nMeta_Full = numpy.prod(MeanShape) ## total amount of classes

BD_Mean = numpy.zeros(MeanShape)+numpy.nan
Por_Mean = numpy.zeros(MeanShape)+numpy.nan
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
    
    Filter = numpy.all(MetaCode[:,:-1] == MetaFilt,axis=1) & (~numpy.isnan(Ks))
    nS = Filter.sum()
    print('%s : %s Samples'%(' | '.join(MetaFilt),nS))
    
    if nS > 0: 
        BD_Mean[tuple(MetaIndex)] = numpy.mean(BulkDens_List[Filter])
        Por_Mean[tuple(MetaIndex)] = numpy.mean(Por[Filter])

        fig=pylab.figure()
        fig.suptitle('_'.join(MetaFilt))
        ax=fig.add_subplot(211)
        for iS in range(nS):
            color = cmap(iS/(nS-0.99))
            ax.plot(10**Data4_pF_All[Filter,:][iS,:] ,Data4_WC_All[Filter,:][iS,:]*1e-2,'.',color = color)
            ax.plot(-Tensionplot ,HypropModel_SWR(Tensionplot,ParHyprop_SWR_List[Filter,:][iS,:]),'-',color = color,label=MetaCode[Filter,-1][iS])
        ax.set_xscale('log')
        ax.set_xlabel('Matric Potential (cm)')
        ax.set_ylabel('Water Content (m3/m3)')
        ax.legend()
        ax=fig.add_subplot(212)
        for iS in range(nS):
            color = cmap(iS/(nS-0.99))
            ax.plot(10**Data5_pF_All[Filter,:][iS,:] ,10**Data5_Kh_All[Filter,:][iS,:],'.',color = color)
            ax.plot(-Tensionplot ,HypropModel_Kh(Tensionplot,ParHyprop_Kh_List[Filter,:][iS,:]),'-',color = color,label=MetaCode[Filter,-1][iS])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Matric Potential (cm)')
        ax.set_ylabel('Hydraulic Conductivity (cm/d)')
        figlist.append(fig)


#==============================================================================
# Plot
#==============================================================================

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.hist(BulkDens_List[~numpy.isnan(BulkDens_List)],bins=20)
ax.set_xlabel('Bulk Density (kg/m3)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution BD')
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
                                 FileNames[:,None], MetaCode,
                                 ParHyprop_SWR_List,
                                 ParHyprop_Kh_List,
                                 BulkDens_List[:,None], Por[:,None],
                                 AC[:,None], AC2[:,None], MP[:,None], MP2[:,None],
                                 AC_P[:,None], AC2_P[:,None], MP_P[:,None], MP2_P[:,None],
                                 SWCPoints,KhPoints,SupplementaryFitPar
                                 ,),axis=1)

Header1L=[]
Header1L = Header1L + list(ID_IDHeader)
Header1L = Header1L + list(ID_MetaHeader)
Header1L.append((Header2[iH_Com]))
Header1L.append('FileName')
Header1L = Header1L + list(ID_MetaHeader)
#Header1L = Header1L + map(lambda x : 'HypropSWR_Par%s'%x,range(ParHyprop_SWR_List.shape[1]))
Header1L = Header1L + ['HypropSWR_Par%s'%x for x in range(ParHyprop_SWR_List.shape[1])]
#Header1L = Header1L + map(lambda x : 'HypropKh_Par%s'%x,range(ParHyprop_Kh_List.shape[1]))
Header1L = Header1L + ['HypropKh_Par%s'%x for x in range(ParHyprop_Kh_List.shape[1])]
Header1L = Header1L + ['Bulk Density (kg/m3)','Porosity (m3/m3)',]
Header1L = Header1L + ['Air Capacity (50cm,thetas) (m3/m3)','Air Capacity (100cm,thetas) (m3/m3)','Macro Porosity (10cm,thetas) (m3/m3)','Macro Porosity (30cm,thetas) (m3/m3)',]
Header1L = Header1L + ['Air Capacity (50cm,porosity) (m3/m3)','Air Capacity (100cm,porosity) (m3/m3)','Macro Porosity (10cm,porosity) (m3/m3)','Macro Porosity (30cm,porosity) (m3/m3)',]
#Header1L = Header1L + map(lambda x : 'SWRCPoint_pF%1.2f'%x,numpy.log10(-TensionPoints))
Header1L = Header1L + ['SWRCPoint_pF%1.2f'%x for x in numpy.log10(-TensionPoints)]
#Header1L = Header1L + map(lambda x : 'KhPoint_pF%1.2f'%x,numpy.log10(-TensionPoints))
Header1L = Header1L + ['KhPoint_pF%1.2f'%x for x in numpy.log10(-TensionPoints)]
#Header1L = Header1L + map(lambda x : 'SupplementaryFit_Par%s'%x,range(SupplementaryFitPar.shape[1]))
Header1L = Header1L + ['SupplementaryFit_Par%s'%x for x in range(SupplementaryFitPar.shape[1])]

Header1 = '\t'.join(Header1L)
fmt1='%s'

resultfile1='%s_%s_1.txt'%(Project,ThisAnalysisType)
resultname1=os.path.join(WriteRoot,resultfile1)
numpy.savetxt(resultname1,FinalResults1, fmt= fmt1, delimiter='\t', newline='\n', header= Header1)



#[BD_Mean,Por_Mean] = map(lambda x : x.flatten()[:,None],[BD_Mean,Por_Mean])
[BD_Mean,Por_Mean] = [x.flatten()[:,None] for x in [BD_Mean,Por_Mean]]
MetaCode_2 = numpy.reshape(MetaCode_2,[-1,nMetaLevel-1])
FinalResults2=numpy.concatenate((MetaCode_2,
                                 BD_Mean,
                                 Por_Mean,),axis=1)

Header2L=[]
Header2L = Header2L + list(ID_MetaHeader[:-1])
Header2L = Header2L + ['Bulk Density (kg/m3)','Porosity (m3/m3)']

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
