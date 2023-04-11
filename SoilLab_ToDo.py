#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
## ------------------------------------------------------------------------- ##
#                        WHAT HAS BEEN DONE AND WHAT NOT?                     #
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

from SoilLab_Funk import *


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


#MainRoot = '/home/jan/git/SoilLab/TestData'
MainRoot = 'C:\\Users\\admin\\Documents\\GitHub\\SoilLab\\TestData'

#==============================================================================
# Find Identification Data
#==============================================================================

IDFile = os.path.join(MainRoot,'Identification.csv')

delimiter = '\t'
skip_header = 6

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

ID_AllData = numpy.genfromtxt(IDFile,delimiter = delimiter, dtype ='str', skip_header = skip_header)
Analysis_AllData = ID_AllData[:,iH_Ana:iH_Com]=='1'
iAnalysis_Requested = Analysis_AllData.sum(axis=0)>0
Analysis_Requested = Header2[iH_Ana:iH_Com][iAnalysis_Requested]

print("Analysis Required: %s\n\t%s"%(Analysis_Requested.size,"\n\t".join(Analysis_Requested)))

#==============================================================================
# KSAT
#==============================================================================

ThisAnalysisType = 'UMS_Ksat'
DataRoot = os.path.join(MainRoot,'KSAT')


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


ID_FullMeta = map(lambda x : '_'.join(x),ID_MetaData)


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


MetaCode = numpy.zeros((nFiles,nMetaLevel),dtype='S256')

for iF in range(nFiles):
    Meta = FileNames[iF].split('_')
    if len(Meta) != nMetaLevel:
        raise ValueError('file %s has a wrong name'%FileNames[iF])
    MetaCode[iF,:] = Meta

    
MetaUnique=[]
for iM in range(nMetaLevel):
    MetaUnique.append(numpy.unique(MetaCode[:,iM]))

nMetaUnique = map(numpy.size,MetaUnique)








##==============================================================================
## HYPROP
##==============================================================================
#
#path=os.path.join(AlbonRoot,'Albon hyprop data\\exported files')
#
#filelist_Hyprop = []
#for root, dirnames, filenames in os.walk(path):
#  for filename in fnmatch.filter(filenames, '*.csv'):
#      filelist_Hyprop.append(os.path.join(root, filename))
#
#nF_H=len(filelist_Hyprop)
#print("HYPROP: %s files found" %nF_H)
#
#filenames_Hyprop=numpy.array(map(lambda x : os.path.basename(x).replace('.csv','').upper(),filelist_Hyprop)) # format : albon_ + fieldcode + center/kop + top/zool/sub + repetition
#
#Hyprop_FC=[]
#Hyprop_SFC=[]
#Hyprop_HC=[]
#Hyprop_TCS=[]
#Hyprop_K0KS=[]
#Hyprop_Rep=[]
#for iF in range(nF_H):
#    Hyprop_FC.append(filenames_Hyprop[iF].replace('ALBON_','').split('_')[0])
#    Hyprop_SFC.append(filenames_Hyprop[iF].replace('ALBON_','').split('_')[1])
#    Hyprop_HC.append(filenames_Hyprop[iF].replace('ALBON_','').split('_')[2])
#    Hyprop_TCS.append(filenames_Hyprop[iF].replace('ALBON_','').split('_')[3])
#    Hyprop_K0KS.append(filenames_Hyprop[iF].replace('ALBON_','').split('_')[4])
#    Hyprop_Rep.append(filenames_Hyprop[iF].replace('ALBON_','').split('_')[5])
#Hyprop_FC=numpy.array(Hyprop_FC)
#Hyprop_SFC=numpy.array(Hyprop_SFC)
#Hyprop_HC=numpy.array(Hyprop_HC)
#Hyprop_TCS=numpy.array(Hyprop_TCS)
#Hyprop_K0KS=numpy.array(Hyprop_K0KS)
#Hyprop_Rep=numpy.array(Hyprop_Rep)
#
##==============================================================================
## Texture
##==============================================================================
#
#TextureFilename=os.path.join(AlbonRoot,'Albon for future research\\Albon_Texture.txt')
#
#fID=open(TextureFilename)
#Tex_Header=numpy.array(fID.readline().replace('\n','').split('\t'))
#fID.close()
#
#Tex_All=numpy.genfromtxt(TextureFilename,delimiter='\t',skip_header=1,dtype='str')
#Tex_FC=Tex_All[:,0]
#Tex_FC=numpy.array([T.upper() for T in Tex_FC])
#Tex_SFC=Tex_All[:,1]
#Tex_SFC=numpy.array([T.upper() for T in Tex_SFC])
#Tex_HC=Tex_All[:,2]
#Tex_HC=numpy.array([T.upper() for T in Tex_HC])
#Tex_TCS=Tex_All[:,3]
#Tex_TCS=numpy.array([T.upper() for T in Tex_TCS])
#
#nTex=Tex_FC.size
#
#print("Texture: %s instances found" %nTex)
#
##Tex_FullCode=[]
##for iH in range(nTex):
##    Tex_FullCode.append('%s_%s_%s_%s'%(Tex_FC[iH],Tex_SFC[iH],Tex_HC[iH],Tex_TCS[iH]))
##Tex_FullCode=numpy.array(Tex_FullCode)
##
##Tex_DataHeader=Tex_Header[4:-1]
##Tex_Data=Tex_All[:,4:-1].astype('str')
##
##nTexData=Tex_DataHeader.size
#
##==============================================================================
## Bulk Density
##==============================================================================
#
#BD_Filename=os.path.join(AlbonRoot,'Albon for future research\\Albon_BD_Ksat.txt')
#
#fID=open(BD_Filename)
#BD_Header=numpy.array(fID.readline().replace('\n','').split('\t'))
#fID.close()
#
#BD_All=numpy.genfromtxt(BD_Filename,delimiter='\t',skip_header=1,dtype='str')
#BD_FC=BD_All[:,0]
#BD_FC=numpy.array([T.upper() for T in BD_FC])
#BD_SFC=BD_All[:,1]
#BD_SFC=numpy.array([T.upper() for T in BD_SFC])
#BD_HC=BD_All[:,2]
#BD_HC=numpy.array([T.upper() for T in BD_HC])
#BD_TCS=BD_All[:,3]
#BD_TCS=numpy.array([T.upper() for T in BD_TCS])
#
#BD_Count=numpy.sum(~numpy.isnan(BD_All[:,4:].astype('float')),axis=1)
#
#print("Bulk Density: %s instances found" %numpy.sum(BD_Count))
#
##==============================================================================
## Wishlist
##==============================================================================
#
#WishlistFilename=os.path.join(AlbonRoot,'Albon for future research\\Albon_Wishlist.txt')
#
#Wish_All=numpy.genfromtxt(WishlistFilename,delimiter='\t',skip_header=0,dtype='str')
#Wish_FC=Wish_All[:,0]
#Wish_SFC=Wish_All[:,1]
#Wish_HC=numpy.array(['H','C'])
#Wish_TCS=numpy.array(['TOP','CSUB','SUB'])
#
#nWish=Wish_FC.size
#nHC=Wish_HC.size
#nTCS=Wish_TCS.size
#
#print("Wishlist: %s instances found" %nWish)
#
##==============================================================================
## MERGE & WRITE
##==============================================================================
#
#import sys
#basename=os.path.basename(sys.argv[0])[:-3]
#resultname=os.path.join(AlbonRoot,'Albon for future research',basename+'.txt')
#Header='FieldCode\tSubFieldCode\tHeadCenter\tFieldCode\tTexture\tHypropK0\tHypropKS\tKsat\tPenetrologger\tBulkDensity\n'
#fID=open(resultname,'w')
#fID.write(Header)
#
#Result=[]
#for iF_SF in range(nWish):
#    FC=Wish_FC[iF_SF].upper()
#    SFC=Wish_SFC[iF_SF].upper()
#    for iHC in range(nHC):
#        HC=Wish_HC[iHC].upper()
#        for iTCS in range(nTCS):
#            TCS=Wish_TCS[iTCS].upper()
#            
#            nTex=numpy.sum((Tex_FC==FC)*(Tex_SFC==SFC)*(Tex_HC==HC)*(Tex_TCS==TCS))
#            nHyprop_K0=numpy.sum((Hyprop_FC==FC)*(Hyprop_SFC==SFC)*(Hyprop_HC==HC)*(Hyprop_TCS==TCS)*(Hyprop_K0KS=='K0'))
#            nHyprop_KS=numpy.sum((Hyprop_FC==FC)*(Hyprop_SFC==SFC)*(Hyprop_HC==HC)*(Hyprop_TCS==TCS)*(Hyprop_K0KS=='KS'))
#            nKsat=numpy.sum((Ksat_FC==FC)*(Ksat_SFC==SFC)*(Ksat_HC==HC)*(Ksat_TCS==TCS))
#            nPen=numpy.sum((Pen_FC==FC)*(Pen_SFC==SFC))
#            nBD=BD_Count[(BD_FC==FC)*(BD_SFC==SFC)*(BD_HC==HC)*(BD_TCS==TCS)][0]
#            
#            Result.append([nTex,nHyprop_K0,nHyprop_KS,nKsat,nPen,nBD])
#            
#            fID.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(FC,SFC,HC,TCS,nTex,nHyprop_K0,nHyprop_KS,nKsat,nPen,nBD))
#
#fID.close()
#
#Result=numpy.array(Result)
#
##==============================================================================
## Checksum
##==============================================================================
#
#ResultSum=Result.sum(axis=0)
#print('\nChecksum')
#print('Penetrologger: %s'%(ResultSum[4]/6.0))
#print('KSAT: %s'%ResultSum[3])
#print('Hyprop_K0: %s'%ResultSum[1])
#print('Hyprop_KS: %s'%ResultSum[2])
#print('Texture: %s'%ResultSum[0])
#print('Bulk Density: %s'%ResultSum[5])
#
#TexCheck=numpy.zeros(nTex)
#PenCheck=numpy.zeros(nF_P)
#
#for iF_SF in range(nWish):
#    FC=Wish_FC[iF_SF].upper()
#    SFC=Wish_SFC[iF_SF].upper()
#    for iHC in range(nHC):
#        HC=Wish_HC[iHC].upper()
#        for iTCS in range(nTCS):
#            TCS=Wish_TCS[iTCS].upper()
#            
#            TexCheck=(Tex_FC==FC)*(Tex_SFC==SFC)*(Tex_HC==HC)*(Tex_TCS==TCS)+TexCheck
#            PenCheck=(Pen_FC==FC)*(Pen_SFC==SFC)+PenCheck
#
#filenames_Pen[PenCheck==0]
#Tex_All[TexCheck==0,:4]
#
#pylab.show()

