# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 17:33:44 2014

@author: jan
"""

import numpy
import datetime

def GenericFileOpening(filename):
    fID=open(filename,'rb')
    line=fID.read()
    fID.close()
    
    # decode special signs, convert to normal string
    line=line.decode()
    
    # Required file format:
    # decimal sign: .
    # delimiter: ;
    # line end: '\n'

    if ';' not in line:
        line=line.replace(',',';')
    else:
        line=line.replace(',','.')
        
    line=line.replace('\r\n','\n')

    return line

def ReadHypropCSV(filename):
#    fID=open(filename,'rb')
#    line=fID.read()
#    fID.close()
#    line=line.replace(',','.')
    line = GenericFileOpening(filename)
    
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
#     indSV2=line[indSV0:].find('\n')
    Sample=float(line[indSN0:][indSN1:indSN2])
    
#    indIWC0=line.find('Initial water content [Vol%]:;')
#    indIWC1=line[indIWC0:].find(';')+1
#    indIWC2=line[indIWC0:].find('\r\n')
#    IWC=float(line[indIWC0:][indIWC1:indIWC2])
    
    indSV0=line.find('Soil volume [cm3]:;')
    indSV1=line[indSV0:].find(';')+1
    indSV2=line[indSV0:].find('\r\n')
#     indSV2=line[indSV0:].find('\n')
    SV=float(line[indSV0:][indSV1:indSV2])
    
    indDSW0=line.find('Dry soil weight [g]:;')
    indDSW1=line[indDSW0:].find(';')+1
    indDSW2=line[indDSW0:].find('\r\n')
#     indDSW2=line[indDSW0:].find('\n')
    DSW=float(line[indDSW0:][indDSW1:indDSW2])
    
    DataStartString='Date / Time;Tension bottom [hPa];Tension top [hPa];Temperature [°C];Gross weight [g];Net weight [g];Weight change [g]\r\n'
    DataEndString='Date / Time;Tension bottom [hPa];Tension top [hPa];Net weight [g]'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data[:-1])))
    Data_Time=Data[:,0]
    Data_Time=numpy.array(list(map(lambda x : datetime.datetime.strptime(x,'%d*%m*%Y %H:%M:%S'),Data_Time)))
    Data_Time_Sec=numpy.array(list(map(lambda x: x.total_seconds(),Data_Time-Data_Time[0])))
    
    Data_TensionB=Data[:,1].astype('float')
    Data_TensionT=Data[:,2].astype('float')
    Data_Temp=Data[:,3].astype('float')
    Data_MassNet=Data[:,5].astype('float')

    return Sample,SV,DSW,Data_Time,Data_Time_Sec,Data_TensionB,Data_TensionT,Data_Temp,Data_MassNet

def ReadHyprop_VERSION(filename):
#    fID=open(filename,'rb')
#    line=fID.read()
#    fID.close()
    line = GenericFileOpening(filename)

    NewVersion=0
    if "Model Description,PDI-variant of the bimodal constrained van Genuchten model " in line:
        NewVersion=1
    return NewVersion

def ReadHypropCSV_V2_RawData(filename):
    line = GenericFileOpening(filename)

#    line=line.replace('Soil surface area [cm2]:','Soil volume [cm3]:')
    
#    indIWC0=line.find('Initial water content [Vol%]:;')
#    indIWC1=line[indIWC0:].find(';')+1
#    indIWC2=line[indIWC0:].find('\r\n')
#    IWC=float(line[indIWC0:][indIWC1:indIWC2])
    
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
    Sample=line[indSN0:][indSN1:indSN2]
    
    ## Old version
#    indSV0=line.find('Soil volume [cm3]:;')
#    indSV1=line[indSV0:].find(';')+1
#    indSV2=line[indSV0:].find('\n')
#    SV=float(line[indSV0:][indSV1:indSV2])
    
    ## New version
    indSA0=line.find('Soil surface area [cm2]:;')
    indSA1=line[indSA0:].find(';')+1
    indSA2=line[indSA0:].find('\n')
    SA=float(line[indSA0:][indSA1:indSA2])
    ## BUG
    SA=50.0
    
    indSH0=line.find('Soil column height [cm]:;')
    indSH1=line[indSH0:].find(';')+1
    indSH2=line[indSH0:].find('\n')
    SH=float(line[indSH0:][indSH1:indSH2])
    
    SV=SA*SH
    
    
    indDSW0=line.find('Dry soil weight [g]:;')
    indDSW1=line[indDSW0:].find(';')+1
    indDSW2=line[indDSW0:].find('\n')
    DSW=float(line[indDSW0:][indDSW1:indDSW2])
    
    DataStartString='Date / Time;Tension bottom [hPa];Tension top [hPa];Temperature ['
    DataEndString='Date / Time;Gross weight [g];Net weight [g];Weight change [g]\n'
    indData1=line.find(DataStartString)
    indData2=line.find(DataEndString)
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()[1:] # skip the celsius thing
    Data=numpy.array(list(map(lambda x: x.split(';'),Data[:-1])))
    Data_Time=Data[:,0]
    Data1_Time=numpy.array(list(map(lambda x : datetime.datetime.strptime(x,'%d*%m*%Y %H:%M:%S'),Data_Time)))
    Data1_Time_Sec=numpy.array(list(map(lambda x: x.total_seconds(),Data1_Time-Data1_Time[0])))
    
    Data1_TensionB=Data[:,1].astype('float')
    Data1_TensionT=Data[:,2].astype('float')
    Data1_Temp=Data[:,3].astype('float')
    
    
    DataStartString='Date / Time;Gross weight [g];Net weight [g];Weight change [g]\n'
    DataEndString='Date / Time;Tension bottom [hPa];Tension top [hPa];Net weight [g]\n'
    indData1=line.find(DataStartString)
    indData2=line.find(DataEndString)
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data[:-1])))
    Data_Time=Data[:,0]
    Data2_Time=numpy.array(list(map(lambda x : datetime.datetime.strptime(x,'%d*%m*%Y %H:%M:%S'),Data_Time)))
    Data2_Time_Sec=numpy.array(list(map(lambda x: x.total_seconds(),Data2_Time-Data2_Time[0])))
    
#    Data2_MassBruto=Data[:,1].astype('float')
    Data2_MassNet=Data[:,2].astype('float')
#    Data2_dMass=Data[:,3].astype('float')

    return Sample,SV,DSW,Data1_Time,Data1_Time_Sec,Data1_TensionB,Data1_TensionT,Data1_Temp,Data2_Time,Data2_Time_Sec,Data2_MassNet


def ReadHypropCSV_V2_TensMass(filename):
    line = GenericFileOpening(filename)

    line=line.replace('Soil surface area [cm2]:','Soil volume [cm3]:')
    
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
    Sample=line[indSN0:][indSN1:indSN2]

    DataStartString='Date / Time;Tension bottom [hPa];Tension top [hPa];Net weight [g]\n'
    DataEndString='pF [-];Water Content [Vol%]\n'
    indData1=line.find(DataStartString)
    indData2=line.find(DataEndString)
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data[:-1])))
    Data_Time=Data[:,0]
    Data3_Time=numpy.array(list(map(lambda x : datetime.datetime.strptime(x,'%d*%m*%Y %H:%M:%S'),Data_Time)))
    Data3_Time_Sec=numpy.array(list(map(lambda x: x.total_seconds(),Data3_Time-Data3_Time[0])))
    
    Data3_TensionB=Data[:,1].astype('float')
    Data3_TensionT=Data[:,2].astype('float')
    Data3_MassNet=Data[:,3].astype('float')

    return Sample,Data3_Time,Data3_Time_Sec,Data3_TensionB,Data3_TensionT,Data3_MassNet


def ReadHypropCSV_V2_NiceData(filename):
    line = GenericFileOpening(filename)
    
    line=line.replace('Soil surface area [cm2]:','Soil volume [cm3]:')
        
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
    Sample=line[indSN0:][indSN1:indSN2]
    
    NewVersion='pF [-];Water Content [Vol%];Weight [-]\n' in line
    if NewVersion:
        DataStartString='pF [-];Water Content [Vol%];Weight [-]\n'
        DataEndString='pF [-];log 10 K [cm/d];Weight [-]\n'
    else:
        DataStartString='pF [-];Water Content [Vol%]\n'
        DataEndString='pF [-];log 10 K [cm/d]\n'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)+indData1
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))

    Data4_pF=Data[:,0].astype('float')
    Data4_WC=Data[:,1].astype('float')
    
#    NewVersion=False
    
    if NewVersion:
        DataStartString='pF [-];log 10 K [cm/d];Weight [-]\n'
        DataEndString='Water Content [Vol%];log 10 K [cm/d]\n'
    else:
        DataStartString='pF [-];log 10 K [cm/d]\n'
        DataEndString='Water Content [Vol%];log 10 K [cm/d]\n'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)+indData1
    
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
    Data5_pF=Data[:,0].astype('float')
    Data5_Kh=Data[:,1].astype('float')
    
    if NewVersion:
        DataStartString='Water Content [Vol%];log 10 K [cm/d]\n'
        DataEndString='pF [-];Water Content [Vol%]\n'
    else:
        DataStartString='Water Content [Vol%];log 10 K [cm/d]\n'
        DataEndString='pF [-];Water Content [Vol%]\n'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)+indData1
    
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
    Data6_WC=Data[:,0].astype('float')
    Data6_Kh=Data[:,1].astype('float')
    
    return Sample,Data4_pF,Data4_WC,Data5_pF,Data5_Kh,Data6_WC,Data6_Kh


def ReadHypropCSV_V2_ModelData(filename):
    line = GenericFileOpening(filename)
    
    line=line.replace('Soil surface area [cm2]:','Soil volume [cm3]:')
        
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
    Sample=line[indSN0:][indSN1:indSN2]
    
    NewVersion='pF [-];Water Content [Vol%];Weight [-]\n' in line
    if NewVersion:
        DataStartString='pF [-];Water Content [Vol%]\n'
        DataEndString='pF [-];log 10 K [cm/d]\n'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)+indData1
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))

    Data4_pF=Data[:,0].astype('float')
    Data4_WC=Data[:,1].astype('float')
    
    if NewVersion:
        DataStartString='pF [-];log 10 K [cm/d]\n'
        DataEndString='Water Content [Vol%];log 10 K [cm/d]\n'
    else:
        DataStartString='pF [-];log 10 K [cm/d]\n'
        DataEndString='Water Content [Vol%];log 10 K [cm/d]\n'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)+indData1
    
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
    Data5_pF=Data[:,0].astype('float')
    Data5_Kh=Data[:,1].astype('float')
    
    if NewVersion:
        DataStartString='Water Content [Vol%];log 10 K [cm/d]\n'
        DataEndString='pF [-];Water Content [Vol%]\n'
    else:
        DataStartString='Water Content [Vol%];log 10 K [cm/d]\n'
        DataEndString='pF [-];Water Content [Vol%]\n'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)+indData1
    
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
    Data6_WC=Data[:,0].astype('float')
    Data6_Kh=Data[:,1].astype('float')
    
    return Sample,Data4_pF,Data4_WC,Data5_pF,Data5_Kh,Data6_WC,Data6_Kh


def ReadHypropCSV_V2_Fit(filename):
    from io import StringIO
    import re
    
    line = GenericFileOpening(filename)
    
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
#     indSV2=line[indSV0:].find('\n')
    Sample=line[indSN0:][indSN1:indSN2]
    
    DataStartString='Is fitted;Parameter;Value;Min;Max;2.5%;97.5%;Unit;Model\n'
#    DataEndString=' ;alpha;n;th_r;th_s;tau;Ks\n'
    DataEndString=' ;alpha;n;th_r;th_s;alpha2;n2;w2;tau;Ks;omega\n'
    indData1=line.find(DataStartString)
    indData2=line.find(DataEndString)
    NewVersion=0
    if indData1 == -1:
        # v3.0 (2017)
        NewVersion=1
        DataStartString='Is fitted;Parameter;Value;Min;Max;2.5%;97.5%;Unit\n'
        DataEndString=' ;alpha1;n1;th_r;th_s;pF_dry;alpha2;n2;w2;Ks;tau;omega;a\n'
        indData1=line.find(DataStartString)
        indData2=line.find(DataEndString)
        
    if indData2 == -1:
        # v4.0 (2019)
        NewVersion=1
        DataStartString='Is fitted;Parameter;Value;Min;Max;2.5%;97.5%;Unit\n'
        DataEndString=' ;alpha1;n1;th_r;th_s;alpha2;n2;w2;Ks;tau;omega\n'
        
        indData1=line.find(DataStartString)
        indData2=line.find(DataEndString)
        
    DataString=line[indData1+len(DataStartString):indData2]
    DataString=re.sub(r'\([^)]*\)', '', DataString)# remove text between brackets to avoid confusion
    DataString=DataString.replace('*','') 
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
    if NewVersion:
        Data=Data[numpy.array([0,1,2,3,5,6,7,9,8,10]),:] # ;alpha1;n1;th_r;th_s;pF_dry;alpha2;n2;w2;Ks;tau;omega;a
        
    if indData1 !=-1:
        Data7_isfit=Data[:,0].astype('float')
        Data7_Par=Data[:,1]
        Data7_ParVal=Data[:,2].astype('float')
        Data7_Model=Data[:,-1]
    else:
        Data7_isfit=numpy.nan
        Data7_Par=numpy.nan
        Data7_ParVal=numpy.nan
        Data7_Model=numpy.nan
        
    return Sample,Data7_isfit,Data7_Par,Data7_ParVal,Data7_Model,NewVersion


def ReadHypropCSV_V3_RawData(filename):
    line = GenericFileOpening(filename)
    
#    line=line.replace('Soil surface area [cm2]:','Soil volume [cm3]:')
    
#    indIWC0=line.find('Initial water content [Vol%]:;')
#    indIWC1=line[indIWC0:].find(';')+1
#    indIWC2=line[indIWC0:].find('\r\n')
#    IWC=float(line[indIWC0:][indIWC1:indIWC2])
    
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
    Sample=line[indSN0:][indSN1:indSN2]
    
    ## Old version
#    indSV0=line.find('Soil volume [cm3]:;')
#    indSV1=line[indSV0:].find(';')+1
#    indSV2=line[indSV0:].find('\n')
#    SV=float(line[indSV0:][indSV1:indSV2])
    
    ## New version
    indSA0=line.find('Soil surface area [cm2]:;')
    indSA1=line[indSA0:].find(';')+1
    indSA2=line[indSA0:].find('\n')
    SA=float(line[indSA0:][indSA1:indSA2])
    ## BUG
    SA=50.0
    
    indSH0=line.find('Soil column height [cm]:;')
    indSH1=line[indSH0:].find(';')+1
    indSH2=line[indSH0:].find('\n')
    SH=float(line[indSH0:][indSH1:indSH2])

    SV=SA*SH
    
    
    indDSW0=line.find('Dry soil weight [g]:;')
    indDSW1=line[indDSW0:].find(';')+1
    indDSW2=line[indDSW0:].find('\n')
    DSW=float(line[indDSW0:][indDSW1:indDSW2])
    
    DataStartString='Date / Time;Tension bottom [hPa];Tension top [hPa];Temperature [\xc2\xb0C]'
    DataEndString='Date / Time;Gross weight [g];Net weight [g];Weight change [g]'
    indData1=line.find(DataStartString)
    indData2=line.find(DataEndString)
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data[1:])))
    Data_Time=Data[:,0]
    # print(Data_Time)
    Data1_Time=numpy.array(list(map(lambda x : datetime.datetime.strptime(x,'%d-%b-%y %I:%M:%S %p'),Data_Time)))
    Data1_Time_Sec=numpy.array(list(map(lambda x: x.total_seconds(),Data1_Time-Data1_Time[0])))
    
    Data1_TensionB=Data[:,1].astype('float')
    Data1_TensionT=Data[:,2].astype('float')
    Data1_Temp=Data[:,3].astype('float')
    
    
    DataStartString='Date / Time;Gross weight [g];Net weight [g];Weight change [g]'
    DataEndString='Date / Time;Tension bottom [hPa];Tension top [hPa];Net weight [g]'
    indData1=line.find(DataStartString)
    indData2=line.find(DataEndString)
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data[1:])))
    Data_Time=Data[:,0]
    Data2_Time=numpy.array(list(map(lambda x : datetime.datetime.strptime(x,'%d-%b-%y %I:%M:%S %p'),Data_Time)))
    Data2_Time_Sec=numpy.array(list(map(lambda x: x.total_seconds(),Data2_Time-Data2_Time[0])))
    
#    Data2_MassBruto=Data[:,1].astype('float')
    Data2_MassNet=Data[:,2].astype('float')
#    Data2_dMass=Data[:,3].astype('float')

    return Sample,SV,DSW,Data1_Time,Data1_Time_Sec,Data1_TensionB,Data1_TensionT,Data1_Temp,Data2_Time,Data2_Time_Sec,Data2_MassNet

def ReadHypropCSV_V3_NiceData(filename):
    line = GenericFileOpening(filename)
    
    line=line.replace('Soil surface area [cm2]:','Soil volume [cm3]:')
        
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
    Sample=line[indSN0:][indSN1:indSN2]
    
    NewVersion='pF [-];Water Content [Vol%];Weight [-]' in line
    if NewVersion:
        DataStartString='pF [-];Water Content [Vol%];Weight [-]'
        DataEndString='pF [-];log 10 K [cm/d];Weight [-]'
    else:
        DataStartString='pF [-];Water Content [Vol%]'
        DataEndString='pF [-];log 10 K [cm/d]'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)+indData1
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()[1:]
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
    Data4_pF=Data[:,0].astype('float')
    Data4_WC=Data[:,1].astype('float')
    
#    NewVersion=False
    
    if NewVersion:
        DataStartString='pF [-];log 10 K [cm/d];Weight [-]'
        DataEndString='Water Content [Vol%];log 10 K [cm/d]'
    else:
        DataStartString='pF [-];log 10 K [cm/d]'
        DataEndString='Water Content [Vol%];log 10 K [cm/d]'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)+indData1
    
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()[1:]
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
    Data5_pF=Data[:,0].astype('float')
    Data5_Kh=Data[:,1].astype('float')
    
    if NewVersion:
        DataStartString='Water Content [Vol%];log 10 K [cm/d]'
        DataEndString='pF [-];Water Content [Vol%]'
    else:
        DataStartString='Water Content [Vol%];log 10 K [cm/d]'
        DataEndString='pF [-];Water Content [Vol%]'
    indData1=line.find(DataStartString)
    indData2=line[indData1:].find(DataEndString)+indData1
    DataString=line[indData1+len(DataStartString):indData2]
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()[1:]
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
    Data6_WC=Data[:,0].astype('float')
    Data6_Kh=Data[:,1].astype('float')
    
    return Sample,Data4_pF,Data4_WC,Data5_pF,Data5_Kh,Data6_WC,Data6_Kh


def ReadHypropCSV_V3_Fit(filename):
    from io import StringIO
    import re

    line = GenericFileOpening(filename)
    
    
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
#     indSV2=line[indSV0:].find('\n')
    Sample=line[indSN0:][indSN1:indSN2]
    
    DataStartString='Is fitted;Parameter;Value;Min;Max;2.5%;97.5%;Unit;Model'
#    DataEndString=' ;alpha;n;th_r;th_s;tau;Ks\n'
    DataEndString=' ;alpha;n;th_r;th_s;alpha2;n2;w2;tau;Ks;omega'
    indData1=line.find(DataStartString)
    indData2=line.find(DataEndString)
    NewVersion=0
    if indData1 == -1:
        NewVersion=1
        DataStartString='Is fitted;Parameter;Value;Min;Max;2.5%;97.5%;Unit'
        DataEndString=' ;alpha1;n1;th_r;th_s;pF_dry;alpha2;n2;w2;Ks;tau;omega;a'
        indData1=line.find(DataStartString)
        indData2=line.find(DataEndString)
        
    DataString=line[indData1+len(DataStartString):indData2]
    DataString=re.sub(r'\([^)]*\)', '', DataString)# remove text between brackets to avoid confusion
    DataString=DataString.replace('*','') 
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()[1:]
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
    if NewVersion:
        Data=Data[numpy.array([0,1,2,3,5,6,7,9,8,10]),:]
        
    if indData1 !=-1:
        Data7_isfit=Data[:,0].astype('float')
        Data7_Par=Data[:,1]
        Data7_ParVal=Data[:,2].astype('float')
        Data7_Model=Data[:,-1]
    else:
        Data7_isfit=numpy.nan
        Data7_Par=numpy.nan
        Data7_ParVal=numpy.nan
        Data7_Model=numpy.nan
        
    return Sample,Data7_isfit,Data7_Par,Data7_ParVal,Data7_Model,NewVersion


#==============================================================================
# AIR ENTRY AND ALL
#==============================================================================

def RunMean_noedgeeffect(X,n):
    ## Make it smoooooooth without edge effects
    import numpy
    nX=X.size
    X_extend=numpy.zeros(nX+2*n)
    X_extend[:n]=X[0]
    X_extend[n:nX+n]=X
    X_extend[nX+n:]=X[-1]
    Xsmooth=numpy.convolve(X_extend,numpy.ones(n)/n,'same')[n:n+nX]
    return Xsmooth
    
def AirEntryExtrapolation(Data1_Time_Sec,Data1_TensionB,Data1_TensionT,AElimit=-3e-7,pAirEntry=8800.0,nForInterpol=20):
    import numpy
    
    dT=Data1_Time_Sec[1:]-Data1_Time_Sec[:-1]
    Data_Time_Sec_mid=(Data1_Time_Sec[1:]+Data1_Time_Sec[:-1])/2
    ddT=(dT[1:]+dT[:-1])/2
    
    ## Smoooooth
    Data1_TensionB_smooth=RunMean_noedgeeffect(Data1_TensionB,10)
    Data1_TensionT_smooth=RunMean_noedgeeffect(Data1_TensionT,10)
    dData_TensionB=(Data1_TensionB_smooth[1:] - Data1_TensionB_smooth[:-1])/dT
    dData_TensionT=(Data1_TensionT_smooth[1:] - Data1_TensionT_smooth[:-1])/dT
    
#    dData_TensionB=(Data1_TensionB[1:] - Data1_TensionB[:-1])/dT
#    dData_TensionT=(Data1_TensionT[1:] - Data1_TensionT[:-1])/dT
    
    ddData_TensionB=(dData_TensionB[1:] - dData_TensionB[:-1])/ddT
    ddData_TensionT=(dData_TensionT[1:] - dData_TensionT[:-1])/ddT
    
    marge=20
    IStopT=numpy.argmax(Data1_TensionT)
    IStopB=numpy.argmax(Data1_TensionB)
    IAirEntryT=numpy.argmin(ddData_TensionT[IStopT+marge:])+IStopT+marge+1
#    IAirEntryB=numpy.argmax(ddData_TensionB[IStopB+marge:])+IStopB+marge+1
    IAirEntryB=numpy.argmin(ddData_TensionB[IAirEntryT-1:])+IAirEntryT+1
    
#    AElimit=-3e-7 # hpa/s2
#    pAirEntry=8800.0 # hpa
#    nForInterpol=200
#    nForInterpol=20
    
    if min(ddData_TensionT[IStopT+marge:])<AElimit:
        
        Data_TensionT_forinterpol=numpy.append(Data1_TensionT[IStopT-nForInterpol:IStopT],pAirEntry)
        TT_forinterpol=numpy.append(Data1_Time_Sec[IStopT-nForInterpol:IStopT],Data1_Time_Sec[IAirEntryT])
        pT=numpy.polyfit(TT_forinterpol,Data_TensionT_forinterpol,3)
        vT=numpy.polyval(pT,Data1_Time_Sec[IStopT:IAirEntryT])
        
        Data_TensionT_TOT=numpy.concatenate((Data1_TensionT[:IStopT],vT,numpy.array([pAirEntry,])))
        
        if (IAirEntryT>IStopB) & (min(ddData_TensionB[IStopB+marge:])<AElimit):
            Data_TensionB_forinterpol=numpy.append(Data1_TensionB[IStopB-nForInterpol:IStopB],pAirEntry)
            TB_forinterpol=numpy.append(Data1_Time_Sec[IStopB-nForInterpol:IStopB],Data1_Time_Sec[IAirEntryB])
            pB=numpy.polyfit(TB_forinterpol,Data_TensionB_forinterpol,3)
            vB=numpy.polyval(pB,Data1_Time_Sec[IStopB:IAirEntryT+1])
        
            Data_TensionB_TOT=numpy.concatenate((Data1_TensionB[:IStopB],vB))
            Data1_Time_Sec_TOT=Data1_Time_Sec[:IAirEntryT+1]
        else:
            iStop=min([IStopB,IAirEntryT])+1
            Data_TensionT_TOT=Data_TensionT_TOT[:iStop]
            Data_TensionB_TOT=Data1_TensionB[:iStop]
            Data1_Time_Sec_TOT=Data1_Time_Sec[:iStop]
    
    else:
        Data_TensionT_TOT=Data1_TensionT[:IStopB]
        Data_TensionB_TOT=Data1_TensionB[:IStopB]
        Data1_Time_Sec_TOT=Data1_Time_Sec[:IStopB]
        
    TensionMean_TOT=(Data_TensionT_TOT+Data_TensionB_TOT)/2 # hPa
    
    return Data1_Time_Sec_TOT,TensionMean_TOT,Data_TensionB_TOT,Data_TensionT_TOT,IStopB,IStopT,IAirEntryB,IAirEntryT

def ReadHypropCSV_V4_Fit(filename):
    from io import StringIO
    import re
    
    line = GenericFileOpening(filename)
        
    indSN0=line.find('Sample name:;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
#     indSV2=line[indSV0:].find('\n')
    Sample=line[indSN0:][indSN1:indSN2]

    indSN0=line.find('Model Description;')
    indSN1=line[indSN0:].find(';')+1
    indSN2=line[indSN0:].find('\n')
#     indSV2=line[indSV0:].find('\n')
    ModelDescription=line[indSN0:][indSN1:indSN2]
    
    DataStartString='Is fitted;Parameter;Value;Min;Max;2.5%;97.5%;Unit;Model\n'
#    DataEndString=' ;alpha;n;th_r;th_s;tau;Ks\n'
    DataEndString=' ;alpha;n;th_r;th_s;alpha2;n2;w2;tau;Ks;omega\n'
    indData1=line.find(DataStartString)
    indData2=line.find(DataEndString)
    NewVersion=0
    
    if indData1 == -1:
        NewVersion=1
        DataStartString='Is fitted;Parameter;Value;Min;Max;2.5%;97.5%;Unit\n'
#        DataEndString=' ;alpha1;n1;th_r;th_s;pF_dry;alpha2;n2;w2;Ks;tau;omega;a\n'
        DataEndString=' ;'
        indData1=line.find(DataStartString)
        indData2=line[indData1:].find(DataEndString)+indData1
        
    DataString=line[indData1+len(DataStartString):indData2]
    DataString=re.sub(r'\([^)]*\)', '', DataString)# remove text between brackets to avoid confusion
    DataString=DataString.replace('*','') 
    #Data=numpy.fromstring(DataString,dtype='str',sep='\r\n')
    Data=DataString.splitlines()
    Data=numpy.array(list(map(lambda x: x.split(';'),Data)))
    
#    if NewVersion:
#        Data=Data[numpy.array([0,1,2,3,5,6,7,9,8,10]),:]

    if indData1 !=-1:
        Data7_isfit=Data[:,0].astype('float')
        Data7_Par=Data[:,1]
        Data7_ParVal=Data[:,2].astype('float')
        Data7_Model=Data[:,-1]
    else:
        Data7_isfit=numpy.nan
        Data7_Par=numpy.nan
        Data7_ParVal=numpy.nan
        Data7_Model=numpy.nan
        
    return Sample,Data7_isfit,Data7_Par,Data7_ParVal,Data7_Model,NewVersion,ModelDescription


#==============================================================================
# Water Retention
#==============================================================================

def Vangenuchten4(h,Par):
    thetar,thetas,alpha,n=Par
    m=1.0-1.0/n
    theta=thetar+(thetas-thetar)/(1.0+(alpha*numpy.abs(h))**n)**m
    return theta

def Vangenuchten5(h,Par):
    thetar,thetas,alpha,n,m=Par
    theta=thetar+(thetas-thetar)/(1.0+(alpha*numpy.abs(h))**n)**m
    return theta

def Durner(h,Par):
    # Double van Genuchten
    thetar,thetas,w1,alpha1,n1,alpha2,n2=Par
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    w2=1.0-w1
#    theta=thetar+(thetas-thetar)*(w1/(1.0+(alpha1*numpy.abs(h))**n1)**m1+w2/(1.0+(alpha2*numpy.abs(h))**n2)**m2)
    theta=thetar+(thetas-thetar)*(w1*Vangenuchten4(h,[0.0,1.0,alpha1,n1])+w2*Vangenuchten4(h,[0.0,1.0,alpha2,n2]))
    return theta

def PDI_AdsorptiveSaturation(h,Par):
    import numpy
    if len(Par)==4:
        thetar,thetas,alpha,n=Par
        pF0=6.8
    elif len(Par)==5:
        thetar,thetas,alpha,n,pF0=Par
#    pFa=numpy.log10(1/alpha)
    pFa=numpy.log10(1/alpha)
    b=0.1+(0.2/n**2)*(1.0-numpy.exp(-(thetar/(thetas-thetar))**2))
    
    Sad=1+(1.0/(pFa-pF0))*(numpy.log10(numpy.abs(h))-pFa+b*numpy.log(1.0+numpy.exp((pFa-numpy.log10(numpy.abs(h)))/b)))
    Sad[Sad<0]=0
    return Sad

def PDI_Ret(h,Par):
    if len(Par)==4:
        thetar,thetas,alpha,n=Par
        pF0=6.8
    elif len(Par)==5:
        thetar,thetas,alpha,n,pF0=Par
    h0=10**pF0
    SadPar=[thetar,thetas,alpha,n,pF0]
    Sad=PDI_AdsorptiveSaturation(h,SadPar)
    
    Scap=(Vangenuchten4(h,[0.0,1.0,alpha,n])-Vangenuchten4(h0,[0.0,1.0,alpha,n])) \
        /(1-Vangenuchten4(h0,[0.0,1.0,alpha,n])) # Scaled to be 0 at pF0
    
    theta=thetar*Sad+(thetas-thetar)*Scap
    theta[theta<0]=0
    return theta

def BimodalPDI_Ret(h,Par):
    if len(Par)==7:
        thetar,thetas,w1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
    elif len(Par)==8:
        thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0=Par
    w2=1.0-w1
    h0=10**pF0
    SadPar=[thetar,thetas,[alpha1,alpha2][numpy.argmax([alpha1,alpha2])],[n1,n2][numpy.argmax([alpha1,alpha2])],pF0]
    Sad=PDI_AdsorptiveSaturation(h,SadPar)
    
    Scap=((w1*Vangenuchten4(h,[0.0,1.0,alpha1,n1])+w2*Vangenuchten4(h,[0.0,1.0,alpha2,n2]))-(w1*Vangenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*Vangenuchten4(h0,[0.0,1.0,alpha2,n2]))) \
        /(1-(w1*Vangenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*Vangenuchten4(h0,[0.0,1.0,alpha2,n2]))) # Scaled to be 0 at pF0
    
    theta=thetar*Sad+(thetas-thetar)*Scap
    theta[theta<0]=0
    return theta
    
#==============================================================================
# Hydraulic Conductivity
#==============================================================================

def Mualem4(h,Par):
    Ks,l,alpha,n=Par
    m=1.0-1.0/n
    K=Ks*Vangenuchten5(h,[0.0,1.0,alpha,n,m])**l * (alpha*(1.0-(1.0-Vangenuchten5(h,[0.0,1.0,alpha,n,m])**(1.0/m))**m))**2/(alpha**2)
    return K
    
def Mualem5(h,Par):
    Ks,l,alpha,n,m=Par
    K=Ks*Vangenuchten5(h,[0.0,1.0,alpha,n,m])**l * (alpha*(1.0-(1.0-Vangenuchten5(h,[0.0,1.0,alpha,n,m])**(1.0/m))**m))**2/(alpha**2)
    return K

def MualemBimodal(h,Par):
    Ks,l,w1,alpha1,n1,alpha2,n2=Par
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    w2=1.0-w1
    K=Ks*(Durner(h,[0.0,1.0,w1,alpha1,n1,alpha2,n2])**l \
        *(w1*alpha1 * (1.0 - (alpha1*numpy.abs(h))**(n1-1) * Vangenuchten5(h,[0.0,1.0,alpha1,n1,m1])) \
        + w2*alpha2 * (1.0 - (alpha2*numpy.abs(h))**(n2-1) * Vangenuchten5(h,[0.0,1.0,alpha2,n2,m2])))**2 \
        /(w1*alpha1+w2*alpha2)**2)
    return K

def PetersDurnerIBimodal(h,Par):
    # old version hyprop
    Ks,l,w1,omega1,alpha1,n1,alpha2,n2=Par
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    w2=1.0-w1
    omega2=1.0-omega1
    K=Ks*(omega2*(Durner(h,[0.0,1.0,w1,alpha1,n1,alpha2,n2])**l \
        *(w1*alpha1 * (1.0 - (alpha1*numpy.abs(h))**(n1-1) * Vangenuchten5(h,[0.0,1.0,alpha1,n1,m1])) \
        + w2*alpha2 * (1.0 - (alpha2*numpy.abs(h))**(n2-1) * Vangenuchten5(h,[0.0,1.0,alpha2,n2,m2])))**2 \
        /(w1*alpha1+w2*alpha2)**2) \
        +omega1*(Durner(h,[0.0,1.0,w1,alpha1,n1,alpha2,n2])**l))
    return K

    
def PDI_Kfilm(h,Par):
    if len(Par)==7:
        thetar,thetas,Ks,l,omega1,alpha,n=Par
        pF0=6.8
        a=-1.5
    elif len(Par)==9:
        thetar,thetas,Ks,l,omega1,alpha,n,pF0,a=Par
    
    ha=1.0/numpy.max(alpha)
    h0=10**pF0
    
    SadPar=[0.0,1.0,alpha,n,pF0]
    Sad=PDI_AdsorptiveSaturation(h,SadPar)
    Kfilm=(h0/ha)**(a*(1.0-Sad))
    
    return Kfilm
    
    
def PDI_K(h,Par,VaporConductivity=False):
    if len(Par)==7:
        thetar,thetas,Ks,l,omega1,alpha,n=Par
        pF0=6.8
        a=-1.5
    elif len(Par)==9:
        thetar,thetas,Ks,l,omega1,alpha,n,pF0,a=Par
    
    m=1.0-1.0/n
    omega2=1.0-omega1
    ha=1.0/alpha
    pFa=numpy.log10(ha)
    h0=10**pF0
    Par_Ret=[thetar,thetas,alpha,n,pF0]
    
    Kfilm=PDI_Kfilm(h,Par)
    
    Scap=(Vangenuchten4(h,[0.0,1.0,alpha,n])-Vangenuchten4(h0,[0.0,1.0,alpha,n])) \
        /(1-Vangenuchten4(h0,[0.0,1.0,alpha,n])) # Scaled to be 0 at pF0
        
    Kcap=Scap**l \
        *(1.0 -
        (1.0 - Vangenuchten5(h,[0.0,1.0,alpha,n,m])**(1/m))**m \
        /(1.0 - Vangenuchten5(h0,[0.0,1.0,alpha,n,m])**(1/m))**m \
        )**2
    
    if VaporConductivity:
        Kvap=PDI_Kvap(h,Par_Ret,SWR_Model=PDI_Ret)
        K=Ks*(omega2*Kcap+omega1*Kfilm)+Kvap
    else:
        K=Ks*(omega2*Kcap+omega1*Kfilm)
    
    return K


def BimodalPDI_Kfilm(h,Par):
    if len(Par)==10:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
        a=-1.5
    elif len(Par)==12:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2,pF0,a=Par
    ha=1.0/numpy.max([alpha1,alpha2])
    h0=10**pF0
    
    SadPar=[0.0,1.0,[alpha1,alpha2][numpy.argmax([alpha1,alpha2])],[n1,n2][numpy.argmax([alpha1,alpha2])],pF0]
    Sad=PDI_AdsorptiveSaturation(h,SadPar)
    Kfilm=(h0/ha)**(a*(1.0-Sad))
    
    return Kfilm

def BimodalPDI_K(h,Par,VaporConductivity=False):
    if len(Par)==10:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
        a=-1.5
    elif len(Par)==12:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2,pF0,a=Par
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    w2=1.0-w1
    omega2=1.0-omega1
    
    h0=10**pF0
    Par_Ret=[thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0]
    
    Kfilm=BimodalPDI_Kfilm(h,Par)
    Scap=((w1*Vangenuchten4(h,[0.0,1.0,alpha1,n1])+w2*Vangenuchten4(h,[0.0,1.0,alpha2,n2]))-(w1*Vangenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*Vangenuchten4(h0,[0.0,1.0,alpha2,n2]))) \
        /(1-(w1*Vangenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*Vangenuchten4(h0,[0.0,1.0,alpha2,n2]))) # Scaled to be 0 at pF0
        
    Kcap=Scap**l \
        *(1.0 -
        (w1*alpha1 * (1.0 - Vangenuchten5(h,[0.0,1.0,alpha1,n1,m1])**(1/m1))**m1 \
        + w2*alpha2 * (1.0 - Vangenuchten5(h,[0.0,1.0,alpha2,n2,m2])**(1/m2))**m2) \
        /(w1*alpha1 * (1.0 - Vangenuchten5(h0,[0.0,1.0,alpha1,n1,m1])**(1/m1))**m1 \
        + w2*alpha2 * (1.0 - Vangenuchten5(h0,[0.0,1.0,alpha2,n2,m2])**(1/m2))**m2) )**2
    
    if VaporConductivity:
        Kvap=PDI_Kvap(h,Par_Ret,SWR_Model=BimodalPDI_Ret)
        K=Ks*(omega2*Kcap+omega1*Kfilm)+Kvap
    else:
        K=Ks*(omega2*Kcap+omega1*Kfilm)
    
    return K

def PDI_Kvap(h,Par,SWR_Model=BimodalPDI_Ret):
    import numpy
    thetas=Par[1]
    rho_w=1e3 # kg/m3
    M=0.018015 # kg/mol
    g=9.81 # m/s2
    R=8.134 #J/mol/kg
    T=300 # K
    
    theta_A=thetas-SWR_Model(h,Par)
    Tort=theta_A**(7.0/3)/thetas**2
    D_A=2.14e-5*(T/273.15)**2
    D=theta_A*D_A*Tort
    
    rho_sv=1e-3*numpy.exp(31.3716-6014.76/T-7.924951e-3*T)/T
    Hr=numpy.exp(-numpy.abs(h*100.0)*M*g/(R*T)) # Convert h to meters!
    
    Kvap=rho_sv/rho_w*Hr*D*M*g/(R*T)
    
    return Kvap

#==============================================================================
# Capacity
#==============================================================================

def PDI_Cap(h,Par):
    if len(Par)==4:
        thetar,thetas,alpha,n=Par
        pF0=6.8
    elif len(Par)==5:
        thetar,thetas,alpha,n,pF0=Par
    h0=10**pF0
    pFa=numpy.log10(1/alpha)
    pF=numpy.log10(numpy.abs(h))
    m=1.0-1.0/n
    b=0.1+(0.2/n**2)*(1.0-numpy.exp(-(thetar/(thetas-thetar))**2))
    
    Cap_AD=(h*numpy.log(10.0)*(pF-pF0))**(-1.0) * (1.0 - numpy.exp((pFa-pF)/b)/(1.0+numpy.exp((pFa-pF)/b)))
    
    dSdh=-alpha*n*m*(alpha*h)**(n-1.0)*(1+(alpha*h)**n)**(-m-1.0)
    
    Cap=(thetas-thetar)/(1.0-Vangenuchten5(h0,[0.0,1.0,alpha,n,m])) * dSdh + thetar*Cap_AD
    
    return Cap

def BimodalPDI_Cap(h,Par):
    if len(Par)==7:
        thetar,thetas,w1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
    elif len(Par)==8:
        thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0=Par
    w2=1.0-w1
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    h0=10**pF0
    
    SadPar=[thetar,thetas,[alpha1,alpha2][numpy.argmax([alpha1,alpha2])],[n1,n2][numpy.argmax([alpha1,alpha2])],pF0]
    
    pFa=numpy.log10(1.0/SadPar[2])
    pF=numpy.log10(numpy.abs(h))
    b=0.1+(0.2/SadPar[3]**2)*(1.0-numpy.exp(-(thetar/(thetas-thetar))**2))
    
    Cap_AD=(h*numpy.log(10.0)*(pF-pF0))**(-1.0) * (1.0 - numpy.exp((pFa-pF)/b)/(1.0+numpy.exp((pFa-pF)/b)))
    
    dSdh= w1 * (-alpha1*n1*m1*(alpha1*h)**(n1-1.0)*(1+(alpha1*h)**n1)**(-m1-1.0)) \
        + w2 * (-alpha2*n2*m2*(alpha2*h)**(n2-1.0)*(1+(alpha2*h)**n2)**(-m2-1.0))
    
    Cap=(thetas-thetar)/(1.0-(w1*Vangenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*Vangenuchten4(h0,[0.0,1.0,alpha2,n2]))) * dSdh + thetar*Cap_AD
    
    return Cap


#==============================================================================
# Fitting
#==============================================================================

def FitRetentionModel0(Par,h,thetareal,model=Vangenuchten4):
    thetavg=model(h,Par)
#    Diff=((thetavg-thetareal)/thetareal)**2
    Diff=((thetavg-thetareal))**2
    return Diff.sum()

def FitRetentionModel(Par,h,thetareal,model=Vangenuchten4):
    thetavg=model(h,Par)
    Diff=((thetavg-thetareal)/thetareal)**2
#    Diff=((thetavg-thetareal))**2
    return Diff.sum()
    
    
def FitRetentionModel2(Par,h,thetareal,model=Vangenuchten4):
    thetavg=model(h,Par)
#    Diff=((numpy.log(thetavg)-numpy.log(thetareal))/(numpy.log(thetareal)))**2
    Diff=(numpy.log(thetavg)-numpy.log(thetareal))**2
    return Diff.sum()


def FitRetentionConductivityModel(Par,
                                  Retention_h,Retention_theta,Retention_ParI,Retention_model,
                                  Conductivity_h,Conductivity_K,Conductivity_ParI,Conductivity_model,
                                  ArbitraryExtraWeight=1e2 ):
    
    Retention_Par=Par[Retention_ParI]
    Retention_Pred=Retention_model(Retention_h,Retention_Par)
    Conductivity_Par=Par[Conductivity_ParI]
    Conductivity_Pred=Conductivity_model(Conductivity_h,Conductivity_Par)
    

#    Retention_Diff=((numpy.log(Retention_Pred)-numpy.log(Retention_theta))/(numpy.log(Retention_theta)))**2
#    Conductivity_Diff=((numpy.log(Conductivity_Pred)-numpy.log(Conductivity_K))/(numpy.log(Conductivity_K)))**2
    
#    Retention_Diff=(Retention_Pred-Retention_theta)**2
#    Conductivity_Diff=(numpy.log(Conductivity_Pred)-numpy.log(Conductivity_K))**2
    
#    Retention_Diff=1e3*((Retention_Pred-Retention_theta)/Retention_theta)**2
#    Conductivity_Diff=(numpy.log(Conductivity_Pred)-numpy.log(Conductivity_K))**2/numpy.std(numpy.log(Conductivity_K))
    
    Retention_Diff=ArbitraryExtraWeight*((Retention_Pred-Retention_theta)/Retention_theta)**2
    Conductivity_Diff=(numpy.log(Conductivity_Pred)-numpy.log(Conductivity_K))**2/numpy.std(numpy.log(Conductivity_K))
#    Conductivity_Diff=((numpy.log(Conductivity_Pred)-numpy.log(Conductivity_K))/(numpy.log(Conductivity_K)))**2
    
    Objective = Retention_Diff.mean() + Conductivity_Diff.mean()
    return Objective
