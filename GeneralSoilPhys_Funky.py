# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 13:11:50 2015

@author: Administrator
"""

#==============================================================================
# General Soil Physics functions
#==============================================================================

import numpy
import pylab

#==============================================================================
# Textuur
#==============================================================================

def DefineTexClassPolygons():
    
    #==============================================================================
    # define polygon
    #==============================================================================
    
    ## Draw polygons: Sand / Clay coordinates
    poly_U = numpy.array([[0.0  ,100.0],
                          [65.0 ,35.0 ],
                          [10.0 ,35.0 ],
                          [0.0  ,45.0 ],
                          [0.0  ,100.0]])
    
    poly_E = numpy.array([[65.0 ,35.0 ],
                          [10.0 ,35.0 ],
                          [0.0  ,45.0 ],
                          [0.0  ,30.0 ],
                          [25.0 ,17.5 ],
                          [82.5 ,17.5 ],
                          [65.0 ,35.0 ]])
    
    poly_A = numpy.array([[0.0  ,30.0 ],
                          [15.0 ,22.5 ],
                          [15.0 ,0.0  ],
                          [0.0  ,0.0  ],
                          [0.0  ,30.0 ]])
    
    poly_L = numpy.array([[15.0 ,22.5 ],
                          [25.0 ,17.5 ],
                          [67.5 ,17.5 ],
                          [67.5 ,12.5 ],
                          [50.0 ,12.5 ],
                          [50.0 ,0.0  ],
                          [15.0 ,0.0  ],
                          [15.0 ,22.5 ]])
    
    poly_P = numpy.array([[50.0 ,0.0  ],
                          [67.5 ,0.0  ],
                          [67.5 ,12.5 ],
                          [50.0 ,12.5 ],
                          [50.0 ,0.0  ]])
    
    poly_S = numpy.array([[67.5 ,0.0  ],
                          [67.5 ,17.5 ],
                          [82.5 ,17.5 ],
                          [92.5 ,7.5  ],
                          [82.5 ,7.5  ],
                          [82.5 ,0.0  ],
                          [67.5 ,0.0  ]])
    
    poly_Z = numpy.array([[82.5 ,0.0  ],
                          [82.5 ,7.5  ],
                          [92.5 ,7.5  ],
                          [100.0,0.0  ],
                          [82.5 ,0.0  ]])
    
    Poly_All_Labels=['U','E','A','L','P','S','Z']
    Poly_All=[poly_U,poly_E,poly_A,poly_L,poly_P,poly_S,poly_Z]
    Poly_All=map(lambda x : x/100.0,Poly_All)
    nP=len(Poly_All)
    
    Poly_Dict={}
    for iP in range(nP):
        Poly_Dict.update({Poly_All_Labels[iP] : Poly_All[iP]})
        
    return Poly_Dict


def DefineTexClassPolygons_USDA():
    
    #==============================================================================
    # define polygon
    #==============================================================================
    
    ## Draw polygons: Sand / Clay coordinates
    poly_Cl = numpy.array([[0.0 ,100.0],
                          [45.0  ,55.0 ],
                          [45.0 ,40.0 ],
                          [20.0  ,40.0 ],
                          [0.0  ,60.0],
                          [0.0 ,100.0],])
    
    poly_SaCl = numpy.array([[45.0 ,55.0 ],
                          [45.0 ,35.0 ],
                          [65.0  ,35.0 ],
                          [45.0 ,55.0 ],])
    
    poly_SiCl = numpy.array([[0.0 ,60.0 ],
                          [ 0.0 , 40.0],
                          [20.0 , 40.0],
                          [ 0.0 , 60.0],])
    
    poly_SaClL = numpy.array([[65.0 ,35.0 ],
                          [45.0 ,35.0 ],
                          [45.0 ,27.5 ],
                          [52.5 ,20.0 ],
                          [80.0 ,20.0 ],
                          [65.0 ,35.0 ],])
    
    poly_ClL = numpy.array([[45.0 ,40.0  ],
                          [20.0 , 40.0  ],
                          [20.0 , 27.5 ],
                          [45.0 , 27.5 ],
                          [45.0 ,40.0  ],])
    
    poly_SiClL = numpy.array([[20.0 ,40.0 ],
                          [ 0.0 , 40.0 ],
                          [ 0.0 , 27.5 ],
                          [20.0 , 27.5 ],
                          [20.0 ,40.0  ],])
    
    poly_SaL = numpy.array([[80.0 ,20.0  ],
                          [52.5 ,20.0  ],
                          [52.5 ,7.5 ],
                          [42.5 ,7.5 ],
                          [50.0 ,0.0  ],
                          [70.0 ,0.0  ],
                          [85.0 ,15.0  ],
                          [80.0 ,20.0  ],])

    poly_LSa = numpy.array([[85.0 ,15.0  ],
                          [70.0 ,0.0  ],
                          [85.0 ,0.0 ],
                          [90.0 ,10.0 ],
                          [85.0 ,15.0  ],])

    poly_Sa = numpy.array([[90.0 ,10.0  ],
                          [85.0 ,0.0  ],
                          [100.0 ,0.0 ],
                          [90.0 ,10.0  ],])

    poly_L = numpy.array([[45.0 ,27.5 ],
                          [22.5 ,27.5 ],
                          [42.5 ,7.5 ],
                          [52.5 ,7.5 ],
                          [52.5 ,20.0 ],
                          [45.0 ,27.5 ],])
    
    poly_SiL = numpy.array([[22.5 ,27.5 ],
                          [0.0 , 27.5 ],
                          [0.0 , 12.5 ],
                          [7.5 ,12.5  ],
                          [20.0 ,0.0  ],
                          [50.0 ,0.0  ],
                          [22.5 ,27.5 ],])
    
    poly_Si = numpy.array([[7.5 ,12.5 ],
                          [0.0 ,12.5 ],
                          [0.0 ,0.0  ],
                          [20.0 ,0.0 ],
                          [7.5 ,12.5 ],])
    
    Poly_All_Labels=['Cl','SaCl','SiCl','SaClL','ClL','SiClL','SaL','LSa','Sa','L','SiL','Si']
    Poly_All=[poly_Cl,poly_SaCl,poly_SiCl,poly_SaClL,poly_ClL,poly_SiClL,poly_SaL,poly_LSa,poly_Sa,poly_L,poly_SiL,poly_Si]
    Poly_All=[x/100.0 for x in Poly_All]
    nP=len(Poly_All)
    
    Poly_Dict={}
    for iP in range(nP):
        Poly_Dict.update({Poly_All_Labels[iP] : Poly_All[iP]})
        
    return Poly_Dict

def GetTexture(Sand,Clay,DoPlot=False,Type='Belgium'):
    
    import matplotlib.path as mplPath
    import matplotlib.patches as mplPatches

    if Type == 'Belgium':
      Poly_Dict=DefineTexClassPolygons()
    elif Type =='USDA':
      Poly_Dict=DefineTexClassPolygons_USDA()
    Poly_All_Labels=list(Poly_Dict.keys())
    nP=len(Poly_Dict)
    
    #==============================================================================
    # Get Texture    
    #==============================================================================
    
    returnfloat = False
    if isinstance(Sand,float):
        Sand=numpy.array([Sand,])
        Clay=numpy.array([Clay,])
        returnfloat = True
    
    nTex=Sand.size
    TexClass=numpy.zeros(nTex,dtype='U10')
    TexClass[:]='nan'
    
    if numpy.size(Sand)>0:
        for iP in range(nP):
            TexLabel=Poly_All_Labels[iP]
            bbPath = mplPath.Path(Poly_Dict[TexLabel])
            TexClass[bbPath.contains_points(numpy.array([Sand,Clay]).transpose())]=TexLabel
    else:
        TexClass=numpy.array(['nan',])

    if returnfloat:
       TexClass = TexClass[0]
    
    #==============================================================================
    # Draw
    #==============================================================================

    if DoPlot:
        
        fig=pylab.figure()
        ax=fig.add_subplot(111)
        
        for iP in range(nP):
            TexLabel=Poly_All_Labels[iP]
            bbPath = mplPath.Path(Poly_Dict[TexLabel])
            patch=mplPatches.PathPatch(bbPath,facecolor='w', lw=2)
        
            ax.add_patch(patch)
            ax.text(Poly_Dict[TexLabel].mean(axis=0)[0],Poly_Dict[TexLabel].mean(axis=0)[1],TexLabel,
                    horizontalalignment='center', verticalalignment='top')
        
        ax.plot(Sand,Clay,'r+')
        
        ax.set_aspect('equal')
        ax.set_xlabel('Sand (kg/kg)')
        ax.set_ylabel('Clay (kg/kg)')

        return TexClass,fig

    else:
        return TexClass

def DrawTexTriangle(Poly_Dict,fig=-1,ax=-1):

    import matplotlib.path as mplPath
    import matplotlib.patches as mplPatches
    
    if ax==-1:
        if fig==-1:
            fig=pylab.figure()
        ax=fig.add_subplot(111)
    
    Poly_All_Labels=Poly_Dict.keys()
    nP=len(Poly_Dict)
    
    for iP in range(nP):
        TexLabel=Poly_All_Labels[iP]
        bbPath = mplPath.Path(Poly_Dict[TexLabel])
        patch=mplPatches.PathPatch(bbPath,facecolor=None,fill=False, lw=2)
    
        ax.add_patch(patch)
        ax.text(Poly_Dict[TexLabel].mean(axis=0)[0],Poly_Dict[TexLabel].mean(axis=0)[1],TexLabel,
                horizontalalignment='center', verticalalignment='top')
    
    ax.set_aspect('equal')
    ax.set_xlabel('Sand (kg/kg)')
    ax.set_ylabel('Clay (kg/kg)')
    
    return fig,ax


def TexCode2Name(TexCode):
    TexCode=numpy.array(TexCode)
    TexName=numpy.zeros_like(TexCode).astype('S16')
    TexCode_unique=numpy.array(['U','E','A','L','P','S','Z'])
    TexName_unique=numpy.array(['Zware Klei','Klei','Leem','Zandleem','Lichte Zandleem','Lemig Zand','Zand'])
    
    nTex=len(TexCode_unique)
    
    for iT in range(nTex):
        TexName[TexCode==[TexCode_unique[iT]]]=TexName_unique[iT]
    
    return TexName
    
def TexCode2Name_USDA(TexCode):
    TexCode=numpy.array(TexCode)
    TexName=numpy.zeros_like(TexCode).astype('S16')
    TexCode_unique=numpy.array(['Cl','SaCl','SiCl','SaClL','ClL','SiClL','SaL','LSa','Sa','L','SiL','Si'])
    TexName_unique=numpy.array(['Clay','Sandy Clay','Silty Clay','Sandy Clay Loam','Clay Loam','Silty Clay Loam','Sandy Loam','Loamy Sand','Sand','Loam','Silt Loam','Silt'])
    
    nTex=len(TexCode_unique)
    
    for iT in range(nTex):
        TexName[TexCode==[TexCode_unique[iT]]]=TexName_unique[iT]
    
    return TexName
    
