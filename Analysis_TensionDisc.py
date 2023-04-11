#==============================================================================
# Tension Disc Infiltrometer
# Logsdon - Jaynes Analysis
#
# Input: tab-delimited file with 4 data columns
# - time (min)
# - water level in reservoir (cm) (should increase + the script can deal with refills)
# - tension head (cm)
# - reservoir area (m2)
# Data needs to be chronologically!
# + Transect area of the disc!
#
# Ghent University (2018)
# Jan De Pue (jan.depue@ugent.be)
# Give me some kick
# Before I loose
# my mind
#==============================================================================

import numpy
import pylab
from scipy import optimize

figlist = []


#==============================================================================
# Input Data
#==============================================================================


root = '/home/jan/git/SoilLab/TestData/TensionDisc'
filename = os.path.join(root,'Testdata_TensionDisc.txt')

Data = numpy.genfromtxt(filename,delimiter='\t',skip_header=1)

Time = Data[:,0]*60.0 # min => sec
Level = Data[:,1]*1e-2 # cm => m
Head = Data[:,2] # cm
Area = Data[:,3] # m2

DiscDiameter = 0.204 # m
DiscArea = numpy.pi * (DiscDiameter/2)**2


#==============================================================================
# Parse Data
#==============================================================================

HeadUnique = numpy.unique(Head)
nH = HeadUnique.size # number of tension heads

Time_L = []
Flux_L = []
SteadyFlux_L = numpy.zeros(nH)
for iH in range(nH):
    Head_h = HeadUnique[iH]
    Time_h = Time[Head == Head_h]
    Time_h = Time_h - Time_h[0]    # reset time
    Time_C = (Time_h[1:]+Time_h[:-1])/2.0
    dT = Time_h[1:] - Time_h[:-1]

    Volume = Level[Head == Head_h] * Area[Head == Head_h]
    dV = Volume[1:]-Volume[:-1]  # should be positive, unless refill happened

    FiltRefill = dV>0   ## remove refills
    dV = dV[FiltRefill]
    dT = dT[FiltRefill]
    Time_C = Time_C[FiltRefill]

    Discharge = dV/dT
    Flux = Discharge/DiscArea

    SteadyFlux = numpy.mean(Flux[-4:]) # mean of last 4 measurements (arbitrary!)

    # plot
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.plot(Time_C,Flux,'o-k')
    ax.axhline(SteadyFlux,color='0.3',ls='--')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Flux (m/s)')
    ax.set_title('Head = %1.2f cm'%Head_h)
    figlist.append(fig)
                
    Time_L.append(Time_C)
    Flux_L.append(Flux)
    SteadyFlux_L[iH] = SteadyFlux



#==============================================================================
# Logsdon and Jaynes (1993)
#==============================================================================

def GardnerExponential(Par,h):
    Ks,alpha = Par
    K = Ks * numpy.exp(alpha*h)
    return K

def LogsdonJaynes1993(Par,h,R):
    # ! Units alpha, h and R !
    Ks,alpha = Par
    q = Ks * numpy.exp(alpha*h) + (4* Ks * numpy.exp(alpha*h))/(numpy.pi*R*alpha)
    return q

def FitLogsdonJaynes(Par,h,R,qMeas):
    qMod = LogsdonJaynes1993(Par,h,R)
    Diff = numpy.mean((numpy.log10(qMod) - numpy.log10(qMeas))**2) ## better objective fuction
    # Diff = numpy.mean((qMod - qMeas)**2)
    return Diff

DiscRadius = DiscDiameter*0.5*1e2 # cm !
Par0=(1e-3,1e-1) # set appropriatly
bnd = ((1e-9,1e-2),(1e-4,1.0))
options = {'disp' : True,'maxiter' : 1000}
Min=optimize.minimize(FitLogsdonJaynes, Par0,
                      args=(HeadUnique,DiscRadius,SteadyFlux_L),
                      bounds = bnd,
                      tol= 1e-4,
                      # tol= 1e-20,
                      method='SLSQP',
                      options = options)
ParOpt = Min.x


hPlot = -numpy.linspace(0,20,100)
QLJ = LogsdonJaynes1993(ParOpt,hPlot,DiscRadius)

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(-HeadUnique,SteadyFlux_L,'ok')
ax.plot(-hPlot,QLJ,'-r')
ax.set_xlabel('Tension Head (cm)')
ax.set_ylabel('Flux (m/s)')
# ax.set_yscale('log')
figlist.append(fig)


hPlot = -numpy.logspace(-1,2,100)
KGardner = GardnerExponential(ParOpt,hPlot)

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(-hPlot,KGardner,'-k')
ax.set_xlabel('Tension Head (cm)')
ax.set_ylabel('Hydraulic Conductivity (m/s)')
ax.set_xscale('log')
figlist.append(fig)



#==============================================================================
# Write
#==============================================================================

print('Writing data'.center(50,'='))

from matplotlib.backends.backend_pdf import PdfPages

import inspect
ScriptBasename = os.path.basename(inspect.getfile(inspect.currentframe()))[:-3]
basename = ScriptBasename
Root = root

# save plots to pdf
pdfname = os.path.join(Root, basename + '' + '.pdf')
pp = PdfPages(pdfname)
for fig in figlist:
    pp.savefig(fig, bbox_inches='tight')
pp.close()

# # save plots to png
# extension='.eps'
# for iF in range(len(figlist)):
#     fig=figlist[iF]
#     figname = os.path.join(Root, basename+'_%s'%iF + extension)
#     fig.savefig(figname, bbox_inches='tight')

# pylab.show()

#os.startfile(pdfname)


