# Texture triangle demo for Lidong

import numpy
import pylab

from GeneralSoilPhys_Funky import *

Sand = numpy.array([0.39,0.32])
Clay = numpy.array([0.21,0.5])

TextureClass = GetTexture(Sand,Clay,DoPlot=True)[0]


pylab.show()
