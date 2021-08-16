import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import scipy as scipy
import matplotlib.lines as mlines
import itertools
from matplotlib.patches import Ellipse, Polygon
import math
import scipy.optimize as optimization
from scipy.interpolate import interp1d
#______________________________



xccam, yccam, eyccam = loadtxt('continuum_CCAM_other_method.dat', unpack=True)
plt.plot(xccam, yccam, label='CCAM - power law')
#plt.errorbar(xccam, yccam, eyccam)

xspit, yspit, eyspit = loadtxt('continuum_spitzer_other_method_ngc1614.dat', unpack=True)
plt.plot(xspit, yspit, label='Spitzer - power law (NGC1614)')
#plt.errorbar(xspit, yspit, eyspit)

xspit0, yspit0, eyspit0 = loadtxt('continuum_spitzer_other_method_ESO244G012.dat', unpack=True)
plt.plot(xspit0, yspit0, label='Spitzer - power law (ESO244-G012)')
#plt.errorbar(xspit0, yspit0, eyspit0)

xspit1, yspit1 = loadtxt('continuum_spitzer_other_method_ngc1614_plus_water.dat', unpack=True)
plt.plot(xspit1, yspit1, label='NGC1614 + ISM + water (5-8.9)um')

xspit2, yspit2 = loadtxt('continuum_spitzer_other_method_eso244G012_plus_water.dat', unpack=True)
plt.plot(xspit2, yspit2, label='ESO244-G012 + ISM + water (5-8.9)um')

xspit3, yspit3 = loadtxt('continuum_spitzer_other_method_ngc1614_plus_water2.dat', unpack=True)
plt.plot(xspit3, yspit3, label='NGC1614 + ISM + water (2.5-22.2)um')

xspit4, yspit4 = loadtxt('continuum_spitzer_other_method_eso244G012_plus_water2.dat', unpack=True)
plt.plot(xspit4, yspit4, label='ESO244-G012 + ISM + water (2.5-22.2)um')

plt.errorbar(8.37, 33., 4., color='brown', fmt='o')

plt.legend(loc = 'best', numpoints=1, prop={'size':20})
plt.xlabel(" Rest-frame Wavelength ($\mu$m)",fontsize='26')
plt.ylabel(r"$log_{10}(f_{\nu})$ (mJy)",fontsize='26')


plt.show()
































plt.show()