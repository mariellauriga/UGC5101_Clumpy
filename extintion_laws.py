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

zred = 0.039367
c=3.e14


x_lamb, y_flux=loadtxt('../spitzer_spectrum/spitzer_spec_completo.dat', unpack=True)
spit_y = y_flux*1000.
spit_x = x_lamb/(1.+zred)



x_ext, A_gc, A_ism = loadtxt('extinction_chiar_tielens2006_original.dat',unpack=True)
tau_ism = A_ism*0.09/(1.086)
tau_gc = A_gc*0.09/(1.086)


plt.figure()
print 'Interpolating extinction law'
exf  = interp1d(x_ext, tau_gc, kind='cubic')
f_gc = exf(spit_x)


exf_2  = interp1d(x_ext, tau_ism, kind='cubic')
f_ism = exf_2(spit_x)

plt.plot(x_ext, tau_gc, label='GC-extinction')
#plt.plot(spit_x, f_gc, label='cubic')


plt.plot(x_ext, tau_ism, label='ISM-extinction')
#plt.plot(spit_x, f_ism, label='cubic')


plt.legend( loc='best')
plt.show()