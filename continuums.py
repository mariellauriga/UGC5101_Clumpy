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

# Compara los continuos
zred = 0.039367

x_ccam, y_ccam, ey_ccam = loadtxt('continuum_CCAM_other_method_ESO244G012_c2_fixed.dat', unpack = True)
x_spitzer, y_spitzer, ey_spitzer = loadtxt('continuum_spitzer_other_method_ESO244_G12_plus_water.dat', unpack = True)
nx,ny,ney=loadtxt('spectrum_CCAM_other_method_ESO244G012_c2_fixed.dat', unpack=True)

plt.figure()

plt.plot(x_ccam, y_ccam, 'k--',color='orange', label="Continumm component (CCAM)")
plt.plot(x_spitzer, y_spitzer, 'k-.',color='green', label="Continumm component (Spitzer)")
plt.errorbar(nx/(1.+zred), ny, ney,fmt='o', label='Spectrum-continuum')

plt.fill_between(x_ccam, y_ccam + ey_ccam, y_ccam, color='Linen')
plt.fill_between(x_ccam, y_ccam - ey_ccam, y_ccam, color='Linen')

plt.fill_between(x_spitzer, y_spitzer + ey_spitzer, y_spitzer, color='LightGreen')
plt.fill_between(x_spitzer, y_spitzer - ey_spitzer, y_spitzer, color='LightGreen')

lccam=8.7
lccam=lccam/(1.+zred)
fccam=36.
efccam= 3.7

plt.errorbar(lccam,fccam, efccam, fmt='o',color='red', label='CCAM (PSF) present work')
plt.xlabel(" Rest-frame Wavelength ($\mu$m)",fontsize='20')
plt.ylabel(r"$f_{\nu}$ (mJy)",fontsize='20')
plt.legend(loc = 'best', numpoints=1)


plt.show()
