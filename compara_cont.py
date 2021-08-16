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
from leastsqbound import leastsqbound


x_nfix, y_nfix, ey_nfix =loadtxt('spectrum_CCAM_other_method_ESO244G012.dat', unpack=True)

x_fix, y_fix, ey_fix = loadtxt('spectrum_CCAM_other_method_ESO244G012_c2_fixed.dat', unpack=True)

plt.errorbar(x_nfix, y_nfix, ey_nfix,color='blue')
plt.plot(x_nfix, y_nfix, color='blue', label='C2 No-fixed')

plt.errorbar(x_fix, y_fix, ey_fix,color='red')
plt.plot(x_fix, y_fix, color='red', label='C2 fixed')


plt.legend(loc = 'best', numpoints=1, prop={'size':16})

plt.show()