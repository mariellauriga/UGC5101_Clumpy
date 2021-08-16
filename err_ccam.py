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
#%%%%%%%%%%Constants%%%%%%%%%%%
#------------------------------
###Ya no es necesario usar este programa
zred = 0.0393671
c=3.e14

ccam_x, ccam_y, ef, flag=loadtxt('Spec_m1_stck_UGC5101obs_215_216.dat',unpack=True) # A este le he quitado el punto negativo
ccam_x =(ccam_x/10000.)/(1.+zred)
ccam_y = (ccam_y*1000.)
ef = (ef*1000.)

def qmean(ccam_y):
    return sqrt(sum(n*n for n in ccam_y[272:275])/len(ccam_y))
err1=qmean(ccam_y)  

ccam_ey = sqrt(err1**(2.) + (0.1*ef[272:275])**(2.))

for i in range(0,4):
 print ccam_ey[i]