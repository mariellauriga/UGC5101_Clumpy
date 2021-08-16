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
#from fitting import *

#______________________________
#%%%%%%%%%%Constants%%%%%%%%%%%
#------------------------------
zred = 0.039367
c=3.0e14

#_________________________________________________
#%%%%%%%%%%% Spectrum to be fitting%%%%%%%%%%%%%%%
#-------------------------------------------------
spit_x, spit_y, spitz_ey=loadtxt('../spitzer_spectrum/spitzer_spec.dat', unpack=True)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#_____________________________________________
#%%%%%%%%%% Sartburts Templates Library%%%%%%%%%%%%%%%
#--------------------------------------------
#x_temp, temp_9_75, temp10, temp10_25, temp10_5, temp10_75, temp11, temp11_25, temp11_5, temp11_75, temp_12, temp12_25, temp12_5, temp12_75, temp_13 =loadtxt('Rieke_template.dat', unpack = True)
#x_temp, temp1,  temp2,  temp3, temp4,lum= loadtxt('template_SBs_shmith.dat', unpack=True)
#x_temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11=loadtxt('Rieke_individuals_templates.dat', unpack=True)
#temp9 = 0.5*temp9
x_temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13=loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/individual_temp_Brandl_etal.dat', unpack=True)
#temp1 = temp1/x_temp2
x_temp_agn, temp1a, temp2a, temp3a, temp4a, temp5a = loadtxt('Assef_atal_template.dat', unpack = True)
x_temp2=c/x_temp
temp_agn2 = temp2a*(1.e-14)/(1.e-23)
#temp_agn2 = 1.e-15*temp_agn2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#_____________________________________________
#%%%%%%%%%% Extinction Laws Library%%%%%%%%%%%%%%%
#--------------------------------------------
x_ext, A_gc, A_ism = loadtxt('extinction_chiar_tielens2006.dat',unpack=True)
tau_gc = A_gc*0.09/(1.086)

#x_ext_d, alb, cos, Cext, Kabs, Cos2 = loadtxt('extinction_law_draine_R_3p1.dat', unpack=True)
#x_ext=range(1063)
#Cext=range(1063)
#x_ext[0:] = x_ext_d[:0]
#Cext[0:] = Cext_d[:0]
#A_gc=Cext[:0]*5.8e21
#tau_gc=A_gc/(1.086)

#_______________________________
#%%%%%% Interpoling Starburst Templates%%%%%
#-------------------------------
plt.figure()
print 'Interpolating Starbust'
f  = interp1d(x_temp, temp3, kind='cubic')
f1 = f(spit_x)
plt.plot(x_temp, temp3,'o', spit_x, f1)
plt.legend(['PAH template', 'cubic'], loc='best')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_______________________________
#%%%%%% Interpoling AGN Templates%%%%%
#-------------------------------
plt.figure()
print 'Interpolating ANG'
f  = interp1d(x_temp_agn, temp_agn2, kind='cubic')
f3 = f(spit_x)
plt.plot(x_temp_agn, temp_agn2,'o', spit_x, f3)
plt.legend(['AGN template', 'cubic'], loc='best')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#__________________________
#%%%%%%Interpoling Extintion Laws%%%%%
#---------------------------
plt.figure()
print 'Interpolating extinction law'
exf  = interp1d(x_ext, tau_gc, kind='cubic')
f2 = exf(spit_x)
plt.plot(x_ext, tau_gc,'o', spit_x, f2)
plt.legend(['Extinction', 'cubic'], loc='best')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#_____________________________________
#%%%%%%%%%%%%Fitting %%%%%%%%%
#----------------------------------
f30 = open('parameters.dat','w+')
print 'Fitting...'
#p0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
sigma = 0.1*spit_y
def func(xmodel, c1, c2, c3, c4):
  return (c1*f3*np.exp(-c2*f2))+(c3*f1*np.exp(-c4*f2))
popt, pcov = optimization.curve_fit(func, spit_x, spit_y, None, sigma)
if sigma is None:
 chi2 = sum(((func(spit_x,*popt)-spit_y))**2)
else:
 chi2 = sum(((func(spit_x,*popt)-spit_y)/sigma)**2)
dof = len(spit_y) - len(popt)
rchi2 = chi2/dof
print 'a1=',popt[0],',', 'a2=',popt[1],',', 'a3=',popt[2] 
print 'a4=',popt[3],',', 'rchi2=', rchi2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#_____________________________________
#%%%%%%%%%%%%Ploting Fitting %%%%%%%%%
#----------------------------------
z = func(spit_x, *popt)
plt.figure()
plt.errorbar(spit_x/(1.+zred), spit_y, 0.2*spit_y, color ='grey')
plt.plot(spit_x/(1.+zred), spit_y, 'ko', label="Spitzer Data UGC5101")
plt.plot(spit_x, z, 'r-', label="Fitted Curve ($\chi^{2}\sim2.7$)")


a1 = popt[0]
a2 = popt[1]
a3 = popt[2]
a4 = popt[3]

y_cont=a1*f3*np.exp(-a2*f2)
y_sb = a3*f1*np.exp(-a4*f2)

plt.yscale('log')
plt.plot(spit_x, y_cont, 'k--', label="Continuum template")
plt.plot(spit_x, y_sb, 'b-', label="PAH template (NGC2623)")
plt.legend()
xz=[5.,14.]
yz=[0.,0.]
plt.plot(xz, yz, color='black')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#chi2 = sum(((func(y_sint,*popt)-y_sint)/ey_sint)**2)
#dof = 180-2
#rchi2 = chi2/dof
#print >>f3, x_sint, popt[0], popt[1], rchi2
#f3.close()
#xpar, a1, a2, a_rchi2 =loadtxt('parameters.dat', unpack=True)
plt.xlabel("Rest-frame Wavelength ($\mu$m)",fontsize='20')
plt.ylabel(r"$f_{\nu}$ (Jy)",fontsize='20')

plt.show()