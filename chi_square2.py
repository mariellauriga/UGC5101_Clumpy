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

#%%%%%%%%%%Constants%%%%%%%%%%%
#------------------------------
zred = 0.039367

#Reading the x-values
ccam_x, ccam_y, ccam_ey = loadtxt('spec_ccam_ug5101_resampled_to_decompose.dat',unpack=True)

#Reading the starburst template
x_temp, ngc1614, ngc2369, ngc3256, ngc4194, iras12112, iras14348, iras17208, iras22491, Arp220, ESO0320, Zw049=loadtxt('Rieke_individuals_templates_original.dat', unpack=True)
temp1 = ngc1614*1000.
#print 'Interpolating starbust'
f  = interp1d(x_temp, temp1, kind='cubic')
f_sb = f(ccam_x)

#Reading the extinction law
x_ext, A_gc, A_ism = loadtxt('extinction_chiar_tielens2006_original_plus_ism.dat',unpack=True)
tau_gc = A_ism*0.09/(1.086)
#print 'Interpolating extinction law'
exf  = interp1d(x_ext, tau_gc, kind='cubic')
f_ext = exf(ccam_x)

#Calculating the model

xmodel = ccam_x

#c1 = 17.
#c2 = 1.4
#c3 = 41.
#c4 = 0.2
#c5 = 18.

numb, d1, d2, d3, d4, d5 = loadtxt('parameters3.dat', unpack=True)

k = 0
for l in range(10):
 c1 = d1[k:10+k]
 c2 = d2[k:10+k]
 c3 = d3[k:10+k]
 c4 = d4[k:10+k]
 c5 = d5[k:10+k]
 
 Fmodel = np.ones((66,10))
 y_cont = np.ones((66,10))
 y_sb = np.ones((66,10))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%Model's grid%%%%%%%%%%%%%%
 for i in range(10):
   for j in range(66):
     Fmodel[j,i] = c1[i]*xmodel[j]**(c2[i])*np.exp(-c3[i]*f_ext[j]) + c4[i]*f_sb[j]*np.exp(-c5[i]*f_ext[j])
     y_cont[j,i] = c1[i]*xmodel[j]**(c2[i])*np.exp(-c3[i]*f_ext[j])
     y_sb[j,i] = c4[i]*f_sb[j]*np.exp(-c5[i]*f_ext[j])
 file100 = open('grid_spectral_sampled_1.dat','w+')
 def printArray(args):
     print >>file100,"\t".join(args)


 for row in Fmodel:
  printArray([str(x) for x in row])
 file100.close() 

 mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10=loadtxt('grid_spectral_sampled_1.dat', unpack=True)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%Continuum's grid%%%%%%%%%%%%%
 file300 = open('grid_contiuum.dat','w+')
 def printArray(args):
     print >>file300,"\t".join(args)

 for row in y_cont:
  printArray([str(x) for x in row])
 file300.close() 

 cont1,cont2,cont3,cont4,cont5,cont6,cont7,cont8,cont9,cont10=loadtxt('grid_contiuum.dat', unpack=True)






#Calculating the chi-square
 dof = len(ccam_y) - 5

 chi2_1 = sum(((mod1-ccam_y)/ccam_ey)**2.)
 rchi2_1 = chi2_1/dof

 chi2_2 = sum(((mod2-ccam_y)/ccam_ey)**2.)
 rchi2_2 = chi2_2/dof

 chi2_3 = sum(((mod3-ccam_y)/ccam_ey)**2.)
 rchi2_3 = chi2_3/dof

 chi2_4 = sum(((mod4-ccam_y)/ccam_ey)**2.)
 rchi2_4 = chi2_4/dof

 chi2_5 = sum(((mod5-ccam_y)/ccam_ey)**2.)
 rchi2_5 = chi2_5/dof

 chi2_6 = sum(((mod6-ccam_y)/ccam_ey)**2.)
 rchi2_6 = chi2_6/dof

 chi2_7 = sum(((mod7-ccam_y)/ccam_ey)**2.)
 rchi2_7 = chi2_7/dof

 chi2_8 = sum(((mod8-ccam_y)/ccam_ey)**2.)
 rchi2_8 = chi2_8/dof

 chi2_9 = sum(((mod9-ccam_y)/ccam_ey)**2.)
 rchi2_9 = chi2_9/dof

 chi2_10 = sum(((mod10-ccam_y)/ccam_ey)**2.)
 rchi2_10 = chi2_10/dof

 file200 = open('grid_chi2_sampled_1.dat','w+')

# print 'chi2=', 'reduce-chi2'
 print chi2_1, rchi2_1
 print chi2_2, rchi2_2
 print chi2_3, rchi2_3
 print chi2_4, rchi2_4
 print chi2_5, rchi2_5
 print chi2_6, rchi2_6
 print chi2_7, rchi2_7
 print chi2_8, rchi2_7
 print chi2_9, rchi2_9
 print chi2_10, rchi2_10
 file200.close()


#Errors propagation
#dfunc_c1 = xmodel**(c2)*np.exp(-c3*f_ext)
#dfunc_c2 = c1*xmodel**(c2)*log(xmodel)*np.exp(-c3*f_ext)
#dfunc_c3 = -f_ext*c1*xmodel**(c2)*np.exp(-c3*f_ext)
#dfunc_c4 = f_sb*np.exp(-c5*f_ext)
#dfunc_c5 = -f_ext*c4*f_sb*np.exp(-c5*f_ext)

#delta_func = sqrt((dfunc_c1**(2)*sig1**(2) + dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)) + (dfunc_c4**(2)*sig4**(2)) + (dfunc_c5**(2)*sig5**(2)))
#delta_cont = sqrt((dfunc_c1**(2)*sig1**(2) + dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)))
#delta_sb   = sqrt((dfunc_c4**(2)*sig4**(2)) + (dfunc_c5**(2)*sig5**(2)))


 plt.figure('Model1')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod1, label='model1')
 plt.plot(xmodel, cont1, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
# plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod1/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')

 plt.figure('Model2')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod2, label='model2')
 plt.plot(xmodel, cont2, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
 #plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod2/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')

 plt.figure('Model3')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod3, label='model3')
 plt.plot(xmodel, cont3, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
 #plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod3/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')

 plt.figure('Model4')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod4, label='model4')
 plt.plot(xmodel, cont4, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
 #plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod4/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')


 plt.figure('Model5')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod5, label='model5')
 plt.plot(xmodel, cont5, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
#plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod5/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')


 plt.figure('Model6')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod6, label='model6')
 plt.plot(xmodel, cont6, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
 plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod6/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')


 plt.figure('Model7')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod7, label='model7')
 plt.plot(xmodel, cont7, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
 plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod7/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')


 plt.figure('Model8')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod8, label='model8')
 plt.plot(xmodel, cont8, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
 plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod8/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')

 plt.figure('Model9')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod9, label='model9')
 plt.plot(xmodel, cont9, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
 plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod9/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')


 plt.figure('Model10')
 plt.subplot(211)
 plt.plot(ccam_x, ccam_y, color='black')
 plt.plot(xmodel, mod10, label='model10')
 plt.plot(xmodel, cont10, 'k:')
#plt.plot(xmodel, y_sb, label='Starburst')
 plt.legend(loc = 'best', numpoints=1, prop={'size':20})
 plt.subplot(212)
 plt.xlabel(" Rest-frame Wavelength ($\mu$m)")
 plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')
 plt.plot(ccam_x, ((mod10/ccam_y)*100.),'o')
 x_ref = [7,12]
 y_ref = [100, 100]
 plt.plot(x_ref, y_ref, 'k--')
 
 k = k + 10

plt.show()