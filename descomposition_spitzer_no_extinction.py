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

lPAH_7=[7.7,7.7]
lPAH_785=[7.85, 7.85]
lPAH_833=[8.33, 8.33]
lPAH_8=[8.61,8.61]
lSIV=[10.5,10.5]
lPAH=[11.3,11.3]
lPAH_12=[12.8,12.8]
lPAH_15=[15.56,15.56]
lPAH_1124=[11.24, 11.24]
lPAH_119=[11.99, 11.9]
lOIV=[25.89, 25.89]
f_ind=[10,20000]

#______________________________
#%%%%%%%%%%Constants%%%%%%%%%%%
#------------------------------
zred = 0.039367
c=3.e14
#_________________________________________________
#%%%%%%%%%%% Spectrum to be fitting%%%%%%%%%%%%%%%
#-------------------------------------------------
x_lamb, y_flux=loadtxt('../spitzer_spectrum/spitzer_spec_completo.dat', unpack=True)
spit_y = y_flux*1000.
spit_x = x_lamb/(1.+zred)

#spit_x, spit_y, spitz_ey=loadtxt('../spitzer_spectrum/spitzer_spec_corr_z.dat', unpack=True)
#spit_y=spit_y*1000
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_____________________________________________
#%%%%%%%%%% Sartburts Templates Library%%%%%%%%%%%%%%%
#--------------------------------------------
#x_temp, ngc1614, ngc2369, ngc3256, ngc4194, iras12112, iras14348, iras17208, iras22491, Arp220, ESO0320, Zw049=loadtxt('Rieke_individuals_templates_original.dat', unpack=True)
#temp1 = ESO0320*1000.

#-----Brand et al. 2006
#x_temp, mrk266 = loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/Mrk266/spec_mrk266_complet.dat', unpack=True)
#x_temp ,ngc520 = loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC520/spec_ngc520_complet.dat', unpack=True)
#x_temp ,ngc1097 =  loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC1097/spec_ngc1097_complet.dat', unpack=True)
#x_temp ,ngc1222 =  loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC1222/spec_ngc1222_complet.dat', unpack=True)
#x_temp ,ngc2623 =  loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC2623/spec_ngc2623_complet.dat', unpack=True)
#x_temp ,ngc3310 =  loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC3310/spec_ngc3310_complet.dat', unpack=True)
#x_temp ,ngc4088 =  loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC4088/spec_ngc4088_complet.dat', unpack=True)
#x_temp ,ngc4676 =  loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC4676/spec_ngc4676_complet.dat', unpack=True)
#x_temp ,ngc4818 =  loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC4818/spec_ngc4818_complet.dat', unpack=True)
#x_temp ,ngc7252 =  loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC7252/spec_ngc7252_complet.dat', unpack=True)
#x_temp ,ngc7714 = loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC7714/spec_ngc7714_complet.dat', unpack=True)
#temp1 = ngc7714*1000.
#z_mrk266 = 0.027863
#z_ngc520 = 0.007609
#z_ngc1097 = 0.004240
#z_ngc1222 = 0.008079
#z_ngc2623 = 0.018509
#z_ngc3310 = 0.003312
#z_ngc4088 = 0.002524 
#z_ngc4676 = 0.022049
#z_ngc4818 = 0.003552
#z_ngc7252 = 0.015984
#z_ngc7714 = 0.009333
#x_temp = x_temp/(1.+z_ngc7714)
#-----
#x_temp, m82 = loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/m82_sws_complete.dat', unpack = True)
#temp1 = m82*1000.
#zm82 = 0.000677
#x_temp = x_temp/(1.+zm82)
x_temp, temp1, etemp1=loadtxt('templates/ESO_244-G_012_NED02.dat', unpack=True)
x_temp=x_temp/(1.+0.022903)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_________________________________________________________________
#_____________________________________________
#%%%%%%%%%% Extinction Laws Library%%%%%%%%%%%%%%%
#--------------------------------------------
#x_ext, A_gc= loadtxt('extinction_chiar_tielens2006_original.dat',unpack=True)
x_ext, A_gc, A_ism= loadtxt('extinction_chiar_tielens2006_original_plus_ism.dat',unpack=True)
tau_gc = A_ism*0.09/(1.086)
#_______________________________
#%%%%%% Interpoling Starburst Templates%%%%%
#-------------------------------
plt.figure()
print 'Interpolating starbust'
f  = interp1d(x_temp, temp1, kind='cubic')
f1 = f(spit_x)
plt.plot(x_temp, temp1,'o', spit_x, f1)
plt.legend(['PAH template', 'cubic'], loc='best')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__________________________
#%%%%%%Interpoling Extintion Laws%%%%%
#---------------------------
#plt.figure()
print 'Interpolating extinction law'
exf  = interp1d(x_ext, tau_gc, kind='cubic')
f2 = exf(spit_x)


#plt.plot(x_ext, tau_gc,'o', spit_x, f2)
#plt.legend(['ISM-extinction', 'cubic'], loc='best')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_____________________________________
#%%%%%%%%%%%%Fitting %%%%%%%%%
#----------------------------------
f3 = open('parameters.dat','w+')
print 'Fitting...'
#p0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

sigma_obs_spec = 0.1*spit_y
#sigma_obs_spec = None
def func(xmodel, c2, c3, c4): 
  return (((xmodel)**(c2)*np.exp(-c3*f2)) + (c4*f1))
popt, pcov = optimization.curve_fit(func, spit_x, spit_y, None, sigma_obs_spec)
if sigma_obs_spec is None:
 chi2 = sum(((func(spit_x,*popt)-spit_y))**2)
else:
 chi2 = sum(((func(spit_x,*popt)-spit_y)/sigma_obs_spec)**2)
dof = len(spit_y) - len(popt)
rchi2 = chi2/dof


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_____________________________________
#%%%%%%%%%%%%Ploting Fitting %%%%%%%%%
#----------------------------------
zf = func(spit_x, popt[0], popt[1], popt[2])

a2 = popt[0]
a3 = popt[1]
a4 = popt[2]

y_cont = spit_x**(a2)*np.exp(-a3*f2) 

y_sb = a4*f1

y_tot = y_cont + y_sb

# now plot the best fit curve and also +- 1 sigma curves
# (the square root of the diagonal covariance matrix  
# element is the uncertianty on the fit parameter.)
sigma = [sqrt(pcov[0,0]), \
         sqrt(pcov[1,1]), \
         sqrt(pcov[2,2])]
#sig1 = sigma[0]
sig1 = 1.
sig2 = sigma[0]
sig3 = sigma[1]
sig4 = sigma[2]

#Errors propagation
dfunc_c2 = a2*spit_x**(a2-1)*np.exp(-a3*f2)
dfunc_c3 = -f2*spit_x**(a2)*np.exp(-a3*f2)
dfunc_c4 = f1

delta_func = sqrt((dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)) + (dfunc_c4**(2)*sig4**(2)))
delta_cont = sqrt((dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)))
delta_sb   = sqrt((dfunc_c4**(2)*sig4**(2)))


# a1 parameter and its error
#a1 = (y_tot-(a4*f1*np.exp(-a5*f2)))*c**(a2)*np.exp(a3*f2)/spit_x**(a2)
#da1 =c**(a2)*np.exp(a3*f2)*((log(c)*spit_x**(-a2))-(a2*spit_x**(-a2-1.)))*(y_tot-(a4*f1*np.exp(a5*f2)))
#da2 =(y_tot-(a4*f1*np.exp(-a5*f2)))*c**(a2)*spit_x**(-a2)*f2*np.exp(a3*f2)
#da3 =-f1*np.exp(-a5*f1)*c**(a2)*spit_x**(-a2)*np.exp(a3*f2)
#da4 =a4*f1*f2*np.exp(-a5*f2)*c**(a2)*spit_x**(-a2)*np.exp(a3*f2)
#da_tot = sqrt((da1*sig2)**(2.)+(da2*sig3)**(2.)+(da3*sig4)**(2.)+(da4*sig5)**(2.))

print 'Degrees of freddom=', dof
#print 'a1=',a1[0],'+/-',da_tot[0]
print 'a2=',popt[0],'+/-',sigma[0]
print 'a3=',popt[1],'+/-',sigma[1] 
print 'a4=',popt[2],'+/-',sigma[2]
print 'rchi2=', rchi2, ',','chi2=',chi2

nspit_x = c/spit_x

#y_cont2 = (a1*nspit_x**(-a2)*np.exp(-a3*f2))
#y_tot2 = y_cont2 + y_sb


#print 'Printing continuum - Spitzer'
#for i in range(0,len(spit_x)):
# print spit_x[i], y_cont[i], delta_cont[i]

#print 'Printing - Spitzer Spectrum in restframe'
#for i in range(0,len(spit_x)):
# print spit_x[i], spit_y[i], 0.2*spit_y[i]

plt.figure()

plt.subplot(211)

#plt.plot(lSIV,f_ind, 'k--',color='black')
#plt.text(9.6,800,"[SIV]", color='black', fontsize='16')
#plt.plot(lPAH,f_ind, 'k--',color='black')
#plt.text(10.5,800,"PAH", color='black', fontsize='16')
#plt.plot(lPAH_8,f_ind, 'k--',color='black')
#plt.text(7.8,800,"PAH", color='black', fontsize='16')
#plt.plot(lPAH_7,f_ind, 'k--',color='black')
#plt.text(6.9,800,"PAH", color='black', fontsize='16')
#plt.plot(lPAH_12, f_ind, 'k--', color='black')
#plt.text(12.0,800,"PAH", color='black', fontsize='16')


#plt.errorbar(spit_x, spit_y, 0.2*spit_y, color ='grey')
plt.fill_between(spit_x, log10(spit_y), log10(spit_y)+(0.1/log(10)), color='grey')
plt.fill_between(spit_x, log10(spit_y), log10(spit_y)-(0.1/log(10)), color='grey')
plt.plot(spit_x, log10(spit_y), 'k-',linewidth = 4, label="UGC5101 IRS")
#plt.errorbar(spit_x, log10(spit_y), 0.2/log(10), color='grey')
plt.plot(spit_x, log10(zf), 'r--', linewidth = 4, label="Power law + Starburst template")
#plt.errorbar(spit_x, log10(zf), delta_func/(zf*log(10)), color='pink')
plt.fill_between(spit_x, log10(zf), log10(zf)+delta_func/(zf*log(10)), color='pink')
plt.fill_between(spit_x, log10(zf), log10(zf)-delta_func/(zf*log(10)), color='pink')
plt.plot(spit_x, log10(y_cont), 'g:',linewidth = 4, label="Power law")
#plt.errorbar(spit_x, log10(y_cont), delta_cont/(y_cont*log(10)), color='LightGreen')
plt.plot(spit_x, log10(y_sb), 'b-.', linewidth = 4, label="Starburst template")
#plt.errorbar(spit_x, log10(y_sb), delta_sb/(y_sb*log(10)), color='LightBlue')

tick_params(axis='x', labelsize=18)
tick_params(axis='y', labelsize=18)

#plt.plot(spit_x,y_cont+y_sb, color='magenta')




#plt.yscale('log')
plt.xlim(5,26)
#plt.legend(loc = (0.01, 0.8), numpoints=1)
plt.legend(loc = 'best', numpoints=1, prop={'size':20})

plt.ylabel(r"$log_{10}(f_{\nu})$ (mJy)",fontsize='26')

setp(gca(), xticks=[])

plt.subplot(212)

#w,h = figaspect(-2.)
plt.xlabel(" Rest-frame Wavelength ($\mu$m)",fontsize='26')
plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')

plt.plot(spit_x, ((zf/spit_y)*100.),'o')

for i in range(0, len(spit_x)):
 print spit_y[i], zf[i]
 
x_ref = [0,30]
y_ref = [100, 100]

plt.plot(x_ref, y_ref, 'k--')

plt.xlim(5,26)
tick_params(axis='x', labelsize=18)
tick_params(axis='y', labelsize=18)
plt.show()