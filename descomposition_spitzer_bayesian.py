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


#______________________________
#%%%%%%%%%%Constants%%%%%%%%%%%
#------------------------------
zred = 0.039367
c=3.e14 # um s^{-1}
#_________________________________________________
#%%%%%%%%%%% Spectrum to be fitting%%%%%%%%%%%%%%%
#-------------------------------------------------
x_lamb, y_flux=loadtxt('../spitzer_spectrum/spitzer_spec_completo_to_ism.dat', unpack=True)
spit_y = y_flux*1000.
spit_x = x_lamb/(1.+zred)

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
x_ext, A_gc, A_ism = loadtxt('extinction_chiar_tielens2006_original_plus_ism.dat',unpack=True)
tau_gc = A_ism*0.09/(1.086)

#x_ext, A_gc= loadtxt('extinction_chiar_tielens2006_original.dat',unpack=True)
#tau_gc = A_gc*0.09/(1.086)

#_______________________________
#%%%%%% Interpoling Starburst Templates%%%%%
#-------------------------------
plt.figure()
print 'Interpolating starbust'
f  = interp1d(x_temp, temp1, kind='cubic')
f1 = f(spit_x)
plt.plot(x_temp, temp1, label='PAH template')
plt.plot(spit_x, f1, label='cubic')
plt.legend( loc='best')
plt.xlim(5,35)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__________________________
#%%%%%%Interpoling Extintion Laws%%%%%
#---------------------------
plt.figure()
print 'Interpolating extinction law'
exf  = interp1d(x_ext, tau_gc, kind='cubic')
f2 = exf(spit_x)

plt.plot(x_ext, tau_gc,'o', spit_x, f2)
plt.legend(['ISM-extinction', 'cubic'], loc='best')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_____________________________________
#%%%%%%%%%%%%Fitting %%%%%%%%%
#----------------------------------

sigma_obs_spec = 0.1*spit_y

#spit_y= spit_y - sigma_obs_spec

def func(p, spit_x): 
  c1, c2, c3, c4, c5 = p
  return c1*spit_x**(c2)*np.exp(-c3*f2) + c4*f1*np.exp(-c5*f2)

def err(p, spit_y, spit_x):
  return (spit_y - func(p, spit_x))/sigma_obs_spec
p0 =[1., 2., 22., 0.2, 18.]
#p0 =[10., 2., 40., 0.3, 10.]
#p0 =[1., 1., 1., 1., 1.]
bounds = [(0,None), (0, None), (0,None), (0,None), (0,None)]

p, cov_x, infodic, mesg, ier = leastsqbound(err, p0, args=(spit_y, spit_x), bounds = bounds, full_output=True)

print "Bounded Least Squares fitting results:"
print "p:", p
print "cov_x:", cov_x
#print "infodic['nfev']:", infodic['nfev']
#print "infodic['fvec']:", infodic['fvec']
#print "infodic['fjac']:", infodic['fjac']
#print "infodic['ipvt']:", infodic['ipvt']
#print "infodic['qtf']:", infodic['qtf']
print "mesg:", mesg
print "ier:", ier


if sigma_obs_spec is None:
 chi2 = sum(((func(p, spit_x)-spit_y))**2)
else:
 chi2 = sum(((func(p, spit_x)-spit_y)/sigma_obs_spec)**2)
dof = len(spit_y) - len(p)
rchi2 = chi2/dof


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_____________________________________
#%%%%%%%%%%%%Ploting Fitting %%%%%%%%%
#----------------------------------

a1 = p[0]*c**(p[1])
a2 = p[1]
a3 = p[2]
a4 = p[3]
a5 = p[4]

nspit_x = c/spit_x

y_cont = a1*nspit_x**(-a2)*np.exp(-a3*f2) 

y_sb = a4*f1*np.exp(-a5*f2)

y_tot = y_cont + y_sb

# now plot the best fit curve and also +- 1 sigma curves
# (the square root of the diagonal covariance matrix  
# element is the uncertianty on the fit parameter.)
sigma = [sqrt(cov_x[0,0]), \
         sqrt(cov_x[1,1]), \
         sqrt(cov_x[2,2]), \
         sqrt(cov_x[3,3]), \
         sqrt(cov_x[4,4])]
#sig1 = sigma[0]
sig1 = c**(p[1])*sigma[0]
sig2 = sigma[1]
sig3 = sigma[2]
sig4 = sigma[3]
sig5 = sigma[4]


#Errors propagation
dfunc_c1 = nspit_x**(-a2)*np.exp(-a3*f2)
dfunc_c2 = -a1*np.exp(-a3*f2)*log(nspit_x)/(nspit_x**(a2))
dfunc_c3 = -f2*a1*nspit_x**(-a2)*np.exp(-a3*f2)
dfunc_c4 = f1*np.exp(-a5*f2)
dfunc_c5 = -f2*a4*f1*np.exp(-a5*f2)

delta_func = sqrt((dfunc_c1**(2)*sig1**(2) + dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)) + (dfunc_c4**(2)*sig4**(2)) + (dfunc_c5**(2)*sig5**(2)))
delta_cont = sqrt((dfunc_c1**(2)*sig1**(2) + dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)))
delta_sb   = sqrt((dfunc_c4**(2)*sig4**(2)) + (dfunc_c5**(2)*sig5**(2)))


print 'Degrees of freddom=', dof
print 'a1=',a1,'+/-',sig1
print 'a2=',p[1],'+/-',sigma[1]
print 'a3=',p[2],'+/-',sigma[2] 
print 'a4=',p[3],'+/-',sigma[3]
print 'a5=',p[4],'+/-',sigma[4]
print 'rchi2=', rchi2, ',','chi2=',chi2

print 'c1=', p[0]*c**(p[1])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file400 = open('continuum_spitzer_other_method_ESO244G012.dat','w+')
print >>file400, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
for i in range(0,len(spit_x)):
  print >>file400, spit_x[i], y_cont[i], delta_cont[i]
file400.close()  

file500 = open('spitzer_components_upper.dat','w+')
print >>file500, '#lambda(um)Rest-Frame', 'f_nu(mJy)Total', 'F_nu(mJy) Cont.', 'F_nu(mJy) SB'
for i in range(0,len(spit_x)):
  print >>file500, spit_x[i], y_tot[i], y_cont[i], y_sb[i]
file500.close()  




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plt.figure()

plt.subplot(211)


plt.plot(spit_x, log10(spit_y), 'k-',linewidth = 4, label="UGC5101 IRS")

plt.plot(spit_x, log10(y_tot), 'r--', linewidth = 4, label="Power law + Starburst template")

#plt.fill_between(spit_x, log10(y_tot), log10(y_tot)+(delta_func/(y_tot*log(10))), color='pink')
#plt.fill_between(spit_x, log10(y_tot), log10(y_tot)-(delta_func/(y_tot*log(10))), color='pink')
#plt.fill_between(spit_x, log10(spit_y), log10(spit_y)+(0.1/log(10)), color='grey')
#plt.fill_between(spit_x, log10(spit_y), log10(spit_y)-(0.1/log(10)), color='grey')

plt.fill_between(spit_x, log10(y_tot), log10(y_tot+delta_func), color='pink')
plt.fill_between(spit_x, log10(y_tot), log10(y_tot-delta_func), color='pink')
plt.fill_between(spit_x, log10(spit_y), log10(spit_y+0.1*spit_y), color='grey')
plt.fill_between(spit_x, log10(spit_y), log10(spit_y-0.1*spit_y), color='grey')



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

plt.plot(spit_x, ((y_tot-spit_y)/spit_y)*100.,'o')

#for i in range(0, len(spit_x)):
# print spit_y[i], zf[i]
 
x_ref = [0,30]
y_ref = [0, 0]

plt.plot(x_ref, y_ref, 'k--')

plt.xlim(5,26)
tick_params(axis='x', labelsize=18)
tick_params(axis='y', labelsize=18)


plt.show()