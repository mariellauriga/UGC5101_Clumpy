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
zred = 0.0393671
c=3.e14
#_________________________________________________
#%%%%%%%%%%% Spectrum to be fitting%%%%%%%%%%%%%%%
#-------------------------------------------------
#%%%%%%%%%%%%%%%%CCAM spectrum%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ccam_x, ccam_y, ef, flag=loadtxt('Spec_m1_stck_UGC5101obs_215_216.dat',unpack=True)
#ccam_x, ccam_y, ef, flag=loadtxt('Spec_m1_stck_UGC5101obs_215_216_test1.dat',unpack=True) # A este le he quitado el punto negativo

ccam_x, ccam_y, ccam_ey = loadtxt('spec_ccam_ug5101_resampled_to_decompose.dat',unpack=True)

#Curva de transmision

#x2,f2cc,ef,flag=loadtxt('Spec_m1_stck_UGC5101obs_215_216.dat',unpack=True)
#fccam=(f1cc+f2cc)/2.
#ef3=((0.2*f1cc)**(2.)+(0.2*f2cc)**(2.))**(0.5)
#xsc, ysc,eysc=loadtxt('ccam_spec_scaled.dat', unpack=True)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%% Extinction Laws Library%%%%%%%%%%%%%%%
#--------------------------------------------
x_ext, A_gc = loadtxt('extinction_chiar_tielens2006_original.dat',unpack=True)
tau_gc = A_gc*0.09/(1.086)
#x_ext, A_gc, A_ism = loadtxt('extinction_chiar_tielens2006_original_plus_ism.dat',unpack=True)
#tau_gc = A_ism*0.09/(1.086)

#_____________________________________________
#%%%%%%%%%% The best Starburts Template%%%%%%%%%%%%%%%
#--------------------------------------------
x_temp_sb, ngc1614, ngc2369, ngc3256, ngc4194, iras12112, iras14348, iras17208, iras22491, Arp220, ESO0320, Zw049=loadtxt('Rieke_individuals_templates_original.dat', unpack=True)
temp1 = ngc4194*1000.

#x_temp, temp1, etemp1=loadtxt('templates/ESO_244-G_012_NED02.dat', unpack=True)
#x_temp_sb=x_temp/(1.+0.022903)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%Interpoling Extintion Laws%%%%%
#---------------------------
#plt.figure()
print 'Interpolating extinction law'
f3  = interp1d(x_ext, tau_gc, kind='cubic')
f_ext = f3(ccam_x)


plt.plot(x_ext, tau_gc,'o', ccam_x, f_ext)
plt.legend(['Extinction', 'cubic'], loc='best')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_______________________________
#%%%%%% Interpoling Starburst Templates%%%%%
#-------------------------------
plt.figure()
print 'Interpolating starbust'
f_1  = interp1d(x_temp_sb, temp1, kind='cubic')
f_sb = f_1(ccam_x)

plt.xlim(0, 20)
plt.legend(['PAH template', 'cubic'], loc='best')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#_____________________________________

#sigma_obs_spec = None
p0 = np.array([10, 2.0, 40., 0.2, 0.1])
#p0 = 20.
#p0 = None
sigma_obs_spec = ccam_ey
#%%%%%%%%Model 1 y 2
def func(xmodel, c1,c2, c3, c4, c5): 
  return (c1*xmodel**(c2)*np.exp(-c3*f_ext)) + (c4*f_sb*np.exp(-c5*f_ext))
popt, pcov = optimization.curve_fit(func, ccam_x, ccam_y, p0, sigma_obs_spec)
if sigma_obs_spec is None:
 chi2 = sum((func(ccam_x,*popt)-ccam_y)**2)
else:
 chi2 = sum(((func(ccam_x,*popt)-ccam_y)/sigma_obs_spec)**2)
dof = len(ccam_y) - len(popt)
rchi2 = chi2/dof

# now plot the best fit curve and also +- 1 sigma curves
# (the square root of the diagonal covariance matrix  
# element is the uncertianty on the fit parameter.)
sigma = [sqrt(pcov[0,0]), \
         sqrt(pcov[1,1]), \
         sqrt(pcov[2,2]), \
         sqrt(pcov[3,3]), \
         sqrt(pcov[4,4])]
                         
zf = func(ccam_x, *popt)

#%%%%%%%%%%%%Ploting Fitting %%%%%%%%%
a1 = popt[0]
a2 = popt[1]
a3 = popt[2]
a4 = popt[3]
a5 = popt[4]

sig1 = sigma[0]
sig2 = sigma[1]
sig3 = sigma[2]
sig4 = sigma[3]
sig5 = sigma[4]

#Errors propagation
dfunc_c1 = ccam_x**(a2)*np.exp(-a3*f_ext)
dfunc_c2 = a1*ccam_x**(a2)*log(ccam_x)*np.exp(-a3*f_ext)
dfunc_c3 = -f_ext*a1*ccam_x**(a2)*np.exp(-a3*f_ext)
dfunc_c4 = f_sb*np.exp(-a5*f_ext)
dfunc_c5 = -f_ext*a4*f_sb*np.exp(-a5*f_ext)

delta_func = sqrt((dfunc_c1**(2)*sig1**(2))  + (dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)) + (dfunc_c4**(2)*sig4**(2)) + (dfunc_c5**(2)*sig5**(2)))
delta_cont = sqrt((dfunc_c1**(2)*sig1**(2))  + (dfunc_c1**(2)*sig1**(2))  +(dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)))
delta_sb   = sqrt((dfunc_c4**(2)*sig4**(2)) + (dfunc_c5**(2)*sig5**(2)))


y_cont = a1*ccam_x**(a2)*np.exp(-a3*f_ext)
y_sb = a4*f_sb*np.exp(-a5*f_ext)
y_tot = y_cont + y_sb

print 'Degrees of freddom=', dof

#a1 = (y_tot-(a4*f_sb*np.exp(-a5*f_ext)))*c**(a2)*np.exp(a3*f_ext)/ccam_x**(a2)
#da1 = abs(c**(a2)*np.exp(a3*f_ext)*((log(c)*ccam_x**(-a2))-(a2*ccam_x**(-a2-1.)))*(y_tot-(a4*f_sb*np.exp(a5*f_ext))))
#da2 = abs(((y_tot)-(a4*f_sb*np.exp(-a5*f_ext)))*c**(a2)*ccam_x**(-a2)*f_ext*np.exp(a3*f_ext))
#da3 = abs(-f_sb*np.exp(-a5*f_sb)*c**(a2)*ccam_x**(-a2)*np.exp(a3*f_ext))
#da4 = abs(a4*f_sb*f_ext*np.exp(-a5*f_ext)*c**(a2)*ccam_x**(-a2)*np.exp(a3*f_ext))
#da_tot = sqrt((da1*sig2)**(2.)+(da2*sig3)**(2.)+(da3*sig4)**(2.)+(da4*sig5)**(2.))


print 'a1=',popt[0],'+/-',sigma[0]
print 'a2=',popt[1],'+/-',sigma[1]
print 'a3=',popt[2],'+/-',sigma[2] 
print 'a4=',popt[3],'+/-',sigma[3]
print 'a5=',popt[4],'+/-',sigma[4]
print 'rchi2=', rchi2, ',','chi2=',chi2
 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plt.figure()
plt.subplot(211)
#plt.fill([8.3/(1+zred),8.86/(1+zred),8.86/(1+zred),8.3/(1+zred)],[-0.2,-0.2,200,200], color='Gainsboro', fill=False, hatch='\\')
#plt.fill([9.36/(1+zred),10.04/(1+zred),10.04/(1+zred),9.36/(1+zred)],[-0.2,-0.2,200,200], color='Gainsboro',fill=False, hatch='\\')
#plt.fill([12.5,12.64,12.64,12.5],[-0.2,-0.2,2,2], color='Gainsboro',fill=False,hatch='\\')

plt.plot(ccam_x, ccam_y, 'k-', linewidth=4, label="UGC5101 CC")

plt.plot(ccam_x, zf, 'r--', linewidth=4, label="Power law + Starburst template")
plt.fill_between(ccam_x, zf, zf+delta_func, color='pink')
plt.fill_between(ccam_x, zf, zf-delta_func, color='pink')

plt.plot(ccam_x, y_cont, 'g:', linewidth=4,label="Power law")
plt.plot(ccam_x, y_sb, 'b-.', linewidth=4,label="Starburst template")
###plt.plot(ccam_x, y_cont+y_sb, color='magenta')
###plt.fill_between(ccam_x, y_sb, y_sb-y_sb_sigma, color='LightBlue')
####plt.fill_between(ccam_x, y_sb, y_sb+y_sb_sigma, color='LightBlue')


plt.fill_between(ccam_x, y_cont, y_cont+delta_cont, color='LightGreen')
plt.fill_between(ccam_x, y_cont, y_cont-delta_cont, color='LightGreen')

plt.fill_between(ccam_x, ccam_y+ ccam_ey, ccam_y, color='grey')
plt.fill_between(ccam_x, ccam_y- ccam_ey, ccam_y, color='grey')

nccam_x = c/ccam_x
y_cont2 = (a1*nccam_x**(-a2)*np.exp(-a3*f_ext))
y_tot2 = y_cont2 + y_sb

#print 'Printing continuum - CCAM'
#for i in range(0,len(ccam_x)):
#  print ccam_x[i], nccam_x[i], y_cont[i], delta_cont[i]
  
#print 'Printing rest-frame - CCAM-spectrum'
#for i in range(0,len(ccam_x)):
#  print ccam_x[i], ccam_y[i]

lccam=8.7
l12_5=12.5
lccam=lccam/(1.+zred)
l12_5 = l12_5/(1.+zred)

fccam=36.
#Flujo de Soifer et al. 2000
f12_5 = 92.5

efccam= 3.7
ef12_5= 5.

plt.errorbar(lccam,fccam, efccam, fmt='o',color='red',markersize=10, label='CCAM (PSF) present work')
plt.errorbar(l12_5, f12_5, ef12_5, fmt='o', color='orange',markersize=10, label='Soifer et al. 2000')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu_ccam_comp = ccam_y - y_sb
nerr = ccam_ey+delta_sb


#plt.errorbar(ccam_x, y_cont, sigm_cont, color='magenta')
#print 'Nuclear component'
#for i in range(0, len(ccam_x)):
#  print  nccam_x[i], ccam_x[i], nu_ccam_comp[i], (nu_ccam_comp[i]-abs(nu_ccam_comp_minus[i])+ nu_ccam_comp_plus[i] -nu_ccam_comp[i])/2

nerr= ccam_ey*nu_ccam_comp/ccam_y

file300 = open('nuclear_component_sampled_to_bayesclumpy.dat','w+')
print >>file300, '#lambda(um)OBS', 'f_nu(mJy)'
for i in range(0,len(ccam_x)):
  print >>file300, ccam_x[i]*(1.+ zred), nu_ccam_comp[i], nerr[i]
file300.close()                                                      

tick_params(axis='x', labelsize=18)
tick_params(axis='y', labelsize=18)

file400 = open('continuum_component_cc.dat','w+')
print >>file400, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
for i in range(0,len(ccam_x)):
  print >>file400, ccam_x[i], y_cont[i], delta_cont[i]
file400.close()

plt.ylim(0,200)
plt.xlim(7.5, 12.3)
#plt.xscale('log')
#plt.xlabel(" Rest-frame Wavelength ($\mu$m)",fontsize='26')
setp(gca(), xticks=[])
plt.ylabel(r"$f_{\nu}$ (mJy)",fontsize='26')

plt.legend(loc = 'best', numpoints=1, prop={'size':16})


plt.subplot(212)

#w,h = figaspect(-2.)
plt.xlabel(" Rest-frame Wavelength ($\mu$m)",fontsize='26')
plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')

plt.plot(ccam_x, ((zf/ccam_y)*100.),'o')

 
x_ref = [0,30]
y_ref = [100, 100]

plt.plot(x_ref, y_ref, 'k--')

plt.xlim(7.5, 12.3)
tick_params(axis='x', labelsize=18)
tick_params(axis='y', labelsize=18)


plt.figure()

plt.errorbar(ccam_x, nu_ccam_comp, nerr, color='pink')
plt.plot(ccam_x, nu_ccam_comp, color='red')
plt.plot(ccam_x, y_cont, 'g:', linewidth=4,label="Continuum")
plt.fill_between(ccam_x, y_cont, y_cont+delta_cont, color='LightGreen')
plt.fill_between(ccam_x, y_cont, y_cont-delta_cont, color='LightGreen')


plt.ylim(0,200)
plt.xlim(7.5, 12.5)
#plt.xscale('log')
plt.xlabel(" Rest-frame Wavelength ($\mu$m)",fontsize='26')
plt.ylabel(r"$f_{\nu}$ (mJy)",fontsize='26')

plt.legend(loc = 'best', numpoints=1, prop={'size':20})


plt.figure()

plt.plot(ccam_x, f_sb, label='before')
plt.plot(ccam_x, y_sb, 'b-.', linewidth=4,label="Starburst template")


plt.show()


