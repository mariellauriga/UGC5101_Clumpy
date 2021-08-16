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
#x_ext, A_gc = loadtxt('extinction_chiar_tielens2006_original.dat',unpack=True)
#tau_gc = A_gc*0.09/(1.086)
x_ext, A_gc, A_ism = loadtxt('extinction_chiar_tielens2006_original_plus_ism.dat',unpack=True)
tau_gc = A_ism*0.09/(1.086)

#_____________________________________________
#%%%%%%%%%% The best Starburts Template%%%%%%%%%%%%%%%
#--------------------------------------------
#x_temp_sb, ngc1614, ngc2369, ngc3256, ngc4194, iras12112, iras14348, iras17208, iras22491, Arp220, ESO0320, Zw049=loadtxt('Rieke_individuals_templates_original.dat', unpack=True)
#temp1 = ngc1614*1000.

x_temp, temp1, etemp1=loadtxt('templates/ESO_244-G_012_NED02.dat', unpack=True)
x_temp_sb=x_temp/(1.+0.022903)

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
#%%%%%%%%%% Water absorbcion %%%%%%%%%%%%%%%
#--------------------------------------------
x_lamb, tau_hielo = loadtxt('h2op30k_new2.dat',unpack=True)
x_lamb, tau_hielo = loadtxt('h2op30k_new2.dat',unpack=True)
tau_hielo = log(10)*tau_hielo
#%%%%%%%%%% Water absorbcion %%%%%%%%%%%%%%%
#--------------------------------------------
#x_lamb, tau_hielo = loadtxt('h2op30k_new2.dat',unpack=True)
#x_lamb, tau_hielo = loadtxt('h2op30k_new2.dat',unpack=True)
#tau_hielo = log(10)*tau_hielo


#%%%%%%Interpoling Water absorbcion%%%%%
#---------------------------
#plt.figure()
#print 'Interpolating Water absorbcion'
#perfil  = interp1d(x_lamb, tau_hielo, kind='cubic')
#t_hielo = perfil(ccam_x)

#file500 = open('h2op30k_new_ccam.dat','w+')
#print >>file500, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
#for i in range(0,len(ccam_x)):
#  print >>file500, ccam_x[i], t_hielo[i]
#file500.close() 
#
x_lamb, t_hielo = loadtxt('h2op30k_new_ccam.dat',unpack=True)
plt.figure()
plt.plot(x_lamb, t_hielo, label='water profile')

#______________________________
#%%%%%%%%Model
def func(p, ccam_x):
  c1, c2, c3, c4, c5, c6, c7 = p
  return c1*ccam_x**(c2)*np.exp(-c3*f_ext)*np.exp(-c6*t_hielo) + c4*f_sb*np.exp(-c5*f_ext)*np.exp(-c7*t_hielo)

def err(p, ccam_y, ccam_x):
  return (ccam_y - func(p, ccam_x))/ccam_ey

#p0 =[17., 1.4, 41., 0.2, 9.]

p0 =[1., 2., 60., 1., 80., 1., 40.]
bounds = [(0,None), (0, 3.5), (0, None), (0,None), (0,None), (0,None), (0,None)]

#p0 =[1., 2., 60., 1., 18., 1., 40.]
#bounds = [(0,None), (0, None), (0, None), (0,None), (0,None), (0,None), (0,None)]


p, cov_x, infodic, mesg, ier = leastsqbound(err, p0, args=(ccam_y, ccam_x), bounds = bounds, full_output=True)

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
            
sigma_obs_spec = ccam_ey
if sigma_obs_spec is None:
 chi2 = sum((func(p, ccam_x)-ccam_y)**2)
else:
 chi2 = sum(((func(p, ccam_x)-ccam_y)/sigma_obs_spec)**2)
dof = len(ccam_y) - len(p)
rchi2 = chi2/dof

# now plot the best fit curve and also +- 1 sigma curves
# (the square root of the diagonal covariance matrix  
# element is the uncertianty on the fit parameter.)

sigma = [sqrt(cov_x[0,0]), \
         sqrt(cov_x[1,1]), \
         sqrt(cov_x[2,2]), \
         sqrt(cov_x[3,3]), \
         sqrt(cov_x[4,4]), \
         sqrt(cov_x[5,5]), \
         sqrt(cov_x[6,6])]

print 'a1=',p[0],'+/-',sigma[0]
print 'a2=',p[1],'+/-',sigma[1]
print 'a3=',p[2],'+/-',sigma[2] 
print 'a4=',p[3],'+/-',sigma[3]
print 'a5=',p[4],'+/-',sigma[4]
print 'a6=',p[5],'+/-',sigma[5]
print 'a7=',p[6],'+/-',sigma[6]

print 'rchi2=', rchi2, ',','chi2=',chi2

print 'c1=', p[0]*c**(p[1])

y_cont = p[0]*ccam_x**(p[1])*np.exp(-p[2]*f_ext)*np.exp(-p[5]*t_hielo)
y_sb = p[3]*f_sb*np.exp(-p[4]*f_ext)*np.exp(-p[6]*t_hielo)
y_tot = y_cont + y_sb

#%%%%%%%%%%%%Ploting Fitting %%%%%%%%%
a1 = p[0]
a2 = p[1]
a3 = p[2]
a4 = p[3]
a5 = p[4]
a6 = p[5]
a7 = p[6]

sig1 = sigma[0]
sig2 = sigma[1]
sig3 = sigma[2]
sig4 = sigma[3]
sig5 = sigma[4]
sig6 = sigma[5]
sig7 = sigma[6]


#Errors propagation
dfunc_c1 = ccam_x**(a2)*np.exp(-a3*f_ext)*np.exp(-a6*t_hielo)
dfunc_c2 = a1*ccam_x**(a2)*log(ccam_x)*np.exp(-a3*f_ext)*np.exp(-a6*t_hielo)
dfunc_c3 = -f_ext*a1*ccam_x**(a2)*np.exp(-a3*f_ext)*np.exp(-a6*t_hielo)
dfunc_c6 = -t_hielo*a1*ccam_x**(a2)*np.exp(-a3*f_ext)*np.exp(-a6*t_hielo)
dfunc_c4 = f_sb*np.exp(-a5*f_ext)*np.exp(-a7*t_hielo)
dfunc_c5 = -f_ext*a4*f_sb*np.exp(-a5*f_ext)*np.exp(-a7*t_hielo)
dfunc_c7 = -t_hielo*a4*f_sb*np.exp(-a5*f_ext)*np.exp(-a7*t_hielo)

delta_func = sqrt((dfunc_c1**(2)*sig1**(2))  + (dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)) + (dfunc_c4**(2)*sig4**(2)) + (dfunc_c5**(2)*sig5**(2)) + (dfunc_c6**(2)*sig6**(2)) + (dfunc_c7**(2)*sig7**(2)))
delta_cont = sqrt((dfunc_c1**(2)*sig1**(2)) + (dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2))+ (dfunc_c6**(2)*sig6**(2)))
delta_sb   = sqrt((dfunc_c4**(2)*sig4**(2)) + (dfunc_c5**(2)*sig5**(2)) + (dfunc_c7**(2)*sig7**(2)))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file400 = open('continuum_CCAM_mod3.dat','w+')
print >>file400, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
for i in range(0,len(ccam_x)):
  print >>file400, ccam_x[i], y_cont[i], delta_cont[i]
file400.close()  

file500 = open('starbust_CCAM_mod3.dat','w+')
print >>file500, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
for i in range(0,len(ccam_x)):
  print >>file500, ccam_x[i], y_sb[i], delta_sb[i]
file500.close() 


file600 = open('composite_CCAM_mod3.dat','w+')
print >>file600, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
for i in range(0,len(ccam_x)):
  print >>file600, ccam_x[i], y_tot[i], delta_func[i]
file600.close() 



file700 = open('spectrum_CCAM_other_method_ESO244G012_plus_water.dat','w+')
print >>file700, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
for i in range(0,len(ccam_x)):
  print >>file700, ccam_x[i]*(1.+zred), ccam_y[i]-y_sb[i], sqrt(delta_sb[i]**(2)+ccam_ey[i]**(2))
file700.close()  

#for i in range(0, len(ccam_x)):
#  print ccam_x[i], ccam_y[i], ccam_ey[i]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plt.figure()
plt.subplot(211)

plt.plot(ccam_x, ccam_y, 'k-', linewidth=4, label="UGC5101 CC")
plt.plot(ccam_x, y_tot, 'r--', linewidth=4, label="Power law + Starburst template")

plt.plot(ccam_x, y_cont, 'g:', linewidth=4,label="Power law")


plt.fill_between(ccam_x, y_cont, y_cont+delta_cont, color='LightGreen')
plt.fill_between(ccam_x, y_cont, y_cont-delta_cont, color='LightGreen')
plt.plot(ccam_x, y_sb, 'b-.', linewidth=4,label="Starburst template")

plt.fill_between(ccam_x, y_tot, y_tot+delta_func, color='pink')
plt.fill_between(ccam_x, y_tot, y_tot-delta_func, color='pink')
plt.fill_between(ccam_x, ccam_y+ ccam_ey, ccam_y, color='grey')
plt.fill_between(ccam_x, ccam_y- ccam_ey, ccam_y, color='grey')

lccam=8.7
l12_5=12.5
lccam=lccam/(1.+zred)
l12_5 = l12_5/(1.+zred)

fccam=36.
#Flujo de Soifer et al. 2000
f12_5 = 92.5

efccam= 4.
ef12_5= 5.

plt.errorbar(lccam,fccam, efccam, fmt='o',color='red',markersize=10, label='CCAM (PSF) present work')
plt.errorbar(l12_5, f12_5, ef12_5, fmt='o', color='orange',markersize=10, label='Soifer et al. 2000')


setp(gca(), xticks=[])
plt.ylabel(r"$f_{\nu}$ (mJy)",fontsize='26')

plt.legend(loc = 'best', numpoints=1, prop={'size':16})

plt.xlim(7.5, 12.3)
plt.subplot(212)

#w,h = figaspect(-2.)
plt.xlabel(" Rest-frame Wavelength ($\mu$m)",fontsize='26')
plt.ylabel(r"$Residual$ $(\%)$",fontsize='26')

plt.plot(ccam_x, ((y_tot-ccam_y)/ccam_y)*100.,'o')

 
x_ref = [0,30]
y_ref = [0, 0]

plt.plot(x_ref, y_ref, 'k--')

plt.xlim(7.5, 12.3)
tick_params(axis='x', labelsize=18)
tick_params(axis='y', labelsize=18)

plt.figure()

plt.plot(ccam_x, ccam_y, 'k-', linewidth=4, label="UGC5101 CC")
plt.plot(ccam_x, y_tot, 'r--', linewidth=4, label="Power law + Starburst template")

plt.plot(ccam_x, y_cont, 'g:', linewidth=4,label="Power law")


plt.fill_between(ccam_x, y_cont, y_cont+delta_cont, color='LightGreen')
plt.fill_between(ccam_x, y_cont, y_cont-delta_cont, color='LightGreen')
plt.plot(ccam_x, y_sb, 'b-.', linewidth=4,label="Starburst template")
plt.fill_between(ccam_x, ccam_y+ ccam_ey, ccam_y, color='grey')
plt.fill_between(ccam_x, ccam_y- ccam_ey, ccam_y, color='grey')
plt.fill_between(ccam_x, y_tot, y_tot+delta_func, color='pink')
plt.fill_between(ccam_x, y_tot, y_tot-delta_func, color='pink')

lccam=8.7
l12_5=12.5
lccam=lccam/(1.+zred)
l12_5 = l12_5/(1.+zred)

fccam=36.
#Flujo de Soifer et al. 2000
f12_5 = 92.5

efccam= 4.
ef12_5= 5.

plt.errorbar(lccam,fccam, efccam, fmt='o',color='red',markersize=10, label='CCAM (PSF) present work')
plt.errorbar(l12_5, f12_5, ef12_5, fmt='o', color='orange',markersize=10, label='Soifer et al. 2000')

setp(gca(), xticks=[])
plt.ylabel(r"$f_{\nu}$ (mJy)",fontsize='26')

plt.legend(loc = 'best', numpoints=1, prop={'size':16})

plt.xlim(7.5, 12.3)






plt.show()

