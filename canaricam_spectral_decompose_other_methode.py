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

#x_temp, temp1, etemp1=loadtxt('templates/ESO_244-G_012_NED02.dat', unpack=True)
#x_temp_sb=x_temp/(1.+0.022903)


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
x_temp ,ngc7714 = loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/NGC7714/spec_ngc7714_complet.dat', unpack=True)
temp1 = ngc7714*1000.
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
z_ngc7714 = 0.009333
x_temp_sb = x_temp/(1.+z_ngc7714)














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
#%%%%%%%%Model
def func(p, ccam_x):
  c1, c2, c3, c4, c5 = p
  return c1*ccam_x**(c2)*np.exp(-c3*f_ext) + c4*f_sb*np.exp(-c5*f_ext)

def err(p, ccam_y, ccam_x):
#   return (ccam_y - func(p, ccam_x))
  return (ccam_y - func(p, ccam_x))/ccam_ey

#p0 =[17., 1.4, 41., 0.2, 9.]

p0 =[10., 2., 40., 15., 18.]

#p0 = [8.9, 1.8, 41.6, 0.5, 9.3]

bounds = [(0,None), (0, 4.), (0, None), (0,None), (0,None)]

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

sigma = [sqrt(cov_x[0,0]), \
         sqrt(cov_x[1,1]), \
         sqrt(cov_x[2,2]), \
         sqrt(cov_x[3,3]), \
         sqrt(cov_x[4,4])]

print 'a1=',p[0],'+/-',sigma[0]
print 'a2=',p[1],'+/-',sigma[1]
print 'a3=',p[2],'+/-',sigma[2] 
print 'a4=',p[3],'+/-',sigma[3]
print 'a5=',p[4],'+/-',sigma[4]

print 'rchi2=', rchi2, ',','chi2=',chi2

print 'c1=', p[0]*c**(p[1])

# now plot the best fit curve and also +- 1 sigma curves
# (the square root of the diagonal covariance matrix  
# element is the uncertianty on the fit parameter.)
sigma = [sqrt(cov_x[0,0]), \
         sqrt(cov_x[1,1]), \
         sqrt(cov_x[2,2]), \
         sqrt(cov_x[3,3]), \
         sqrt(cov_x[4,4])]                       

y_cont = p[0]*ccam_x**(p[1])*np.exp(-p[2]*f_ext)
y_sb = p[3]*f_sb*np.exp(-p[4]*f_ext)
y_tot = y_cont + y_sb

#%%%%%%%%%%%%Ploting Fitting %%%%%%%%%
a1 = p[0]
a2 = p[1]
a3 = p[2]
a4 = p[3]
a5 = p[4]

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
delta_cont = sqrt((dfunc_c1**(2)*sig1**(2)) + (dfunc_c2**(2)*sig2**(2)) + (dfunc_c3**(2)*sig3**(2)))
delta_sb   = sqrt((dfunc_c4**(2)*sig4**(2)) + (dfunc_c5**(2)*sig5**(2)))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#file400 = open('continuum_CCAM_other_method_ESO244G012.dat','w+')
#print >>file400, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
#for i in range(0,len(ccam_x)):
#  print >>file400, ccam_x[i], y_cont[i], delta_cont[i]
#file400.close()  


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#file400 = open('continuum_CCAM_mod8.dat','w+')
#print >>file400, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
#for i in range(0,len(ccam_x)):
#  print >>file400, ccam_x[i], y_cont[i], delta_cont[i]
#file400.close()  

#file500 = open('starbust_CCAM_mod8.dat','w+')
#print >>file500, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
#for i in range(0,len(ccam_x)):
#  print >>file500, ccam_x[i], y_sb[i], delta_sb[i]
#file500.close() 


#file600 = open('composite_CCAM_mod8.dat','w+')
#print >>file600, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
#for i in range(0,len(ccam_x)):
#  print >>file600, ccam_x[i], y_tot[i], delta_func[i]
#file600.close() 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Errores
x_stad, std_comp = loadtxt('std_composite.dat', unpack=True)
std_comp = sqrt(std_comp**(2)+(ccam_ey)**(2))
x_stad, std_cont = loadtxt('std_continuo.dat', unpack=True)
std_cont = sqrt(std_cont**(2)+(ccam_ey**(2)))
x_stad, std_starb = loadtxt('std_starburst.dat', unpack=True)



file500 = open('spectrum_CCAM_other_method_ESO244G012_all_0_to_2p5.dat','w+')
print >>file500, '#lambda(um)Rest-Frame', 'f_nu(mJy)'
for i in range(0,len(ccam_x)):
  print >>file500, ccam_x[i]*(1.+zred), ccam_y[i]-y_sb[i], sqrt(std_starb[i]**(2)+ccam_ey[i]**(2))
file500.close()  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#for i in range(0, len(ccam_x)):
#  print ccam_x[i], ccam_y[i], ccam_ey[i]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plt.figure()
plt.subplot(211)

plt.plot(ccam_x, ccam_y, 'k-', linewidth=4, label="UGC5101 CC")
plt.plot(ccam_x, y_tot, 'r--', linewidth=4, label="Power law + Starburst template")

plt.plot(ccam_x, y_cont, 'g:', linewidth=4,label="Power law")


plt.fill_between(ccam_x, y_cont, y_cont+std_cont, color='LightGreen')
plt.fill_between(ccam_x, y_cont, y_cont-std_cont, color='LightGreen')
plt.plot(ccam_x, y_sb, 'b-.', linewidth=4,label="Starburst template")
plt.fill_between(ccam_x, y_tot, y_tot+std_comp, color='pink')
plt.fill_between(ccam_x, y_tot, y_tot-std_comp, color='pink')
plt.fill_between(ccam_x, ccam_y+ ccam_ey, ccam_y, color='grey')
plt.fill_between(ccam_x, ccam_y- ccam_ey, ccam_y, color='grey')
plt.fill_between(ccam_x, y_tot, y_tot+std_comp, color='pink')
plt.fill_between(ccam_x, y_tot, y_tot-std_comp, color='pink')

lccam=8.7
l12_5=12.5
lccam=lccam/(1.+zred)
l12_5 = l12_5/(1.+zred)

fccam=33.
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

nx,ny,ney=loadtxt('spectrum_CCAM_other_method_ESO244G012.dat', unpack=True)
plt.errorbar(nx/(1.+zred), ny, ney,fmt='o')

setp(gca(), xticks=[])
plt.ylabel(r"$f_{\nu}$ (mJy)",fontsize='26')

plt.legend(loc = 'best', numpoints=1, prop={'size':16})

plt.xlim(7.5, 12.3)

plt.figure()


plt.plot(ccam_x, ccam_y, 'k-', linewidth=4, label="UGC5101 CC")
plt.plot(ccam_x, y_tot, 'r--', linewidth=4, label="Power law + Starburst template")

plt.plot(ccam_x, y_cont, 'g:', linewidth=4,label="Power law")


plt.fill_between(ccam_x, y_cont, y_cont+std_cont, color='LightGreen')
plt.fill_between(ccam_x, y_cont, y_cont-std_cont, color='LightGreen')
plt.plot(ccam_x, y_sb, 'b-.', linewidth=4,label="Starburst template")
plt.fill_between(ccam_x, y_tot, y_tot+std_comp, color='pink')
plt.fill_between(ccam_x, y_tot, y_tot-std_comp, color='pink')
plt.fill_between(ccam_x, ccam_y+ ccam_ey, ccam_y, color='grey')
plt.fill_between(ccam_x, ccam_y- ccam_ey, ccam_y, color='grey')
plt.fill_between(ccam_x, y_tot, y_tot+std_comp, color='pink')
plt.fill_between(ccam_x, y_tot, y_tot-std_comp, color='pink')

lccam=8.7
l12_5=12.5
lccam=lccam/(1.+zred)
l12_5 = l12_5/(1.+zred)

fccam=33.
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

