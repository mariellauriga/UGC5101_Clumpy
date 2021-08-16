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
f_ind=[-100,800]

#______________________________
#%%%%%%%%%%Constants%%%%%%%%%%%
#------------------------------
zred = 0.039367
c=3.e14
#_________________________________________________
#%%%%%%%%%%% Spectrum to be fitting%%%%%%%%%%%%%%%
#-------------------------------------------------
spit_x, spit_y, spitz_ey=loadtxt('../spitzer_spectrum/spitzer_spec_corr_z.dat', unpack=True)
spit_y=spit_y*1000


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_____________________________________________
#%%%%%%%%%% Sartburts Templates Library%%%%%%%%%%%%%%%
#--------------------------------------------
#x_temp, temp_9_75, temp10, temp10_25, temp10_5, temp10_75, temp11, temp11_25, temp11_5, temp11_75, temp_12, temp12_25, temp12_5, temp12_75, temp_13 =loadtxt('Rieke_template.dat', unpack = True)
#temp_11=temp_12*1000
#x_temp, temp1,  temp2,  temp3, temp4,lum= loadtxt('template_SBs_shmith.dat', unpack=True)
#lambda 
x_temp, ngc1614, ngc2369, ngc3256, ngc4194, iras12112, iras14348, iras17208, iras22491, Arp220, ESO0320, Zw049=loadtxt('Rieke_individuals_templates.dat', unpack=True)
temp1 = Zw049*1000.
#x_temp, ic342, mrk266, ngc520, ngc1097, ngc1222, ngc1365, ngc2623, ngc3310, ngc4088, ngc4676, ngc4818, ngc7252, ngc7714 = loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/starburst_template/individual_temp_Brandl_etal.dat', unpack=True)
#temp1 = ngc4088*1000.
#z_ic342 = 0.000103
#z_mrk266 = 0.027863
#z_ngc520 = 0.007609
#z_ngc1097 = 0.004240
#z_ngc1222 = 0.008079
#z_ngc1365 = 0.005457
#z_ngc2623 = 0.018509
#z_ngc3310 = 0.003312
#z_ngc4088 = 0.002524 
#z_ngc4676 = 0.022049
#z_ngc4818 = 0.003552
#z_ngc7252 = 0.015984
#z_ngc7714 = 0.009333
#x_temp = x_temp/(1.+z_ngc4088)
#print x_temp
#x_temp, m82 = loadtxt('/home/mariela/piratas-project/QSOs/SB-templates/m82_sws.dat', unpack = True)
#temp1 = m82*1000.
#zm82 = 0.000677
#x_temp = x_temp/(1.+zm82)
#x_temp, temp1, etemp1=loadtxt('templates/ESO_244-G_012_NED02.dat', unpack=True)
#x_temp=x_temp/(1.+0.022903)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_________________________________________________________________
#%%%%%%%%%%%%%%%%CCAM spectrum%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#x1,f1cc,ef1,flag=loadtxt('Spec_m1_stck_UGC5101obs113_114.dat',unpack=True)
x2,f2cc,ef,flag=loadtxt('Spec_m1_stck_UGC5101obs_215_216.dat',unpack=True)
#fccam=(f1cc+f2cc)/2.
#ef3=((0.2*f1cc)**(2.)+(0.2*f2cc)**(2.))**(0.5)
#xsc, ysc,eysc=loadtxt('ccam_spec_scaled.dat', unpack=True)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_____________________________________________
#%%%%%%%%%% Extinction Laws Library%%%%%%%%%%%%%%%
#--------------------------------------------
x_ext, A_gc = loadtxt('extinction_chiar_tielens2006.dat',unpack=True)
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
#plt.figure()
print 'Interpolating starbust'
f  = interp1d(x_temp, temp1, kind='cubic')
f1 = f(spit_x)
#plt.plot(x_temp, temp1,'o', spit_x, f1)
#plt.legend(['PAH template', 'cubic'], loc='best')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__________________________
#%%%%%%Interpoling Extintion Laws%%%%%
#---------------------------
#plt.figure()
print 'Interpolating extinction law'
exf  = interp1d(x_ext, tau_gc, kind='cubic')
f2 = exf(spit_x)
plt.plot(x_ext, tau_gc,'o', spit_x, f2)
plt.legend(['Extinction', 'cubic'], loc='best')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#_____________________________________
#%%%%%%%%%%%%Fitting %%%%%%%%%
#----------------------------------
f3 = open('parameters.dat','w+')
print 'Fitting...'
#p0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

sigma_obs_spec = 0.1*spit_y
#sigma_obs_spec = None
def func(xmodel, c2, c3, c4, c5): 
  return (((xmodel)**(c2)*np.exp(-c3*f2)) + (c4*f1*np.exp(-c5*f2)))
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
zf = func(spit_x, popt[0], popt[1], popt[2], popt[3])

a2 = popt[0]
a3 = popt[1]
a4 = popt[2]
a5 = popt[3]

y_cont = spit_x**(a2)*np.exp(-a3*f2) 

y_sb = a4*f1*np.exp(-a5*f2)

y_tot = y_cont + y_sb

# now plot the best fit curve and also +- 1 sigma curves
# (the square root of the diagonal covariance matrix  
# element is the uncertianty on the fit parameter.)
sigma = [sqrt(pcov[0,0]), \
         sqrt(pcov[1,1]), \
         sqrt(pcov[2,2]), \
         sqrt(pcov[3,3])]
#sig1 = sigma[0]
sig1 = 1.
sig2 = sigma[0]
sig3 = sigma[1]
sig4 = sigma[2]
sig5 = sigma[3]

tau_tot = (log(a4*f1) + a2*log(spit_x) - log(y_tot))/(a3+a5)
derv1 = abs(log(spit_x)/(a3+a5))
derv2 = abs((-log(a4*f1)-(a2*log(spit_x))+log(y_tot))/(a3+a5)**(2.))
derv3 = abs((1./((a3+a5)*a4)) + log(f1)/(a3+a5))
derv4 = abs((-log(a4*f1)-(a2*log(spit_x))+log(y_tot))/(a3+a5)**(2.))

delta_tau_tot = derv1*sig2 + derv2*sig3 + derv3*sig4 + derv4*sig5

#print 'Printing extinction - Spitzer'
#for i in range(0,175):
#print spit_x[i], tau_tot[i], delta_tau_tot[i]

a1 = (y_tot-(a4*f1*np.exp(-a5*f2)))*c**(a2)*np.exp(a3*f2)/spit_x**(a2)
da1 =c**(a2)*np.exp(a3*f2)*((log(c)*spit_x**(-a2))-(a2*spit_x**(-a2-1.)))*(y_tot-(a4*f1*np.exp(a5*f2)))
da2 =(y_tot-(a4*f1*np.exp(-a5*f2)))*c**(a2)*spit_x**(-a2)*f2*np.exp(a3*f2)
da3 =-f1*np.exp(-a5*f1)*c**(a2)*spit_x**(-a2)*np.exp(a3*f2)
da4 =a4*f1*f2*np.exp(-a5*f2)*c**(a2)*spit_x**(-a2)*np.exp(a3*f2)
da_tot = sqrt((da1*sig2)**(2.)+(da2*sig3)**(2.)+(da3*sig4)**(2.)+(da4*sig5)**(2.))

print 'Degrees of freddom=', dof
print 'a1=',a1[0],'+/-',da_tot[0]
print 'a2=',popt[0],'+/-',sigma[0]
print 'a3=',popt[1],'+/-',sigma[1] 
print 'a4=',popt[2],'+/-',sigma[2]
print 'a5=',popt[3],'+/-',sigma[3]
print 'rchi2=', rchi2, ',','chi2=',chi2

nspit_x=c/spit_x
y_cont2 = (a1*nspit_x**(-a2)*np.exp(-a3*f2))+y_sb

zf_sigma = sqrt((a2*np.exp(-a3*f2)*sig2)**(2) + (-f2*spit_x**(a2)*np.exp(-a3*f2)*sig3)**(2) + (f1*np.exp(-a5*f2)*sig4)**(2) + (-f2*a4*f1*np.exp(-a5*f2)*sig5)**(2))
y_cont_sigma =sqrt((a2*np.exp(-a3*f2)*sig2)**(2) + (-f2*spit_x**(a2)*np.exp(-a3*f2)*sig3)**(2))
y_sb_sigma =sqrt((f1*np.exp(-a5*f2)*sig4)**(2) + (-f2*a4*f1*np.exp(-a5*f2)*sig5)**(2))

print 'Printing continuum flux - Spitzer'
for i in range(0,175):
 print spit_x[i], y_cont[i], y_cont_sigma[i]

file400 = open('continuum_spitzer.dat','w+')
print >>file400, '#lambda(um)', 'f_nu(mJy)'
for i in range(0,175):
  print >>file400, spit_x[i]*(1.+ zred), y_cont[i], y_cont_sigma[i]
file400.close()


plt.figure()

plt.plot(lSIV,f_ind, 'k--',color='black')
plt.text(10.6,400,"[SIV]", color='black')
plt.plot(lPAH,f_ind, 'k--',color='black')
plt.text(11.4,400,"PAH", color='black')
plt.plot(lPAH_8,f_ind, 'k--',color='black')
plt.text(8.7,400,"PAH", color='black')
plt.plot(lPAH_7,f_ind, 'k--',color='black')
plt.text(7.8,400,"PAH", color='black')
plt.plot(lPAH_12, f_ind, 'k--', color='black')
plt.text(12.9,400,"PAH", color='black')

#plt.errorbar(spit_x, spit_y, 0.2*spit_y, color ='grey')
plt.fill_between(spit_x, spit_y, spit_y+(0.2*spit_y), color='grey')
plt.fill_between(spit_x, spit_y, spit_y-(0.2*spit_y), color='grey')
plt.plot(spit_x, spit_y, 'k-', label="UGC5101 Spitzer-composite")
plt.plot(spit_x, zf, 'r--', linewidth = 2, label="Continuum + PAH template")

plt.plot(spit_x, y_cont, 'g:',linewidth = 2, label="Continuum")
plt.plot(spit_x, y_sb, 'b-.', linewidth = 2, label="PAH template")


plt.plot(spit_x,y_cont+y_sb, color='magenta')

#Scaled to nuclear spectrum
#plt.plot(spit_x, y_sb/2.3, 'b-', label="PAH template (Zw049)")

#plt.plot((x1/10000.)/(1.+zred),f1cc*1000., color='green',linewidth=2.0, label='CCAM-spec')
#plt.plot((x1/10000.)/(1.+zred),f2cc*1000., color='orange',linewidth=2.0, label='CCAM-spec')
#plt.plot(xsc/(1.+zred),ysc*1000., color='brown',linewidth=2.0, label='CCAM-spec scaled')

plt.fill_between(spit_x, zf-zf_sigma, zf, color='pink')
plt.fill_between(spit_x, zf+zf_sigma, zf, color='pink')


plt.legend(loc = (0.01, 0.8), numpoints=1)
xz=[5.,14.]
yz=[0.,0.]
#plt.plot(xz, yz, color='black')
plt.ylim(-17,550)
plt.xlabel(" Rest-frame Wavelength ($\mu$m)",fontsize='20')
plt.ylabel(r"$f_{\nu}$ (mJy)",fontsize='20')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%Pritting components %%%%%%%%%
#----------------------------------
#file100 = open('starburst_component_Zw049.dat','w+')
#print >>file100, '#lambda(um)', 'f_nu(mJy)'
#for i in range(0,175):
#  print >>file100, spit_x[i], y_sb[i]
#file100.close()  
  
#file200 = open('nuclear_component_using_Zw049.dat','w+')
#print >>file200, '#lambda(um)', 'f_nu(mJy)'
#for i in range(0,175):
#  print >>file200, spit_x[i], y_cont[i]  
#file200.close()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%Nuclear photometry%%%%%%%%%%%%%%%%%%%%
lf110w=1.10
lf160w=1.60
lf222w=2.22
lccam=8.7
lk=2.2
lLp=3.8
lN=8.75
lq=17.69

lf110w=lf110w/(1.+zred)
lf160w=lf160w/(1.+zred)
lf222w=lf222w/(1.+zred)
lccam=lccam/(1.+zred)
lk=lk/(1.+zred)
lLp=lLp/(1.+zred)
# Ya estan en rest-frame
lN=lN
lq=lq

f110w_SN=1.18
f160w_SN=3.30
f222w_SN=10.20
fccam=36.
f160w_mm=1.27
f222w_mm=5.38
fk=10.223
fLp=38.830
fN=229.4
fq=157

efccam= 3.
ef160w_mm =0.2
ef222w_mm =0.3
efq=9.
efLp=7.14


#%%%%%%%%%%Ploting%%%%%%%%%%%%%%%%%%%%
#plt.errorbar(lf110w,f110w_SN,0.2*f110w_SN, fmt='o',color='magenta', label='HST (1.1") Scoville+02')
#plt.errorbar(lk, fk, 0.2*fk, fmt='o', color='pink', label='Imanishi+14')
#plt.errorbar(lf160w,f160w_mm,0.2*f160w_mm, fmt='o',color='green', label='HST (PSF) present work')
#plt.errorbar(lccam,fccam,0.2*fccam, fmt='o',color='red', label='CCAM (PSF) present work')
#plt.errorbar(lq,fq,0.2*fq,fmt='o', color='brown', label='Asmus+14')
#plt.errorbar(lf110w,f110w_SN,0.2*f110w_SN, fmt='o',color='magenta')
#plt.errorbar(lf222w,f222w_SN,0.2*f110w_SN, fmt='o',color='magenta')
#plt.errorbar(lf160w,f160w_SN,0.2*f110w_SN, fmt='o',color='magenta')
#plt.errorbar(lf160w,f160w_mm,0.2*f160w_mm, fmt='o',color='green')
#plt.errorbar(lf222w,f222w_mm,0.2*f222w_mm, fmt='o',color='green')
#plt.errorbar(lccam,fccam,0.2*fccam, fmt='o',color='red')
#plt.errorbar(lLp,fLp,0.2*fLp, fmt='o', color='pink')

#plt.errorbar((x2/10000.)/(1.+zred),f2cc*1000., ef*1000.,color='Moccasin')
#plt.plot((x2/10000.)/(1.+zred),f2cc*1000., color='orange',linewidth=2.0, label='CCAM-spec')
#plt.plot(spit_x, y_cont, 'k--', label="Continuum")
plt.xlim(5.2,13.5)
plt.legend(loc='best', numpoints=1)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Otra figura%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plt.figure('test')
#sca_fac=2.5
#sca_fac=3.
sca_fac=7.6

#plt.plot(spit_x, y_sb/3., label='SB template', color='red')
x_ccam = (x2/10000.)/(1.+zred)
y_ccam = f2cc*1000./2.4
ey_ccam=ef*1000./2.4

y_intr_ccam  = interp1d(x_ccam, y_ccam, kind='cubic')
y_ccam_intr = y_intr_ccam(spit_x[73:159])


nx_ccam = spit_x[73:159]
ny_ccam = y_sb[73:159]/sca_fac
ny_ccam_sigma = y_sb_sigma[73:159]/sca_fac

y_ccam_clean_1 = (y_ccam_intr - ny_ccam)
y_ccam_clean_1_sigma=y_ccam_intr-ny_ccam_sigma

#plt.plot(nx_ccam, y_agn_ccam, color='purple')

plt.errorbar(x_ccam,y_ccam, ey_ccam,color='Moccasin')
plt.plot(x_ccam,y_ccam, color='orange')


plt.plot(nx_ccam, y_ccam_intr, color='brown', label='Interpolated CCAM-spec')
#plt.plot(spit_x, y_cont/3.2, label='Continuum', color='black')

#plt.plot(nx_ccam, y_ccam_clean_1, 'g--',label='sin valor abs')



#plt.errorbar(lccam,fccam,0.2*fccam, fmt='o',color='red', label='CCAM (PSF) present work')
#plt.plot(xsc/(1.+zred),ysc*1000., color='brown',linewidth=2.0, label='CCAM-spec scaled')
print 'longitud del espectro nuclear',len(nx_ccam)
file300 = open('nuclear_nuclear_ccam_2.dat','w+')
print >>file300, '#lambda(um)', 'f_nu(mJy)'
for i in range(0,86):
  print >>file300, nx_ccam[i]*(1.+ zred), y_ccam_clean_1[i], 0.3*y_ccam_clean_1[i] 
file300.close()

#__________________________
#%%%%%%Interpoling Extintion Law for better Rieke template%%%%%
#---------------------------
#print x_temp[115]

print 'Interpolating extinction law to NIR'
exf2  = interp1d(x_ext, tau_gc, kind='cubic')
x_rieke=range(359)
x_rieke=x_temp[115:]
f2_rieke = exf2(x_rieke)
#plt.plot(x_ext, tau_gc,'o', x_rieke, f2_rieke)
#plt.legend(['Extinction law for better Rieke template', 'cubic'], loc='best')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%Making Rieke's component%%%%%%%%%%%%%
print  'Making Riekes component'
f1_rieke=range(359)
f1_rieke=temp1[115:]
y_sb_rieke = a4*f1_rieke*np.exp(-a5*f2_rieke)
y_sb_sigma_r =sqrt((f1_rieke*np.exp(-a5*f2_rieke)*sig4)**(2) + (-f1_rieke*a4*f2_rieke*np.exp(-a5*f2_rieke)*sig5)**(2))


file1 = open('starburst_component_Zw049_to_NIR.dat','w+')
print >>file1, '#lambda(um)', 'f_nu(mJy)', 'fnu_plus', 'fnu_minus'
for i in range(0,358):
#  print >>file1, x_rieke[i], y_sb_rieke[i]/3., y_sb_plus_r[i]/3., y_sb_minus_r[i]/3.
  print >>file1, x_rieke[i], y_sb_rieke[i]/sca_fac, y_sb_sigma_r[i]/sca_fac

file1.close()


#__________________________________________________
#%%%%%%%%%%%%%%%%%%Printing templates%%%%%%%%%%%%%
print 'plotting...'
x_comp1, ycomp1, eycomp1_sigma = loadtxt('starburst_component_Zw049_to_NIR.dat', unpack=True)

plt.plot(x_comp1, ycomp1, color='blue', linewidth =2, label='Template SB scaled')
plt.fill_between(x_comp1, ycomp1, ycomp1 - eycomp1_sigma, color='LightBlue')
plt.fill_between(x_comp1, ycomp1, ycomp1 + eycomp1_sigma, color='LightBlue')


plt.plot(nx_ccam, y_ccam_clean_1, color='black', linewidth = 2,label='CCAM-spec minus template SB scaled')
plt.fill_between(nx_ccam, y_ccam_clean_1, y_ccam_clean_1_sigma - y_ccam_clean_1, color='grey')
plt.fill_between(nx_ccam, y_ccam_clean_1, y_ccam_clean_1_sigma + y_ccam_clean_1, color='grey')


plt.errorbar(lf110w,f110w_SN,0.2*f110w_SN, fmt='o',color='magenta', label='HST (1.1") Scoville+02')
plt.errorbar(lk, fk, 0.2*fk, fmt='o', color='pink', label='Imanishi+14')
plt.errorbar(lf160w,f160w_mm, ef160w_mm, fmt='o',color='green', label='HST (PSF) present work')
plt.errorbar(lccam,fccam, efccam, fmt='o',color='red', label='CCAM (PSF) present work')
plt.errorbar(lq,fq, efq,fmt='o', color='brown', label='Asmus+14')
plt.errorbar(lf110w,f110w_SN,0.2*f110w_SN, fmt='o',color='magenta')
plt.errorbar(lf222w,f222w_SN,0.2*f110w_SN, fmt='o',color='magenta')
plt.errorbar(lf160w,f160w_SN,0.2*f110w_SN, fmt='o',color='magenta')

plt.errorbar(lf222w, f222w_mm, ef222w_mm, fmt='o',color='green')


plt.errorbar(lLp,fLp,0.2*fLp, fmt='o', color='black')


plt.legend(loc='best', numpoints=1)
plt.xlabel(" Rest-frame Wavelength ($\mu$m)",fontsize='20')
plt.ylabel(r"$f_{\nu}$ (mJy)",fontsize='20')
plt.ylim(-10,170)
plt.xlim(0.9,18)


plt.show()