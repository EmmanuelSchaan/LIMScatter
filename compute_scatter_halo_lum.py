from headers import *

import universe
reload(universe)
from universe import *

import mass_function
reload(mass_function)
from mass_function import *

import profile
reload(profile)
from profile import *


nProc = 4


#####################################################################################

pathFig = "./figures/halo_luminosities/"
if not os.path.exists(pathFig):
   os.mkdir(pathFig)

#####################################################################################

u = UnivPlanck15()

#massFunc = MassFuncPS(u, nProc=nProc, save=False)
massFunc = MassFuncST(u, nProc=nProc, save=False)
#massFunc = MassFuncTinker(u, nProc=nProc, save=False)







#####################################################################################
#####################################################################################
# SFR, SFRD and effective halo mass cutoff


#####################################################################################
# SFR

def sfr(m, z):
   ''' SFR [Msun/yr] as a function of halo mass [Msun] and redshift
   from Fonseca+16 (1607.05288v2), Eq 11.
   I inferred the units from Eq 9 and 11.
   '''

   # bounds for the fitting function
   if z<0. or z>5.:
      return 0.

   # Table I 
   Z = np.array([0., 1., 2., 3., 4., 5.])
   M0 = np.array([3.0e-10, 1.7e-9, 4.0e-9, 1.1e-8, 6.6e-8, 7.0e-7])  # [Msun/yr]
   Mb = np.array([6.0e10, 9.0e10, 7.0e10, 5.0e10, 5.0e10, 6.0e10])   # [Msun]
   Mc = np.array([1.0e12, 2.0e12, 2.0e12, 3.0e12, 2.0e12, 2.0e12])   # [Msun]
   a = np.array([3.15, 2.9, 3.1, 3.1, 2.9, 2.5])
   b = np.array([-1.7, -1.4, -2.0, -2.1, -2.0, -1.6])
   c = np.array([-1.7, -2.1, -1.5, -1.5, -1.0, -1.0])
   
   # interpolate values
   fM0 = interp1d(Z, M0, kind='linear', bounds_error=False, fill_value=0.) # [Msun/yr]
   fMb = interp1d(Z, Mb, kind='linear', bounds_error=False, fill_value=0.) # [Msun]
   fMc = interp1d(Z, Mc, kind='linear', bounds_error=False, fill_value=0.) # [Msun]
   fa = interp1d(Z, a, kind='linear', bounds_error=False, fill_value=0.)
   fb = interp1d(Z, b, kind='linear', bounds_error=False, fill_value=0.)
   fc = interp1d(Z, c, kind='linear', bounds_error=False, fill_value=0.)
   
   # below Eq 11
   Ma = 1.e8   # [Msun]

   # Eq 11
   result = fM0(z) * (m/Ma)**fa(z) * (1.+m/fMb(z))**fb(z) * (1.+m/fMc(z))**fc(z)
   return result




#####################################################################################
# SFRD

def sfrd(z):
   '''SFR density:
   \int dm dn/dm SFR(m) [Msun / yr / Mpc^3]
   '''
   def integrand(lnm):
      '''the mass in integral is in [Msun/h]
      '''
      m = np.exp(lnm)
      return m * massFunc.fmassfunc(m, 1./(1.+z)) * sfr(m / u.bg.h, z)
   # [Msun /yr / (Mpc/h)^3]
   result = integrate.quad(integrand, np.log(massFunc.mMin), np.log(massFunc.mMax), epsabs=0., epsrel=1.e-3)[0]
   # [Msun /yr / Mpc^3]
   result *= u.bg.h**3 
   return result

   
def sfrdBehroozi13(z):
  return 0.180 / (10.**(-0.997*(z-1.243)) + 10.**(0.241*(z-1.243))) 

'''
# SFRD
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# SMill fitting function
z = np.linspace(0., 5., 6)
Sfrd = np.array(map(sfrd, z))
ax.plot(z, np.log10(Sfrd), 'b--', label=r'Fonseca+16')
#ax.plot(z, Sfrd)
#
# Behroozi+13 parametrization
z = np.linspace(0., 5., 101)
ax.plot(z, np.log10(sfrdBehroozi13(z)), 'r-', label=r'Behroozi+13')
#
ax.legend(loc='lower center')
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'log$_10$(SFRD$(z)$/[M$_\odot$ / yr / Mpc$^3$])')
#
fig.savefig(pathFig + "sfrd.pdf", bbox_inches='tight')
plt.show()
fig.clf()
'''


#####################################################################################
#####################################################################################
# Effective mean halo number density, and number per voxel


def haloMassIntegral(z, f=lambda m: 1., mMin=None, mMax=None):
   '''Compute int dm dn/dm f(m)  [(f unit) / (Mpc/h)^3]
   m [Msun/h]
   '''
   if mMin is None:
      mMin = massFunc.mMin
   if mMax is None:
      mMax = massFunc.mMax

   def integrand(lnm):
      m = np.exp(lnm)
      return m * massFunc.fmassfunc(m, 1./(1.+z)) * f(m)
   result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0., epsrel=1.e-3)[0]
   return result

def nHaloEff(z, mMin=None, a1=1., a2=1.):
   '''
   computes \bar{n}^{h eff}_{12} [1/(Mpc/h)^3]
   for lines 1 and 2, such that 
   L1 propto SFR^a1, L2 propto SFR^a2.
   m [Msun/h]
   '''
   result = haloMassIntegral(z, lambda m: sfr(m * u.bg.h, z)**a1, mMin=mMin)
   result *= haloMassIntegral(z, lambda m: sfr(m * u.bg.h, z)**a2, mMin=mMin)
   result /= haloMassIntegral(z, lambda m: sfr(m * u.bg.h, z)**(a1+a2), mMin=mMin)
   return result



"""
Z = np.array([0.001, 1., 2., 3., 4., 5.])

# SPHEREx voxel size
# the spectral resolution power is R=40 for the lower z, and 150 for the high z
R = 40.
# hence the redshift size of the voxel
dz = (1. + Z) / R
# and the comoving depth of the voxel
dChi = dz * 3.e5 / u.hubble(Z)   # [Mpc/h]
# angular pixel size: 6.2 arcsec
thetaPix = 6.2 * np.pi/(180.*3600.)
# hence the voxel comoving volume
vVoxSpherex = (u.bg.comoving_distance(Z) * thetaPix)**2 * dChi  # [(Mpc/h)^3]



# range of reasonable scalings for L = SFR^alpha
Alpha = [0.8, 1., 1.1]



fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for alpha in Alpha:
   f = lambda z: nHaloEff(z, a1=alpha, a2=alpha)
   NHaloEff = np.array(map(f, Z))
   ax.plot(Z, NHaloEff, label=r'$\alpha_i=\alpha_j=$'+str(alpha))
#
ax.legend(loc=2, labelspacing=0.2)
#ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\bar{n}^\text{h eff}_{ij}$ [(Mpc/h)$^3$]')
#
fig.savefig(pathFig+"nheff.pdf", bbox_inches='tight')
plt.show()
fig.clf()



fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for alpha in Alpha:
   f = lambda z: nHaloEff(z, a1=alpha, a2=alpha)
   NHaloEff = np.array(map(f, Z))
   #ax.plot(Z, vVoxSpherex)
   #ax.plot(Z, NHaloEff)
   ax.plot(Z, NHaloEff * vVoxSpherex, label=r'$\alpha_i=\alpha_j=$'+str(alpha))
#
ax.legend(loc=2, labelspacing=0.2)
#ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\bar{N}^\text{h eff}_{ij}$ per SPHEREx voxel')
#
fig.savefig(pathFig+"halo_sparsity_spherex.pdf", bbox_inches='tight')
plt.show()
fig.clf()
"""


#####################################################################################
# Mmax dependence of the low k 1-halo term

profNFW = ProfNFW(u)

def fdP1hdlnm(m, z, mMin=None):
   result =  m * massFunc.fmassfunc(m, 1./(1.+z))
   result *= sfr(m * u.bg.h, z)**2
   return result

'''
M = np.logspace(np.log10(1.e10), np.log10(2.e14), 101, 10.) # masses in h^-1 solarM

Z = np.array([0.001, 1., 2., 3., 4., 5.])


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for iZ in range(len(Z)):
   z = Z[iZ]
   f = lambda m: fdP1hdlnm(m, z)
   y = np.array(map(f, M))
   ax.plot(M, y, c=plt.cm.autumn_r(iZ/(len(Z)-1.)), label=r'$z=$'+str(int(z)))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.2)
ax.set_xscale('log', nonposx='clip')
ax.set_xlabel(r'Halo mass $m$ [$M_\odot/h$]')
ax.set_ylabel(r'$\frac{dP^\text{1-halo}(k=0, \mu=0, z)}{d\text{ln} m} = m\ n(m)\ \text{SFR}(m)^2 $  [arbitrary unit]', fontsize=14)
#
fig.savefig(pathFig+"dp1hdm.pdf", bbox_inches='tight')
plt.show()
'''

#####################################################################################
# Mmax and k dependences of the 1-halo term



def fP1h(k, z, mMin=None, mMax=None):
   '''Assumes mu=0, ie no RSD.
   '''
   f = lambda m: sfr(m * u.bg.h, z)**2 * profNFW.nfw(k, m, z)**2
   result = haloMassIntegral(z, f, mMin=mMin, mMax=mMax)
   return result

'''
# Contributions of each halo mass to the 1-halo term
MMax = np.logspace(np.log10(1.e12), np.log10(5.e13), 5, 10.) # masses in h^-1 solarM
K = np.logspace(np.log10(1.e-3), np.log10(10.), 51, 10.)
Z = [1., 2.]

for z in Z:
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   for iMMax in range(len(MMax))[::-1]:
      mMax = MMax[iMMax]
      f = lambda k: fP1h(k, z, mMax=mMax)
      P1h = np.array(map(f, K))
      print P1h
      ax.plot(K, P1h, c=plt.cm.winter_r(iMMax/(len(MMax)-1.)), label=r'$M_\text{max}=$'+intExpForm(mMax, round=1)+r'$M_\odot$/h')
   #
   ax.legend(loc=3, fontsize='x-small', labelspacing=0.2)
   ax.set_xscale('log', nonposx='clip')
   ax.set_yscale('log', nonposy='clip')
   #ax.set_xlim((1.e10, 4.e14))
   ax.set_xlabel(r'$k$ [h/Mpc]')
   ax.set_ylabel(r'$P^\text{1h}_{ij}(k, \mu=0, z='+str(int(z))+r')$ [arbitrary unit]')
   #
   fig.savefig(pathFig+"p1h_mmax_z"+str(int(z))+".pdf", bbox_inches='tight')
   plt.show()
   fig.clf()
'''


#####################################################################################
# z dependence of the 1-halo term

'''
K = np.logspace(np.log10(1.e-3), np.log10(10.), 51, 10.)
Z = np.array([0.001, 1., 2., 3., 4., 5.])


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for iZ in range(len(Z)):
   z = Z[iZ]
   f = lambda k: fP1h(k, z)
   y = np.array(map(f, K))
   ax.plot(K, y, c=plt.cm.autumn_r(iZ/(len(Z)-1.)), label=r'$z=$'+str(int(z)))
#
ax.legend(loc=3, fontsize='x-small', labelspacing=0.2)
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$k$ [h/Mpc]')
ax.set_ylabel(r'$P^\text{1h}_{ij}(k, \mu=0, z)$ [arbitrary unit]')
#
fig.savefig(pathFig+"p1h.pdf", bbox_inches='tight')
plt.show()
'''


#####################################################################################
#####################################################################################







#####################################################################################
#####################################################################################
# SPHEREx lines:
#Ha, Hb, Lya, OIII 5007A



#####################################################################################
# Lyman-alpha 121.6 nm


def kLyA(z):
   '''Kennicut-Schmidt constant for Lyman-alpha luminosity
   from Fonseca+16 (1607.05288v2).
   '''

   # Eq 16: 
   # for case B recombination and optically thick ISM
   rLyARec = 1.1e42  # [erg/s]
   #rLyAExc = 4.0e41  # [erg/s], ignored because subdominant

   # Eq 17
   fUVEsc = 0.2
   # dust extinction of 1 mag
   #EUVMin = 0.8
   #EUVMax = 1.2
   EUVTypical = 1.
   fUVDust = 10.**(-EUVTypical/2.5)

   # line escape fraction, Eq 18
   C = 1.67e-3 # [dimless]
   xi = 2.57   # [dimless]
   fLyAEsc = lambda z: C * (1.+z)**xi

   # Kennicut-Schmidt constant, Eq 15
   kLyA = lambda z: (fUVDust - fUVEsc) * fLyAEsc(z) * rLyARec
   return kLyA(z)


def luminosityLyA(m, z):
   '''Lyman-alpha luminosity
   from Fonseca+16 (1607.05288v2).
   '''
   gLyA = 1. 
   result = kLyA(z) * sfr(m, z)**gLyA
   return result


#####################################################################################
# H-alpha 656.281 nm

def kHA(z):
   '''Kennicut-Schmidt constant for H-alpha luminosity
   from Fonseca+16 (1607.05288v2).
   '''

   # Eq 19: 
   # for case B recombination and optically thick ISM
   rHA = 1.3e41  # [erg/s]

   # Eq 17
   fUVEsc = 0.2
   # dust extinction of 1 mag
   #EUVMin = 0.8
   #EUVMax = 1.2
   EUVTypical = 1.
   fUVDust = 10.**(-EUVTypical/2.5)

   # line escape fraction, Eq 17
   EHATypical = 1.
   fHAEsc = 10.**(-EHATypical/2.5)

   # Kennicut-Schmidt constant, Eq 15
   kHA = lambda z: (fUVDust - fUVEsc) * fHAEsc * rHA
   return kHA(z)


def luminosityHA(m, z):
   '''H-alpha luminosity
   from Fonseca+16 (1607.05288v2).
   '''
   gHA = 1. 
   result = kHA(z) * sfr(m, z)**gHA
   return result


#####################################################################################
# H-beta 486.1 nm

def kHB(z):
   '''Kennicut-Schmidt constant for H-beta luminosity
   from Fonseca+16 (1607.05288v2).
   '''

   # Eq 20: 
   # for case B recombination and optically thick ISM
   rHB = 4.45e40  # [erg/s]

   # Eq 17
   fUVEsc = 0.2
   # dust extinction of 1 mag
   #EUVMin = 0.8
   #EUVMax = 1.2
   EUVTypical = 1.
   fUVEsc = 10.**(-EUVTypical/2.5)

   # line escape fraction, Eq 18
   EHBTypical = 1.38
   fHBEsc = 10.**(-EHBTypical/2.5)

   # Kennicut-Schmidt constant, Eq 15
   kHB = lambda z: (fUVDust - fUVEsc) * fHBEsc * rHB
   return kHA(z)


def luminosityHB(m, z):
   '''H-beta luminosity
   from Fonseca+16 (1607.05288v2).
   '''
   gHB = 1. 
   result = kHB(z) * sfr(m, z)**gHB
   return result


#####################################################################################
# Doublet [OIII] 500.7 nm and [OIII] 495.9 nm
# The luminosity below really represents the sum of the two lines in the doublet.


def kOIII5007(z):
   '''Kennicut-Schmidt constant for [OIII] doublet (500.7nm and 495.9 nm) luminosity
   from Fonseca+16 (1607.05288v2).
   '''

   # Eq 202: 
   # for case B recombination and optically thick ISM
   rHB = 1.3e41  # [erg/s]

   # Eq 17
   fUVEsc = 0.2
   # dust extinction of 1 mag
   #EUVMin = 0.8
   #EUVMax = 1.2
   EUVTypical = 1.
   fUVEsc = 10.**(-EUVTypical/2.5)

   # line escape fraction, Eq 18
   EHBTypical = 1.35
   fHBEsc = 10.**(-EHBTypical/2.5)

   # Kennicut-Schmidt constant, Eq 15
   kHB = lambda z: (fUVDust - fUVEsc) * fHBEsc * rHB
   return kHA(z)


def luminosityOIII5007(m, z):
   '''[OIII] doublet (500.7 nm and 494.9 nm) luminosity
   from Fonseca+16 (1607.05288v2).
   '''
   gOIII5007 = 1. 
   result = kOIII5007(z) * sfr(m, z)**gOIII5007
   return result


#####################################################################################
#####################################################################################


#####################################################################################
# Correlation of the halo line noise
# if two lines scale as SFR to the same power, then the halo line noises
# are 100% correlated.
# In the literature, slightly different scalings are given for the line luminosities
# as a function of SFR.
# Here, I give the correlation coefficients for several line combinations.

'''
# NII 122mu
# Fonseca+16, from Spinoglio+12
gNII = 1.01 # +/- 0.04

# NIII 58mu
# Fonseca+16, from Spinoglio+12
gNIII = 0.78 # +/- 0.10

# CII 158mu
# Finseca+16, from de Looze+11
gCII = 1.02

# CO transitions lower than 4-3
# Fonseca+16, from Sargent+14
gCOLow = 0.81
# CO transitions 4-3 and higher
# Finseca+16, from Liu+15
gCOHigh = 1.



#lines = ['Lya', 'Ha', 'Hb', 'OIII5007', 'NII122', 'NIII58', 'CII158', 'COLow', 'COHigh']
#gamma = [1., 1., 1., 1., 1.01, 0.78, 1.02, 0.81, 1.]
lines = ['NII122', 'NIII58', 'CII158', 'COLow', 'COHigh']
gamma = [1.01, 0.78, 1.02, 0.81, 1.]

z = np.linspace(0., 5., 6)
s2ijh = np.zeros((len(lines), len(lines), len(Z)))
rijh = np.zeros((len(lines), len(lines), len(Z)))



# compute the covariance matrix
# of the halo line luminosities
for iLine1 in range(len(lines)):
   g1 = gamma[iLine1]
   for iLine2 in range(len(lines)):
      g2 = gamma[iLine2]

      f = lambda z: meanSfr(z, alpha=g1+g2) / meanSfr(z, alpha=g1) / meanSfr(z, alpha=g2) - 1.
      s2ijh[iLine1, iLine2, :] = np.array(map(f, Z))


# compute correlation coefficient
# of the halo line luminosities
for iZ in range(len(Z)):
   s = np.sqrt(np.diag(s2ijh[:,:,iZ]))
   rijh[:,:,iZ] = s2ijh[:,:,iZ] / np.outer(s,s)


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for iLine1 in range(len(lines)):
   line1 = lines[iLine1]
   for iLine2 in range(iLine1,len(lines)):
      line2 = lines[iLine2]
      ax.plot(Z, rijh[iLine1, iLine2], label=line1+', '+line2)
#
ax.legend(loc=3, fontsize='x-small', labelspacing=0.1, handlelength=0.1, ncol=2)
ax.set_ylim((0.92, 1.001))
fig.savefig(pathFig+"rijh.pdf", bbox_inches='tight')

plt.show()
'''





