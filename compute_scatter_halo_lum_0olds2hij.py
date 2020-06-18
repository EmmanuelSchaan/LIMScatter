from headers import *

import universe
reload(universe)
from universe import *

import mass_function
reload(mass_function)
from mass_function import *


nProc = 4


#####################################################################################

pathFig = "./figures/halo_luminosities/"
if not os.path.exists(pathFig):
   os.mkdir(pathFig)

#####################################################################################

u = UnivPlanck15()

#massFunc = MassFuncPS(u, nProc=nProc, save=True)
massFunc = MassFuncST(u, nProc=nProc, save=True)
#massFunc = MassFuncTinker(u, nProc=nProc, save=True)







#####################################################################################
#####################################################################################
# SPHEREx lines:
#Ha, Hb, Lya, OIII 5007A


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



#''' SFR [Msun/yr] as a function of halo mass [Msun] and redshift
#from Fonseca+16 (1607.05288v2), Eq 11.
#I inferred the units from Eq 9 and 11.
#'''
## Table I 
#Z = np.array([0., 1., 2., 3., 4., 5.])
#M0 = np.array([3.0e-10, 1.7e-9, 4.0e-9, 1.1e-8, 6.6e-8, 7.0e-7])  # [Msun/yr]
#Mb = np.array([6.0e10, 9.0e10, 7.0e10, 5.0e10, 5.0e10, 6.0e10])   # [Msun]
#Mc = np.array([1.0e12, 2.0e12, 2.0e12, 3.0e12, 2.0e12, 2.0e12])   # [Msun]
#a = np.array([3.15, 2.9, 3.1, 3.1, 2.9, 2.5])
#b = np.array([-1.7, -1.4, -2.0, -2.1, -2.0, -1.6])
#c = np.array([-1.7, -2.1, -1.5, -1.5, -1.0, -1.0])
## below Eq 11
#Ma = 1.e8   # [Msun]
## Eq 11
#Sfr = M0 * (m/Ma)**a * (1.+m/Mb)**b * (1.+m/Mc)**c
#
## interpolate
#sfr = interp1d(Z, Sfr, kind='linear', bounds_error=False, fill_value=0.)


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
# Halo mass and SFR scatter as a function of redshift


def nHalo(z, mMin=None):
   '''Mean number density of halos of all masses, at redshift z
   nHalo(z) = int dm dn/dm [(Mpc/h)^{-3}]
   '''
   if mMin is None:
      mMin = massFunc.mMin
   
   def f(lnM):
      m = np.exp(lnM)
      return m * massFunc.fmassfunc(m, 1./(1.+z))

   result = integrate.quad(f, np.log(mMin), np.log(massFunc.mMax), epsabs=0., epsrel=1.e-3)[0]
   return result


def haloAverage(z, f=lambda m: 1., mMin=None, norm=True):
   '''Compute int dm dn/dm f(m) / int dm dn/dm  [f unit]
   '''
   if mMin is None:
      mMin = massFunc.mMin

   def integrand(lnm):
      m = np.exp(lnm)
      return m * massFunc.fmassfunc(m, 1./(1.+z)) * f(m)
   result = integrate.quad(integrand, np.log(mMin), np.log(massFunc.mMax), epsabs=0., epsrel=1.e-3)[0]
   if norm:
      result /= nHalo(z)
   return result

def meanHaloMass(z, mMin=None):
   '''Mean mass of all halos
   meanHaloMass(z) = <m> = int dm dn/dm m  /  nHalo(z)   [Msun/h]
   '''
   return haloAverage(z, lambda m: m, mMin=mMin)

def meanHaloMassSquared(z, mMin=mMin):
   '''Mean mass of all halos
   meanHaloMassSquared(z) = <m^2> = int dm dn/dm m^2  /  nHalo(z)   [Msun/h]
   '''
   return haloAverage(z, lambda m: m**2, mMin=mMin)

def relativeVarHaloMass(z, mMin=None):
   '''Relative mass variance across halos
   varHaloMass(z) = <m^2> / <m>^2 - 1. [dimless]
   '''
   result = meanHaloMassSquared(z, mMin) / meanHaloMass(z, mMin)**2 - 1.
   return result

def meanSfr(z, alpha=1, mMin=None):
   '''Computes <SFR^alpha>,
   where the average is over the halo mass function.
   '''
   return haloAverage(z, lambda m: sfr(m * u.bg.h, z)**alpha, mMin=mMin)


def meanSfrSquared(z, mMin):
   return haloAverage(z, lambda m: sfr(m * u.bg.h, z)**2, mMin=mMin)

def relativeVarHaloSfr(z, mMin=None):
   '''<SFR^2> / <SFR>^2, dimensionless.
   '''
   result = meanSfr(z, alpha=2, mMin=mMin) / meanSfr(z, mMin=mMin)**2 - 1.
   return result



def nHaloEff(z, mMin=None):
   '''Effective mean number density of halos for the 1-halo term,
   i.e. P1h = I_1I_2 / nHaloEff
   '''
   if mMin is None:
      mMin = massFunc.mMin

  result = meanSfr(z, alpha=2, mMin=mMin)
  result /= meanSfr(z, alpha=1, mMin=mMin)**2
  result /= nHalo(z, mMin=mMin)
  return result






#####################################################################################
# Sparsity of halos


# Evaluate all
z = np.linspace(0., 5., 6)
NHalo = np.array(map(nHalo, z))
NHaloEff = np.array(map(nHaloEff, z))
mMean = np.array(map(meanHaloMass, z))
m2Mean = np.array(map(meanHaloMassSquared, z))
sigma2HaloMass = np.array(map(relativeVarHaloMass, z))

sfrMean = np.array(map(meanSfr, z))
sfr2Mean = np.array(map(meanSfrSquared, z))
sigma2HaloSfr = np.array(map(relativeVarHaloSfr, z))








# Effective mass cutoff from SFR,
# which determines which halo masses contribute
# to the effective number density of halos
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(z, sigma2HaloMass, 'k')
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$<\text{SFR}>^2 / <\text{SFR}^2>$')
ax.set_yscale('log', nonposy='clip')
#
fig.savefig(pathFig + "mass_cutoff_sfr.pdf", bbox_inches='tight')
plt.show()
fig.clf()





# Effective mean number density of halos
# which contribute to the 1-halo term
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(z, NHaloEff, 'k', label=r'$\bar{n}_{h\text{ eff}}$')
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\bar{n}_{h\text{ eff}} (z)$ [(Mpc/h)$^{-3}$]')
ax.set_yscale('log', nonposy='clip')
#
fig.savefig(pathFig + "nheff.pdf", bbox_inches='tight')
plt.show()
fig.clf()







# mean halo number density
# as a function of the minimum mass cut
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# mean number density of halos weighted by SFR
f = lambda z: nHaloEff(z, mMin=None)
y = array(map(f, z))
ax.plot(z, y, 'k', label=r'$\bar{n}_h$')

#
# Mean number density of halos for each mMin cutoff
nMMin = 5
MMin = np.logspace(np.log10(1.e8), np.log10(1.e12), nMMin, 10.)
for iMMin in range(nMMin):
   mMin = MMin[iMMin]
   f = lambda z: nHalo(z, mMin)
   y = np.array(map(f, z))
   ax.plot(z, y, c=plt.cm.autumn(iMMin/(nMMin-1.)), label=floatExpForm(mMin))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\bar{n}_h (z)$ [(Mpc/h)$^{-3}$]')
ax.set_yscale('log', nonposy='clip')
#
fig.savefig(pathFig + "nhalo.pdf", bbox_inches='tight')
plt.show()
fig.clf()





# SPHEREx voxel size
z = np.linspace(0., 5., 6)
# the spectral resolution power is R=40 for the lower z, and 150 for the high z
R = 40.
# hence the redshift size of the voxel
dz = (1. + z) / R
# and the comoving depth of the voxel
dChi = dz * 3.e5 / u.hubble(z)   # [Mpc/h]
# angular pixel size: 6.2 arcsec
thetaPix = 6.2 * np.pi/(180.*3600.)
# hence the voxel comoving volume
vVoxSpherex = (u.bg.comoving_distance(z) * thetaPix)**2 * dChi  # [(Mpc/h)^3]



# SPHEREx: mean halo number per voxel
# as a function of the minimum mass cut
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
#
nMMin = 5
MMin = np.logspace(np.log10(1.e8), np.log10(1.e12), nMMin, 10.)
for iMMin in range(nMMin):
   mMin = MMin[iMMin]
   f = lambda z: nHalo(z, mMin)
   y = np.array(map(f, z))
   ax.plot(z, y * vVoxSpherex, c=plt.cm.autumn(iMMin/(nMMin-1.)), label=floatExpForm(mMin))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\bar{n}_h (z)$ [(Mpc/h)$^{-3}$]')
ax.set_yscale('log', nonposy='clip')
#
fig.savefig(pathFig + "nhalo_voxel_spherex.pdf", bbox_inches='tight')
plt.show()
fig.clf()


















'''
# mean halo mass
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(z, mMean)
#
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\langle m_h \rangle (z)$ [M$_\odot$/h]')
ax.set_yscale('log', nonposy='clip')
#
fig.savefig(pathFig + "m_mean.pdf", bbox_inches='tight')
plt.show()
fig.clf()



# mean squared halo mass
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(z, m2Mean)
#
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\langle m_h^2 \rangle (z)$ [(M$_\odot$/h)$^2$]')
ax.set_yscale('log', nonposy='clip')
#
fig.savefig(pathFig + "m2_mean.pdf", bbox_inches='tight')
plt.show()
fig.clf()



# relative variance of halo mass
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(z, sigma2HaloMass)
#
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\sigma^2_\text{halo mass} \equiv \langle m_h^2 \rangle / \langle m_h \rangle^2 - 1$')
ax.set_yscale('log', nonposy='clip')
#
fig.savefig(pathFig + "sigma2_halo_mass.pdf", bbox_inches='tight')
plt.show()
fig.clf()


# mean SFR
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(z, sfrMean)
#
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\langle \text{SFR} \rangle (z)$ [M$_\odot$ / yr]')
#
fig.savefig(pathFig + "sfr_mean.pdf", bbox_inches='tight')
plt.show()
fig.clf()

# mean squared SFR
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(z, sfr2Mean)
#
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\langle \text{SFR}^2 \rangle (z)$ [(M$_\odot$ / yr)$^2$]')
#
fig.savefig(pathFig + "sfr2_mean.pdf", bbox_inches='tight')
plt.show()
fig.clf()

# Relative variance of halo SFR
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(z, sigma2HaloSfr)
#
ax.set_ylim((0., 1.1 * np.max(sigma2HaloSfr)))
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\sigma^2_\text{halo SFR} \equiv \langle \text{SFR}^2 \rangle / \langle \text{SFR} \rangle^2 - 1$')
#
fig.savefig(pathFig + "sigma2_halo_sfr.pdf", bbox_inches='tight')
plt.show()
fig.clf()
'''


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





