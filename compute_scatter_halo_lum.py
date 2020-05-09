from headers import *

import universe
reload(universe)
from universe import *

import mass_function
reload(mass_function)
from mass_function import *


nProc = 3


u = UnivPlanck15()

#massFunc = MassFuncPS(u, nProc=nProc, save=True)
#massFunc = MassFuncST(u, nProc=nProc, save=True)
massFunc = MassFuncTinker(u, nProc=nProc, save=False)


# SPHEREx lines:
#Ha, Hb, Lya, OIII 5007A


def sfr(m, z):
   ''' SFR [Msun/yr] as a function of halo mass [Msun] and redshift
   from Fonseca+16 (1607.05288v2), Eq 11.
   I inferred the units from Eq 9 and 11.
   '''
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


def luminosityHA(m, z):
   '''Ha luminosity
   '''
   pass


def luminosityLyA(m, z):
   '''Lyman-alpha luminosity
   '''

   # Eq 16: 
   # for case B recombination and optically thick ISM
   rLyARec = 1.1e42  # [erg/s]
   rLyaExc = 4.0e41  # [erg/s]

   # Eq 17
   fUVEsc = 0.2
   # dust extinction of 1 mag
   EUVMin = 0.8
   EUVMax = 1.2
   EUVTypical = 1.
   fUVDust = 10.**(-EUVTypical/2.5)

   # Lyman alpha escape fraction
   fLyaEsc = lambda z: C * (1.+z)**xi
