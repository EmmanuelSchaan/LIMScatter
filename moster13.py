
# I believe that Minf is M_{200crit} in the Moster+13 notation

# According to Moster+13, sec 3.3.2,
# The stellar mass in the central galaxies is obtained 
# when inputing M_{200crit} for the halo mass.
# The stellar mass of the satellite galaxies is obtained
# when inputing the infall mass and infall redshift for the subhalo

def moster(Minf,z):
    """
    moster(Minf,z):
    Returns the stellar mass (M*) given Minf and z from Table 1 and
    Eq. (2,11-14) of Moster++13 [1205.5807].
    This version works in terms of Msun units, not Msun/h.
    To get "true" stellar mass, add 0.15 dex of lognormal scatter.
    To get "observed" stellar mass, add between 0.1-0.45 dex extra scatter.
    """
    zzp1  = z/(1+z)
    M1    = 10.0**(11.590+1.195*zzp1)
    mM    = 0.0351 - 0.0247*zzp1
    beta  = 1.376  - 0.826*zzp1
    gamma = 0.608  + 0.329*zzp1
    Mstar = 2*mM / ( (Minf/M1)**(-beta) + (Minf/M1)**gamma )
    Mstar*= Minf
    return(Mstar)
    #

