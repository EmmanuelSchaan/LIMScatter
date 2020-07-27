from headers import *
import pycolfits
from universe import *

#####################################################################################

u = UnivPlanck15()


#####################################################################################
plot = False


pathCat = "./output/EGG/catalog_EGG_LIM.fits"
pathFig = "./figures/EGG/"
pathOut = "./output/EGG_results/"

if not os.path.exists(pathFig):
   os.makedirs(pathFig)
if not os.path.exists(pathOut):
   os.makedirs(pathOut)


# Read the catalog in memory
cat = pycolfits.readfrom(pathCat, lower_case=True)
#cat['id'] # galaxy id
#cat['m'] # log10( stellar mass / msun )
#cat['z'] # redshift
#cat['sfr'], cat['sfruv'], cat['sfrir'] # various SFR?
#cat['lines']  # line names


nLines = len(cat['lines'])
nBins = 51


#####################################################################################

# Save line names
np.savetxt(pathOut + "line_names.txt", cat['lines'], delimiter=" ", newline = "\n", fmt="%s")

# save line wavelengths [mu]
np.savetxt(pathOut + "line_lambda_microns.txt", cat['line_lambda'])


#####################################################################################
# Explore basic properties
'''
# Redshifts
myHistogram(cat['z'], nBins=101, lim=None, S2Theory=[], path=pathFig+'z.pdf', plot=plot, nameLatex=r'$z$', semilogx=False, semilogy=True, doGauss=False)

# Stellar masses
myHistogram(cat['m'], nBins=101, lim=None, S2Theory=[], path=pathFig+'log10mstellar.pdf', plot=plot, nameLatex=r'$\text{log}_{10}\left(M_\star/M_\odot\right)$', semilogx=False, semilogy=True, doGauss=False)
'''

#####################################################################################

# Define redshift bins
nZE = 11
nZC = nZE-1
# z bin edges
zE = np.linspace(np.min(cat['z']), np.max(cat['z']), nZE)
dZ = zE[1:] - zE[:-1]
# get effective bin centers and bin counts
zC, zE, zBinIndices = stats.binned_statistic(cat['z'], cat['z'], statistic='mean', bins=zE)
zBinCounts, _, _ = stats.binned_statistic(cat['z'], cat['z'], statistic='count', bins=zE)

## Comoving volume per redshift bin
#skyArea = 1. * (np.pi/180.)**2   # [sr]
#chi = u.bg.comoving_distance(zC) # [Mpc/h]
#hubble = u.hubble(zC)   # [km/s/(Mpc/h)]
#comovVolume = chi**2 * skyArea * (3.e5/hubble) * dZ   #[(Mpc/h)^3]

## Comoving volume per redshift bin
#skyArea = 1. * (np.pi/180.)**2   # [sr]
#chi = u.bg.comoving_distance(zC) # [Mpc/h]
#dChi = u.bg.comoving_distance(zE[1:]) - u.bg.comoving_distance(zE[:-1])
##hubble = u.hubble(zC)   # [km/s/(Mpc/h)]
#comovVolume = chi**2 * skyArea * dChi  #(3.e5/hubble) * dZ   #[(Mpc/h)^3]

# Comoving volume per redshift bin
skyArea = 1. * (np.pi/180.)**2   # [sr]
def integrand(z):
   chi = u.bg.comoving_distance(z)
   result = skyArea * chi**2  # [(Mpc/h)^2]
   result *= 3.e5/u.hubble(z) # [(Mpc/h)]
   return result

comovVolume = np.zeros(nZC)
for iZ in range(nZC):
   comovVolume[iZ] = integrate.quad(integrand, zE[iZ], zE[iZ+1], epsabs=0., epsrel=1.e-3)[0]






'''
# Check my comoving volume
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.semilogy(zC, 4.*np.pi / skyArea * comovVolume.cumsum() / 1.e9)
#
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'Comoving volume of Universe [(Gpc/h)$^3$]')
plt.show()
'''


#####################################################################################
# galaxy number density

nGal = zBinCounts / comovVolume
# Poisson uncertainty on the galaxy number density
snGal = np.sqrt(zBinCounts) / comovVolume


# save to file
data = np.zeros((nZC, 3))
data[:,0] = zC
data[:,1] = nGal
data[:,2] = snGal
np.savetxt(pathOut + "ngal.txt", data)


'''
# Plot galaxy number density
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.semilogy(zC, nGal, 'b')
#
ax.fill_between(zC, nGal-snGal, nGal+snGal, edgecolor='', facecolor='b')
#
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\bar{n}_\text{gal}$ [(Mpc/h)$^{-3}$]')
#
fig.savefig(pathFig+'ngal_z.pdf', bbox_inches='tight')
plt.show()
plt.clf()
'''


#####################################################################################
# Mean galaxy luminosity and
# total galaxy luminosity density per unit comoving volume

meanLum = np.zeros((nZC, nLines))
sMeanLum = np.zeros((nZC, nLines))
lumDensity = np.zeros((nZC, nLines))
sLumDensity = np.zeros((nZC, nLines))
for iLine in range(nLines):
   meanLum[:,iLine], _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine], statistic='mean', bins=zE)
   #
   sMeanLum[:,iLine], _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine], statistic='std', bins=zE)
   sMeanLum[:,iLine] /= np.sqrt(zBinCounts)
   
   lumDensity[:,iLine], _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine], statistic='sum', bins=zE)
   #
   sLumDensity[:,iLine], _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine], statistic='std', bins=zE)
   sLumDensity[:,iLine] *= np.sqrt(zBinCounts)


lumDensity /= comovVolume[:,np.newaxis]
sLumDensity /= comovVolume[:,np.newaxis]

# save to file
np.savetxt(pathOut + "mean_gal_lum.txt", meanLum)
np.savetxt(pathOut + "total_gal_lum_density.txt", lumDensity)

'''
# Plot mean galaxy luminosity
fig=plt.figure(0, figsize=(10,8))
ax=fig.add_subplot(111)
#
for iLine in range(nLines):
   #ax.semilogy(zC, meanLum[iLine,:], label=cat['lines'][iLine].replace('_', ' '))
   ax.fill_between(zC, meanLum[:,iLine]-sMeanLum[:,iLine], meanLum[:,iLine]+sMeanLum[:,iLine], label=cat['lines'][iLine].replace('_', ' '))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\bar{L}_j^\text{gal}$ [$L_\odot$]')
ax.set_title(r'Mean galaxy luminosity')
#plt.tight_layout()
#
fig.savefig(pathFig+'mean_gal_lum_z.pdf', bbox_inches='tight')
plt.show()
plt.clf()




# Plot luminosity density
fig=plt.figure(0, figsize=(10,8))
ax=fig.add_subplot(111)
#
for iLine in range(nLines):
   #ax.semilogy(zC, lumDensity[iLine,:], label=cat['lines'][iLine].replace('_', ' '))
   ax.fill_between(zC, lumDensity[:,iLine]-sLumDensity[:,iLine], lumDensity[:,iLine]+sLumDensity[:,iLine], label=cat['lines'][iLine].replace('_', ' '))
   #
   # check that lumDensity is meanLum * nGal
   #ax.semilogy(zC, meanLum[iLine,:] * nGal, label=cat['lines'][iLine].replace('_', ' '))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\mathcal{L}_j$ [$L_\odot$ (Mpc/h)$^{-3}$]')
ax.set_title(r'Luminosity density')
#
fig.savefig(pathFig+'luminosity_density_z.pdf', bbox_inches='tight')
#plt.show()
plt.clf()
'''


#####################################################################################
# Covariance of the different line luminosities
# at each redshift

meanLineLum = np.zeros((nZC, nLines))
s2ij = np.zeros((nZC, nLines, nLines))   # s2ij = cov(Li,Lj)/Li/Lj
rij = np.zeros((nZC, nLines, nLines)) # correlation coefficient


for iZ in range(nZC):
   z = zC[iZ]

   for iLine1 in range(nLines):
      line1 = cat['lines'][iLine1]
      line1Lum = cat['line_lum'][:,iLine1]
      line1Lum = line1Lum[np.where((zE[iZ]<=cat['z'])*(cat['z']<zE[iZ+1]))]

      meanLineLum[iZ, iLine1] = np.mean(line1Lum)

      for iLine2 in range(0, iLine1+1):
         line2 = cat['lines'][iLine2]
         line2Lum = cat['line_lum'][:,iLine2]
         line2Lum = line2Lum[np.where((zE[iZ]<=cat['z'])*(cat['z']<zE[iZ+1]))]

         # compute cov and corr coeff
         s2ij[iZ, iLine1, iLine2] = np.cov(np.vstack((line1Lum, line2Lum)))[0,1] / np.mean(line1Lum) / np.mean(line2Lum)
         rij[iZ, iLine1, iLine2] = np.corrcoef(line1Lum, line2Lum)[0,1]


'''
   # plot correlation matrix
   fig=plt.figure(0, figsize=(18,12))
   ax=fig.add_subplot(111)
   #
   mask = np.triu(np.ones((nLines, nLines)), k=1)
   sns.heatmap(rij[iZ,:,:], annot=True, mask=mask, cbar=False)
   ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
   ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
   #
   path = pathFig + "rij_z"+floatExpForm(zC[iZ], round=2)+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()

   # plot relative cov matrix
   fig=plt.figure(0, figsize=(26,14))
   ax=fig.add_subplot(111)
   #
   mask = np.triu(np.ones((nLines, nLines)), k=1)
   sns.heatmap(s2ij[iZ,:,:], annot=True, mask=mask, cbar=False)
   ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
   ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
   #
   path = pathFig + "s2ij_z"+floatExpForm(zC[iZ], round=2)+".pdf"
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
   #plt.show()
'''


# convert to 2d array then save to file
np.savetxt(pathOut + "s2ij.txt", s2ij.flatten())
np.savetxt(pathOut + "rij.txt", rij.flatten())


#####################################################################################
# stellar mass function in each redshift bin


# bins
nME = 31
ME = np.logspace(np.log10(1.e8), np.log10(1.e12), nME, 10.)
nMC = nME - 1

# bin widths
dM = ME[1:] - ME[:-1]
dLog10M = np.log10(ME[1:]) - np.log10(ME[:-1])

# mean mass and lin/log mass function
MC = np.zeros((nZC, nMC))
LinMF = np.zeros((nZC, nMC))
LogMF = np.zeros((nZC, nMC))

for iZ in range(nZC):
   z = zC[iZ]

   # stellar mass
   m = cat['m']
   m = m[np.where((zE[iZ]<=cat['z'])*(cat['z']<zE[iZ+1]))]
   m = 10.**m  # get masses in Msun

   # get central mass in the bin
   MC[iZ,:], _, _ = stats.binned_statistic(m, m, statistic='mean', bins=ME)

   # get histogram
   LinMF[iZ,:], _, _ = stats.binned_statistic(m, m, statistic='count', bins=ME)
   # convert to [(Mpc/h)^{-3}]
   LinMF[iZ,:] /= comovVolume[iZ]
   LogMF[iZ,:] = LinMF[iZ,:].copy()

   # normalize the lin and log mass functions
   LinMF[iZ,:] /= dM
   LogMF[iZ,:] /= dLog10M

'''
# Plot log stellar mass functions
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for iZ in range(nZC):
   z = zC[iZ]
   plt.loglog(MC[iZ,:], LogMF[iZ,:], label=r'$z=$'+str(round(z, 1)))
#
ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'Stellar mass $M$ $[M_\odot]$')
ax.set_ylabel(r'$\frac{dN}{dVd\text{log}_{10}(M/M_\odot)}$ [$(\text{Mpc}/h)^{-3} \text{dex}^{-1}$]')
ax.set_title(r'EGG Stellar mass functions')

plt.show()

# Same, without the little h unit
# to compare with Schreiber+15 fig 3,
# and Schreiber+17 fig 1
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for iZ in range(nZC):
   z = zC[iZ]
   plt.loglog(MC[iZ,:], LogMF[iZ,:] * u.bg.h**3, label=r'$z=$'+str(round(z, 1)))
#
ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
ax.set_xlim((1.e8, 1.e12))
ax.set_ylim((1.e-5, 1.e-1))
ax.set_xlabel(r'Stellar mass $M$ $[M_\odot]$')
ax.set_ylabel(r'$\frac{dN}{dVd\text{log}_{10}(M/M_\odot)}$ [$\text{Mpc}^{-3} \text{dex}^{-1}$]')
ax.set_title(r'EGG Stellar mass functions')

plt.show()
'''



#####################################################################################
# Stellar mass function


# bins
nLE = 31
LE = np.logspace(np.log10(1.e5), np.log10(1.e11), nLE, 10.)
nLC = nLE - 1

# bin widths
dL = LE[1:] - LE[:-1]
dLog10L = np.log10(LE[1:]) - np.log10(LE[:-1])

# mean mass and lin/log mass function
LC = np.zeros((nLines, nZC, nLC))
LinLF = np.zeros((nLines, nZC, nLC))
LogLF = np.zeros((nLines, nZC, nLC))



for iZ in range(nZC):
   z = zC[iZ]

   for iLine1 in range(nLines):
      line1 = cat['lines'][iLine1]

      # line luminosity
      line1Lum = cat['line_lum'][:,iLine1]
      line1Lum = line1Lum[np.where((zE[iZ]<=cat['z'])*(cat['z']<zE[iZ+1]))]

      # get central luminosity in the bin
      LC[iLine1,iZ,:], _, _ = stats.binned_statistic(line1Lum, line1Lum, statistic='mean', bins=LE)

      # get histogram
      LinLF[iLine1,iZ,:], _, _ = stats.binned_statistic(line1Lum, line1Lum, statistic='count', bins=LE)
      # convert to [(Mpc/h)^{-3}]
      LinLF[iLine1,iZ,:] /= comovVolume[iZ]
      LogLF[iLine1,iZ,:] = LinLF[iLine1,iZ,:].copy()

      # normalize the lin and log mass functions
      LinLF[iLine1,iZ,:] /= dL
      LogLF[iLine1,iZ,:] /= dLog10L



# Plot line luminosity functions
#for iLine1 in range(nLines):
for iLine1 in [10]:

   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   for iZ in range(nZC):
      z = zC[iZ]
      plt.loglog(LC[iLine1,iZ,:], LogLF[iLine1,iZ,:], label=r'$z=$'+str(round(z, 1)))
   #
   ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'Line luminosity $[L_\odot]$')
   ax.set_ylabel(r'$\frac{dN}{dVd\text{log}_{10}(L/L_\odot)}$ [$(\text{Mpc}/h)^{-3} \text{dex}^{-1}$]')
   ax.set_title(r'EGG '+cat['lines'][iLine1]+' line luminosity functions')

   plt.show()


# Same, but comparing with Schreiber+in prep, fig A4
#for iLine1 in range(nLines):
for iLine1 in [10]:

   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   for iZ in range(nZC):
      z = zC[iZ]
      plt.loglog(LC[iLine1,iZ,:], LogLF[iLine1,iZ,:] * u.bg.h**3, label=r'$z=$'+str(round(z, 1)))
   #
   ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   ax.set_xlim((10.**7.5, 10.**10.5))
   ax.set_ylim((3.e-6, 3.e-2))
   ax.set_xlabel(r'Line luminosity $[L_\odot]$')
   ax.set_ylabel(r'$\frac{dN}{dVd\text{log}_{10}(L/L_\odot)}$ [$\text{Mpc}^{-3} \text{dex}^{-1}$]')
   ax.set_title(r'EGG '+cat['lines'][iLine1]+' line luminosity functions')

   plt.show()



# Same, but comparing with fig1 in Gong+17
#for iLine1 in range(nLines):
for iLine1 in [10]:

   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   for iZ in range(nZC):
      z = zC[iZ]
      plt.loglog(LC[iLine1,iZ,:] * 3.839e33,  # [Lsun] to [erg/sec] 
            LogLF[iLine1,iZ,:] * u.bg.h**3,  # no change to log-luminosity function
            label=r'$z=$'+str(round(z, 1)))
   #
   ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   ax.set_xlim((1.e39, 1.e44))
   ax.set_ylim((10.**(-6.5), 10.**(-1.)))
   ax.set_xlabel(r'Line luminosity [erg/sec]')
   ax.set_ylabel(r'$\frac{dN}{dVd\text{log}_{10}(L/[\text{erg/sec}])}$ [$\text{Mpc}^{-3} \text{dex}^{-1}$]')
   ax.set_title(r'EGG '+cat['lines'][iLine1]+' line luminosity functions')

   plt.show()






























######################################################################################
## Exact ngaleff for the shot noise auto-spectra
#'''
#nGalEffAuto = np.zeros((nLines, nZC))
#for iLine in range(nLines):
#   summedSquaredLum, _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine]**2, statistic='sum', bins=zE)
#   summedLum, _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine], statistic='sum', bins=zE)
#   nGalEffAuto[iLine,:] = summedLum**2 / summedSquaredLum / comovVolume
#'''
#'''
#fig=plt.figure(0)
#ax=fig.add_subplot(111)
##
#for iLine in range(nLines):
#   ax.semilogy(zC, nGalEffAuto[iLine,:], label=cat['lines'][iLine].replace('_', ' '))
##
#ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
#ax.set_yscale('log', nonposy='clip')
#ax.set_xlabel(r'$z$')
#ax.set_ylabel(r'$\bar{n}^\text{gal eff}_{ii}$ [(Mpc/h)$^{-3}$]')
##
#fig.savefig(pathFig+'ngeff_auto_z.pdf', bbox_inches='tight')
##plt.show()
#plt.clf()
#'''
#
#
#
#
######################################################################################
######################################################################################
## Average over all redshifts:
## line mean, cov and corr coeff
#
######################################################################################
## look at scatter: log-normal?
#
#'''
#for iLine1 in range(nLines):
#   line1 = cat['lines'][iLine1]
#   line1Lum = cat['line_lum'][:,iLine1]
#
#   print 'std dev of '+line1+': '+str(np.std(line1Lum))
#   print 'std dev of log10('+line1+'): '+str(np.std(np.log10(line1Lum)))
#   
#   if plot:
#      #myHistogram(lineFlux[line1][mask], nBins=nBins, lim=None, S2Theory=[], path=None, plot=True, nameLatex=line1, semilogx=True, semilogy=False, doGauss=True)
#      myHistogram(np.log10(line1Lum), nBins=nBins, lim=None, S2Theory=[], path=pathFig+'hist_log10_'+line1+'.pdf', plot=False, nameLatex=r'log$_{10}$('+line1.replace('_',' ')+r')', semilogx=False, semilogy=False, doGauss=True)
#'''
#
######################################################################################
#'''
#meanLineLum = np.zeros(nLines)
#s2ij = np.zeros((nLines, nLines))   # s2ij = cov(Li,Lj)/Li/Lj
#rij = np.zeros((nLines, nLines)) # correlation coefficient
#s2ijlog = np.zeros((nLines, nLines))   # s2ij = cov(log10Li,log10Lj)/ log10Li / log10Lj
#rijlog = np.zeros((nLines, nLines)) # correlation coefficient of log10Li
#
## scatter plot: correlation
#for iLine1 in range(nLines):
#   line1 = cat['lines'][iLine1]
#   line1Lum = cat['line_lum'][:,iLine1]
#
#   meanLineLum[iLine1] = np.mean(line1Lum)
#   print(r'\bar{L}_{'+line1+'} = '+floatExpForm(meanLineLum[iLine1], round=2)+r' L_\odot \\')
#
#   for iLine2 in range(0, iLine1+1):
#      line2 = cat['lines'][iLine2]
#      line2Lum = cat['line_lum'][:,iLine2]
#
#      # compute cov and corr coeff
#      s2ij[iLine1, iLine2] = np.cov(np.vstack((line1Lum, line2Lum)))[0,1] / np.mean(line1Lum) / np.mean(line2Lum)
#      rij[iLine1, iLine2] = np.corrcoef(line1Lum, line2Lum)[0,1]
#      s2ijlog[iLine1, iLine2] = np.cov(np.vstack((np.log10(line1Lum), np.log10(line2Lum))))[0,1] / np.mean(np.log10(line1Lum)) / np.mean(np.log10(line2Lum))
#      rijlog[iLine1, iLine2] = np.corrcoef(np.log10(line1Lum), np.log10(line2Lum))[0,1]
#      
#      
#      if (iLine1<>iLine2):
#         fig=plt.figure(0)
#         ax=fig.add_subplot(111)
#         #
#         ax.errorbar(line1Lum, line2Lum, fmt='.', alpha=0.01)
#         #
#         ax.set_xscale('log', nonposx='clip')
#         ax.set_yscale('log', nonposy='clip')
#         ax.set_xlabel(line1.replace('_',''))
#         ax.set_ylabel(line2.replace('_',' '))
#         #
#         fig.savefig(pathFig+"scatterplot_"+line1+"_"+line2+".pdf", bbox_inches='tight')
#         fig.clf()
#         #plt.show()
#      
#
## Save matrix of cov and corr coeff
#path = pathFig + "s2ij.txt"
#header = "s2ij = cov(Li, Lj)/Li/Lj [dimless]\n" + ', '.join(cat['lines'])
#np.savetxt(path, np.round(s2ij.T, 2), fmt='%.2f', header=header)
##
#path = pathFig + "rij.txt"
#header = "rij, correlation coefficient of line luminosities\n" + ', '.join(cat['lines'])
#np.savetxt(path, np.round(rij.T, 2), fmt='%.2f', header=header)
##
#path = pathFig + "s2ijlog10.txt"
#header = "s2ijlog = cov(log10Li, log10Lj) / log10Li / log10Lj [dimless]\n" + ', '.join(cat['lines'])
#np.savetxt(path, np.round(s2ijlog.T, 2), fmt='%.2f', header=header)
##
#path = pathFig + "rijlog10.txt"
#header = "rijlog, correlation coefficient of the log10 of line luminosities\n" + ', '.join(cat['lines'])
#np.savetxt(path, np.round(rijlog.T, 2), fmt='%.2f', header=header)
#'''
#
#'''
## plot correlation matrix
#fig=plt.figure(0, figsize=(18,12))
#ax=fig.add_subplot(111)
##
#mask = np.triu(np.ones((nLines, nLines)), k=1)
#sns.heatmap(rij, annot=True, mask=mask, cbar=False)
#ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
#ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
##
#path = pathFig + "rij.pdf"
#fig.savefig(path, bbox_inches='tight')
#fig.clf()
##plt.show()
#
## plot correlation matrix in log10
#fig=plt.figure(0, figsize=(18,12))
#ax=fig.add_subplot(111)
##
#mask = np.triu(np.ones((nLines, nLines)), k=1)
#sns.heatmap(rijlog, annot=True, mask=mask, cbar=False)
#ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
#ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
##
#path = pathFig + "rijlog10.pdf"
#fig.savefig(path, bbox_inches='tight')
#fig.clf()
##plt.show()
#
## plot relative cov matrix
#fig=plt.figure(0, figsize=(18,12))
#ax=fig.add_subplot(111)
##
#mask = np.triu(np.ones((nLines, nLines)), k=1)
#sns.heatmap(s2ij, annot=True, mask=mask, cbar=False)
#ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
#ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
##
#path = pathFig + "s2ij.pdf"
#fig.savefig(path, bbox_inches='tight')
#fig.clf()
##plt.show()
#
## plot relative cov matrix in log10
#fig=plt.figure(0, figsize=(18,12))
#ax=fig.add_subplot(111)
##
#mask = np.triu(np.ones((nLines, nLines)), k=1)
#sns.heatmap(s2ijlog, annot=True, mask=mask, cbar=False)
#ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
#ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
##
#path = pathFig + "s2ijlog10.pdf"
#fig.savefig(path, bbox_inches='tight')
#fig.clf()
##plt.show()
#'''
#
######################################################################################
## Approximate effective number density of galaxies
## for auto-correlations
#'''
#nGalEffAutoApprox = np.zeros((nLines, nZC))
#for iLine in range(nLines):
#   nGalEffAutoApprox[iLine,:] = nGal / (1. + s2ij[iLine, iLine])
#'''
#'''
#fig=plt.figure(0)
#ax=fig.add_subplot(111)
##
#for iLine in range(nLines):
#   ax.semilogy(zC, nGalEffAutoApprox[iLine,:], label=cat['lines'][iLine].replace('_', ' '))
##
#ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
#ax.set_yscale('log', nonposy='clip')
#ax.set_xlabel(r'$z$')
#ax.set_ylabel(r'$\bar{n}_\text{gal}(z) / (1 + \sigma^2_{\text{gal }ii})$ [(Mpc/h)$^{-3}$]')
##
#fig.savefig(pathFig+'approx_ngeff_auto_z.pdf', bbox_inches='tight')
#plt.show()
#plt.clf()
#'''
#



