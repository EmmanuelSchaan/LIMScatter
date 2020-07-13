from headers import *
import pycolfits
from universe import *

#####################################################################################

u = UnivPlanck15()


#####################################################################################
plot = False


pathCat = "./output/EGG/catalog_EGG_LIM.fits"
pathFig = "./figures/EGG/"

if not os.path.exists(pathFig):
   osmakedirs(pathFig)


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
# Explore basic properties
'''
# Redshifts
myHistogram(cat['z'], nBins=101, lim=None, S2Theory=[], path=pathFig+'z.pdf', plot=plot, nameLatex=r'$z$', semilogx=False, semilogy=True, doGauss=False)

# Stellar masses
myHistogram(cat['m'], nBins=101, lim=None, S2Theory=[], path=pathFig+'log10mstellar.pdf', plot=plot, nameLatex=r'$\text{log}_{10}\left(M_\star/M_\odot\right)$', semilogx=False, semilogy=True, doGauss=False)
'''

#####################################################################################

# Define redshift bins
nZE = 101
nZC = nZE-1
# z bin edges
zE = np.linspace(np.min(cat['z']), np.max(cat['z']), nZE)
dZ = zE[1:] - zE[:-1]
# get effective bin centers and bin counts
zC, zE, zBinIndices = stats.binned_statistic(cat['z'], cat['z'], statistic='mean', bins=zE)
zBinCounts, _, _ = stats.binned_statistic(cat['z'], cat['z'], statistic='count', bins=zE)

# Comoving volume per redshift bin
skyArea = 1. * (np.pi/180.)**2   # [sr]
chi = u.bg.comoving_distance(zC) # [Mpc/h]
hubble = u.hubble(zC)   # [km/s/(Mpc/h)]
comovVolume = chi**2 * skyArea * (3.e5/hubble) * dZ

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
#plt.show()
plt.clf()
'''


#####################################################################################
# Total galaxy luminosity per unit comoving volume

meanLum = np.zeros((nLines, nZC))
sMeanLum = np.zeros((nLines, nZC))
lumDensity = np.zeros((nLines, nZC))
sLumDensity = np.zeros((nLines, nZC))
for iLine in range(nLines):
   meanLum[iLine,:], _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine], statistic='mean', bins=zE)
   #
   sMeanLum[iLine,:], _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine], statistic='std', bins=zE)
   sMeanLum[iLine,:] /= np.sqrt(zBinCounts)
   
   lumDensity[iLine,:], _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine], statistic='sum', bins=zE)
   #
   sLumDensity[iLine,:], _, _ = stats.binned_statistic(cat['z'], cat['line_lum'][:,iLine], statistic='std', bins=zE)
   sLumDensity[iLine,:] *= np.sqrt(zBinCounts)

lumDensity /= comovVolume
sLumDensity /= comovVolume
   
'''
# Plot luminosity density
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for iLine in range(nLines):
   #ax.semilogy(zC, meanLum[iLine,:], label=cat['lines'][iLine].replace('_', ' '))
   ax.fill_between(zC, meanLum[iLine,:]-sMeanLum[iLine,:], meanLum[iLine,:]+sMeanLum[iLine,:], label=cat['lines'][iLine].replace('_', ' '))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'Mean gal. lum. $\bar{L}_j$ [$L_\odot$]')
#
fig.savefig(pathFig+'mean_gal_lum_z.pdf', bbox_inches='tight')
#plt.show()
plt.clf()
'''


'''
# Plot luminosity density
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for iLine in range(nLines):
   #ax.semilogy(zC, lumDensity[iLine,:], label=cat['lines'][iLine].replace('_', ' '))
   ax.fill_between(zC, lumDensity[iLine,:]-sLumDensity[iLine,:], lumDensity[iLine,:]+sLumDensity[iLine,:], label=cat['lines'][iLine].replace('_', ' '))
   #
   # check that lumDensity is meanLum * nGal
   #ax.semilogy(zC, meanLum[iLine,:] * nGal, label=cat['lines'][iLine].replace('_', ' '))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'Luminosity density $\mathcal{L}_j$ [$L_\odot$ (Mpc/h)$^{-3}$]')
#
fig.savefig(pathFig+'luminosity_density_z.pdf', bbox_inches='tight')
#plt.show()
plt.clf()
'''

#####################################################################################
# look at scatter: log-normal?

'''
for iLine1 in range(nLines):
   line1 = cat['lines'][iLine1]
   line1Lum = cat['line_lum'][:,iLine1]

   print 'std dev of '+line1+': '+str(np.std(line1Lum))
   print 'std dev of log10('+line1+'): '+str(np.std(np.log10(line1Lum)))
   
   if plot:
      #myHistogram(lineFlux[line1][mask], nBins=nBins, lim=None, S2Theory=[], path=None, plot=True, nameLatex=line1, semilogx=True, semilogy=False, doGauss=True)
      myHistogram(np.log10(line1Lum), nBins=nBins, lim=None, S2Theory=[], path=pathFig+'hist_log10_'+line1+'.pdf', plot=False, nameLatex=r'log$_{10}$('+line1.replace('_',' ')+r')', semilogx=False, semilogy=False, doGauss=True)
'''


#####################################################################################
#####################################################################################
# Average over all redshifts:
# line mean, cov and corr coeff


meanLineLum = np.zeros(nLines)
s2ij = np.zeros((nLines, nLines))   # s2ij = cov(Li,Lj)/Li/Lj
rij = np.zeros((nLines, nLines)) # correlation coefficient
s2ijlog = np.zeros((nLines, nLines))   # s2ij = cov(log10Li,log10Lj)/ log10Li / log10Lj
rijlog = np.zeros((nLines, nLines)) # correlation coefficient of log10Li

# scatter plot: correlation
for iLine1 in range(nLines):
   line1 = cat['lines'][iLine1]
   line1Lum = cat['line_lum'][:,iLine1]

   meanLineLum[iLine1] = np.mean(line1Lum)
   print(r'\bar{L}_{'+line1+'} = '+floatExpForm(meanLineLum[iLine1], round=2)+r' L_\odot \\')

   for iLine2 in range(0, iLine1+1):
      line2 = cat['lines'][iLine2]
      line2Lum = cat['line_lum'][:,iLine2]

      # compute cov and corr coeff
      s2ij[iLine1, iLine2] = np.cov(np.vstack((line1Lum, line2Lum)))[0,1] / np.mean(line1Lum) / np.mean(line2Lum)
      rij[iLine1, iLine2] = np.corrcoef(line1Lum, line2Lum)[0,1]
      s2ijlog[iLine1, iLine2] = np.cov(np.vstack((np.log10(line1Lum), np.log10(line2Lum))))[0,1] / np.mean(np.log10(line1Lum)) / np.mean(np.log10(line2Lum))
      rijlog[iLine1, iLine2] = np.corrcoef(np.log10(line1Lum), np.log10(line2Lum))[0,1]
      
      '''
      if (iLine1<>iLine2):
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.errorbar(line1Lum, line2Lum, fmt='.', alpha=0.01)
         #
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlabel(line1.replace('_',''))
         ax.set_ylabel(line2.replace('_',' '))
         #
         fig.savefig(pathFig+"scatterplot_"+line1+"_"+line2+".pdf", bbox_inches='tight')
         fig.clf()
         #plt.show()
      '''

# Save matrix of cov and corr coeff
path = pathFig + "s2ij.txt"
header = "s2ij = cov(Li, Lj)/Li/Lj [dimless]\n" + ', '.join(cat['lines'])
np.savetxt(path, np.round(s2ij.T, 2), fmt='%.2f', header=header)
#
path = pathFig + "rij.txt"
header = "rij, correlation coefficient of line luminosities\n" + ', '.join(cat['lines'])
np.savetxt(path, np.round(rij.T, 2), fmt='%.2f', header=header)
#
path = pathFig + "s2ijlog10.txt"
header = "s2ijlog = cov(log10Li, log10Lj) / log10Li / log10Lj [dimless]\n" + ', '.join(cat['lines'])
np.savetxt(path, np.round(s2ijlog.T, 2), fmt='%.2f', header=header)
#
path = pathFig + "rijlog10.txt"
header = "rijlog, correlation coefficient of the log10 of line luminosities\n" + ', '.join(cat['lines'])
np.savetxt(path, np.round(rijlog.T, 2), fmt='%.2f', header=header)

'''
# plot correlation matrix
fig=plt.figure(0, figsize=(18,12))
ax=fig.add_subplot(111)
#
mask = np.triu(np.ones((nLines, nLines)), k=1)
sns.heatmap(rij, annot=True, mask=mask, cbar=False)
ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
#
path = pathFig + "rij.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()
#plt.show()

# plot correlation matrix in log10
fig=plt.figure(0, figsize=(18,12))
ax=fig.add_subplot(111)
#
mask = np.triu(np.ones((nLines, nLines)), k=1)
sns.heatmap(rijlog, annot=True, mask=mask, cbar=False)
ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
#
path = pathFig + "rijlog10.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()
#plt.show()

# plot relative cov matrix
fig=plt.figure(0, figsize=(18,12))
ax=fig.add_subplot(111)
#
mask = np.triu(np.ones((nLines, nLines)), k=1)
sns.heatmap(s2ij, annot=True, mask=mask, cbar=False)
ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
#
path = pathFig + "s2ij.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()
#plt.show()

# plot relative cov matrix in log10
fig=plt.figure(0, figsize=(18,12))
ax=fig.add_subplot(111)
#
mask = np.triu(np.ones((nLines, nLines)), k=1)
sns.heatmap(s2ijlog, annot=True, mask=mask, cbar=False)
ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
#
path = pathFig + "s2ijlog10.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()
#plt.show()
'''

#####################################################################################
#####################################################################################
# Approximate effective number density of galaxies
# for auto-correlations


nGalEffAuto = np.zeros((nLines, nZC))
for iLine in range(nLines):
   nGalEffAuto[iLine,:] = nGal / (1. + s2ij[iLine, iLine])

fig=plt.figure(0)
ax=fig.add_subplot(111)
#
for iLine in range(nLines):
   ax.semilogy(zC, nGalEffAuto[iLine,:], label=cat['lines'][iLine].replace('_', ' '))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$\bar{n}_\text{gal}(z) / (1 + \sigma^2_{\text{gal }ij})$')
#
fig.savefig(pathFig+'approx_ngeff_auto_z.pdf', bbox_inches='tight')
#plt.show()
plt.clf()





