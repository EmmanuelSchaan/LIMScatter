from headers import *
from astropy.cosmology import Planck15 as cosmo

plot = True

nBins = 15

pathInput = "../../data/GOALS/"
pathFig = "./figures/GOALS/"
if not os.path.exists(pathFig):
   os.makedirs(pathFig)


#######################################################################

lines = np.array(['CII158', 'NII122', 'OI63', 'OIII88'])
nLines = len(lines)

nGal = 300  # needs to be at least as large as the number of galaxies
galMask = {}
lineFlux = {}
sLineFlux = {}
contFlux = {}
sContFlux = {}
raString = {}
decString = {}
ra = {}
dec = {}
z = {}


# Read all line data
for line in lines:
   lineData = np.genfromtxt(pathInput + "goals_herschel_pacs_"+line+"_linecatalog_HIPEv13.dat", skip_header=2, dtype='str')
   # find out which galaxies have a measured flux
   galId = lineData[:,0].astype(np.int)
   print str(len(galId))+" galaxies have a measured "+line+" line"

   # intialize arrays
   galMask[line] = np.zeros(nGal, dtype='bool')
   lineFlux[line] = np.zeros(nGal)
   sLineFlux[line] = np.zeros(nGal)
   contFlux[line] = np.zeros(nGal)
   sContFlux[line] = np.zeros(nGal)
   raString[line] = np.zeros(nGal, dtype=np.dtype('U12'))
   decString[line] = np.zeros(nGal, dtype=np.dtype('U12'))
   ra[line] = np.zeros(nGal)
   dec[line] = np.zeros(nGal)

   # read data from file
   galMask[line][galId] = True
   lineFlux[line][galId] = lineData[:,17].astype(np.float) * 1.e-17  # [W/m^2]
   sLineFlux[line][galId] = lineData[:,18].astype(np.float) * 1.e-17 # [W/m^2]
   contFlux[line][galId] = lineData[:,19].astype(np.float)  # [mJy]
   sContFlux[line][galId] = lineData[:,20].astype(np.float) # mJy
   raString[line][galId] = lineData[:,6]  # format '00h09m53.32s'
   decString[line][galId] = lineData[:,7]  # format '+25d55m26.2s'

   # Convert RA, dec to degrees
   for iGal in range(nGal):
      if galMask[line][iGal]:
         # remove the h, m, s (format '00 09 53.32')
         raString[line][iGal] = raString[line][iGal].replace('h', ' ').replace('m', ' ').replace('s', '')
         h, m, s = np.array(raString[line][iGal].split(' '), dtype=np.float)
         ra[line][iGal] = (360./24.) * (h + m/60. + s/3600.)   # [deg]

         # remove the h, m, s (format '00 09 53.32')
         decString[line][iGal] = decString[line][iGal].replace('d', ' ').replace('m', ' ').replace('s', '')
         if decString[line][iGal][0]=='-':
            sgn = -1.
         else:
            sgn = 1.
         d, m, s = np.array(decString[line][iGal].split(' '), dtype=np.float)
         #sgn = lambda x: 1.*(x>=0.) - 1.*(x<0.)
         #dec[line][iGal] = sgn(d) * (np.abs(d) + m/60. + s/3600.)   # [deg]
         dec[line][iGal] = sgn * (np.abs(d) + m/60. + s/3600.)   # [deg]


# Combine the ra, dec from all files into one list
ra['union'] = np.zeros(nGal)
dec['union'] = np.zeros(nGal)
galMask['union'] = np.zeros(nGal, dtype=bool)
for iGal in range(nGal):
   counter = 0
   for line in lines:
      if galMask[line][iGal]:
         counter += 1
         ra['union'][iGal] += ra[line][iGal]
         dec['union'][iGal] += dec[line][iGal]
   if counter:
      galMask['union'][iGal] = True
      ra['union'][iGal] /= counter
      dec['union'][iGal] /= counter



#######################################################################
# Read redshifts

refData = np.genfromtxt(pathInput + "goals_galaxies_ned.txt", skip_header=5, dtype='str')
refZ = refData[:,6]
refNGal = len(refZ)
# replace the '...' redshifts with 'nan'
# and convert to float
refZ = np.array([refZ[i].replace('...', 'nan') for i in range(refNGal)], dtype=np.float)
refRaString = refData[:,2]
refDecString = refData[:,3]


# Convert RA, dec to degrees
refRa = np.zeros(refNGal)
refDec = np.zeros(refNGal)
for iGal in range(refNGal):
      # remove the h, m, s (format '00 09 53.32')
      refRaString[iGal] = refRaString[iGal].replace('h', ' ').replace('m', ' ').replace('s', '')
      h, m, s = np.array(refRaString[iGal].split(' '), dtype=np.float)
      refRa[iGal] = (360./24.) * (h + m/60. + s/3600.)   # [deg]

      # remove the h, m, s (format '00 09 53.32')
      refDecString[iGal] = refDecString[iGal].replace('d', ' ').replace('m', ' ').replace('s', '')
      if refDecString[iGal][0]=='-':
         sgn = -1.
      else:
         sgn = 1.
      d, m, s = np.array(refDecString[iGal].split(' '), dtype=np.float)
      #sgn = lambda x: 1.*(x>=0.) - 1.*(x<0.)
      #dec[line][iGal] = sgn(d) * (np.abs(d) + m/60. + s/3600.)   # [deg]
      refDec[iGal] = sgn * (np.abs(d) + m/60. + s/3600.)   # [deg]


# Check that the RA and dec match
galMask['redshift'] = np.zeros(nGal, dtype=bool)
z['redshift'] = np.zeros(nGal)
for iGal in range(nGal):
   if galMask['union'][iGal]:
      diff = (ra['union'][iGal] - refRa)**2 / (1./3600.)**2 # [arcsec^2]
      diff += (dec['union'][iGal] - refDec)**2 / (1./3600.)**2 # [arcsec^2]
      diff = np.sqrt(diff)
      minDiff = np.min(diff)
      iMin = np.where(diff==minDiff)[0][0]
      #print iGal- iMin, minDiff[iGal]
      # Consider you have a match if coordinate distance < 25 arcsec
      # this is because the IFU has 5x5 pixels with 9,4'' on the side of a pixel
      # so the line fluxes are from the center of one of the pixels, chosen somehow.
      # So the match is not expected to be perfect
      if (minDiff < 25.) and (iGal==iMin) and (np.isfinite(refZ[iMin])):   # [arcsec]
         galMask['redshift'][iGal] = True
         z['redshift'][iGal] = refZ[iMin]




########################################################################
# Inspect the redshifts and luminosity distances

# Inspect redshifts
myHistogram(z['redshift'][galMask['redshift']], nBins=nBins, lim=None, S2Theory=[], path=pathFig+'hist_z.pdf', plot=False, nameLatex=r'$z$', semilogx=False, semilogy=False, doGauss=False)
# Inspect luminosity distances
dist = np.array(cosmo.luminosity_distance(z['redshift'][galMask['redshift']]))
myHistogram(dist, nBins=nBins, lim=None, S2Theory=[], path=pathFig+'hist_luminosity_distances.pdf', plot=False, nameLatex=r'$D_L$ [Mpc]', semilogx=False, semilogy=False, doGauss=False)



########################################################################
## Match with redshifts from IRAS RBGS
#
#irasHdu = fits.open(pathInput + "iras_rbgs_vizier.fits")
#irasZ = irasHdu[1].data['cz']
#irasNGal = len(irasZ)
#irasRa = irasHdu[1].data['_RAJ2000'] # [deg]
#irasDec = irasHdu[1].data['_DEJ2000']  # [deg]
#
#
#irasRaString = irasHdu[1].data['RAJ2000']
#irasDecString = irasHdu[1].data['DEJ2000']
## Convert RA, dec to degrees
#MyirasRa = np.zeros(irasNGal)
#MyirasDec = np.zeros(irasNGal)
#for iGal in range(irasNGal):
#   # RA
#   h, m, s = np.array(irasRaString[iGal].split(' '), dtype=np.float)
#   MyirasRa[iGal] = (360./24.) * (h + m/60. + s/3600.)   # [deg]
#   # dec
#   if irasDecString[iGal][0]=='-':
#      sgn = -1.
#   else:
#      sgn = 1.
#   d, m, s = np.array(irasDecString[iGal].split(' '), dtype=np.float)
#   MyirasDec[iGal] = sgn * (np.abs(d) + m/60. + s/3600.)   # [deg]
#
#
#
#
## Match GOALS to IRAS RBGS, to find the galaxy redshifts
#minDiff = -np.ones(nGal)
#for iGal in range(nGal):
#   if galMask['union'][iGal]:
#      diff = (ra['union'][iGal] - irasRa)**2 / (1./3600.)**2
#      diff += (dec['union'][iGal] - irasDec)**2 / (1./3600.)**2
#      #diff = (dec['union'][iGal] - irasDec)**2 / (1./3600.)**2
#      diff = np.sqrt(diff)
#      minDiff[iGal] = np.min(diff)
#
#plt.hist(minDiff, 1001)
#plt.yscale('log')
#plt.show()


#######################################################################
# Convert line fluxes to luminosities


lineLum = {}
sLineLum = {}

# Read all line data
for line in lines:
   lineLum[line] = np.zeros(nGal)
   sLineLum[line] = np.zeros(nGal)

   # keep only objects with a line flux and redshift
   mask = galMask[line]*galMask['redshift']

   # units
   MpcInm = 3.085678e+22   # [Mpc / m]
   LsunInW = 3.827e+26  # [Lsun / W]
   # convert line flux to line lumi [W]
   lineLum[line][mask] = lineFlux[line][mask] * 4. * np.pi * (cosmo.luminosity_distance(z['redshift'][mask]) * MpcInm)**2 # [W]
   lineLum[line][mask] /= 3.827e+26 # [Lsun]
   #
   sLineLum[line][mask] = sLineFlux[line][mask] * 4. * np.pi * (cosmo.luminosity_distance(z['redshift'][mask]) * MpcInm)**2 # [W]
   sLineLum[line][mask] /= 3.827e+26 # [Lsun]




#######################################################################
# Inspect line fluxes


# Plot histograms: is it log normal?
for line1 in lines:
   mask = galMask[line1] * (lineFlux[line1]>0.)
   flux = lineFlux[line1][mask]

   print 'std dev of '+line1+': '+str(np.std(flux))+' W/m^2'
   print 'std dev of log10('+line1+'): '+str(np.std(np.log10(flux)))
   print 'sigma_i^2 = var[Li]/Li = '+str( np.var(flux) / np.mean(flux)**2 )
   
   if plot:
      #myHistogram(lineFlux[line1][mask], nBins=nBins, lim=None, S2Theory=[], path=None, plot=True, nameLatex=line1, semilogx=True, semilogy=False, doGauss=True)
      myHistogram(np.log10(lineFlux[line1][mask]), nBins=nBins, lim=None, S2Theory=[], path=pathFig+'hist_log10_'+line1+'.pdf', plot=False, nameLatex=r'log$_{10}$('+line1+r'/[W/m$^2$])', semilogx=False, semilogy=False, doGauss=True)



# line cov and corr coeff
s2ij = np.zeros((nLines, nLines))   # s2ij = cov(Li,Lj)/Li/Lj
rij = np.zeros((nLines, nLines)) # correlation coefficient

for iLine1 in range(nLines):
   line1 = lines[iLine1]
   for iLine2 in range(iLine1+1):
      line2 = lines[iLine2]
# loop over distinct pairs of lines
#for line1, line2 in list(itertools.combinations(lines, r=2)):
      print "Scatter plot for", line1, line2

      # keep the galaxies that have both lines measured
      # negative line fluxes mean upper limits: discard for now
      mask = galMask[line1] * galMask[line2] * (lineFlux[line1]>0) * (lineFlux[line2]>0)

      flux1 = lineFlux[line1][mask]
      flux2 = lineFlux[line2][mask]
      sFlux1 = sLineFlux[line1][mask]
      sFlux2 = sLineFlux[line2][mask]

      # compute cov and corr coeff
      s2ij[iLine1, iLine2] = np.cov(np.vstack((flux1, flux2)))[0,1] / np.mean(flux1) / np.mean(flux2)
      rij[iLine1, iLine2] = np.corrcoef(flux1, flux2)[0,1]

      print 'corr coeff of '+line1+' and '+line2+': '+str(rij[iLine1, iLine2])
      print 'sigma_[i,j}^2 = cov[Li, Lj]/LiLj ='+str(s2ij[iLine1, iLine2])
      
      if plot and (iLine1<>iLine2):
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.errorbar(flux2, flux1, xerr=sFlux2, yerr=sFlux1, fmt='bo')
         #
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlabel(line2 + r' [W/m$^2$]')
         ax.set_ylabel(line1 + r' [W/m$^2$]')
         #
         fig.savefig(pathFig+"scatterplot_"+line1+"_"+line2+".pdf", bbox_inches='tight')
         fig.clf()
         #plt.show()


#plot = True

# Save matrix of cov and corr coeff
path = pathFig + "s2ij.txt"
header = "s2ij = cov(Li, Lj)/Li/Lj [dimless]\n" + ', '.join(lines)
np.savetxt(path, np.round(s2ij.T, 2), fmt='%.2f', header=header)
#
path = pathFig + "rij.txt"
header = "rij, correlation coefficient of line luminosities\n" + ', '.join(lines)
np.savetxt(path, np.round(rij.T, 2), fmt='%.2f', header=header)


if plot:
   # plot correlation matrix
   fig=plt.figure(0, figsize=(14,10))
   ax=fig.add_subplot(111)
   #
   sns.heatmap(rij.T, annot=True)
   ax.set_xticklabels(lines, rotation=45)
   ax.set_yticklabels(lines)
   #
   path = pathFig + "rij.pdf"
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
   #plt.show()

   # plot relative cov matrix
   fig=plt.figure(0, figsize=(14,10))
   ax=fig.add_subplot(111)
   #
   sns.heatmap(s2ij.T, annot=True)
   ax.set_xticklabels(lines, rotation=45)
   ax.set_yticklabels(lines)
   #
   path = pathFig + "s2ij.pdf"
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
   #plt.show()




#######################################################################
# Inspect line luminosities




# Plot histograms: is it log normal?
for line1 in lines:
   mask = galMask[line1] * (lineLum[line1]>0.) * galMask['redshift']
   lum = lineLum[line1][mask]

   print 'std dev of '+line1+': '+str(np.std(lum))+' W/m^2'
   print 'std dev of log10('+line1+'): '+str(np.std(np.log10(lum)))
   print 'sigma_i^2 = var[Li]/Li = '+str( np.var(lum) / np.mean(lum)**2 )

   if plot:
      #myHistogram(lineLum[line1][mask], nBins=nBins, lim=None, S2Theory=[], path=None, plot=True, nameLatex=line1, semilogx=True, semilogy=False, doGauss=True)
      myHistogram(np.log10(lineLum[line1][mask]), nBins=nBins, lim=None, S2Theory=[], path=pathFig+'hist_log10_lum_'+line1+'.pdf', plot=False, nameLatex=r'log$_{10}$('+line1+r'/[$L_\odot$])', semilogx=False, semilogy=False, doGauss=True)



# line cov and corr coeff
s2ij = np.zeros((nLines, nLines))   # s2ij = cov(Li,Lj)/Li/Lj
rij = np.zeros((nLines, nLines)) # correlation coefficient

for iLine1 in range(nLines):
   line1 = lines[iLine1]
   for iLine2 in range(iLine1+1):
      line2 = lines[iLine2]
# loop over distinct pairs of lines
#for line1, line2 in list(itertools.combinations(lines, r=2)):
      print "Scatter plot for", line1, line2

      # keep the galaxies that have both lines measured
      # negative line lumes mean upper limits: discard for now
      mask = galMask[line1] * galMask[line2] * (lineLum[line1]>0) * (lineLum[line2]>0)
      mask *= galMask['redshift']

      lum1 = lineLum[line1][mask]
      lum2 = lineLum[line2][mask]
      sLum1 = sLineLum[line1][mask]
      sLum2 = sLineLum[line2][mask]

      # compute cov and corr coeff
      s2ij[iLine1, iLine2] = np.cov(np.vstack((lum1, lum2)))[0,1] / np.mean(lum1) / np.mean(lum2)
      rij[iLine1, iLine2] = np.corrcoef(lum1, lum2)[0,1]

      print 'corr coeff of '+line1+' and '+line2+': '+str(rij[iLine1, iLine2])
      print 'sigma_[i,j}^2 = cov[Li, Lj]/LiLj ='+str(s2ij[iLine1, iLine2])

      if plot and (iLine1<>iLine2):
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.errorbar(lum2, lum1, xerr=sLum2, yerr=sLum1, fmt='bo')
         #
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlabel(line2 + r' [$L_\odot$]')
         ax.set_ylabel(line1 + r' [$L_\odot$]')
         #
         fig.savefig(pathFig+"scatterplot_lum_"+line1+"_"+line2+".pdf", bbox_inches='tight')
         fig.clf()
         #plt.show()


#plot = True

# Save matrix of cov and corr coeff
path = pathFig + "s2ij_lum.txt"
header = "s2ij = cov(Li, Lj)/Li/Lj [dimless]\n" + ', '.join(lines)
np.savetxt(path, np.round(s2ij.T, 2), fmt='%.2f', header=header)
#
path = pathFig + "rij_lum.txt"
header = "rij, correlation coefficient of line luminosities\n" + ', '.join(lines)
np.savetxt(path, np.round(rij.T, 2), fmt='%.2f', header=header)


if plot:
   # plot correlation matrix
   fig=plt.figure(0, figsize=(14,10))
   ax=fig.add_subplot(111)
   #
   mask = np.triu(np.ones((nLines, nLines)), k=1)
   sns.heatmap(rij, annot=True, mask=mask, cbar=False)
   ax.set_xticklabels(lines, rotation=45)
   ax.set_yticklabels(lines)
   #
   path = pathFig + "rij_lum.pdf"
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
   #plt.show()

   # plot relative cov matrix
   fig=plt.figure(0, figsize=(14,10))
   ax=fig.add_subplot(111)
   #
   mask = np.triu(np.ones((nLines, nLines)), k=1)
   sns.heatmap(s2ij, annot=True, mask=mask, cbar=False)
   ax.set_xticklabels(lines, rotation=45)
   ax.set_yticklabels(lines)
   #
   path = pathFig + "s2ij_lum.pdf"
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
   #plt.show()

