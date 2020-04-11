from headers import *



pathInput = "../../data/GOALS/"

pathFig = "./figures/GOALS/"
if not os.path.exists(pathFig):
   os.makedirs(pathFig)




lines = np.array(['CII158', 'NII122', 'OI63', 'OIII88'])

nGal = 300  # needs to be at least as large as the number of galaxies
galMask = {}
lineFlux = {}
sLineFlux = {}
contFlux = {}
sContFlux = {}


# Read all data
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

   # read data from file
   galMask[line][galId] = True
   lineFlux[line][galId] = lineData[:,17].astype(np.float)  # [W/m^2]
   sLineFlux[line][galId] = lineData[:,18].astype(np.float) # [W/m^2]
   contFlux[line][galId] = lineData[:,19].astype(np.float)  # [mJy]
   sContFlux[line][galId] = lineData[:,20].astype(np.float) # mJy






# Plot histograms: is it log normal?
for line1 in lines:
   nBins = 15
   mask = galMask[line1] * (lineFlux[line1]>0.)
   flux = lineFlux[line1][mask]

   print 'std dev of '+line1+': '+str(np.std(flux))+' W/m^2'
   print 'std dev of log10('+line1+'): '+str(np.std(np.log10(flux)))
   print 'sigma_i^2 = var[Li]/Li = '+str( np.std(flux) / np.mean(flux) )

   #myHistogram(lineFlux[line1][mask], nBins=nBins, lim=None, S2Theory=[], path=None, plot=True, nameLatex=line1, semilogx=True, semilogy=False, doGauss=True)
   myHistogram(np.log10(lineFlux[line1][mask]), nBins=nBins, lim=None, S2Theory=[], path=pathFig+'hist_log10_'+line1+'.pdf', plot=False, nameLatex=r'log$_{10}$('+line1+r'/[W/m$^2$])', semilogx=False, semilogy=False, doGauss=True)


# loop over distinct pairs of lines
for line1, line2 in list(itertools.combinations(lines, r=2)):
   print "Scatter plot for", line1, line2

   # keep the galaxies that have both lines measured
   # negative line fluxes mean upper limits: discard for now
   mask = galMask[line1] * galMask[line2] * (lineFlux[line1]>0) * (lineFlux[line2]>0)

   flux1 = lineFlux[line1][mask]
   flux2 = lineFlux[line2][mask]
   sFlux1 = sLineFlux[line1][mask]
   sFlux2 = sLineFlux[line2][mask]


   print 'corr coeff of '+line1+' and '+line2+': '+str(np.corrcoef(flux1, flux2)[0,1])
   #print 'corr coeff of log10('+line1+') and log10('+line2+'): '+str(np.corrcoef(np.log10(flux1), np.log10(flux2))[0,1])
   print 'sigma_[i,j}^2 = cov[Li, Lj]/LiLj ='+str( np.cov(np.vstack((flux1, flux2)))[0,1] / np.mean(flux1) / np.mean(flux2) )
   
   #my2dHistogram(flux1, flux2, nBins=(nBins, nBins), limx=None, limy=None, limc=None, fTheory=[], path=None, plot=True, nameLatexX=line1, nameLatexY=line2, logx=True, logy=True, logColor=False, cmap=cmaps.viridis)

   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   ax.errorbar(flux1, flux2, xerr=sFlux1, yerr=sFlux2, fmt='bo')
   #
   ax.set_xscale('log', nonposx='clip')
   ax.set_yscale('log', nonposy='clip')
   ax.set_xlabel(line1 + r' [W/m$^2$]')
   ax.set_ylabel(line2 + r' [W/m$^2$]')
   #
   fig.savefig(pathFig+"scatterplot_"+line1+"_"+line2+".pdf", bbox_inches='tight')
   fig.clf()
   #plt.show()



