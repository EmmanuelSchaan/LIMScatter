from headers import *

import pycolfits

pathCat = "./output/EGG/catalog_EGG.fits"
pathFig = "./figures/EGG/"

if not os.path.exists(pathFig):
   osmakedirs(pathFig)


# Read the catalog in memory
cat = pycolfits.readfrom(pathCat, lower_case=True)

# galaxy IDs
#cat['id']

# galaxy lines
#print cat['lines']

nLines = len(cat['lines'])
nBins = 51


# look at scatter: log-normal?
for iLine1 in range(nLines):
   line1 = cat['lines'][iLine1]
   line1Lum = cat['line_lum'][:,iLine1]

   print 'std dev of '+line1+': '+str(np.std(line1Lum))
   print 'std dev of log10('+line1+'): '+str(np.std(np.log10(line1Lum)))

   #myHistogram(lineFlux[line1][mask], nBins=nBins, lim=None, S2Theory=[], path=None, plot=True, nameLatex=line1, semilogx=True, semilogy=False, doGauss=True)
   myHistogram(np.log10(line1Lum), nBins=nBins, lim=None, S2Theory=[], path=pathFig+'hist_log10_'+line1+'.pdf', plot=False, nameLatex=r'log$_{10}$('+line1.replace('_',' ')+r')', semilogx=False, semilogy=False, doGauss=True)




# scatter plot: correlation
for iLine1 in range(nLines):
   line1 = cat['lines'][iLine1]
   line1Lum = cat['line_lum'][:,iLine1]

   for iLine2 in range(0, iLine1):
      line2 = cat['lines'][iLine2]
      line2Lum = cat['line_lum'][:,iLine2]

      
      print 'corr coeff of '+line1+' and '+line2+': '+str(np.corrcoef(line1Lum, line2Lum)[0,1])
      print 'corr coeff of log10('+line1+') and log10('+line2+'): '+str(np.corrcoef(np.log10(line1Lum), np.log10(line2Lum))[0,1])
      
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







