from headers import *
import pycolfits



plot = False


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
   
   if plot:
      #myHistogram(lineFlux[line1][mask], nBins=nBins, lim=None, S2Theory=[], path=None, plot=True, nameLatex=line1, semilogx=True, semilogy=False, doGauss=True)
      myHistogram(np.log10(line1Lum), nBins=nBins, lim=None, S2Theory=[], path=pathFig+'hist_log10_'+line1+'.pdf', plot=False, nameLatex=r'log$_{10}$('+line1.replace('_',' ')+r')', semilogx=False, semilogy=False, doGauss=True)



# line cov and corr coeff
s2ij = np.zeros((nLines, nLines))   # s2ij = cov(Li,Lj)/Li/Lj
rij = np.zeros((nLines, nLines)) # correlation coefficient

# scatter plot: correlation
for iLine1 in range(nLines):
   line1 = cat['lines'][iLine1]
   line1Lum = cat['line_lum'][:,iLine1]

   for iLine2 in range(0, iLine1+1):
      line2 = cat['lines'][iLine2]
      line2Lum = cat['line_lum'][:,iLine2]

      # compute cov and corr coeff
      s2ij[iLine1, iLine2] = np.cov(np.vstack((line1Lum, line2Lum)))[0,1] / np.mean(line1Lum) / np.mean(line2Lum)
      rij[iLine1, iLine2] = np.corrcoef(line1Lum, line2Lum)[0,1]
      
      if plot and (iLine1<>iLine2):
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


# Save matrix of cov and corr coeff
path = pathFig + "s2ij.txt"
header = "s2ij = cov(Li, Lj)/Li/Lj [dimless]\n" + ', '.join(cat['lines'])
np.savetxt(path, np.round(s2ij.T, 2), fmt='%.2f', header=header)
#
path = pathFig + "rij.txt"
header = "rij, correlation coefficient of line luminosities\n" + ', '.join(cat['lines'])
np.savetxt(path, np.round(rij.T, 2), fmt='%.2f', header=header)

plot=True

if plot:
   # plot correlation matrix
   fig=plt.figure(0, figsize=(18,12))
   ax=fig.add_subplot(111)
   #
   sns.heatmap(rij.T, annot=True)
   ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
   ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
   #
   path = pathFig + "rij.pdf"
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
   #plt.show()


   # plot relative cov matrix
   fig=plt.figure(0, figsize=(18,12))
   ax=fig.add_subplot(111)
   #
   sns.heatmap(s2ij.T, annot=True)
   ax.set_xticklabels(cat['lines'].replace('_', ' '), rotation=45)
   ax.set_yticklabels(cat['lines'].replace('_', ' '), rotation=0)
   #
   path = pathFig + "s2ij.pdf"
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
   #plt.show()


#if plot:
#   # plot correlation matrix
#   fig=plt.figure(0, figsize=(14,10))
#   ax=fig.add_subplot(111)
#   #
#   sns.heatmap(rij.T, annot=True)
#   ax.set_xticklabels(lines, rotation=45)
#   ax.set_yticklabels(lines)
#   #
#   path = pathFig + "rij.pdf"
#   fig.savefig(path, bbox_inches='tight')
#   fig.clf()
#   #plt.show()
#
#   # plot relative cov matrix
#   fig=plt.figure(0, figsize=(14,10))
#   ax=fig.add_subplot(111)
#   #
#   sns.heatmap(s2ij.T, annot=True)
#   ax.set_xticklabels(lines, rotation=45)
#   ax.set_yticklabels(lines)
#   #
#   path = pathFig + "s2ij.pdf"
#   fig.savefig(path, bbox_inches='tight')
#   fig.clf()
#   #plt.show()
#
