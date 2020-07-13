# Script to generate the mock catalog using EGG


# generate the catalog using EGG:
#egg-gencat verbose maglim=27 selection_band=hst-f160w area=0.08 out=./output/EGG/catalog_EGG.fits
# get the output in an ASCII file
#egg-gencat verbose maglim=27 selection_band=hst-f160w area=0.08 out=./output/EGG/catalog_EGG.fits ascii

# Selection relevant for LIM
# remove unwanted output positions to reduce file size
# area=2 sq deg is probably sufficient to have enough galaxies,
# since all we want is mean and variance of the luminosities
# keep the default redshift range z=0.05-10.5
egg-gencat verbose area=1. mmin=8 mmax=14 zmin=0.05 zmax=10. no_pos no_clust out=./output/EGG/catalog_EGG_LIM.fits


## For visualization with SkyMaker
## convert the output to a SkyMaker input
#egg-2skymaker cat=./output/EGG/catalog_EGG.fits verbose band=hst-f160w template=goodss-hst-f160w.conf
## run SkyMaker to generate a mock image
#sky ./output/EGG/catalog_EGG-hst-f160w.cat -c ./output/EGG/catalog_EGG-hst-f160w-sky.conf




# Input parameters for egg-gencat
#area # [sq deg]
#maglim
#selection_band
#mmin # to be specified as log10(Mmin/Msun). This disables the limiting magnitude
#mmax # to be specified as log10(Mmin/Msun)
#zmin # default is z=0.05
#zmax # default is z=10.5
#save_sed to save the full SED of each galaxy
#ascii

# line info available:
#(c2_157, n2_205, c1_609, co10, co21, co32, co43, co54, co65, co76, halpha, hbeta, hgamma, hdelta, n2_6583, n2_6548, o3_5007, o3_4959, o2_3727, lyalpha)

# To list the available bands for galaxy selection:
#egg-gencat list_bands
