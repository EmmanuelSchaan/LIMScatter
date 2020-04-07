

# Script to generate the mock catalog using EGG


# generate the catalog using EGG:
egg-gencat verbose maglim=27 selection_band=hst-f160w area=0.08 out=./output/EGG/catalog_EGG.fits
# get the output in an ASCII file
#egg-gencat verbose maglim=27 selection_band=hst-f160w area=0.08 out=./output/EGG/catalog_EGG.fits ascii





## For visualization with SkyMaker
## convert the output to a SkyMaker input
#egg-2skymaker cat=./output/EGG/catalog_EGG.fits verbose band=hst-f160w template=goodss-hst-f160w.conf
## run SkyMaker to generate a mock image
#sky ./output/EGG/catalog_EGG-hst-f160w.cat -c ./output/EGG/catalog_EGG-hst-f160w-sky.conf




# Input parameters for egg-gencat
#area
#maglim
#selection_band
#mmin
#mmax
#zmin
#zmax
#save_sed to save the full SED of each galaxy
#ascii

# line info available:
#(c2_157, n2_205, c1_609, co10, co21, co32, co43, co54, co65, co76, halpha, hbeta, hgamma, hdelta, n2_6583, n2_6548, o3_5007, o3_4959, o2_3727, lyalpha)

# To list the available bands for galaxy selection:
#egg-gencat list_bands
