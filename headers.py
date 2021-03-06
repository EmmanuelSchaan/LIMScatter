#from pylab import *


import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import os
import scipy
from scipy import special, optimize, integrate, stats
from scipy.interpolate import UnivariateSpline, RectBivariateSpline, interp1d, interp2d, BarycentricInterpolator
from time import time
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
import matplotlib.colors as mc
from matplotlib.mlab import bivariate_normal
import colorsys
from timeit import timeit
from time import time
from copy import copy
import sys
import seaborn as sns

# parallelizing "map"
# version that works when the function is a class module
from pathos.multiprocessing import ProcessingPool as Pool
#from multiprocess import Pool
# import sharedmem   # library from Yu Feng to fork rather than pickle

# Yu Feng's version of multiprocessing, relying on forking rather than pickling
import sharedmem

import vegas   # for Monte Carlo integration, for CMB lens reconstruction
#import gvar as gv
from astropy.io import fits   # for saving/reeading maps
import colormaps as cmaps  # for viridis and plasma colormaps

# to save wav files
from scipy.io import wavfile

import itertools

# for faster FFT
#import pyfftw
#pyfftw.interfaces.cache.enable() # so subsequent FFTs use the wisdom from the first one
## however, this wisdom is only kept for 0.1sec, or for x seconds if using:
##pyfftw.interfaces.cache.set_keepalive_time(x)


# neutrino masses: https://github.com/ThomasTram/iCLASS/blob/master/neutrinohierarchy.ipynb
# cosmo doc: https://lesgourg.github.io/class-tour/Narbonne.pdf
#from classy import Class

# http://classylss.readthedocs.io/en/stable/index.html
import classylss
import classylss.binding as CLASS

##################################################################################
# for pretty plots

from matplotlib import rc
#rc('font',**{'size':'20','family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
rc('font',**{'size':'22','family':'serif','serif':['CMU serif']})
rc('mathtext', **{'fontset':'cm'})
rc('text', usetex=True)
rc('text.latex', preamble='\usepackage{amsmath}, \usepackage{amssymb}')
#rc('font', size=20)
rc('legend',**{'fontsize':'18'})

# fonty stuffs
#font.serif: CMU Serif
#font.family: serif
#mathtext.fontset: cm
#text.usetex: False
#text.latex.preamble: \usepackage{amsmath}




##################################################################################

import utils
reload(utils)
from utils import *

