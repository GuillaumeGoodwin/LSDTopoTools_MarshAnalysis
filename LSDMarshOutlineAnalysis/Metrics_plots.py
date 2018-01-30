"""
This is a Python script plots the results of the MarshFinder


"""


#Set up display environment in putty
import matplotlib
matplotlib.use('Agg')

#----------------------------------------------------------------
#1. Load useful Python packages

import os
import sys


import numpy as np
import functools
import math as mt
import cmath
import scipy as sp
import scipy.stats as stats
from datetime import datetime
import cPickle

from pylab import *
import functools

import itertools as itt
from osgeo import gdal, osr
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
from copy import copy


from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import scipy.ndimage as scim

import cPickle as pickle
#------------------------------------------------------------------
# Import the marsh-finding functions
from LSDMarshPlatform_functions import Open_tide_stats
from LSDMarshPlatform_functions import ENVI_raster_binary_to_2d_array
from LSDMarshPlatform_functions import ENVI_raster_binary_from_2d_array
from LSDMarshPlatform_functions import Distribution
from LSDMarshPlatform_functions import Outline
from LSDMarshPlatform_functions import kernel



from New_functions import *




#############################################################################
#------------------------------------------------------------------
# These are some variables and arrays that need to be set up at the beginning

Nodata_value = -9999

Gauges = []
for i in range(1,16):
    if i != 3:
        Gauges.append("WOR_domain_%g" % (i))
                
# This is temporary of course
Gauges = ['WOR_domain_1', 'WOR_domain_2', 'WOR_domain_4', 'WOR_domain_13', 'WOR_domain_14', 'WOR_domain_15']




#domain = pickle.load( open(write_path+write_name+'_domain.pkl', "rb" ) )

##############################################################################    
# Select the plots to draw

Make_simple_outline = True

Draw_combined_pdf = False
Draw_Tidal_prism = False
Draw_scarp_stats = True


tiddir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Tide/"
srcdir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshAnalysis/Input/"
dstdir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshAnalysis/Output/Marsh_metrics/"

if Make_simple_outline == True:
    for gauge in Gauges:
        Platform, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array (srcdir+"%s_Marsh.bil" % (gauge), gauge)
        
        #Platform = Platform [800:1200, 50:450]
        Platform = Platform [250:500, 950:1150]
        

        Outline = True_Outline (Platform, 1, 1, Nodata_value)

        Labeled_outline = label_connected (Outline, Nodata_value)
        
        Simple_outline = select_few_longest (Labeled_outline, Nodata_value, 2)
        
        Creekless_outline = Delete_creeks (Platform, Simple_outline, 11, Nodata_value)

        Creekless_outline[Creekless_outline == 0] = -10
        
        fig=plt.figure("Combined PdF", facecolor='White',figsize=[50,50])

        ax_raw = plt.subplot2grid((1,2),(0,0), colspan=1, rowspan=1, axisbg='white')
        ax_processed = plt.subplot2grid((1,2),(0,1), colspan=1, rowspan=1, axisbg='white')

        ax_raw.imshow (Platform, cmap=plt.cm.gist_earth, interpolation = 'None')
        #ax_processed.imshow (Platform, cmap=plt.cm.gist_earth, interpolation = 'None')
        ax_raw.imshow (Outline, interpolation = 'None', alpha = 0.6)
        ax_processed.imshow (Creekless_outline, cmap=plt.cm.jet, interpolation = 'None', alpha = 1)

        
        
        plt.savefig(dstdir+'TEST2.png')
        
        STOP

    
STOP

if Draw_combined_pdf == True:
    Combined_pdf (sourcedir, destindir, Gauges, tiddir)
if Draw_Tidal_prism == True:
    Tidal_prism (sourcedir, destindir, Gauges, tiddir)
if Draw_scarp_stats == True:
    Scarp_stats (sourcedir, destindir, Gauges, tiddir)
    

