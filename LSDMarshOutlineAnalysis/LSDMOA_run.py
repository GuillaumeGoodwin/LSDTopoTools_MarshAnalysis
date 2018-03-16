"""
LSDMOA_run.py

This file loads the inputs, launches the analysis and saves the outputs of the MOA.

Please read the README and the instructions in this script before you run it.

Authors: Guillaume CH Goodwin and Simon Marius Mudd

"""

#------------------------------------------------------------------
#0. Set up display environment if you are working on a terminal with no GUI.
import matplotlib
matplotlib.use('Agg')

#------------------------------------------------------------------

# Useful Python packages
import pandas as bb
import numpy as np
import cPickle as pickle
import timeit
import os

from LSDMOA_classes import *
from LSDMOA_functions import *




def MarshOutlineAnalysis(Input_dir =  "/Example_Data/",
            Output_dir = "/Example_Data/Output/",
            Site = ["FEL_DEM_clip"], opt1 = -2.0, opt2 = 0.85, opt3 = 8.0):
    """
    This function wraps all the marsh ID scripts in one location

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites
        opt1 (flt): first optimisation
        opt2 (flt): 2nd optimisation
        opt3 (flt): 3rd optimisation
        compare_with_digitised_marsh (bool): If true, this will compare the data with a digitised marsh platform

    Author:
        GCHG, Modified by SMM 02/10/2017
    """
    #------------------------------------------------------------------
    # Timing the run
    Start = timeit.default_timer()
    print("\nWelcome to the MOA programme!\n")
    print("I am opening the input files in: "+Input_dir)

    # Set the value for empty DEM cells
    Nodata_value = -1000


    for site in Site:
        print("Loading input data from site: "+site)
        # NB: When loading input data, please make sure the naming convention shown here is respected.

        print(" Loading DEM")
        DEM_fname = site+"_DEM.bil"
        DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (Input_dir+DEM_fname)
        print DEM.shape

        print " Loading Marsh Platform"
        # Make sure we have the right slope file name
        marsh_fname = site+"_Marsh.bil"
        Marsh, post_Marsh, envidata_Marsh =  ENVI_raster_binary_to_2d_array (Input_dir+marsh_fname)
        print Marsh.shape



        # Here begins the actual MOA
        print "\nProducing Lines"
        #Make a proper object for the marsh
        Marsh_object = Marsh_platform(Marsh.shape[0], Marsh.shape[1])
        Marsh_object = Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = True)
        Marsh_DEM = Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = False)
        Marsh_labels = Marsh_object.label_connected (Nodata_value)

        if not os.path.isfile(Output_dir+'Marsh_metrics/'+str(site)+'_Out.pkl'):
            print 'Extracting outlines from array'
            Outlines = Marsh_labels.extract_outlines()
            Outlines.Save_as_pickle(Output_dir+'Marsh_metrics/',str(site)+'_Out.pkl')
        else:
            print 'Loading outlines from Pickle'
            Outlines = pickle.load( open( Output_dir+'Marsh_metrics/'+str(site)+'_Out.pkl', "rb" ) )
        Outlines.save_to_shp (envidata_DEM, DEM, Output_dir+'Shapefiles/', site)


        print "\nProducing Transects"
        All_transects = Outlines.Polyline_transects(10,20, envidata_DEM, DEM, Output_dir, site)
        All_transects = All_transects.get_attribute_from_basemap (20, DEM, 'Z', Nodata_value)
        All_transects = All_transects.get_attribute_from_basemap (1, Marsh_object, 'Marsh', Nodata_value)
        All_transects_mean, All_transects_stdev, Big_mean, Big_stdev, Bigtransect = All_transects.Polyline_stats()

        print "\nSaving outputs"

        All_transects.Save_as_pickle(Output_dir+'Marsh_metrics/',str(site)+'_Tr.pkl')
        All_transects_mean.Save_as_pickle(Output_dir+'Marsh_metrics/',str(site)+'_mean.pkl')
        All_transects_stdev.Save_as_pickle(Output_dir+'Marsh_metrics/',str(site)+'_stdev.pkl')

        Bigtransect.Save_as_pickle(Output_dir+'Marsh_metrics/',str(site)+'_BigTr.pkl')
        Big_mean.Save_as_pickle(Output_dir+'Marsh_metrics/',str(site)+'_Bigmean.pkl')
        Big_stdev.Save_as_pickle(Output_dir+'Marsh_metrics/',str(site)+'_Bigstdev.pkl')





        Stop = timeit.default_timer()
        print '\nAnalysis runtime was: ', Stop - Start , 's'


























        """Marsh = np.ones((50,60), dtype = np.float)
        from random import randint

        DEM = np.ones((50,60), dtype = np.float)
        for i in range(len(DEM)):
            for j in range(len(DEM[0])):
                DEM[i,j] = DEM[i,j] * randint(1,5)

        # Modify the marsh
        Marsh[0:30,:] = 0
        Marsh[0:20,0:15] = 1
        Marsh[30:,20] = 0
        Marsh[4:26,25:48] = 1
        Marsh[12:14,30:48] = 0
        Marsh[38:40,40:45] = 0
        Marsh[31, 6:10] = 0
        Marsh[8, 14:15] = 0
        Marsh[7:9, 13] = 0
        Marsh[6:11, 5:13] = 0
        Marsh[30:32, 19] = 0
        Marsh[30, 18:20] = 0
        Marsh[3:5, 47:49] = 1
        Marsh[25:27, 26] = 1
        Marsh[28, 29:35] = 1
        Marsh[29, 29] = 1
        Marsh[29, 35] = 1
        Marsh[18, 5:9] = 0
        Marsh[20:27, 41] = 0
        Marsh[2:7, 25:30] = 1
        Marsh[30:35, 0] = 0"""


        """
        Okay, so far we have lots of single lines stored each in a shapefile, and for each line a set of perpendicular transects stored in another shapefile.

        The next steps are:
        - Clean up the code that got us there, i.e. make sensible comments and documentation. Almost Ok
        - Sort the issue of having redundant lines.
        - Find a way to extract profiles. OK
        - Identify sites to work on and get tidal range and wind/wave data and regional bathymetry.
        - Think about point-based forcing indicators
        - Start writing up the article and poster to get a headstart.
        - Send an email to Jaap.

        """
