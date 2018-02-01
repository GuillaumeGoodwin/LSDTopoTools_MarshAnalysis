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
import numpy as np
import cPickle
import timeit
import os

# A very useful package

from LSDMOA_classes import *

from LSDMOA_functions import ENVI_raster_binary_to_2d_array
from LSDMOA_functions import ENVI_raster_binary_from_2d_array
#from LSDMOA_functions import MARSH_ID


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

    
    
    ######
    #TEST ZONE STARTS
    ######


    ######
    #TEST ZONE ENDS
    ######  
    
    
    

    for site in Site:
        print("Loading input data from site: "+site)
        # NB: When loading input data, please make sure the naming convention shown here is respected.
        
        print(" Loading DEM")
        DEM_fname = site+"_DEM.bil"
        DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (Input_dir+DEM_fname)
        print DEM.shape
        
        print " Loading Slopes"
        # Make sure we have the right slope file name
        slope_fname = site+"_slope.bil"
        if not os.path.isfile(Input_dir+slope_fname):
            slope_fname = site+"_SLOPE.bil"   
        Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array (Input_dir+slope_fname)
        print Slope.shape
        
        print " Loading Curvature"
        # Make sure we have the right slope file name
        curv_fname = site+"_curv.bil"
        if not os.path.isfile(Input_dir+curv_fname):
            curv_fname = site+"_CURV.bil"   
        Curv, post_Curv, envidata_Curv =  ENVI_raster_binary_to_2d_array (Input_dir+curv_fname)
        print Curv.shape
        
        print " Loading Marsh Platform"
        # Make sure we have the right slope file name
        marsh_fname = site+"_Marsh.bil"
        Marsh, post_Marsh, envidata_Marsh =  ENVI_raster_binary_to_2d_array (Input_dir+marsh_fname)
        print Marsh.shape
        
  
        #DEM = DEM [700:1300,100:600]
        #Marsh = Marsh [700:1300,100:600]

        
        DEM = DEM [700:800,100:200]
        Marsh = Marsh [700:800,100:200]

        # Here begins the actual MOA
        print "\nRunning the analysis"
        
        #Make a proper object for the marsh
        Marsh_object = 0 * Marsh_platform(Marsh.shape[0], Marsh.shape[1])
        Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = True)    
        Marsh_object.plot_map(Output_dir+'Figures/', '01_Marsh', 'Sous-fifre', Nodata_value)       

        #Marsh_labelled = Marsh_object.label_connected (Nodata_value)
        #Marsh_labelled.plot_map(Output_dir+'Figures/', '02_Marsh_lab', 'Sous-fifre', Nodata_value)       

        #Make the Outline
        Outline_object = 0 * Marsh_outline(Marsh.shape[0], Marsh.shape[1])
        Outline_object.extract_outline_from_array (Marsh_object)
        Outline_object.plot_map(Output_dir+'Figures/', '02_Line', 'Sous-fifre', Nodata_value)

        #Label the connected components
        Outline_labels = Outline_object.label_connected (Nodata_value)
        #Outline_labels.plot_map(Output_dir+'Figures/', '03_Labeled_Line', 'Sous-fifre', Nodata_value)

        Outline_lengths = Outline_labels.calc_outline_length ()
        Outline_lengths.plot_map(Output_dir+'Figures/', '04_Lengths_line', 'Sous-fifre', Nodata_value)
        
        
        
        
        
        
        
        STOP
        
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        DEM_work = np.copy(DEM)
        Search_space, Scarps, Platform = MARSH_ID(DEM, Slope, Nodata_value, opt1, opt2, opt3)
        Platform_work = np.copy(Platform)
        Scarps[Scarps == 0] = Nodata_value

        # Here is where you save your output files for use in a GIS software
        print "Saving marsh features"
        new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, Output_dir+"%s_Search_space.bil" % (site), post_DEM, Search_space)
        new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, Output_dir+"%s_Scarps.bil" % (site), post_DEM, Scarps)
        new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, Output_dir+"%s_Marsh.bil" % (site), post_DEM, Platform)





    # Comment these 2 lines if you don't want to know how long the script run for.
    Stop = timeit.default_timer()
    print 'Analysis runtime was: ', Stop - Start , 's'
