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
from LSDMOA_functions import plot_lines_on_basemap
from LSDMOA_functions import plot_transects_on_basemap
from LSDMOA_functions import Select_few_longest
from LSDMOA_functions import Lines_to_shp
from LSDMOA_functions import Generate_transects
from LSDMOA_functions import Shp_to_lines






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







    def Transect_values (Lines_row, Lines_col, Envidata, Enviarray):
        """
        This function takes a basemap array and a line (or several lines) and returns the values of the cells that include part of the line. These values should be ordered by their distance to the centre of the transect.

        Surely someone has made a similar function before.

        You could also make a function that changes the line into an array selection and operate cell replacement. But this already exists and has been made by the guy who made the transect function (I think).

        Args:
            Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
            Outline_array (2D numpy array): a 2-D array destined to store the outline.
            Outline_value (float): The value to be given to outline cells
            Nodata_value (float): The value for empty cells

        Returns:
            Nothing. It just saves a bunch of shapefiles

        Author: GCHG
        """


        return Lines_dist, Lines_value







    """basemap_array = np.ones((3,3),dtype=np.float)
    from random import randint
    for i in range(len(basemap_array)):
        for j in range(len(basemap_array[0])):
            basemap_array[i,j] = basemap_array[i,j] * randint(1,5)
    print 'This is the array'
    print basemap_array


    Lines_row = [[0,1]]
    Lines_col = [[0,1]]
    print 'This is the line'
    print Lines_row, Lines_col


    print

    quit()"""

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


        #DEM = DEM [850:1050,230:430]
        #Marsh = Marsh [850:1050,230:430]

        DEM = DEM [700:1300,100:600]
        Marsh = Marsh [700:1300,100:600]









        """Marsh = np.ones((50,50), dtype = np.float)
        from random import randint

        DEM = np.ones((50,50), dtype = np.float)
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




        # Here begins the actual MOA
        print "\nRunning the analysis"

        # STEP 1: Make the lines
        #Make a proper object for the marsh
        Marsh_object = Marsh_platform(Marsh.shape[0], Marsh.shape[1])
        Marsh_object = Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = True)
        Marsh_DEM = Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = False)
        Marsh_labels = Marsh_object.label_connected (Nodata_value)


        #Make a proper object for the outline
        Outline_object = 0*Marsh_outline(Marsh.shape[0], Marsh.shape[1])
        Outline_object, Outline_value = Outline_object.complete_outline_from_array (Marsh_object)
        #Tightrope_object, Tightrope_value = Outline_object.tightrope_outline_from_array (Marsh_object)
        Outline_labels = Outline_object.label_connected (Nodata_value)

        #Now get rid of some useless lines
        Outline_simple = Outline_labels.reduce_to_marsh_labels (Marsh_labels, Nodata_value)

        #Now get rid of useless bits of every line
        Outline_trimmed = Outline_simple.trim_to_main_stem (Marsh_object, Nodata_value)

        # Now calculate the lengths
        Outline_length, Lines_row, Lines_col, Lines_dist, Lines_code = Outline_trimmed.calc_outline_length (1)

        # Now select the longest lines
        nLines_row, nLines_col, nLines_dist, nLines_code = Select_few_longest (Lines_row, Lines_col, Lines_dist, Lines_code)

        # Show results of a circular kernel
        #Outlines, Outlines_row, Outlines_col = Outline_labels.vectorise(Nodata_value)
        #Outline_buffer = Outline_object.swath (Outline_value, Nodata_value)

        # Plot the things
        #Marsh_object.plot_map(Output_dir+'Figures/', '00_Marsh_object', 'Sous-fifre', Nodata_value)
        #Marsh_DEM.plot_map(Output_dir+'Figures/', '01_Marsh_DEM', 'Sous-fifre', Nodata_value)
        #Marsh_labels.plot_map(Output_dir+'Figures/', '02_Marsh_Labels', 'Sous-fifre', Nodata_value)
        #Outline_object.plot_map(Output_dir+'Figures/', '03_Outline_object', 'Sous-fifre', Nodata_value)
        #Outline_labels.plot_map(Output_dir+'Figures/', '04_Outline_labels', 'Sous-fifre', Nodata_value)
        #Outline_simple.plot_map(Output_dir+'Figures/', '05_Outline_simple', 'Sous-fifre', Nodata_value)
        #Outline_trimmed.plot_map(Output_dir+'Figures/', '06_Outline_trimmed', 'Sous-fifre', Nodata_value)
        #Outline_length.plot_map(Output_dir+'Figures/', '07_Outline_length', 'Sous-fifre', Nodata_value)
        #plot_lines_on_basemap(Lines_row, Lines_col, Lines_dist, Lines_code, Outline_length, Output_dir+'Figures/', '08_Lines_DivAll', Nodata_value)
        #plot_lines_on_basemap(nLines_row, nLines_col, nLines_dist, nLines_code, Outline_length, Output_dir+'Figures/', '09_Lines_DivAll_simple', Nodata_value)


        # STEP 2: Store the lines into shapefiles
        Lines_to_shp(nLines_row, nLines_col, envidata_DEM, DEM, Output_dir+'Shapefiles/', site)



        # STEP 3: Make the transects with that weird function
        Transects_row = []; Transects_col = []
        # For each label
        for label in range(len(nLines_row)):
            print " \nThe label is:", label+1
            if len(nLines_row[label]) > 0:
                for code in range(len(nLines_row[label])):
                    print "  Saving code number ", code+1
                    Generate_transects(Output_dir+'Shapefiles/'+'%s_%s_%s.shp' % (site,label, code),Output_dir+'Shapefiles/'+'%s_%s_%s_Tr.shp' % (site,label, code), 10, 20)
                    # STEP 4: Put the transect lines into our array reference system
                    Transect_row, Transect_col = Shp_to_lines (Output_dir+"Shapefiles/", "WOR_domain_1_0_0_Tr", envidata_DEM, DEM)

        print Transect_row

        plot_transects_on_basemap(Transect_row, Transect_col, Outline_length, Output_dir+'Figures/', '09_Lines_DivAll_simple', Nodata_value)

        quit()


        """
        Okay, so far we have lots of single lines stored each in a shapefile, and for each line a set of perpendicular transects stored in another shapefile.

        The next steps are:
        - Clean up the code that got us there, i.e. make sensible comments and documentation.
        - Sort the issue of having redundant lines.
        - Find a way to extract profiles.
        - Identify sites to work on and get tidal range and wind/wave data and regional bathymetry.
        - Think about point-based forcing indicators
        - Start writing up the article and poster to get a headstart.
        - Send an email to Jaap.

        """


        #new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Marsh.bil", post_DEM, DEM)

        #ENVI_raster_binary_from_2d_array((envidata, file_out, post, DEM))





        # Right, we have the polylines, but now we need to make them into shapefiles.
        #https://gis.stackexchange.com/questions/85448/creating-polygon-shapefile-from-list-of-x-y-coordinates-using-python

        # FOR SWATHS, USE THiS

        # https://gis.stackexchange.com/questions/50108/elevation-profile-10-km-each-side-of-a-line



        sys.exit()






















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
