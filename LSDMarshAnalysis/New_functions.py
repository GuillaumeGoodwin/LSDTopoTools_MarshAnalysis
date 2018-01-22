# This is a new functions file


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




def True_Outline (Raster, Raster_val, Outline_value, Nodata_value):

    #Find a way to get a 1-pixel-wide outline
    
    Raster[Raster > 0 ] = 1
    Raster_2 = np.zeros(Raster.shape, dtype = np.float)
    Inside = np.where(Raster == 1)
    for i in range(len(Inside[0])):
        x = Inside[0][i]; y = Inside[1][i]
        K = kernel (Raster, 3, x, y)    
        if np.count_nonzero(K) <=  K.size-1 :
            Raster_2[x, y] = Outline_value
    
    Raster_3 = np.copy(Raster_2)       
    Outline = np.where (Raster_2 > 0)
    for i in range(len(Outline[0])):
        x = Outline[0][i]; y = Outline[1][i]
        K_r = kernel (Raster, 3, x, y)
        K = kernel (Raster_2, 3, x, y)
        if np.sum(K_r) == 0:
            print 'By the sea, Mr. Todd'
            Raster_3[x, y] = 0

    return Raster_3



def label_connected (array, Nodata_value):
    array, numfeat = scim.label(array)

    for value in np.arange(1, np.amax(array)+1, 1):
        line = np.where(array == value)

        for n in range(len(line[0])):
            x=line[0][n]; y=line[1][n]
            array_kernel = kernel (array, 3, x, y)
            neighbour = np.where(np.logical_and(array_kernel > 0, array_kernel != value))

            if len(neighbour[0]) > 0 :
                for X in neighbour[0]:
                    for Y in neighbour[1]:
                        array[x+X-1, y+Y-1] = value

    return array




def select_few_longest (array, Nodata_value, num):
    # num is the number of lines you select
    array_2 = np.zeros(array.shape, dtype = np.float)
    
    values = range (np.amin(array[array>0]), np.amax(array), 1)
    line_lengths = []
    for value in values:
        line_lengths.append(len(np.where(array == value)[0]))

    line_lengths = np.asarray(line_lengths)
    Longest = np.where(line_lengths == np.amax(line_lengths))
    print values[Longest[0][0]], line_lengths[Longest[0][0]]
    array_2[array == values[Longest[0][0]]] = values[Longest[0][0]]
    
    if num > 0:
        for i in range(num):
            line_lengths[Longest[0][0]] = 0
            Longest = np.where(line_lengths == np.amax(line_lengths))
            print Longest[0][0]
            array_2[array == values[Longest[0][0]]] = values[Longest[0][0]]

    return array_2




def Measure_polyline_length (array, Nodata_value):
    array_2 = np.zeros(array.shape, dtype = np.float)   
    values = range (int(np.amin(array[array>0])), int(np.amax(array))+1, 1)
    
    for val in range(len(values)):
        print 'Value ', val+1, '/', len(values), ' (',values[val],')'
        
        Lines_x = []
        Lines_y = []
        Lines_dist = []
        
        #Measure all the lines
        num_filled = 0
        num_elements = len (np.where(array == values[val])[0])
        while num_filled < num_elements:
            array_2, Line_x, Line_y, Line_dist = Measure_line_length (array, array_2, values[val], Nodata_value)
            
            Lines_x.append(Line_x)
            Lines_y.append(Line_y)
            Lines_dist.append(Line_dist)
            
            num_filled = len (np.where(np.logical_and(array == values[val], array_2 !=0))[0])
            print 'Number of filled elements = ', num_filled, '/', num_elements
        
        #Stitch the lines
        Bigline_x = []
        Bigline_y = []
        Bigline_dist = []
        
        # Find the line number of the closest starting point to the x-edge (top wall)
        min_x_dist = []
        for i in range(len(Lines_x)):
            min_x_dist.append(Lines_x[i][0]**2)
            min_x_dist.append(Lines_x[i][-1]**2)
        
        A = np.where (min_x_dist == min(min_x_dist)) [0][0]
        print 'Start - line number:', A, 'distance:', np.sqrt(min(min_x_dist)) 
        # Add this line to the Bigline
        if (-1)**A >0:
            print 'Start - coordinates: x = ', Lines_x[A/2][0], ', y = ', Lines_y[A/2][0]
            for i in range(len(Lines_x[A/2])):
                Bigline_x.append(Lines_x[A/2][i])
                Bigline_y.append(Lines_y[A/2][i])
                Bigline_dist.append(Lines_dist[A/2][i])
            Lines_x.remove(Lines_x[A/2])
            Lines_y.remove(Lines_y[A/2])
            Lines_dist.remove(Lines_dist[A/2])
        else:
            # Be careful to reorder the Bigline by inverting the distances
            A = A+1
            print 'Start - coordinates: x = ', Lines_x[A/2][-1], ', y = ', Lines_y[A/2][-1]
            for i in range(len(Lines_x[(A)/2])-1, 0, -1):
                Bigline_x.append(Lines_x[A/2][i])
                Bigline_y.append(Lines_y[A/2][i])
                Bigline_dist.append(Lines_dist[A/2][len(Lines_x[(A)/2])-1-i])
            Lines_x.remove(Lines_x[A/2])
            Lines_y.remove(Lines_y[A/2])
            Lines_dist.remove(Lines_dist[A/2])
        
        print 'End - coordinates: x = ', Bigline_x[-1], ', y = ', Bigline_y[-1]
        print 'End - distance: d = ', Bigline_dist[-1]
        
        
        #for all the next bits:
        while len(Bigline_x) < num_elements:
            print 'Bigline length = ', len(Bigline_x), '/', num_elements
            # Find the closest starting point to the origin
            min_square_dist = []
            for i in range(len(Lines_x)):
                x_prev = Bigline_x[-1]
                y_prev = Bigline_y[-1]
                dist_prev = Bigline_dist[-1]
                            
                head_dist = (Lines_x[i][0]-x_prev)**2+(Lines_y[i][0]-y_prev)**2
                tail_dist = (Lines_x[i][-1]-x_prev)**2+(Lines_y[i][-1]-y_prev)**2
                min_square_dist.append(head_dist)
                min_square_dist.append(tail_dist) 
                
            A = np.where (min_square_dist == min(min_square_dist)) [0][0]
            print 'Next start - line number:', A, 'distance:', np.sqrt(min(min_square_dist)) 
            print 'Next start - distance: d = ', Bigline_dist[-1]
            # Add this line to the Bigline
            if (-1)**A >0:
                print 'Next start - coordinates: x = ', Lines_x[A/2][0], ', y = ', Lines_y[A/2][0]
                
                
                """figure out why they don't save the same thing...s"""
                print len(Lines_x[A/2])
                print len(Lines_y[A/2])
                print len(Lines_dist[A/2])
                
                for i in range(len(Lines_x[A/2])):
                    Bigline_x.append(Lines_x[A/2][i])
                    Bigline_y.append(Lines_y[A/2][i])
                    Bigline_dist.append(Lines_dist[A/2][i] + dist_prev + min(min_square_dist))
                Lines_x.remove(Lines_x[A/2])
                Lines_y.remove(Lines_y[A/2])
                Lines_dist.remove(Lines_dist[A/2])
            else:
                # Be careful to reorder the Bigline by inverting the distances
                A = A+1
                print 'Next start - coordinates: x = ', Lines_x[A/2][-1], ', y = ', Lines_y[A/2][-1]
                for i in range(len(Lines_x[(A)/2])-1, 0, -1):
                    Bigline_x.append(Lines_x[A/2][i])
                    Bigline_y.append(Lines_y[A/2][i])
                    Bigline_dist.append(Lines_dist[A/2][len(Lines_x[(A)/2])-1-i]+dist_prev+min(min_square_dist))
                Lines_x.remove(Lines_x[A/2])
                Lines_y.remove(Lines_y[A/2])
                Lines_dist.remove(Lines_dist[A/2])
            print 'End - coordinates: x = ', Bigline_x[-1], ', y = ', Bigline_y[-1]
            print 'End - distance: d = ', Bigline_dist[-1]
            
        for i in range(len(Bigline_x)):
            array_2[Bigline_x[i], Bigline_y[i]] = Bigline_dist[i]
        
        break
     
     
    
    return array_2









def Measure_line_length (array, array_2, val, Nodata_value):
    num_elements = len (np.where(array == val)[0])
    if num_elements > 0:

        # initialise
        Line_x = []
        Line_y = []
        Line_dist = []
        
        #first point
        Line_x.append(np.where(np.logical_and(array == val, array_2==0))[0][0])
        Line_y.append(np.where(np.logical_and(array == val, array_2==0))[1][0])
        Line_dist.append(0.01)
        array_2 [Line_x[-1], Line_y[-1]] = Line_dist[-1]
        

        #second point
        x = Line_x[0]; y = Line_y[0]
        kernel_size = 3
        K = kernel (array, kernel_size, x, y)
        K_2 = kernel (array_2, kernel_size, x, y)

        if K.size == kernel_size**2:
            New_points = np.where(np.logical_and(K == val, K_2 == 0))
 
            Sea_points = np.where(K == 0)
            if len(New_points[0]) > 0:
                for i in range(len(New_points[0])):
                    for j in range(len(Sea_points[0])):
                        X=New_points[0][i]; Y=New_points[1][i]
                        Xs=Sea_points[0][j]; Ys=Sea_points[1][j]
                        
                        Dist_to_sea = np.sqrt((X-Xs)**2+(Y-Ys)**2)   
                        if Dist_to_sea <= np.sqrt(2):
                            
                            Distance = np.sqrt((X-1)**2+(Y-1)**2)
                            
                            Line_x.append(x+X-1)
                            Line_y.append(y+Y-1)
                            Line_dist.append(Line_dist[-1]+Distance)
                            
                            array_2 [Line_x[-1], Line_y[-1]] = Line_dist[-1]
            
        for t in range(1,1000):           
            x = Line_x[-1]; y = Line_y[-1]
            kernel_size = 3
            K = kernel (array, kernel_size, x, y)
            K_2 = kernel (array_2, kernel_size, x, y)

            if K.size == kernel_size**2:
                New_points = np.where(np.logical_and(K == val, K_2 == 0))

                Sea_points = np.where(K == 0)
                if len(New_points[0]) > 0:
                    for i in range(len(New_points[0])):
                        for j in range(len(Sea_points[0])):
                            X=New_points[0][i]; Y=New_points[1][i]
                            Xs=Sea_points[0][j]; Ys=Sea_points[1][j]

                            Dist_to_sea = np.sqrt((X-Xs)**2+(Y-Ys)**2)   
                            if Dist_to_sea <= np.sqrt(2):
                                Distance = np.sqrt((X-1)**2+(Y-1)**2)

                                Line_x.append(x+X-1)
                                Line_y.append(y+Y-1)
                                Line_dist.append(Line_dist[-1]+Distance)
                                              
                                array_2 [Line_x[-1], Line_y[-1]] = Line_dist[-1]

    return array_2, Line_x, Line_y, Line_dist

















def Kernel_circle (arr, radius, x_centre, y_centre, inclusive = True):
    arr_2 = np.copy(arr)
    
    X_to_0 = x_centre
    X_to_End = len(arr)-x_centre
    Y_to_0 = y_centre
    Y_to_End = len(arr[0])-y_centre

    Lim_left = x_centre - min(np.floor(radius), X_to_0)
    Lim_right = x_centre + min(np.floor(radius)+1, X_to_End)
    Lim_top = y_centre - min(np.floor(radius), Y_to_0)
    Lim_bottom = y_centre + min(np.floor(radius)+1, Y_to_End)
    
    kernel = arr_2[int(Lim_left):int(Lim_right), int(Lim_top):int(Lim_bottom)]

    for x in range(len(kernel)):
        for y in range(len(kernel[0])):
            X = x_centre - radius + x # coordinates in the array
            Y = y_centre - radius + y # coordinates in the array
            Distance = np.sqrt((x_centre-X)**2 + (y_centre-Y)**2)
            
            if Distance <= radius:
                kernel[x,y] = Distance

    return kernel
    

def Kernel_value (arr, radius, x_centre, y_centre, inclusive = True):
    arr_2 = np.copy(arr)
    
    X_to_0 = x_centre
    X_to_End = len(arr)-x_centre
    Y_to_0 = y_centre
    Y_to_End = len(arr[0])-y_centre

    Lim_left = x_centre - min(np.floor(radius), X_to_0)
    Lim_right = x_centre + min(np.floor(radius)+1, X_to_End)
    Lim_top = y_centre - min(np.floor(radius), Y_to_0)
    Lim_bottom = y_centre + min(np.floor(radius)+1, Y_to_End)
    
    kernel = arr_2[int(Lim_left):int(Lim_right), int(Lim_top):int(Lim_bottom)]

    for x in range(len(kernel)):
        for y in range(len(kernel[0])):
            X = x_centre - radius + x # coordinates in the array
            Y = y_centre - radius + y # coordinates in the array
            Distance = np.sqrt((x_centre-X)**2 + (y_centre-Y)**2)
            
            if Distance <= radius:
                kernel[x,y] = arr[X,Y]

    return kernel




def Delete_creeks (marsh, array, window_size, Nodata_value):
    
    array_2 = Measure_polyline_length (array, Nodata_value)
    
    values = range (int(np.amin(array[array>0])), int(np.amax(array))+1, 1)
    for val in range(len(values)):
        print 'Value ', val+1, '/', len(values), ' (',values[val],')'
        
        """lengths = range(int(np.amin(array_2[array>1])), int(np.amax(array_2))+1, 1)
        
        for lgth in range(len(lengths)):
            print 'Length ', lgth+1, '/', len(lengths), ' (',lengths[lgth],')'
            
            Current_position = np.where(np.logical_and(array == values[val], array_2 == lengths[lgth]))

            for i in range(len(Current_position[0])):
                x=Current_position[0][i]; y=Current_position[1][i]
                dist_kernel = Kernel_circle (array, window_size, x, y, inclusive = True)
                object_kernel = Kernel_value (array, window_size, x, y, inclusive = True)
                
                for X in range(len(object_kernel)):
                    for Y in range(len(object_kernel[0])):
                        pass
                        if dist_kernel[X,Y] > 0 and object_kernel [X,Y] != array [x,y]:
                          
                            array_2[X,Y] = 0""""""if dist_kernel[X,Y] > 0 and object_kernel [X,Y] != array [x,y]:
                          
                            array_2[X,Y] = 0"""


        
        
    return array_2




################################################################
def Combined_pdf (Gauges, Metrix_gauges):
    fig=plt.figure("Combined PdF", facecolor='White',figsize=[8,8])

    # Set up the fonts
    matplotlib.rc('xtick', labelsize=9) 
    matplotlib.rc('ytick', labelsize=9)

    ax_raw = plt.subplot2grid((2,1),(0,0), colspan=1, rowspan=1, axisbg='white')
    ax_raw.set_xlabel('Elevation (m)', fontsize = 11)   
    ax_raw.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax_raw.set_ylabel('PDF', fontsize = 11)

    ax_normalised = plt.subplot2grid((2,1),(1,0), colspan=1, rowspan=1, axisbg='white')

    ax_normalised.set_xlabel('Normalised elevation (0-1)', fontsize = 11)   
    ax_normalised.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax_normalised.set_ylabel('PDF', fontsize = 11)

    i = 0
    for gauge in Gauges:

        # Load tidal data
        tidedir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/"
        sourcedir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshAnalysis/"
        Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats (tidedir+"Input/Tide/WOR/WOR_", gauge)
        Metrix_gauges[i,0] = np.mean (Metric2_tide[3])-np.mean (Metric2_tide[0])


        # Load the elevation data
        DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (sourcedir+"Input/%s_DEM.bil" % (gauge), gauge) 

        Platform, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (sourcedir+"Input/%s_Marsh.bil" % (gauge), gauge)
        
        Tidal_flat = np.copy(DEM)
        Tidal_flat [Platform>0] = Nodata_value
        
        Platform [Platform>0] = DEM [Platform>0]
        Platform [Platform <=0] = Nodata_value

        
        
        DEM [DEM < -100] = Nodata_value
        
        bins_DEM, hist_DEM = Distribution (DEM, Nodata_value)
        bins_Pla, hist_Pla = Distribution (Platform, Nodata_value)
        bins_Tid, hist_Tid = Distribution (Tidal_flat, Nodata_value)
        

        Tide_max = np.mean (Metric2_tide[3])
        Tide_neap_max = np.mean (Metric2_tide[2])
        Tide_neap_min = np.mean (Metric2_tide[1])
        Tide_min = np.mean (Metric2_tide[0])
        
        TR = Tide_max - Tide_min
  
        #ax_raw.plot( bins_DEM, hist_DEM, '-', color=plt.cm.jet(i*20), linewidth = 1.5)   
        ax_raw.plot( bins_Pla, hist_Pla, '-', color=plt.cm.jet(i*40), linewidth = 1.5)  
        ax_raw.plot( bins_Tid, hist_Tid, '-', color=plt.cm.jet(0.1*Metrix_gauges[i,0]), linewidth = 0.75)   
        
    
        ax_raw.hlines([0.1 + 0.008*i, 0.1+ 0.008*i], [0], [Tide_neap_min, Tide_neap_max], color=plt.cm.jet(0.1*Metrix_gauges[i,0]), lw=2)
        ax_raw.hlines([0.1 + 0.008*i, 0.1+ 0.008*i], [0], [Tide_min, Tide_max], color=plt.cm.jet(0.1*Metrix_gauges[i,0]), lw=1)

        bins_DEM = np.linspace ((np.amin (bins_DEM)-Tide_min)/TR, (np.amax (bins_DEM)-Tide_min)/TR, len(bins_DEM))
        bins_Pla = np.linspace ((np.amin (bins_Pla)-Tide_min)/TR, (np.amax (bins_Pla)-Tide_min)/TR, len(bins_Pla))
        bins_Tid = np.linspace ((np.amin (bins_Tid)-Tide_min)/TR, (np.amax (bins_Tid)-Tide_min)/TR, len(bins_Tid))

        
        ax_normalised.plot( bins_Pla, hist_Pla, '-', color=plt.cm.jet(i*40), linewidth = 1.5)   
        ax_normalised.plot( bins_Tid, hist_Tid, '-', color=plt.cm.jet(0.1*Metrix_gauges[i,0]), linewidth = 0.75)
        
        ax_raw.set_ylim (0, 0.15)
        ax_raw.set_xlim (-4, 8)
        ax_normalised.set_ylim (0, 0.15)
        
        i = i+1


    plt.savefig(sourcedir+'Output/Marsh_metrics/Combined_pdf_Solway.png')


    
################################################################
def Tidal_prism (Gauges, Metrix_gauges):
    fig=plt.figure("Combined PdF", facecolor='White',figsize=[8,8])

    # Set up the fonts
    matplotlib.rc('xtick', labelsize=9) 
    matplotlib.rc('ytick', labelsize=9)

    ax_raw = plt.subplot2grid((2,1),(0,0), colspan=1, rowspan=1, axisbg='white')
    ax_raw.set_xlabel('Tidal Elevation (m)', fontsize = 11)   
    ax_raw.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax_raw.set_ylabel('Tidal_prism (m3)', fontsize = 11)

    ax_normalised = plt.subplot2grid((2,1),(1,0), colspan=1, rowspan=1, axisbg='white')

    ax_normalised.set_xlabel('Normalised elevation (0-1)', fontsize = 11)   
    ax_normalised.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax_normalised.set_ylabel('Tidal Prism', fontsize = 11)

    i = 0
    for gauge in Gauges:
        # Load tidal data
        Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats ("Input/Tide/%s/%s_" % (gauge,gauge), gauge)
        Tide_max = np.mean (Metric2_tide[3])
        Tide_neap_max = np.mean (Metric2_tide[2])
        Tide_neap_min = np.mean (Metric2_tide[1])
        Tide_min = np.mean (Metric2_tide[0])
        TR = Tide_max - Tide_min
        
        # Load the elevation data
        DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_DEM_clip.bil" % (gauge,gauge), gauge) 
        Platform, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_O1_-2.0_O2_0.85_O3_8.0_Marsh_nofilter.bil" % (gauge,gauge), gauge)
        
        #Tidal_flat = np.copy(DEM)
        #Tidal_flat [Platform>0] = Nodata_value
        
        Platform [Platform>0] = DEM [Platform>0]
        Platform [Platform <=0] = Nodata_value
        

        tidal_elevation = np.linspace (Tide_min, Tide_max, 50)
        tidal_prism = np.zeros (len(tidal_elevation), dtype = np.float)
        tidal_prism_Pla = np.zeros (len(tidal_elevation), dtype = np.float)
        
        for z in range(1, len(tidal_elevation)):
            Flooded = tidal_elevation[z] - DEM; Flooded [Flooded <= 0] = 0; Flooded [DEM < Nodata_value +1] = 0
            Flooded_Pla = tidal_elevation[z] - Platform; Flooded_Pla [Flooded_Pla <= 0] = 0; Flooded_Pla [Platform < Nodata_value+1] = 0
            
            tidal_prism[z] = (np.sum(Flooded) + tidal_prism[z-1]) / (DEM.size) # *  tidal_elevation[z]
            tidal_prism_Pla[z] = (np.sum(Flooded_Pla) + tidal_prism_Pla[z-1]) / (Platform.size)
  
        
    
        ax_raw.plot( tidal_elevation, tidal_prism, '-', color=plt.cm.jet(0.1*Metrix_gauges[i,0]), linewidth = 0.75)   
        ax_raw.plot( tidal_elevation, tidal_prism_Pla, '-', color=plt.cm.jet(0.1*Metrix_gauges[i,0]), linewidth = 1.5)   
   
        #ax_raw.hlines([0.1 + 0.008*i, 0.1+ 0.008*i], [0], [Tide_neap_min, Tide_neap_max], color=plt.cm.jet(0.1*Metrix_gauges[i,0]), lw=2)
        #ax_raw.hlines([0.1 + 0.008*i, 0.1+ 0.008*i], [0], [Tide_min, Tide_max], color=plt.cm.jet(0.1*Metrix_gauges[i,0]), lw=1)
        
        
        
        #ax_normalised.plot( bins_Pla, hist_Pla, '-', color=plt.cm.jet(0.1*Metrix_gauges[i,0]), linewidth = 1.5)   
        #ax_normalised.plot( bins_Tid, hist_Tid, '-', color=plt.cm.jet(0.1*Metrix_gauges[i,0]), linewidth = 0.75)
        
        ax_raw.set_ylim (ymin = 0)
        #ax_normalised.set_ylim (0, 0.15)
        
        i = i+1


    plt.savefig('Output/Marsh_metrics/Tidal_prism.png')


    
    
    
def Scarp_stats (srcdir, dstdir, Gauges, tiddir):
    
    i = 0
    Metrix_gauges = np.zeros((len(Gauges),2), dtype = np.float)    
    for gauge in Gauges:
        # Load tidal data
        Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats (tiddir+"WOR/WOR_", gauge)
        Metrix_gauges[i,0] = np.mean (Metric2_tide[3])-np.mean (Metric2_tide[0])

        # Load topography and marshes
        DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array (srcdir+"%s_DEM.bil" % (gauge), gauge) 
        Platform, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array (srcdir+"%s_Marsh.bil" % (gauge), gauge)
        
        
        #DEM = DEM[1000:1200, 50:250]
        #Platform = Platform[1000:1200, 50:250]
        
        Platform [Platform > 0] = 1
  
        Platform_outline = Outline (Platform, 2, Nodata_value)
        
        Platform_outline[Platform_outline == 1] = Nodata_value
        
        Edge = np.where (Platform_outline == 2)

        for n in range(len(Edge[0])):
            x=Edge[0][n]; y=Edge[1][n]
            Platform_outline_kernel = kernel (Platform_outline, 15, x, y)
            DEM_kernel = kernel (DEM, 15, x, y)
            
            for X in range(Platform_outline_kernel.shape[0]):
                for Y in range(Platform_outline_kernel.shape[1]):
                    
                    if np.amin(DEM_kernel) != Nodata_value:
                    
                        if x+X-1 > 0 and x+X-1 < Platform_outline.shape[0]:
                            if y+Y-1 > 0 and y+Y-1 < Platform_outline.shape[1]:

                                #Platform_outline[x+X-1, y+Y-1] = DEM[x+X-1, y+Y-1]
                                Platform_outline[x+X-1, y+Y-1] = np.amax(DEM_kernel) - np.amin(DEM_kernel[DEM_kernel>-10])
    
    
        Platform_outline[Platform_outline<=0] = Nodata_value
        Platform_outline = np.ma.masked_where(Platform_outline == Nodata_value, Platform_outline)
        
        bins_Pla, hist_Pla = Distribution (Platform, Nodata_value)
        
        pickle.dump(bins_Pla, open(sourcedir+'Output/Marsh_metrics/bins_'+gauge+'.pkl', "wb" ) )
        pickle.dump(hist_Pla, open(sourcedir+'Output/Marsh_metrics/hist_'+gauge+'.pkl', "wb" ) )
        
        
        
        #Make the figure
        fig=plt.figure("Combined PdF", facecolor='White',figsize=[8,8])

        # Set up the fonts
        matplotlib.rc('xtick', labelsize=9) 
        matplotlib.rc('ytick', labelsize=9)

        ax_raw = plt.subplot2grid((2,1),(0,0), colspan=1, rowspan=1, axisbg='white')
        #ax_raw.set_xlabel('Elevation (m)', fontsize = 11)   
        #ax_raw.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        #ax_raw.set_ylabel('PDF', fontsize = 11)

        ax_normalised = plt.subplot2grid((2,1),(1,0), colspan=1, rowspan=1, axisbg='white')

        ax_normalised.set_xlabel('elevation (m)', fontsize = 11)   
        #ax_normalised.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax_normalised.set_ylabel('PDF', fontsize = 11)

        Map = ax_raw.imshow(Platform_outline, interpolation='None', cmap=plt.cm.gist_earth)#, vmin=0, vmax=3)  
        
        ax_normalised.plot(bins_Pla, hist_Pla)
    
        plt.savefig(sourcedir + 'Output/Marsh_metrics/Scarp_stats_%g.png' % (i))

    
        i += 1


