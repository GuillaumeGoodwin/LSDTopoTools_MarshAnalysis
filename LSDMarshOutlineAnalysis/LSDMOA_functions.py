"""
    This package processes topographic data in order to extract marsh platforms
    Guillaume C.H. Goodwin
    Released unedr GPL3
"""

# Load useful Python packages
import os
import sys
import numpy as np

import cPickle

import gdal



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
#from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import scipy.ndimage as scim

import cPickle as pickle


from osgeo import ogr
from shapely.geometry import LineString

import fiona

import pandas as bb


import numpy as np
from fiona import collection
from shapely.geometry import LineString, MultiLineString



import LSDMOA_classes as cl


##########################################################################################################
##########################################################################################################
def ENVI_raster_binary_to_2d_array(file_name):
    """
    This function transforms a raster into a numpy array.

    Args:
        file_name (ENVI raster): the raster you want to work on.
        gauge (string): a name for your file

    Returns:
        image_array (2-D numpy array): the array corresponding to the raster you loaded
        pixelWidth (geotransform, inDs) (float): the size of the pixel corresponding to an element in the output array.

    Source: http://chris35wills.github.io/python-gdal-raster-io/
    """

    driver = gdal.GetDriverByName('ENVI')

    driver.Register()

    inDs = gdal.Open(file_name)

    if inDs is None:
        print "Couldn't open this file: " + file_name
        print "Perhaps you need an ENVI .hdr file? "
        sys.exit("Try again!")
    else:
        print "%s opened successfully" %file_name

        print '~~~~~~~~~~~~~~'
        print 'Get image size'
        print '~~~~~~~~~~~~~~'
        cols = inDs.RasterXSize
        rows = inDs.RasterYSize
        bands = inDs.RasterCount

        print "columns: %i" %cols
        print "rows: %i" %rows
        print "bands: %i" %bands

        print '~~~~~~~~~~~~~~'
        print 'Get georeference information'
        print '~~~~~~~~~~~~~~'
        geotransform = inDs.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]

        print "origin x: %i" %originX
        print "origin y: %i" %originY
        print "width: %2.2f" %pixelWidth
        print "height: %2.2f" %pixelHeight

        # Set pixel offset.....
        print '~~~~~~~~~~~~~~'
        print 'Convert image to 2D array'
        print '~~~~~~~~~~~~~~'
        band = inDs.GetRasterBand(1)
        image_array = band.ReadAsArray(0, 0, cols, rows)
        image_array_name = file_name
        print type(image_array)
        print image_array.shape

        return image_array, pixelWidth, (geotransform, inDs)



##########################################################################################################
##########################################################################################################
def ENVI_raster_binary_from_2d_array(envidata, file_out, post, image_array):
    """
    This function transforms a numpy array into a raster.

    Args:
        envidata: the geospatial data needed to create your raster
        file_out (string): the name of the output file
        post: coordinates for the goegraphical transformation
        image_array (2-D numpy array): the input raster

    Returns:
        new_geotransform
        new_projection: the projection in which the raster
        file_out (ENVI raster): the raster you wanted

    Source: http://chris35wills.github.io/python-gdal-raster-io/
    """
    driver = gdal.GetDriverByName('ENVI')

    original_geotransform, inDs = envidata

    rows, cols = image_array.shape
    bands = 1

    # Creates a new raster data source
    outDs = driver.Create(file_out, cols, rows, bands, gdal.GDT_Float32)

    # Write metadata
    originX = original_geotransform[0]
    originY = original_geotransform[3]

    outDs.SetGeoTransform([originX, post, 0.0, originY, 0.0, -post])
    outDs.SetProjection(inDs.GetProjection())

    #Write raster datasets
    outBand = outDs.GetRasterBand(1)
    outBand.WriteArray(image_array)

    new_geotransform = outDs.GetGeoTransform()
    new_projection = outDs.GetProjection()

    print "Output binary saved: ", file_out

    return new_geotransform,new_projection,file_out



##########################################################################################################
##########################################################################################################
def kernel (array, kernel_size, row_centre, col_centre):
    """
    This function defines a square kernel within an array (array), centred on (x_centre, y_centre). The is of a width of kernel_size.
    Args:
        array (2D numpy array): a 2-D array.
        kernel_size (float): the width of the square defining the size of the kernel. kernel_size MUST be an ODD number to account for the central element.
        x_centre (int): The index of the element in the 1st dimension.
        y_centre (int): The index of the element in the 2nd dimension.

    Returns:
        kernel (2D numpy array): The kernel of selected elements.
        kernel_row (2D numpy array): The kernel of row indices of the selected elements.
        kernel_col (2D numpy array): The kernel of column indices of the selected elements.

    Author: GCHG
    """

    if (-1)**kernel_size < 0:
        row_to_0 = row_centre
        row_to_End = array.shape[0]-row_centre
        col_to_0 = col_centre
        col_to_End = array.shape[1]-col_centre

        width = np.floor(kernel_size/2)

        Lim_top = row_centre - min(width, row_to_0)
        Lim_bottom = row_centre + min(width+1, row_to_End)
        Lim_left = col_centre - min(width, col_to_0)
        Lim_right = col_centre + min(width+1, col_to_End)

        kernel = array [int(Lim_top):int(Lim_bottom), int(Lim_left):int(Lim_right)]
        kernel_row = np.arange(int(Lim_top),int(Lim_bottom))
        kernel_col = np.arange(int(Lim_left),int(Lim_right))

    else:
        print
        print " ... WARNING: you need to choose an odd kernel size, buddy"
        print
        pass

    return kernel, kernel_row, kernel_col


##########################################################################################################
##########################################################################################################
def Surface_outline (Outline_array, Surface_array, Outline_value, Nodata_value):
    """
    This function calculates the inner outline of a surface within an array (Surface_array) where the element to outline has a value of 1, and stores that outline in a second array (Outline_array) under the value Outline_value.
    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Outline_array (2D numpy array): a 2-D array populated with the outline cells.

    Author: GCHG
    """

    Surface_array[Surface_array > 0. ] = 1
    Inside = np.where(Surface_array == 1)

    for i in range(len(Inside[1])):
        x = Inside[0][i]; y = Inside[1][i]
        K, Kx, Ky = kernel (Surface_array, 3, x, y)

        if np.count_nonzero(K) <=  K.size-1 :
            Outline_array[x, y] = Outline_value

    """ THIS PART IS CURRENTLY USELESS.
    for i in range(len(Outline[0])):
        x = Outline[0][i]; y = Outline[1][i]
        K_r, Kx_r, Ky_r = kernel (Surface_array, 3, x, y)
        K, Kx, Ky = kernel (Twin, 3, x, y)
        if np.sum(K_r) == 0:
            Outline_array[x, y] = 0"""

    return Outline_array, Outline_value


##########################################################################################################
##########################################################################################################
def stitch_segments (segments, M_code, reset_length = False, select_longest = True):
    """This function stitches tiny segments into something that makes sense"""

    #Initiate the objects
    segments[np.isnan(segments)] = 0
    L_code = 1
    # This is the list of outlines
    M_outlines = cl.Polyline()

    # loop until all the segments are stitched
    while np.amax(segments) != 0:
        # Setup the pandas outline
        M_outline = cl.Line()

        #Find the first non-null element in the array
        nonzero_x = np.where(segments[:,0]!=0)[0][0]
        nonzero_y = np.where(segments[:,1]!=0)[0][0]

        # choose a this segment to be the first element in the Pandas
        M_outline.set_first_point(M_code, L_code, segments[nonzero_y,1], segments[nonzero_x,0], 1)
        M_outline = M_outline.add_element(1,M_code, L_code, segments[nonzero_y+1,1], segments[nonzero_x+1,0])

        # Now stitch the connected segments together in both directions. Because we're 2D, bruv!
        for d in [-1,0]:
            for i in range(3,len(segments),3):
                # Find a segment that starts at the end of the previous segment
                last_point = M_outline['rowcol'].iloc[d]
                find_next = np.where(np.logical_and(segments[:,1] == last_point.row(), segments[:,0] == last_point.col()))[0]

                # Add the segment. Invert it if need be.
                for j in range(len(find_next)):
                    if find_next[j]/3. == find_next[j]//3:
                        M_outline = M_outline.add_element(-d*len(M_outline['rowcol']),M_code, L_code, segments[find_next[j]+1,1], segments[find_next[j]+1,0])
                        segments[find_next[j]:find_next[j]+3,:] = 0
                    elif (find_next[j]-1)/3. == (find_next[j]-1)//3:
                        M_outline = M_outline.add_element(-d*len(M_outline['rowcol']),M_code, L_code, segments[find_next[j]-1,1], segments[find_next[j]-1,0])
                        segments[find_next[j]-1:find_next[j]+2,:] = 0

        #Recalculate the distances to make it nicer. Takes a loooong time
        if reset_length is True:
            M_outline = M_outline.recalc_length(1)

        #add this outline to the list of outlines
        M_outlines.append(M_outline)

        # update the L_code counter
        L_code+=1

    #This part is to group everything in one Pandas instead of a silly list
    Final_outlines = M_outlines[0]
    for i in range(1,len(M_outlines)):
        Final_outlines = bb.concat([Final_outlines,M_outlines[i]])

    if select_longest is True:
        Lengths = []
        value_range = [min(Final_outlines['L_code']),max(Final_outlines['L_code'])]
        for L in value_range:
            Lengths.append(len(Final_outlines.loc[Final_outlines['L_code'] == L]))
        Longest = max(Lengths)
        for L in value_range:
            if len(Final_outlines.loc[Final_outlines['L_code'] == L]) < 0.5 * Longest:
                Final_outlines = Final_outlines[Final_outlines.L_code != L]

    return Final_outlines


##########################################################################################################
##########################################################################################################
def Pandas_outline (Surface_array, M_code, scale):
    """
    This magic function from the internet (https://stackoverflow.com/questions/24539296/outline-a-region-in-a-graph) takes an array of positives v. negatives and outlines the border between the two. It gives you a nice Pandas Dataframe at the end.

    Args:

    Returns:

    Author: GCHG
    """

    image = Surface_array
    maskimg = np.zeros(Surface_array.shape, dtype='int')
    maskimg[image == M_code] = 3

    x0 = 0; x1 = Surface_array.shape[1]
    y0 = 0; y1 = Surface_array.shape[0]

    # our image with the numbers 1-3 is in array maskimg
    # create a boolean image map which has trues only where maskimg[x,y] == 3
    mapimg = (maskimg == 3)

    # a vertical line segment is needed, when the pixels next to each other horizontally
    #   belong to diffferent groups (one is part of the mask, the other isn't)
    # after this ver_seg has two arrays, one for row coordinates, the other for column coordinates
    ver_seg = np.where(mapimg[:,1:] != mapimg[:,:-1])
    # the same is repeated for horizontal segments
    hor_seg = np.where(mapimg[1:,:] != mapimg[:-1,:])

    # if we have a horizontal segment at 7,2, it means that it must be drawn between pixels
    #   (2,7) and (2,8), i.e. from (2,8)..(3,8)
    # in order to draw a discountinuous line, we add Nones in between segments
    l = []
    for p in zip(*hor_seg):
        l.append((p[1], p[0]+1))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan,np.nan))

    # and the same for vertical segments
    for p in zip(*ver_seg):
        l.append((p[1]+1, p[0]))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan, np.nan))

    # now we transform the list into a numpy array of Nx2 shape
    segments = np.array(l)
    # now we need to know something about the image which is shown
    #   at this point let's assume it has extents (x0, y0)..(x1,y1) on the axis
    #   drawn with origin='lower'
    # with this information we can rescale our points
    segments[:,0] = (x0 + (x1-x0) * segments[:,0] / mapimg.shape[1]) - scale/2.
    segments[:,1] = (y0 + (y1-y0) * segments[:,1] / mapimg.shape[0]) - scale/2.

    #This is me now ^^
    M_outline = stitch_segments(segments,M_code)

    return M_outline








##########################################################################################################
##########################################################################################################
def Line_to_shp (line, Envidata, Enviarray, save_dir, file_name):
    """
    This function takes all these masses of riggly lines and turns each of them into a shapefile. This is going to take a lot of space, so we should probably think of deleting the shapefiles when we're done with them.

    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Nothing. It just saves a bunch of shapefiles

    Author: GCHG
    """
    X_origin = Envidata[0][0]; X_cell_width = Envidata[0][1]
    Y_origin = Envidata[0][3]; Y_cell_width = Envidata[0][5]

    Coordinates = []
    for i in range(len(line['rowcol'])):
        x = line['rowcol'].iloc[i].col()
        y = line['rowcol'].iloc[i].row()
        coord_pair = (X_cell_width * x + X_origin, Y_cell_width * y + Y_origin)
        Coordinates.append(coord_pair)
    # Here's an example Shapely geometry
    poly = LineString(Coordinates)
    # Now convert it to a shapefile with OGR
    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource('%s/%s.shp' % (save_dir,file_name))
    #ds = driver.CreateDataSource(save_dir+str(label)+'_'+str(code)+'.shp')
    layer = ds.CreateLayer('', None, ogr.wkbLineString)
    # Add one attribute
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    defn = layer.GetLayerDefn()
    ## If there are multiple geometries, put the "for" loop here
    # Create a new feature (attribute and geometry)
    feat = ogr.Feature(defn)
    feat.SetField('id', 123)
    # Make a geometry, from Shapely object
    geom = ogr.CreateGeometryFromWkb(poly.wkb)
    feat.SetGeometry(geom)
    layer.CreateFeature(feat)
    feat = geom = None  # destroy these
    # Save and close everything
    ds = layer = feat = geom = None

    print "file saved:" '%s/%s.shp' % (save_dir,file_name)


##########################################################################################################
##########################################################################################################
def Shp_to_transects (file_dir, file_name, M_code, L_code, Envidata, Enviarray):
    """

    Args:

    Returns:

    Author: GCHG
    """
    # Load the file
    shape = fiona.open(file_dir+file_name)

    # This s the environment of the raster
    X_origin = Envidata[0][0]; X_cell_width = Envidata[0][1]
    Y_origin = Envidata[0][3]; Y_cell_width = Envidata[0][5]

    # initialise storage list
    Transects = cl.Transect()
    Transects.set_first_point(0, 0, 0, 0)

    for i in range(len(shape)):
        attr = shape.next()
        geometry = attr['geometry']
        coordinates = geometry['coordinates']

        Start_row = (coordinates[0][1]-Y_origin)/Y_cell_width
        End_row = (coordinates[1][1]-Y_origin)/Y_cell_width
        Start_col = (coordinates[0][0]-X_origin)/X_cell_width
        End_col = (coordinates[1][0]-X_origin)/X_cell_width

        Transects = Transects.add_element(len(Transects['rowcol']), M_code, L_code, Start_row, Start_col)
        Transects = Transects.add_element(len(Transects['rowcol']), M_code, L_code, End_row, End_col)

    Transects = Transects.iloc[1:]
    Transects = Transects.assign_transect_code()

    return Transects


##########################################################################################################
##########################################################################################################
def Make_transects (in_shp, out_shp, spc, sect_len):
    #-------------------------------------------------------------------------------
    # Name:        perp_lines.py
    # Purpose:     Generates multiple profile lines perpendicular to an input line
    #
    # Author:      JamesS
    #
    # Created:     13/02/2013
    #-------------------------------------------------------------------------------
    """ Takes a shapefile containing a single line as input. Generates lines
        perpendicular to the original with the specified length and spacing and
        writes them to a new shapefile.

        The data should be in a projected co-ordinate system.
    """

    # ##############################################################################
    # User input

    # Input shapefile. Must be a single, simple line, in projected co-ordinates


    # The shapefile to which the perpendicular lines will be written


    # Profile spacing. The distance at which to space the perpendicular profiles
    # In the same units as the original shapefile (e.g. metres)


    # Length of cross-sections to calculate either side of central line
    # i.e. the total length will be twice the value entered here.
    # In the same co-ordinates as the original shapefile

    #https://gis.stackexchange.com/questions/50108/elevation-profile-10-km-each-side-of-a-line

    # ##############################################################################

    # Open the shapefile and get the data
    source = collection(in_shp, "r")
    data = source.next()['geometry']
    line = LineString(data['coordinates'])

    # Define a schema for the output features. Add a new field called 'Dist'
    # to uniquely identify each profile
    schema = source.schema.copy()
    schema['properties']['Dist'] = 'float'

    # Open a new sink for the output features, using the same format driver
    # and coordinate reference system as the source.
    sink = collection(out_shp, "w", driver=source.driver, schema=schema, crs=source.crs)

    # Calculate the number of profiles to generate
    n_prof = int(line.length/spc)

    # Start iterating along the line
    for prof in range(1, n_prof+1):
        # Get the start, mid and end points for this segment
        seg_st = line.interpolate((prof-1)*spc)
        seg_mid = line.interpolate((prof-0.5)*spc)
        seg_end = line.interpolate(prof*spc)

        # Get a displacement vector for this segment
        vec = np.array([[seg_end.x - seg_st.x,], [seg_end.y - seg_st.y,]])

        # Rotate the vector 90 deg clockwise and 90 deg counter clockwise
        rot_anti = np.array([[0, -1], [1, 0]])
        rot_clock = np.array([[0, 1], [-1, 0]])
        vec_anti = np.dot(rot_anti, vec)
        vec_clock = np.dot(rot_clock, vec)

        # Normalise the perpendicular vectors
        len_anti = ((vec_anti**2).sum())**0.5
        vec_anti = vec_anti/len_anti
        len_clock = ((vec_clock**2).sum())**0.5
        vec_clock = vec_clock/len_clock

        # Scale them up to the profile length
        vec_anti = vec_anti*sect_len
        vec_clock = vec_clock*sect_len

        # Calculate displacements from midpoint
        prof_st = (seg_mid.x + float(vec_anti[0]), seg_mid.y + float(vec_anti[1]))
        prof_end = (seg_mid.x + float(vec_clock[0]), seg_mid.y + float(vec_clock[1]))

        # Write to output
        rec = {'geometry':{'type':'LineString', 'coordinates':(prof_st, prof_end)},
               'properties':{'Dist':(prof-0.5)*spc, u'id':0}}

        sink.write(rec)

    # Tidy up
    source.close()
    sink.close()




##########################################################################################################
##########################################################################################################
def plot_transects_on_basemap (Lines_row, Lines_col, Basemap, save_dir, fig_name, Nodata_value):
    """
    This function calculates the length of several connected lines of connected cells in an Outline_array and returns lists of x,y coordinates (Line_x, Line_y) and length values (Line_dist) of the cells along the line.


    There's a labeled array and a storage array

    of a surface within an array (Surface_array) where the element to outline has a value of 1, and stores that outline in a second array (Outline_array) under the value Outline_value.
    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Outline_array (2D numpy array): a 2-D array populated with the outline cells.

    Author: GCHG
    """
    from matplotlib.lines import Line2D
    print 'Plotting a coup'
    Nodata_value = -1000


    twin  = Basemap.copy()

    fig_height = min(np.floor(twin.shape[1])/5, 50)
    fig_width = min(np.floor(twin.shape[1])/5, 50)

    fig=plt.figure('Some Title', facecolor='White',figsize=[fig_height,fig_width])
    ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
    ax1.tick_params(axis='x', colors='black')
    ax1.tick_params(axis='y', colors='black')

    Vmin = min(np.amin(twin[twin!=Nodata_value])*0.95, np.amin(twin[twin!=Nodata_value])*1.05)
    Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)

    Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
    ax2 = fig.add_axes([0.1, 0.98, 0.85, 0.02])
    scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
    cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)

    for i in range(len(Lines_row)):
        #for j in range(len(Lines_row[i])):

        Scatt = ax1.scatter(Lines_col[i][0], Lines_row[i][0], marker  = '+', color = 'b')
        Scatt2 = ax1.scatter(Lines_col[i][-1], Lines_row[i][-1], marker = 'x', color = 'r')

        Line = Line2D(Lines_col[i], Lines_row[i], color = plt.cm.gist_earth(50), alpha = 0.5)
        ax1.add_line(Line)

    #ax1.set_xlim(0, len(Basemap)-1)
    #ax1.set_ylim(len(Basemap[0])-1, 0)

    plt.savefig(save_dir+fig_name+'.png')



##########################################################################################################
##########################################################################################################
def Line_length (Length_array, Labels_array, Code_array, Label_value, Scale, Code):
    """
    This function calculates the length of a single line of connected cells in an Outline_array and returns lists of x,y coordinates (Line_x, Line_y) and length values (Line_dist) of the cells along the line. It also returns a coded array (Code_array) for the line in case you need to stitch several lines. It uses the scale of the raster to calculate distance in a tricky way.

    Args:
        Length_array (2D numpy array): a 2-D array containing the Length array to fill.
        Labels_array (2D numpy array): a 2-D array containing the Labels of the lines.
        Code_array (2D numpy array): a 2-D array containing the Code array to fill.
        Label_value (float): A float telling us which label to look at.
        Scale (float): A float telling us the scale of the original raster.
        Code (float): A float telling us what code to attribute to the measured line
        Short (Bool): Decide whether you want the line to be the longest or shortest possible
    Returns:
        Length_array (2D numpy array): a 2-D array containing the filled Length array.
        Line_row (list): a list of row coordinates for the filled cells.
        Line_col (list): a list of column coordinates for the filled cells.
        Line_dist (list): a list of distances in the filled cells.
        Code_array (2D numpy array): a 2-D array containing the code value in the filled cells.
        Filled_elements:

    Author: GCHG
    """
    kernel_size = 3

    Elements = np.where(Length_array == Label_value); num_elements = len(Elements[0])
    # initialise line
    This_line = cl.Line()
    This_dist = []

    if num_elements > 0:
        #Add the the first point's coordinates to the lists (first = closest to the origin)
        First_point = cl.Point(Elements[0][0],Elements[1][0])
        This_line.append(First_point); This_dist.append(0.001*Scale)
        #Change the lengths array to include the new distance
        Length_array [First_point.row(), First_point.col()] = This_dist[0]
        #Update the Code array
        Code_array [First_point.row(), First_point.col()] = Code

        for t in range(1, num_elements):
            #Find the coordinates of the last point added
            row = This_line.end_point().row(); col = This_line.end_point().col()

            #Make a kernel around that point
            K_lab, Kr_lab, Kc_lab = kernel (Labels_array, kernel_size, row, col)
            K_len, Kr_len, Kc_len = kernel (Length_array, kernel_size, row, col)
            #This is where the Kernel contains un-distanced but labelled line values
            Next_points = np.where(np.logical_and(K_lab == Label_value, K_len == Label_value))

            if len(Next_points[0]) > 0:
                for i in range(len(Next_points[0])):
                    #These are the coordinats of the candidates in the big array
                    ROW = Kr_lab[Next_points[0]]; COL = Kc_lab[Next_points[1]]
                    #This is the distance between the candidates and the centre
                    Dist_to_previous = np.sqrt((ROW-row)**2+(COL-col)**2)
                    #We select the closest one to make sure we get all the points
                    Selected_next = np.where (Dist_to_previous == np.amin(Dist_to_previous[Dist_to_previous>0]))
                    #Now we add these values to our list
                    Next_point = cl.Point(ROW[Selected_next[0][0]],COL[Selected_next[0][0]])
                    This_line.append(Next_point)
                    # We also increase the distance
                    Next_dist = This_dist[-1] + Dist_to_previous[Selected_next[0][0]]; This_dist.append(Next_dist)
                    #Change the lengths array to include the new distance
                    Length_array [Next_point.row(), Next_point.col()] = Next_dist
                    #Record the code for that line
                    Code_array [Next_point.row(), Next_point.col()] = Code
                    break

        Filled_elements = len(This_dist)

    This_line.append(This_dist)

    return This_line, Length_array, Code_array, Filled_elements




##########################################################################################################
##########################################################################################################
def Measure_all_lines (Labels_array, Label_value, Scale):
    """
    This function calculates the length of several lines an Outline_array and returns lists of x,y coordinates (Line_x, Line_y) and length values (Line_dist) of the cells along the line.


    There's a labeled array and a storage array

    of a surface within an array (Surface_array) where the element to outline has a value of 1, and stores that outline in a second array (Outline_array) under the value Outline_value.
    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Outline_array (2D numpy array): a 2-D array populated with the outline cells.

    Author: GCHG
    """
    Length_array = Labels_array.copy(); Code_array = 0 * Labels_array.copy()

    Elements = np.where(Labels_array == Label_value); num_elements = len (Elements[0])
    print ' This label has ', num_elements, ' elements'

    #Initialise the lists to store point object in the lines
    This_polyline = cl.Polyline()
    print type(This_polyline)
    # Initialise the number of elements filled and the counter for the line code
    Filled_total = 0; code = 1
    #THIS IS THE BIT WHERE YOU CALCULATE THE LENGTH OF EACH INDIVIDUAL LINE
    while Filled_total < num_elements:
        # Calculate the length of each squiggly line (which may stop abruptly)
        This_line, Length_array, Code_array, Filled_elements = Line_length (Length_array, Labels_array, Code_array, Label_value, Scale, code)
        This_line.append(code)
        #add to the lists of all the line coordinates, lengths, and code
        This_polyline.append(This_line)
        # Keep track of filled elements
        Filled_total = Filled_total + Filled_elements
        print '  Number of filled elements: ', Filled_total, '/', num_elements, ' (code:' , code ,')'
        #Update the line code
        code += 1
    # Reduce to 0 the elements that were not filled
    Length_array[Length_array == Label_value] = 0

    return This_polyline, Length_array, Code_array



##########################################################################################################
##########################################################################################################
def Stitch_diverging_starts (Length_array, Labels_array, Label_value, This_polyline, Code_array, Scale):
    """
    This function stitches lines that start at neighbouring points

    There's a labeled array and a storage array

    of a surface within an array (Surface_array) where the element to outline has a value of 1, and stores that outline in a second array (Outline_array) under the value Outline_value.
    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Outline_array (2D numpy array): a 2-D array populated with the outline cells.

    Author: GCHG
    """

    print "########################################"
    print "We are looking for divergent line starts"

    kernel_size = 3
    Indices_to_clean = []

    # For each line in the polyline
    for i in range(len(This_polyline)):
        This_line = This_polyline[i]
        #See how long each line is is.
        Counter = 0
        for j in range(len(This_line)):
            if type(This_line[j]) is cl.Point:
                Counter += 1

        #If this is not a 1-pixel-long line:
        if Counter>1:
            # Find the starting point
            Startpoint_row = This_line.start_point()[0]; Startpoint_col = This_line.start_point()[1]
            # Make kernels around the starting point
            sK_lab, sKr_lab, sKc_lab = kernel (Labels_array, kernel_size, Startpoint_row, Startpoint_col)
            sK_len, sKr_len, sKc_len = kernel (Length_array, kernel_size, Startpoint_row, Startpoint_col)
            sK_cod, sKr_cod, sKc_cod = kernel (Code_array, kernel_size, Startpoint_row, Startpoint_col)
            #Also find its end point
            Endpoint_row = This_line.end_point()[0]; Endpoint_col = This_line.end_point()[1]
            # Make kernels around the end point
            eK_lab, eKr_lab, eKc_lab = kernel (Labels_array, kernel_size, Endpoint_row, Endpoint_col)
            eK_len, eKr_len, eKc_len = kernel (Length_array, kernel_size, Endpoint_row, Endpoint_col)
            eK_cod, eKr_cod, eKc_cod = kernel (Code_array, kernel_size, Endpoint_row, Endpoint_col)
            # Also retrieve the line's code
            This_line_code = Code_array[Startpoint_row, Startpoint_col]
            This_line_code_index = This_line_code - 1

            #Tell something to the confused user
            print '\nLine starts (', Startpoint_row, Startpoint_col, ').', 'Code is :' , int(This_line_code)
            print '....Line ends (', Endpoint_row, Endpoint_col, '). Length is :' , Length_array[Endpoint_row, Endpoint_col]

            # Record where the length kernel is non-null AND has a different line code.
            # This is a connecting line
            Touches_start = np.where(np.logical_and(sK_len > 0, sK_cod.astype(int) != int(This_line_code)))
            Touches_end = np.where(np.logical_and(eK_len > 0, eK_cod.astype(int) != int(This_line_code)))

            # If it's connected at the start
            if len(sK_len[Touches_start[:][0]]) > 0 and np.count_nonzero(sK_cod[Touches_start[:][0]]) > 0 :
                print '.....This line connects to another at its starting point.'
                # Find the smallest distance in that subset and where it is in the array
                Smallest_dist = np.amin (sK_len[Touches_start]); Smallest_dist_index = np.where(sK_len == Smallest_dist)
                # If the shortest distance is the beginning of a line AND you have two of those, apply Rule X
                if Smallest_dist <= 0.001*Scale and len(Smallest_dist_index[0]) == 2 :
                    print '......Tis a case of diverging starts'
                    # Find the coordinates and distance of these points.
                    First_line_row = sKr_len[Smallest_dist_index[0][0]]; First_line_col = sKc_len[Smallest_dist_index[1][0]]
                    Second_line_row = sKr_len[Smallest_dist_index[0][1]]; Second_line_col = sKc_len[Smallest_dist_index[1][1]]
                    #Find out how far aprt they are
                    Distance = np.sqrt((First_line_row-Second_line_row)**2 + (First_line_col-Second_line_col)**2)
                    # Find the code associated to those two points
                    First_line_code = Code_array[First_line_row, First_line_col]; Second_line_code = Code_array[Second_line_row, Second_line_col]
                    print ".......connected code are ", First_line_code, Second_line_code
                    print sK_cod
                    # Select the two lines that go with these codes
                    First_line_index = int(First_line_code)-1; Second_line_index = int(Second_line_code)-1

                    Second_Counter = 0
                    for j in range(len(This_polyline[Second_line_index])):
                        if type(This_polyline[Second_line_index][j]) is cl.Point:
                            Second_Counter += 1
                    # Reverse the second line (arbitrary)
                    This_polyline[Second_line_index][:Second_Counter] = reversed(This_polyline[Second_line_index][:Second_Counter])


                    # What if we could stop there and let the second function do the job?

                    """#Make sure it's not one of those short lines
                    if len(This_polyline[Second_line_index][:Second_Counter]) > 1:
                        # Add the extra distance values to the first line
                        This_polyline[First_line_index][-2] = (np.asarray(This_polyline[First_line_index][-2]) + This_polyline[Second_line_index][-2][-1] + Distance).tolist()
                        # Stitch the lines
                        This_polyline[Second_line_index][-2] = This_polyline[Second_line_index][-2] + This_polyline[First_line_index][-2]
                        This_polyline[Second_line_index][-1] = First_line_code
                        This_polyline[Second_line_index][:-2] = This_polyline[Second_line_index][:-2] + This_polyline[First_line_index][:-2]

                        # Update the length array
                        Line_row = []; Line_col = []
                        for i in range(len(This_polyline[Second_line_index][:Second_Counter+Counter])):
                            Line_row.append(This_polyline[Second_line_index][i][0])
                            Line_col.append(This_polyline[Second_line_index][i][1])

                        Length_array[Line_row, Line_col] = This_polyline[Second_line_index][-2]
                        Code_array[Line_row, Line_col] = First_line_code
                        Indices_to_clean.append(First_line_index)"""

    """#Now rid of the lines that have become redundant
    for i in reversed(Indices_to_clean):
        del This_polyline[i]

    # And reset the codes for the remaining lines.
    for i in range(len(This_polyline)):
        This_polyline[i][-1] = i+1
        Line_row = []; Line_col = []
        for j in range(len(This_polyline[i][:-2])):
            Line_row.append(This_polyline[i][j][0])
            Line_col.append(This_polyline[i][j][1])

        Code_array[Line_row, Line_col] = i+1"""


    return This_polyline, Length_array, Code_array



##########################################################################################################
##########################################################################################################
def Graft_diverging_branch (Length_array, Labels_array, Label_value, This_polyline, Code_array, Scale):
    """
    This function stitches a line that has its starting point in contact with another line that is already on its way.

    There's a labeled array and a storage array

    of a surface within an array (Surface_array) where the element to outline has a value of 1, and stores that outline in a second array (Outline_array) under the value Outline_value.
    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Outline_array (2D numpy array): a 2-D array populated with the outline cells.

    Author: GCHG
    """

    print "########################################"
    print "We are looking for diverging branches"

    kernel_size = 3

    # For each line in the polyline
    for i in range(len(This_polyline)):
        This_line = This_polyline[i]
        #See how long each line is is.
        Counter = 0
        for j in range(len(This_line)):
            if type(This_line[j]) is cl.Point:
                Counter += 1

        #If this is not a 1-pixel-long line:
        if Counter>1:
            # Find the starting point
            Startpoint_row = This_line.start_point()[0]; Startpoint_col = This_line.start_point()[1]
            # Make kernels around the starting point
            sK_lab, sKr_lab, sKc_lab = kernel (Labels_array, kernel_size, Startpoint_row, Startpoint_col)
            sK_len, sKr_len, sKc_len = kernel (Length_array, kernel_size, Startpoint_row, Startpoint_col)
            sK_cod, sKr_cod, sKc_cod = kernel (Code_array, kernel_size, Startpoint_row, Startpoint_col)
            #Also find its end point
            Endpoint_row = This_line.end_point()[0]; Endpoint_col = This_line.end_point()[1]
            # Make kernels around the end point
            eK_lab, eKr_lab, eKc_lab = kernel (Labels_array, kernel_size, Endpoint_row, Endpoint_col)
            eK_len, eKr_len, eKc_len = kernel (Length_array, kernel_size, Endpoint_row, Endpoint_col)
            eK_cod, eKr_cod, eKc_cod = kernel (Code_array, kernel_size, Endpoint_row, Endpoint_col)
            # Also retrieve the line's code
            This_line_code = Code_array[Startpoint_row, Startpoint_col]
            This_line_code_index = This_line_code-1

            #Tell something to the confused user
            print '\nLine starts (', Startpoint_row, Startpoint_col, ').', 'Code is :' , int(This_line_code)
            print '....Line ends (', Endpoint_row, Endpoint_col, '). Length is :' , Length_array[Endpoint_row, Endpoint_col]

            # Record where the length kernel is non-null AND has a different line code.
            # This is a connecting line
            Touches_start = np.where(np.logical_and(sK_len > 0, sK_cod.astype(int) != int(This_line_code)))
            Touches_end = np.where(np.logical_and(eK_len > 0, eK_cod.astype(int) != int(This_line_code)))

            # If it's connected at the start
            if len(sK_len[Touches_start[:][0]]) > 0 and np.count_nonzero(sK_cod[Touches_start[:][0]]) > 0 :
                print '.....This line connects to another at its starting point.'
                # Find the (ironically named) biggest distance in that subset and where it is in the array
                Smallest_dist = np.amax (sK_len[Touches_start]); Smallest_dist_index = np.where(sK_len == Smallest_dist)

                # If the shortest distance is NOT  beginning of a line:
                if Smallest_dist > 0.001*Scale:
                    print '......this line was born in the wake of another'
                    # Find the coordinates and distance of the connecting point on the other line
                    This_line_row = Startpoint_row; This_line_col = Startpoint_col
                    Other_line_row = sKr_len[Smallest_dist_index[0][0]]; Other_line_col = sKc_len[Smallest_dist_index[1][0]]

                    #Find out how far apart they are
                    Distance = np.sqrt((Other_line_row-This_line_row)**2 + (Other_line_col-This_line_col)**2)
                    # Find the codes
                    This_line_code = Code_array[This_line_row, This_line_col]
                    Other_line_code = Code_array[Other_line_row, Other_line_col]

                    # Now find the line indices
                    This_line_index = int(This_line_code)-1; Other_line_index = int(Other_line_code)-1

                    # Find the right distance along the other line

                    print sK_len
                    print sK_lab
                    print sK_cod

                    print This_line_index, Other_line_index
                    print This_polyline[This_line_index][-2]
                    print This_polyline[Other_line_index][-2]
                    print Smallest_dist

                    Distance_index = np.where(This_polyline[Other_line_index][-2] == Smallest_dist)[0][0]

                    #Prepare to append the other line
                    Dist_to_append = This_polyline[Other_line_index][-2][0:Distance_index+1]
                    Coord_to_append = This_polyline[Other_line_index][0:Distance_index+1]

                    #If there's something to append
                    if len(Dist_to_append) > 1:
                        #Add the correct distances
                        This_polyline[This_line_index][-2] = (np.asarray(This_polyline[This_line_index][-2]) + Dist_to_append[-1] + Distance).tolist()

                        # Stitch the lines
                        This_polyline[This_line_index][-2] = Dist_to_append + This_polyline[This_line_index][-2]

                        This_polyline[This_line_index][:-2] = Coord_to_append + This_polyline[This_line_index][:-2]

                        #Update the arrays
                        Line_row = []; Line_col = []
                        for i in range(len(This_polyline[This_line_index][:-2])):
                            Line_row.append(This_polyline[This_line_index][i][0])
                            Line_col.append(This_polyline[This_line_index][i][1])

                        Length_array[Line_row, Line_col] = This_polyline[This_line_index][-2]
                        #Code_array[Lines_row[This_line_index], Lines_col[This_line_index]] = This_line_code



    return This_polyline, Length_array



##########################################################################################################
##########################################################################################################
def Graft_converging_branch (Length_array, Labels_array, Label_value, Lines_row, Lines_col, Lines_dist, Lines_code, Code_array, Scale):
    """
    This function stitches a line that has its starting point in contact with another line that is already on its way.

    There's a labeled array and a storage array

    of a surface within an array (Surface_array) where the element to outline has a value of 1, and stores that outline in a second array (Outline_array) under the value Outline_value.
    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Outline_array (2D numpy array): a 2-D array populated with the outline cells.

    Author: GCHG
    """

    print "\n########################################"
    print "We are looking for converging branches"

    kernel_size = 3

    # For each line
    for i in range(len(Lines_row)):
        #If this is not a 1-pixel-long line:
        if len(Lines_row[i])>1:
            # Find the starting point
            Startpoint_row = Lines_row[i][0]; Startpoint_col = Lines_col[i][0]
            # Make kernels around the starting point
            sK_lab, sKr_lab, sKc_lab = kernel (Labels_array, kernel_size, Startpoint_row, Startpoint_col)
            sK_len, sKr_len, sKc_len = kernel (Length_array, kernel_size, Startpoint_row, Startpoint_col)
            sK_cod, sKr_cod, sKc_cod = kernel (Code_array, kernel_size, Startpoint_row, Startpoint_col)
            #Also find its end point
            Endpoint_row = Lines_row[i][-1]; Endpoint_col = Lines_col[i][-1]
            # Make kernels around the end point
            eK_lab, eKr_lab, eKc_lab = kernel (Labels_array, kernel_size, Endpoint_row, Endpoint_col)
            eK_len, eKr_len, eKc_len = kernel (Length_array, kernel_size, Endpoint_row, Endpoint_col)
            eK_cod, eKr_cod, eKc_cod = kernel (Code_array, kernel_size, Endpoint_row, Endpoint_col)
            # Also retrieve the line's code
            This_line_code = Code_array[Startpoint_row, Startpoint_col]; This_line_code_index = np.where(Lines_code == This_line_code)[0][0]

            #Tell something to the confused user
            print '\nLine starts (', Startpoint_row, Startpoint_col, ').', 'Code is :' , int(This_line_code)
            print '....Line ends (', Endpoint_row, Endpoint_col, '). Length is :' , Length_array[Endpoint_row, Endpoint_col]

            # Record where the length kernel is non-null AND has a different line code.
            # This is a connecting line
            Touches_start = np.where(np.logical_and(sK_len > 0, sK_cod.astype(int) != int(This_line_code)))
            Touches_end = np.where(np.logical_and(eK_len > 0, eK_cod.astype(int) != int(This_line_code)))


            #If it's connected at the end
            if len(eK_len[Touches_end[:][0]]) > 0 and np.count_nonzero(eK_cod[Touches_end[:][0]]) > 0 :
                print '.....This line connects to another at its ending point.'
                # Find the smallest distance in that subset and where it is in the array
                Smallest_dist = np.amin (eK_len[Touches_end]); Smallest_dist_index = np.where(eK_len == Smallest_dist)
                # Find the connecting line code
                #Smallest_dist_code = np.amin(eK_cod[np.where(eK_len == Smallest_dist)])

                #If this distance is shorter than the one of the line we are looking at
                if Smallest_dist < Length_array[Endpoint_row, Endpoint_col]:
                    print '......The other line is shorter'
                    #The first line is the one we are looking at, and the second is the other one
                    Long_line_row = Endpoint_row; Long_line_col = Endpoint_col
                    Short_line_row = eKr_len[Smallest_dist_index[0][0]]; Short_line_col = eKc_len[Smallest_dist_index[1][0]]
                else:
                    #Well it's the other way around, innit?
                    print '......The other line is longer'
                    Short_line_row = Endpoint_row; Short_line_col = Endpoint_col
                    Long_line_row = eKr_len[Smallest_dist_index[0][0]]; Long_line_col = eKc_len[Smallest_dist_index[1][0]]


                #Tis the distance spearating these lines
                Distance = np.sqrt((Short_line_row-Long_line_row)**2 + (Short_line_col-Long_line_col)**2)
                # Find the code associated to those two points
                Short_line_code = Code_array[Short_line_row, Short_line_col]; Long_line_code = Code_array[Long_line_row, Long_line_col]
                # Select the two lines that go with these codes
                Short_line_index = np.where(Lines_code == Short_line_code); Long_line_index = np.where(Lines_code == Long_line_code)

                #Tell something to the confused user
                print 'Here are the starting points:'
                print Lines_row[Short_line_index[0]][0], Lines_col[Short_line_index[0]][0] , '/', Lines_row[Long_line_index[0]][0], Lines_col[Long_line_index[0]][0]

                #Check that these lines don't have a common origin
                if Lines_row[Short_line_index[0]][0] != Lines_row[Long_line_index[0]][0] or  Lines_col[Short_line_index[0]][0] != Lines_col[Long_line_index[0]][0] :

                    print "WHOOOEEEE"

                    #sys.exit()

                    # Reverse the short line (arbitrary)
                    Lines_row[Short_line_index[0]] = list(reversed(np.asarray(Lines_row[Short_line_index[0]])))
                    Lines_col[Short_line_index[0]] = list(reversed(np.asarray(Lines_col[Short_line_index[0]])))

                    #Make sure it's not one of those short lines
                    if len(Lines_dist[Short_line_index[0]]) > 1:
                    # Add the extra distance values to the first line
                        Lines_dist[Long_line_index[0]] = (np.asarray(Lines_dist[Long_line_index[0]]) + Lines_dist[Short_line_index[0]][-1] + Distance).tolist()

                        #Where to stop stitching
                        Limit_dist_index = np.where(np.floor(Lines_dist[Long_line_index[0]]) == np.floor(Smallest_dist))
                        #print Smallest_dist
                        #print Lines_dist[Long_line_index[0]]
                        #print Limit_dist_index

                        if len(Limit_dist_index[0]) > 0:

                            # Stitch the lines
                            Lines_dist[Short_line_index[0]] =  Lines_dist[Long_line_index[0]][0:Limit_dist_index[0]] + Lines_dist[Short_line_index[0]]
                            Lines_row[Short_line_index[0]] = Lines_row[Long_line_index[0]][0:Limit_dist_index[0]] + Lines_row[Short_line_index[0]]
                            Lines_col[Short_line_index[0]] = Lines_col[Long_line_index[0]][0:Limit_dist_index[0]] + Lines_col[Short_line_index[0]]
                            # Update the length array
                            Length_array[Lines_row[Short_line_index[0]], Lines_col[Short_line_index[0]]] = Lines_dist[Short_line_index[0]]


    return Length_array, Lines_row, Lines_col, Lines_dist, Lines_code



##########################################################################################################
##########################################################################################################
def Select_few_longest (Polylines):
    """
    This function selects the longest lines in a set of polylines. The longest line within a set sharing the same origin is selected.
    The structure for each list is: 1st index = label; 2nd index = code; 3rd index = position on line

    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Outline_array (2D numpy array): a 2-D array populated with the outline cells.

    Author: GCHG
    """

    print "\n########################################"
    print "Let's find the longest branch"


    # Instead of deleting, store the selected things in a new array. it's easier.
    # Initialise the new lines
    new_Polylines = cl.Polyline()
    # For each label
    for label in range(len(Polylines)):
        print " \nThe label is:", label+1
        # Initialise a list of indices to keep
        # Its first version is a dummy.
        This_new_polyline = cl.Polyline()
        This_new_polyline.append([cl.Point(0,0),[0], 0])

        if len(Polylines[label]) > 0:
            for code in range(len(Polylines[label])):
                # This is the origin of this line in the label
                Origin = Polylines[label][code][0]
                #print "   The code is:", Polylines[label][code][-1]
                #print "    Line origin:", Origin,"; Line length:", Polylines[label][code][-2][-1]

                # Make a list of all the origin points in the new storage polyline for this label
                Origins_list = []
                for i in range(len(This_new_polyline)):
                    Origins_list.append(This_new_polyline[i][0])
                #Find origins in the Origins list that match that of the examined line
                Matching_origin = [item for item in Origins_list if item == Origin]

                # If our origin matches that storage polyline
                if len(Matching_origin) > 0 :
                    print "   The code is:", Polylines[label][code][-1]
                    print "    Line origin:", Origin,"; Line length:", Polylines[label][code][-2][-1]
                    print "     matching origins:", Matching_origin
                    Matching_origin_id = Origins_list.index(Matching_origin[0])
                    # Compare lengths between the line examined and the stored one:
                    if Polylines[label][code][-2][-1] > This_new_polyline[Matching_origin_id][-2][-1]:
                        # The longest line takes the spot in the list
                        This_new_polyline[Matching_origin_id] = Polylines[label][code]

                # if these origins are not in the set:
                else:
                    # We're not adding really small lines
                    if len(Polylines[label][code][-2]) > 1:
                        print "     tis a new line"
                        # Add the line to the selected set.
                        This_new_polyline.append(Polylines[label][code])

        # Get rid of the dummy
        del This_new_polyline[0]

        # Store the result in the new big polylines list
        new_Polylines.append(This_new_polyline)



    return new_Polylines




##########################################################################################################
##########################################################################################################
def porthole (array, radius, row_centre, col_centre):
    """
    This function defines a circular kernel within an array (array), centred on (x_centre, y_centre). The is of a width of kernel_size.
    Args:
        array (2D numpy array): a 2-D array.
        radius (float): the raidus of the circle defining the.
        x_centre (int): The index of the element in the 1st dimension.
        y_centre (int): The index of the element in the 2nd dimension.

    Returns:
        kernel (2D numpy array): The kernel of selected elements.
        kernel_row (2D numpy array): The kernel of row indices of the selected elements.
        kernel_col (2D numpy array): The kernel of column indices of the selected elements.

    Author: GCHG
    """

    porthole_row = []; porthole_col = []; porthole_dist = []

    row_to_0 = row_centre
    row_to_End = len(array[:,0])-row_centre
    col_to_0 = col_centre
    col_to_End = len(array)-col_centre

    Lim_top = row_centre - min(np.floor(radius), row_to_0)
    Lim_bottom = row_centre + min(np.floor(radius)+1, row_to_End)
    Lim_left = col_centre - min(np.floor(radius), col_to_0)
    Lim_right = col_centre + min(np.floor(radius)+1, col_to_End)


    kernel = array [int(Lim_top):int(Lim_bottom), int(Lim_left):int(Lim_right)]

    for row in range(len(kernel)):
        for col in range(len(kernel[0])):
            ROW = row_centre - radius + row # coordinates in the array
            COL = col_centre - radius + col # coordinates in the array
            Distance = np.sqrt((row_centre-ROW)**2 + (col_centre-COL)**2)

            if Distance < radius:
                porthole_dist.append(Distance)
                porthole_row.append(ROW)
                porthole_col.append(COL)

    return porthole_row, porthole_col, porthole_dist










##########################################################################################################
##########################################################################################################
def Billhook (Outline_labels_array, Marsh_array, Label_value):
    """
    This function takes

    of a surface within an array (Surface_array) where the element to outline has a value of 1, and stores that outline in a second array (Outline_array) under the value Outline_value.
    Args:
        Surface_array (2D numpy array): a 2-D array containing the surface to outline with the value 1. Undesirable elements have the value 0 or Nodata_value.
        Outline_array (2D numpy array): a 2-D array destined to store the outline.
        Outline_value (float): The value to be given to outline cells
        Nodata_value (float): The value for empty cells

    Returns:
        Outline_array (2D numpy array): a 2-D array populated with the outline cells.

    Author: GCHG
    """
    kernel_size = 3

    #Select a line to work on
    Active_line = np.where (Outline_labels_array == Label_value)
    num_elements = len (Active_line[0])

    #If it's not empty
    if num_elements > 0:
        #For each cell in the line
        for t in range(num_elements):
            row = Active_line [0][t]
            col = Active_line [1][t]
            #Make a kernel of labels and of marsh values
            K_lab, Kr_lab, Kc_lab = kernel (Outline_labels_array, kernel_size, row, col)
            K_Mar, Kr_Mar, Kc_Mar = kernel (Marsh_array, kernel_size, row, col)

            #If the Kernel has 5 or more labelled cells:
            if np.sum(K_lab) >= 5 * Label_value and np.sum(K_Mar) >= 8:
                #Check that they don't have a full empty row or col:
                if np.sum(K_lab[0]) > Label_value and np.sum(K_lab[2]) > Label_value and np.sum(K_lab[:,0]) > Label_value and np.sum(K_lab[:,2]) > Label_value:
                    Outline_labels_array[row, col] = 0

    return Outline_labels_array














    """for j in range(len(Previous_points[0])):
        X=New_points[0][i]; Y=New_points[1][i]
        Xp=Previous_points[0][j]; Yp=Previous_points[1][j]

        Dist_to_previous = np.sqrt((X-Xp)**2+(Y-Yp)**2)



        if Dist_to_previous <= np.sqrt(2):

            Distance = np.sqrt((X-1)**2+(Y-1)**2)

            Line_x.append(x+X-1)
            Line_y.append(y+Y-1)
            Line_dist.append(Line_dist[-1]+Distance)

            array_2 [Line_x[-1], Line_y[-1]] = Line_dist[-1]


    print 'so far so good'

    #Add the other points' coordinates to the lists
    for t in range(1, num_elements):
        print 'iterating', t
        # Latest coordinate stored
        x = Line_x[-1]; y = Line_y[-1]
        K_lab = kernel (Labels_array, kernel_size, x, y)
        K_len = kernel (Length_array, kernel_size, x, y)

        print K_lab

        #if K_lab.size == kernel_size**2:


        New_points = np.where(K_lab == Label_value)
        Old_points = np.where(K_lab == Line_dist[-2])

        if len(New_points[0]) > 0:

            print 'New points can be added'

            for i in range(len(New_points[0])):
                for j in range(len(Sea_points[0])):
                    X=New_points[0][i]; Y=New_points[1][i]
                    #Xs=Sea_points[0][j]; Ys=Sea_points[1][j]

                    #Dist_to_sea = np.sqrt((X-Xs)**2+(Y-Ys)**2)
                    #if Dist_to_sea <= np.sqrt(2):

                        #print 'Distance is right'

                    Distance = np.sqrt((X-1)**2+(Y-1)**2)

                    Line_x.append(x+X-1)
                    Line_y.append(y+Y-1)
                    Line_dist.append(Line_dist[-1]+Distance)

                    Length_array [Line_x[-1], Line_y[-1]] = Line_dist[-1]

                    print 'appended'"""





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

























#---------------------------------------------------------------
# This function opens tidal stats data
def Open_tide_stats (file_location, gauge):
    print 'Loading tidal statistics of' + file_location + ' %s' % (gauge)


    with open (file_location + "Metric1.pkl", 'rb') as input_file:
        Metric1 = cPickle.load(input_file)
    with open (file_location + "Metric2.pkl", 'rb') as input_file:
        Metric2 = cPickle.load(input_file)
    with open (file_location + "Metric3.pkl", 'rb') as input_file:
        Metric3 = cPickle.load(input_file)
    with open (file_location + "Subsample.pkl", 'rb') as input_file:
        Subsample = cPickle.load(input_file)

    return Metric1, Metric2, Metric3, Subsample








#-----------------------------------------------------------------------------------------------------------
def Distribution(Data2D, Nodata_value):
    """
    This simple function takes a 2-D array (Data2D) and makes a probability distribution of its values. It is set to ignore elements with a specific value (Nodata_value).

    Args:
        Data2D (2D numpy array): the 2D array you want a distribution for
        Nodata_value (float): The value for ignored elements

    Returns:
        bins [1D numpy array]: the value bins
        hist [1D numpy array]: the probability associated to the bins

    Author: GCHG
    """

    Data1D = Data2D.ravel()

    Max_distribution = max(Data1D)
    if len(Data1D[Data1D>Nodata_value]) == 0:
        Min_distribution = -1
    else:
        Min_distribution = min(Data1D[Data1D>Nodata_value])

    bin_size = (Max_distribution - Min_distribution) / 100

    X_values = np.arange(Min_distribution, Max_distribution, bin_size)


    hist, bins = np.histogram (Data1D, X_values, density=True)
    hist=hist/sum(hist)
    bins=bins[:-1]


    return bins,hist



#-----------------------------------------------------------------------------------------------------------


def Outline (Raster, Outline_value, Nodata_value):
    """
    This simple function takes a 2-D array (Raster) and attributes a specific value (Outline value) to elements at the limit of a bloc of elements with identical values. Effectively, it draws an outline around a group of elements with the same value. It is set to ignore elements with a specific value (Nodata_value).

    Args:
        Raster (2D numpy array): the 2-D array
        Outline_value (float): The value associated to the outline. Be smart and select a different value from those already in your 2-D array.
        Nodata_value (float): The value for ignored elements

    Returns:
        Raster (2D numpy array): the 2-D array, with the outlines given their own value.

    Author: GCHG
    """
    P1 = np.where(Raster[:,1:] != Raster[:,:-1])
    Raster[P1] = Outline_value

    P2 = np.where(Raster[1:,:] != Raster[:-1,:])
    Raster[P2] = Outline_value

    for i in range(len(Raster)):
        for j in range(len(Raster[0,:])):
            if Raster[i,j] == Outline_value:
                K = kernel (Raster, 3, i, j)
                if np.mean(K) < 0:
                    Raster[i,j] = Nodata_value

    return Raster




#-----------------------------------------------------------------------------------------------------
def define_search_space (DEM, Slope, Nodata_value, opt):
    """
   This function defines a search space (Search_space) within a 2-D array, based on the combined values of 2 2-D arrays (DEM and Slope) of the same dimensions. It defines the threshold for the selection of the search space according to a threshold value (opt). It is set to ignore elements with a specific value (Nodata_value).
    Args:
        DEM (2D numpy array): a 2-D array (here a DEM) used as a first condition for the definition of the search space
        Slope (2D numpy array): a 2-D array (here a DEM) used as a second condition for the definition of the search space
        Nodata_value (float): The value for ignored elements
        opt (float): the value of the threshold for the selection of the search space

    Returns:
        Search_space (2D numpy array): The resulting search space array. Search_space has a value of 0 for non-selected elements and 1 for selected elements.
        Crossover (2D numpy array): The array resulting of the multiplication of relative slope and relative relief.
        bins (1D array): the value bins for the Crossover array
        hist (1D array): the value hist for the Crossover array
        Inflecion_point(float): the value of the threshold for the search space selection.

    Author: GCHG
    """


    print 'Choosing a holiday destination ...'
    Height = len(DEM); Width = len(DEM[0,:])
    Search_space = np.zeros((Height,Width), dtype=np.float)

    # We calculate the relative relief of the DEM to have values of elevation between 0 and 1
    Relief = DEM-np.amin(DEM[DEM > Nodata_value])
    Rel_relief = Relief/np.amax(Relief)
    Rel_relief[DEM == Nodata_value] = Nodata_value

    # We then do the same thing for slope
    Rel_slope = Slope/np.amax(Slope)
    Rel_slope[Slope == Nodata_value] = Nodata_value

    # We then multiply these new relative relief and slope arrays and biologically name them "Crossover"
    Crossover = Rel_relief * Rel_slope
    Crossover[DEM == Nodata_value] = Nodata_value

    # We make a curve of the frequency of values in this Crossover
    # That curve should look like a decreasing exponential function
    data = Crossover.ravel(); data = data[data>0]
    step = (max(data) - min(data)) / 100
    value = np.arange(min(data), max(data), step)
    hist, bins = np.histogram (data, value, density=True)
    hist=hist/sum(hist); bins=bins[:-1]

    # We now find the slope of that curve
    hist_der = np.zeros(len(hist), dtype = np.float)
    for j in range(1, len(hist), 1):
        hist_der[j] = (hist[j]-hist[j-1])/step

    # If the slope gets above the -1 threshold, now that we have hit the closest point to the origin.
    # We call it the inflexion point even though it's not really an inflexion point.
    for j in range(1, len(hist)-1, 1):
        if hist_der[j] < opt and hist_der[j+1] >= opt:
            Inflexion_point = bins[j]

    # Points within the search space should have a Crossover value above the inflexion point
    Search = np.where(Crossover > Inflexion_point)
    Search_space[Search] = 1

    # We get rid of the borders of the DEM because otherwise it will be difficult to work with the smaller slope array
    Search_space[0,:] = 0; Search_space[Height-1,:] = 0; Search_space[:,0] = 0; Search_space[:,Width-1] = 0

    # And update the search locations for the shaved edges
    Search = np.where(Search_space == 1)

    # If this happens, your landscape is weird
    if np.amax(Search_space) == 0:
        print
        print " ... Your search space is empty! Are you sure there's a marsh platform here?"
        print
        STOP

    return Search_space, Crossover, bins, hist, Inflexion_point


#-----------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
def peak_flag (Slope, Search_space, Order):
    """
    This function is the first stage of a routing process used to identify lines of maximum slopes.
    This function identifies multiple local maxima in an array (Slope), within a predefined search space (Search_space). The identified maxima are given a value of Order.

    Args:
        Slope (2D numpy array): the input 2-D array, here issued from a slope raster.
        Search_space (2D numpy array): the search space array in which to look for local maxima.
        Order (int): the value given to the local maxima points.

    Returns:
        Peaks (2D numpy array): a 2-D array where the local maxima have a value of Order and other elements are null.
        Slope_copy (2D numpy array): a copy of the input array where the value of the selected local maxima has been set to 0.

    Author: GCHG
    """

    print 'Finding local slope maxima ...'
    Slope_copy = np.copy(Slope) # the copy of the initial data array
    Search = np.where(Search_space == 1) # the searched locations
    Peaks = np.zeros((len(Slope),len(Slope[0,:])),dtype = np.float)

    for i in range(len(Search[0])):
        x=Search[0][i]; y=Search[1][i] # coordinates of the kernel's centre
        Kernel_slope = kernel (Slope, 3, x, y)
        Kernel_search = kernel(Search_space, 3, x, y)

        # if the centre of the kernel is its maximum and is not an isolated point
        if Kernel_slope[1,1] == np.amax(Kernel_slope) and np.amax(Kernel_search[Kernel_search<=Kernel_search[1,1]] > 0):
            Peaks[x,y] = Order # The kernel centre becomes a local peak
            Slope_copy[x,y] = 0 # The slope of the modified data array drops to 0

    return Peaks, Slope_copy



#-----------------------------------------------------------------------------------------------------
def initiate_ridge (Slope, Search_space, Peaks, Order):
    """
    This function is the second stage of a routing process used to identify lines of maximum slopes.
    This function identifies multiple duplets of elements in an array (Slope), within a predefined search space (Search_space) and within the neighbourhood of the local maxima identified in a second input array (Peaks). The identified elements are given a value of Order. To make this function work, the input array Slope should be the output array Slope_copy of the function peak_flag.

    Args:
        Slope (2D numpy array): the input 2-D array, here issued from a slope raster where the local maximal values have been replaced by 0.
        Search_space (2D numpy array): the search space array.
        Peaks (2D numpy array): A 2-D array containing elements with a value of 1. These elements have the same indices as the elements with a value of 0 in Slope.
        Order (int): the value given to the identified elements. it should be superior by 1 to the value of Order in the function peak_flag.

    Returns:
        Ridges (2D numpy array): a 2-D array where the identified elements have a value of Order. This array is modified from the Peaks array and therefore also contains elements of a value equal to the Order in the function peak_flag.
        Slope_copy (2D numpy array): a copy of the input array where the value of the selected elements has been set to 0.

    Author: GCHG
    """

    print ' ... Starting ridges ...'
    Slope_copy = np.copy(Slope) # the copy of the initial data array
    Search = np.where(Search_space == 1) # the searched locations
    Search_peaks = np.where(Peaks == Order-1) # the searched locations where the peaks are
    Ridges = np.copy(Peaks)

    # Define Kernels
    for i in range(len(Search_peaks[0])):
        x=Search_peaks[0][i]; y=Search_peaks[1][i] # coordinates of the kernel's centre
        Kernel_slope = kernel (Slope, 3, x, y)
        Kernel_slope_copy = kernel (Slope_copy, 3, x, y)
        Kernel_ridges = kernel (Ridges, 3, x, y)
        Kernel_search = kernel (Search_space, 3, x, y)

        # 1/ If there are no other peaks, we have two ridge starters
        if np.count_nonzero(Kernel_ridges) == 1:
            Ridge_starter1 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
            X1=Ridge_starter1[0][0]; Y1=Ridge_starter1[1][0]

            # if it is within the initial search space
            if Search_space[x+X1-1, y+Y1-1] != 0:
                Ridges[x+X1-1, y+Y1-1] = Order
                Slope_copy[x+X1-1, y+Y1-1] = 0

                # Look for a second ridge starter
                Ridge_starter2 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
                X2=Ridge_starter2[0][0]; Y2=Ridge_starter2[1][0]
                Distance = np.sqrt((X2-X1)**2+(Y2-Y1)**2)

                # if it is within the initial search space AND not next to the first ridge starter
                if Search_space[x+X2-1, y+Y2-1] != 0 and Distance > np.sqrt(2):
                    Ridges[x+X2-1, y+Y2-1] = Order
                    Slope_copy[x+X2-1, y+Y2-1] = 0

                # Otherwise, look for second ridge starter elsewhere in the kernel
                elif Search_space[x+X2-1, y+Y2-1] != 0 and Distance <= np.sqrt(2):
                    for j in np.arange(0,9,1):
                        Kernel_slope_copy[X2, Y2] = 0

                        Ridge_starter2 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
                        X2=Ridge_starter2[0][0]; Y2=Ridge_starter2[1][0]
                        Distance = np.sqrt((X2-X1)**2+(Y2-Y1)**2)

                        if Search_space[x+X2-1, y+Y2-1] != 0 and Distance > np.sqrt(2):
                            Ridges[x+X2-1, y+Y2-1] = Order
                            Slope_copy[x+X2-1, y+Y2-1] = 0
                            break


        # 2/ If there are two peaks, we have one ridge starter
        elif np.count_nonzero(Kernel_ridges) == 2:
            Ridge_starter1 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
            X1=Ridge_starter1[0][0]; Y1=Ridge_starter1[1][0]

            # if it is within the initial search space
            if Search_space[x+X1-1, y+Y1-1] != 0:
                Ridges[x+X1-1, y+Y1-1] = Order
                Slope_copy[x+X1-1, y+Y1-1] = 0

    return Ridges, Slope_copy






#-----------------------------------------------------------------------------------------------------
def Continue_ridge (Slope, Search_space, Peaks, Order):
    """
    This function is the third and final stage of a routing process used to identify lines of maximum slopes.
    IMPORTANT: this function is meant to be run several times! It requires the incrementation of the Order value with each iteration.
    This function identifies multiple elements in an array (Slope), within a predefined search space (Search_space) and within the neighbourhood of the local maxima identified in a second input array (Peaks).  The identified elements are given a value of Order. To make this function work, the input array Slope should be the output array Slope_copy of the function initiate_ridge.

    Args:
        Slope (2D numpy array): the input 2-D array, here issued from a slope raster where the elements selected in the initiate_ridge function have been replaced by 0.
        Search_space (2D numpy array): the search space array.
        Peaks (2D numpy array): A 2-D array containing elements with a value of 1. These elements have the same indices as the elements with a value of 0 in Slope.
        Order (int): the value given to the identified elements. On the first iteration it should be superior by 1 to the value of Order in the function initiate_ridge. the value of Order then needs to be incremented with every iteration.

    Returns:
        Ridges (2D numpy array): a 2-D array where the identified elements have a value of Order. This array is modified from the Peaks array and therefore also contains elements of a value equal to the Order in the functions peak_flag and initiate_ridge.
        Slope_copy (2D numpy array): a copy of the input array where the value of the selected elements has been set to 0.

    Author: GCHG
    """

    print ' ... Prolongating ridges ...'
    Slope_copy = np.copy(Slope) # the copy of the initial slope array
    Search = np.where(Search_space == 1) # the searched locations
    Search_peaks = np.where(Peaks == Order-1) # the searched locations where the peaks are
    Ridges = np.copy(Peaks)

    # Define Kernels
    for i in range(len(Search_peaks[0])):
        x=Search_peaks[0][i]; y=Search_peaks[1][i] # coordinates of the kernel's centre

        Kernel_slope = kernel (Slope, 3, x, y)
        Kernel_slope_copy = kernel (Slope_copy, 3, x, y)
        Kernel_ridges = kernel (Ridges, 3, x, y)
        Kernel_search = kernel (Search_space, 3, x, y)

        # Count the number of nonzero points in the kernel of the ridge array
        Ridge_count = np.count_nonzero(Kernel_ridges)

        # If there are only the 2 previous ridge points, draw a third point that is far enough from the previous point
        if Ridge_count == 2:
            New_point = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
            X=New_point[0][0]; Y=New_point[1][0]
            Grandad_point = np.where (Kernel_ridges == Order-2)
            Xgd=Grandad_point[0][0]; Ygd=Grandad_point[1][0]
            Distance = np.sqrt((X-Xgd)**2+(Y-Ygd)**2)

            if Search_space[x+X-1, y+Y-1] != 0 and Distance > np.sqrt(2):
                Ridges[x+X-1, y+Y-1] = Order
                Slope_copy[x+X-1, y+Y-1] = 0

            elif Search_space[x+X-1, y+Y-1] != 0 and Distance <= np.sqrt(2):
                for j in np.arange(0,9,1):
                    Kernel_slope_copy[X, Y] = 0

                    New_point = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
                    X=New_point[0][0]; Y=New_point[1][0]
                    Distance = np.sqrt((X-Xgd)**2+(Y-Ygd)**2)

                    if Search_space[x+X-1, y+Y-1] != 0 and Distance > np.sqrt(2):
                        Ridges[x+X-1, y+Y-1] = Order
                        Slope_copy[x+X-1, y+Y-1] = 0
                        break

    return Ridges, Slope_copy




#-----------------------------------------------------------------------------------------------------
def Clean_ridges (Peaks, DEM, Nodata_value, opt):
    """
    This function eliminates some of the ridges (Peaks) identified by the trio of functions (peak_flag, initiate_ridge and continue_ridge). The elimination process depends on local relief, which uses a DEM (DEM) and a threshold value (opt). It is set to ignore elements with a value of  Nodata_value.

    Args:
        Peaks (2D numpy array): the input 2-D arraym which is the output of the ridge identification process.
        DEM (2D numpy array): the DEM array used as a base for the elimination of unnecessary ridges.
        Nodata_value (float): The value for ignored elements.
        opt (float): The value of the threshold to eliminate unnecessary ridges.

    Returns:
        Peaks (2D numpy array): a 2-D array much like the input Peaks array, but the unnecessary elemets have been reset to 0.

    Author: GCHG
    """

    print "Cleaning up ridges ..."
    DEM_copy = np.copy(DEM)
    DEM_copy[DEM_copy==Nodata_value] = 0
    Search_ridge = np.where (Peaks != 0)

    Cutoff = np.percentile(DEM_copy,75)
    Threshold = np.amax(DEM_copy[DEM_copy<Cutoff])
    DEM_copy[DEM_copy>Threshold]=Threshold

    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
        Kernel_DEM = kernel (DEM_copy, 9, x, y)
        Kernel_DEM[Kernel_DEM==Nodata_value]=0

        if np.amax(Kernel_DEM)/Threshold < opt:
            Peaks[x,y] = 0

    Search_ridge = np.where (Peaks != 0)
    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
        Kernel_ridges = kernel (Peaks, 9, x, y)
        # If there aren't at least 8 ridge points in the neighbourhood of 10 by 10
        if np.count_nonzero(Kernel_ridges) < 8:
            Peaks[x,y] = 0

    return Peaks




#-----------------------------------------------------------------------------------------------------
def Fill_marsh (DEM, Peaks, Nodata_value, opt):
    """
    This function builds a marsh platform array by using the Peaks array as a starting point. It uses the DEM array to establish conditions on the elements to select. the opt parameter sets a threshold value to eliminate superfluous elements. It is set to ignore elements with a value of Nodata_value.

    Args:
        DEM (2D numpy array): the DEM array.
        Peaks (2D numpy array): the 2-D array of ridge elements, which is the output of the ridge identification and cleaning process.
        Nodata_value (float): The value for ignored elements.
        opt (float): The value of the threshold to eliminate unnecessary elements.

    Returns:
        Marsh (2D numpy array): a 2-D array where the marsh platform elements are identified by strictly positive values. Other elements have a valuof 0 or Nodata_value.

    Author: GCHG
    """

    print "Initiate platform ..."
    DEM_copy = np.copy(DEM)
    Marsh = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)

    Counter = 1
    Search_ridges = np.where (Peaks > 0)
    for i in range(len(Search_ridges[0])):
        x=Search_ridges[0][i]; y=Search_ridges[1][i]
        Kernel_ridges = kernel (Peaks, 3, x, y)
        Kernel_DEM = kernel (DEM, 3, x, y)

        Marsh_point = np.where (np.logical_and (Kernel_DEM >= Kernel_DEM[1,1], Kernel_ridges == 0))
        for j in range(len(Marsh_point[0])):
            X=Marsh_point[0][j]; Y=Marsh_point[1][j]
            Marsh[x+X-1, y+Y-1] = Counter

    Search_marsh_start = np.where (Marsh == 1)
    for i in range(len(Search_marsh_start[0])):
        x=Search_marsh_start[0][i]; y=Search_marsh_start[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        Kernel_ridges = kernel (Peaks, 3, x, y)
        if np.count_nonzero(Kernel_marsh) <=2:
            Marsh[x,y] = 0

    print ' ... Build platform ...'
    while Counter < 100:
        Counter = Counter+1
        Search_marsh = np.where (Marsh == Counter-1)
        for i in range(len(Search_marsh[0])):
            x = Search_marsh[0][i]; y = Search_marsh[1][i]
            Kernel_DEM = kernel (DEM, 3, x, y)
            Kernel_DEM_copy = kernel (DEM_copy, 3, x, y)
            Kernel_ridges = kernel (Peaks, 3, x, y)
            Kernel_marsh = kernel (Marsh, 3, x, y)
            Big_Kernel_DEM = kernel (DEM, 11, x, y)
            Big_Kernel_DEM_copy = kernel (DEM_copy, 11, x, y)


            Conditions = np.zeros((len(Kernel_DEM), len(Kernel_DEM[0,:])), dtype = np.float)
            # 1: free space
            Condition_1 = np.where (np.logical_and(Kernel_ridges == 0, Kernel_marsh == 0)); Conditions[Condition_1] = 1
            # 2: not topped
            Condition_2 = np.where (np.logical_and(Kernel_DEM_copy > np.amax(Big_Kernel_DEM_copy)-0.2, Conditions == 1)); Conditions[Condition_2] = 2


            #This is a distance thing to make sure you don't cross the ridges agin
            Here_be_ridges = np.where (Kernel_ridges != 0)
            Here_be_parents = np.where (Kernel_marsh == Counter-1)

            for j in range(len(Condition_2[0])):
                X=Condition_2[0][j]; Y=Condition_2[1][j]
                Distance_to_ridges = []
                Distance_to_parents = []

                for k in range(len(Here_be_ridges[0])):
                    Xr=Here_be_ridges[0][k]; Yr=Here_be_ridges[1][k]
                    Distance = np.sqrt((X-Xr)**2+(Y-Yr)**2)
                    Distance_to_ridges.append(Distance)

                for k in range(len(Here_be_parents[0])):
                    Xp=Here_be_parents[0][k]; Yp=Here_be_parents[1][k]
                    Distance = np.sqrt((X-Xp)**2+(Y-Yp)**2)
                    Distance_to_parents.append(Distance)

                if len(Distance_to_ridges)>0:
                    if min(Distance_to_ridges) > min(Distance_to_parents):
                        Marsh[x+X-1, y+Y-1] = Counter
                else:
                    Marsh[x+X-1, y+Y-1] = Counter
                    DEM_copy[x+X-1, y+Y-1] = 0


    print ' ... defining the elimination of low platforms ...'
    Platform = np.copy(Marsh)
    Platform[Platform > 0] = DEM [Platform > 0]
    Platform_bins, Platform_hist = Distribution(Platform,0)

    #1. Find the highest and biggest local maximum of frequency distribution
    # Initialize Index
    Index = len(Platform_hist)-1
    # Initiate Cutoff_Z value
    Cutoff_Z = 0

    for j in range(1,len(Platform_hist)-1):
        if Platform_hist[j]>0.9*max(Platform_hist) and Platform_hist[j]>Platform_hist[j-1] and Platform_hist[j]>Platform_hist[j+1]:
            Index  = j

    #2. Now run a loop from there toward lower elevations.
    Counter = 0
    for j in range(Index,0,-1):
        # See if you cross the mean value of frequency. Count for how many indices you are under.
        if Platform_hist[j] < np.mean(Platform_hist):
            Counter = Counter + 1
        # Reset the counter value if you go above average again
        else:
            Counter = 0

        #If you stay long enough under (10 is arbitrary for now), initiate cutoff and stop the search
        if Counter > opt:
            Cutoff = j
            Cutoff_Z = Platform_bins[Cutoff]
            break

    # If you stay under for more than 5, set a Cutoff_Z value but keep searching
    if Counter > opt/2:
        Cutoff = j
        Cutoff_Z = Platform_bins[Cutoff]

    Marsh[Platform<Cutoff_Z] = 0


    print " ... Fill high areas left blank ..."
    Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_marsh = np.where (DEM >= Platform_bins[Index])
    Search_marsh_condition [Search_marsh] = 1
    Search_marsh_2 = np.where (np.logical_and(Marsh == 0, Search_marsh_condition == 1))
    Marsh[Search_marsh_2] = 3

    print ' ... Fill the interior of pools ...'
    for Iteration in np.arange(0,10,1):
        Counter = 100
        while Counter > 3:
            #print Counter
            Counter = Counter-1
            Search_marsh = np.where (Marsh == Counter+1)
            Non_filled = 0
            for i in range(len(Search_marsh[0])):
                x = Search_marsh[0][i]; y = Search_marsh[1][i]
                Kernel_DEM = kernel (DEM, 3, x, y)
                Kernel_ridges = kernel (Peaks, 3, x, y)
                Kernel_marsh = kernel (Marsh, 3, x, y)

                if Non_filled <len(Search_marsh[0]):
                    if np.count_nonzero(Kernel_marsh) > 6:
                        Condition = np.where (np.logical_and(Kernel_marsh == 0, Kernel_ridges == 0))
                        for j in range(len(Condition[0])):
                            X=Condition[0][j]; Y=Condition[1][j]
                            Marsh[x+X-1, y+Y-1] = Counter
                    else:
                        Non_filled = Non_filled + 1

    print 'done with pools'
    # Reapply the cutoff because the straight line thing is ugly
    Platform = np.copy(Marsh)
    Platform[Platform > 0] = DEM [Platform > 0]
    Marsh[Platform<Cutoff_Z] = 0



    # We fill in the wee holes
    Search_marsh = np.where (np.logical_and(Marsh == 0, Peaks == 0))
    for i in range(len(Search_marsh[0])):
        #print 'hole No ', i
        x = Search_marsh[0][i]; y = Search_marsh[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.count_nonzero(Kernel_marsh) == 8:
            Marsh[x,y] = 105



    print ' ... Adding the ridges'
    # We get rid of scarps that do not have a marsh next to them
    Search_false_scarp = np.where (Peaks > 0)
    for i in range(len(Search_false_scarp[0])):
        x = Search_false_scarp[0][i]; y = Search_false_scarp[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.count_nonzero (Kernel_marsh) == 0:
            Peaks[x, y] = 0

    # We get rid of the sticky-outy bits
    Search_ridge = np.where (Peaks > 0)
    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i]
        Kernel_ridges = kernel (Peaks, 9, x, y)
        if np.count_nonzero(Kernel_ridges) < 8:
            Peaks[x,y] = 0

    # We put the scarps in the platform
    Search_side = np.where (Peaks > 0)
    Marsh[Search_side] = 110

    print " ... eliminate patches of empty elements ..."
    Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_marsh = np.where (DEM >= Platform_bins[Index])
    Search_marsh_condition [Search_marsh] = 1
    Search_marsh_2 = np.where (np.logical_and(Marsh == 0, Search_marsh_condition == 1))
    Marsh[Search_marsh_2] = 3

    print ' ... Fill the interior of pools ...'
    for Iteration in np.arange(0,10,1):
        Counter = 110
        while Counter > 3:
            #print Counter
            Counter = Counter-1
            Search_marsh = np.where (Marsh == Counter+1)
            Non_filled = 0
            for i in range(len(Search_marsh[0])):
                x = Search_marsh[0][i]; y = Search_marsh[1][i]
                Kernel_DEM = kernel (DEM, 3, x, y)
                Kernel_ridges = kernel (Peaks, 3, x, y)
                Kernel_marsh = kernel (Marsh, 3, x, y)

                if Non_filled <len(Search_marsh[0]):
                    if np.count_nonzero(Kernel_marsh) > 6:
                        Condition = np.where (np.logical_and(Kernel_marsh == 0, Kernel_ridges == 0))
                        for j in range(len(Condition[0])):
                            X=Condition[0][j]; Y=Condition[1][j]
                            Marsh[x+X-1, y+Y-1] = Counter
                    else:
                        Non_filled = Non_filled + 1

    print ' ... defining the elimination of low platforms ...'
    Platform = np.copy(Marsh)
    Platform[Platform > 0] = DEM [Platform > 0]
    Marsh[Platform<Cutoff_Z] = 0

    Marsh[DEM == Nodata_value] = Nodata_value

    return Marsh





#---------------------------------------------------------------
def MARSH_ID (DEM, Slope, Nodata_value, opt1, opt2, opt3):
    """
    This is the master function for marsh identification. It defines in which order the functions define_search_space, peak_flag, initiate_ridge, Continue_ridge, Clean_ridges, Fill_marsh are executed. It is set to repeat the iteration of the Continue_ridge function 50 times.

    Args:
        DEM (2D numpy array): the input DEM array.
        Slope (2D numpy array): the input Slope array.
        Nodata_value (float): The value for ignored elements.
        opt1 (float): The value of the threshold used in the define_search_space function.
        opt2 (float): The value of the threshold used in the Clean_ridges function.
        opt3 (float): The value of the threshold used in the Fill_marsh function.

    Returns:
        Search_space (2D numpy array): The output search space of the define_search_space function.
        Ridge (2D numpy array): The output ridges of the peak_flag, initiate_ridge, Continue_ridge, Clean_ridges functions.
        Marsh (2D numpy array): The output marsh platform of the Fill_marsh function.

    Author: GCHG
    """

    DEM_work = np.copy(DEM); Slope_work = np.copy(Slope);

    Platform = np.copy(DEM_work)
    Ridge = np.copy(DEM_work)
    Marsh = np.copy(DEM_work)

    Platform[Platform != Nodata_value] = 0
    Summit = np.where (Platform==np.amax(Platform))
    Platform[Summit] = 1


    Search_space, Crossover, bins, hist, Inflexion_point = define_search_space (DEM_work, Slope_work, Nodata_value,opt1)

    Order = 1
    Ridge, Slope_temp = peak_flag (Slope_work, Search_space, Order)

    Order = Order+1
    Ridge, Slope_temp = initiate_ridge (Slope_temp, Search_space, Ridge, Order)

    while Order < 50:
        Order = Order+1
        Ridge, Slope_temp = Continue_ridge (Slope_temp, Search_space, Ridge, Order)

    Ridge = Clean_ridges (Ridge, DEM_work, Nodata_value, opt2)

    Marsh = Fill_marsh (DEM_work, Ridge, Nodata_value, opt3)


    print "My hovercraft is full of eels!"
    print


    return Search_space, Ridge, Marsh




#-----------------------------------------------------------------------------------------------------
def Confusion (Subject, Reference, Nodata_value):
    """
    This function compares a Subject 2-D array to a Reference 2-D array and returns an array of differences, which we call a confusion array or confusion map if it look like a map. It then calculates a number of metrics relative to the adequation between the subject and the reference. It is set to ignore elements with a value of Nodata_value.

    To learn more about confusion matrices and their associated metrics, please visit the Wikipedia page: https://en.wikipedia.org/wiki/Confusion_matrix

    Args:
        Subject (2D numpy array): the input array. This is the one you want to test
        Reference (2D numpy array): the reference array. This one is supposed to contain correct information
        Nodata_value (float): The value for ignored elements.

    Returns:
        Confusion_matrix (2D numpy array): an array containing the values 1 (True Positive), 2 (True Negative), -1 (False Positive) and -2 (False Negative).
        Performance (1D numpy array): the number of (respectively) True Positives, True Negatives, False Positives and False Negatives in Confusion_matrix.
        Metrix (1D numpy array): The values of (respectively) Accuracy, Reliability, Sensitivity, F1 derived from the Performance array.

    Author: GCHG
    """

    Height = len(Subject[:,0]); Width = len(Subject[0,:])
    Height_R = len(Reference[:,0]); Width_R = len(Reference[0,:])

    print Height, Width
    print Height_R, Width_R

    H = min (Height, Height_R)
    W = min (Width, Width_R)

    Confusion_matrix = Nodata_value*np.ones((Height, Width), dtype = np.float)

    Subject_marsh = np.where (np.logical_and(Subject != 0, Subject != Nodata_value))
    Reference_marsh = np.where (np.logical_and(Reference != 0, Reference != Nodata_value))

    Subject[Subject_marsh] = 1.
    Reference[Reference_marsh] = 1.

    for i in range (H):
        for j in range (W):
            if Subject[i,j] == 1 and Reference[i,j] == 1: # TRUE POSITIVE
                Confusion_matrix[i,j] = 1
            elif Subject[i,j] == 0 and Reference[i,j] == 0: # TRUE NEGATIVE
                Confusion_matrix[i,j] = 2
            elif Subject[i,j] == 1 and Reference[i,j] == 0: # FALSE POSITIVE
                Confusion_matrix[i,j] = -1
            elif Subject[i,j] == 0 and Reference[i,j] == 1: # FALSE NEGATIVE
                Confusion_matrix[i,j] = -2

    True_positive = np.sum(Confusion_matrix[Confusion_matrix == 1])
    True_negative = np.sum(Confusion_matrix[Confusion_matrix == 2])/2
    False_positive = -np.sum(Confusion_matrix[Confusion_matrix == -1])
    False_negative = -np.sum(Confusion_matrix[Confusion_matrix == -2])/2

    Reliability = True_positive / (True_positive+False_positive)
    Sensitivity = True_positive / (True_positive+False_negative)
    Accuracy = (True_positive+True_negative) / (True_positive+True_negative+False_positive+False_negative)
    F1 = 2*True_positive/(2*True_positive+False_positive+False_negative)

    Performance = np.array([True_positive,True_negative,False_positive,False_negative])
    Metrix = np.array([Accuracy, Reliability, Sensitivity, F1])


    return Confusion_matrix, Performance, Metrix



"""def True_Outline (Raster, Raster_val, Outline_value, Nodata_value):

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

    return Raster_3"""
