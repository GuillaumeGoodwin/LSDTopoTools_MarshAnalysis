"""
    This package processes topographic data in order to extract marsh platforms
    Guillaume C.H. Goodwin
    Released unedr GPL3
"""

# Load useful Python packages
import os
import sys
import numpy as np

#from osgeo import gdal, osr, gdalconst
import matplotlib.pyplot as plt
#from osgeo.gdalconst import *
import cPickle

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D

import timeit

import LSDMOA_functions as fct

import copy

##########################################################################################################
##########################################################################################################
class Land_surface (np.ndarray):
    def __new__ (Land_surface, x_length, y_length):
        print 'In __new__ with class %s' % Land_surface
        return np.ndarray.__new__(Land_surface, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        self[np.isnan(self)] = 0
        self = 0 * self

    # We have no copy method because the ndarray.copy method is good enough

    def set_attribute (self, select_array, select_value, attribute_array, Nodata_value, classification = False):
        new_array = self.copy()

        Select = np.where (select_array >= select_value)
        Discard = np.where (np.logical_and(select_array < select_value, select_array != Nodata_value))
        Nodata = np.where (attribute_array == Nodata_value)

        new_array[Discard] = 0
        new_array[Nodata] = Nodata_value

        if classification == False:
            new_array[Select] = attribute_array [Select]
        else:
            new_array[Select] = select_value

        return new_array


    def label_connected (self, Nodata_value):
        new_array = self.copy()
        Ref = self.copy()

        import scipy.ndimage as scim
        array, numfeat = scim.label(self)

        for value in np.arange(1, np.amax(array)+1, 1):
            line = np.where(array == value)
            new_array[line] = value
            new_array[Ref == Nodata_value] = Nodata_value

        return new_array


    def plot_map (self, save_dir, figname, title, Nodata_value):
        print ' Plotplotplotplot'
        twin  = self.copy()

        fig_height = min(np.floor(twin.shape[1])/5, 50)
        fig_width = min(np.floor(twin.shape[1])/5, 50)

        fig=plt.figure(title, facecolor='White',figsize=[fig_height,fig_width])

        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        Vmin = min(np.amin(twin[twin!=Nodata_value])*0.95, np.amin(twin[twin!=Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)

        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        ax2 = fig.add_axes([0.1, 0.98, 0.85, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)

        plt.savefig(save_dir+figname+'.png')



##########################################################################################################
##########################################################################################################
class Marsh_platform (Land_surface):
    def __new__ (Marsh_platform, x_length, y_length):
        print 'In __new__ with class %s' % Marsh_platform
        return np.ndarray.__new__(Marsh_platform, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        self[np.isnan(self)] = 0
        self = 0 * self


    def calc_labelled_area (self, label_value):
        Selected = np.where(self == label_value)
        Labelled_area = len(Selected[0])

        return Labelled_area




##########################################################################################################
##########################################################################################################
class Marsh_outline (Land_surface):
    def __new__ (Marsh_outline, x_length, y_length):
        print 'In __new__ with class %s' % Marsh_outline
        return np.ndarray.__new__(Marsh_outline, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        self[np.isnan(self)] = 0
        self = 0 * self

    def complete_outline_from_array (self, array):
        print '\nExtracting outline...'
        Start = timeit.default_timer()
        new_array = self.copy()
        new_array, Outline_value = fct.Surface_outline (new_array, array, 2)
        Stop = timeit.default_timer()
        print '  Extraction runtime : ', Stop - Start , 's'

        return new_array, Outline_value

    def tightrope_outline_from_array (self, array):
        print '\nExtracting outline...'
        Start = timeit.default_timer()
        new_array = self.copy()
        new_array, Outline_value = fct.Tightrope_outline (new_array, array, 2)
        Stop = timeit.default_timer()
        print '  Extraction runtime : ', Stop - Start , 's'

        return new_array, Outline_value



    def reduce_to_marsh_labels (self, Marsh_labels, Nodata_value):
        # This one relies on the fact that the longest continuous line should be the marsh outline for each label
        print '\nReducing outline...'
        Start = timeit.default_timer()

        new_array = 0 * self.copy()

        # Find out how many marsh labels you have
        M_Labels = range (int(np.amin(Marsh_labels[Marsh_labels>0])), int(np.amax(Marsh_labels[Marsh_labels>0]))+1, 1)

        # Find out how many outline labels you have
        L_Labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)

        # Make a list counting the elements for each label. The 0 index is for the label value, the 1 index is for the number of elements
        Num_Elements = [[],[]]
        for label in L_Labels:
            Elements = np.where(self == label)
            num_elements = len (Elements[0])
            Num_Elements[0].append(label); Num_Elements[1].append(num_elements)

        for i in range(len(M_Labels)):
            # If the label has a non-negligible area
            Labelled_area = Marsh_labels.calc_labelled_area (M_Labels[i])

            if Labelled_area >1 :
                Most_pop_index = np.where(np.asarray(Num_Elements[1])==max(np.asarray(Num_Elements[1])))[0][0]

                Most_pop_label = Num_Elements[0][Most_pop_index]
                new_array[self == Most_pop_label] = Most_pop_label

                # Mow remove it.
                del Num_Elements[0][Most_pop_index]; del Num_Elements[1][Most_pop_index]

        Stop = timeit.default_timer()
        print '  Reduction runtime : ', Stop - Start , 's'

        return new_array



    def trim_to_main_stem (self, Marsh_array, Nodata_value):
        print '\nTrimming the messy outlines...'
        Start = timeit.default_timer()
        new_array = self.copy()
        # Find out how many labels you have
        Labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)
        #And loop through the labels
        for lab in range(len(Labels)):
            print '  This is the label: ', lab+1, '/', len(Labels), ' (', Labels[lab], ')'
            new_array = fct.Billhook (new_array, Marsh_array, Labels[lab])
        Stop = timeit.default_timer()
        print '  Gardening runtime : ', Stop - Start , 's'

        return new_array



    def calc_outline_length (self, Scale, Short = False):
        print '\nCalculating outline length...'
        Start = timeit.default_timer()
        #Setup the environmnents
        Length_array = self.copy()
        Labels_array = self.copy()
        #Initiate the vector
        Polylines = Polyline()
        # Find out how many labels you have
        Labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)
        #And loop through thev labels
        for lab in range(len(Labels)):
            print '\nThis is the label: ', lab+1, '/', len(Labels), ' (', Labels[lab], ')'
            # Measure the length of the stitched line for this label value
            This_polyline, Length_array, Code_array = fct.Measure_all_lines (Length_array, Labels[lab], Scale)
            #Stitch the diverging starts
            print
            # Rehabilitate this in a different manner
            This_polyline, Length_array, Code_array = fct.Stitch_diverging_starts (Length_array, Labels_array, Labels[lab], This_polyline, Code_array, Scale)
            print
            #Stitch outward going branches
            This_polyline, Length_array = fct.Graft_diverging_branch (Length_array, Labels_array, Labels[lab], This_polyline, Code_array, Scale)
            print

            #Stitch inward going branches, but not yet.
            #new_array, Line_row, Line_col, Line_
            #Line_dist, Line_code = fct.Graft_converging_branch (new_array, Labels_array, Labels[lab], Line_row, Line_col, Line_dist, Line_code, Code_array, Scale)

            Polylines.append(This_polyline)

        Stop = timeit.default_timer()
        print '  Surveying runtime : ', Stop - Start , 's'

        return Length_array, Polylines




##########################################################################################################
##########################################################################################################
class Point (tuple):
    """With this class we make a point that has row and column properties and that can also hold other properties"""

    def __new__(self, row, col):
       return tuple.__new__(Point, (row, col))

    def row (self):
        row = self[0]
        return row

    def col (self):
        col = self[1]
        return col

##########################################################################################################
##########################################################################################################
class Line (list):
    def __init__(self):
        list.__init__(self)

    def start_point (self):
        return self[0]

    def end_point (self):
        if len(self)>1:
            for i in range(1,len(self)):
                if type(self[-i]) is Point:
                    A = self[-i]
                    break
        else:
            A = self[0]
        return A


    def add_attribute_list (self, attribute_list):
        if len(self[start_point:end_point]) == len(attribute_list):
            self.append(attribute_list)
        return self


    def save_to_shp (self, Envidata, Enviarray, save_dir, file_name):
        fct.Line_to_shp(self, Envidata, Enviarray, save_dir, file_name)


    def prepare_for_plotting(self,colour):
        Line_row = []; Line_col = []
        for i in range(self.index(self.start_point()), self.index(self.end_point())+1):
            Line_row.append(self[i][0]); Line_col.append(self[i][1])
        Line = Line2D(Line_col, Line_row, color = plt.cm.gist_earth(50*colour), alpha = 0.5)

        return Line



    def select_line (self, condition):

        if condition == 1:
            halfway = int(len(self[-1])/2)
            # If the marsh bit is at the beginning
            if np.count_nonzero(self[-1][0:halfway]) >= 0.9 * len(self[-1][0:halfway]):
                # If there are no other marsh bits after
                if np.count_nonzero(self[-1][halfway:]) < 0.2 * len(self[-1][halfway:]):
                    self.append(True)
                else:
                    self.append(False)
            # The same in reverse
            elif np.count_nonzero(self[-1][halfway:]) >= 0.9 * len(self[-1][halfway:]):
                # If there are no other marsh bits after
                if np.count_nonzero(self[-1][0:halfway]) < 0.2 * len(self[-1][0:halfway]):
                    self.append(True)
                else:
                    self.append(False)
            else:
                self.append(False)

        return self



    def property_indices (self):
        indices = []
        for i in range(len(self)):
            if type(self[i]) is List:
                indices.append(i)
        return indices



    def extract_values (self, basemap):
        prop_list = []
        for i in range(self.index(self.start_point()), self.index(self.end_point())+1):

            if int(self[i][0]) < basemap.shape[0] and int(self[i][1]) < basemap.shape[1]:
                value = basemap[int(self[i][0]),int(self[i][1])]
            else:
                value = 0
            prop_list.append(value)

        if type(self[-1]) is Point:
            self.append(prop_list)
        elif type(self[-1]) is bool:
            self.insert(-1, prop_list)

        return self


    def subdivide (self,sub_number):
        #Only works if the initial line only has two points.
        if self.index(self.end_point()) <2:
            sub_row = []; sub_col = []
            for sub in range(1,sub_number):
                sub_row.append(self[0][0]+sub*float(self[1][0]-self[0][0])/sub_number)
                sub_col.append(self[0][1]+sub*float(self[1][1]-self[0][1])/sub_number)

            for sub in range(0,sub_number-1):
                subpoint = Point(sub_row[sub], sub_col[sub])
                self.insert(-1,subpoint)

        return self


##########################################################################################################
##########################################################################################################
class Polyline (list):
    """A polyline object can have line objects or polyline objects inside."""
    def __init__(self):
        list.__init__(self)


    def structure (self):
        #This lists the types of objects in the polyline
        Structure_list = [type(self)]
        This_level = self[0]
        for i in range(0,10):
            Structure_list.append(type(This_level))
            if type(This_level) is Point:
                break
            else:
                This_level = This_level[0]
        return Structure_list


    def select_few_longest (self):
        New_polyline = fct.Select_few_longest(self)
        return New_polyline


    def save_to_shp (self, Envidata, Enviarray, save_dir, site_name):
        #Structure_list = self.structure()
        for label in range(len(self)):
            if len(self[label]) > 0:
                for code in range(len(self[label])):
                    if type(self[label]) is Line:
                        self[label].save_to_shp (self, Envidata, Enviarray, save_dir, site_name+"_"+str(label)+"_"+str(code))



    def transects (self, spacing, length, Envidata, Enviarray, save_dir, site_name):
        #first you must save the polylines as a shapefile shapefiles
        self.save_to_shp (Envidata, Enviarray, save_dir, site_name)
        #then use the outside function to make the tr
        All_transects = Polyline()
        # For each label
        for label in range(len(self)):
            print " \nThe label is:", label+1
            Label_transects = Polyline()
            if len(self[label]) > 0:
                for code in range(len(self[label])):
                    print "  Saving code number ", code+1
                    fct.Save_transects(save_dir+'Shapefiles/'+'%s_%s_%s.shp' % (site_name,label, code),save_dir+'Shapefiles/'+'%s_%s_%s_Tr.shp' % (site_name,label, code), spacing, length)
                    # STEP 4: Put the transect lines into our array reference system
                    Code_transects = fct.Shp_to_lines (save_dir+"Shapefiles/", "%s_%s_%s_Tr" % (site_name,label, code), Envidata, Enviarray)

                    Label_transects.append(Code_transects)

            All_transects.append(Label_transects)

        return All_transects


    def transect_properties (self, refinement, basemap):
        Structure_list = self.structure()

        if len(Structure_list) == 4 and Structure_list[-1] is Point:
            for i in range(len(self)):
                for j in range(len(self[i])):
                    self[i][j] = self[i][j].subdivide(refinement)
                    self[i][j] = self[i][j].extract_values(basemap)

        elif len(Structure_list) == 5 and Structure_list[-1] is Point:
            for h in range(len(self)):
                    for i in range(len(self[h])):
                        for j in range(len(self[h][i])):
                            self[h][i][j] = self[h][i][j].subdivide(refinement)
                            self[h][i][j] = self[h][i][j].extract_values(basemap)

        return self




    def select_transects_from_property (self, condition):
        Structure_list = self.structure()

        if len(Structure_list) == 4 and Structure_list[-1] is Point:
            for i in range(len(self)):
                for j in range(len(self[i])):
                    self[i][j].select_line(condition)

        elif len(Structure_list) == 5 and Structure_list[-1] is Point:
            for h in range(len(self)):
                    for i in range(len(self[h])):
                        for j in range(len(self[h][i])):
                            self[h][i][j].select_line(condition)

        return self




    def transect_stats (self, refinement):
        Structure_list = self.structure()

        if len(Structure_list) == 4 and Structure_list[-1] is Point:
            for i in range(len(self)):
                for j in range(len(self[i])):
                    values = self[i][j][-2]


        elif len(Structure_list) == 5 and Structure_list[-1] is Point:
            for h in range(len(self)):
                    for i in range(len(self[h])):

                        # this contains values to consider for stats
                        Polyline_values = np.zeros((len(self[h][i]), refinement+1), dtype = np.float)

                        for j in range(len(self[h][i])): # j is the index of each line in the polyline
                            for k in range(len(self[h][i][j][-2])): # k is the index of each point in the line
                                if self[h][i][j][-1] == True:
                                    Polyline_values[j,k] = self[h][i][j][-2][k]

                        Polyline_mean = np.zeros(refinement+1, dtype = np.float)
                        Polyline_stdev = np.zeros(refinement+1, dtype = np.float)


                        if len(Polyline_values) > 1:
                            if len(Polyline_values[0]) > 1:

                                for col in range(len(Polyline_values[0])):
                                    Polyline_mean[col] = np.sum(Polyline_values[:,col])/np.count_nonzero(Polyline_values[:,col])

                                    for row in range(len(Polyline_values)):
                                        if Polyline_values[row,col] != 0:
                                            Polyline_stdev[col] += (Polyline_values[row,col]-Polyline_mean[col])**2

                                    Polyline_stdev[col] = np.sqrt(Polyline_stdev[col]/ np.count_nonzero(Polyline_values[:,col]))


                        Polyline_mean[np.isnan(Polyline_mean)] = 0
                        Polyline_stdev[np.isnan(Polyline_stdev)] = 0

                        self[h][i].append(Polyline_mean)
                        self[h][i].append(Polyline_stdev)

                        #print self[h][i]
                        #print

        return self




    def plot_on_basemap(self,basemap, save_dir, figname, Nodata_value):
        twin  = basemap.copy()
        #Make the canvas
        fig_height = min(np.floor(twin.shape[1])/5, 50); fig_width = min(np.floor(twin.shape[1])/5, 50)
        fig=plt.figure(figname, facecolor='White',figsize=[fig_height,fig_width])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        # Make the basemap
        Vmin = min(np.amin(twin[twin!=Nodata_value])*0.95, np.amin(twin[twin!=Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)

        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        ax2 = fig.add_axes([0.1, 0.98, 0.85, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)

        # Draw the lines.
        Structure_list = self.structure()
        print Structure_list

        if len(Structure_list) == 4 and Structure_list[-1] is Point:
            for i in range(len(self)):
                for j in range(len(self[i])):
                    Line1 = self[i][j]
                    # This is what you plot
                    #ax1.scatter(Line1[0][0], Line1[0][1], marker  = 'o', color = 'b')
                    Line1 = Line1.prepare_for_plotting(5*j+1)
                    ax1.add_line(Line1)

        elif len(Structure_list) == 5 and Structure_list[-1] is Point:
            for h in range(len(self)):
                    for i in range(len(self[h])):
                        for j in range(len(self[h][i])):
                            Line1 = self[h][i][j]
                            # This is what you plot
                            #ax1.scatter(Line1[0][0], Line1[0][1], marker  = 'o', color = 'b')
                            Line1 = Line1.prepare_for_plotting(5*j+1)
                            ax1.add_line(Line1)

        plt.savefig(save_dir+figname+'.png')




    def plot_property(self, save_dir, figname, draw):
        fig=plt.figure(figname, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')


        # Draw the lines.
        Structure_list = self.structure()
        print Structure_list

        if len(Structure_list) == 4 and Structure_list[-1] is Point:
            for i in range(len(self)):
                for j in range(len(self[i])):
                    Line1 = self[i][j]
                    if Line1[-1] == True:
                        # This is what you plot
                        ax1.plot(Line1[draw], color = plt.cm.jet(i*20), alpha = 0.5+j/100)

        elif len(Structure_list) == 5 and Structure_list[-1] is Point:
            for h in range(len(self)):
                    for i in range(len(self[h])):
                        for j in range(len(self[h][i])):
                            Line1 = self[h][i][j]
                            if Line1[-1] == True:
                                # This is what you plot
                                ax1.plot(Line1[draw], color = plt.cm.jet(i*20), alpha = 0.5+j/100)

        plt.savefig(save_dir+figname+'.png')


    def plot_property_stats(self, save_dir, figname, draw):
        fig=plt.figure(figname, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        # Draw the lines.
        Structure_list = self.structure()

        if len(Structure_list) == 4 and Structure_list[-1] is Point:
            for i in range(len(self)):
                ax1.plot(self[i][-2], color = plt.cm.jet(i*20), alpha = 0.5/100)

        elif len(Structure_list) == 5 and Structure_list[-1] is Point:
            for h in range(len(self)):
                    for i in range(len(self[h])):
                        ax1.plot(self[h][i][-2], color = plt.cm.jet(i*20), alpha = 1)
                        if max( self[h][i][-2]) > 0:
                            ax1.fill_between(range(0,21), self[h][i][-2]-self[h][i][-1], self[h][i][-2]+self[h][i][-1], color = plt.cm.jet(i*20), alpha = 0.3)

                        print self[h][i][-2]
                        print self[h][i][-1]
                        print

        plt.savefig(save_dir+figname+'.png')







    """def select_few_longest (self, Nodata_value, num):
        # num is the number of lines you select
        empty = np.zeros(self.shape, dtype = np.float)

        twin = np.copy(self)

        values = range (np.amin(self[self>0]), np.amax(self), 1)
        line_lengths = []
        for value in values:
            line_lengths.append(len(np.where(self == value)[0]))

        line_lengths = np.asarray(line_lengths)
        Longest = np.where(line_lengths == np.amax(line_lengths))
        print values[Longest[0][0]], line_lengths[Longest[0][0]]
        array_2[array == values[Longest[0][0]]] = values[Longest[0][0]]

        if num > 0:
            for i in range(num):
                line_lengths[Longest[0][0]] = 0
                Longest = np.where(line_lengths == np.amax(line_lengths))
                print Longest[0][0]
                self[twin == values[Longest[0][0]]] = values[Longest[0][0]]

        return self"""



    """def select_few_longest (array, Nodata_value, num):
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

    return array_2"""

    def vectorise (self, Nodata_value):
        from matplotlib.lines import Line2D

        Outlines = []; Outlines_row = []; Outlines_col = []

        Labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)

        for lab in range(len(Labels)):
            Select = np.where (self == lab)

            Line_row = Select[0]; Line_col = Select[1]
            Line = Line2D(Line_col, Line_row)

            Outlines_row.append(Line_row); Outlines_col.append(Line_col)
            Outlines.append(Line)



        twin  = self.copy()

        fig=plt.figure('title', facecolor='White',figsize=[np.floor(twin.shape[1])/5,np.floor(twin.shape[0])/5])

        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')

        Vmin = min(np.amin(twin[twin!=Nodata_value])*0.95, np.amin(twin[twin!=Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)

        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        ax2 = fig.add_axes([0.1, 0.98, 0.85, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)


        for Line in Outlines:
            ax1.add_line(Line)

        plt.savefig('/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshAnalysis/Example_Data/Output/Figures/'+'TEST'+'.png')


        return Outlines, Outlines_row, Outlines_col

######################################################

######################################################
class Tidal_flat (np.ndarray):
    def __new__ (Tidal_flat, x_length, y_length):
        print 'In __new__ with class %s' % Tidal_flat
        return np.ndarray.__new__(Tidal_flat, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
######################################################




######################################################
class Grid (np.ndarray):
    def __new__ (Grid, x_length, y_length):
        print 'In __new__ with class %s' % Grid
        return np.ndarray.__new__(Grid, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length

    def add_terrace (self, mudflat_width, mudflat_depth):
        self [:, self.shape[1]-mudflat_width:] = self [:, self.shape[1]-mudflat_width:] + mudflat_depth
        return self

    def add_creek (self,creek,x_start):
        y_start = self.shape[1] - creek.shape[1]
        self[x_start:x_start+creek.shape[0], y_start:] = creek
        return self

    def add_patch (self, patch, x_start, y_start):
        self[x_start:x_start+patch.shape[0], y_start:y_start+patch.shape[1]] =  self[x_start:x_start+patch.shape[0], y_start:y_start+patch.shape[1]] + patch
        return self
