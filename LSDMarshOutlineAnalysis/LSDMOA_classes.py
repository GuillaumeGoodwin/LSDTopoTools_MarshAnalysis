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

import timeit

import LSDMOA_functions as fct


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

                #print Most_pop_index

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

        new_array = self.copy()
        Labels_array = self.copy()
        #Initiate the vectors
        Lines_row = []; Lines_col = []; Lines_dist = []; Lines_code = []
        # Find out how many labels you have
        Labels = range (int(np.amin(self[self>0])), int(np.amax(self[self>0]))+1, 1)

        #And loop through thev labels
        for lab in range(len(Labels)):
            print '\nThis is the label: ', lab+1, '/', len(Labels), ' (', Labels[lab], ')'

            # Measure the length of the stitched line for this label value
            new_array, Line_row, Line_col, Line_dist, Line_code, Code_array = fct.Measure_all_lines (new_array, Labels[lab], Scale)
            #Stitch the diverging starts
            print
            new_array, Line_row, Line_col, Line_dist, Line_code = fct.Stitch_diverging_starts (new_array, Labels_array, Labels[lab], Line_row, Line_col, Line_dist, Line_code, Code_array, Scale)
            print
            #Stitch outward going branches
            new_array, Line_row, Line_col, Line_dist, Line_code = fct.Graft_diverging_branch (new_array, Labels_array, Labels[lab], Line_row, Line_col, Line_dist, Line_code, Code_array, Scale)

            #Stitch inward going branches
            #new_array, Line_row, Line_col, Line_dist, Line_code = fct.Graft_converging_branch (new_array, Labels_array, Labels[lab], Line_row, Line_col, Line_dist, Line_code, Code_array, Scale)

            Lines_row.append (Line_row); Lines_col.append (Line_col); Lines_dist.append (Line_dist); Lines_code.append (Line_code)


        Stop = timeit.default_timer()
        print '  Surveying runtime : ', Stop - Start , 's'

        return new_array, Lines_row, Lines_col, Lines_dist, Lines_code




##########################################################################################################
##########################################################################################################
class Point (tuple):
    def __new__(self, x, y):
       return tuple.__new__(Point, (x, y))

##########################################################################################################
##########################################################################################################
class Line (list):



##########################################################################################################
##########################################################################################################
class Land_vector (Land_surface):
    def __new__ (Marsh_outline, x_length, y_length):
        print 'In __new__ with class %s' % Marsh_outline
        return np.ndarray.__new__(Marsh_outline, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        self[np.isnan(self)] = 0
        self = 0 * self






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
