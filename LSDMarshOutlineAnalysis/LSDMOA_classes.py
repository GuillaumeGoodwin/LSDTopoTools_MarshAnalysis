"""
    This package processes topographic data in order to extract marsh platforms
    Guillaume C.H. Goodwin
    Released unedr GPL3
"""

# Load useful Python packages
import os
import sys
import numpy as np

from osgeo import gdal, osr, gdalconst
import matplotlib.pyplot as plt
from osgeo.gdalconst import *
import cPickle

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import LSDMOA_functions as fct

######################################################
class Land_surface (np.ndarray):
    def __new__ (Land_surface, x_length, y_length):
        print 'In __new__ with class %s' % Land_surface
        return np.ndarray.__new__(Land_surface, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        
    def set_attribute (self, select_array, select_value, array, Nodata_value, classification = False):
        Select = np.where (select_array >= select_value)
        Discard = np.where (np.logical_and(select_array < select_value, select_array != Nodata_value))
        Nodata = np.where (array == Nodata_value)
        
        self[Discard] = 0
        self[Nodata] = Nodata_value
        
        if classification == False:
            self[Select] = array [Select]
        else:
            self[Select] = 1

        return self
        
        
    def save_plot(self, save_dir, figname, title, Nodata_value):
        fig=plt.figure(title, facecolor='White',figsize=[np.floor(self.shape[1])/50,np.floor(self.shape[0])/50])

        ax1 = plt.subplot2grid((2,2),(0,0),colspan=1, rowspan=2)    
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')
        
        Vmin = min(np.amin(self[self!=Nodata_value])*0.95, np.amin(self[self!=Nodata_value])*1.05)
        Vmax = max(np.amax(self)*0.95, np.amax(self)*1.05)

        Map = ax1.imshow(self, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        ax2 = fig.add_axes([0.1, 0.92, 0.35, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)
        
                
        ax3 = plt.subplot2grid((2,2),(0,1),colspan=1, rowspan=1)
        ax3.plot (self[np.floor(self.shape[0]/2),:])
        
        ax3 = plt.subplot2grid((2,2),(1,1),colspan=1, rowspan=1)
        ax3.plot (self[:,np.floor(self.shape[1]/2)])

        plt.savefig(save_dir+figname+'.png')

        
######################################################
######################################################
class Marsh_platform (Land_surface):
    def __new__ (Marsh_platform, x_length, y_length):
        print 'In __new__ with class %s' % Marsh_platform
        return np.ndarray.__new__(Marsh_platform, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
    
    def calc_area (self, Nodata_value):
        return np.count_nonzero(self[self!=Nodata_value])
    
    def extract_scarp (self, Nodata_value):
        
        
######################################################

######################################################
class Marsh_scarp (Land_surface):
    def __new__ (Marsh_scarp, x_length, y_length):
        print 'In __new__ with class %s' % Marsh_scarp
        return np.ndarray.__new__(Marsh_scarp, shape=(x_length, y_length), dtype =np.float)
    def __init__ (self, x_length, y_length):
        self.X_length = x_length
        self.Y_length = y_length
        
    #def calc_length (self, Nodata_value):

        #return
        
        
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
        
    
        



