"""
LSDMOA_plots.py

This file plots stuff.

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
import cPickle as pickle
import timeit
import os
from random import randint
import pandas as bb
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


from LSDMOA_classes import *
from LSDMOA_functions import *


def Plot_transects_hs_and_profiles (Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Site = ["FEL"]):
    """
    This plots the extracted marsh platform on a hillshade

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites

    Author: GCHG
    """
    #Plot 1: Draw the platform on a DEM, superimposed on a hillshade
    for site in Site:

        Nodata_value = -1000

        basedir = Input_dir
        outdir = basedir + 'Output/'

        tidedir = basedir + 'Tide/'
        figdir = outdir + 'Figures/'
        metricsdir = outdir + 'Marsh_metrics/'


        print(" Loading DEM")
        DEM_fname = site+"_DEM.bil"
        DEM, post_DEM, envidata_DEM =  fct.ENVI_raster_binary_to_2d_array (basedir+DEM_fname)

        print " Loading Marsh Platform"
        # Make sure we have the right slope file name
        marsh_fname = site+"_Marsh.bil"
        Marsh, post_Marsh, envidata_Marsh =  fct.ENVI_raster_binary_to_2d_array (basedir+marsh_fname)

        print " Loading Hillshade"
                # Make sure we have the right slope file name
        hs_fname = site+"_hs.bil"
        hs, post_hs, envidata_hs =  ENVI_raster_binary_to_2d_array (Input_dir+hs_fname)

        Marsh_object = Marsh_platform(Marsh.shape[0], Marsh.shape[1])
        Marsh_object = Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = True)
        Marsh_DEM = Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = False)
        Marsh_labels = Marsh_object.label_connected (Nodata_value)

        print 'Loading outlines from Pickle'
        Outlines = pickle.load( open(outdir+'Marsh_metrics/'+str(site)+'_Out.pkl', "rb" ) )

        print 'Loading transects from Pickle'
        Transects = pickle.load( open(outdir+'Marsh_metrics/'+str(site)+'_Tr.pkl', "rb" ) )
        Transects_mean = pickle.load( open(outdir+'Marsh_metrics/'+str(site)+'_mean.pkl', "rb" ) )
        Transects_stdev = pickle.load( open(outdir+'Marsh_metrics/'+str(site)+'_stdev.pkl', "rb" ) )

        #Transects.plot_transects_basemap_and_profiles(DEM, Marsh_object, hs, figdir, site + '_combined', Nodata_value)

        twin = DEM.copy()
        twinmask = Marsh_object.copy()

        fig_height = min(np.floor(twin.shape[1])/5, 25)
        fig_width = min(np.floor(twin.shape[1])/5, 10)
        fig=plt.figure('figname', facecolor='White',figsize=[fig_height,fig_width])

        ax1 = plt.subplot2grid((1,2),(0,0),colspan=1, rowspan=1)
        ax1.tick_params(axis='x', colors='black')
        ax1.tick_params(axis='y', colors='black')
        plt.xlabel('x (m)', fontsize=18)
        plt.ylabel('y (m)', fontsize=18)

        ax1.annotate('a.', fontsize = 14,
             xy=(0.95, 0.05),
             xycoords='axes fraction',color='white')

        # configure the basemap
        Vmin = min(np.amin(twin[twin>Nodata_value])*0.95, np.amin(twin[twin>Nodata_value])*1.05)
        Vmax = max(np.amax(twin)*0.95, np.amax(twin)*1.05)
        twinmask = np.ma.masked_where(twinmask == 1, twinmask)

        Map = ax1.imshow(hs, interpolation='None', cmap=plt.cm.gray, vmin=0, vmax=210, alpha = 1.0)
        Map = ax1.imshow(twin, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax, alpha = 0.6)
        Map = ax1.imshow(twinmask, interpolation='None', cmap= plt.cm.gray, vmin=Vmin, vmax=Vmax, alpha = 0.4)

        ax2 = fig.add_axes([0.13, 0.98, 0.345, 0.02])
        scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
        cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)
        cb1.set_label('Elevation (m a.m.s.l.)', fontsize=18)


        # Draw the lines, panda style
        colour = 0
        for i in range(len(Transects)):
            Pandaline = Transects[i]
            if Pandaline.size > 0:
                colour+=1
                L_labels = range (1,max(Pandaline['L_code'])+1)
                for L in L_labels:
                    Pandaline_slice = Pandaline.loc[Pandaline['L_code'] == L]
                    #print Pandaline_slice
                    if Pandaline_slice.size > 0:
                        T_labels = range (1,max(Pandaline_slice['T_code'])+1)
                        for T in T_labels:
                            Pandaline_zest = Pandaline_slice.loc[Pandaline_slice['T_code'] == T]
                            To_draw = Pandaline_zest.prepare_for_plotting(plt.cm.jet(colour*50),opaque = Pandaline_zest['select'].iloc[0])
                            ax1.add_line(To_draw)



        # Make the other plot
        ax3 = plt.subplot2grid((1,2),(0,1),colspan=1, rowspan=1)
        ax3.tick_params(axis='x', colors='black')
        ax3.tick_params(axis='y', colors='black')
        plt.ylabel('Elevation (m a.m.s.l.)', fontsize=18)
        plt.xlabel('Transect length (m)', fontsize=18)

        ax3.annotate('b.', fontsize = 14,
             xy=(0.05, 0.05),
             xycoords='axes fraction')

        ax3.set_xlim(0,20)
        #ax3.set_ylim(ymin=-2, ymax = 10)

        colour = 0
        for i in range(len(Transects)):
            Pandaline = Transects[i]
            Pandaline_m = Transects_mean[i]
            Pandaline_s = Transects_stdev[i]
            if Pandaline.size > 0:
                colour+=1
                L_labels = range (1,max(Pandaline['L_code'])+1)
                for L in L_labels:
                    Pandaline_slice = Pandaline.loc[Pandaline['L_code'] == L]
                    #print Pandaline_slice
                    if Pandaline_slice.size > 0:
                        T_labels = range (1,max(Pandaline_slice['T_code'])+1)
                        for T in T_labels:
                            Pandaline_zest = Pandaline_slice.loc[Pandaline_slice['T_code'] == T]
                            if Pandaline_zest['select'].iloc[0] == True:
                                ax3.plot(Pandaline_zest['Z'], color = plt.cm.jet(colour*50),alpha = 0.15)

                if Pandaline_m['select'].iloc[0] == True:
                    ax3.plot(Pandaline_m['Z'], color = plt.cm.jet(colour*50),alpha = 1.0,linewidth = 4.0)

        ax3.axvline(10, color='black', lw=1.0, alpha=0.8)

        ax4 = fig.add_axes([0.55, 0.98, 0.345, 0.02])
        scheme = plt.cm.jet; norm = colors.Normalize(vmin=0, vmax=colour)
        cb2 = matplotlib.colorbar.ColorbarBase(ax4, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)
        cb2.set_label('Platform number', fontsize=18)

        plt.savefig(figdir+site+'_Combined_Tr.png', bbox_inches='tight')







def Plot_all_site_stats (Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Site = ["FEL"]):

    #Hack your inputs!!!!
    Map_Site = ['BOU','FEL','CRO','SHE', 'WOR_domain_4', 'HEY','HIN']
    Tide_Site = ['BOU','FEL','CRO','SHE', 'WOR','HEY','HIN']

    Map_Site = ['BOU','FEL','CRO','SHE','HEY','HIN']
    Tide_Site = ['BOU','FEL','CRO','SHE','HEY','HIN']

    fig_height = 20
    fig_width = 10

    #fig=plt.figure('figname', facecolor='White',figsize=[fig_height,fig_width])

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(fig_height, fig_width))

    axes[0,0] = plt.subplot2grid((2,3),(0,0),colspan=1, rowspan=1)

    #axes[0,0] = plt.subplot2grid((1,2),(0,0),colspan=1, rowspan=1)
    axes[0,0].tick_params(axis='x', colors='black')
    axes[0,0].tick_params(axis='y', colors='black')
    plt.ylabel('Elevation / Spring Tidal Range', fontsize=18)
    plt.xlabel('Transect length (m)', fontsize=18)
    axes[0,0].annotate('a.', fontsize = 14, xy=(0.95, 0.05), xycoords='axes fraction')

    axes[0,1] = plt.subplot2grid((6,2),(0,1),colspan=1, rowspan=5)
    plt.ylabel('Elevation (m a.m.s.l.)', fontsize=14)
    plt.xlabel('Spring Tidal Range (m)', fontsize=14)
    axes[0,1].annotate('b.', fontsize = 14, xy=(0.05, 0.95), xycoords='axes fraction')

    #ax1.set_ylim(ymin = -1, ymax = 8)

    #ax2 = fig.add_axes([0.13, 0.98, 0.345, 0.02])
    #scheme = plt.cm.gist_earth; norm = colors.Normalize(vmin=Vmin, vmax=Vmax)
    #cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)
    #cb1.set_label('Elevation (m a.m.s.l.)', fontsize=18)

    """ax3 = plt.subplot2grid((1,2),(0,1),colspan=1, rowspan=1)
    ax3.tick_params(axis='x', colors='black')
    ax3.tick_params(axis='y', colors='black')
    plt.ylabel('Local slope (n/a)', fontsize=18)
    plt.xlabel('Transect length (m)', fontsize=18)
    ax3.annotate('b.', fontsize = 14, xy=(0.95, 0.05), xycoords='axes fraction')"""

    i = 0
    for site in Map_Site:
        axes[0,0] = plt.subplot2grid((6,2),(i,0),colspan=1, rowspan=1)
        axes[0,0].tick_params(axis='x', colors='black')
        axes[0,0].tick_params(axis='y', colors='black')
        if i == 3:
            plt.ylabel('Elevation (m a.m.s.l.)', fontsize=14)
        if i == 5:
            plt.xlabel('Transect length (m)', fontsize=14)
        axes[0,0].annotate('a.'+str(i+1), fontsize = 14, xy=(0.95, 0.8), xycoords='axes fraction')

        if i<5:
            axes[0,0].set_xticks([])

        Nodata_value = -1000
        basedir = Input_dir
        outdir = basedir + 'Output/'
        tidedir = basedir + 'Tide/'
        figdir = outdir + 'Figures/'
        metricsdir = outdir + 'Marsh_metrics/'

        print 'Loading transect stats from Pickle'
        Transect_mean = pickle.load( open(metricsdir+str(site)+'_mean.pkl', "rb" ) )
        Transect_stdev = pickle.load( open(metricsdir+str(site)+'_stdev.pkl', "rb" ) )
        Big_Tr = pickle.load( open(metricsdir+str(site)+'_BigTr.pkl', "rb" ) )
        Big_mean = pickle.load( open(metricsdir+str(site)+'_Bigmean.pkl', "rb" ) )
        Big_stdev = pickle.load( open(metricsdir+str(site)+'_Bigstdev.pkl', "rb" ) )
        print 'Loading tidal stats from Pickle'
        Tidemetrix = pickle.load( open(tidedir+str(Tide_Site[i])+'_Metric2.pkl', "rb" ) )
        Tide = bb.DataFrame()
        Tide['LWS'] = [np.mean(Tidemetrix[0])]; Tide['LWN'] = [np.mean(Tidemetrix[1])]
        Tide['HWN'] = [np.mean(Tidemetrix[2])]; Tide['HWS'] = [np.mean(Tidemetrix[3])]
        TR = Tide['HWS']-Tide['LWS']; TR = TR.iloc[0]; print TR
        Big_Tr = Big_Tr.add_attribute_list('TR',TR)
        Draw = Big_Tr.loc[Big_Tr['select']==True]


        colour = plt.cm.jet(TR/10)

        r = axes[0,1].violinplot(np.asarray(Draw['Z']), [TR], points=20, widths=0.3, showmeans=True, showextrema=False, showmedians=False)#, colour = plt.cm.jet(TR/10))
        r['cmeans'].set_color('k')
        for vp in r['bodies']:
            vp.set_facecolor(colour)
            #vp.set_edgecolor(rrred)
            #vp.set_linewidth(1)
            vp.set_alpha(0.3)

        rs = axes[0,1].violinplot(np.asarray(Draw['Z'].loc[0]), [TR], points=20, widths=0.3, showmeans=True, showextrema=False, showmedians=False)
        rs['cmeans'].set_color('k')
        for vp in rs['bodies']:
            vp.set_facecolor(colour)
            #vp.set_edgecolor(rrred)
            #vp.set_linewidth(1)
            vp.set_alpha(0.6)
        re = axes[0,1].violinplot(np.asarray(Draw['Z'].loc[20]), [TR], points=20, widths=0.3, showmeans=True, showextrema=False, showmedians=False)
        re['cmeans'].set_color('k')
        for vp in re['bodies']:
            vp.set_facecolor(colour)
            #vp.set_edgecolor(rrred)
            #vp.set_linewidth(1)
            vp.set_alpha(0.6)

        """rd = ax3.violinplot(np.asarray(Draw['Z'].loc[0])-Draw['Z'].loc[20], [TR], points=20, widths=0.3, showmeans=True, showextrema=False, showmedians=False)
        rd['cmeans'].set_color(colour)
        for vp in rd['bodies']:
            vp.set_facecolor(colour)
            #vp.set_edgecolor(rrred)
            #vp.set_linewidth(1)
            vp.set_alpha(0.6)"""



        M_counter = 0
        for M in range(len(Transect_mean)):
            if Transect_mean[M].size > 0:
                to_draw = Transect_mean[M]
                #to_fill = Transect_stdev[M]
                #for i in range(len(max(to_draw.index))):
                if min(to_draw['Z']) > -10:
                    if to_draw['select'].iloc[0] == True:
                        #ax3.plot(Pandaline_m['Z'], color = plt.cm.jet(colour*50),alpha = 1.0,linewidth = 4.0)
                        axes[0,0].plot(to_draw['Z'], color = plt.cm.jet(TR/10+M/20), alpha = 0.20)
                        #ax1.plot(to_draw['Z'], color = plt.cm.jet(i*40))
                        #ax3.plot(to_draw['dZ'], color = plt.cm.jet(TR/10), alpha = 0.15)
                    #ax1.fill_between(range(0,21), min(to_draw['Z']), to_draw['Z'],color = plt.cm.jet(i*50), alpha = 0.3)
                    #ax1.fill_between(range(0,21), to_draw['Z']-to_fill['Z'], to_draw['Z']+to_fill['Z'],color = plt.cm.jet(i*50), alpha = 0.3)
                M_counter+=1

        axes[0,0].plot(Big_mean['Z'], color = plt.cm.jet(TR/10), alpha = 0.9,linewidth = 2.5)
        #ax3.plot(Big_mean['dZ'], color = plt.cm.jet(TR/10), alpha = 0.9,linewidth = 3.0)

        axes[0,0].fill_between([9.95, 10.05], Tide['LWS'], Tide['HWS'],color='k', lw=0.0, alpha=0.8)
        axes[0,0].fill_between([9.9, 10.1], Tide['LWN'], Tide['HWN'],color='k', lw=0.0, alpha=0.8)

        axes[0,0].set_ylim(ymin = 0, ymax = 1.3*max(Tide['HWS'].iloc[0], max(to_draw['Z']) ))
        majorFormatter = FormatStrFormatter('%i')
        axes[0,0].yaxis.set_major_formatter(majorFormatter)

        i+=1

    ax4 = fig.add_axes([0.55, 0.13, 0.345, 0.02])
    scheme = plt.cm.jet; norm = colors.Normalize(vmin=0, vmax=11.7)
    cb2 = matplotlib.colorbar.ColorbarBase(ax4, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)
    cb2.set_label('Spring Tidal Range (m)', fontsize=18)





    plt.savefig(figdir+'BIGSTUFF.png', bbox_inches='tight')







def Plot_all_site_slopes (Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Site = ["FEL"]):

    #Hack your inputs!!!!
    Map_Site = ['BOU','FEL','CRO','SHE', 'WOR_domain_4', 'HEY','HIN']
    Tide_Site = ['BOU','FEL','CRO','SHE', 'WOR','HEY','HIN']

    Map_Site = ['BOU','FEL','CRO','SHE','HEY','HIN']
    Tide_Site = ['BOU','FEL','CRO','SHE','HEY','HIN']

    fig_height = 10
    fig_width = 8

    #fig=plt.figure('figname', facecolor='White',figsize=[fig_height,fig_width])

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(fig_height, fig_width))


    axes[0,0] = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=5)
    plt.ylabel('Local slope ', fontsize=14)
    plt.xlabel('Transect length (m)', fontsize=14)


    i = 0
    for site in Map_Site:

        Nodata_value = -1000
        basedir = Input_dir
        outdir = basedir + 'Output/'
        tidedir = basedir + 'Tide/'
        figdir = outdir + 'Figures/'
        metricsdir = outdir + 'Marsh_metrics/'

        print 'Loading transect stats from Pickle'
        Transect_mean = pickle.load( open(metricsdir+str(site)+'_mean.pkl', "rb" ) )
        Transect_stdev = pickle.load( open(metricsdir+str(site)+'_stdev.pkl', "rb" ) )
        Big_Tr = pickle.load( open(metricsdir+str(site)+'_BigTr.pkl', "rb" ) )
        Big_mean = pickle.load( open(metricsdir+str(site)+'_Bigmean.pkl', "rb" ) )
        Big_stdev = pickle.load( open(metricsdir+str(site)+'_Bigstdev.pkl', "rb" ) )
        print 'Loading tidal stats from Pickle'
        Tidemetrix = pickle.load( open(tidedir+str(Tide_Site[i])+'_Metric2.pkl', "rb" ) )
        Tide = bb.DataFrame()
        Tide['LWS'] = [np.mean(Tidemetrix[0])]; Tide['LWN'] = [np.mean(Tidemetrix[1])]
        Tide['HWN'] = [np.mean(Tidemetrix[2])]; Tide['HWS'] = [np.mean(Tidemetrix[3])]
        TR = Tide['HWS']-Tide['LWS']; TR = TR.iloc[0]; print TR
        Big_Tr = Big_Tr.add_attribute_list('TR',TR)
        Draw = Big_Tr.loc[Big_Tr['select']==True]


        colour = plt.cm.jet(TR/10)

        M_counter = 0
        for M in range(len(Transect_mean)):
            if Transect_mean[M].size > 0:
                to_draw = Transect_mean[M]
                #to_fill = Transect_stdev[M]
                #for i in range(len(max(to_draw.index))):
                if min(to_draw['Z']) > -10:
                    if to_draw['select'].iloc[0] == True:
                        #ax3.plot(Pandaline_m['Z'], color = plt.cm.jet(colour*50),alpha = 1.0,linewidth = 4.0)
                        axes[0,0].plot(to_draw['dZ'], color = plt.cm.jet(TR/10), alpha = 0.10)
                        #ax1.plot(to_draw['Z'], color = plt.cm.jet(i*40))
                        #ax3.plot(to_draw['dZ'], color = plt.cm.jet(TR/10), alpha = 0.15)
                    #ax1.fill_between(range(0,21), min(to_draw['Z']), to_draw['Z'],color = plt.cm.jet(i*50), alpha = 0.3)
                    #ax1.fill_between(range(0,21), to_draw['Z']-to_fill['Z'], to_draw['Z']+to_fill['Z'],color = plt.cm.jet(i*50), alpha = 0.3)
                M_counter+=1

        axes[0,0].plot(Big_mean['dZ'], color = plt.cm.jet(TR/10), alpha = 0.9,linewidth = 2.5)
        #axes[0,0].plot(Draw['dZ'], color = plt.cm.jet(TR/10), alpha = 0.1)#,linewidth = 2.5)
        #ax3.plot(Big_mean['dZ'], color = plt.cm.jet(TR/10), alpha = 0.9,linewidth = 3.0)

        #axes[0,0].fill_between([9.95, 10.05], Tide['LWS'], Tide['HWS'],color='k', lw=0.0, alpha=0.8)
        #axes[0,0].fill_between([9.9, 10.1], Tide['LWN'], Tide['HWN'],color='k', lw=0.0, alpha=0.8)

        #axes[0,0].set_ylim(ymin = 0, ymax = 1.3*max(Tide['HWS'].iloc[0], max(to_draw['Z']) ))
        #majorFormatter = FormatStrFormatter('%i')
        #axes[0,0].yaxis.set_major_formatter(majorFormatter)

        i+=1

    axes[0,0].axvline(10, color='black', lw=1.0, alpha=0.8)

    axes[0,0].set_ylim(ymin = -0.8,ymax = 0.4)

    ax4 = fig.add_axes([0.15, 0.83, 0.345, 0.02])
    scheme = plt.cm.jet; norm = colors.Normalize(vmin=0, vmax=11.7)
    cb2 = matplotlib.colorbar.ColorbarBase(ax4, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)
    cb2.set_label('Spring Tidal Range (m)', fontsize=14)





    plt.savefig(figdir+'BIGSTUFF2.png', bbox_inches='tight')
    quit()






    # Draw the lines, panda style
    colour = 0
    for i in range(len(self)):
        Pandaline = self[i]
        if Pandaline.size > 0:
            colour+=1
            L_labels = range (1,max(Pandaline['L_code'])+1)
            for L in L_labels:
                Pandaline_slice = Pandaline.loc[Pandaline['L_code'] == L]
                #print Pandaline_slice
                if Pandaline_slice.size > 0:
                    T_labels = range (1,max(Pandaline_slice['T_code'])+1)
                    for T in T_labels:
                        Pandaline_zest = Pandaline_slice.loc[Pandaline_slice['T_code'] == T]
                        To_draw = Pandaline_zest.prepare_for_plotting(plt.cm.jet(colour*50),opaque = Pandaline_zest['select'].iloc[0])
                        ax1.add_line(To_draw)

    ax1.fill_between(range(0,21), self[h][i][-2]-self[h][i][-1], self[h][i][-2]+self[h][i][-1], alpha = 0.3)

    # Make the other plot
    ax3 = plt.subplot2grid((1,2),(0,1),colspan=1, rowspan=1)
    ax3.tick_params(axis='x', colors='black')
    ax3.tick_params(axis='y', colors='black')
    plt.ylabel('Elevation (m a.m.s.l.)', fontsize=16)
    plt.xlabel('Transect length (m)', fontsize=16)

    ax3.annotate('b.', fontsize = 14,
         xy=(0.05, 0.05),
         xycoords='axes fraction')

    ax3.set_xlim(0,20)
    #ax3.set_ylim(ymin=-2, ymax = 10)

    colour = 0
    for i in range(len(self)):
        Pandaline = self[i]
        if Pandaline.size > 0:
            colour+=1
            L_labels = range (1,max(Pandaline['L_code'])+1)
            for L in L_labels:
                Pandaline_slice = Pandaline.loc[Pandaline['L_code'] == L]
                #print Pandaline_slice
                if Pandaline_slice.size > 0:
                    T_labels = range (1,max(Pandaline_slice['T_code'])+1)
                    for T in T_labels:
                        Pandaline_zest = Pandaline_slice.loc[Pandaline_slice['T_code'] == T]
                        if Pandaline_zest['select'].iloc[0] == True:
                            ax3.plot(Pandaline_zest['Z'], color = plt.cm.jet(colour*50))

    ax3.axvline(10, color='black', lw=1.0, alpha=0.8)

    ax4 = fig.add_axes([0.55, 0.98, 0.345, 0.02])
    scheme = plt.cm.jet; norm = colors.Normalize(vmin=0, vmax=colour)
    cb2 = matplotlib.colorbar.ColorbarBase(ax4, cmap=scheme, norm=norm, orientation='horizontal', alpha = 0.6)
    cb2.set_label('Platform number', fontsize=18)

























    for i in range(10):

        fig=plt.figure(1, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)

        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)


        # Make a mask!
        Platform_mask = np.ma.masked_where(Platform <=0, Platform)
        Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]

        # Make a map!
        #Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
        #Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)

        plt.savefig(Output_dir+'Platform_DEM_%s.png' % (site))



        #quit()

        """# Set the value for empty DEM cells
        Nodata_value = -1000

        Sites = ['BOU','FEL','CRO','HEY','HIN','SHE']

        basedir = '/home/willgoodwin/Software/LSDCoastal/LSD_MOA/Example_Data/'
        outdir = basedir + 'Output/'

        tidedir = basedir + 'Tide/'
        figdir = outdir + 'Figures/'
        metricsdir = outdir + 'Marsh_metrics/'


        Outlines = pickle.load( open(outdir+'Marsh_metrics/'+'BOU'+'_Out.pkl', "rb" ) )

        quit()"""




    for site in Sites:
        print(" Loading DEM")
        DEM_fname = site+"_DEM.bil"
        DEM, post_DEM, envidata_DEM =  fct.ENVI_raster_binary_to_2d_array (basedir+DEM_fname)

        print " Loading Marsh Platform"
        # Make sure we have the right slope file name
        marsh_fname = site+"_Marsh.bil"
        Marsh, post_Marsh, envidata_Marsh =  fct.ENVI_raster_binary_to_2d_array (basedir+marsh_fname)


        Marsh_object = cl.Marsh_platform(Marsh.shape[0], Marsh.shape[1])
        Marsh_object = Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = True)
        Marsh_DEM = Marsh_object.set_attribute (Marsh, 1, DEM, Nodata_value, classification = False)

        print 'Loading outlines from Pickle'
        Outlines = pickle.load( open(outdir+'Marsh_metrics/'+str(site)+'_Out.pkl', "rb" ) )

        print 'Loading outlines from Pickle'
        cl.Transects = pickle.load( open(outdir+'Marsh_metrics/'+str(site)+'_Tr.pkl', "rb" ) )


        Transects.plot_transects_basemap_and_profiles(DEM, Marsh_object, figdir, site + '_combined', Nodata_value)


    quit()





"""def MarshOutlineAnalysis(Input_dir =  "/Example_Data/",
            Output_dir = "/Example_Data/Output/",
            Site = ["FEL_DEM_clip"], opt1 = -2.0, opt2 = 0.85, opt3 = 8.0):

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
        print Marsh.shape"""



def Plot_platform_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This plots the extracted marsh platform on a hillshade

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites

    Author: GCHG
    """
    #Plot 1: Draw the platform on a DEM, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(1, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)

        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)


        # Make a mask!
        Platform_mask = np.ma.masked_where(Platform <=0, Platform)
        Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]

        # Make a map!
        #Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
        #Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)

        plt.savefig(Output_dir+'Platform_DEM_%s.png' % (site))





def Plot_Elevation_PDF(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This plots the extracted marsh platform on a hillshade

    Args:
        Input_dir (str): Name your data input directory
        Output_dir (str): Name your results output directory
        Sites (str list): A list of strings. The file names are modified based on these sites

    Author: GCHG
    """
    #Plot 1: Draw the platform on a DEM, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(1, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)

        # Name the axes
        ax1.set_xlabel('Elevation (m)', fontsize = 12)
        ax1.set_ylabel('Probability Distribution (m)', fontsize = 12)

        # Load the relevant data
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)
        Platform [Platform>0] = DEM [Platform>0]

        bins_z, hist_z = Distribution (DEM, Nodata_value)
        bins_M, hist_M = Distribution (Platform, Nodata_value)

        Elevation_range_z = np.arange(min(bins_z[bins_z!=Nodata_value]), max(bins_z), 0.1)
        Elevation_range_M = np.arange(min(bins_z[bins_z!=Nodata_value]), max(bins_M), 0.1)
        Ratio = (max(hist_z[hist_z < 0.2])/max(hist_M[hist_M < 0.2]))
        hist_z_copy = hist_z / Ratio
        hist_M[0] = 0


        ax1.fill_between( bins_z, -5, hist_z_copy, color=plt.cm.gray(0), alpha = 0.5, linewidth = 0.0)
        ax1.plot( bins_M, hist_M, '-r', linewidth = 2.0)


        # Set the ticks
        A = 0.01
        #for x in range(len(hist_M)-1):
            #if hist_M[x]==0 and hist_M[x+1]>0:
                #A = bins_M[x]
                #break
        #xmin = max(-5,A)
        ymax = max(hist_M[hist_M<0.2])

        #ax1.set_xlim (xmin = xmin)
        ax1.set_ylim (ymin = 0, ymax = ymax*1.05)

        #majorLocator1 = MultipleLocator(np.floor(100*ymax)/200)
        #ax1.yaxis.set_major_locator(majorLocator1)
        #majorLocator2 = MultipleLocator(1)
        #ax1.xaxis.set_major_locator(majorLocator2)


        plt.savefig(Output_dir+'Elevation_PDF_%s.png' % (site))






def Plot_marsh_outline_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This draws the marsh outline on a DEM, superimposed on a hillshade.

    Args:
        Input_dir (str): Name your data input directory.
        Output_dir (str): Name your results output directory.
        Sites (str list): A list of strings. The file names are modified based on these sites.

    Author: GCHG
    """

    for site in Sites:
        fig=plt.figure(2, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)

        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)

        # Outline the marsh
        Platform[Platform > 0] = 1
        Marsh_outline = Outline (Platform, 2, Nodata_value)


        # Make a mask!
        Outline_mask = np.ma.masked_where(Marsh_outline <=1, Marsh_outline)

        # Make a map!
        Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)

        Map_Marsh = ax1.imshow(Outline_mask, interpolation='None', cmap=plt.cm.Reds, vmin = 0, vmax = 2, alpha = 1)


    plt.savefig(Output_dir+'Platform_outline_%s.png' % (site))


def Plot_marsh_reference_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This draws the marsh reference on a hillshade

    Args:
        Input_dir (str): Name your data input directory.
        Output_dir (str): Name your results output directory.
        Sites (str list): A list of strings. The file names are modified based on these sites.

    Author: GCHG
    """
    #Plot 3: Draw the marsh map and reference outline, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(3, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)


        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array (Input_dir+"%s.bil" % (site), site)
        Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Marsh.bil" % (site), site)
        Reference, post_Reference, envidata_Reference = ENVI_raster_binary_to_2d_array (Input_dir+"%s_ref.bil" % (site), site)

        # Outline the reference
        Reference[Reference > 0] = 1
        Ref_outline = Outline (Reference, 2, Nodata_value)


        # Make a mask!
        Outline_mask = np.ma.masked_where(Ref_outline <=1, Ref_outline)


        # Make a map!
        Platform_mask = np.ma.masked_where(Platform <=0, Platform)
        Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]

        Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
        Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)

        Map_Marsh = ax1.imshow(Outline_mask, interpolation='None', cmap=plt.cm.Reds, vmin = 0, vmax = 2, alpha = 1)

    plt.savefig(Output_dir+'Reference_map_%s.png' % (site))



def Plot_confusion_map_on_hillshade(Input_dir =  "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Output_dir = "/LSDTopoTools/LSDTopoTools_MarshPlatform/Example_data/",
            Sites = ["FEL"]):
    """
    This draws the marsh reference on a hillshade

    Args:
        Input_dir (str): Name your data input directory.
        Output_dir (str): Name your results output directory.
        Sites (str list): A list of strings. The file names are modified based on these sites.

    Author: GCHG
    """
    #Plot 4: Draw the confusion map, superimposed on a hillshade
    for site in Sites:
        fig=plt.figure(4, facecolor='White',figsize=[10,10])
        ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)


        # Name the axes
        ax1.set_xlabel('x (m)', fontsize = 12)
        ax1.set_ylabel('y (m)', fontsize = 12)

        # Load the relevant data
        HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array (Input_dir+"%s_hs.bil" % (site), site)
        Confusion, post_Confusion, envidata_Confusion = ENVI_raster_binary_to_2d_array (Output_dir+"%s_Confusion.bil" % (site), site)


        # Make a mask!
        Confusion_mask = np.ma.masked_where(Confusion == Nodata_value, Confusion)


        # Make a map!
        Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
        Map_DEM = ax1.imshow(Confusion_mask, interpolation='None', cmap=plt.cm.Spectral, vmin = -2, vmax = 2, alpha = 0.5)


    plt.savefig(Output_dir+'Confusion_%s.png' % (site))
