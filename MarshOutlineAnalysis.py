"""
MarshOutlineAnalysis.py

This script drives marsh outline analysis.

Authors: Guillaume CH Goodwin, Simon M. Mudd, University of Edinburgh
         Fiona J. Clubb, Institute of something in Potsdam
         Jaap Nienhuis, Florida State University


"""

# First import the mecessary modules
import os
import sys
import LSDMarshOutlineAnalysis as MOA


#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to run the LSD Marsh Outline Analysis (MOA).\n")
    print("Use the -dir flag to define the working directory.")
    print("If you don't I will assume the data is in the directory '/Example_data'.")
    print("Use the -driver flag to define the driver file name.")
    print("If you don't I will use the most recent driver file in the working directory.\n")
    print("For help type:")
    print("   python MarshOutlineAnalysis.py -h\n")
    print("=======================================================================\n\n ")

#=============================================================================
# This is the main function that runs the whole thing
#=============================================================================
def main(argv):

    # print("On some windows systems you need to set an environment variable GDAL_DATA")
    # print("If the code crashes here it means the environment variable is not set")
    # print("Let me check gdal enviroment for you. Currently is is:")
    # print(os.environ['GDAL_DATA'])
    #os.environ['GDAL_DATA'] = os.popen('gdal-config --datadir').read().rstrip()
    #print("Now I am going to get the updated version:")
    #print(os.environ['GDAL_DATA'])

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()
    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, default = "/Example_Data/", help="The base directory with the inputs for MOA. If this isn't defined I'll assume it's the directory '/Example_Data/'.")
    parser.add_argument("-site", "--site_name",type=str,default = "", help = "This is the prefix of the site files that selects the site for MOA. This file shoud be stored in the base directory. Default = no site")

    # Do you want to run the analysis?
    parser.add_argument("-MOA", "--MarshOutlineAnalysis", type=bool, default=True, help="If this is true, I will run the MOA algorithm")

    # Do you want plots?
    parser.add_argument("-Plots", "--MOA_plots", type=bool, default=False, help="If this is true I'll plot the results of MOA")

    args = parser.parse_args()

    site = []
    if not args.site_name:
        print("WARNING! You haven't supplied your driver file name. Please specify this with the flag '-site'")
        sys.exit()
    else:
        print("The site you want to analyse is: ")
        site = [str(item) for item in args.site_name.split(',')]
        print(site)

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
        print("You gave me the base directory:")
        print(this_dir)
    else:
        this_dir = os.getcwd()
        print("You didn't give me a directory. I am using the directory:")
        print(this_dir)

    # Run the analysis if you want it
    if args.MarshOutlineAnalysis:
        MOA.MarshOutlineAnalysis(Input_dir = this_dir, Output_dir = this_dir+'/Output/', Site=site)

    # make the plots depending on your choices
    if args.MOA_plots:
        MOA.Plot_platform_on_hillshade(Input_dir = this_dir, Output_dir = this_dir+'/Output/Figures/', Site=site)
        MOA.Plot_marsh_outline_on_hillshade(Input_dir = this_dir, Output_dir = this_dir+'/Output/Figures/', Site=site)
        MOA.Plot_Elevation_PDF(Input_dir = this_dir, Output_dir = this_dir+'/Output/Figures/', Site=site)


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
    