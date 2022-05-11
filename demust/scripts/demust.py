import sys
import argparse
from demust import __version__ as cp_vers
from demust.scripts.maps import mapsApp
from demust.scripts.plots import plotsApp
#from demust.scripts.compare import compareApp

#TODO:

def usage_main():
    """
    Show how to use this program!
    """
    print("""
Example usage:
demust -h
Demust contains three apps:
 - maps
 - plots
 - compare
You can get more information about each individual app as follows:
demust maps -h
demust plots -h
demust compare -h
""")


def main():

    print("""
|------------------------------------------demust---------------------------------------------------|
                                                                            
 demust       :  A Python toolkit to visualize and analyze deep mutational scanning data of proteins.                                
                                                                             
 Copyright   (C) Mustafa Tekpinar, 2022                         
 Address      :  UMR 7238 CNRS - LCQB, Sorbonne Université, 75005 Paris, France                                
 Email        :  tekpinar@buffalo.edu                          
 Licence      :  GNU LGPL V3                               
                                                                             
 Documentation:                                           
 Citation     : ....................................................................................|      
 Version      : {0}                                
|---------------------------------------------------------------------------------------------------|
""".format(cp_vers))

    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(dest='command')

    maps_parser  = subparsers.add_parser('maps')

    maps_parser.add_argument('-i', '--inputfile', dest='inputfile', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)

    maps_parser.add_argument('-d', '--datatype', dest='datatype', type=str, \
        help='gemme, rhapsody, foldx or evmutation', \
        required=False, default='gemme')

    maps_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Name of the output png file', \
        required=False, default='output.png')

    maps_parser.add_argument('--offset', dest='offset', type=int, \
        help='An integer value to offset the xlabels for incomplete sequences',
        required=False, default=0)

    maps_parser.add_argument('--colormap', dest='colormap', type=str, \
        help='A colormap as defined in matplotlib',
        required=False, default='coolwarm_r')    

    maps_parser.add_argument('--field', dest='field', type=int, \
        help='An integer value starting from 0 for reading the rhapsody output',
        required=False, default=10)
    
    maps_parser.add_argument('-b', '--beginning', dest='beginning', type=int, \
        help='An integer to indicate the first residue index.',
        required=False, default=1)

    maps_parser.add_argument('-e', '--end', dest='end', type=int, \
        help='An integer to indicate the final residue index.',
        required=False, default=None)

    maps_parser.add_argument('--ranknorm', dest='ranknorm', type=bool, \
        help='An True or False value to apply rank normalization to data matrix',
        required=False, default=False)

    plots_parser = subparsers.add_parser('plots')
    plots_parser.add_argument('-i', '--inputfile', dest='inputfile', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    plots_parser.add_argument('-d', '--datatype', dest='datatype', type=str, \
        help='gemme, rhapsody, foldx or evmutation', \
        required=False, default='gemme')
    plots_parser.add_argument('-t', '--type', dest='type', type=str, \
        help='Type of the 2D data that you want to extract. \n It can be min or max.', \
        required=False, default=False)
    plots_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Name of the output file.', \
        required=False, default='output.txt')

    
    args = main_parser.parse_args()

    if args.command == "maps":
        mapsApp(args)
    elif args.command == "plots":
        plotsApp(args)
#    elif setup.py == "compare":
        #     compareApp()
    # elif sys.argv[1] == "-h" or sys.argv[1] == "--help":
    #     usage_main()
    else:
        usage_main()
        sys.exit(-1)
    
if __name__ == "__main__":
    main()
