import sys
import argparse
from demust import __version__ as cp_vers
from demust.scripts.maps import mapsApp
from demust.scripts.plots import plotsApp
from demust.scripts.compare import compareApp
from demust.scripts.convert import convertApp
from demust.scripts.removegaps import removeGapsApp
from demust.scripts.riesselman import riesselmanApp

#TODO:

def usage_main():
    """
    Show how to use this program!
    """
    print("""
Example usage:
demust -h
Demust contains five apps:
 - maps
 - plots
 - compare
 - convert
 - removegaps
You can get more information about each individual app as follows:
demust maps -h
demust plots -h
demust compare -h
demust convert -h
demust removegaps -h
""")


def main():

    print("""
|-------------------------------------------demust----------------------------------------------------|
|                                                                                                     |
| demust       :  A Python toolkit to modify, visualize and analyze deep mutational scanning (DMS)    |                                
|                 data of proteins.                                                                   |
| Copyright   (C) Mustafa Tekpinar, 2022                                                              |
| Address      :  UMR 7238 CNRS - LCQB, Sorbonne Université, 75005 Paris, France                      |          
| Email        :  tekpinar@buffalo.edu                                                                |
| Licence      :  GNU LGPL V3                                                                         |
|                                                                                                     |
| Documentation:                                                                                      |
| Citation     : .....................................................................................|      
| Version      : {0}                                                                                |
|-----------------------------------------------------------------------------------------------------|
""".format(cp_vers))

    main_parser = argparse.ArgumentParser(description=\
        "A Python toolkit to modify, visualize and analyze deep mutational scanning (DMS) data of proteins.")
    subparsers = main_parser.add_subparsers(dest='command')

    #maps script argument parsing
    maps_parser  = subparsers.add_parser('maps')

    maps_parser.add_argument('-i', '--inputfile', dest='inputfile', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)

    maps_parser.add_argument('-d', '--datatype', dest='datatype', type=str, \
        help='gemme, rhapsody, foldx, evmutation or proteingym', \
        required=False, default='gemme')

    maps_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Name of the output file without file extension. Default extension is png.', \
        required=False, default='output')

    maps_parser.add_argument('--offset', dest='offset', type=int, \
        help='An integer value to offset the xlabels for incomplete sequences',
        required=False, default=0)

    maps_parser.add_argument('--colormap', dest='colormap', type=str, \
        help='A colormap as defined in matplotlib',
        required=False, default='coolwarm_r')

    maps_parser.add_argument('--paginate', dest='paginate', type=str, \
        help='A true or false value',
        required=False, default='false')    

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
        help='A True or False value to apply rank normalization to data matrix',
        required=False, default=False)

    maps_parser.add_argument('-s', '--sequence', dest='sequence', type=str, \
        help='Input sequence file in fasta format', \
        required=False, default=None)

    maps_parser.add_argument('--aaorder', dest='aaorder', type=str, \
        help='Amino acid order as a single string. Default is alphabetical: \"ACDEFGHIKLMNPQRSTVWY\"', \
        required=False, default='ACDEFGHIKLMNPQRSTVWY')

    #plots script argument parsing
    plots_parser = subparsers.add_parser('plots')
    plots_parser.add_argument('-i', '--inputfile', dest='inputfile', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    plots_parser.add_argument('-d', '--datatype', dest='datatype', type=str, \
        help='gemme, rhapsody, foldx, evmutation, riesselman, jet or pdb.', \
        required=False, default='gemme')
    plots_parser.add_argument('-b', '--beginning', dest='beginning', type=int, \
        help='An integer to indicate the first residue index.',
        required=False, default=1)
    plots_parser.add_argument('-e', '--end', dest='end', type=int, \
        help='An integer to indicate the final residue index.',
        required=False, default=None)
    plots_parser.add_argument('-t', '--type', dest='type', type=str, \
        help='Type of the 2D data that you want to extract. \n It can be min or max for DMS maps. It can be trace, pc or cv for JET files.', \
        required=False, default=False)
    plots_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Name of the output file.', \
        required=False, default='output.txt')


    #compare script argument parsing
    compare_parser = subparsers.add_parser('compare')
    compare_parser.add_argument('-i', '--inputfile1', dest='inputfile1', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    compare_parser.add_argument('--itype', dest='itype', type=str, \
        help='gemme, rhapsody, foldx or evmutation. Default is gemme.', \
        required=False, default='gemme')
    compare_parser.add_argument('-j', '--inputfile2', dest='inputfile2', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    compare_parser.add_argument('--jtype', dest='jtype', type=str, \
        help='gemme, rhapsody, foldx or evmutation. Default is gemme.', \
        required=False, default='gemme')
    compare_parser.add_argument('-m', '--metric', dest='metric', type=str, \
        help='Comparison metric.\n It can be spearman or pearson. Default is spearman', \
        required=False, default="spearman")
    compare_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Name of the output file. Default is output.txt', \
        required=False, default='output.txt')

    #convert script argument parsing
    convert_parser = subparsers.add_parser('convert')
    convert_parser.add_argument('-i', '--inputfile', dest='inputfile', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    convert_parser.add_argument('--itype', dest='itype', type=str, \
        help='gemme, rhapsody, foldx or evmutation. Default is gemme.', \
        required=False, default='foldx')
    convert_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Default is gemme', \
        required=True, default='output.txt')
    convert_parser.add_argument('--otype', dest='otype', type=str, \
        help='gemme or dms (standard format). Default is gemme.', \
        required=False, default='gemme')
    convert_parser.add_argument('--aaorder', dest='aaorder', type=str, \
        help='Amino acid order as a single string. Default is alphabetical: \"ACDEFGHIKLMNPQRSTVWY\"', \
        required=False, default='ACDEFGHIKLMNPQRSTVWY')

    #removegaps script argument parser
    removegaps_parser = subparsers.add_parser('removegaps')

    removegaps_parser.add_argument('-i', '--inputgappedmsa', dest='inputgappedmsa', type=str, \
        help="Name of the input gapped Multiple Sequence Alignment file in fasta format.",
        required=True, default=None)

    removegaps_parser.add_argument('-o', '--outputungappedmsa', dest='outputungappedmsa', type=str, \
        help="Name of the output gapped Multiple Sequence Alignment file in fasta format. Default is outputungappedmsa.fasta",
        required=False, default="outputungappedmsa.fasta")

    # riesselman script argument parsing
    riesselman_parser = subparsers.add_parser('riesselman', description=\
        "This script will parse and plot experimental DMS maps from Riesselman, 2016.")

    riesselman_parser.add_argument('-d', '--dataset', dest='dataset', type=str, \
        help='Name of the dataset that you want to retrieve.', \
        required=True, default=None)
    
    riesselman_parser.add_argument('-e', '--experiment', dest='experiment', type=str, \
        help='Experiment type field such as fitness, abundance, screenscore.', \
        required=True, default="DMS_score")

    riesselman_parser.add_argument('--colormap', dest='colormap', type=str, \
        help='A colormap as defined in matplotlib',
        required=False, default='coolwarm_r')

    riesselman_parser.add_argument('-s', '--sequence', dest='sequence', type=str, \
        help='One of the input sequence file in fasta format', \
        required=False, default=None)

    riesselman_parser.add_argument('-m', '--metric', dest='metric', type=str, \
        help='Comparison metric.\n It can be spearman or pearson. Default is spearman', \
        required=False, default="spearman")
    
    riesselman_parser.add_argument('-o', '--output', dest='output', type=str, \
        help='Name of the output file prefix for the csv and the png file. Default is output', \
        required=False, default='output')

    args = main_parser.parse_args()

    if args.command == "maps":
        mapsApp(args)
    elif args.command == "plots":
        plotsApp(args)
    elif args.command == "compare":
       compareApp(args)
    elif args.command == "convert":
       convertApp(args)
    elif args.command == "removegaps":
       removeGapsApp(args)
    elif args.command == "riesselman":
       riesselmanApp(args)
    elif args.command == "-h" or args.command == "--help":
        usage_main()
    else:
        usage_main()
        sys.exit(-1)
    
if __name__ == "__main__":
    main()
