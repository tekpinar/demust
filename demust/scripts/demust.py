import sys
import argparse
from demust import __version__ as cp_vers
from demust.scripts.maps import mapsApp
from demust.scripts.plots import plotsApp
from demust.scripts.compare import compareApp
from demust.scripts.convert import convertApp
from demust.scripts.extract import extractApp
from demust.scripts.removegaps import removeGapsApp
from demust.scripts.riesselman import riesselmanApp
from demust.scripts.randsubset import randsubsetApp
from demust.scripts.diffmap import diffmapApp
from demust.scripts.inputgenerator import inputgeneratorApp
from demust.scripts.postprocess import postprocessApp
#TODO:

def usage_main():
    """
    Show how to use this program!
    """
    print("""
Example usage:
demust -h
Demust contains following apps:
 - maps
 - plots
 - compare
 - convert
 - extract
 - removegaps
 - randsubset
 - diffmap
 - inputgenerator
You can get more information about each individual app as follows:
demust maps -h
demust plots -h
demust compare -h
demust convert -h
demust extract -h
demust removegaps -h
demust randsubset -h
demust diffmap -h
demust inputgenerator -h
""")


def main():

    print("""
                                                                                             
| demust       :  A Python toolkit to modify, visualize and analyze deep mutational scanning (DMS)                                    
|                 data of proteins.                                                                   
| Copyright   (C) Mustafa Tekpinar, 2022-2023
| Address      :  UMR 7238 CNRS - LCQB, Sorbonne University, 75005 Paris, France                                
| Email        :  tekpinar@buffalo.edu                                                                
| Licence      :  GNU LGPL V3                                                                                                                                                                 
| Documentation:                                                                                      
| Citation     : .............................................     
| Version      : {0}

""".format(cp_vers))

    main_parser = argparse.ArgumentParser(description=\
        "A Python toolkit to prepare, modify, visualize and analyze deep mutational scanning (DMS) data of proteins.")
    subparsers = main_parser.add_subparsers(dest='command')

    #maps script argument parsing
    maps_parser  = subparsers.add_parser('maps', description="Plots colored DMS maps.")

    maps_parser.add_argument('-i', '--inputfile', dest='inputfile', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)

    maps_parser.add_argument('-d', '--datatype', dest='datatype', type=str, \
        help='gemme, rhapsody, foldx, evmutation, proteingym, singleline', \
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

    maps_parser.add_argument('--paginate', dest='paginate', type=int, \
        help='An integer value. Default is 0, which means to plot without pagination.',
        required=False, default=0)    

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
    
    maps_parser.add_argument('--iscolorbaron', dest='iscolorbaron', type=bool, \
        help='A True or False value to put color bar on the map or not. Default is False.',
        required=False, default=False)

    maps_parser.add_argument('-s', '--sequence', dest='sequence', type=str, \
        help='Input sequence file in fasta format', \
        required=False, default=None)

    maps_parser.add_argument('--aaorder', dest='aaorder', type=str, \
        help='Amino acid order as a single string. Default is alphabetical: \"ACDEFGHIKLMNPQRSTVWY\"', \
        required=False, default='ACDEFGHIKLMNPQRSTVWY')
    
    maps_parser.add_argument('--onlydnaaccessible', dest='onlydnaaccessible', type=bool, \
        help='A True or False value to mask mutations that are more than one nucleotide away. Default is False.',
        required=False, default=False)

    #plots script argument parsing
    plots_parser = subparsers.add_parser('plots', \
        description="Can plot 2D data: 'min or max from DMS maps' or 'trace, pc or cv from JET files.")
    plots_parser.add_argument('-i', '--inputfile', dest='inputfile', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    plots_parser.add_argument('-d', '--datatype', dest='datatype', type=str, \
        help='gemme, rhapsody, foldx, evmutation, riesselman, , fasta, jet or pdb.', \
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
    plots_parser.add_argument('--paginate', dest='paginate', type=int, \
        help='An integer value. Default is 0, which means to plot without pagination.',
        required=False, default=0)    

    plots_parser.add_argument('--offset', dest='offset', type=int, \
        help='An integer value to offset the xlabels for incomplete sequences',
        required=False, default=0) 


    #compare script argument parsing
    compare_parser = subparsers.add_parser('compare', \
        description="Compares two DMS results and reports Spearman correlation.")
    compare_parser.add_argument('-i', '--inputfile1', dest='inputfile1', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    compare_parser.add_argument('--itype', dest='itype', type=str, \
        help='gemme, rhapsody, foldx, evmutation or singleline. Default is gemme.', \
        required=False, default='gemme')
    compare_parser.add_argument('-j', '--inputfile2', dest='inputfile2', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    compare_parser.add_argument('--jtype', dest='jtype', type=str, \
        help='gemme, rhapsody, foldx, evmutation or singleline. Default is gemme.', \
        required=False, default='gemme')
    compare_parser.add_argument('-m', '--metric', dest='metric', type=str, \
        help='Comparison metric.\n It can be spearman or pearson. Default is spearman', \
        required=False, default="spearman")
    compare_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Name of the output file. Default is output.txt', \
        required=False, default='output.txt')
    compare_parser.add_argument('--msafile', dest='msafile', type=str, \
        help='An aligned multiple sequence alignment file in fasta format.', \
        required=False, default=None)
    compare_parser.add_argument('--ssfile', dest='ssfile', type=str, \
        help='An secondary structure file in plain text format.', \
        required=False, default=None)

    #convert script argument parsing
    convert_parser = subparsers.add_parser('convert', description="Converts outputs of various formats.")
    convert_parser.add_argument('-i', '--inputfile', dest='inputfile', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    convert_parser.add_argument('--itype', dest='itype', type=str, \
        help='gemme, rhapsody, foldx or evmutation. Default is gemme.', \
        required=False, default='gemme')
    convert_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Default is output.txt', \
        required=False, default='output.txt')
    convert_parser.add_argument('--otype', dest='otype', type=str, \
        help='gemme or dms (standard format). Default is gemme.', \
        required=False, default='gemme')
    convert_parser.add_argument('--aaorder', dest='aaorder', type=str, \
        help='Amino acid order as a single string. Default is alphabetical: \"ACDEFGHIKLMNPQRSTVWY\"', \
        required=False, default='ACDEFGHIKLMNPQRSTVWY')
    convert_parser.add_argument('-f', '--fastafile', dest='fastafile', type=str, \
        help='A fasta file of a single gene to deduce amino acid order of the reference sequence.', \
        required=False, default=None)
    convert_parser.add_argument('-b', '--beginning', dest='beginning', type=int, \
        help='An integer indicating the amount to add to the first residue index (1)."+\
            " If not specified, nothing(0) will be added.',
        required=False, default=0)

    #removegaps script argument parser
    removegaps_parser = subparsers.add_parser('removegaps', \
        description="Remove columns containing gaps in the query sequence.")

    removegaps_parser.add_argument('-i', '--inputgappedmsa', dest='inputgappedmsa', type=str, \
        help="Name of the input gapped Multiple Sequence Alignment file in fasta format.",
        required=True, default=None)

    removegaps_parser.add_argument('-o', '--outputungappedmsa', dest='outputungappedmsa', type=str, \
        help="Name of the output ungapped Multiple Sequence Alignment file in fasta format. Default is outputungappedmsa.fasta",
        required=False, default="outputungappedmsa.fasta")

    # riesselman script argument parsing
    riesselman_parser = subparsers.add_parser('riesselman', description=\
        "This script will parse and plot experimental DMS maps from Riesselman, 2016 and ProteinGym datasets.")

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

    riesselman_parser.add_argument('--source', dest='source', type=str, \
        help='Dataset source. It has two options: riesselman or proteingym. Default value is proteingym',
        required=False, default='proteingym')

    riesselman_parser.add_argument('-m', '--metric', dest='metric', type=str, \
        help='Comparison metric.\n It can be spearman or pearson. Default is spearman', \
        required=False, default="spearman")
    
    riesselman_parser.add_argument('-o', '--output', dest='output', type=str, \
        help='Name of the output file prefix for the csv and the png file. Default is output', \
        required=False, default='output')
    riesselman_parser.add_argument('--otype', dest='otype', type=str, \
        help='gemme or dat (each mutation in a single line, separated by single space). Default is gemme.', \
        required=False, default='gemme')
    riesselman_parser.add_argument('--shift', dest='shift', type=int, \
        help='An integer indicating the amount of the shift for the mutation file.',
        required=False, default=None)

        #extract script argument parsing
    extract_parser = subparsers.add_parser('extract', \
        description="Extracts single point mutation (M212A), average or alanine scanning data from a GEMME output!")
    extract_parser.add_argument('-i', '--inputfile', dest='inputfile', type=str, \
        help='One of the output files of gemme, gemme_singleline', \
        required=True, default=None)
    extract_parser.add_argument('--itype', dest='itype', type=str, \
        help='gemme, gemme_singleline. Default is gemme.', \
        required=False, default='gemme')
    extract_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Default is output.txt', \
        required=False, default='output.txt')
    extract_parser.add_argument('--otype', dest='otype', type=str, \
        help='average or mutA (which means all alanine mutations). '+ \
            'You can also extract result of a single mutation (such as M694V)'+ \
            'from a file in singleline format. Default is average.', \
        required=False, default='average')
    # extract_parser.add_argument('--aaorder', dest='aaorder', type=str, \
    #     help='Amino acid order as a single string. Default is alphabetical: \"ACDEFGHIKLMNPQRSTVWY\"', \
    #     required=False, default='ACDEFGHIKLMNPQRSTVWY')
    # extract_parser.add_argument('-f', '--fastafile', dest='fastafile', type=str, \
    #     help='A fasta file of a single gene to deduce amino acid order of the reference sequence.', \
    #     required=False, default=None)
    # extract_parser.add_argument('-b', '--beginning', dest='beginning', type=int, \
    #     help='An integer indicating the amount to add to the first residue index (1)."+\
    #         " If not specified, nothing(0) will be added.',
    #     required=False, default=0)

    #randsubset script argument parser
    randsubset_parser = subparsers.add_parser('randsubset', \
        description="Selects 'numsel' sequences randomly from an input MSA and writes them to an outputmsa.fasta file!")

    randsubset_parser.add_argument('-i', '--inputmsa', dest='inputmsa', type=str, \
        help="Name of the input Multiple Sequence Alignment file in fasta format.",
        required=True, default=None)

    randsubset_parser.add_argument('-o', '--outputmsa', dest='outputmsa', type=str, \
        help="Name of the output Multiple Sequence Alignment file in fasta format. \
              Default is outputmsa.fasta",
        required=False, default="outputmsa.fasta")

    randsubset_parser.add_argument('-n', '--numsel', dest='numsel', type=int, \
        help='An integer number showing number of randomly selected sequences from the input MSA.',
        required=True, default=None)

        #compare script argument parsing
    diffmap_parser = subparsers.add_parser('diffmap', description="Obtain a difference map between two predictions.")
    diffmap_parser.add_argument('-i', '--inputfile1', dest='inputfile1', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    diffmap_parser.add_argument('--itype', dest='itype', type=str, \
        help='gemme, rhapsody, foldx, evmutation or singleline. Default is gemme.', \
        required=False, default='gemme')
    diffmap_parser.add_argument('-j', '--inputfile2', dest='inputfile2', type=str, \
        help='One of the output files of gemme, rhapsody or evmutation', \
        required=True, default=None)
    diffmap_parser.add_argument('--jtype', dest='jtype', type=str, \
        help='gemme, rhapsody, foldx, evmutation or singleline. Default is gemme.', \
        required=False, default='gemme')
    diffmap_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help='Name of the output file. Default is output.txt', \
        required=False, default='output')
    diffmap_parser.add_argument('--paginate', dest='paginate', type=int, \
        help='An integer value. Default is 0, which means to plot without pagination.',
        required=False, default=0)    

    diffmap_parser.add_argument('-b', '--beginning', dest='beginning', type=int, \
        help='An integer to indicate the first residue index.',
        required=False, default=1)

    diffmap_parser.add_argument('-e', '--end', dest='end', type=int, \
        help='An integer to indicate the final residue index.',
        required=False, default=None)

    diffmap_parser.add_argument('--colormap', dest='colormap', type=str, \
        help='A colormap as defined in matplotlib',
        required=False, default='coolwarm_r')
    # diffmap_parser.add_argument('--ranknorm', dest='ranknorm', type=bool, \
    #     help='A True or False value to apply rank normalization to data matrix',
    #     required=False, default=False)

    diffmap_parser.add_argument('--iscolorbaron', dest='iscolorbaron', type=bool, \
        help='A True or False value to put color bar on the map or not. Default is False.',
        required=False, default=False)
    
    diffmap_parser.add_argument('--offset', dest='offset', type=int, \
        help='An integer value to offset the xlabels for incomplete sequences',
        required=False, default=0)
    
    diffmap_parser.add_argument('-s', '--sequence', dest='sequence', type=str, \
        help='Input sequence file in fasta format', \
        required=False, default=None)

    diffmap_parser.add_argument('--aaorder', dest='aaorder', type=str, \
        help='Amino acid order as a single string. Default is alphabetical: \"ACDEFGHIKLMNPQRSTVWY\"', \
        required=False, default='ACDEFGHIKLMNPQRSTVWY')
    # diffmap_parser.add_argument('--msafile', dest='msafile', type=str, \
    #     help='An aligned multiple sequence alignment file in fasta format.', \
    #     required=False, default=None)
    # diffmap_parser.add_argument('--ssfile', dest='ssfile', type=str, \
    #     help='An secondary structure file in plain text format.', \
    #     required=False, default=None)


    #inputgenerator script argument parser
    inputgenerator_parser = subparsers.add_parser('inputgenerator', \
        description="Generates inputs for polyphen2 and (hopefully) VESPA using a single sequence file in fasta format!")

    inputgenerator_parser.add_argument('-i', '--inputsequence', dest='inputsequence', type=str, \
        help="Name of the input sequence file in fasta format.",
        required=True, default=None)

    inputgenerator_parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, \
        help="Name of the output file. \
              Default is output.txt",
        required=False, default="output.txt")
    
    inputgenerator_parser.add_argument('--otype', dest='otype', type=str, \
        help='polyphen or VESPA. Default is polyphen.', \
        required=False, default='polyphen')
    
    #postprocess script argument parsing; copied as it is from sedy. 
    postprocess_parser  = subparsers.add_parser('postprocess')
    postprocess_parser.add_argument("-i", "--input", \
        help="gemme file in a two columns format.", \
        dest="input", required=True, default=None)
    postprocess_parser.add_argument("--itype", \
        help="It is gemme, generic or esm1b.", \
        dest="itype", required=False, default="gemme")
    postprocess_parser.add_argument("-c", "--column", 
        help="An integer value for the data column to process.",
        dest="column", type=int,  required=False, default=2)
    postprocess_parser.add_argument("--process", \
        help="One of these options: ranksort, 1-ranksort, minmax, 1-minmax, 1-values, zscore, negzscore, attenuate", 
        dest="process", type=str,  required=False, default=None)
    postprocess_parser.add_argument("-o", "--outfile", 
        help="Name of your output file.", 
        dest="outfile", default="outfile.dat")
    args = main_parser.parse_args()

    if args.command == "maps":
        mapsApp(args)
    elif args.command == "plots":
        plotsApp(args)
    elif args.command == "compare":
       compareApp(args)
    elif args.command == "convert":
       convertApp(args)
    elif args.command == "extract":
       extractApp(args)
    elif args.command == "removegaps":
       removeGapsApp(args)
    elif args.command == "riesselman":
       riesselmanApp(args)
    elif args.command == "randsubset":
       randsubsetApp(args)
    elif args.command == "diffmap":
       diffmapApp(args)
    elif args.command == "inputgenerator":
       inputgeneratorApp(args)
    elif args.command == "postprocess":
        postprocessApp(args)
    elif args.command == "-h" or args.command == "--help":
        usage_main()
    else:
        usage_main()
        sys.exit(-1)
    
if __name__ == "__main__":
    main()
