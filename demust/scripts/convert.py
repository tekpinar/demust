"""
    This script is written to convert FoldX output to GEMME format. 
    It can also be used to convert GEMME standard output, which is in 
    alphabetical one-letter amino acid (aa) order, to another aa order. 
"""
import argparse
from demust.io import *

def convertApp(args):
    if (args.inputfile == None):
        print('Usage: demust convert [-h] [-i INPUTFILE] [--itype GEMME] [-o OUTPUTFILE] [--otype GEMME]')
        print('\nError: missing arguments: Please provide --inputfile, --itype, --outputfile, --otype ')
        sys.exit(-1)


    ###########
    print("\nRunning 'demust convert' app...\n")
    
    if (args.itype.lower()=='gemme'):
        scanningMatrix = parseGEMMEoutput(args.inputfile, verbose=False)
        
    elif (args.itype.lower()==('rhapsody' or 'rapsody')):
        scanningMatrix = parseRHAPSODYoutput(args.inputfile, field=args.field)
        
    elif (args.itype.lower()=='foldx'):
        scanningMatrix = parseFOLDXoutput(args.inputfile, colorThreshhold=7.5, colorCorrect=True)
    else:
        print("\nError: Unknown itype!")
        print("         Input data types can be gemme, rhapsody or foldx!")

    if (args.otype.lower()=='gemme'):
        writeGEMMEmatrix(scanningMatrix, args.outputfile, beg=0, end=None, \
                        aaOrder = alphabeticalAminoAcidsList, \
                        residueList = None,
                        offSet=0)

    # elif (args.otype.lower()==('rhapsody' or 'rapsody')):
    #     scanningMatrix2 = parseRHAPSODYoutput(args.inputfil2, field=args.field)
        
    # elif (args.otype.lower()=='foldx'):
    #     scanningMatrix2 = parseFOLDXoutput(args.inputfile2, colorThreshhold=7.5, colorCorrect=True)
    else:
        print("\nError: Unknown otype!")
        print("         Output data type can only be gemme")

    print("\nFinished running 'demust convert' app!\n")


