"""
    This script is written to convert FoldX output to GEMME format. 
    It can also be used to convert GEMME standard output, which is in 
    alphabetical one-letter amino acid (aa) order, to another aa order. 
"""
import argparse
import sys
from demust.io import *
from Bio import SeqIO
#import pandas as pd

def convertApp(args):
    if (args.inputfile == None):
        print('Usage: demust convert [-h] [-i INPUTFILE] [--itype GEMME] [-o OUTPUTFILE] [--otype GEMME]')
        print('\nError: missing arguments: Please provide --inputfile, --itype, --outputfile, --otype ')

        sys.exit(-1)


    
    # An amino acid order list from 
    # https://www.embopress.org/doi/full/10.15252/msb.20177908
    aaOrderV2 = 'AVLIMFYWRHKDESTNQGCP'
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


    aaOrderList = list(args.aaorder)
    print("\nWriting the residues in the following order:\n")
    print(aaOrderList)
    localResidueList = None
    if(args.fastafile != None):
        referenceSeq = SeqIO.read(args.fastafile, 'fasta')
        localResidueList = list(referenceSeq.seq)
        #print(residueList)

    if (args.otype.lower()=='gemme'):
        writeGEMMEmatrix(scanningMatrix, args.outputfile, beg=0, end=None, \
                        aaOrder = aaOrderList, \
                        residueList = None,
                        offSet=0)
    elif (args.otype.lower()=='dms'):
        print("@> Protein= {:}, Average={:.2f}, Min={:.2f}, Max={:.2f} ".format(args.inputfile, \
                                                                        np.mean(scanningMatrix.flatten()), \
                                                                        np.amin(scanningMatrix.flatten()), \
                                                                        np.amax(scanningMatrix.flatten())))
        writeDMSformat(scanningMatrix, args.outputfile, residueList = localResidueList,\
                        beg=args.beginning, end=None, aaOrder = aaOrderList, \
                        offSet=0)

    # elif (args.otype.lower()==('rhapsody' or 'rapsody')):
    #     scanningMatrix2 = parseRHAPSODYoutput(args.inputfil2, field=args.field)
        
    # elif (args.otype.lower()=='foldx'):
    #     scanningMatrix2 = parseFOLDXoutput(args.inputfile2, colorThreshhold=7.5, colorCorrect=True)
    else:
        print("\nError: Unknown otype!")
        print("         Output data type can only be gemme")

    print("\nFinished running 'demust convert' app!\n")


