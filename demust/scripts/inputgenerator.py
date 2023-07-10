"""
    This script is written to extract data from GEMME output. 
    It can extract column averages. Here, a column is basically 
    residue index. 
    Moreover, it can extract effects of mutations to 
    only a single amino acids like alanine. This means to extract a row
    from GEMME full matrix output. You have to specify --otype mutA as 
    output type. 
"""
import sys
import argparse
from demust.io import *
from Bio import SeqIO
alphabeticalAminoAcidsList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# def getGEMMEmatrixRow(inputfile, mutationsTo="a"):
#     """
#         This method gets full single point mutation scanning matrix produced by 
#         GEMME and extracts all mutations to a sinle residue such as alanine.
#     """
#     import pandas as pd
#     #print("Here I am!")
#     df = pd.read_table(inputfile, sep=" ")
#     df_transposed = df.T
#     print(df_transposed[mutationsTo])
#     df_transposed.fillna(0, inplace=True)
#     return (df_transposed[mutationsTo].to_list())

def inputgeneratorApp(args):
    if (args.inputsequence == None):
        print('Usage: demust extract [-h] [-i INPUTSEQUENCE] [-o OUTPUTFILE] [--otype OTYPE]')
        print('\nError: missing arguments: Please provide --inputsequence, --itype, --outputfile, --otype ')

        sys.exit(-1)
    print("@> Running inputgenerator app...\n\n")
    if (args.otype.lower()=='polyphen'):
        mySequence = SeqIO.read(args.inputsequence, "fasta")
        print("Preparing polyphen2 input file from fasta file!")
        #Converts all lowercase letters into uppercase
        #because polyphen gives error with lowercase letters. 

        mySequence.seq = mySequence.seq.upper()
        #print(mySequence.name)
        #print(mySequence.seq)
        #mySequenceList = list(mySequence.seq)
        with open(args.outputfile, 'w') as my_file:
            for i in range (len(mySequence)):
                for aa in alphabeticalAminoAcidsList:
                    if (aa != mySequence[i].upper()):
                        my_file.write("{} {} {} {}\n".format(mySequence.id, \
                                                    (i+1), mySequence[i], aa))
        print("Created {} polyphen2 input file successfully!".\
              format(args.outputfile))
    elif (args.otype.lower()=='vespa'):
        print("Well, I think this part is not implemented yet!")
    else:
        print("@ ERROR: Unknown output type: It can only be polyphen or vespa!")
