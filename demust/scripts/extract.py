"""
    This script is written to extract data from GEMME output. 
    It can extract column averages. Here, a column is basically 
    residue index. 
    Moreover, it can extract effects of mutations to 
    only a single amino acids like alanine. This means to extract a row
    from GEMME full matrix output. You have to specify --otype mutA as 
    output type. 
"""
import argparse
from demust.io import *
alphabeticalAminoAcidsList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def getGEMMEmatrixRow(inputfile, mutationsTo="a"):
    import pandas as pd
    print("Here I am!")
    df = pd.read_table(inputfile, sep=" ")
    df_transposed = df.T
    print(df_transposed[mutationsTo])
    df_transposed.fillna(0, inplace=True)
    return (df_transposed[mutationsTo].to_list())

def extractApp(args):
    if (args.inputfile == None):
        print('Usage: demust extract [-h] [-i INPUTFILE] [--itype GEMME] [-o OUTPUTFILE] [--otype GEMME]')
        print('\nError: missing arguments: Please provide --inputfile, --itype, --outputfile, --otype ')

        sys.exit(-1)

    probableMutationsList = []

    for item in alphabeticalAminoAcidsList:
        print(item)
        probableMutationsList.append(("mut"+item).lower())

    print(probableMutationsList)

    if (args.itype.lower()=='gemme'):
        scanningMatrix = parseGEMMEoutput(args.inputfile, verbose=False)
        if (args.otype.lower()=='average'):
            dataArray = getGEMMEAverages(scanningMatrix.T)
        
        elif (args.otype.lower() in probableMutationsList):
            dataArray = getGEMMEmatrixRow(args.inputfile, mutationsTo=args.otype.lower()[-1])

        with open(args.outputfile, 'w', encoding='utf-8') as datafile:
            for i in range (len(dataArray)):
                datafile.write("{} {}\n".format(i+1, dataArray[i]))

            print("@> demust extract wrote the data to {} succesfully!".format(args.outputfile))