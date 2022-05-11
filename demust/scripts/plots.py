from demust.io import *
import argparse

def getMinMaxData(scanningMatrix, type, outfile, printDetails=False):
    """
        This function gets min and max for each column (namely,
        over all variants at a certain position). 

    Parameters
    ----------
    scanningMatrix: numpy array of arrays
        Data matrix to rank normalize
    type: string
        type can only be min or max.    
    outfile: string
        Prefix for the output file.
    printDetails: a boolean value
        If True, it will print details for debugging (Default: False)
    Returns
    -------
    Nothing
    """


    FILE = open(outfile, "w")
    if(type.lower()=='max'):
        for col in range(len(scanningMatrix[0])):
            X_max = (scanningMatrix.T[col].max())            
            FILE.write("{}\t{}\n".format(col, X_max))
            
            if(printDetails):
                print(scanningMatrix.T[col])
    elif(type.lower()=='min'):
        for col in range(len(scanningMatrix[0])):
            X_min = (scanningMatrix.T[col].min())
            FILE.write("{}\t{}\n".format(col, X_min))

            if(printDetails):
                print(scanningMatrix.T[col])
    else:
        print("ERROR: Unknown type!")
        print("       Type can only be min or max!")
    
    FILE.close()


def plotsApp(args):
    if args.inputfile == None:
        print('Usage: python demust.py [-h] [-i INPUTFILE] [-d GEMME]')
        print('\nError: missing arguments: Please provide --inputfile and/or --datatype')
        sys.exit(-1)

    
    gemmeData = demust.parseGEMMEoutput(args.inputfile, verbose=False)
    if((args.type.lower()=='min') or (args.type.lower()=='max')):
        getMinMaxData(gemmeData, args.type, outfile=args.outputfile, printDetails=False)

    
