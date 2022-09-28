import os
import sys
from scipy import stats
from demust.io import *
import numpy as np

def compareMapsSpearman(scanningMatrix1, scanningMatrix2):
    """
        A function to compare similarity of two deep mutational scanning matrices. 
  
    Parameters
    ----------
    scanningMatrix1: numpy array of arrays
        The first data matrix.

    scanningMatrix2: numpy array of arrays
        The second data matrix.

    outFile: string
        Name of the output png image
    
    Returns
    -------
    correlation: float or ndarray (2-D square)
        Spearman correlation coefficient 

    pvalue: float
        The p-value for a hypothesis test whose null hypotheisis is that two 
        sets of data are uncorrelated.

    """
    # print(scanningMatrix1)
    # print(scanningMatrix2)
    correlation, pvalue = stats.spearmanr(scanningMatrix1.flatten(), scanningMatrix2.flatten())
    return(correlation, pvalue)

def compareApp(args):
    if (args.inputfile1 == None or args.inputfile2 == None):
        print('Usage: demust compare [-h] [-i INPUTFILE1] [--itype GEMME] [-j INPUTFILE2] [--jtype GEMME]')
        print('\nError: missing arguments: Please provide --inputfile1, --itype, --inputfile2, --jtype ')
        sys.exit(-1)

    if (os.path.isfile(args.inputfile1)==False  or os.path.isfile(args.inputfile2)==False):
        print("ERROR: {} or {} files do not exist in the folder!".format(args.inputfile1, args.inputfile2))
        print("       How can I calculate Spearman without them?")
        sys.exit(-1)
    
        

    debug = False
    print("\nRunning 'demust compare' app...\n")
    if((args.itype == "gemme") and (args.jtype == "gemme")):
        dataSet1 = parseGEMMEoutput(args.inputfile1, verbose=False)
        dataSet2 = parseGEMMEoutput(args.inputfile2, verbose=False)

    elif((args.itype == "singleline") and (args.jtype == "singleline")):
        with open(args.inputfile1, 'r') as file:
            allExpLines = file.readlines()

        if(debug):
            print(allExpLines)

        with open(args.inputfile2, 'r') as file:
            allCompLines = file.readlines()

        if(debug):
            print(allCompLines)

        dataSet1 = []
        dataSet2 = []
        for line in allExpLines:
            mutationE = line.split()[0]
            for compline in allCompLines:
                mutationC = compline.split()[0]
                if(mutationE == mutationC):
                    if(1):
                        print(mutationE, line.split()[1], compline.split()[1])
                    dataSet1.append(float(line.split()[1]))
                    dataSet2.append(float(compline.split()[1]))
        
        dataSet1 = np.array(dataSet1)
        dataSet2 = np.array(dataSet2)
    else:
        print("@> ERROR: Unknown --itype or --jtype specified.")
        print("@>        It can be gemme or singleline!")
        sys.exit(-1)
    if(args.metric.lower() == 'spearman'):
        #if((args.itype == "gemme") and (args.jtype == "gemme")):
        correlation, pvalue = compareMapsSpearman(dataSet1, dataSet2)
        print("\nSpearman comparison of {} and {}: correlation={:.3f} - pvalue={:.3E}\n".format(args.inputfile1, args.inputfile2, correlation, pvalue))
