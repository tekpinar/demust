import argparse
from scipy import stats
from demust.io import *

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

    print("\nRunning 'demust compare' app...\n")
    gemmeData1 = parseGEMMEoutput(args.inputfile1, verbose=False)
    gemmeData2 = parseGEMMEoutput(args.inputfile2, verbose=False)

    if(args.metric.lower() == 'spearman'):
        correlation, pvalue = compareMapsSpearman(gemmeData1, gemmeData2)
        print("\nSpearman comparison of {} and {}: correlation={:.3f} - pvalue={:.3E}\n".format(args.inputfile1, args.inputfile2, correlation, pvalue))
