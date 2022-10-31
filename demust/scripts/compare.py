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

def innerLoopV4(allExpTypedList, allCompTypedList, debug):
    """
        This is my inner loop for numba. 
        This function works on the assumption that for each 
        experimental value, there is a non-nan (float) value
        and line number for both of them are same. 
        I had to do this assumption. Otherwise, computation
        of Spearman correlation for multiple mutations can become
        a pain in the ass. 
    """
    dataSet1 = []
    dataSet2 = []
    # allExpTypedList = List()
    # allCompTypedList = List()
    # [allExpTypedList.append(x) for x in allExpLines]
    # [allCompTypedList.append(x) for x in allCompLines]
    i = 0
    for i in range (len(allExpTypedList)):
        mutationE = allExpTypedList[i].split()[0]
        mutationC = allCompTypedList[i].split()[0]

        if(mutationE == mutationC) and (allCompTypedList[i].split()[1] != 'NA'):
            # if(debug):
            #     print(mutationE, line.split()[1], compline.split()[1])
            dataSet1.append((allExpTypedList[i].split()[1]))
            dataSet2.append((allCompTypedList[i].split()[1]))
        if(i%1000==0):
            print(i)

    return dataSet1, dataSet2

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

        # allExpTypedList = List()
        # allCompTypedList = List()
        # [allExpTypedList.append(x) for x in allExpLines]
        # [allCompTypedList.append(x) for x in allCompLines]  
        #dataSet1, dataSet2 = innerLoop(allExpTypedList, allCompTypedList, debug=True)
        dataSet1, dataSet2 = innerLoopV4(allExpLines, allCompLines, debug=True)
        #dataSet1, dataSet2 = innerLoopV2(allExpTypedList, allCompTypedList, debug=True)
        
        # for line in allExpLines:
        #     mutationE = line.split()[0]
        #     for compline in allCompLines:
        #         mutationC = compline.split()[0]
        #         if(mutationE == mutationC) and (compline.split()[1] != 'NA'):
        #             if(debug):
        #                 print(mutationE, line.split()[1], compline.split()[1])
        #             dataSet1.append(float(line.split()[1]))
        #             dataSet2.append(float(compline.split()[1]))
        
        dataSet1 = np.array(dataSet1, dtype=float)
        dataSet2 = np.array(dataSet2, dtype=float)
    else:
        print("@> ERROR: Unknown --itype or --jtype specified.")
        print("@>        It can be gemme or singleline!")
        sys.exit(-1)
    if(args.metric.lower() == 'spearman'):
        #if((args.itype == "gemme") and (args.jtype == "gemme")):
        correlation, pvalue = compareMapsSpearman(dataSet1, dataSet2)
        print("\nSpearman comparison of {} and {}: correlation={:.3f} ; pvalue={:.3E}\n".format(args.inputfile1, args.inputfile2, correlation, pvalue))
