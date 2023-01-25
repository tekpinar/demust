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
alphabeticalAminoAcidsList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def getGEMMEmatrixRow(inputfile, mutationsTo="a"):
    """
        This method gets full single point mutation scanning matrix produced by 
        GEMME and extracts all mutations to a sinle residue such as alanine.
    """
    import pandas as pd
    print("Here I am!")
    df = pd.read_table(inputfile, sep=" ")
    df_transposed = df.T
    print(df_transposed[mutationsTo])
    df_transposed.fillna(0, inplace=True)
    return (df_transposed[mutationsTo].to_list())

def getSingleLineAverages(inputfile, outputfile):
    """
        This method gets a mutations file in singleline format and return the average
        for a particular residue.

        Singleline format is as follows:
        L1A 0.0
        L1C 0.3
        C2A 0.9
        .......
        The output may not be uniform because there may be ten mutations tested for residue
        1 but only three mutations tested for residue 2. It is what it is, unfortunately.
    """
    #Read the inputfile
    singleFile = open(inputfile, "r")
    allLines = singleFile.readlines()
    singleFile.close()

    mutationsList = []
    for line in allLines:
        data = line.split()
        if(":" in data[0]) or ("," in data[0]):
            continue
        else:
            mutationsList.append(data[0][0:-1])

    # insert the list to the set
    listSet = sorted(set(mutationsList), key=mutationsList.index)
    # convert the set to the list
    uniqueMutations = (list(listSet))

    #print(uniqueMutations)

    averagesList = []
    for mut in uniqueMutations:
        tempSum = 0.0
        tempCtr = 0
        for line in allLines:
            data = line.split()
            if(":" in data[0]) or ("," in data[0]):
                continue
            else:
                if (mut == (data[0][0:-1])):
                    tempSum = tempSum + float(data[1])
                    tempCtr += 1

        averagesList.append(tempSum/tempCtr)

    #print(averagesList)

    if(len(uniqueMutations)==len(averagesList)):
        with open(outputfile, 'w', encoding='utf-8') as datafile:
            for i in range (len(uniqueMutations)):
                datafile.write("{} {:.4f}\n".format(uniqueMutations[i][1:], averagesList[i]))

            print("@> 'demust extract' wrote the data to {} succesfully!".format(outputfile))
    else:
        print("@> ERROR: Lenghts of the arrays do not match!")
        sys.exit(-1)

def getSingleLineSingleMutations(inputfile, outputfile, mutationsTo="a"):
    """
        This method gets a mutation result for a single amino acid in singleline format.
        It is useful to obtain results like alanine scanning experiments from a mixed dataset. 
    """
    #Read the inputfile
    singleFile = open(inputfile, "r")
    allLines = singleFile.readlines()
    singleFile.close()

    mutationsList = []
    for line in allLines:
        data = line.split()
        if(":" in data[0]) or ("," in data[0]):
            continue
        else:
            mutationsList.append(data[0][0:-1])

    # insert the list to the set
    listSet = sorted(set(mutationsList), key=mutationsList.index)
    # convert the set to the list
    uniqueMutations = (list(listSet))

    #print(uniqueMutations)

    mutationList = []
    residueIDList = []

    for line in allLines:
        data = line.split()
        if(":" in data[0]) or ("," in data[0]):
            continue
        else:
            if (mutationsTo.lower() == (data[0][-1]).lower()):
                mutationList.append(float(data[1]))
                residueIDList.append(int(data[0][1:-1]))

    #return(mutationList)

    if(len(mutationList)==len(residueIDList)):
        with open(outputfile, 'w', encoding='utf-8') as datafile:
            for i in range (len(mutationList)):
                datafile.write("{} {:.4f}\n".format(residueIDList[i], mutationList[i]))

            print("@> 'demust extract' wrote the data to {} succesfully!".format(outputfile))
    else:
        print("@> ERROR: Lenghts of the arrays do not match!")
        sys.exit(-1)

def extractApp(args):
    if (args.inputfile == None):
        print('Usage: demust extract [-h] [-i INPUTFILE] [--itype GEMME] [-o OUTPUTFILE] [--otype GEMME]')
        print('\nError: missing arguments: Please provide --inputfile, --itype, --outputfile, --otype ')

        sys.exit(-1)

    probableMutationsList = []
    for item in alphabeticalAminoAcidsList:
        # print(item)
        probableMutationsList.append(("mut"+item).lower())

    #print(probableMutationsList)

    if (args.itype.lower()=='gemme'):
        scanningMatrix = parseGEMMEoutput(args.inputfile, verbose=False)
        if (args.otype.lower()=='average'):
            dataArray = getGEMMEAverages(scanningMatrix.T)
        
        elif (args.otype.lower() in probableMutationsList):
            dataArray = getGEMMEmatrixRow(args.inputfile, mutationsTo=args.otype.lower()[-1])

        with open(args.outputfile, 'w', encoding='utf-8') as datafile:
            for i in range (len(dataArray)):
                datafile.write("{} {}\n".format(i+1, dataArray[i]))

            print("@> 'demust extract' wrote the data to {} succesfully!".format(args.outputfile))
    
    elif (args.itype.lower()=='singleline'):
        if (args.otype.lower()=='average'):
            getSingleLineAverages(args.inputfile, args.outputfile)
        elif (args.otype.lower() in probableMutationsList):
            getSingleLineSingleMutations(args.inputfile, args.outputfile, mutationsTo=args.otype.lower()[-1])
            # with open(args.outputfile, 'w', encoding='utf-8') as datafile:
            #     for i in range (len(dataArray)):
            #         datafile.write("{} {}\n".format(i+1, dataArray[i]))

            #     print("@> 'demust extract' wrote the data to {} succesfully!".format(args.outputfile))
    
    else:
        print("@ ERROR: Unknown input type: It can only be gemme or singleline!")
