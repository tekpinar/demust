import sys
# from prody import *
import numpy as np
# from matplotlib.pylab import *
# from sedy.postprocess import *
import pandas as pd
#This script was part of sedy package.
from scipy.stats import rankdata,zscore
import matplotlib.pyplot as plt
# from demust.scripts import extract

alphabeticalAminoAcidsList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
#Part of Sedy Package
def rankSortData(dataArray):
    """
        This function ranksorts protein data and converts it
        to values between [0,1.0].
    
    Parameters:
    ----------
    dataArray: numpy array of arrays
               data read by numpy.loadtxt

    Returns:
    -------
    normalizedRankedDataArray: numpy array
    """
    normalizedRankedDataArray = rankdata(dataArray)/float(len(dataArray))

    return (normalizedRankedDataArray)

def minMaxNormalization(data):
    """
        Min-max normalization of a data array.

    Parameters:
    ----------
    data: numpy array
          

    Returns:
    -------
    normalizeddata: numpy array
    """
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def zscore(data):
    """
        Z-score normalization

    Parameters:
    ----------
    data: numpy array
          

    Returns:
    -------
    normalizeddata: numpy array

    """
    return (data - np.mean(data))/np.std(data)

def attenuateEndPoints(data):                                                                                                    
    """                                                                                                                         
        Reduce signal strength of the signal on N and C terminals
        with a sigmoid function.                                                                
    """                                                                                                                         
                                                                                                                                
    debug = 0                                                                  
                                                                               
    n = len(data)                                                              
    if(debug):                                                                 
        #Check array lengths                                                   
        print(n/2)                                                             
        print(len(data))

    if(n%2 == 0):                                                              
        x1 = np.linspace(-0, 10, int(n/2))                                     
        x2 = np.linspace(-10, 0, int(n/2))                                     
    else:                                                                      
        x1 = np.linspace(-0, 10, int((n-1)/2))                                 
        x2 = np.linspace(-10, 0, int((n+1)/2))                                 
                                                                               
    z1 = 1/(1 + np.exp(-2.50*x1))                                              
    z2 = 1.0 - (1/(1 + np.exp(-2.50*x2)))                                      
                                                                               
                                                                               
    z = np.append(z1, z2)
    if(debug):                                                                 
        #Check array lengths                                                   
        print(len(z1))                                                     
        print(len(z2))                                                         
        print(len(z))                                                          
        #Plot the arrays                                                       
        plt.plot(z)                                                            
        plt.plot(data)                                                         
        plt.plot(data*z)                                                       
        plt.xlabel("x")                                                        
        plt.ylabel("Sigmoid(X)")                                               
        plt.show()

    if(len(data) == len(z)):                                                   
        return (data*z)                                                                                                          
    else:                                                                       
        print("ERROR: Can not attenuate N and C terminal data signals!")        
        sys.exit(-1)

def colRescaleWithAverageRaw(inputfile, outputfile):
    """
        Trying to rescale values in each column (namely all mutations for a
        position) with the average. Basically, we are distancing away the 
        values from the average. This function uses the raw data produced 
        by ESCOTT or iGEMME.  
    """
    #Read the inputfile
    singleFile = open(inputfile, "r")
    allLines = singleFile.readlines()
    singleFile.close()

    positionsList = []
    mutationsList = []
    for line in allLines:
        data = line.split()
        if(":" in data[0]) or ("," in data[0]):
            continue
        else:
            positionsList.append(data[0][0:-1])
            mutationsList.append(data[0])

    # insert the list to the set
    listSet = sorted(set(positionsList), key=positionsList.index)
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
    averagesDict = dict(map(lambda i,j : (i,j) , uniqueMutations,averagesList))
    rescaledValuesList = []
    # for mut in uniqueMutations:
    #     tempSum = 0.0
    #     tempCtr = 0
    for line in allLines:
        data = line.split()
        if(":" in data[0]) or ("," in data[0]):
            continue
        else:
            position = (data[0][0:-1])
            rescaledValuesList.append(2*float(data[1]) - float(averagesDict[position]))

    if(len(mutationsList)==len(rescaledValuesList)):
        with open(outputfile, 'w', encoding='utf-8') as datafile:
            for i in range (len(rescaledValuesList)):
                datafile.write("{} {:.4f}\n".format(mutationsList[i], rescaledValuesList[i]))

            print("@> 'demust postprocess' wrote the column data rescaled with the averages to {} succesfully!".format(outputfile))
        return rescaledValuesList
    else:
        print("@> ERROR: Lenghts of the arrays do not match!")
        sys.exit(-1)

def colRescaleWithAverageRanked(inputfile, outputfile):
    """
        Trying to rescale values in each column (namely all mutations for a
        position) with the average. Basically, we are distancing away the 
        values from the average. This function uses the raw data produced 
        by ESCOTT or iGEMME.  
    """
    #Read the inputfile
    singleFile = open(inputfile, "r")
    allLines = singleFile.readlines()
    singleFile.close()

    positionsList = []
    mutationsList = []
    for line in allLines:
        data = line.split()
        if(":" in data[0]) or ("," in data[0]):
            continue
        else:
            positionsList.append(data[0][0:-1])
            mutationsList.append(data[0])

    # insert the list to the set
    listSet = sorted(set(positionsList), key=positionsList.index)
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
    averagesDict = dict(map(lambda i,j : (i,j) , uniqueMutations,averagesList))
    rescaledValuesList = []
    # for mut in uniqueMutations:
    #     tempSum = 0.0
    #     tempCtr = 0
    for line in allLines:
        data = line.split()
        if(":" in data[0]) or ("," in data[0]):
            continue
        else:
            position = (data[0][0:-1])
            tempValue = 2*float(data[1]) - float(averagesDict[position])
            if (tempValue>1.0):
                tempValue = 1.0
            if(tempValue <=0.0):
                tempValue = 0.0

            rescaledValuesList.append(tempValue)

    if(len(mutationsList)==len(rescaledValuesList)):
        with open(outputfile, 'w', encoding='utf-8') as datafile:
            for i in range (len(rescaledValuesList)):
                datafile.write("{} {:.4f}\n".format(mutationsList[i], rescaledValuesList[i]))

            print("@> 'demust postprocess' wrote the column data rescaled with the averages to {} succesfully!".format(outputfile))
        return rescaledValuesList
    else:
        print("@> ERROR: Lenghts of the arrays do not match!")
        sys.exit(-1)
def colRescaleWithAverageRankedMinmax(inputfile, outputfile):
    """
        Trying to rescale values in each column (namely all mutations for a
        position) with the average. Basically, we are distancing away the 
        values from the average. This function uses the ranked data produced 
        by ESCOTT or iGEMME. There is a min-max per column.   
    """
    #Read the inputfile
    singleFile = open(inputfile, "r")
    allLines = singleFile.readlines()
    singleFile.close()

    positionsList = []
    mutationsList = []
    for line in allLines:
        data = line.split()
        if(":" in data[0]) or ("," in data[0]):
            continue
        else:
            positionsList.append(data[0][0:-1])
            mutationsList.append(data[0])

    # insert the list to the set
    listSet = sorted(set(positionsList), key=positionsList.index)
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
    averagesDict = dict(map(lambda i,j : (i,j) , uniqueMutations,averagesList))
    rescaledValuesList = []
    for mut in uniqueMutations:
        tempSum = 0.0
        tempCtr = 0
        tempColumnList = []
        for line in allLines:
            data = line.split()

            if(":" in data[0]) or ("," in data[0]):
                continue
            else:
                position = (data[0][0:-1])
                if(position == mut):
                    tempValue = 2*float(data[1]) - float(averagesDict[position])
                    
                    tempColumnList.append(tempValue)
        resultPerCol = list((np.array(tempColumnList) - np.min(tempColumnList)) / (np.max(tempColumnList) - np.min(tempColumnList)))
        rescaledValuesList.extend(resultPerCol)
        print(resultPerCol)
        # sys.exit(-1)

    if(len(mutationsList)==len(rescaledValuesList)):
        with open(outputfile, 'w', encoding='utf-8') as datafile:
            for i in range (len(rescaledValuesList)):
                datafile.write("{} {:.4f}\n".format(mutationsList[i], rescaledValuesList[i]))

            print("@> 'demust postprocess' wrote the column data rescaled with the averages to {} succesfully!".format(outputfile))
        return rescaledValuesList
    else:
        print("@> ERROR: Lenghts of the arrays do not match!")
        sys.exit(-1)


def postprocessApp(args):

    if(args.itype=='gemme'):
        # #Mostyl, I am using normPred_Combi_singleline as input file and it doesn't have a header.
        # df = pd.read_table(args.input, sep="\s+", header=None)

        # #data = np.genfromtxt(args.input,dtype=None)
        # data = df.to_numpy()

        # # print(data)
        # rawData = data.T[args.column-1]
        # # print(rawData)
        # # processedData = rankdata(rawData.T[args.column - 1])/float(len(rawData.T[args.column - 1]))
        if (args.fasta == None):
            print("@> ERROR: Fasta file is mandatory for gemme format input type!")
            print("@> Please provide a fasta file with -f option!")
            sys.exit(-1)
        else:
            from Bio import SeqIO
            seq = SeqIO.read(args.fasta, "fasta").seq
            # print(seq)
        gemmeDF = pd.read_table(args.input, sep="\s+")
        #matrix = gemmeDF.to_numpy()
        maxValue =np.nanmax(gemmeDF.to_numpy())
        # print(maxValue)
        matrix = gemmeDF.to_numpy(na_value=maxValue)

        rawData = matrix.T.flatten()

        #offset = 0
        # # print(gemmeDF)
        # gemmeDFtrans.columns = gemmeDFtrans.columns.str.upper()
        aaAndPosition = []
        # oldNamesList = gemmeDF.columns.tolist()
        for i in range(len(list(seq))):
            for item in alphabeticalAminoAcidsList:
                aaAndPosition.append(list(seq)[i]+str(i+1+args.offset)+item)

        data = np.array(list(zip(aaAndPosition, rawData)))

    elif(args.itype=='singleline'):
        #Mostyl, I am using normPred_Combi_singleline as input file and it doesn't have a header.
        df = pd.read_table(args.input, sep="\s+", header=None)

        #data = np.genfromtxt(args.input,dtype=None)
        data = df.to_numpy()

        # print(data)
        rawData = data.T[args.column-1]
        # print(rawData)
        # processedData = rankdata(rawData.T[args.column - 1])/float(len(rawData.T[args.column - 1]))
    elif(args.itype=='esm1b'):
        dfESMvariants = pd.read_csv(args.input)
        # dfESMvariants.rename(columns = {'mut_name':'mutant'}, inplace = True)
        rawData = dfESMvariants['esm_score'].to_numpy()

    else:
        print("Unknown input type: {}".format(args.itype))
        sys.exit(-1)

    if(args.process == 'ranksort'):
        processedData = rankSortData(rawData)
    elif(args.process == '1-ranksort'):
        processedData = 1.0 - rankSortData(rawData)
        #print("@> Average pathogenicity for {:}={:.4f}".format(args.input, np.mean(processedData)))
    elif(args.process == 'minmax'):
        processedData = minMaxNormalization(rawData)
    elif(args.process == '1-minmax'):
        processedData = 1.0 - minMaxNormalization(rawData)
    elif(args.process == '1-values'):
        processedData = 1.0 - rawData
    elif(args.process == 'zscore'):
        processedData = zscore(rawData)
    elif(args.process == 'negzscore'):
        processedData = -1.0*zscore(rawData)
    elif(args.process == 'attenuate'):
        processedData = attenuateEndPoints(rawData)
    elif(args.process == 'colrescaleraw'):
        processedData = colRescaleWithAverageRaw(args.input, args.outfile)
    elif(args.process == 'colrescaleranked'):
        processedData = colRescaleWithAverageRanked(args.input, args.outfile)
    elif(args.process == 'colrescalerankedminmax'):
        processedData = colRescaleWithAverageRankedMinmax(args.input, args.outfile)
    elif(args.process == None):
        print("@> No postprocessing applied!")
    else:
        print("@> Unknown postprocessing option!")
        sys.exit(-1)

    print("@> Processing the data started!")

    if(args.itype=='gemme'):
        with open(args.outfile, 'w') as f:
            #f.write("#Resid Value\n")
            for i in range (len(processedData)):
                f.write("{:} {:6.2f}\n".format(data.T[0][i], processedData[i]))
    
    if(args.itype=='singleline'):
        with open(args.outfile, 'w') as f:
            #f.write("#Resid Value\n")
            for i in range (len(processedData)):
                f.write("{:} {:6.2f}\n".format(data.T[0][i], processedData[i]))
    if(args.itype=='esm1b'):
        dfESMvariants['esm_score'] = processedData
        dfESMvariants.to_csv(args.outfile, index=False)

    print("@> Processing the data finished successfully!")
    sys.exit()

