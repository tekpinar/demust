import os
import sys
from scipy import stats
from demust.io import *
import numpy as np
import matplotlib.pyplot as plt

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

def innerLoopV4(allExpTypedList, allCompTypedList, writeOutput=True, outfile="exp-vs-comp.txt"):
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
    spearmanFile = open(outfile, "w")
    for i in range (len(allExpTypedList)):
        mutationE = allExpTypedList[i].split()[0]
        mutationC = allCompTypedList[i].split()[0]

        if(mutationE == mutationC) and (allCompTypedList[i].split()[1] != 'NA'):
            # if(debug):
            #     print(mutationE, line.split()[1], compline.split()[1])
            dataSet1.append((allExpTypedList[i].split()[1]))
            dataSet2.append((allCompTypedList[i].split()[1]))
            if(writeOutput):
                spearmanFile.write("{} {:.4f} {:.4f}\n".format(mutationE, \
                                    float(allExpTypedList[i].split()[1]), \
                                    float(allCompTypedList[i].split()[1])))

        if(i%1000==0):
            print(i)

    return dataSet1, dataSet2
def calcMutationalState(align_array, mutation):
    """
        Given a multiple sequence alignment file as a matrix (align_array) 
        and a mutation I26F, calculate mutational state. 

        Mutational State (MS) is a parameter indicating if a mutation is from a highly
        conserved state to a lowly represented state or vice versa.
    """

    mutPosition=int(mutation[1:-1])-1
    Istate = str(mutation[0])
    Fstate = str(mutation[-1])
    
    #Calculate mutational state
    iCount = list(align_array[mutPosition]).count(Istate)
    fCount = list(align_array[mutPosition]).count(Fstate)

    MS = float(iCount-fCount)/len(align_array[mutPosition])

    return MS

def calcMutationalStateV2(align_array, mutation):
    """
        Given a multiple sequence alignment file as a matrix (align_array) 
        and a mutation I26F, calculate mutational state. 

        Mutational State (MS) is a parameter indicating if a mutation is from a highly
        conserved state to a lowly represented state or vice versa.
    """

    mutPosition=int(mutation[1:-1])-1
    gapPerColumnCount=list(align_array[mutPosition]).count("-") 
    Istate = str(mutation[0])
    Fstate = str(mutation[-1])
    
    #Calculate mutational state
    iCount = list(align_array[mutPosition]).count(Istate)
    fCount = list(align_array[mutPosition]).count(Fstate)

    MS = float(iCount-fCount)/(len(align_array[mutPosition])-gapPerColumnCount)

    return MS

def innerLoopV5(allExpTypedList, allCompTypedList, msafile=None, writeOutput=True, outfile="exp-vs-comp.txt"):
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

    if(msafile!=None):
        from Bio import AlignIO
        dataArray = []
        alignment = AlignIO.read(msafile, "fasta")
        proteinSize = (len(alignment[0].seq))
        numberOfRows = (len(alignment))
        align_array = np.array([list(rec) for rec in alignment], np.unicode_).transpose()

        for i in range(0, proteinSize):
            #print(i)
            gapPerColumnCount=list(align_array[i]).count("-")
            dataArray.append(1.0 - (gapPerColumnCount/float(numberOfRows)))
    # allExpTypedList = List()
    # allCompTypedList = List()
    # [allExpTypedList.append(x) for x in allExpLines]
    # [allCompTypedList.append(x) for x in allCompLines]
    i = 0
    spearmanFile = open(outfile, "w")
    for i in range (len(allExpTypedList)):
        mutationE = allExpTypedList[i].split()[0]
        mutationC = allCompTypedList[i].split()[0]

        if(mutationE == mutationC) and (allCompTypedList[i].split()[1] != 'NA'):
            if(False):
                print(mutationE, (allExpTypedList[i].split()[1]), (allCompTypedList[i].split()[1]))
            dataSet1.append((allExpTypedList[i].split()[1]))
            dataSet2.append((allCompTypedList[i].split()[1]))
            if(writeOutput):
                if(msafile!=None):
                    multipleMutation = mutationE.split(":")
                    # print(multipleMutation)
                    tempReliability = 0.0
                    tempMutationalState = 0.0
                    for item in multipleMutation:
                        mutPosition=int(item[1:-1])-1

                        # print(mutPosition)
                        #Calculate reliability even for multiple point mutations 
                        gapPerColumnCount=list(align_array[mutPosition]).count("-")
                        valueR = (1.0 - (gapPerColumnCount/float(numberOfRows)))
                        tempReliability = tempReliability + valueR

                        valueMS = calcMutationalState(align_array, item)
                        tempMutationalState = tempMutationalState + valueMS


                    tempReliability = tempReliability / len(multipleMutation)
                    tempMutationalState = tempMutationalState / len(multipleMutation)


                    spearmanFile.write("{} {:.4f} {:.4f} {:.2f} {:.2f}\n".format(mutationE, \
                                        float(allExpTypedList[i].split()[1]), 
                                        float(allCompTypedList[i].split()[1]), tempReliability, tempMutationalState))
                else:
                    spearmanFile.write("{} {:.4f} {:.4f}\n".format(mutationE, \
                                        float(allExpTypedList[i].split()[1]), 
                                        float(allCompTypedList[i].split()[1])))

        if(i%1000==0):
            print(i)

    return dataSet1, dataSet2    

def innerLoopV6(allExpTypedList, allCompTypedList, msafile=None, ssfile=None, writeOutput=True, outfile="exp-vs-comp.txt"):
    """
        This is my inner loop for numba. 
        This function works on the assumption that for each 
        experimental value, there is a non-nan (float) value
        and line number for both of them are same. 
        I had to do this assumption. Otherwise, computation
        of Spearman correlation for multiple mutations can become
        a pain in the ass. 
        The difference of this function to innerLoopV5:
        It calculates mutational state by removing the number of gaps from the 
        columns.

        ssfile = Secondary structure file
                 It is something like this:
                 1,H
                 2,E
                 3,C
                 ...
                 159,E
        The letters are dssp secondary structure letters. 
    """
    dataSet1 = []
    dataSet2 = []

    if(msafile!=None):
        from Bio import AlignIO
        dataArray = []
        alignment = AlignIO.read(msafile, "fasta")
        proteinSize = (len(alignment[0].seq))
        numberOfRows = (len(alignment))
        align_array = np.array([list(rec) for rec in alignment], np.unicode_).transpose()

        for i in range(0, proteinSize):
            #print(i)
            gapPerColumnCount=list(align_array[i]).count("-")
            dataArray.append(1.0 - (gapPerColumnCount/float(numberOfRows)))
    # allExpTypedList = List()
    # allCompTypedList = List()
    # [allExpTypedList.append(x) for x in allExpLines]
    # [allCompTypedList.append(x) for x in allCompLines]
    i = 0
    spearmanFile = open(outfile, "w")
    for i in range (len(allExpTypedList)):
        mutationE = allExpTypedList[i].split()[0]
        mutationC = allCompTypedList[i].split()[0]

        if(mutationE == mutationC) and (allCompTypedList[i].split()[1] != 'NA' and (allCompTypedList[i].split()[1] != 'NaN')):
            if(False):
                print(mutationE, (allExpTypedList[i].split()[1]), (allCompTypedList[i].split()[1]))
            dataSet1.append((allExpTypedList[i].split()[1]))
            dataSet2.append((allCompTypedList[i].split()[1]))
            if(writeOutput):
                if(msafile!=None):
                    multipleMutation = mutationE.split(":")
                    # print(multipleMutation)
                    tempReliability = 0.0
                    tempMutationalState = 0.0
                    for item in multipleMutation:
                        mutPosition=int(item[1:-1])-1

                        # print(mutPosition)
                        #Calculate reliability even for multiple point mutations 
                        gapPerColumnCount=list(align_array[mutPosition]).count("-")
                        valueR = (1.0 - (gapPerColumnCount/float(numberOfRows)))
                        tempReliability = tempReliability + valueR

                        valueMS = calcMutationalStateV2(align_array, item)
                        tempMutationalState = tempMutationalState + valueMS


                    tempReliability = tempReliability / len(multipleMutation)
                    tempMutationalState = tempMutationalState / len(multipleMutation)


                    spearmanFile.write("{} {:.4f} {:.4f} {:.2f} {:.2f}\n".format(mutationE, \
                                        float(allExpTypedList[i].split()[1]), 
                                        float(allCompTypedList[i].split()[1]), tempReliability, tempMutationalState))
                # ssfile is secondary structure file.
                elif(ssfile!=None):
                    #print("We have a secondary structure file!")
                    #Read secondary structure file
                    mydic = {}
                    with open(ssfile) as f:
                        for line in f:
                            (key, val) = line.split(",")
                            mydic[int(key)] = val
                    multipleMutation = mutationE.split(":")
                    #print(mydic)
                    secStructureState=""
                    for item in multipleMutation:
                        mutPosition=int(item[1:-1])

                        #print(mutPosition)
                        #Add secondary structure information even for multiple point mutations
                        #print(mydic[mutPosition])
                        secStructureState +=mydic[mutPosition].strip()
                    #     gapPerColumnCount=list(align_array[mutPosition]).count("-")
                    #     valueR = (1.0 - (gapPerColumnCount/float(numberOfRows)))
                    #     tempReliability = tempReliability + valueR

                    #     valueMS = calcMutationalStateV2(align_array, item)
                    #     tempMutationalState = tempMutationalState + valueMS


                    # tempReliability = tempReliability / len(multipleMutation)
                    # tempMutationalState = tempMutationalState / len(multipleMutation)


                    spearmanFile.write("{} {:.4f} {:.4f} {}\n".format(mutationE, \
                                        float(allExpTypedList[i].split()[1]), 
                                        float(allCompTypedList[i].split()[1]), secStructureState))
                    
                else:
                    spearmanFile.write("{} {:.4f} {:.4f}\n".format(mutationE, \
                                        float(allExpTypedList[i].split()[1]), 
                                        float(allCompTypedList[i].split()[1])))

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
    
        

    debug = True
    print("\nRunning 'demust compare' app...\n")
    if((args.itype == "gemme") and (args.jtype == "gemme")):
        dataSet1 = parseGEMMEoutput(args.inputfile1, verbose=False)
        dataSet2 = parseGEMMEoutput(args.inputfile2, verbose=False)

    # elif((args.itype == "singleline") and (args.jtype == "singleline")):
    #     with open(args.inputfile1, 'r') as file:
    #         allExpLines = file.readlines()

    #     if(debug):
    #         print(allExpLines)

    #     with open(args.inputfile2, 'r') as file:
    #         allCompLines = file.readlines()

    #     if(debug):
    #         print(allCompLines)

    #     # allExpTypedList = List()
    #     # allCompTypedList = List()
    #     # [allExpTypedList.append(x) for x in allExpLines]
    #     # [allCompTypedList.append(x) for x in allCompLines]  
    #     #dataSet1, dataSet2 = innerLoop(allExpTypedList, allCompTypedList, debug=True)
    #     #dataSet1, dataSet2 = innerLoopV4(allExpLines, allCompLines, writeOutput=True, outfile="spearman.dat")
    #     #dataSet1, dataSet2 = innerLoopV5(allExpLines, allCompLines, msafile=args.msafile, writeOutput=True, outfile="exp-vs-comp.txt")
    #     dataSet1, dataSet2 = innerLoopV6(allExpLines, allCompLines, msafile=args.msafile, \
    #                                     ssfile=args.ssfile, writeOutput=True, outfile=args.outputfile+"-exp-vs-comp.txt")
    #     #dataSet1, dataSet2 = innerLoopV2(allExpTypedList, allCompTypedList, debug=True)
        
    #     # for line in allExpLines:
    #     #     mutationE = line.split()[0]
    #     #     for compline in allCompLines:
    #     #         mutationC = compline.split()[0]
    #     #         if(mutationE == mutationC) and (compline.split()[1] != 'NA'):
    #     #             if(debug):
    #     #                 print(mutationE, line.split()[1], compline.split()[1])
    #     #             dataSet1.append(float(line.split()[1]))
    #     #             dataSet2.append(float(compline.split()[1]))
    #     print(len(dataSet1))
    #     dataSet1 = np.array(dataSet1, dtype=float)
    #     dataSet2 = np.array(dataSet2, dtype=float)
    elif((args.itype == "singleline") and (args.jtype == "singleline")):
        import pandas as pd
        # with open(args.inputfile1, 'r') as file:
        #     allExpLines = file.readlines()
        allExpDF = pd.read_csv(args.inputfile1, sep='\s+', header=None)
        allExpDF.columns = ['variant', 'dms_score']

        if(debug):
            print(allExpDF)

        # with open(args.inputfile2, 'r') as file:
        #     allCompLines = file.readlines()
        allCompDF = pd.read_csv(args.inputfile2, sep='\s+', header=None)
        allCompDF.columns = ['variant', 'comp_score'] 


        if(debug):
            print(allCompDF)

        mergedDF = pd.merge(allExpDF, allCompDF, on='variant')
                # ssfile is secondary structure file.
        if(args.ssfile!=None):
            #print("We have a secondary structure file!")
            #Read secondary structure file
            mydic = {}
            with open(args.ssfile) as f:
                for line in f:
                    (key, val) = line.split(",")
                    mydic[int(key)] = val

            mergedDF['ss'] = ""
            for index, row in mergedDF.iterrows():
                if(':' in row['variant']):
                    multipleMutation = row['variant'].split(":")
                    #print(mydic)
                    secStructureState=""
                    for item in multipleMutation:
                        mutPosition=int(item[1:-1])

                        #print(mutPosition)
                        #Add secondary structure information even for multiple point mutations
                        #print(mydic[mutPosition])
                        secStructureState +=mydic[mutPosition].strip()
                    mergedDF.at[index, 'ss'] = secStructureState
                else:
                    mergedDF.at[index, 'ss'] = mydic[int(row['variant'][1:-1])].strip()

        # allExpTypedList = List()
        # allCompTypedList = List()
        # [allExpTypedList.append(x) for x in allExpLines]
        # [allCompTypedList.append(x) for x in allCompLines]  
        #dataSet1, dataSet2 = innerLoop(allExpTypedList, allCompTypedList, debug=True)
        #dataSet1, dataSet2 = innerLoopV4(allExpLines, allCompLines, writeOutput=True, outfile="spearman.dat")
        #dataSet1, dataSet2 = innerLoopV5(allExpLines, allCompLines, msafile=args.msafile, writeOutput=True, outfile="exp-vs-comp.txt")
        # dataSet1, dataSet2 = innerLoopV6(allExpLines, allCompLines, msafile=args.msafile, \
        #                                 ssfile=args.ssfile, writeOutput=True, outfile=args.outputfile+"-exp-vs-comp.txt")
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
        
        dataSet1 = mergedDF['dms_score'].to_numpy()
        dataSet2 = mergedDF['comp_score'].to_numpy()
        mergedDF.to_csv(args.outputfile+"-exp-vs-comp.txt", float_format='%.8f', header=None, index=None, sep=' ')
        print(len(dataSet1))
    else:
        print("@> ERROR: Unknown --itype or --jtype specified.")
        print("@>        It can be gemme or singleline!")
        sys.exit(-1)
    if(args.metric.lower() == 'spearman'):
        #if((args.itype == "gemme") and (args.jtype == "gemme")):
        correlation, pvalue = compareMapsSpearman(dataSet1, dataSet2)
        fig = plt.figure()
        plt.rcParams.update({'font.size': 16})
        plt.scatter(dataSet1, dataSet2, s=10, color="tab:orange")
        plt.title("Spearman Correlation: {:.3f}".format(correlation))
        plt.xlabel(args.xlabel)
        plt.ylabel(args.ylabel)
        plt.axis('square')
        #plt.set_aspect('equal', adjustable='box')
        plt.savefig(args.outputfile+"-spearman2d.png")

        print("\nSpearman comparison of {} and {}: correlation={:.3f} ; pvalue={:.3E}\n".format(args.inputfile1, args.inputfile2, correlation, pvalue))
