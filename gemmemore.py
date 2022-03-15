import os
import sys
import numpy as np
import matplotlib.pylab as plt

def parseGEMMEoutput(inputFile):
    """
        Parse normalized (I don't know how?!) GEMME output files: 
        Independent, Epistatic and Combined
        Return: Data File as a Numpy array, where each col contains 
        deleteriousness score for an amino acid.
    """
    gemmeDataFile =open(inputFile, 'r')
    allLines = gemmeDataFile.readlines()
    headerLine = allLines[0]
    print(headerLine)
    gemmeData = []
    aaColumn = []
    for line in allLines[1:]:
        tempList = []
        data = line.split()
        for item in data:
            if item[0] == "\"":
                aaColumn.append(item)
            elif (item == 'NA'):
                tempList.append(0.0000000)
            else:
                tempList.append(float(item))
        gemmeData.append(tempList)

    gemmeDataFile.close()
    print(gemmeData[0])
    mutationsData = np.array(gemmeData)
    # mutatedTransposed = mutationsData.T
    # return mutatedTransposed
    return mutationsData

def plotGEMMEmatrix(gemmeData):
    """
        A test to plot a better GEMME matrix
    """
    fig = plt.figure(figsize=(17, 7))
  
    #print(len(gemmeData[0]))

    # import seaborn as snb
    # snb.heatmap(gemmeData)
    plt.imshow(gemmeData, cmap='bwr')

    plt.tight_layout()
    plt.show()
    plt.close()

def getGEMMEAverages(gemmeData):
    """
    Get GEMME columnwise averages from the normalized (?) data.

    These averages can be compared with 2D dynamical quantities like 
    PRS, MSF, stiffness, DCI etc.

    return A list which has the length of N, where N is the number of residues. 
    """
    gemmeAverages = []
    for column in gemmeData:

        gemmeAverages.append(np.mean(column))
    
    return gemmeAverages

def parsePRSdata(prsDataFile):
    """
        Parses Perturbation Response Scanning (PRS) data file obtained from 
        http://enm.pitt.edu/. The data file is csv format. 
        Returns a numpy array with two columns.
        The first column contains residue IDs obtained from the pdb file.  
        The second column contains the PRS values
    """
    prsData = np.loadtxt(prsDataFile, delimiter=',', comments='\"')

    return prsData.transpose()

def plotGEMME2DvsOther(gemme2DData, other2DData, \
                        customColor="purple", \
                        figName="outfile.png", saveFig=True):
    """
        Saving 2D plots of comparison with GEMME 2D data (like averages)
        with other 2D data like PRS, MSF etc.
        Return Nothing
    """
    #Plot the figure
    plt.rcParams.update({'font.size': 20})
    fig, ax1 = plt.subplots()

    fig.set_size_inches(17, 5)
    ax1.set_xlabel("Residue IDs")
    ax1.set_ylabel(r"$(GEMME Score)^2$")

    ax1.plot(np.arange(1, len(gemme2DData)+1), \
        np.multiply(gemme2DData, gemme2DData),\
        linestyle='-', color="black", marker='o',\
        label="GEMME: Deleterious Mutations")

    ax1.legend(loc='upper left')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel("PRS")
    ax2.plot(np.arange(1, len(gemme2DData)+1), other2DData[1], \
            linestyle='--', marker='*', color=customColor, label="PRS")
    
    plt.legend()
    fig.tight_layout()
    if(saveFig==True):
        plt.savefig(figName)
    plt.show()

if (__name__ == '__main__'):

    gemmeData = parseGEMMEoutput(sys.argv[1])

    plotGEMMEmatrix(gemmeData)

    gemmeAverages = getGEMMEAverages(gemmeData.T)
    
    #Read PRS data
    prsData = parsePRSdata(sys.argv[2])

    plotGEMME2DvsOther(gemmeAverages, prsData, \
                        customColor=sys.argv[3], \
                        figName=sys.argv[4])
