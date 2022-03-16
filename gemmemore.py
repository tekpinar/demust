from ctypes import alignment
from lib2to3.pgen2.token import EQUAL
import os
import sys
import argparse
import numpy as np
import matplotlib.pylab as plt

alphabeticalAminoAcidsList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
def parseGEMMEoutput(inputFile):
    """
        Parse normalized (I don't know how?!) GEMME output files: 
        Independent, Epistatic and Combined
        Return: Data File as a Numpy array, where each col contains 
        deleteriousness score for an amino acid.
    """
    #gemmeDataFile =open(inputFile, 'r')
    allLines = inputFile.readlines()
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

    inputFile.close()
    print(gemmeData[0])
    mutationsData = np.array(gemmeData)
    # mutatedTransposed = mutationsData.T
    # return mutatedTransposed
    return mutationsData

def plotGEMMEmatrix(gemmeData, outFile, offSet=0):
    """
        A test to plot a better GEMME matrix
    """
  
    nres_shown = len(gemmeData[0])
    fig_height=8
    # figure proportions
    fig_width = fig_height/2  # inches
    fig_width *= nres_shown/20
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
  
    ##########################################################################
    # Set plotting parameters
    majorTics = 50
    major_nums_x = np.arange(majorTics, len(gemmeData[0]), majorTics, dtype=int)
    major_nums_x = major_nums_x -1
    major_nums_x = np.insert(major_nums_x, 0, 0)
    # print(major_nums_x)

    minor_nums_x = np.arange(10, len(gemmeData[0]), 10, dtype=int)
    minor_nums_x = minor_nums_x - 1
    minor_nums_x = np.insert(minor_nums_x, 0, 0)
    print(minor_nums_x)

    major_labels_x = major_nums_x + 1 + offSet

    major_nums_y = np.arange(0, 20, 1, dtype=int)
    major_labels_y = alphabeticalAminoAcidsList

    plt.xticks(major_nums_x, major_labels_x, size=28)
    ax.set_xticks(minor_nums_x, minor=True)
    
    plt.yticks(major_nums_y, major_labels_y, size=28)
    ax.set_yticklabels(major_labels_y, ha='left')
    ax.tick_params(axis='y', which='major', pad=30)

    #############################################################################
    plt.imshow(gemmeData, cmap='coolwarm_r', aspect=3.0)

    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.tight_layout()
    plt.savefig(outFile)
    plt.show()
    

    #plt.imsave('output.png', gemmeData)
    plt.close()

def parseRHAPSODYoutput(inputFile):
    """
        Parse RHAPSODY output files

    Parameters
    ----------
    inputFile: string
        Name of the Rhapsody output file 

    Returns
    -------
    Data as a Numpy array, where each col contains 
        deleteriousness score for an amino acid mutation.
    """
    #gemmeDataFile =open(inputFile, 'r')
    allLines = inputFile.readlines()
    headerLine = allLines[0]
    print(headerLine)
    gemmeData = []
    aaColumn = []
    tempList = []
    for line in allLines[1:]:

        data = line.split()
        if(data[0] != '#'):
            print(data[10])
            if (data[10] == 'nan'):
                tempList.append(0.0000000)
            else:
                tempList.append(float(data[10]))
        
    gemmeData.append(np.reshape(tempList, (20, 781)))

    inputFile.close()
    print(gemmeData[0])
    mutationsData = np.array(gemmeData)
    # mutatedTransposed = mutationsData.T
    # return mutatedTransposed
    return mutationsData

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

    parser = argparse.ArgumentParser(description='GEMMEMORE: A Python package to modify and annotate GEMME results')
    parser.add_argument('-i', '--inputfile', dest='inputfile', type=argparse.FileType('r'), help='One of the output files of gemme, rhapsody or evmutation', required=True, default=None)
    parser.add_argument('-d', '--datatype', dest='datatype', type=str, help='gemme, rhapsody or evmutation', required=False, default='gemme')
    parser.add_argument('-o', '--outputfile', dest='outputfile', type=str, help='Name of the output png file', required=False, default='output.png')
   
    args = parser.parse_args()
    if args.inputfile == None:
        print('Usage: python gemmemore.py [-h] [-i INPUTFILE] [-d GEMME]')
        print('\nError: missing arguments: Please provide --inputfile and/or --datatype')
        sys.exit(1)

    if (args.datatype.lower()=='gemme'):
        gemmeData = parseGEMMEoutput(args.inputfile)
        plotGEMMEmatrix(gemmeData, args.outputfile, offSet=413)

    elif (args.datatype.lower()=='rhapsody' or 'rapsody'):
        rhapsodyData = parseRHAPSODYoutput(args.inputfile)
        #plotRHAPSODYmatrix(rhapsodyData, args.outputfile, offSet=413)

    print(rhapsodyData)

    # gemmeAverages = getGEMMEAverages(gemmeData.T)
    
    # #Read PRS data
    # prsData = parsePRSdata(sys.argv[2])

    # plotGEMME2DvsOther(gemmeAverages, prsData, \
    #                     customColor=sys.argv[3], \
    #                     figName=sys.argv[4])
