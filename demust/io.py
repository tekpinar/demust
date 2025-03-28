import os
import string
import sys
import argparse
import numpy as np
import matplotlib.pylab as plt
import matplotlib
from matplotlib.gridspec import GridSpec
from Bio import SeqIO
import pandas as pd
from scipy.stats import rankdata
from collections import OrderedDict


alphabeticalAminoAcidsList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aaSingleNucleotideNeighborDict = OrderedDict()
aaSingleNucleotideNeighborDict['A'] = ['V', 'D', 'E', 'G', 'T', 'P', 'S']
aaSingleNucleotideNeighborDict['C'] = ['F', 'S', 'Y', 'W', 'R', 'G']
aaSingleNucleotideNeighborDict['D'] = ['V', 'A', 'G', 'E', 'N', 'H', 'Y']
aaSingleNucleotideNeighborDict['E'] = ['V', 'A', 'G', 'D', 'K', 'Q']
aaSingleNucleotideNeighborDict['F'] = ['S', 'C', 'Y', 'L', 'I', 'V']
aaSingleNucleotideNeighborDict['G'] = ['V', 'A', 'D', 'E', 'R', 'S', 'W', 'C']
aaSingleNucleotideNeighborDict['H'] = ['L', 'P', 'Q', 'R', 'Y', 'N', 'D']
aaSingleNucleotideNeighborDict['I'] = ['F', 'L', 'V', 'M', 'T', 'N', 'S']
aaSingleNucleotideNeighborDict['K'] = ['N', 'T', 'M', 'R', 'Q', 'E']
aaSingleNucleotideNeighborDict['L'] = ['F', 'I', 'V', 'M', 'P', 'H', 'Q', 'R', 'S', 'W']
aaSingleNucleotideNeighborDict['M'] = ['T', 'K', 'R', 'I', 'L', 'V']
aaSingleNucleotideNeighborDict['N'] = ['I', 'T', 'S', 'K', 'H', 'D', 'Y']
aaSingleNucleotideNeighborDict['P'] = ['L', 'H', 'Q', 'R', 'S', 'T', 'A']
aaSingleNucleotideNeighborDict['Q'] = ['L', 'P', 'R', 'H', 'K', 'E']
aaSingleNucleotideNeighborDict['R'] = ['L', 'P', 'Q', 'H', 'C', 'W', 'S', 'G', 'M', 'T', 'K']
aaSingleNucleotideNeighborDict['S'] = ['F', 'L', 'Y', 'C', 'W', 'P', 'T', 'A', 'I', 'N', 'R', 'G']
aaSingleNucleotideNeighborDict['T'] = ['I', 'M', 'N', 'K', 'S', 'R', 'A', 'P']
aaSingleNucleotideNeighborDict['V'] = ['A', 'D', 'E', 'G', 'F', 'L', 'I', 'M']
aaSingleNucleotideNeighborDict['W'] = ['L', 'S', 'C', 'R', 'G']
aaSingleNucleotideNeighborDict['Y'] = ['F', 'S', 'C', 'H', 'N', 'D']

def parseDNAsequence(inputDNAsequenceFile, inputProteinSequenceFile):
    """
        Parse DNA sequence of the protein, given in fasta format.
    """
    myDNAsequence = SeqIO.read(inputDNAsequenceFile, 'fasta')

    myDNAsequenceList = list(myDNAsequence.seq.upper())
    myProteinSequence = SeqIO.read(inputProteinSequenceFile, 'fasta')
    myProteinSequenceList = list(myProteinSequence.seq.upper())

    #print(myProteinSequenceList)
    # print(len(myDNAsequenceList))
    prot2DNAdictionary = OrderedDict()
    if 3*len(myProteinSequenceList) != len(myDNAsequenceList):
        print("Lengths of protein sequence and DNA sequence do not match!")
        sys.exit(-1)
    else:
        for i in range(0, len(myProteinSequenceList)):
            # print(i)
            prot2DNAdictionary[i] = myDNAsequenceList[3*i]+myDNAsequenceList[3*i+1]+myDNAsequenceList[3*i+2]

    return prot2DNAdictionary

DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}


def getPossible1NucleotideNeighbors(codon):
    dnaNucleotides = ['A', 'T', 'G', 'C']

    #Find a set that doesnt contain codon[i], i=1,3
    #Do this 3 times for each position separately
    changeCodonList = []

    diffPosition1 = list(set(dnaNucleotides) - set(list(codon[0])))
    for i in diffPosition1:
        changedCodon = i+codon[1]+codon[2]
        changeCodonList.append(changedCodon)
    
    diffPosition2 = list(set(dnaNucleotides) - set(list(codon[1])))
    for j in diffPosition2:
        changedCodon = codon[0]+j+codon[2]
        changeCodonList.append(changedCodon)

    diffPosition3 = list(set(dnaNucleotides) - set(list(codon[2])))
    for k in diffPosition3:
        changedCodon = codon[0]+codon[1]+k
        changeCodonList.append(changedCodon)      
        
    changedAAList = []
    for item in changeCodonList:
        changedAAList.append(DNA_Codons[item])
        # print(item)

    #Remove stop codons because we are looking for only possible missense mutations here!
    if "_" in changedAAList:
        changedAAList.remove("_")

    #Remove duplicate elements. 
    return list(set(changedAAList))

def parseGEMMEoutput(inputFile, verbose):
    """
        Parse normalized (I don't know how?!) GEMME output files: 
        Independent, Epistatic and Combined
        Return: Data File as a Numpy array, where each col contains 
        deleteriousness score for an amino acid.
        verbose is a boolean value
    """
    gemmeDataFile =open(inputFile, 'r')
    allLines = gemmeDataFile.readlines()
    gemmeDataFile.close()

    headerLine = allLines[0]
    if(verbose):
        print(headerLine)
    matrixData = []
    aaColumn = []
    for line in allLines[1:]:
        tempList = []
        data = line.split()
        aaColumn.append(data[0])
        for item in data[1:]:
            if (item == 'NA'):
                tempList.append(0.0000000)
            else:
                tempList.append(float(item))
        matrixData.append(tempList)
    mutationsData = np.array(matrixData)
    # mutatedTransposed = mutationsData.T
    # return mutatedTransposed
    return mutationsData

def plotGEMMEmatrix(scanningMatrix, outFile, beg, end, \
                    colorMap = 'turbo_r', \
                    offSet=0, pixelType='square',
                    aaOrder="ACDEFGHIKLMNPQRSTVWY", \
                    sequence=None,\
                    interactive=False,\
                    isColorBarOn=False,\
                    onlyDNAaccessible=None,\
                    climmin=None,\
                    climmax=None):
    """
        A function to plot deep mutational scanning matrices. 
  
    Parameters
    ----------
    scanningMatrix: numpy array of arrays
        Data matrix to plot

    outFile: string
        Name of the output image without file extension.
        Default file extension (png) is added by the program
    
    beg: int
        The first residue to use. It is used to select a subrange 
        of amino acids. It starts from 1.
    
    end: int
        The last residue to use. It is used to select a subrange 
        of amino acids.

    colorMap: matplotlib cmap
        Any colormap existing in matplotlib.
        Default is coolwarm. 

    offSet: int
        It is used to match the residue IDs in your PDB
        file with 0 based indices read from scanningMatrix matrix

    pixelType: string
        It can only have 'square' or 'rectangle' values.
        It is a matter of taste but I added it as an option.
        Default is 'square'

    sequence: string
        A fasta file of one letter amino acid codes from N terminal to C terminal. 

    interactive: bool
        If True, it will plot the map interactively. 

    isColorBarOn: bool
        If True, it will show a colorbar to show the numerical scale of
        the colors. Default is False. 

    onlyDNAaccessible: string
        If it is not None, you must provide a DNA sequence of the protein. 
        It will calculate mutations that are accessible from 'only' one 
        nucleotide change. The other locations will be shown as white. 

    climmin: float
        It is used to put a limit to minimum value of the colorbar. 

    climmax: float
        It is used to put an upper limit to maximum value of the colorbar. 

    Returns
    -------
    Nothing

    """

    #We subtract 1 from beg bc matrix indices starts from 0
    if(end == None):
        end = len(scanningMatrix[0])
        
    print("Beginning: "+str(beg))
    print("End      : "+str(end))
    #print(scanningMatrix)
    #np.reshape(scanningMatrix.flatten(), (20, 286))
    #print(len(scanningMatrix[0]))

    subMatrix = scanningMatrix[:, (beg-1):end]
    #print(subMatrix)

    ##########################################################################
    # Set plotting parameters
    nres_shown = len(subMatrix[0])
    print(nres_shown)
    fig_height=8
    # figure proportions
    fig_width = fig_height/2  # inches
    fig_width *= nres_shown/20
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    if (nres_shown >150):
        majorTics = 50
    else:
        majorTics = 20


    major_nums_x = np.arange(majorTics, len(subMatrix[0]), majorTics, dtype=int)
    #major_nums_x = major_nums_x -1
    major_nums_x = np.insert(major_nums_x, 0, 0)
    #print(major_nums_x)

    minor_nums_x = np.arange(10, len(subMatrix[0]), 10, dtype=int)
    #minor_nums_x = minor_nums_x - 1
    minor_nums_x = np.insert(minor_nums_x, 0, 0)
    #print(minor_nums_x)

    major_labels_x = major_nums_x + 1 + offSet

    major_nums_y = np.arange(0, 20, 1, dtype=int)
    major_labels_y = list(aaOrder)

    plt.xticks(major_nums_x, major_labels_x, size=28)
    ax.set_xticks(minor_nums_x, minor=True)
    
    plt.yticks(major_nums_y, major_labels_y, size=16)
    ax.set_yticklabels(major_labels_y, ha='left')
    ax.tick_params(axis='y', which='major', pad=30)

    # print(subMatrix)
    #############################################################################
    onlyDNAaccessibleList = []
    current_cmap = matplotlib.cm.get_cmap()
    current_cmap.set_bad(color='white')
    if(sequence!=None):
        mySeqFile = SeqIO.read(sequence, 'fasta')

        # onylDNAaccessible = True
        if(onlyDNAaccessible==True):
            mySeqFile = SeqIO.read(sequence, 'fasta')
            # print("I am plotting only DNA accessible mutations.")
            #print(mySeqFile)
            #Convert aaOrder to a list.
            aaOrderList = list(aaOrder)
            for i in range (len(subMatrix[0])):
                j = beg-1+i
                # print(i, aaOrderList.index(sequence[i]))
                notInTheList = list(set(alphabeticalAminoAcidsList) - set(aaSingleNucleotideNeighborDict[mySeqFile.seq[j].upper()]))
                notInTheList.remove(mySeqFile.seq[j].upper())
                for item in notInTheList:
                    # plt.scatter(i, aaOrderList.index(item), s=75, c='white', marker='x')
                    #plt.scatter(i, aaOrderList.index(item), s=200, c='white', marker='s')
                    onlyDNAaccessibleList.append([i, aaOrderList.index(item)])
            for item in onlyDNAaccessibleList:
                # print(item)
                subMatrix[item[1]][item[0]] = np.nan
    # print(subMatrix)
    if(pixelType=='square'):
        #For plotting square pixels
        img = plt.imshow(subMatrix, cmap=colorMap)
    elif(pixelType=='rectangle'):
        #For plotting rectangular pixels
        img = plt.imshow(subMatrix, cmap=colorMap, aspect=3.0)
    else:
        print("\nERROR: Unknown pixelType specified!\n")
        sys.exit(-1)
    
    #To make the colors consistent if there are submatrices.
    #plt.clim(np.min(scanningMatrix), np.max(scanningMatrix))
    if((climmin!=None) and (climmax!=None)):
        plt.clim([climmin, climmax])

    if(sequence!=None):
        mySeqFile = SeqIO.read(sequence, 'fasta')
        #print(mySeqFile)
        #Convert aaOrder to a list.
        aaOrderList = list(aaOrder)
        for i in range (len(subMatrix[0])):
            j = beg-1+i
            # print(i, aaOrderList.index(sequence[i]))
            plt.scatter(i, aaOrderList.index(mySeqFile.seq[j].upper()), s=5, c='black', marker='o')


        # if(onlyDNAaccessible==True):
        #     mySeqFile = SeqIO.read(sequence, 'fasta')
        #     onlyDNAaccessibleList = []
        #     #print(mySeqFile)
        #     #Convert aaOrder to a list.
        #     aaOrderList = list(aaOrder)
        #     for i in range (len(subMatrix[0])):
        #         j = beg-1+i
        #         # print(i, aaOrderList.index(sequence[i]))
        #         notInTheList = list(set(alphabeticalAminoAcidsList) - set(aaSingleNucleotideNeighborDict[mySeqFile.seq[j].upper()]))
        #         notInTheList.remove(mySeqFile.seq[j].upper())
        #         for item in notInTheList:
        #             # plt.scatter(i, aaOrderList.index(item), s=75, c='white', marker='x')
        #             plt.scatter(i, aaOrderList.index(item), s=216, c='white', marker='s')
        #             # onlyDNAaccessibleList.append([i, aaOrderList.index(item)])
        # #     for item in onlyDNAaccessibleList:
        # #         subMatrix[item[0]][item[1]] = np.nan


    if(isColorBarOn):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=fig_width/nres_shown, pad=0.2)
        # plt.clim([-10.0, 2.0])
        plt.colorbar(img, cax=cax)

    plt.tight_layout()
    plt.savefig(outFile+".png")
    if(interactive):
        plt.show()
    

    #plt.imsave('output.png', subMatrix)
    plt.close()

def plotDecoratedMatrix2(scanningMatrix, outFile, beg, end, \
                    colorMap = 'coolwarm', \
                    offSet=0, pixelType='square',
                    aaOrder="ACDEFGHIKLMNPQRSTVWY", \
                    sequence=None,\
                    interactive=True,\
                    isColorBarOn=False):
    """
        A function to plot deep mutational scanning matrices. 
  
    Parameters
    ----------
    scanningMatrix: numpy array of arrays
        Data matrix to plot

    outFile: string
        Name of the output image without file extension.
        Default file extension (png) is added by the program
    
    beg: int
        The first residue to use. It is used to select a subrange 
        of amino acids. It starts from 1.
    
    end: int
        The last residue to use. It is used to select a subrange 
        of amino acids.

    colorMap: matplotlib cmap
        Any colormap existing in matplotlib.
        Default is coolwarm. 

    offSet: int
        It is used to match the residue IDs in your PDB
        file with 0 based indices read from scanningMatrix matrix

    pixelType: string
        It can only have 'square' or 'rectangle' values.
        It is a matter of taste but I added it as an option.
        Default is 'square'

    sequence: string
        A fasta file of one letter amino acid codes from N terminal to C terminal. 

    interactive: bool
        If True, it will plot the map interactively. 

    isColorBarOn: bool
        If True, it will show a colorbar to show the numerical scale of
        the colors. Default is False. 

    Returns
    -------
    Nothing

    """

    #We subtract 1 from beg bc matrix indices starts from 0
    if(end == None):
        end = len(scanningMatrix[0])

    print("Beginning: "+str(beg))
    print("End      : "+str(end))
    
    print(len(scanningMatrix[0]))
    subMatrix = scanningMatrix[:, (beg-1):end]
    #print(subMatrix)

    ##########################################################################
    # Set plotting parameters
    nres_shown = len(subMatrix[0])
    fig_height=8
    # figure proportions
    fig_width = fig_height/2  # inches
    fig_width *= nres_shown/20
    
    fig, ax = plt.subplots(2, 1, figsize=(fig_width, fig_height))

    wspace = 0.5  # inches
    plt.subplots_adjust(wspace=wspace/fig_width, hspace=0.15)

    # figure structure
    gs = GridSpec(2, 1, height_ratios=[1, 5])
    ax0 = plt.subplot(gs[0])  # secondary structure strip
    ax1 = plt.subplot(gs[1])  # mutagenesis table
    # axcb = plt.subplot(gs[1, 1])  # colorbar
    # ax2 = plt.subplot(gs[2, 0])  # average profile

    if (nres_shown >150):
        majorTics = 50
    else:
        majorTics = 20


    major_nums_x = np.arange(majorTics, len(subMatrix[0]), majorTics, dtype=int)
    #major_nums_x = major_nums_x -1
    major_nums_x = np.insert(major_nums_x, 0, 0)
    #print(major_nums_x)

    minor_nums_x = np.arange(10, len(subMatrix[0]), 10, dtype=int)
    #minor_nums_x = minor_nums_x - 1
    minor_nums_x = np.insert(minor_nums_x, 0, 0)
    #print(minor_nums_x)

    major_labels_x = major_nums_x + 1 + offSet

    major_nums_y = np.arange(0, 20, 1, dtype=int)
    major_labels_y = list(aaOrder)

    plt.xticks(major_nums_x, major_labels_x, size=28)
    ax1.set_xticks(minor_nums_x, minor=True)
    
    plt.yticks(major_nums_y, major_labels_y, size=16)
    ax1.set_yticklabels(major_labels_y, ha='left')
    ax1.tick_params(axis='y', which='major', pad=30)

    
    #############################################################################
    if(pixelType=='square'):
        #For plotting square pixels
        img = plt.imshow(subMatrix, cmap=colorMap)
    elif(pixelType=='rectangle'):
        #For plotting rectangular pixels
        img = plt.imshow(subMatrix, cmap=colorMap, aspect=3.0)
    else:
        print("\nERROR: Unknown pixelType specified!\n")
        sys.exit(-1)
    
    #To make the colors consistent if there are submatrices.
    plt.clim(np.min(scanningMatrix), np.max(scanningMatrix)) 

    if(sequence!=None):
        mySeqFile = SeqIO.read(sequence, 'fasta')
        #Convert aaOrder to a list.
        aaOrderList = list(aaOrder)
        for i in range (len(subMatrix[0])):
            j = beg-1+i
            # print(i, aaOrderList.index(sequence[i]))
            plt.scatter(i, aaOrderList.index(mySeqFile.seq[j]), s=5, c='black', marker='o')
    
    if(isColorBarOn):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="1%", pad=0.2)
        plt.colorbar(img, cax=cax)

    decorated = True
    if(decorated == True):
        ax0.bar(np.arange(len(scanningMatrix[0])), np.random.rand(len(scanningMatrix[0])))


    plt.tight_layout()
    plt.savefig(outFile+".png")
    if(interactive):
        plt.show()
    

    #plt.imsave('output.png', subMatrix)
    plt.close()
def plotDecoratedMatrix(scanningMatrix, outFile, beg, end, \
                    colorMap = 'coolwarm', \
                    offSet=0, pixelType='square',
                    aaOrder="ACDEFGHIKLMNPQRSTVWY", \
                    sequence=None,\
                    interactive=True,\
                    isColorBarOn=False):
    """
        A function to plot deep mutational scanning matrices. 
  
    Parameters
    ----------
    scanningMatrix: numpy array of arrays
        Data matrix to plot

    outFile: string
        Name of the output image without file extension.
        Default file extension (png) is added by the program
    
    beg: int
        The first residue to use. It is used to select a subrange 
        of amino acids. It starts from 1.
    
    end: int
        The last residue to use. It is used to select a subrange 
        of amino acids.

    colorMap: matplotlib cmap
        Any colormap existing in matplotlib.
        Default is coolwarm. 

    offSet: int
        It is used to match the residue IDs in your PDB
        file with 0 based indices read from scanningMatrix matrix

    pixelType: string
        It can only have 'square' or 'rectangle' values.
        It is a matter of taste but I added it as an option.
        Default is 'square'

    sequence: string
        A fasta file of one letter amino acid codes from N terminal to C terminal. 

    interactive: bool
        If True, it will plot the map interactively. 

    isColorBarOn: bool
        If True, it will show a colorbar to show the numerical scale of
        the colors. Default is False. 

    Returns
    -------
    Nothing

    """

    #We subtract 1 from beg bc matrix indices starts from 0
    if(end == None):
        end = len(scanningMatrix[0])

    print("Beginning: "+str(beg))
    print("End      : "+str(end))
    
    print(len(scanningMatrix[0]))
    subMatrix = scanningMatrix[:, (beg-1):end]
    #print(subMatrix)

    ##########################################################################
    # Set plotting parameters
    nres_shown = len(subMatrix[0])
    fig_height=8
    # figure proportions
    fig_width = fig_height/2  # inches
    fig_width *= nres_shown/20
    
    fig, ax = plt.subplots(2, figsize=(fig_width, fig_height), \
                            sharex=True)
                            #, gridspec_kw={'height_ratios': [0.1, 1], 'width_ratios': [1]}
                            #, subplot_kw={'aspect':'equal'})
    wspace = 0.5  # inches
    plt.subplots_adjust(wspace=0)
    # adjust subplot sizes
    # gs = GridSpec(2, 1, width_ratios=[1])

    if (nres_shown >150):
        majorTics = 50
    else:
        majorTics = 20


    major_nums_x = np.arange(majorTics, len(subMatrix[0]), majorTics, dtype=int)
    #major_nums_x = major_nums_x -1
    major_nums_x = np.insert(major_nums_x, 0, 0)
    #print(major_nums_x)

    minor_nums_x = np.arange(10, len(subMatrix[0]), 10, dtype=int)
    #minor_nums_x = minor_nums_x - 1
    minor_nums_x = np.insert(minor_nums_x, 0, 0)
    #print(minor_nums_x)

    major_labels_x = major_nums_x + 1 + offSet

    major_nums_y = np.arange(0, 20, 1, dtype=int)
    major_labels_y = list(aaOrder)

    plt.xticks(major_nums_x, major_labels_x, size=28)
    ax[1].set_xticks(minor_nums_x, minor=True)
    
    plt.yticks(major_nums_y, major_labels_y, size=6)
    ax[1].set_yticklabels(major_labels_y, ha='left')
    ax[1].tick_params(axis='y', which='major', pad=30)

    
    #############################################################################
    if(pixelType=='square'):
        #For plotting square pixels
        img = plt.imshow(subMatrix, cmap=colorMap)
    elif(pixelType=='rectangle'):
        #For plotting rectangular pixels
        img = plt.imshow(subMatrix, cmap=colorMap, aspect=3.0)
    else:
        print("\nERROR: Unknown pixelType specified!\n")
        sys.exit(-1)
    
    #To make the colors consistent if there are submatrices.
    plt.clim(np.min(scanningMatrix), np.max(scanningMatrix)) 

    if(sequence!=None):
        mySeqFile = SeqIO.read(sequence, 'fasta')
        #Convert aaOrder to a list.
        aaOrderList = list(aaOrder)
        for i in range (len(subMatrix[0])):
            j = beg-1+i
            # print(i, aaOrderList.index(sequence[i]))
            plt.scatter(i, aaOrderList.index(mySeqFile.seq[j]), s=5, c='black', marker='o')
    
    if(isColorBarOn):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes("right", size="1%", pad=0.2)
        plt.colorbar(img, cax=cax)

    decorated = True
    if(decorated == True):
        ax[0].set_aspect(10)
        dataArray = getMinMaxData(scanningMatrix, 'min', outfile=None, printDetails=False)
        dataArray = rankdata((-1.0)*np.array(dataArray))/float(len(dataArray))
        ax[0].bar(np.arange(len(dataArray)), dataArray, color='red')


    
    
    if(interactive):
        plt.show()
    fig.tight_layout()
    fig.savefig(outFile+".png")

    #plt.imsave('output.png', subMatrix)
    plt.close()

def getMinMaxData(scanningMatrix, type, outfile, printDetails=False, ranksort=True):
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

    data = []
    
    if(type.lower()=='max'):
        for col in range(len(scanningMatrix[0])):
            X_max = (scanningMatrix.T[col].max())            
            # FILE.write("{}\t{}\n".format(col, X_max))
            data.append(X_max)
            
            if(printDetails):
                print(scanningMatrix.T[col])
        if(ranksort == True):
            data = rankdata(np.array(data))/float(len(data))

    elif(type.lower()=='min'):
        for col in range(len(scanningMatrix[0])):
            X_min = (scanningMatrix.T[col].min())
            # FILE.write("{}\t{}\n".format(col, X_min))
            data.append(X_min)

            if(printDetails):
                print(scanningMatrix.T[col])
        if(ranksort == True):
            data = rankdata((-1.0)*np.array(data))/float(len(data))

    else:
        print("ERROR: Unknown type!")
        print("       Type can only be min or max!")
        sys.exit(-1)
    if(outfile != None):
        FILE = open(outfile, "w")
        for i in range (len(data)):
            FILE.write("{}\t{}\n".format(i, data[i]))
        FILE.close()
    return data

def plotExperimentalMatrix(scanningMatrix, outFile, beg, end, \
                    colorMap = 'coolwarm', \
                    offSet=0, pixelType='square',
                    aaOrder="ACDEFGHIKLMNPQRSTVWY", \
                    sequence=None,\
                    interactive=True,\
                    isColorBarOn=False,\
                    onlyDNAaccessible=None,\
                    climmin=None,\
                    climmax=None):
    """
        A function to plot deep mutational scanning data from
        the Rieselamn dataset. 
  
    Parameters
    ----------
    scanningMatrix: numpy array of arrays
        Data matrix to plot

    outFile: string
        Name of the output image without file extension.
        Default file extension (png) is added by the program
    
    beg: int
        The first residue to use. It is used to select a subrange 
        of amino acids. It starts from 1.
    
    end: int
        The last residue to use. It is used to select a subrange 
        of amino acids.

    colorMap: matplotlib cmap
        Any colormap existing in matplotlib.
        Default is coolwarm. 

    offSet: int
        It is used to match the residue IDs in your PDB
        file with 0 based indices read from scanningMatrix matrix

    pixelType: string
        It can only have 'square' or 'rectangle' values.
        It is a matter of taste but I added it as an option.
        Default is 'square'

    sequence: string
        A fasta file of one letter amino acid codes from N terminal to C terminal. 

    interactive: bool
        If True, it will plot the map interactively. 

    isColorBarOn: bool
        If True, it will show a colorbar to show the numerical scale of
        the colors. Default is False.
        
    onlyDNAaccessible: string
        If it is not None, you must provide a DNA sequence of the protein. 
        It will calculate mutations that are accessible from 'only' one 
        nucleotide change. The other locations will be shown as white. 

    climmin: float
        It is used to put a limit to minimum value of the colorbar. 

    climmax: float
        It is used to put an upper limit to maximum value of the colorbar. 

    Returns
    -------
    Nothing

    """

    scanningMatrix = np.ma.array (scanningMatrix, mask=np.isnan(scanningMatrix))

    #We subtract 1 from beg bc matrix indices starts from 0
    if(end == None):
        end = len(scanningMatrix[0])

    print("Beginning: "+str(beg))
    print("End      : "+str(end))
    
    print(len(scanningMatrix[0]))
    subMatrix = scanningMatrix[:, (beg-1):end]
    #print(subMatrix)

    ##########################################################################
    # Set plotting parameters
    nres_shown = len(subMatrix[0])
    fig_height=8
    # figure proportions
    fig_width = fig_height/2  # inches
    fig_width *= nres_shown/20
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    if (nres_shown >150):
        majorTics = 50
    else:
        majorTics = 20


    major_nums_x = np.arange(majorTics, len(subMatrix[0]), majorTics, dtype=int)
    #major_nums_x = major_nums_x -1
    major_nums_x = np.insert(major_nums_x, 0, 0)
    #print(major_nums_x)

    minor_nums_x = np.arange(10, len(subMatrix[0]), 10, dtype=int)
    #minor_nums_x = minor_nums_x - 1
    minor_nums_x = np.insert(minor_nums_x, 0, 0)
    #print(minor_nums_x)

    major_labels_x = major_nums_x + 1 + offSet

    major_nums_y = np.arange(0, 20, 1, dtype=int)
    major_labels_y = list(aaOrder)

    plt.xticks(major_nums_x, major_labels_x, size=28)
    ax.set_xticks(minor_nums_x, minor=True)
    
    plt.yticks(major_nums_y, major_labels_y, size=16)
    ax.set_yticklabels(major_labels_y, ha='left')
    ax.tick_params(axis='y', which='major', pad=30)

    #############################################################################
    onlyDNAaccessibleList = []
    current_cmap = matplotlib.cm.get_cmap()
    current_cmap.set_bad(color='white')
    if(sequence!=None):
        mySeqFile = SeqIO.read(sequence, 'fasta')

        # onylDNAaccessible = True
        if(onlyDNAaccessible!=None):
            mySeqFile = SeqIO.read(sequence, 'fasta')
            seqDictionary = parseDNAsequence(onlyDNAaccessible, sequence)
            # print(seqDictionary)

            print("I am plotting only DNA accessible mutations.")
            #print(mySeqFile)
            #Convert aaOrder to a list.
            aaOrderList = list(aaOrder)
            for i in range (len(subMatrix[0])):
                j = beg-1+i
                # print(i, aaOrderList.index(sequence[i]))

                dnaAccessibleList = getPossible1NucleotideNeighbors(seqDictionary[j])
                # print(dnaAccessibleList)
                notInTheList = list(set(alphabeticalAminoAcidsList) - set(dnaAccessibleList))
                if mySeqFile.seq[j].upper() in notInTheList:
                    notInTheList.remove(mySeqFile.seq[j].upper())
                for item in notInTheList:
                    # plt.scatter(i, aaOrderList.index(item), s=75, c='white', marker='x')
                    #plt.scatter(i, aaOrderList.index(item), s=200, c='white', marker='s')
                    onlyDNAaccessibleList.append([i, aaOrderList.index(item)])
            for item in onlyDNAaccessibleList:
                # print(item)
                subMatrix[item[1]][item[0]] = np.nan
    #############################################################################
    if(pixelType=='square'):
        #For plotting square pixels
        img = plt.imshow(subMatrix, cmap=colorMap)
    elif(pixelType=='rectangle'):
        #For plotting rectangular pixels
        img = plt.imshow(subMatrix, cmap=colorMap, aspect=3.0)
    else:
        print("\nERROR: Unknown pixelType specified!\n")
        sys.exit(-1)
    
    #To make the colors consistent if there are submatrices.
    # plt.clim(np.min(scanningMatrix), np.max(scanningMatrix))

    #This value here does not make sense if you are really plotting experimental data.
    #I had to do this to make the colors consistent to plot only DNA accessible mutations.
    #The source data for them is frequency modifed ESGEMME scores and they are typically 
    #between -10 and 2.
    if((climmin!=None) and (climmax!=None)):
        plt.clim([climmin, climmax])

    if(sequence!=None):
        mySeqFile = SeqIO.read(sequence, 'fasta')
        #Convert aaOrder to a list.
        aaOrderList = list(aaOrder)
        for i in range (len(subMatrix[0])):
            j = beg-1+i
            # print(i, aaOrderList.index(sequence[i]))
            plt.scatter(i, aaOrderList.index(mySeqFile.seq[j].upper()), s=5, c='black', marker='o')
    
    if(isColorBarOn):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=fig_width/nres_shown, pad=0.2)
        # plt.clim([-10.0, 2.0])
        plt.colorbar(img, cax=cax)

    plt.tight_layout()
    plt.savefig(outFile+".png")
    if(interactive):
        plt.show()
    

    #plt.imsave('output.png', subMatrix)
    plt.close()


def parseRHAPSODYoutput(inputFile, field=10):
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
    rhapsodyDataFile=open(inputFile, 'r')
    allLines = rhapsodyDataFile.readlines()
    rhapsodyDataFile.close()
    matrixData = []

    if (allLines[0][0]=='#'):
        print(allLines[0])
        del allLines[0]

    seqLength = int(len(allLines)/19)
    for i in range (0, seqLength):

        #Determine the missing amino acid and its index
        missingAminoAcidsList = []
        for j in range (0, 19):
            #print(i*19+j)
            data = allLines[i*19+j].split()
            #print(data)
            missingAminoAcidsList.append(data[3])
        
        # print(missingAminoAcidsList)
        missing = set(alphabeticalAminoAcidsList).difference(missingAminoAcidsList)
        # print(missing)
        missingIndex = alphabeticalAminoAcidsList.index(str(list(missing)[0]))
        # print(missingIndex)


        #Read the data and append it to a temporary list
        tempList = []
        for j in range (0, 19):
            data = allLines[i*19+j].split()
            if (data[field] == 'nan'):
                tempList.append(0.000)
            else:
                tempList.append(float(data[field]))
        
        
        
        #Fill the index location for the missing amino acid. 
        tempList.insert(missingIndex, 0.000)
        matrixData.append(tempList)            

    mutationsData = np.array(matrixData).T
    # mutatedTransposed = mutationsData.T
    # return mutatedTransposed
    return mutationsData

def parseFOLDXoutput(inputFile, colorThreshhold=10.0, colorCorrect=True):
    """
        Parse FoldX output files

    Parameters
    ----------
    inputFile: string
        Name of the FoldX output file

    colorThreshhold: float
        If colorCorrect True, it will reduce values greater
        than colorThreshhold to make the deep mutational scanning
        matrices look better.
    colorCorrect: bool
        If True, it will change values the values greater than
        the colorThreshold to  colorThreshold 

    Returns
    -------
    Data as a Numpy array, where each col contains 
        delta delta G for an amino acid mutation.
    """
    foldxDataFile=open(inputFile, 'r')
    allLines = foldxDataFile.readlines()
    foldxDataFile.close()
    matrixData = []

    if (allLines[0][0]=='#'):
        print(allLines[0])
        del allLines[0]
    A_list = []
    C_list = []
    D_list = []
    E_list = []
    F_list = []
    G_list = []
    H_list = []
    I_list = []
    K_list = []
    L_list = []
    M_list = []
    N_list = []
    P_list = []
    Q_list = []
    R_list = []
    S_list = []
    T_list = []
    V_list = []
    W_list = []
    Y_list = []
    chunkSize = 21
    seqLength = int(len(allLines)/chunkSize)
    for i in range (0, seqLength):

        #Append the data to corresponding line 
        #Please note that the foldx output is in a 
        #diffent order than the alphabetical single 
        #letter amino acids order. 

        G_list.append(float(allLines[i*chunkSize+1].split()[1]))
        A_list.append(float(allLines[i*chunkSize+2].split()[1]))
        L_list.append(float(allLines[i*chunkSize+3].split()[1]))
        V_list.append(float(allLines[i*chunkSize+4].split()[1]))
        I_list.append(float(allLines[i*chunkSize+5].split()[1]))
        P_list.append(float(allLines[i*chunkSize+6].split()[1]))
        R_list.append(float(allLines[i*chunkSize+7].split()[1]))
        T_list.append(float(allLines[i*chunkSize+8].split()[1]))
        S_list.append(float(allLines[i*chunkSize+9].split()[1]))
        C_list.append(float(allLines[i*chunkSize+10].split()[1]))
        M_list.append(float(allLines[i*chunkSize+11].split()[1]))
        K_list.append(float(allLines[i*chunkSize+12].split()[1]))
        E_list.append(float(allLines[i*chunkSize+13].split()[1]))
        Q_list.append(float(allLines[i*chunkSize+14].split()[1]))
        D_list.append(float(allLines[i*chunkSize+15].split()[1]))
        N_list.append(float(allLines[i*chunkSize+16].split()[1]))
        W_list.append(float(allLines[i*chunkSize+17].split()[1]))
        Y_list.append(float(allLines[i*chunkSize+18].split()[1]))
        F_list.append(float(allLines[i*chunkSize+19].split()[1]))
        H_list.append(float(allLines[i*chunkSize+20].split()[1]))
        
    matrixData.append(A_list)
    matrixData.append(C_list)
    matrixData.append(D_list)
    matrixData.append(E_list)
    matrixData.append(F_list)
    matrixData.append(G_list)
    matrixData.append(H_list)
    matrixData.append(I_list)
    matrixData.append(K_list)
    matrixData.append(L_list)
    matrixData.append(M_list)
    matrixData.append(N_list)
    matrixData.append(P_list)
    matrixData.append(Q_list)
    matrixData.append(R_list)
    matrixData.append(S_list)
    matrixData.append(T_list)
    matrixData.append(V_list)
    matrixData.append(W_list)
    matrixData.append(Y_list)
    
    mutationsData = np.array(matrixData)
    if(colorCorrect):
        mutationsData[mutationsData > colorThreshhold] = colorThreshhold
    # mutatedTransposed = mutationsData.T
    # return mutatedTransposed
    # print(mutationsData[0])
    return mutationsData

def writeGEMMEmatrix(scanningMatrix, outFile, beg, end, \
                    aaOrder = alphabeticalAminoAcidsList, \
                    residueList = None,
                    offSet=0):
    """
        A function to write deep mutational scanning matrices. 

    The format, even thought weird, is as follows:
        X1      X2      ...     XN   
    "A" 0.5     0.1     ...     0.2
    "C" 0.0     0.9     ...     0.4
    .
    .
    .
    "Y" 0.0     0.8     ...     0.7

    It can be read as a dataframe with pandas or any other library.
    X is the unmutated residues and 1,2....N are residues in the protein.
  
    Parameters
    ----------
    scanningMatrix: numpy array of arrays
        Data matrix to plot

    outFile: string
        Name of the output png image
    
    # beg: int
    #     The first residue to use. It is used to select a subrange 
    #     of amino acids.
    
    # end: int
    #     The last residue to use. It is used to select a subrange 
    #     of amino acids.

    aaOrder: Python list
        A list of 20 amino acid characters such as ['A', 'C',...., 'Y'].  
        Default is alphabetical list.

    residueList: Python list
        A list of the residues in the protein such as ['X1', 'X2', ..., 'XN']. 
        X is one-letter amino acid code and N is the total number of residues. 
        Default is None.
        If not provided, the program will write them starting from 1 till the end. 

    offSet: int
        It is used to match the residue IDs in your PDB
        file with 0 based indices read from scanningMatrix matrix

    Returns
    -------
    Nothing

    """
    debug = 0

    #We subtract 1 from beg bc matrix indices starts from 0
    if(end == None):
        end = len(scanningMatrix[0])

    if(debug):
        print("Beginning: "+str(beg))
        print("End      : "+str(end))
    
        print(len(scanningMatrix[0]))

    # I don't want to write the submatrices for now!
    # Let's keep it simple for now! 
    # subMatrix = scanningMatrix[:, (beg-1):end]
    realResidueList = []
    if(residueList == None):
        for i in range(0,len(scanningMatrix[0])):
            tempString = 'X'+str(i+1)
            realResidueList.append(tempString)
    else:
        realResidueList = residueList

    with open(outFile, 'w') as file:
        
        for res in range(len(realResidueList)):
            #file.write(' ')
            file.write(realResidueList[res])
            if(res != (len(realResidueList)-1)):
                file.write(' ')
        file.write('\n')

        for i in range(20):
            file.write("\""+aaOrder[i]+"\"")
            file.write(' ')
            for res in range (len(scanningMatrix[i])):
                file.write(str(scanningMatrix[alphabeticalAminoAcidsList.index(aaOrder[i])][res]))
                if(res != (len(realResidueList)-1)):
                    file.write(' ')
            file.write('\n')

def writeDMSformat(scanningMatrix, outFile, residueList, 
                    beg, end, \
                    aaOrder = alphabeticalAminoAcidsList, \
                    offSet=0):
    """
        A function to convert GEMME matrix format data to experimental 
        deep mutational scanning format. 

    The format has only two columns. In the first column, one-letter wild
    type amino acid, residue number and mutated amino acid are combined without any
    gap or any other character, such as A19C. The second column is typically a float
    value of the predicted effect.
    
    Residue numbering starts with 1.
    "A1A 0.211
    "A1C 0.124
    .
    .
    .
    "T262 0.434

    Parameters
    ----------
    scanningMatrix: numpy array of arrays
        Data matrix to plot

    outFile: string
        Name of the output png image
    
    residueList: Python list
        A list of the residues in the protein such as ['X1', 'X2', ..., 'XN']. 
        X is one-letter amino acid code and N is the total number of residues. 
        Default is None.
        If not provided, the program will write them starting from 1 till the end. 

    beg: int
        The first residue in the fasta file may not necessarily be the first residue
        in the protein. This beg parameter gives us a chance to adjust the first residue
        other than one. It only effects the DMS format text output. Default value is 0.
    
    end: int
        The last residue to use. It is used to select a subrange 
        of amino acids.

    aaOrder: Python list
        A list of 20 amino acid characters such as ['A', 'C',...., 'Y'].  
        Default is alphabetical list.

    offSet: int
        It is used to match the residue IDs in your PDB
        file with 0 based indices read from scanningMatrix matrix

    Returns
    -------
    Nothing

    """
    debug = 0

    #We subtract 1 from beg bc matrix indices starts from 0
    if(end == None):
        end = len(scanningMatrix[0])

    if(debug):
        print("Beginning: "+str(beg))
        print("End      : "+str(end))
    
        print(len(scanningMatrix[0]))

    # I don't want to write the submatrices for now!
    # Let's keep it simple for now! 
    # subMatrix = scanningMatrix[:, (beg-1):end]
    realResidueList = []
    if(residueList == None):
        for i in range(0,len(scanningMatrix[0])):
            tempString = 'X'+str(i+1+beg)
            realResidueList.append(tempString)
    else:
        for i in range(0,len(scanningMatrix[0])):
            tempString = residueList[i].upper()+str(i+1+beg)
            realResidueList.append(tempString)

    with open(outFile, 'w') as file:
        for res in range (len(scanningMatrix[0])):
            for i in range(20):
                file.write(realResidueList[res]+aaOrder[i]+" ")
                file.write(str(scanningMatrix[alphabeticalAminoAcidsList.index(aaOrder[i])][res])+"\n")

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

def rankNormalization(scanningMatrix, printDetails=False):
    """
        This function normalizes mutations for each column (namely,
        over all variants at a certain position) and gives them a 
        score between 0 and 1. Obviously, 0 means low impact and 
        1 means the highest impact of the mutation.  


    Parameters
    ----------
    scanningMatrix: numpy array of arrays
        Data matrix to rank normalize
    printDetails: a boolean value
        If True, it will print details for debugging (Default: False)
    Returns
    -------
    rankNormalizedData: numpy array of arrays
        Rank normalized data matrix
    """

    rankNormalizedData = []
    for col in range(len(scanningMatrix[0])):

        X_min = (scanningMatrix.T[col].min())
        X_max = (scanningMatrix.T[col].max())

        temp = []
        if(X_min != X_max):
            for x in scanningMatrix.T[col]:
                y = (x - X_min)/(X_max-X_min)
                temp.append(y)
        else:
            for x in scanningMatrix.T[col]:
                temp.append(0.0)

        
        if(printDetails):
            print(scanningMatrix.T[col])
            print(temp)

        rankNormalizedData.append(temp)

    return np.array(rankNormalizedData).T

def parseExperimentalData(inputcsv, experiment='DMS_score', 
                            outputcsv="convertedMatrix.csv",
                            debug = False, datasource="proteingym"):
    """
        Parse Riesselman or ProteinGym experimental data files from the
        https://www.nature.com/articles/s41592-018-0138-4 and 
        https://doi.org/10.48550/arXiv.2205.13760


    Parameters
    ----------
    inputcsv: string
        Name of the experimental input data file in csv format
    
    experiment: string
        Name of the experiment such as score, screenscore, fitness,
        abundance. To parse ProteinGym data, set this value to 
        'DMS_score'.

    outputcsv: string
        Name of the output data file in csv format.

    debug: bool
        If True, it will print some extra data to stdout. 
    
    datasource: string
        It should take 'proteingym' value if you want to parse
        ProteinGym data. Otherwise, don't use it. In this case,
        it assumes that the datasource is Riesseman, 2016 csv files.

    Returns
    -------
    Data as a Numpy array, where each col contains 
        an experimental parameter for an amino acid mutation.
        Rows of the matrix are amino acids in alphabetical order. 
    """
    if(datasource == "proteingym"):
        df = pd.read_csv(inputcsv)
    else:
        df = pd.read_csv(inputcsv, sep=";", decimal=",")

    #Add three new columns.
    df['wt'] = ""
    df['mutation']=""
    df['resid']=""

    if(debug):
        print (df.columns)

    # Fill the 
    for index, row in df.iterrows():
        df['wt'].iloc[index] = row['mutant'][0]
        df['mutation'].iloc[index] = row['mutant'][-1]
        df['resid'].iloc[index] = int(row['mutant'][1:-1])
        #print(index, row['mutant'], row['DMS_score'])

    #Create a new dataframe from only the necessary columns
    new_df = df[['wt','resid', 'mutation', experiment]].copy()
    if(debug):
        print(list(new_df['resid']))

    minResid = new_df['resid'].min()
    maxResid = new_df['resid'].max()
    
    numCols = maxResid - minResid + 1

    print(numCols)

    # dmsMap  = pd.DataFrame(index=range(numRows), columns=['resid']+alphabeticalAminoAcidsList, dtype=float)

    colsList = [str(x) for x in np.arange(1, numCols+1)] 
    #dmsMap  = pd.DataFrame(index=range(20), columns=colsList, dtype=float)
    dmsMap  = pd.DataFrame(index=alphabeticalAminoAcidsList, columns=colsList, dtype=float)
    
    #print(dmsMap.columns)
    #print(['resid']+alphabeticalAminoAcidsList)

    for index, row in new_df.iterrows():
        # print(index, row)
        for aa in alphabeticalAminoAcidsList:
            if(row['mutation'] == aa):
                dmsMap[str(row['resid']-minResid + 1)].loc[aa] = row[experiment]
    
    if(outputcsv != None):
        dmsMap.to_csv(outputcsv, index=True)

    return dmsMap.to_numpy()

def parseSingleLineData(inputcsv, experiment='DMS_score', 
                            outputcsv="convertedMatrix.csv",
                            debug = False):
    """
        Parse Singleline Format data


    Parameters
    ----------
    inputcsv: string
        Name of the experimental input data file in csv format
    
    experiment: string
        Name of the experiment such as score, screenscore, fitness,
        abundance. To parse ProteinGym data, set this value to 
        'DMS_score'.

    outputcsv: string
        Name of the output data file in csv format.

    debug: bool
        If True, it will print some extra data to stdout. 
    
    datasource: string
        It should take 'proteingym' value if you want to parse
        ProteinGym data. Otherwise, don't use it. In this case,
        it assumes that the datasource is Riesseman, 2016 csv files.

    Returns
    -------
    Data as a Numpy array, where each col contains 
        an experimental parameter for an amino acid mutation.
        Rows of the matrix are amino acids in alphabetical order. 
    """

    df = pd.read_csv(inputcsv, header=None, sep='\s+')
    df.columns=['mutant', experiment]
    print(df)
    #Add three new columns.
    df['wt'] = ""
    df['mutation']=""
    df['resid']=""

    if(debug):
        print (df.columns)

    # Fill the 
    for index, row in df.iterrows():
        # df['wt'].iloc[index] = row['mutant'][0]
        # df['mutation'].iloc[index] = row['mutant'][-1]
        # df['resid'].iloc[index] = int(row['mutant'][1:-1])
        #print(index, row['mutant'], row['DMS_score'])

        df.loc[index, 'wt'] = row['mutant'][0]
        df.loc[index, 'mutation'] = row['mutant'][-1]
        df.loc[index, 'resid'] = int(row['mutant'][1:-1])

    #Create a new dataframe from only the necessary columns
    new_df = df[['wt','resid', 'mutation', experiment]].copy()
    if(debug):
        print(list(new_df['resid']))

    minResid = new_df['resid'].min()
    maxResid = new_df['resid'].max()
    
    numCols = maxResid - minResid + 1

    print(numCols)

    # dmsMap  = pd.DataFrame(index=range(numRows), columns=['resid']+alphabeticalAminoAcidsList, dtype=float)

    colsList = [str(x) for x in np.arange(1, numCols+1)] 
    #dmsMap  = pd.DataFrame(index=range(20), columns=colsList, dtype=float)
    dmsMap  = pd.DataFrame(index=alphabeticalAminoAcidsList, columns=colsList, dtype=float)
    
    #print(dmsMap.columns)
    #print(['resid']+alphabeticalAminoAcidsList)

    for index, row in new_df.iterrows():
        # print(index, row)
        for aa in alphabeticalAminoAcidsList:
            if(row['mutation'] == aa):
                dmsMap[str(row['resid']-minResid + 1)].loc[aa] = row[experiment]
    
    if(outputcsv != None):
        dmsMap.to_csv(outputcsv, index=True)

    return dmsMap.to_numpy()

