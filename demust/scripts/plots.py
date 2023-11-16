from typing import OrderedDict
from demust.io import *
import argparse
from scipy.stats import rankdata
import pandas as pd
import prody
import numpy as np
import sys

# def getMinMaxData(scanningMatrix, type, outfile, printDetails=False):
#     """
#         This function gets min and max for each column (namely,
#         over all variants at a certain position). 

#     Parameters
#     ----------
#     scanningMatrix: numpy array of arrays
#         Data matrix to rank normalize
#     type: string
#         type can only be min or max.    
#     outfile: string
#         Prefix for the output file.
#     printDetails: a boolean value
#         If True, it will print details for debugging (Default: False)
#     Returns
#     -------
#     Nothing
#     """

#     data = []
#     FILE = open(outfile, "w")
#     if(type.lower()=='max'):
#         for col in range(len(scanningMatrix[0])):
#             X_max = (scanningMatrix.T[col].max())            
#             FILE.write("{}\t{}\n".format(col, X_max))
#             data.append(X_max)
            
#             if(printDetails):
#                 print(scanningMatrix.T[col])
        
#     elif(type.lower()=='min'):
#         for col in range(len(scanningMatrix[0])):
#             X_min = (scanningMatrix.T[col].min())
#             FILE.write("{}\t{}\n".format(col, X_min))
#             data.append(X_min)

#             if(printDetails):
#                 print(scanningMatrix.T[col])
#     else:
#         print("ERROR: Unknown type!")
#         print("       Type can only be min or max!")
#         sys.exit(-1)
    
#     FILE.close()
#     return data

def plot1DHeatMap(dataArray, outFile, beg, end, \
                colorMap = 'Reds', \
                offSet=0, pixelType='square',
                sequence=None,\
                interactive=True, isColorBarOn=False):
    """
        A function to plot deep mutational scanning matrices. 
  
    Parameters
    ----------
    dataArray: numpy array
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

    colorMap: matplotlib color
        Any color existing in matplotlib.
        Default is red. 

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
    # print(dataArray)
    # print(np.arange(len(dataArray)))
    #dataArray = np.ma.array (dataArray, mask=np.isnan(dataArray))
    #We subtract 1 from beg bc matrix indices starts from 0
    if(end == None):
        end = len(dataArray)


    print("Beginning: "+str(beg))
    print("End      : "+str(end))
    
    # print(len(dataArray))
    subArray = dataArray[(beg-1):end]
    #print(subArray)

    ##########################################################################
    # Set plotting parameters
    nres_shown = len(subArray)
    fig_height=8
    # figure proportions
    fig_width = fig_height/2  # inches
    fig_width *= nres_shown/20
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height/4))

    if (nres_shown >150):
        majorTics = 50
    else:
        majorTics = 20


    major_nums_x = np.arange(majorTics, len(subArray), majorTics, dtype=int)
    #major_nums_x = major_nums_x -1
    major_nums_x = np.insert(major_nums_x, 0, 0)
    #print(major_nums_x)
    minor_nums_x = np.arange(10, len(subArray), 10, dtype=int)
    #minor_nums_x = minor_nums_x - 1
    minor_nums_x = np.insert(minor_nums_x, 0, 0)
    #print(minor_nums_x)
    major_labels_x = major_nums_x + 1 + offSet
    plt.xticks(major_nums_x, major_labels_x, size=28)
    ax.set_xticks(minor_nums_x, minor=True)
    
    major_nums_y = np.arange(0, 1, 0.5)
    # major_labels_y = list(major_nums_y)
    major_labels_y = [" ", " "]
    plt.yticks(major_nums_y, major_labels_y, size=16)
    # ax.set_yticklabels([str(round(float(label), 2)) for label in major_labels_y], ha='left')
    ax.tick_params(axis='y', which='major', pad=30)

    
    #############################################################################
    if(pixelType=='square'):
        #For plotting square pixels
        #plt.xlim([-0.1,len(subArray)])
        #img = plt.bar(np.arange(len(subArray)), subArray, color="red")
        img = plt.imshow([subArray], cmap=colorMap)
    elif(pixelType=='rectangle'):
        # For plotting rectangular pixels
        img = plt.imshow([subArray], cmap=colorMap, aspect=3.0)
    else:
        print("\nERROR: Unknown pixelType specified!\n")
        sys.exit(-1)
    
    #To make the colors consistent if there are submatrices.
    #plt.clim(np.min(scanningMatrix), np.max(scanningMatrix)) 
    if(isColorBarOn):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=fig_width/nres_shown, pad=0.2)
        plt.clim([0.0, 1.0])
        plt.colorbar(img, cax=cax)

    plt.tight_layout()
    plt.savefig(outFile+".png")
    if(interactive):
        plt.show()
    
    plt.close()


def option6(a, order='descending'):
    b = -a if order=='descending' else a        
    idx = b.argsort(0,'stable')
    n = idx.shape
    out = np.empty(n, dtype=float)
    np.put_along_axis(out, idx, np.arange(1,n+1), axis=0)
    return np.where(np.isnan(a), np.nan, out)

def calculateReliabilityScore(inputfile):
    from Bio import AlignIO
    dataArray = []
    alignment = AlignIO.read(inputfile, "fasta")
    proteinSize = (len(alignment[0].seq))
    numberOfRows = (len(alignment))
    align_array = np.array([list(rec) for rec in alignment], np.unicode_).transpose()

    for i in range(0, proteinSize):
        #print(i)
        gapPerColumnCount=list(align_array[i]).count("-")
        dataArray.append(1.0 - (gapPerColumnCount/float(numberOfRows)))
    
    return dataArray

def plotsApp(args):
    # if args.inputfile == None:
    #     print('Usage: python demust.py [-h] [-i INPUTFILE] [-d GEMME]')
    #     print('\nError: missing arguments: Please provide --inputfile and/or --datatype')
    #     sys.exit(-1)
    debug = True
    print("\nRunning 'demust plots' app...\n")
    if (args.datatype.lower()=='gemme'): 
        scanningMatrix = parseGEMMEoutput(args.inputfile, verbose=False)

        if((args.type.lower()=='min')):
            dataArray = getMinMaxData(scanningMatrix, args.type, \
                                      outfile=args.outputfile+".dat",\
                                      printDetails=False, ranksort=True)
            #dataArray = rankdata((-1.0)*np.array(dataArray))/float(len(dataArray))
            plot1DHeatMap(dataArray, args.outputfile, beg=args.beginning, end=args.end, \
                         colorMap = 'Reds', \
                         offSet=0, pixelType='square',
                         interactive=False)
        else:
            print("ERROR: Unknown plot type!{}".format(args.type.lower()))
            sys.exit(-1)
    elif (args.datatype.lower()=='riesselman'):

        # Please note that we are not using Riesselman, 2016 csv files here directly. 
        # We are using the newly created csv files produced by 
        # demust riesselman utility. The converted format looks like old GEMME format.
        # Namely, each column contains the results of 20 mutations of a certain aa position. 

        df = pd.read_csv(args.inputfile)
        #Drop the aa names column
        df = df.iloc[: , 1:]

        scanningMatrix = df.to_numpy()
        # print(scanningMatrix)   
        scanningMatrix = np.ma.array(scanningMatrix, mask=np.isnan(scanningMatrix))
        if((args.type.lower()=='max') or (args.type.lower()=='min')):
            # Please note that max value may/might represent wild-type like behaviour. 
            # Here, we are trying to highlight the most deleterious locations.
            # The following four lines of data processing are performed for this purpose. 
            dataArray = getMinMaxData(scanningMatrix, args.type, \
                                      outfile=args.outputfile+".dat",\
                                      printDetails=False, ranksort=False)

            # print(dataArray)
            data = np.nan_to_num(np.array(dataArray), nan=np.array(dataArray).max())
            # np.ma.array(dataArray, mask=np.isnan(dataArray))    
            # print(data)
            # #dataArray = rankdata((-1.0)*np.array(dataArray))/float(len(dataArray))
            temp = []
            nanLocations = []
            j = 0
            nonNanLocations = OrderedDict()
            for i in range(len(data)):
                if(np.isnan(data[i])):
                    nanLocations.append(i)
                else:
                    temp.append(data[i])
                    nonNanLocations[i]=j
                    j = j+1
            if((args.type.lower()=='min')):
                temp = 1.0 - (rankdata(np.array(temp))/float(len(temp)))
            else:
                temp = (rankdata(np.array(temp))/float(len(temp)))
            for key, value in nonNanLocations.items():
                # print(key, value)
                dataArray[key] = temp[nonNanLocations[key]]
            np.savetxt(args.outputfile+".dat", dataArray)
            plot1DHeatMap(np.array(dataArray), args.outputfile, beg=args.beginning, end=args.end, \
                          colorMap = 'Reds', \
                          offSet=0, pixelType='square',\
                          interactive=False)

        else:
            print("ERROR: Unknown plot type!{}".format(args.type.lower()))
            sys.exit(-1)
    elif (args.datatype.lower()=='jet'):
        df = pd.read_csv(args.inputfile, delimiter=r"\s+")
        
        if((args.type.lower()=='trace')):
            if 'trace' in df.columns:
                dataArray = (df['trace'].to_numpy())
                if(debug):
                    print(dataArray)
                plot1DHeatMap(dataArray, args.outputfile, beg=args.beginning, end=args.end, \
                            colorMap = 'Greens', \
                            offSet=0, pixelType='square',\
                            interactive=False, isColorBarOn=True)
            else:
                print("ERROR: The file does not contain a column called 'trace'!")
                sys.exit(-1)
        elif((args.type.lower()=='pc')):
            if 'pc' in df.columns:
                dataArray = (df['pc'].to_numpy())
                if(debug):
                    print(dataArray)
                plot1DHeatMap(dataArray, args.outputfile, beg=args.beginning, end=args.end, \
                            colorMap = 'Blues', \
                            offSet=0, pixelType='square',\
                            interactive=False, isColorBarOn=True)
            else:
                print("ERROR: The file does not contain a column called 'pc'!")
                sys.exit(-1)
        elif((args.type.lower()=='cv')):
            if 'cv' in df.columns:
                dataArray = (df['cv'].to_numpy())
                if(debug):
                    print(dataArray)
                plot1DHeatMap(dataArray, args.outputfile, beg=args.beginning, end=args.end, \
                            colorMap = 'Oranges', \
                            offSet=0, pixelType='square',\
                            interactive=False, isColorBarOn=True)
            else:
                print("ERROR: The file does not contain a column called 'cv'!")
                sys.exit(-1)
        else:
            print("ERROR: You can only obtain trace, pc or cv from a jet results file!")
            sys.exit(-1)
    
    elif (args.datatype.lower()=='pdb'):
        performDSSP(args.inputfile, parseall=False, stderr=True) 
        structure = parsePDB(args.inputfile)

    elif (args.datatype.lower()=='fasta'):
        from Bio import AlignIO
        dataArray = []
        alignment = AlignIO.read(args.inputfile, "fasta")
        proteinSize = (len(alignment[0].seq))
        numberOfRows = (len(alignment))
        align_array = np.array([list(rec) for rec in alignment], np.unicode_).transpose()

        for i in range(0, proteinSize):
            #print(i)
            #Reliability calculation
            gapPerColumnCount=list(align_array[i]).count("-")
            dataArray.append(1.0 - (gapPerColumnCount/float(numberOfRows)))
        
        if(args.end == None):
            args.end = len(dataArray)
            print("Length of data array"+str(args.end))
        if(args.paginate!=0):
            sequenceLength = (args.end - args.beginning - 1) # -1 is for starting the count from 0. 

            rowLength = args.paginate
            numberOfImageChunks = int(int(sequenceLength)/int(rowLength))
            
            for i in range(numberOfImageChunks):
                # plotGEMMEmatrix(gemmeData, args.outputfile+"_part_"+str(i+1), \
                #     i*rowLength + args.beginning, \
                #     (i+1)*rowLength + args.beginning -1,\
                #     colorMap=args.colormap, offSet=i*rowLength + args.offset, pixelType='square',\
                #     aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
                plot1DHeatMap(dataArray, args.outputfile+"_part_"+str(i+1), beg=i*rowLength + args.beginning, \
                    end=(i+1)*rowLength + args.beginning -1, \
                    colorMap = 'Greys', \
                    offSet=i*rowLength + args.offset, pixelType='square',\
                    interactive=False, isColorBarOn=True)
            if(sequenceLength%rowLength != 0):
                # plotGEMMEmatrix(gemmeData, args.outputfile+"_part_"+str(i+2), \
                #     (i+1)*rowLength + args.beginning, \
                #     args.end,\
                #     colorMap=args.colormap, offSet=(i+1)*rowLength + args.offset, pixelType='square',\
                #     aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
                plot1DHeatMap(dataArray, args.outputfile+"_part_"+str(i+2), beg=(i+1)*rowLength + args.beginning, \
                                end=args.end, \
                                colorMap = 'Greys', \
                                offSet=(i+1)*rowLength + args.offset, pixelType='square',\
                                interactive=False, isColorBarOn=True)
        else:
            plot1DHeatMap(dataArray, args.outputfile, beg=args.beginning, end=args.end, \
                        colorMap = 'Greys', \
                        offSet=0, pixelType='square',\
                        interactive=False, isColorBarOn=True)
        
        averageReliability = np.sum(dataArray)/len(dataArray)
        print("Average Reliability Score={:.2f}".format(averageReliability))

    elif (args.datatype.lower()=='average'):
        df = pd.read_csv(args.inputfile, delimiter=r"\s+", header=None)
        print(df)
        dataArray = (df[1].to_numpy())
        if(debug):
            print(dataArray)
        plot1DHeatMap(dataArray, args.outputfile, beg=args.beginning, end=args.end, \
                    colorMap = 'rainbow', \
                    offSet=args.offset, pixelType='rectangle',\
                    interactive=False, isColorBarOn=True)
    else:
        print("ERROR: Unknown data type: {}.".format(args.datatype))
