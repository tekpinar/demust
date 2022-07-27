from demust.io import *
import argparse

def getMinMaxData(scanningMatrix, type, outfile, printDetails=False):
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
    FILE = open(outfile, "w")
    if(type.lower()=='max'):
        for col in range(len(scanningMatrix[0])):
            X_max = (scanningMatrix.T[col].max())            
            FILE.write("{}\t{}\n".format(col, X_max))
            data.append(X_max)
            
            if(printDetails):
                print(scanningMatrix.T[col])
        
    elif(type.lower()=='min'):
        for col in range(len(dataArray)):
            X_min = (scanningMatrix.T[col].min())
            FILE.write("{}\t{}\n".format(col, X_min))
            data.append(X_min)

            if(printDetails):
                print(scanningMatrix.T[col])
    else:
        print("ERROR: Unknown type!")
        print("       Type can only be min or max!")
        sys.exit(-1)
    
    FILE.close()
    return data

def plotMinMax(dataArray, outFile, beg, end, \
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
        end = len(dataArray)

    print("Beginning: "+str(beg))
    print("End      : "+str(end))
    
    print(len(dataArray))
    subArray = dataArray[(beg-1):end]
    #print(subArray)

    ##########################################################################
    # Set plotting parameters
    nres_shown = len(subArray)
    fig_height=8
    # figure proportions
    fig_width = fig_height/2  # inches
    fig_width *= nres_shown/20
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

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

    major_nums_y = np.arange(0, 20, 1, dtype=int)
    major_labels_y = list(aaOrder)

    plt.xticks(major_nums_x, major_labels_x, size=28)
    ax.set_xticks(minor_nums_x, minor=True)
    
    plt.yticks(major_nums_y, major_labels_y, size=16)
    ax.set_yticklabels(major_labels_y, ha='left')
    ax.tick_params(axis='y', which='major', pad=30)

    
    #############################################################################
    if(pixelType=='square'):
        #For plotting square pixels
        img = plt.bar(subArray)
    #elif(pixelType=='rectangle'):
        #For plotting rectangular pixels
        #img = plt.imshow(subArray, cmap=colorMap, aspect=3.0)
    else:
        print("\nERROR: Unknown pixelType specified!\n")
        sys.exit(-1)
    
    #To make the colors consistent if there are submatrices.
    #plt.clim(np.min(scanningMatrix), np.max(scanningMatrix)) 

    # if(sequence!=None):
    #     mySeqFile = SeqIO.read(sequence, 'fasta')
    #     #Convert aaOrder to a list.
    #     aaOrderList = list(aaOrder)
    #     for i in range (len(subArray)):
    #         j = beg-1+i
    #         # print(i, aaOrderList.index(sequence[i]))
    #         plt.scatter(i, aaOrderList.index(mySeqFile.seq[j]), s=5, c='black', marker='o')
    
    # if(isColorBarOn):
    #     from mpl_toolkits.axes_grid1 import make_axes_locatable
    #     divider = make_axes_locatable(ax)
    #     cax = divider.append_axes("right", size="1%", pad=0.2)
    #     plt.colorbar(img, cax=cax)

    plt.tight_layout()
    plt.savefig(outFile+".png")
    if(interactive):
        plt.show()
    

    #plt.imsave('output.png', subArray)
    plt.close()


def plotsApp(args):
    if args.inputfile == None:
        print('Usage: python demust.py [-h] [-i INPUTFILE] [-d GEMME]')
        print('\nError: missing arguments: Please provide --inputfile and/or --datatype')
        sys.exit(-1)

    print("\nRunning 'demust plots' app...\n")
    if (args.datatype.lower()=='gemme'): 
        scanningMatrix = parseGEMMEoutput(args.inputfile, verbose=False)

    if((args.type.lower()=='min') or (args.type.lower()=='max')):
        dataArray = getMinMaxData(scanningMatrix, args.type, outfile=args.outputfile, printDetails=False)

    else:
        print("ERROR: Unknown plot type!{}".format(args.type.lower()))
        sys.exit(-1)
    
