"""
    Parse experimental deep mutational scanning maps.
"""
import sys
import argparse
import pandas as pd
import numpy as np
from demust.io import *
from demust.scripts.compare import compareMapsSpearman
import matplotlib.pylab as plt
#from Bio import SeqIO
#import seaborn as sns

# alphabeticalAminoAcidsList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
#                               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def experimental2singleline(inputcsv, experiment='DMS_score', 
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
        

    header = ["mutant", experiment]
    df.to_csv(outputcsv, columns = header, index=None, header=None, sep=' ')
    print(df[['mutant', experiment]])

def experimental2mutfile(inputcsv, experiment='DMS_score', 
                            output="converted.mut",
                            debug = False, datasource="proteingym"):
    """
        Parse Riesselman or ProteinGym experimental data files from the
        https://www.nature.com/articles/s41592-018-0138-4 and 
        https://doi.org/10.48550/arXiv.2205.13760
        are converted to mutation files where each single point mutation is
        given in a single line. If there are multiple mutations separated by
        columns (:), they also can be handled. The output 
        mut file is used as input in the GEMME. 

    Parameters
    ----------
    inputcsv: string
        Name of the experimental input data file in csv format
    
    experiment: string
        Name of the experiment such as score, screenscore, fitness,
        abundance. To parse ProteinGym data, set this value to 
        'DMS_score'.

    output: string
        Name of the output mut file.

    debug: bool
        If True, it will print some extra data to stdout. 
    
    datasource: string
        It should take 'proteingym' value if you want to parse
        ProteinGym data. Otherwise, don't use it. In this case,
        it assumes that the datasource is Riesseman, 2016 csv files.

    Returns
    -------
    Nothing: a mut file is written.  
    """
    if(datasource == "proteingym"):
        df = pd.read_csv(inputcsv)
    else:
        df = pd.read_csv(inputcsv, sep=";", decimal=",")
        

    header = ["mutant", experiment]
    print(df.columns)
    df.to_csv(output, columns=["mutant"], index=None, header=None, sep=' ')
    #print(df[['mutant', experiment]])

def experimental2mutfileV2(inputcsv, shift=None, experiment='DMS_score', 
                            output="converted.mut",
                            debug = False, datasource="proteingym"):
    """
        Parse Riesselman or ProteinGym experimental data files from the
        https://www.nature.com/articles/s41592-018-0138-4 and 
        https://doi.org/10.48550/arXiv.2205.13760
        are converted to mutation files where each single point mutation is
        given in a single line. If there are multiple mutations separated by
        columns (:), they also can be handled. The output 
        mut file is used as input in the GEMME. 

    Parameters
    ----------
    inputcsv: string
        Name of the experimental input data file in csv format
    shift: integer
        An integer to shift residue numbering if number of amino acids in
        the fasta file and experimental mutation list does not match.
    experiment: string
        Name of the experiment such as score, screenscore, fitness,
        abundance. To parse ProteinGym data, set this value to 
        'DMS_score'.
    output: string
        Name of the output mut file.
    debug: bool
        If True, it will print some extra data to stdout. 
    datasource: string
        It should take 'proteingym' value if you want to parse
        ProteinGym data. Otherwise, don't use it. In this case,
        it assumes that the datasource is Riesseman, 2016 csv files.

    Returns
    -------
    Nothing: a mut file is written.  
    """

    print("Creating a mutations file from the experimental data.")
    if(datasource == "proteingym"):
        df = pd.read_csv(inputcsv)
    else:
        df = pd.read_csv(inputcsv, sep=";", decimal=",")
        
###############################################################################
    # Fill the
    if(shift == None):
        if(debug):
            print(df.columns)
        df.to_csv(output, columns=["mutant"], index=None, header=None, sep=' ')
        #print(df[['mutant', experiment]])
        print("The mutations file was produced successfully!")
    else:
        df['mutantModif']=""
        for index, row in df.iterrows():
            tempString = row['mutant'].split(',|:')
            if(len(tempString)==0):
                print("ERROR: Your input file seems to be empty!")
                sys.exit(-1)
            elif(len(tempString)==1):
                temp2 = int(tempString[0][1:-1]) - shift + 1
                #print(index, row['mutant'], row['DMS_score'])
                mutantModif = tempString[0][0] + str(temp2) + tempString[0][-1]
                df.at[index, 'mutantModif'] = mutantModif
            else:
                print("ERROR: Your input file contains multiple point mutations!")
                sys.exit(-1)
                mutantModif = ""
                for i in range (len(tempString)):
                    temp2 = int(tempString[i][1:-1]) - shift + 1
                    #print(index, row['mutant'], row['DMS_score'])
                    if(i == (len(tempString) - 1)):
                        mutantModif = mutantModif + tempString[i][0] + str(temp2) + tempString[i][-1]
                    else:
                        mutantModif = mutantModif + tempString[i][0] + str(temp2) + tempString[i][-1] + ":"
                
                df.at[index, 'mutantModif'] = mutantModif
        if(debug):
            print(df.columns)
        df.to_csv(output, columns=["mutantModif"], index=None, header=None, sep=' ')
        #print(df[['mutant', experiment]])
        print("The mutations file was produced successfully!")

def experimental2singlelineV2(inputcsv, shift=None, experiment='DMS_score', 
                            output="converted.mut",
                            debug = False, datasource="proteingym"):
    """
        Parse Riesselman or ProteinGym experimental data files from the
        https://www.nature.com/articles/s41592-018-0138-4 and 
        https://doi.org/10.48550/arXiv.2205.13760
        are converted to data files where each single point mutation is
        written on a single line. If there are multiple mutations separated by
        columns (:), they also can be handled. The output 
        dat file is used as input in the demust compare for Spearman correlation
        comparison. 

    Parameters
    ----------
    inputcsv: string
        Name of the experimental input data file in csv format
    shift: integer
        An integer to shift residue numbering if number of amino acids in
        the fasta file and experimental mutation list does not match.
    experiment: string
        Name of the experiment such as score, screenscore, fitness,
        abundance. To parse ProteinGym data, set this value to 
        'DMS_score'.
    output: string
        Name of the output mut file.
    debug: bool
        If True, it will print some extra data to stdout. 
    datasource: string
        It should take 'proteingym' value if you want to parse
        ProteinGym data. Otherwise, don't use it. In this case,
        it assumes that the datasource is Riesseman, 2016 csv files.

    Returns
    -------
    Nothing: a mut file is written.  
    """

    print("Creating a single column data file from the experimental data.")
    if(datasource == "proteingym"):
        df = pd.read_csv(inputcsv)
    else:
        df = pd.read_csv(inputcsv, sep=";", decimal=",")
        
###############################################################################
#Add a new column to the dataframe.
    # Fill the
    if(shift == None):
        if(debug):
            print(df.columns)
        df.to_csv(output, columns=["mutant", experiment], index=None, header=None, sep=' ')
        #print(df[['mutant', experiment]])
        print("The single column experimental data file was produced successfully!")
    else:
        df['mutantModif']=""
        for index, row in df.iterrows():
            tempString = row['mutant'].split(',|:')
            if(len(tempString)==0):
                print("ERROR: Your input file seems to be empty!")
                sys.exit(-1)
            elif(len(tempString)==1):
                temp2 = int(tempString[0][1:-1]) - shift + 1
                #print(index, row['mutant'], row['DMS_score'])
                mutantModif = tempString[0][0] + str(temp2) + tempString[0][-1]
                df.at[index, 'mutantModif'] = mutantModif
            else:
                print("ERROR: Your input file contains multiple point mutations!")
                sys.exit(-1)
                mutantModif = ""
                for i in range (len(tempString)):
                    temp2 = int(tempString[i][1:-1]) - shift + 1
                    #print(index, row['mutant'], row['DMS_score'])
                    if(i == (len(tempString) - 1)):
                        mutantModif = mutantModif + tempString[i][0] + str(temp2) + tempString[i][-1]
                    else:
                        mutantModif = mutantModif + tempString[i][0] + str(temp2) + tempString[i][-1] + ":"
                
                df.at[index, 'mutantModif'] = mutantModif
        if(debug):
            print(df.columns)
        df.to_csv(output, columns=["mutantModif", experiment], index=None, header=None, sep=' ')
        #print(df[['mutant', experiment]])
        print("The single column experimental data file was produced successfully!")


def riesselmanApp(args):

    # riesselman_parser = argparse.ArgumentParser(description=\
    #     "This script will parse and plot experimental DMS maps from Riesselman, 2016.")

    # # # main script argument parsing
    # riesselman_parser.add_argument('-d', '--dataset', dest='dataset', type=str, \
    #     help='Name of the dataset that you want to retrieve.', \
    #     required=True, default=None)
    
    # riesselman_parser.add_argument('-e', '--experiment', dest='experiment', type=str, \
    #     help='Experiment type field such as fitness, abundance, screenscore.', \
    #     required=True, default="DMS_score")

    # riesselman_parser.add_argument('--colormap', dest='colormap', type=str, \
    #     help='A colormap as defined in matplotlib',
    #     required=False, default='coolwarm_r')

    # riesselman_parser.add_argument('-s', '--sequence', dest='sequence', type=str, \
    #     help='One of the input sequence file in fasta format', \
    #     required=False, default=None)

    # riesselman_parser.add_argument('-m', '--metric', dest='metric', type=str, \
    #     help='Comparison metric.\n It can be spearman or pearson. Default is spearman', \
    #     required=False, default="spearman")
    
    # riesselman_parser.add_argument('-o', '--output', dest='output', type=str, \
    #     help='Name of the output file prefix for the csv and the png file. Default is output', \
    #     required=False, default='output')

    # args = riesselman_parser.parse_args()



    # #Read the main file that contains information about all datasets. 
    # I decided to skip that part! 


    

    if((args.otype == 'singleline') or (args.otype == 'dat')):
        # Read the experimental DMS map    
        experimental2singlelineV2(args.dataset, shift=args.shift, \
                                        experiment=args.experiment, \
                                        output=args.output+".csv", \
                                        debug=False, \
                                        datasource=args.source)
    elif(args.otype == 'mut'):
        # Read the experimental DMS map    
        experimental2mutfileV2(args.dataset, shift=args.shift, \
                                        experiment=args.experiment, \
                                        output=args.output+".mut", \
                                        debug=False, \
                                        datasource=args.source)

    else:
        # Read the experimental DMS map
        scanningMatrixExp = parseExperimentalData(args.dataset,\
                                            experiment=args.experiment, \
                                            outputcsv=args.output+".csv",\
                                            debug=False, \
                                            datasource="riesselman")
        print(len(scanningMatrixExp[1]))
        # Plot the experimental DMS map with wild-type residues annotated with dots.
        plotExperimentalMatrix(scanningMatrixExp, outFile=args.output, beg=1, \
                                end=len(scanningMatrixExp[1]), \
                                colorMap = args.colormap, \
                                sequence=args.sequence, pixelType='square',
                                interactive=False)
    sys.exit()
    # Read the GEMME DMS map
    scanningMatrixGEMME = parseGEMMEoutput(path+args.protein.upper()+"/CALM1_normPred_evolCombi.txt", verbose=False)

    # Remove the first column in the GEMME map because the experimental maps starts from the 
    # second amino acid. 
    # print("Length of the 0th line before:")
    # print(len(scanningMatrixGEMME[0]))
    scanningMatrixGEMMEPruned = np.delete(scanningMatrixGEMME, 0, axis=1)
    # print("Length of the 0th line after:")
    # print(len(scanningMatrixGEMMEPruned[0]))

    # Plot the pruned GEMME DMS map
    plotGEMMEmatrix(scanningMatrixGEMMEPruned, outFile=args.protein.upper()+"-gemme", beg=1, end=len(scanningMatrixGEMMEPruned[1]), \
                   colorMap = 'coolwarm_r', offSet=1, pixelType='square', sequence=expSequence.seq)

    correlation, pvalue = compareMapsSpearman(scanningMatrixExp, scanningMatrixGEMMEPruned)
    print("\nSpearman correlation={:.3f} - pvalue={:.3E}\n".format(correlation, pvalue))

