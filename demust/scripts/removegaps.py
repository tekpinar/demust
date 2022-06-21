import argparse

import numpy as np
import pandas as pd

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq

def removeGaps(inputMSAfile, outputMSAfile):
    """
        Remove gaps in the query sequence!

    Parameters
    ----------
    inputMSAfile: str
        Name of the input multiple sequence alignment file in fasta format.
    
    outputMSAfile: str
        Name of the output multiple sequence alignment file in fasta format.  
    Returns
    -------
    Nothing
    """
    # Build a dictionary of original sequence index vs aligned sequence index. 
    # Let's read the MSA file first
    alignment = AlignIO.read(inputMSAfile, 'fasta')
    print(alignment[0].seq)
    numSeqInMSA = len(alignment)

    df = pd.DataFrame()
    for i in range (numSeqInMSA):
        df[alignment[i].id] = list(alignment[i].seq)

    dfTransposed = (df.T)
    print(len(dfTransposed.columns))
    
    myEmptyList = []
    for j in range(len(dfTransposed.columns)):
        if(dfTransposed.iat[0, j] != '-'):
            myEmptyList.append(list(dfTransposed.iloc[:, j]))

    print(len(myEmptyList))
    print(len(myEmptyList[0]))   

    numpy_array = np.array(myEmptyList)
    transpose = numpy_array.T
    transpose_list = transpose.tolist()

    for i in range (numSeqInMSA):
        alignment[i].seq = Seq("".join(transpose_list[i]))

    print(alignment[1].seq)
    with open(outputMSAfile, "w") as handle:
        count = SeqIO.write(alignment, handle, "fasta")

def removeGapsApp(args):
    removegaps_parser = argparse.ArgumentParser()

    removegaps_parser.add_argument('-i', '--inputgappedmsa', dest='inputgappedmsa', type=str, \
        help="Name of the input gapped Multiple Sequence Alignment file in fasta format.",
        required=True, default=None)

    removegaps_parser.add_argument('-o', '--outputungappedmsa', dest='outputungappedmsa', type=str, \
        help="Name of the output gapped Multiple Sequence Alignment file in fasta format. Default is outputungappedmsa.fasta",
        required=False, default="outputungappedmsa.fasta")

    args = removegaps_parser.parse_args() 
    # Check flag arguments
    if args.inputgappedmsa is None:
    	print("ERROR: The program requires an input MSA in fasta format")

    removeGaps(args.inputgappedmsa, args.outputungappedmsa)