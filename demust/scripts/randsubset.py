from Bio import SeqIO
import sys
import itertools
import random


def randsubset(inputMSAfile, outputMSAfile, numSel):
    """
        Select numSel sequences randomly from an inputMSAfile
        and write them to outputMSAfile. 
    """

    if(numSel == None):
        print("ERROR: Number of selected sequences can not be None!")
        print("       You have to specify number of selected sequences")
        print("       like this '-n 10' to randomly select 10 sequences!")
        sys.exit(-1)
    elif (numSel == 0):
        print("ERROR: Number of selected sequences can not be zero!")
        sys.exit(-1)
    else:
        N = int(numSel)
        print("Selecting "+str(N)+" to sequences from the input MSA!")\
            


    allRecords = SeqIO.parse(inputMSAfile, "fasta")

    totalSeqs = 500
    myRandList = random.sample(range(0, totalSeqs), numSel)

    outfile = outputMSAfile
    numRecords = 0
    my_file = open(outfile, "w")
    for record in allRecords:
        if(numRecords in myRandList):
            my_file.write(">"+str(record.id)+"\n")
            my_file.write(str(record.seq)+"\n")
        
        numRecords += 1

    my_file.close()

def randsubsetApp(args):
    """
        Select numSel sequences randomly from an inputMSAfile
        and write them to outputMSAfile. 
    """

    randsubset(args.inputmsa, args.outputmsa, args.numsel)