from Bio import SeqIO
import sys
import itertools
import random


def randsubset(inputMSAfile, outputMSAfile, numSel):
    """
        Select numSel sequences randomly from an inputMSAfile
        in fasta format and write them to outputMSAfile in fasta
        format. 
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

    allRecords = list(SeqIO.parse(inputMSAfile, "fasta"))

    totalSeqs = len(allRecords)
    if(totalSeqs > N):
        print("There are "+str(totalSeqs)+" sequences in the input MSA!")
        print("Selecting "+str(N)+" sequences from the input MSA!")
    else:
        print("ERROR: There are "+str(totalSeqs)+" sequences in the input MSA")
        print("       but you are trying to select " + str(N) + " sequences!!")
        print("Please make numsel smaller than "+str(totalSeqs)+"!")
        sys.exit(-1)
        
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