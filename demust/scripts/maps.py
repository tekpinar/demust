import os
import string
import sys
import argparse
import numpy as np
import matplotlib.pylab as plt
from demust.io import *

def mapsApp(args):

    if args.inputfile == None:
        print('Usage: demust maps [-h] [-i INPUTFILE] [-d GEMME]')
        print('\nError: missing arguments: Please provide --inputfile and/or --datatype')
        sys.exit(-1)

    print("\nRunning 'demust maps' app...\n")
        
    if (args.datatype.lower()=='gemme'):
        gemmeData = parseGEMMEoutput(args.inputfile, verbose=False)
        
        if(args.ranknorm):
            gemmeData = rankNormalization(gemmeData)
            plotGEMMEmatrix(gemmeData, args.outputfile, args.beginning, args.end,\
                colorMap=args.colormap, offSet=args.offset, pixelType='square',\
                aaOrder=args.aaorder)
        else:
            plotGEMMEmatrix(gemmeData, args.outputfile, args.beginning, args.end,\
                colorMap=args.colormap, offSet=args.offset, pixelType='square',\
                aaOrder=args.aaorder)

    elif (args.datatype.lower()==('rhapsody' or 'rapsody')):
        rhapsodyData = parseRHAPSODYoutput(args.inputfile, field=args.field)
        plotGEMMEmatrix(rhapsodyData, args.outputfile, args.beginning, args.end,\
            colorMap='coolwarm', offSet=args.offset, pixelType='square',\
            aaOrder=args.aaorder)
        
    elif (args.datatype.lower()=='foldx'):
        foldxData = parseFOLDXoutput(args.inputfile, colorThreshhold=7.5, colorCorrect=True)
        plotGEMMEmatrix(foldxData, args.outputfile, args.beginning, args.end,\
            colorMap='coolwarm', offSet=args.offset, pixelType='square',\
            aaOrder=args.aaorder)
    else:
        print("\nError: Unknown data type!")
        print("         Data types can only be gemme, rhapsody or foldx!")
    # gemmeAverages = getGEMMEAverages(gemmeData.T)
    
    # #Read PRS data
    # prsData = parsePRSdata(sys.argv[2])

    # plotGEMME2DvsOther(gemmeAverages, prsData, \
    #                     customColor=sys.argv[3], \
    #                     figName=sys.argv[4])
