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
        if(args.end == None):
            args.end = len(gemmeData[0])
        
        if(args.ranknorm):
            gemmeData = rankNormalization(gemmeData)
            plotGEMMEmatrix(gemmeData, args.outputfile, args.beginning, args.end,\
                colorMap=args.colormap, offSet=args.offset, pixelType='square',\
                aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
        else:
            if(args.paginate!=0):
                sequenceLength = (args.end - args.beginning - 1) # -1 is for starting the count from 0. 

                rowLength = args.paginate
                numberOfImageChunks = int(int(sequenceLength)/int(rowLength))
                
                for i in range(numberOfImageChunks):
                    plotGEMMEmatrix(gemmeData, args.outputfile+"_part_"+str(i+1), \
                        i*rowLength + args.beginning, \
                        (i+1)*rowLength + args.beginning -1,\
                        colorMap=args.colormap, offSet=i*rowLength + args.offset, pixelType='square',\
                        aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
                if(sequenceLength%rowLength != 0):
                    plotGEMMEmatrix(gemmeData, args.outputfile+"_part_"+str(i+2), \
                        (i+1)*rowLength + args.beginning, \
                        args.end,\
                        colorMap=args.colormap, offSet=(i+1)*rowLength + args.offset, pixelType='square',\
                        aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
            else:
                plotGEMMEmatrix(gemmeData, args.outputfile, args.beginning, args.end,\
                    colorMap=args.colormap, offSet=args.offset, pixelType='square',\
                    aaOrder=args.aaorder, sequence=args.sequence, interactive=False, isColorBarOn=args.iscolorbaron)
                # plotDecoratedMatrix2(gemmeData, args.outputfile, args.beginning, args.end,\
                #     colorMap=args.colormap, offSet=args.offset, pixelType='square',\
                #     aaOrder=args.aaorder, sequence=args.sequence, interactive=False, isColorBarOn=args.iscolorbaron)

    elif (args.datatype.lower()==('rhapsody' or 'rapsody')):
        rhapsodyData = parseRHAPSODYoutput(args.inputfile, field=args.field)
        if(args.end == None):
            args.end = len(rhapsodyData[0])

        if(args.paginate!=0):
            sequenceLength = (args.end - args.beginning - 1) # -1 is for starting the count from 0. 

            rowLength = args.paginate
            numberOfImageChunks = int(int(sequenceLength)/int(rowLength))
            
            for i in range(numberOfImageChunks):
                plotGEMMEmatrix(rhapsodyData, args.outputfile+"_part_"+str(i+1), \
                    i*rowLength + args.beginning, \
                    (i+1)*rowLength + args.beginning -1,\
                    colorMap=args.colormap, offSet=i*rowLength + args.offset, pixelType='square',\
                    aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
            if(sequenceLength%rowLength != 0):
                plotGEMMEmatrix(rhapsodyData, args.outputfile+"_part_"+str(i+2), \
                    (i+1)*rowLength + args.beginning, \
                    args.end,\
                    colorMap=args.colormap, offSet=(i+1)*rowLength + args.offset, pixelType='square',\
                    aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
        else:        
            plotGEMMEmatrix(rhapsodyData, args.outputfile, args.beginning, args.end,\
                colorMap=args.colormap, offSet=args.offset, pixelType='square',\
                aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
        
    elif (args.datatype.lower()=='foldx'):
        foldxData = parseFOLDXoutput(args.inputfile, colorThreshhold=7.5, colorCorrect=True)
        if(args.end == None):
            args.end = len(foldxData[0])

        if(args.paginate!=0):
            sequenceLength = (args.end - args.beginning - 1) # -1 is for starting the count from 0. 

            rowLength = args.paginate
            numberOfImageChunks = int(int(sequenceLength)/int(rowLength))
            
            for i in range(numberOfImageChunks):
                plotGEMMEmatrix(foldxData, args.outputfile+"_part_"+str(i+1), \
                    i*rowLength + args.beginning, \
                    (i+1)*rowLength + args.beginning -1,\
                    colorMap=args.colormap, offSet=i*rowLength + args.offset, pixelType='square',\
                    aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
            if(sequenceLength%rowLength != 0):
                plotGEMMEmatrix(foldxData, args.outputfile+"_part_"+str(i+2), \
                    (i+1)*rowLength + args.beginning, \
                    args.end,\
                    colorMap=args.colormap, offSet=(i+1)*rowLength + args.offset, pixelType='square',\
                    aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
        else:        
            plotGEMMEmatrix(foldxData, args.outputfile, args.beginning, args.end,\
                colorMap=args.colormap, offSet=args.offset, pixelType='square',\
                aaOrder=args.aaorder, sequence=args.sequence, isColorBarOn=args.iscolorbaron)
    elif (args.datatype.lower()=='proteingym'):
        expDataMatrix = parseExperimentalData(args.inputfile, \
                                            experiment="DMS_score",\
                                            outputcsv=None,\
                                            debug = False,\
                                            datasource="proteingym")
        # Plot the experimental DMS map with wild-type residues annotated with dots.
        plotExperimentalMatrix(expDataMatrix, args.outputfile, \
                               args.beginning, args.end, \
                               colorMap = args.colormap, offSet=args.offset, \
                               pixelType='square', aaOrder=args.aaorder, \
                               sequence=args.sequence, interactive=False, isColorBarOn=args.iscolorbaron)
    elif (args.datatype.lower()=='singleline'):

        expDataMatrix = parseSingleLineData(args.inputfile, \
                                            experiment="DMS_score",\
                                            outputcsv=None,\
                                            debug = False)
        if(args.end == None):
            args.end = len(expDataMatrix[0])
        
        if(args.paginate!=0):
            sequenceLength = (args.end - args.beginning - 1) # -1 is for starting the count from 0. 

            rowLength = args.paginate
            numberOfImageChunks = int(int(sequenceLength)/int(rowLength))
            for i in range(numberOfImageChunks):
                plotExperimentalMatrix(expDataMatrix, args.outputfile+"_part_"+str(i+1), \
                    i*rowLength + args.beginning, \
                    (i+1)*rowLength + args.beginning -1,\
                    colorMap=args.colormap, offSet=i*rowLength + args.offset, pixelType='square',\
                    aaOrder=args.aaorder, sequence=args.sequence, interactive=False, isColorBarOn=args.iscolorbaron)
            if(sequenceLength%rowLength != 0):
                plotExperimentalMatrix(expDataMatrix, args.outputfile+"_part_"+str(i+2), \
                    (i+1)*rowLength + args.beginning, \
                    args.end,\
                    colorMap=args.colormap, offSet=(i+1)*rowLength + args.offset, pixelType='square',\
                    aaOrder=args.aaorder, sequence=args.sequence, interactive=False, isColorBarOn=args.iscolorbaron)


        else:    
            # Plot the experimental DMS map with wild-type residues annotated with dots.
            plotExperimentalMatrix(expDataMatrix, args.outputfile, \
                                args.beginning, args.end, \
                                colorMap = args.colormap, offSet=args.offset, \
                                pixelType='square', aaOrder=args.aaorder, \
                                sequence=args.sequence, interactive=False, isColorBarOn=args.iscolorbaron)

    else:
        print("\nError: Unknown data type!")
        print("         Data types can only be gemme, rhapsody or foldx!")
    # gemmeAverages = getGEMMEAverages(gemmeData.T)
    
    # #Read PRS data
    # prsData = parsePRSdata(sys.argv[2])

    # plotGEMME2DvsOther(gemmeAverages, prsData, \
    #                     customColor=sys.argv[3], \
    #                     figName=sys.argv[4])
