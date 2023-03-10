from demust.io import *

def diffmapApp(args):
    if (args.inputfile1 == None or args.inputfile2 == None):
        print('Usage: demust diffmap [-h] [-i INPUTFILE1] [--itype GEMME] [-j INPUTFILE2] [--jtype GEMME]')
        print('\nError: missing arguments: Please provide --inputfile1, --itype, --inputfile2, --jtype ')
        sys.exit(-1)

    if (os.path.isfile(args.inputfile1)==False  or os.path.isfile(args.inputfile2)==False):
        print("ERROR: {} or {} files do not exist in the folder!".format(args.inputfile1, args.inputfile2))
        sys.exit(-1)
    
    debug = True
    print("\nRunning 'demust diffmap' app...\n")
    if((args.itype == "gemme") and (args.jtype == "gemme")):
        dataSet1 = parseGEMMEoutput(args.inputfile1, verbose=False)
        dataSet2 = parseGEMMEoutput(args.inputfile2, verbose=False)
        if len(dataSet1)==len(dataSet2):
            diffmap = dataSet1 - dataSet2
            if(args.end == None):
                args.end = len(dataSet1[0])
        else:
            print("ERROR: Lengths of two datasets are not equal!")
            sys.exit(-1)

        if(args.paginate!=0):
            sequenceLength = (args.end - args.beginning - 1) # -1 is for starting the count from 0. 

            rowLength = args.paginate
            numberOfImageChunks = int(int(sequenceLength)/int(rowLength))
            
            for i in range(numberOfImageChunks):
                plotGEMMEmatrix(diffmap, args.outputfile+"_part_"+str(i+1), \
                    i*rowLength + args.beginning, \
                    (i+1)*rowLength + args.beginning -1,\
                    colorMap=args.colormap, offSet=i*rowLength + args.offset, pixelType='square',\
                    aaOrder=args.aaorder, sequence=args.sequence, interactive=False, isColorBarOn=args.iscolorbaron)
            if(sequenceLength%rowLength != 0):
                plotGEMMEmatrix(diffmap, args.outputfile+"_part_"+str(i+2), \
                    (i+1)*rowLength + args.beginning, \
                    args.end,\
                    colorMap=args.colormap, offSet=(i+1)*rowLength + args.offset, pixelType='square',\
                    aaOrder=args.aaorder, sequence=args.sequence, interactive=False, isColorBarOn=args.iscolorbaron)
        else:
            plotGEMMEmatrix(diffmap, args.outputfile, args.beginning, args.end,\
                colorMap=args.colormap, offSet=args.offset, pixelType='square',\
                aaOrder=args.aaorder, sequence=args.sequence, interactive=False, isColorBarOn=args.iscolorbaron)

    elif((args.itype == "singleline") and (args.jtype == "singleline")):
        with open(args.inputfile1, 'r') as file:
            allExpLines = file.readlines()

        if(debug):
            print(allExpLines)

        with open(args.inputfile2, 'r') as file:
            allCompLines = file.readlines()

        if(debug):
            print(allCompLines)

        #Do something here!
        dataSet1, dataSet2 = innerLoopV6(allExpLines, allCompLines, msafile=args.msafile, \
                                        ssfile=args.ssfile, writeOutput=True, outfile="exp-vs-comp.txt")
        #dataSet1, dataSet2 = innerLoopV2(allExpTypedList, allCompTypedList, debug=True)
        
        dataSet1 = np.array(dataSet1, dtype=float)
        dataSet2 = np.array(dataSet2, dtype=float)
    else:
        print("@> ERROR: Unknown --itype or --jtype specified.")
        print("@>        It can be gemme or singleline!")
        sys.exit(-1)

