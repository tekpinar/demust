import sys
# from prody import *
import numpy as np
# from matplotlib.pylab import *
# from sedy.postprocess import *
import pandas as pd
#This script was part of sedy package.
from scipy.stats import rankdata,zscore
import matplotlib.pyplot as plt
#Part of Sedy Package
def rankSortData(dataArray):
    """
        This function ranksorts protein data and converts it
        to values between [0,1.0].
    
    Parameters:
    ----------
    dataArray: numpy array of arrays
               data read by numpy.loadtxt

    Returns:
    -------
    normalizedRankedDataArray: numpy array
    """
    normalizedRankedDataArray = rankdata(dataArray)/float(len(dataArray))

    return (normalizedRankedDataArray)

def minMaxNormalization(data):
    """
        Min-max normalization of a data array.

    Parameters:
    ----------
    data: numpy array
          

    Returns:
    -------
    normalizeddata: numpy array
    """
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def zscore(data):
    """
        Z-score normalization

    Parameters:
    ----------
    data: numpy array
          

    Returns:
    -------
    normalizeddata: numpy array

    """
    return (data - np.mean(data))/np.std(data)

def attenuateEndPoints(data):                                                                                                    
    """                                                                                                                         
        Reduce signal strength of the signal on N and C terminals
        with a sigmoid function.                                                                
    """                                                                                                                         
                                                                                                                                
    debug = 0                                                                  
                                                                               
    n = len(data)                                                              
    if(debug):                                                                 
        #Check array lengths                                                   
        print(n/2)                                                             
        print(len(data))

    if(n%2 == 0):                                                              
        x1 = np.linspace(-0, 10, int(n/2))                                     
        x2 = np.linspace(-10, 0, int(n/2))                                     
    else:                                                                      
        x1 = np.linspace(-0, 10, int((n-1)/2))                                 
        x2 = np.linspace(-10, 0, int((n+1)/2))                                 
                                                                               
    z1 = 1/(1 + np.exp(-2.50*x1))                                              
    z2 = 1.0 - (1/(1 + np.exp(-2.50*x2)))                                      
                                                                               
                                                                               
    z = np.append(z1, z2)
    if(debug):                                                                 
        #Check array lengths                                                   
        print(len(z1))                                                     
        print(len(z2))                                                         
        print(len(z))                                                          
        #Plot the arrays                                                       
        plt.plot(z)                                                            
        plt.plot(data)                                                         
        plt.plot(data*z)                                                       
        plt.xlabel("x")                                                        
        plt.ylabel("Sigmoid(X)")                                               
        plt.show()

    if(len(data) == len(z)):                                                   
        return (data*z)                                                                                                          
    else:                                                                       
        print("ERROR: Can not attenuate N and C terminal data signals!")        
        sys.exit(-1)

def postprocessApp(args):
    #Mostyl, I am using normPred_Combi_singleline as input file and it doesn't have a header.
    df = pd.read_table(args.input, sep="\s+", header=None)

    #data = np.genfromtxt(args.input,dtype=None)
    data = df.to_numpy()

    # print(data)
    rawData = data.T[args.column-1]
    # print(rawData)
    # processedData = rankdata(rawData.T[args.column - 1])/float(len(rawData.T[args.column - 1]))

    if(args.process == 'ranksort'):
        processedData = rankSortData(rawData)
    elif(args.process == '1-ranksort'):
        processedData = 1.0 - rankSortData(rawData)
    elif(args.process == 'minmax'):
        processedData = minMaxNormalization(rawData)
    elif(args.process == '1-minmax'):
        processedData = 1.0 - minMaxNormalization(rawData)
    elif(args.process == '1-values'):
        processedData = 1.0 - rawData
    elif(args.process == 'zscore'):
        processedData = zscore(rawData)
    elif(args.process == 'negzscore'):
        processedData = -1.0*zscore(rawData)
    elif(args.process == 'attenuate'):
        processedData = attenuateEndPoints(rawData)
    elif(args.process == None):
        print("@> No postprocessing applied!")
    else:
        print("@> Unknown postprocessing option!")
        sys.exit(-1)

    print("@> Processing the data started!")
    with open(args.outfile, 'w') as f:
        #f.write("#Resid Value\n")
        for i in range (len(processedData)):
            f.write("{:} {:6.2f}\n".format(data.T[0][i], processedData[i]))
    print("@> Processing the data finished successfully!")
    sys.exit()

