import gemmemore
import seaborn
import matplotlib.pyplot as plt
aaList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
def plot2DMethodComparisons(dataMatrix1, dataMatrix2, legend1, legend2, aa='all'):
    plt.figure()
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    plt.xlabel(legend1)
    plt.ylabel(legend2)
    if(aa.upper()=='ALL'):   
        plt.scatter(dataMatrix1[0], dataMatrix2[0], label="A")
        plt.scatter(dataMatrix1[1], dataMatrix2[1], label="C")
        plt.scatter(dataMatrix1[2], dataMatrix2[2], label="D")
        plt.scatter(dataMatrix1[3], dataMatrix2[3], label="E")
        plt.scatter(dataMatrix1[4], dataMatrix2[4], label="F")
        plt.scatter(dataMatrix1[5], dataMatrix2[5], label="G")
        plt.scatter(dataMatrix1[6], dataMatrix2[6], label="H")
        plt.scatter(dataMatrix1[7], dataMatrix2[7], label="I")
        plt.scatter(dataMatrix1[8], dataMatrix2[8], label="K")
        plt.scatter(dataMatrix1[9], dataMatrix2[9], label="L")
        plt.scatter(dataMatrix1[10], dataMatrix2[10], label="M")
        plt.scatter(dataMatrix1[11], dataMatrix2[11], label="N")
        plt.scatter(dataMatrix1[12], dataMatrix2[12], label="P")
        plt.scatter(dataMatrix1[13], dataMatrix2[13], label="Q")
        plt.scatter(dataMatrix1[14], dataMatrix2[14], label="R")
        plt.scatter(dataMatrix1[15], dataMatrix2[15], label="S")
        plt.scatter(dataMatrix1[16], dataMatrix2[16], label="T")
        plt.scatter(dataMatrix1[17], dataMatrix2[17], label="V")
        plt.scatter(dataMatrix1[18], dataMatrix2[18], label="W")
        plt.scatter(dataMatrix1[19], dataMatrix2[19], label="Y")

    else:
        plt.scatter(dataMatrix1[aaList.index(aa.upper())], dataMatrix2[aaList.index(aa.upper())], label=aa.upper())



    plt.legend()
    plt.show()

def plot3DMethodComparisons(dataMatrix1, dataMatrix2, dataMatrix3, legend1, legend2, legend3, aa='all'):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #plt.axhline(y=0, color='k')
    #plt.axvline(x=0, color='k')
    
    ax.set_xlabel(legend1)
    ax.set_ylabel(legend2)
    ax.set_zlabel(legend3)
    if(aa.upper()=='ALL'):      
        ax.scatter(dataMatrix1[0], dataMatrix2[0], dataMatrix2[0], label="A")
        ax.scatter(dataMatrix1[1], dataMatrix2[1], dataMatrix3[1], label="C")
        ax.scatter(dataMatrix1[2], dataMatrix2[2], dataMatrix3[2], label="D")
        ax.scatter(dataMatrix1[3], dataMatrix2[3], dataMatrix3[3], label="E")
        ax.scatter(dataMatrix1[4], dataMatrix2[4], dataMatrix3[4], label="F")
        ax.scatter(dataMatrix1[5], dataMatrix2[5], dataMatrix3[5], label="G")
        ax.scatter(dataMatrix1[6], dataMatrix2[6], dataMatrix3[6], label="H")
        ax.scatter(dataMatrix1[7], dataMatrix2[7], dataMatrix3[7], label="I")
        ax.scatter(dataMatrix1[8], dataMatrix2[8], dataMatrix3[8], label="K")
        ax.scatter(dataMatrix1[9], dataMatrix2[9], dataMatrix3[9], label="L")
        ax.scatter(dataMatrix1[10], dataMatrix2[10], dataMatrix3[10], label="M")
        ax.scatter(dataMatrix1[11], dataMatrix2[11], dataMatrix3[11], label="N")
        ax.scatter(dataMatrix1[12], dataMatrix2[12], dataMatrix3[12], label="P")
        ax.scatter(dataMatrix1[13], dataMatrix2[13], dataMatrix3[13], label="Q")
        ax.scatter(dataMatrix1[14], dataMatrix2[14], dataMatrix3[14], label="R")
        ax.scatter(dataMatrix1[15], dataMatrix2[15], dataMatrix3[15], label="S")
        ax.scatter(dataMatrix1[16], dataMatrix2[16], dataMatrix3[16], label="T")
        ax.scatter(dataMatrix1[17], dataMatrix2[17], dataMatrix3[17], label="V")
        ax.scatter(dataMatrix1[18], dataMatrix2[18], dataMatrix3[18], label="W")
        ax.scatter(dataMatrix1[19], dataMatrix2[19], dataMatrix3[19], label="Y")
    else:
        ax.scatter(dataMatrix1[aaList.index(aa.upper())], \
                   dataMatrix2[aaList.index(aa.upper())], \
                   dataMatrix3[aaList.index(aa.upper())], label=aa.upper())

    ax.legend()
    plt.show()
if (__name__ == '__main__'):
    gemmeData = gemmemore.parseGEMMEoutput("../pyrin/sequences/gemme/resid414-781/resid-414-781/normPred_evolCombi.txt")
    #print(len(gemmeData[0]))
    
    foldxData = gemmemore.parseFOLDXoutput("../pyrin/other-methods/foldX/data/PS_4cg4-bio-assembly1-chainA_just_protein-resid414-781_scanning_output.txt", colorThreshhold=7.5, colorCorrect=True)
    #print(len(foldxData[0]))


    sasaData = gemmemore.parseGEMMEoutput("/home/tekpinar/research/pyrin/other-methods/sasa/sasa-scanning-monomer-modeller-dataset.dat")
    #print(len(sasaData[0]))

    #sasaData = gemmemore.parseGEMMEoutput("/home/tekpinar/research/pyrin/other-methods/sasa/sasa-scanning-dimer-modeller-dataset.dat")
    #print(len(sasaData[0]))

    dmsCGData = gemmemore.parseGEMMEoutput("/home/tekpinar/research/pyrin/other-methods/dms-cg/entropy-scanning-1104-modes.dat")
    print((dmsCGData[0]))


    plot2DMethodComparisons(gemmeData, foldxData, "GEMME", "FoldX", aa='all')

    plot2DMethodComparisons(gemmeData, sasaData, "GEMME", "SASA-Monomer", aa='all')

    plot2DMethodComparisons(gemmeData, dmsCGData, "GEMME", "dmsCG-Monomer", aa='all')

    plot2DMethodComparisons(foldxData, dmsCGData, "FoldX", "dmsCG-Monomer", aa='all')

    plot2DMethodComparisons(foldxData, dmsCGData, "FoldX", "dmsCG-Monomer", aa='A')

    plot3DMethodComparisons(gemmeData, foldxData, sasaData, "GEMME", "FoldX", "SASA-Monomer", aa='all')
    
    #There is an extra information source which shows a linear relation between FoldX and dmsCG-Monomer
    #but I don't know the source of this information
    plot3DMethodComparisons(gemmeData, foldxData, dmsCGData, "GEMME", "FoldX", "dmsCG-Monomer", aa='alL')
