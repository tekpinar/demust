import os
import pandas as pd
import matplotlib.pyplot as plt

path="/home/tekpinar/research/datasets/"

def getGEMMEPaperFigure2Data(path):
    """
    Parse GEMME paper (https://doi.org/10.1093/molbev/msz179)
    official results for Figure 2&4 from the Excel file
    provided as supplementary information.

    Parameters:
    path: String    
          Paths of the reference file and the supplementary excel file.

    Returns:
    df: a pandas data frame.
        The data frame contains ncessary columns to plot Figure 2.

    """

    excelFile="supTab_allResults.xlsx"
    df = pd.read_excel(path+excelFile)
    # print(df)
    # print(df.columns)


    #Read the reference file. 
    #The reference file contains short names, reference experiments, papers, 
    #gene names etc for the Rieselman dataset. 
    referenceFile="cleaned_pdbs-v5/"+"readme_main_single_gemme_fig2_order-v9"

    #Read all the lines in the reference file
    f = open(path+referenceFile,"r")
    lines = f.readlines()
    f.close()

    rho_Combi_GEMME = []
    rho_DEEP = []
    rho_Epi_EVmut = []
    prot_nums = []
    prot_names = []
    i = 0
    for line in lines:
        words=line[:-1].split(",")
        #print(words)
        prot=words[0].strip()
        print(prot)
        expname=words[1].strip()
        print(expname)

        prot_names.append(words[4].strip().strip("\""))
        experiment=words[5].strip()
        print(experiment)

        # Get all rows, namely all experiments, performed on a particular protein. 
        allExpOfProt =df.loc[df.dataset_shortname==prot]
        
        #Append results to the lists
        singleExpOfProtGEMME = allExpOfProt.loc[df.experiment_name==experiment, 'rho_Best_GEMME'].values
        rho_Combi_GEMME.append(singleExpOfProtGEMME[0])

        singleExpOfProtDEEP = allExpOfProt.loc[df.experiment_name==experiment, 'rho_DEEP'].values
        rho_DEEP.append(singleExpOfProtDEEP[0])

        singleExpOfProtEVMut = allExpOfProt.loc[df.experiment_name==experiment, 'rho_Epi_EVmut'].values
        rho_Epi_EVmut.append(singleExpOfProtEVMut[0])
        prot_nums.append(i)
        i +=1 

    dict ={'prot_names':prot_names, 'rho_Combi_GEMME':rho_Combi_GEMME, 'rho_DEEP':rho_DEEP, 'rho_Epi_EVmut':rho_Epi_EVmut}
    dataFrame = pd.DataFrame(dict, index=prot_nums)
    print(dataFrame)
    return (dataFrame)

def plotGEMMEPaperFigure2Data(prot_nums, prot_names, \
                            rho_Combi_GEMME, rho_DEEP, rho_Epi_EVmut, mySpearmanData,\
                            title_rho_Combi_GEMME="Global (GEMME-JET (From the GEMME Paper))) ", \
                            title_rho_DEEP="Latent (DeepSequence)", \
                            title_rho_Epi_EVmut="Pairwise (EVmutation)",\
                            title_mySpearmanData="",\
                            outfile="fig2-gemme-official-with-viral.png"):

    plt.rcParams['font.size'] = '16'
    fig, ax = plt.subplots(figsize=[16, 9])

    plt.ylabel("Spearman rank correlation coefficient")

    plt.scatter(prot_nums, rho_Combi_GEMME, s=50, c="red", marker="^", zorder=1)
    plt.scatter(prot_nums, rho_DEEP, s=20, c="blue", marker="o", zorder=0)
    plt.scatter(prot_nums, rho_Epi_EVmut, s=20, c="green", marker="D", zorder=-1)
    if(mySpearmanData!=[]):
        plt.scatter(prot_nums, mySpearmanData, s=50, c="black", marker="^", zorder=2)

    if(mySpearmanData!=[]):
        plt.legend([title_rho_Combi_GEMME, title_rho_DEEP, title_rho_Epi_EVmut, title_mySpearmanData], loc='lower left')
    else:
        plt.legend([title_rho_Combi_GEMME, title_rho_DEEP, title_rho_Epi_EVmut], loc='lower left')
    plt.axvline(x=28.5, color='k', linestyle='--')
    #ax.set_xticklabels(prot_names)
    ax.set_xticks(prot_nums)

    ax.set_ylim(-0.05, 0.85)
    ax.set_xticklabels(prot_names, rotation = 45, ha='right')
    ax.text(10, 0.80, "Prokaryotic or eukaryotic")

    #Part to modify for including-excluding viral section
    plt.grid()

    plt.tight_layout()
    ax.set_xlim(-0.5, 28.5)
    #plt.savefig("fig2-alphabet-lw-i-7-wout-viral.png")
    ax.set_xlim(-0.5, 33.5)
    ax.text(29, 0.80, "Viral")
    plt.savefig(outfile)
    
if __name__ == '__main__':

    df = getGEMMEPaperFigure2Data(path)

    plotGEMMEPaperFigure2Data(df.index, df['prot_names'], \
                                df['rho_Combi_GEMME'], \
                                df['rho_DEEP'],\
                                df['rho_Epi_EVmut'],\
                                [])

