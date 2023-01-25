'''
    Program:

    Date: 

    Purpose: Produce persistent homology barcodes (using Ripser) on simulated data given a folder of cell-cell distances with the naming conventions below

    Input:
    Folder with .txt files containing cell-cell distances
    Number indicating maximum homology dimension for Ripser computation (usually 1)
    Folder in which to place the persistent homology barcodes from Ripser
    
    Output: Persistent homology barcodes from Ripser

    Example run in terminal:
    

    Authors:
    Melissa R. McGuirl, Brown University, 2019
    Alexandria Volkening, Purdue University, 2021
    Daniel Tolosa, Purdue University, 2022
'''
from ripser import ripser
import argparse,re
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import glob
def main(sim_num):
    descriptor = "calls ripser to compute persistent homology"
    parser = argparse.ArgumentParser(description = descriptor)
    distancesData = glob.glob("Distance_Matrices/WT/sim"+ str(sim_num) + "/*/*.txt")
    numFiles = len(distancesData)
    os.mkdir("barcodes/sim" + str(sim_num))
    directory_mel = "barcodes/sim" + str(sim_num) + "/Mel"
    os.mkdir(directory_mel)
    directory_xanc = "barcodes/sim" + str(sim_num) + "/XanC"
    os.mkdir(directory_xanc)
    directory_xansn = "barcodes/sim" + str(sim_num) + "/XanSn"
    os.mkdir(directory_xansn)
    directory_irid = "barcodes/sim" + str(sim_num) + "/IriD"
    os.mkdir(directory_irid)
    directory_iril = "barcodes/sim" + str(sim_num) + "/IriL"
    os.mkdir(directory_iril)

    for dataj in range(numFiles):
        inFile = distancesData[dataj]
        D = np.loadtxt(inFile, dtype='float', delimiter = ',')
        dim = 1
        try: 
            if D == 0:
                continue
        except: pass
        bars = ripser(D,  distance_matrix=True, maxdim=dim)
        cell_type = distancesData[dataj].split("_")[3]
        time = re.split('_|\.', distancesData[dataj])[5]
        for i in range(dim + 1): 
            PD = bars['dgms'][i]
            np.savetxt("barcodes/sim" + str(sim_num) + "/" + cell_type + "/BC_" + cell_type + "sim" + str(sim_num)+ time +'_dim' +  str(i), PD)    
        
if __name__=="__main__":
    main()