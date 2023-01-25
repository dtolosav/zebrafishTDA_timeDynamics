'''
    Program: get_barcodes_wrapper2021.py

    Date: 31 December 2021

    Purpose: Produce persistent homology barcodes (using Ripser) on simulated data given a folder of cell-cell distances with the naming conventions below

    Input:
    Folder with .txt files containing cell-cell distances
    Number indicating maximum homology dimension for Ripser computation (usually 1)
    Folder in which to place the persistent homology barcodes from Ripser
    
    Output: Persistent homology barcodes from Ripser

    Example run in terminal:
    python3 get_barcodes_wrapper2021.py -i InputFolder -d 1 -o OutputFolder

    Authors:
    Melissa R. McGuirl, Brown University, 2019
    Alexandria Volkening, Purdue University, 2021

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
    #parser.add_argument('-i', '--indir',action = 'store',required = True,  help = '''provide path to directory containing distance matrix input''')
    #parser.add_argument('-d', '--dim', action = 'store', required = True, help = '''provide maximum homology dimension for ripser computation''')
    #parser.add_argument('-o', '--outdir',action = 'store',required = True,  help = '''provide path to directory for output''')

    #distancesData = [f for f in listdir(args.indir) if isfile(join(args.indir, f))]
    distancesData = glob.glob("Distance_Matrices\sim"+ str(sim_num) + "\*\*.txt")
    #print(distancesData)
#    distancesData = []


    numFiles = len(distancesData)

    #args = parser.parse_args()
    
    os.mkdir("PD/sim" + str(sim_num))

    directory_mel = "PD/sim" + str(sim_num) + "/Mel"
    os.mkdir(directory_mel)
    directory_xanc = "PD/sim" + str(sim_num) + "/XanC"
    os.mkdir(directory_xanc)
    directory_xansn = "PD/sim" + str(sim_num) + "/XanSn"
    os.mkdir(directory_xansn)
    directory_irid = "PD/sim" + str(sim_num) + "/IriD"
    os.mkdir(directory_irid)
    directory_iril = "PD/sim" + str(sim_num) + "/IriL"
    os.mkdir(directory_iril)

    for dataj in range(numFiles):
        inFile = distancesData[dataj]
        #args.indir
        #outFile = args.outdir
        print(distancesData[dataj])
        D = np.loadtxt(inFile, dtype='float', delimiter = ',')
        #dim = int(args.dim)
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
            np.savetxt("PD/sim" + str(sim_num) + "/" + cell_type + "/PD_" + cell_type + "sim" + str(sim_num)+ time +'_dim' +  str(i), PD)    
        
if __name__=="__main__":
    main()

