# This is a short scrypt to plot birth/death diagrams, having as input the persistence data that ripser returns. For example on 
# the files that the function get_barcodes_wrapper2021.py creates from the distance matrices.

import numpy as np
import matplotlib.pyplot as plt
from ripser import ripser
from persim import plot_diagrams

# read the file
f = np.genfromtxt(fname='OutputTest1\mel_distances_1_time51_dim1')

# plot

#plot_diagrams(f, show=True)
#plot_diagrams(f, show=True, lifetime = True)
fig, ax = plt.subplots()

ax.stem(f[], y)

ax.set(xlim=(0, 8), xticks=np.arange(1, 8),
       ylim=(0, 8), yticks=np.arange(1, 8))