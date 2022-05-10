#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys

plt.rcParams['font.family'   ] = 'Times New Roman'
plt.rcParams['axes.linewidth'] = 0.5

def load_input_file( filename ):
    infile   = open( filename, 'r' )
    dataline = infile.readlines()
    infile.close()

    matA = []
    ncol = 0
    nrow = 0

    for eachline in dataline:
        eachline.strip()
        item = eachline.split()
        
        for i in range( len(item) ):
            if nrow==1:
                ncol += 1
            matA.append( float(item[i]) )

        nrow += 1

    matA = np.asarray( matA )
    matA = np.reshape( matA, (nrow,ncol) )

    return matA

def plot_spymat( data ):
    fig, ax = plt.subplots( 1, 1, figsize=(6,2) )
    ax.imshow( data, interpolation='none', cmap='jet' )

    plt.tight_layout()
    plt.show()


###########################################
if __name__=='__main__':
    data = load_input_file( sys.argv[1] )
    plot_spymat( data )

