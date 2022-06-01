#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import sys

###################################################
plt.style.use('classic')
plt.rcParams['font.family'   ] = 'Times New Roman'
plt.rcParams['axes.linewidth'] = 0.5
###################################################

def load_result_file( filename ):
    infile    = open( filename, 'r' )
    datalines = infile.readlines()
    infile.close()

    Bdat  = {'x':[], 'y':[]}
    Idat  = {'flag':[], 'x':[], 'y':[], 'phi':[]}
    nmode = -99
    bline = -99
    bpts  = -99
    iline = -99
    ipts  = -99
    cnt   = 0

    for eachline in datalines:
        eachline.strip()
        item = eachline.split()

        # setup number of mode
        if len(item)>3 and item[0]=='MAXIMUM' and item[1]=='MODE' and item[2]=='NUMBER:':
            nmode = int(item[-1])
            for i in range( nmode ):
                Bdat[i] = []
                Idat[i] = []

        # search for the mode of field line
        if len(item)>3 and item[0]=='MODE' and item[1]=='OF' and item[2]=='MAGNETIC':
            bline = cnt
            bpts  = int(item[-1])

        # search for the mode of current
        if len(item)>3 and item[0]=='MODE' and item[1]=='OF' and item[2]=='CURRENT':
            iline = cnt
            ipts  = int(item[-1])

        # fill the field data
        if cnt>=bline+2 and cnt<bline+bpts+2:
            Bdat['x'].append( float(item[0]) )
            Bdat['y'].append( float(item[1]) )
            for i in range( nmode ):
                Bdat[i].append( float(item[i+2]) )

        # fill the current data
        if cnt>=iline+2 and cnt<iline+ipts+2:
            Idat['flag'].append( int  (item[0]) )
            Idat['x'   ].append( float(item[1]) )
            Idat['y'   ].append( float(item[2]) )
            Idat['phi' ].append( np.arctan2(float(item[2]), float(item[1])) )
            for i in range( nmode ):
                Idat[i].append( float(item[i+3]) )

        cnt += 1

    for eachkey in Bdat.keys():
        Bdat[eachkey] = np.asarray( Bdat[eachkey] )
    for eachkey in Idat.keys():
        Idat[eachkey] = np.asarray( Idat[eachkey] )

    return Bdat, Idat

########################################
def fit_func(x,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7):
    y = 0.
    a = [a1,a2,a3,a4,a5]
    b = [b1,b2,b3,b4,b5]
    for i in range( len(a) ):
        y += a[i]*np.cos((i+1)*x) + b[i]*np.sin((i+1)*x)
    return y

########################################
def plot_current_mode( data, layer ):
    fig, ax = plt.subplots( 4, 3, figsize=(14,10) )

    cnt = 0
    
    for i in range( len(ax) ):
        for j in range( len(ax[i]) ):
            phi = data['phi'][ data['flag']==layer ]
            iii = data[cnt  ][ data['flag']==layer ]

            par  = (1,1,1,1,1,1,1,1,1,1,1,1,1,1)
            popt, pcov = scipy.optimize.curve_fit( fit_func, phi, iii, p0=par ) 

            pp  = np.linspace( -np.pi, np.pi, 80 )
            yy  = fit_func(pp, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], \
                               popt[7], popt[8], popt[9], popt[10], popt[11], popt[12], popt[13])

            print('MODE',cnt,popt) 

            #ax[i][j].plot( pp, yy, 'r-' )
            #ax[i][j].plot( phi, iii, 's', mec='b', mfc='none' )
            ax[i][j].plot( phi, iii, '-b' )
            ax[i][j].set_title( 'Mode:%4i' %cnt )
            ax[i][j].set_xlabel( '$\phi$ [rad]' )
            ax[i][j].set_ylabel( 'Current [A]' )
            cnt += 1

    plt.tight_layout()
    plt.savefig('current_mode%i.pdf' %layer)
    plt.show()

########################################
def plot_current_distribution( data ):
    fig, ax = plt.subplots( 4, 3, figsize=(14,12) )

    cnt = 0
    
    for i in range( len(ax) ):
        for j in range( len(ax[i]) ):
            if cnt!=len(ax)*len(ax[0])-1:
                mapp = ax[i][j].scatter( data['x'], data['y'], c=data[cnt], marker='o', s=8, edgecolors='none', cmap='jet' )
                fig.colorbar( mapp, ax=ax[i][j] )
                ax[i][j].set_title( 'Mode:%4i' %cnt )
                ax[i][j].set_xlabel( '$x$ [m]' )
                ax[i][j].set_ylabel( '$y$ [m]' )
                ax[i][j].set_aspect('equal')
            cnt += 1

    fin_mode = data[0]
    for i in range(1,15):
        fin_mode += data[i]

    mapp = ax[-1][-1].scatter( data['x'], data['y'], c=fin_mode, marker='o', s=8, edgecolors='none', cmap='jet' )
    fig.colorbar( mapp, ax=ax[-1][-1] )
    ax[-1][-1].set_title( r'Sum mode from 0 upto 10' )
    ax[-1][-1].set_xlabel( '$x$ [m]' )
    ax[-1][-1].set_ylabel( '$y$ [m]' )
    ax[-1][-1].set_aspect('equal')

    plt.tight_layout()
    plt.savefig('current_mode_distribution.pdf')
    plt.show()

########################################
def plot_stream_distribution( data ):
    fig, ax = plt.subplots( 4, 3, figsize=(14,12) )

    cnt = 0
    
    for i in range( len(ax) ):
        for j in range( len(ax[i]) ):
            if cnt!=len(ax)*len(ax[0])-1:
                mapp = ax[i][j].scatter( data['x'], data['y'], c=data[cnt], marker='o', s=8, edgecolors='none', cmap='jet' )
                fig.colorbar( mapp, ax=ax[i][j] )
                ax[i][j].set_title( 'Mode:%4i' %cnt )
                ax[i][j].set_xlabel( '$x$ [m]' )
                ax[i][j].set_ylabel( '$y$ [m]' )
                ax[i][j].set_aspect('equal')
            cnt += 1

    mapp = ax[-1][-1].scatter( data['x'], data['y'], c=data[0]+data[1]+data[2]+data[3]+data[4]+data[5]+data[6], marker='o', s=8, edgecolors='none', cmap='jet' )
    fig.colorbar( mapp, ax=ax[-1][-1] )
    ax[-1][-1].set_title( r'Sum mode from 0 upto 6' )
    ax[-1][-1].set_xlabel( '$x$ [m]' )
    ax[-1][-1].set_ylabel( '$y$ [m]' )
    ax[-1][-1].set_aspect('equal')

    plt.tight_layout()
    plt.savefig('stream_mode_distribution.pdf')
    plt.show()

########################################
def plot_field_distribution( data ):
    fig, ax = plt.subplots( 4, 3, figsize=(14,12) )

    cnt = 0
    
    for i in range( len(ax) ):
        for j in range( len(ax[i]) ):
            if cnt!=len(ax)*len(ax[0])-1:
                mapp = ax[i][j].scatter( data['x'], data['y'], c=data[cnt], marker='s', s=15, edgecolors='none', cmap='jet' )
                fig.colorbar( mapp, ax=ax[i][j] )
                ax[i][j].set_title( 'Mode:%4i' %cnt )
                ax[i][j].set_xlabel( '$x$ [m]' )
                ax[i][j].set_ylabel( '$y$ [m]' )
                ax[i][j].set_aspect('equal')
            cnt += 1

    fin_mode = data[0]
    for i in range(1,15):
        fin_mode += data[i]

    mapp = ax[-1][-1].scatter( data['x'], data['y'], c=(fin_mode-3.5)/3.5, marker='s', s=15, edgecolors='none', cmap='jet' )
    fig.colorbar( mapp, ax=ax[-1][-1] )
    ax[-1][-1].set_title( r'Sum mode from 0 upto 10' )
    ax[-1][-1].set_xlabel( '$x$ [m]' )
    ax[-1][-1].set_ylabel( '$y$ [m]' )
    ax[-1][-1].set_aspect('equal')

    plt.tight_layout()
    plt.savefig('field_mode_distribution.pdf')
    plt.show()

########################################
if __name__=='__main__':
    bdat, idat = load_result_file( sys.argv[1] )

    plot_current_mode( idat, 1 ) 
    plot_current_distribution( idat ) 
    #plot_stream_distribution( phidat ) 
    plot_field_distribution( bdat )


