#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import sys

####################################################
plt.style.use('classic')
plt.rcParams['font.family'   ] = 'Times New Roman'
plt.rcParams['axes.linewidth'] = 0.8
####################################################

A0    = 0.075
B0    = 0.055
PHI0  = -np.pi
PHI1  =  np.pi
NPHI  = 181
PITCH = 0.002
NLAYER= 16

####################################################
def load_result_file( filename ):
    infile    = open( filename, 'r' )
    datalines = infile.readlines()
    infile.close()

    Bdat  = {'x':[], 'y':[]}
    Idat  = {'flag':[], 'x':[], 'y':[]}
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
            for i in range( nmode ):
                Idat[i].append( float(item[i+3]) )

        cnt += 1

    for eachkey in Bdat.keys():
        Bdat[eachkey] = np.asarray( Bdat[eachkey] )
    for eachkey in Idat.keys():
        Idat[eachkey] = np.asarray( Idat[eachkey] )

    return Bdat, Idat

####################################################
def calc_streamfunc( data, layer, mode0, mode1 ):
    x  = data['x'  ][ data['flag']==layer ]
    y  = data['y'  ][ data['flag']==layer ]
    Is = data[mode0][ data['flag']==layer ]
    phi= np.array([])
    sf = np.array([])
    
    Phi_s = 0.

    if mode1>mode0+1:
        for i in range( mode0+1, mode1 ):
            Is += data[i][ data['flag']==layer ]

    for i in range( len(x) ):
        phi = np.append( phi, PHI0 + i*(PHI1-PHI0)/(NPHI-1) )

        if i==0:
            Phi_s += (Is[-1] + Is[i]) * 0.5
        else:
            Phi_s += (Is[i-1] + Is[i]) * 0.5
        '''
        if i==0:
            Phi_s += (Is[-1]/3.+Is[0]*4/3.+Is[1]/3.)
        elif i==len(x)-1:
            Phi_s += (Is[0]/3.+Is[-1]*4/3.+Is[-2]/3.)
        else:
            Phi_s += (Is[i-1]/3. + Is[i]*4/3. + Is[i+1]/3.)
        '''

        sf = np.append( sf, Phi_s )

    return x, y, phi, Is, sf

####################################################
def search_for_contour( icontour, phi, sf ):
    phi_new = np.array([])
    con_new = np.array([])

    for i in range( 1,len(phi) ):
        if (icontour>sf[i-1] and icontour<=sf[i]) or (icontour>sf[i] and icontour<=sf[i-1]):
            this_phi = phi[i-1] + (phi[i]-phi[i-1])/(sf[i]-sf[i-1]) * (icontour-sf[i-1])
            phi_new = np.append( phi_new, this_phi )
            con_new = np.append( con_new, icontour )

    return phi_new, con_new

####################################################
def calc_contour( layer, phi, sf, sfmin, sfmax, dsf ):
    phi_new = np.array([])
    con_new = np.array([])

    sline = sfmin

    if sfmax>np.max(sf):
        sfmax = np.max(sf)
    if sfmin<np.min(sf):
        sfmin = np.min(sf)

    while sline < sfmax:
        tphi, tcon = search_for_contour(sline, phi, sf)
        phi_new = np.append( phi_new, tphi )
        con_new = np.append( con_new, tcon )
        sline += dsf

    print( 'NUMBER OF CONTOUR FOUND:', len(phi_new) )

    a = A0 + (layer-1)*PITCH
    b = B0 + (layer-1)*PITCH
    e = np.sqrt( a*a - b*b )
    eta0 = np.arctanh( b/a )

    x = e * np.cosh(eta0) * np.cos(phi_new)
    y = e * np.sinh(eta0) * np.sin(phi_new)

    pos = {'x':x, 'y':y, 'phi':phi_new, 'sf':con_new}

    return pos

####################################################
def plot_streamfunc_map( data ):
    fig, ax = plt.subplots( 1, 2, figsize=(12,5) )
    
    MODE0 = 0
    MODE1 = 6
    
    for i in range( NLAYER ):
        x, y, phi, Is, sf = calc_streamfunc( data, i+1, MODE0, MODE1 )
        ax[0].scatter( x, y, c=Is, s=20, marker='o', cmap='jet', edgecolors='none' )
        ax[1].scatter( x, y, c=sf, s=20, marker='o', cmap='jet', edgecolors='none' )

        for j in range( len(ax) ):
            ax[j].set_xlabel(r'$x$ [m]')
            ax[j].set_ylabel(r'$y$ [m]')
            ax[j].set_aspect('equal')

    plt.tight_layout()
    plt.show()

####################################################
def generate_wire_position(data, filename, mode0, mode1, d_sf):
    x, y, phi, Is, sf = calc_streamfunc( data, 1, mode0, mode1 )
    sf_max = np.max(sf)
    sf_min = np.min(sf)
    turns  = (sf_max - sf_min) / d_sf

    print( 'Maximum stream function:%8.2f A' %sf_max )
    print( 'Minimum stream function:%8.2f A' %sf_min )
    print( 'Number of contour:%8i' %turns )
    print( 'Interval of stream function:%8.2f A' %d_sf )

    fig, ax = plt.subplots( 1, 2, figsize=(12,4) )

    infile = open( filename, 'w' )

    for i in range( NLAYER ):
        x,y,phi,Is,sf = calc_streamfunc( data, i+1, mode0, mode1 )

        if i==0:
            sf_min = np.min(sf)

        sf_min = -800000
        sf_max =  800000

        pp = phi[ (phi>=-np.pi*1.0) & (phi<=np.pi*1.0) ]
        ss = sf [ (phi>=-np.pi*1.0) & (phi<=np.pi*1.0) ]

        pos = calc_contour( i+1, pp, ss, sf_min, sf_max, d_sf )
        pos['I'] = -np.ones(len(pos['x']))

        ax[0].plot( phi, sf, 'k-' )
        ax[0].plot( pos['phi'], pos['sf'], 'xr' )

        ax[1].scatter( pos['x'], pos['y'], c= pos['I'], s=6, cmap='bwr', edgecolors='none', vmin=-1, vmax=1)
        #ax[1].scatter( pos['x'],-pos['y'], c= pos['I'], s=6, cmap='bwr', edgecolors='none', vmin=-1, vmax=1)
        #ax[1].scatter(-pos['x'], pos['y'], c=-pos['I'], s=6, cmap='bwr', edgecolors='none', vmin=-1, vmax=1)
        #ax[1].scatter(-pos['x'],-pos['y'], c=-pos['I'], s=6, cmap='bwr', edgecolors='none', vmin=-1, vmax=1)

        for j in range( len(pos['x']) ):
            infile.write( '%4i%15.6e%15.6e%15.6e%5i\n' %(i+1,pos['phi'][j], pos['x'][j], pos['y'][j],-1) )
            # flip
            #infile.write( '%4i%15.6e%15.6e%15.6e%5i\n' %(i+1,pos['phi'][j], pos['x'][j],-pos['y'][j],-1) )
            #infile.write( '%4i%15.6e%15.6e%15.6e%5i\n' %(i+1,pos['phi'][j],-pos['x'][j], pos['y'][j], 1) )
            #infile.write( '%4i%15.6e%15.6e%15.6e%5i\n' %(i+1,pos['phi'][j],-pos['x'][j],-pos['y'][j], 1) )

    infile.close()

    ax[1].set_aspect('equal')
    ax[1].set_xlabel('$x$ [m]')
    ax[1].set_ylabel('$y$ [m]')

    plt.tight_layout()
    plt.savefig('coil_pattern.pdf')
    plt.show()

####################################################
def plot_streamfunc( data ):
    fig, ax = plt.subplots( 1, 2, figsize=(12,5) )
    
    MODE0 = 3
    MODE1 = 3
    
    print( '%6s%6s%14s%14s' %('MODE','LAYER','SIS','SSF') )
    for i in range( 11, 13 ):
        x, y, phi, Is, sf = calc_streamfunc( data, i+1, MODE0, MODE1 )

        print( '%6i%6i%14.6e%14.6e%14.6e%14.6e' %(MODE,i+1,np.sum(Is), np.sum(sf),np.min(Is),np.max(Is)) )

        ax[0].plot( phi*180/np.pi, Is )
        ax[1].plot( phi*180/np.pi, sf )

        for j in range( len(ax) ):
            ax[j].set_xlabel(r'$\phi$ [deg]')

    plt.tight_layout()
    plt.show()

####################################################
if __name__=='__main__':
    bdat, idat = load_result_file( sys.argv[1] ) 
    #plot_streamfunc_map( idat )
    #plot_streamfunc( idat )
    generate_wire_position(idat, 'wire_position.txt', 0, 4, 500)
