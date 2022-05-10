#!/usr/bin/env python

import numpy as np
import sys

rho0  = 1.89
phi0  = -np.pi/4
phi1  =  np.pi/4
nphi  = 100
rref  = 1.89
theta0= 0.
theta1= 2.*np.pi
ntheta= 30
r0    = 0.001
r1    = 0.030
nr    = 20

def generate_points():
    phi   = np.linspace(   phi0,   phi1,   nphi )
    theta = np.linspace( theta0, theta1, ntheta )
    r     = np.linspace(     r0,     r1,     nr )

    data = {'x':[],'y':[],'z':[]}

    for i in range( nphi ):
        for j in range( nr ):
            for k in range( ntheta ):
                x = (rho0 + r[j]*np.cos(theta[k])) * np.cos(phi[i])
                y = r[j]*np.sin(theta[k])
                z = (rho0 + r[j]*np.cos(theta[k])) * np.sin(phi[i])

                data['x'].append( x )
                data['y'].append( y )
                data['z'].append( z )

    for eachkey in data.keys():
        data[eachkey] = np.asarray( data[eachkey] )

    return data

def write_opera_table( data, filename ):
    outfile = open( filename, 'w' )

    outfile.write( ' %i 1 1 2\n' %(len(data['x'])) )
    outfile.write( ' 1 X [METRE]\n' )
    outfile.write( ' 2 Y [METRE]\n' )
    outfile.write( ' 3 Z [METRE]\n' )
    outfile.write( ' 0\n' )

    for i in range( len(data['x']) ):
        outfile.write( '%20.9e%20.9e%20.9e\n' %(data['x'][i], data['y'][i], data['z'][i]) )

    outfile.close()


if __name__=='__main__':
    data = generate_points()
    write_opera_table( data, 'bore_points.table' )



