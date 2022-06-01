import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'   ] = 'Times New Roman'
plt.rcParams['axes.linewidth'] = 0.5

def get_elliptical_coord(a, b, phi0, phi1, nphi):
    phi = np.linspace( phi0, phi1, nphi )
    eta = np.arctanh( b/a )
    e   = np.sqrt(a**2 - b**2)
    x   = e * np.cosh(eta) * np.cos(phi)
    y   = e * np.sinh(eta) * np.sin(phi)
    dis = np.sqrt( (x[1]-x[0])**2 + (y[1]-y[0])**2 )
    #print( 'DISTANCE:%.2f[mm]' %(dis*1e+3) )
    return x[1:-1], y[1:-1]

def generate_multilayer(a0, b0, nlayer, d, phi0=0., phi1=np.pi*0.5, nphi=102):
    layer = {'n':np.array([]), 'x':np.array([]), 'y':np.array([])}
    
    for i in range(nlayer):
        a = a0 + i*d
        b = b0 + i*d
        x, y = get_elliptical_coord( a, b, phi0, phi1, nphi )
        layer['n'] = np.append( layer['n'], i+1 )
        layer['x'] = np.append( layer['x'], x )
        layer['y'] = np.append( layer['y'], y )

    layer['x'] = np.reshape( layer['x'], (nphi-2,nlayer) )
    layer['y'] = np.reshape( layer['y'], (nphi-2,nlayer) )

    return layer

def generate_input_file(filename, x, y, curr=1.0):
    outfile = open( filename, 'w' )

    if x.ndim==2:
        for i in range( len(x) ):
            for j in range( len(x[i]) ):
                outfile.write('%18.9e%18.9e%10.1f\n' %(x[i][j],y[i][j],curr))

    if x.ndim==1:
        for i in range( len(x) ):
            outfile.write('%18.9e%18.9e%10.1f\n' %(x[i],y[i],curr))

    outfile.close()

def plot_points(x, y):
    fig, ax = plt.subplots( 1, 1, figsize=(6,5) )
    
    if x.ndim==1:
        ax.plot( x, y, 'xy' )

    if x.ndim==2:
        color = plt.cm.viridis(np.linspace(0,1,len(x)))
        for i in range( len(x) ):
            ax.plot( x[i], y[i], '.', c=color[i] )

    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.show()
