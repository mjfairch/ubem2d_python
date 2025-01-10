import numpy as np
from ubem2d.aerodynamics.Airfoil import Airfoil
from ubem2d.geometry.Orientation import Orientation

__all__ = ['naca4']

def parse_naca4_code(code):
    if (len(code) != 4):
        raise ValueError('NACA code must be four digits')
    m = .01*float(code[0])    # max camber
    p = .1*float(code[1])     # position of max camber along chord
    th = .01*float(code[2:])  # thickness
    if (m == 0 and p != 0):
        raise ValueError('m==0 and p != 0 are inconsistent')
    if (m != 0 and (p == 0 or p == 1)):
        raise ValueError('p must lie in the interval (0,1) when m != 0')
    return m, p, th

def naca4(code, npan, clamp = True, uniform = False):
    '''
    This method returns an Airfoil object which represents a NACA airfoil of
    the given type and number of panels.  For more details, refer to the
    literature, or to: https://en.wikipedia.org/wiki/NACA_airfoil
    '''
    (m,p,th) = parse_naca4_code(code)
    if (npan != int(npan)):
        raise ValueError('Number of panels must be an integer')
    c = 1.
    # Check input
    if (npan < 4):
        raise Valueerror('Must specify at least four panels')
    # Polynomial coefficients for NACA 4-digit series
    coefs = np.array([-.1015, .2843, -.3516, -.1260, 0.])
    if (clamp):
        coefs[0] = -.1036
    sqrtCoef = .2969
    # Corner counts
    n = npan + 1         # number of corners along entire airfoil
    nu = int(n/2) + n%2  # number of corners along upper surface
    nl = n-nu            # number of corners along lower surface
    # Abcissa along upper surface
    if (uniform):
        xu = np.linspace(c, 0, nu)
    else:
        t = np.linspace(.5*np.pi, 0, nu)
        xu = c*(1-np.cos(t))
    # Abcissa along lower surface determined by those along upper surface
    xl = np.zeros(nl)
    if (npan % 2 == 0):
        for i in range(nl):
            xl[i] = xu[-(i+2)]
    else:
        dxu = np.diff(xu)
        for i in range(nl-1):
            xl[i] = xu[-(i+1)]-.5*dxu[-(i+1)]
        xl[-1] = xu[0]
    # Construct thickness distribution along upper and lower surfaces
    ytu = 5*th*c*(sqrtCoef*np.sqrt(xu/c) + np.polyval(coefs, xu/c))
    ytl = 5*th*c*(sqrtCoef*np.sqrt(xl/c) + np.polyval(coefs, xl/c))
    # Compute upper and lower surfaces
    if (m == 0):
        x = np.concatenate((xu, xl))
        y = np.concatenate((ytu, -1.*ytl))
    else:
        ycu = np.zeros(len(xu))
        dycudx = np.zeros(len(xu))
        ycl = np.zeros(len(xl))
        dycldx = np.zeros(len(xl))
        for i in range(len(ycu)):
            if (xu[i] <= p*c):
                ycu[i] = m*xu[i]/(p*p)*(2*p-xu[i]/c)
                dycudx[i] = 2*m/(p*p)*(p-xu[i]/c)
            else:
                ycu[i] = m*(c-xu[i])/((1-p)**2)*(1+xu[i]/c-2*p)
                dycudx[i] = 2*m/((1-p)**2)*(p-xu[i]/c)
        for i in range(len(ycl)):
            if (xl[i] <= p*c):
                ycl[i] = m*xl[i]/(p*p)*(2*p-xl[i]/c)
                dycldx[i] = 2*m/(p*p)*(p-xl[i]/c)
            else:
                ycl[i] = m*(c-xl[i])/((1-p)**2)*(1+xl[i]/c-2*p)
                dycldx[i] = 2*m/((1-p)**2)*(p-xl[i]/c)
        thu = np.arctan(dycudx)
        thl = np.arctan(dycldx)
        xU = xu - ytu*np.sin(thu)
        yU = ycu + ytu*np.cos(thu)
        xL = xl + ytl*np.sin(thl)
        yL = ycl - ytl*np.cos(thl)
        x = np.concatenate((xU, xL))
        y = np.concatenate((yU, yL))

    return Airfoil(x, y, nu-1, Orientation.CW)

if __name__  == '__main__':
    import sys
    import matplotlib.pyplot as plt
    if (len(sys.argv) < 3):
        sys.exit('Usage: {} naca_code num_panels'.format(sys.argv[0]))
    code = sys.argv[1]
    try:
        npan = int(sys.argv[2])
    except ValueError:
        sys.exit('Bad code: {}'.format(sys.argv[2]))

    foil = naca4(code,npan)
    for i in range(foil.nedge+1):
        print('{:<20.16g}\t{:<20.16g}'.format(foil.x[i],foil.y[i]))
    print('\nClosed:',foil.closed)
    print('Perimeter:',foil.perimeter)
    print('Chord:',foil.chord)
    print('Diameter:',foil.diameter)
    print('Centroid:',foil.centroid)

    plt.plot(foil.x/foil.chord,foil.y/foil.chord,'.-')
    plt.xlabel('x/chord')
    plt.ylabel('y/chord')
    plt.title('NACA {} with {} panels'.format(code,npan))
    plt.axis('equal')
    plt.grid(True)
    plt.show()
