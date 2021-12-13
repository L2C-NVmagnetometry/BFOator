import numpy as np
mu0 = np.pi*4e-7  # vacuum permeability, T.m/A
muB = 9.274e-24  # Bohr magnetron, J/T or A.m^2
a = 3.96e-10
V = a**3

# We assume k along x, and cycloid type 1

def prefactor(zNV, k, ms, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    k = cycloid wavevector, 2pi/period
    ms = magnetic moment of the Fe atoms
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB.
    """
    fraction = (1-(-np.exp(-k*a/np.sqrt(3)))**N)/(1+np.exp(-k*a/np.sqrt(3)))                
    return 2*mu0*ms*muB*np.exp(-k*zNV)*np.sinh(k*a*np.sqrt(3)/6)*fraction/V
                                              
                                              
def Bx(x, y, zNV, period, ms, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    return prefactor(zNV, k, ms, N)*np.cos(k*x)


def By(x, y, zNV, period, ms, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    return np.zeros(len(x))


def Bz(x, y, zNV, period, ms, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    return -prefactor(zNV, k, ms, N)*np.sin(k*x)                

                                              
