import numpy as np
mu0 = np.pi*4e-7  # vacuum permeability, T.m/A
muB = 9.274e-24  # Bohr magnetron, J/T or A.m^2a
a = 3.96e-10
#a = 3.99e-10
V = a**3
#V = a*3.95e-10*3.95e-10

## Fonctions to compute the magnetic field from 001 BFO, cycloid type 1
# We assume P along 111

# cycloid k1 (k1 along 1-10)

def prefactor_k1(zNV, k, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    k = cycloid wavevector, 2pi/period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB.
    """
    return mu0*mDM*muB*np.exp(-k*zNV)*np.sinh(k*a/2)*\
        (1-np.exp(-k*a*N))/(np.sqrt(3)*V*(1-np.exp(-k*a)))


def Bx_k1(x, y, zNV, period, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    return prefactor_k1(zNV, k, mDM, N)*np.sin(k*(x-y)/np.sqrt(2))


def By_k1(x, y, zNV, period, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    return -Bx_k1(x, y, zNV, period, mDM, N)


def Bz_k1(x, y, zNV, period, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    return np.sqrt(2)*prefactor_k1(zNV, k, mDM, N)*np.cos(k*(x-y)/np.sqrt(2))



# cycloid k2 (k2 along -101)

def prefactor_k2(k, mDM):
    """
    Parameters:
    -----------
    k = cycloid wavevector, 2pi/period
    mDM = max magnetic moment in the SDW
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB.
    """
    return mu0*mDM*muB*np.sinh(k*a/(2*np.sqrt(2)))/(np.sqrt(3)*V)

def sum_k2(x, k, zNV, N):
     """
     Parameters:
     -----------
     k = cycloid wavevector, 2pi/period
     zNV = NV flying height
     N = number of BFO layers
     --------
     Everything in SI units !!! Use lengths in m.
     """
     factor = np.exp(-k*zNV/np.sqrt(2))*np.exp(1j*k*(x-zNV)/np.sqrt(2))
     fraction = (1-np.exp(-k*a*N*(1+1j)/np.sqrt(2)))/(1-np.exp(-k*a*(1+1j)/np.sqrt(2)))
     return factor*fraction


def Bx_k2(x, y, zNV, period, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    A = prefactor_k2(k, mDM)
    S = sum_k2(x, k, zNV, N)
    return -A*(np.real(S)-np.imag(S))/np.sqrt(2)


def By_k2(x, y, zNV, period, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    return 0


def Bz_k2(x, y, zNV, period, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    A = prefactor_k2(k, mDM)
    S = sum_k2(x, k, zNV, N)
    return np.sqrt(2)*A*np.real(S)


# cycloid k3 (k3 along 0-11)

def prefactor_k3(k, mDM):
    """
    Parameters:
    -----------
    k = cycloid wavevector, 2pi/period
    mDM = max magnetic moment in the SDW
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB.
    """
    return prefactor_k2(k, mDM)


def sum_k3(y, k, zNV, N):
     """
     Parameters:
     -----------
     k = cycloid wavevector, 2pi/period
     zNV = NV flying height
     N = number of BFO layers
     --------
     Everything in SI units !!! Use lengths in m.
     """
     return sum_k2(y, k, zNV, N)


def Bx_k3(x, y, zNV, period, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    return 0

 
def By_k3(x, y, zNV, period, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    A = prefactor_k3(k, mDM)
    S = sum_k3(y, k, zNV, N)
    return -A*(np.real(S)-np.imag(S))/np.sqrt(2)


def Bz_k3(x, y, zNV, period, mDM, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    A = prefactor_k3(k, mDM)
    S = sum_k3(y, k, zNV, N)
    return np.sqrt(2)*A*np.real(S)
