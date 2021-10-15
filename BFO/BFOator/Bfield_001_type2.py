import numpy as np
mu0 = np.pi*4e-7  # vacuum permeability, T.m/A
muB = 9.274e-24  # Bohr magnetron, J/T or A.m^2
a = 3.96e-10
V = a**3

## Fonctions to compute the magnetic field from 001 BFO, cycloid type 2
# We assume P along 111

# cycloid k1 (k1 along -211)

def prefactor_k1(k, mDM):
    """
    Parameters:
    -----------
    k = cycloid wavevector, 2pi/period
    mDM = max magnetic moment in the SDW
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB.
    """
    return mu0*mDM*muB*np.sinh(k*a*np.sqrt(5/6)/2)/V


def sum_k1(x, y, k, zNV, N):
     """
     Parameters:
     -----------
     k = cycloid wavevector, 2pi/period
     zNV = NV flying height
     N = number of BFO layers
     --------
     Everything in SI units !!! Use lengths in m.
     """
     factor = np.exp(-k*zNV*np.sqrt(5/6))*np.exp(1j*k*(-2*x+y+zNV)/np.sqrt(6))
     fraction = (1-np.exp(-k*a*N*(np.sqrt(5)-1j)/np.sqrt(6)))/\
         (1-np.exp(-k*a*(np.sqrt(5)-1j)/np.sqrt(6)))
     return factor*fraction


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
    A = prefactor_k1(k, mDM)
    S = sum_k1(x, y, k, zNV, N)
    return np.sqrt(2/5)*A*(np.real(S)/np.sqrt(5)+np.imag(S))


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
    k = 2*np.pi/period
    A = prefactor_k1(k, mDM)
    S = sum_k1(x, y, k, zNV, N)
    return A*(np.real(S)/np.sqrt(5)+np.imag(S))/np.sqrt(10)


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
    A = prefactor_k1(k, mDM)
    S = sum_k1(x, y, k, zNV, N)
    return A*((np.sqrt(1/5)+np.sqrt(5))*np.real(S))/np.sqrt(10)


# cycloid k2 (k2 along 1-21)

def prefactor_k2(k, mDM):
    """
    Parameters:
    -----------
    k = cycloid wavevector, 2pi/period
    mDM = max magnetic moment in the SDW
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB.
    """
    return prefactor_k1(k, mDM)


def sum_k2(x, y, k, zNV, N):
     """
     Parameters:
     -----------
     k = cycloid wavevector, 2pi/period
     zNV = NV flying height
     N = number of BFO layers
     --------
     Everything in SI units !!! Use lengths in m.
     """
     factor = np.exp(-k*zNV*np.sqrt(5/6))*np.exp(1j*k*(x-2*y+zNV)/np.sqrt(6))
     fraction = (1-np.exp(-k*a*N*(np.sqrt(5)-1j)/np.sqrt(6)))/\
         (1-np.exp(-k*a*(np.sqrt(5)-1j)/np.sqrt(6)))
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
    S = sum_k2(x, y, k, zNV, N)
    return A*(np.real(S)/np.sqrt(5)+np.imag(S))/np.sqrt(10)


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
    k = 2*np.pi/period
    A = prefactor_k2(k, mDM)
    S = sum_k2(x, y, k, zNV, N)
    return np.sqrt(2/5)*A*(np.real(S)/np.sqrt(5)+np.imag(S))


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
    S = sum_k2(x, y, k, zNV, N)
    return -A*((np.sqrt(1/5)+np.sqrt(5))*np.real(S))/np.sqrt(10)


# cycloid along k3 (k3 along (11-2)

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
    return 0


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
    return 0


# cycloid along k3 with a small rotation alpha in the plane perp to P

def prefactor_k3_alpha(k, mDM, alpha):
    """
    Parameters:
    -----------
    k = cycloid wavevector, 2pi/period
    mDM = max magnetic moment in the SDW
    --------
    Everything in SI units !!! Use lengths in m and mDM in muB.
    """
    sqrt_term = np.sqrt(np.cos(alpha)**2 + 3*np.sin(alpha)**2)
    return 2*ut.mu0*mDM*ut.muB*np.sinh(k*a*sqrt_term/(2*np.sqrt(6)))/(np.sqrt(6)*V)


def sum_k3_alpha(x, y, k, zNV, N, alpha):
    """
    Parameters:
    -----------
    k = cycloid wavevector, 2pi/period
    zNV = NV flying height
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m.
    """
    sqrt_term = np.sqrt(np.cos(alpha)**2 + 3*np.sin(alpha)**2)
    factor = np.exp(-k*zNV*sqrt_term/np.sqrt(6))*\
         np.exp(1j*k*(np.cos(alpha)*(x+y-2*zNV) - np.sqrt(3)*np.sin(alpha)*(x-y))/np.sqrt(6))
    fraction = (1-np.exp(-k*a*N*(sqrt_term+2*1j*np.cos(alpha))/np.sqrt(6)))/\
         (1-np.exp(-k*a*(sqrt_term+2*1j*np.cos(alpha))/np.sqrt(6)))
    return factor*fraction


def Bx_k3_alpha(x, y, zNV, period, mDM, N, alpha):
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
    sqrt_term = np.sqrt(np.cos(alpha)**2 + 3*np.sin(alpha)**2)
    A = prefactor_k3_alpha(k, mDM, alpha)
    S = sum_k3_alpha(x, y, k, zNV, N, alpha)
    return A*(np.sin(alpha)*(np.cos(alpha)-np.sqrt(3)*np.sin(alpha)))*np.imag(S)/sqrt_term


def By_k3_alpha(x, y, zNV, period, mDM, N, alpha):
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
    sqrt_term = np.sqrt(np.cos(alpha)**2 + 3*np.sin(alpha)**2)
    A = prefactor_k3_alpha(k, mDM, alpha)
    S = sum_k3_alpha(x, y, k, zNV, N, alpha)
    return A*(np.sin(alpha)*(np.cos(alpha)+np.sqrt(3)*np.sin(alpha)))*np.imag(S)/sqrt_term


def Bz_k3_alpha(x, y, zNV, period, mDM, N, alpha):
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
    sqrt_term = np.sqrt(np.cos(alpha)**2 + 3*np.sin(alpha)**2)
    A = prefactor_k3_alpha(k, mDM, alpha)
    S = sum_k3_alpha(x, y, k, zNV, N, alpha)
    return -A*np.sin(alpha)*(-np.real(S)+2*np.cos(alpha)*np.imag(S)/sqrt_term)
