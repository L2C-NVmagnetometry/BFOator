import numpy as np
mu0 = np.pi*4e-7  # vacuum permeability, T.m/A
muB = 9.274e-24  # Bohr magnetron, J/T or A.m^2
a = 3.96e-10
V = a**3

# We assume k1 along x, and cycloid type 1, P is along [11-1]

def prefactor_k1(zNV, k):
    """
    Parameters:
    -----------
    zNV = NV flying height
    k = cycloid wavevector, 2pi/period
    --------
    Everything in SI units !!! Use lengths in m.
    """ 
    return mu0*np.exp(-k*zNV)*np.sinh(k*a*np.sqrt(3)/6)/V

def sum_term_cycloid_k1(k, N):
    """
    Parameters:
    -----------
    k = cycloid wavevector, 2pi/period
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m.
    """
    return (1-(-np.exp(-k*a/np.sqrt(3)))**N)/(1+np.exp(-k*a/np.sqrt(3)))                

def sum_term_sdw_k1(k, N):
    """
    Parameters:
    -----------
    k = cycloid wavevector, 2pi/period
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m.
    """
    return (1-np.exp(-k*a*N/np.sqrt(3)))/(1-np.exp(-k*a/np.sqrt(3)))                
                                              
def Bx_k1(x, y, zNV, period, ms, mDM, N, eps, phase):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    eps = +/- 1, rotational sense
    phase = in deg, dephasing between the sdw and the cycloid
    --------
    Everything in SI units !!! Use lengths in m and ms and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    s_cycl = sum_term_cycloid_k1(k, N)
    s_sdw = sum_term_sdw_k1(k, N)
    pref = prefactor_k1(zNV, k)
    phi = phase*np.pi/180
    return pref*(ms*muB*(1-eps/3)*s_cycl*np.cos(k*x)+2*np.sqrt(2)*mDM*muB*s_sdw*np.sin(k*x+phi)/3)


def By_k1(x, y, zNV, period, ms, mDM, N, eps, phase):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    eps = +/- 1, rotational sense
    phase = in deg, dephasing between the sdw and the cycloid
    --------
    Everything in SI units !!! Use lengths in m and ms and mDM in muB, the output is in T.
    """
    return np.zeros(len(x))


def Bz_k1(x, y, zNV, period, ms, mDM, N, eps, phase):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    eps = +/- 1, rotational sense
    phase = in deg, dephasing between the sdw and the cycloid
    --------
    Everything in SI units !!! Use lengths in m and ms and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    s_cycl = sum_term_cycloid_k1(k, N)
    s_sdw = sum_term_sdw_k1(k, N)
    pref = prefactor_k1(zNV, k)
    phi = phase*np.pi/180
    return -pref*(ms*muB*(1-eps/3)*s_cycl*np.sin(k*x)-2*np.sqrt(2)*mDM*muB*s_sdw*np.cos(k*x+phi)/3)        

                                              
def prefactor_k23(zNV, k):
    """
    Parameters:
    -----------
    zNV = NV flying height
    k = cycloid wavevector, 2pi/period
    --------
    Everything in SI units !!! Use lengths in m.
    """ 
    return mu0*np.exp(-k*zNV/np.sqrt(3))*np.sinh(k*a/6)/(3*V)


def sum_term_cycloid_k2(x, y, zNV, k, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    k = cycloid wavevector, 2pi/period
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m.
    """
    kr = -k*x/2 - k*y/(2*np.sqrt(3)) + np.sqrt(2)*k*zNV/np.sqrt(3)
    expo = np.exp(1j*kr)
    num = 1-(-np.exp(-k*a*(1-1j*np.sqrt(2))/3))**N
    den = 1+np.exp(-k*a*(1-1j*np.sqrt(2))/3) 
    return expo*num/den


def sum_term_sdw_k2(x, y, zNV, k, N, phi):
    """
    Parameters:
    -----------
    zNV = NV flying height
    k = cycloid wavevector, 2pi/period
    N = number of BFO layers
    phi = dephasing between the cycloid and the SDW, in rad.
    --------
    Everything in SI units !!! Use lengths in m.
    """
    kr = -k*x/2 - k*y/(2*np.sqrt(3)) + np.sqrt(2)*k*zNV/np.sqrt(3)
    expo = np.exp(1j*(kr+phi))
    num = 1-np.exp(-k*a*N*(1-1j*np.sqrt(2))/3)
    den = 1-np.exp(-k*a*(1-1j*np.sqrt(2))/3) 
    return expo*num/den


def Bx_k2(x, y, zNV, period, ms, mDM, N, eps, phase):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    eps = +/- 1, rotational sense
    phase = in deg, dephasing between the sdw and the cycloid
    --------
    Everything in SI units !!! Use lengths in m and ms and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    pref = np.sqrt(3)*prefactor_k23(zNV, k)/2
    phi = phase*np.pi/180
    s_cycl = sum_term_cycloid_k2(x, y, zNV, k, N)
    s_sdw = sum_term_sdw_k2(x, y, zNV, k, N, phi)
    cycl = ms*muB*(np.sqrt(3)+eps)*(np.real(s_cycl) - np.sqrt(2)*np.imag(s_cycl))
    sdw = mDM*muB*np.sqrt(2)*( np.sqrt(2)*np.real(s_sdw) +np.imag(s_sdw))
    return pref*(cycl+sdw)


def By_k2(x, y, zNV, period, ms, mDM, N, eps, phase):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    eps = +/- 1, rotational sense
    phase = in deg, dephasing between the sdw and the cycloid
    --------
    Everything in SI units !!! Use lengths in m and ms and mDM in muB, the output is in T.
    """
    return Bx_k2(x, y, zNV, period, ms, mDM, N, eps, phase)/np.sqrt(3)


def Bz_k2(x, y, zNV, period, ms, mDM, N, eps, phase):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    eps = +/- 1, rotational sense
    phase = in deg, dephasing between the sdw and the cycloid
    --------
    Everything in SI units !!! Use lengths in m and ms and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    pref = prefactor_k23(zNV, k)
    phi = phase*np.pi/180
    s_cycl = sum_term_cycloid_k2(x, y, zNV, k, N)
    s_sdw = sum_term_sdw_k2(x, y, zNV, k, N, phi)
    cycl = ms*muB*3*(eps+np.sqrt(3))*np.imag(s_cycl)
    sdw  = -np.sqrt(2)*mDM*muB*(2*np.sqrt(2)*np.imag(s_sdw)+ np.real(s_sdw))
    return pref*(cycl + sdw)



def sum_term_cycloid_k3(x, y, zNV, k, N):
    """
    Parameters:
    -----------
    zNV = NV flying height
    k = cycloid wavevector, 2pi/period
    N = number of BFO layers
    --------
    Everything in SI units !!! Use lengths in m.
    """
    kr = -k*x/2 + k*y/(2*np.sqrt(3)) - np.sqrt(2)*k*zNV/np.sqrt(3)
    expo = np.exp(1j*kr)
    num = 1-(-np.exp(-k*a*(1+1j*np.sqrt(2))/3))**N
    den = 1+np.exp(-k*a*(1+1j*np.sqrt(2))/3) 
    return expo*num/den


def sum_term_sdw_k3(x, y, zNV, k, N, phi):
    """
    Parameters:
    -----------
    zNV = NV flying height
    k = cycloid wavevector, 2pi/period
    N = number of BFO layers
    phi = dephasing between the cycloid and the SDW, in rad.
    --------
    Everything in SI units !!! Use lengths in m.
    """
    kr = -k*x/2 + k*y/(2*np.sqrt(3)) - np.sqrt(2)*k*zNV/np.sqrt(3)
    expo = np.exp(1j*(kr+phi))
    num = 1-np.exp(-k*a*N*(1+1j*np.sqrt(2))/3)
    den = 1-np.exp(-k*a*(1+1j*np.sqrt(2))/3) 
    return expo*num/den


def Bx_k3(x, y, zNV, period, ms, mDM, N, eps, phase):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    eps = +/- 1, rotational sense
    phase = in deg, dephasing between the sdw and the cycloid
    --------
    Everything in SI units !!! Use lengths in m and ms and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    pref = np.sqrt(3)*prefactor_k23(zNV, k)/2
    phi = phase*np.pi/180
    s_cycl = sum_term_cycloid_k3(x, y, zNV, k, N)
    s_sdw = sum_term_sdw_k3(x, y, zNV, k, N, phi)
    cycl = ms*muB*(np.sqrt(3)+eps)*(np.real(s_cycl) + np.sqrt(2)*np.imag(s_cycl))
    sdw = -np.sqrt(2)*mDM*muB*(np.sqrt(2)*np.real(s_sdw) - np.imag(s_sdw))
    return pref*(cycl+sdw)


def By_k3(x, y, zNV, period, ms, mDM, N, eps, phase):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    eps = +/- 1, rotational sense
    phase = in deg, dephasing between the sdw and the cycloid
    --------
    Everything in SI units !!! Use lengths in m and ms and mDM in muB, the output is in T.
    """
    return -Bx_k3(x, y, zNV, period, ms, mDM, N, eps, phase)/np.sqrt(3)


def Bz_k3(x, y, zNV, period, ms, mDM, N, eps, phase):
    """
    Parameters:
    -----------
    zNV = NV flying height
    period = cycloid period
    ms = magnetic moment of the Fe atoms
    mDM = max magnetic moment in the SDW
    N = number of BFO layers
    eps = +/- 1, rotational sense
    phase = in deg, dephasing between the sdw and the cycloid
    --------
    Everything in SI units !!! Use lengths in m and ms and mDM in muB, the output is in T.
    """
    k = 2*np.pi/period
    pref = prefactor_k23(zNV, k)
    phi = phase*np.pi/180
    s_cycl = sum_term_cycloid_k3(x, y, zNV, k, N)
    s_sdw = sum_term_sdw_k3(x, y, zNV, k, N, phi)
    cycl = 3*ms*muB*(np.sqrt(3)+eps)*np.imag(s_cycl)
    sdw  = - np.sqrt(2)*mDM*muB*(2*np.sqrt(2)*np.imag(s_sdw) - np.real(s_sdw))
    return pref*(cycl + sdw)
