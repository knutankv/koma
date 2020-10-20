"""
##############################################
Modal post-processing module
##############################################
All functions related to the processing of results from OMA methods or in general.
"""


import numpy as np

def xmacmat(phi1, phi2=None, conjugates=True):
    """
    Modal assurance criterion numbers, cross-matrix between two modal transformation matrices (modes stacked as columns).

    Arguments
    ---------------------------
    phi1 : double
        reference modes
    phi2 : double, optional
        modes to compare with, if not given (i.e., equal default value None), phi1 vs phi1 is assumed
    conjugates : True, optional
        check the complex conjugates of all modes as well (should normally be True)

    Returns
    ---------------------------
    macs : double
        matrix of MAC numbers
    """
    # If no phi2 is given, assign value of phi1
    if phi2 is None:
        phi2 = 1.0*phi1
        
    if len(np.shape(phi1))==1:
        phi1 = np.expand_dims(phi1, axis=0).T
        
    if len(np.shape(phi2))==1:
        phi2 = np.expand_dims(phi2, axis=0).T

    # norms1 = np.dot(np.expand_dims(np.sum(phi1.T * phi1.T.conj(), axis=1), axis=0), np.expand_dims(np.sum(phi2.T * phi2.T.conj(),axis=1), axis=1))
    norms = np.real(np.sum(phi1.T * np.conj(phi1.T), axis=1))[:,np.newaxis] @ np.real(np.sum(phi2.T * np.conj(phi2.T),axis=1))[np.newaxis,:]


    if conjugates:
        macs1 = np.divide(abs(np.dot(phi1.T, phi2))**2, norms)
        macs2 = np.divide(abs(np.dot(phi1.T, phi2.conj()))**2, norms)     
        macs = np.maximum(macs1, macs2)
    else:
        macs = np.divide(abs(np.dot(phi1.T, phi2))**2, norms)

    macs = np.real(macs)
    
    if np.size(macs) == 1:
        macs = macs[0,0]
    
    return macs
    

def xmacmat_alt(phi1, phi2=None, conjugates=True):
    """
    Alternative implementation. Modal assurance criterion numbers, cross-matrix between two modal transformation matrices (modes stacked as columns).

    Arguments
    ---------------------------
    phi1 : double
        reference modes
    phi2 : double, optional
        modes to compare with, if not given (i.e., equal default value None), phi1 vs phi1 is assumed
    conjugates : True, optional
        check the complex conjugates of all modes as well (should normally be True)

    Returns
    ---------------------------
    macs : boolean
        matrix of MAC numbers
    """

    if phi2 is None:
        phi2 = phi1
    
    if phi1.ndim==1:
        phi1 = phi1[:, np.newaxis]
    
    if phi2.ndim==1:
        phi2 = phi2[:, np.newaxis]
    
    A = np.sum(phi1.T * np.conj(phi1.T), axis=1)[:,np.newaxis]
    B = np.sum(phi2.T * np.conj(phi2.T), axis=1)[:, np.newaxis]
    norms = np.abs(A @ B.T)

    if conjugates:
        macs = np.maximum(abs(phi1.T @ phi2)**2/norms, abs(phi1.T @ phi2)**2/norms)     #maximum = element-wise max
    else:
        macs = abs(phi1.T @ phi2)**2/norms
    
    macs = np.real(macs)
    
    return macs


def mac(phi1, phi2):
    """
    Modal assurance criterion number from comparison of two mode shapes.

    Arguments
    ---------------------------
    phi1 : double
        first mode
    phi2 : double
        second mode

    Returns
    ---------------------------
    mac_value : boolean
        MAC number
    """

    mac_value = np.real(np.abs(np.dot(phi1.T,phi2))**2 / np.abs((np.dot(phi1.T, phi1) * np.dot(phi2.T, phi2))))
    return mac_value


def maxreal(phi):
    """
    Rotate complex vectors (stacked column-wise) such that the absolute values of the real parts are maximized.

    Arguments
    ---------------------------
    phi : double
        complex-valued modal transformation matrix (column-wise stacked mode shapes)

    Returns
    ---------------------------
    phi_max_real : boolean
        complex-valued modal transformation matrix, with vectors rotated to have maximum real parts
    """   

    angles = np.expand_dims(np.arange(0,np.pi/2, 0.01), axis=0)
    phi_max_real = np.zeros(np.shape(phi)).astype('complex')
    for mode in range(0,np.shape(phi)[1]):
        rot_mode = np.dot(np.expand_dims(phi[:, mode], axis=1), np.exp(angles*1j))
        max_angle_ix = np.argmax(np.sum(np.real(rot_mode)**2,axis=0), axis=0)

        phi_max_real[:, mode] = phi[:, mode] * np.exp(angles[0, max_angle_ix]*1j)*np.sign(sum(np.real(phi[:, mode])))

    return phi_max_real


def align_modes(phi):
    """
    Flip complex-valued or real-valued mode shapes such that similar modes have the same sign.

    Arguments
    ---------------------------
    phi : double
        complex-valued (or real-valued) modal transformation matrix (column-wise stacked mode shapes)

    Returns
    ---------------------------
    phi_aligned : boolean
        aligned complex-valued (or real-valued) modal transformation matrix
    """   

    prod = (phi.T @ phi)
    most_equal_ix = np.argmax(np.abs(np.sum(np.real(prod), axis=1)))
    modes_to_flip = np.real(prod)[:,most_equal_ix]<0
    phi_aligned = phi*1.0
    phi_aligned[:, modes_to_flip] = -phi[:,modes_to_flip]  # flip modes
    return phi_aligned


def normalize_phi(phi):
    """
    Normalize all complex-valued (or real-valued) mode shapes in modal transformation matrix.

    Arguments
    ---------------------------
    phi : double
        complex-valued (or real-valued) modal transformation matrix (column-wise stacked mode shapes)

    Returns
    ---------------------------
    phi_n : boolean
        modal transformation matrix, with normalized (absolute value of) mode shapes
    mode_scaling : 
        the corresponding scaling factors used to normalize, i.e., phi_n[:,n] * mode_scaling[n] = phi[n]

    """       
    phi_n = phi*0
    n_modes = np.shape(phi)[1]
    mode_scaling = np.zeros([n_modes])
    for mode in range(0, n_modes):
        mode_scaling[mode] = max(abs(phi[:, mode]))
        sign = np.sign(phi[np.argmax(abs(phi[:, mode])), mode])
        phi_n[:, mode] = phi[:, mode]/mode_scaling[mode]*sign

    return phi_n, mode_scaling


def mpc(phi):
    # Based on the current paper:
    # Pappa, R. S., Elliott, K. B., & Schenk, A. (1993). 
    # Consistent-mode indicator for the eigensystem realization algorithm. 
    # Journal of Guidance, Control, and Dynamics, 16(5), 852–858.

    # Ensure on matrix format
    if phi.ndim == 1:
        phi = phi[:,np.newaxis]

    n_modes = np.shape(phi)[1]
    mpc_val = [None]*n_modes

    for mode in range(0,n_modes):
        phin = phi[:, mode]
        Sxx = np.dot(np.real(phin), np.real(phin))
        Syy = np.dot(np.imag(phin), np.imag(phin))
        Sxy = np.dot(np.real(phin), np.imag(phin))

        eta = (Syy-Sxx)/(2*Sxy)

        lambda1 = (Sxx+Syy)/2 + Sxy*np.sqrt(eta**2+1)
        lambda2 = (Sxx+Syy)/2 - Sxy*np.sqrt(eta**2+1)

        mpc_val[mode] = ((lambda1-lambda2)/(lambda1+lambda2))**2

    mpc_val = np.array(mpc_val)
    return mpc_val