from scipy.signal import csd
import numpy as np

def xwelch(x, **kwargs):
    '''
    Compute Welch cross-spectral density matrix for given inputs (stacked column-wise).

    Arguments
    -----------
    x : float
        data matrix (column-wise data)
    **kwargs 
        passed on to scipy's csd function for cross-spectral density estimation

    Returns
    -----------
    f : float
        numpy array of frequencies
    cpsd : float
        3d numpy array with CPSD matrix (frequency component on third index)
    
    '''
    f, __ = csd(x[:,0], x[:,0], **kwargs)
    cpsd = np.zeros([x.shape[1], x.shape[1], len(f)]).astype('complex')

    for i, xi in enumerate(x.T):
        for j, xj in enumerate(x.T):
            f, cpsd[i,j,:] = csd(xi, xj, **kwargs)

    return f, cpsd