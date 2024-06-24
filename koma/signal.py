from scipy.signal import csd
import numpy as np
from scipy.signal import correlate, correlation_lags
from scipy.signal import resample_poly, resample
from scipy.interpolate import interp1d

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


def estimate_lags(data, ref=0, upsample=None, fs=1.0, upsample_method='fourier'):
    '''
    Compute lags between all channels of data (channels stacked column-wise),
    based on max of correlation.

    Arguments
    -----------
    data : float
        data matrix (column-wise data)
    ref : int, default=0
        reference channel (to define resulting lags, will not affect analysis itself)

    upsample : int, optional
        upsample factor; if not specified the lag will be restricted to a resolution
        given by the original sampling factor
    fs : float, default=1.0
        sampling factor of data; if not given
    upsample_method : str, default='fourier'
        method used to conduct the upsampling {'fourier', 'linear', 'poly'}

    Returns
    -----------
    lags : float
        numpy 1d array with lags between all channels and reference channel;
        if fs is not given (i.e. set to its default value of 1.0), 
        the output defines the sample lag
    
    '''   
    if upsample is not None:
        if upsample_method == 'fourier':
            data = resample_poly(data, upsample, 1.0)
        elif upsample_method == 'linear':
            t = np.arange(0, (data.shape[0])*1/fs, 1/fs)
            fsi = fs*upsample
            ti = np.arange(0, (data.shape[0])*upsample*1/fsi, 1/fsi)
            data = interp1d(t, data, axis=0, fill_value='extrapolate')(ti)
        elif upsample_method == 'poly':
            data = resample(data, upsample)
    else:
        upsample = 1.0
    
    l = np.shape(data)[1]        
    R0 = np.zeros([l, l, 2*data.shape[0]-1])
    
    for dof1 in range(0, l):
        for dof2 in range(0, l):
            R0[dof1, dof2, :] = correlate(data[:, dof2], data[:,dof1], mode='full', method='auto')
   
    lag_length = correlation_lags(data[:, dof1].size, data[:, dof2].size, mode="full")   
    lags = lag_length[np.argmax(R0[ref, :, :], axis=-1)]
    
    return lags/(upsample*fs)


def shift_data(data, lags, cut=True):
    '''
    Shift data by specified lags.

    Arguments
    -----------
    data : float
        data matrix (column-wise data)
    lags : float
        list or numpy 1d array with sample lags (not time lag) to apply to all channels
    cut : boolean, default=True
        whether or not to cut the data to common valid range; if not,
        nans are used

    Returns
    -----------
    data_shifted : float
        shifted data matrix
    '''
    data_shifted = data*np.nan
    for ix, lag in enumerate(lags):
        x = np.arange(0, data.shape[0], 1.0)
        data_shifted[:, ix] = interp1d(x, data[:, ix], axis=0, bounds_error=False)(x+lag)
        
    if cut:
        nan_start = int(np.ceil(np.max([0,-np.min(lags)])))
        nan_end = int(-np.ceil(np.max(lags)))
        
        if nan_end == 0:
            nan_end = None 
        
        data_shifted = data_shifted[nan_start:nan_end, :]
        
    return data_shifted