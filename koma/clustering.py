"""
##############################################
Clustering analysis module
##############################################
All functions related to the clustering of poles for automatic OMA.
"""


import numpy as np

def crossdiff(arr, relative=False, allow_negatives=False):
    """
    Establish cross difference matrix used for clustering analysis.

    Arguments
    ---------------------------
    arr : double
        array to provide cross difference of (n_points long)
    relative : False, optional
        output relative difference or absolute difference
    allow_negatives : False, optional
        whether or not to allow negative difference

    Returns
    ---------------------------
    diff : double
        cross-difference matrix (n_points-by-n_points)

    References
    ---------------------------
    Kvåle and Øiseth :cite:`Kvale2020`
    """

    arr1, arr2 = np.meshgrid(arr, arr)

    if relative:
        scaling = arr1
    else:
        scaling = arr1*0+1.0

    if allow_negatives:
        diff = np.minimum(np.real((arr1-arr2)/scaling), np.real((arr1+arr2)/scaling)) + np.minimum(np.imag((arr1-arr2)/scaling), np.imag((arr1+arr2)/scaling))*1j
    else:
        diff = (arr1-arr2)/scaling

    # Remove imaginary 0 if input is float (real)
    if ~np.any(np.iscomplex(arr)):
        diff = np.real(diff).astype('float')

    return diff

