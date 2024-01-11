"""
##############################################
Operational modal analysis module
##############################################
All functions related to operational modal analysis.
"""


import numpy as np
from scipy.signal import correlate
from scipy.linalg import matrix_balance, cholesky
from .modal import xmacmat_alt, mpc


def freq_svd(cpsd):
    """
    TODO: doc
    """
    D = cpsd*0
    U = cpsd*0
    
    for k in range(cpsd.shape[2]):
        U[:,:,k], D[:,:,k], __ = np.linalg.svd(cpsd[:,:,k], compute_uv=True)

    return U, D

def is_pos_def(x):
    """
    Check if matrix is positive definite.

    Arguments
    ---------------------------
    x : double
        matrix to be checked

    Returns
    ---------------------------
    pos_def : boolean
        positive definiteness of matrix
    """
    
    return np.all(np.linalg.eigvals(x) > 0)


def gaxpy_chol(M):
    """
    Gaxpy Cholesky decomposition. Required for Cov-SSI.

    Arguments
    ---------------------------
    M : double
        matrix to be decomposed

    Returns
    ---------------------------
    G : double
        lower-triangular Cholesky factor of M 
    
    Notes
    ---------------------------
    Conducts the decomposition:

    .. math::
        [M] = [G][G]^T
    """

    n = len(M)
    v = np.zeros([n]).astype('complex')
    A = M.astype('complex')

    for j in range(0, n):
        v[j:n] = A[j:n,j] - np.dot(A[j:n, 0:j], A[j, 0:j].T)
        A[j:n,j] = v[j:n]/np.sqrt(v[j])
        
    G = np.tril(A)
    
    return G


def gaxpy_chol_alt(M):
    """
    Alternative Gaxpy Cholesky decomposition.

    Arguments
    ---------------------------
    M : double
        matrix to be decomposed

    Returns
    ---------------------------
    G : double
        lower-triangular Cholesky factor of M 
    
    Notes
    ---------------------------
    Conducts the decomposition:

    .. math::
        [M] = [G][G]^T
    """

    n = np.shape(M)[0]
    G = np.zeros([n,n])
    s = np.zeros([n])

    for j in range(0, n):
        if j==0:
            s[j:n] = M[j:n, j]
        else:
            s[j:n] = M[j:n, j] - np.dot(G[j:n, 0:j], G[j, 0:j].T)
        
        G[j:n, j] = s[j:n]/np.sqrt(s[j])

    return G
    

def xcorr_lag(data, maxlag):
    """
    Computes cross-correlation matrices for multiple time lags based on input data.

    Arguments
    ---------------------------
    data : double
        matrix with data, n_samples-by-n_channels (channels are column-wise)
    maxlag : int
        maximum number of sample lags

    Returns
    ---------------------------
    R : double
        n_channels-by-n_channels-by-n_lags large array, each slice in third dimension 
        corresponds to the cross-correlation matrix for a given time lag 

    """

    l = np.shape(data)[1]
    nsamples = np.shape(data)[0]

    if maxlag>(nsamples-1):
        print('Maximum sample lag is larger than total sample length! Reducing to correspond to one below sample length.')
        maxlag = nsamples-1
    
    R0 = np.zeros([l, l, maxlag+1])
    
    for dof1 in range(0, l):
        for dof2 in range(0, l):
            R0[dof1, dof2, :] = correlate(data[:, dof2], data[:,dof1], mode='full', method='auto')[nsamples-1:nsamples+maxlag]    #matches MATLAB implementation - verify. dof1, dof2 and [(nsamples+1):(nsamples+1-maxlag):-1] looks more correct.
            
    unbiasing_scaling = (nsamples - np.linspace(0, maxlag, maxlag+1))[np.newaxis,np.newaxis,:]
    R = R0/unbiasing_scaling    # ensure unbiased estimate

    return R


def cva_weights(R, balancing_H0=None, balance=False):
    """
    Computes the weights for CVA.

    Arguments
    ---------------------------
    R : double
        n_channels-by-n_channels-by-n_lags large array, each slice in third dimension 
        corresponds to the cross-correlation matrix for a given time lag 
    balancing_H0 : boolean, optional
    balance : boolean, optional

    Returns
    ---------------------------
    W1 : double
        lower triangular Cholesky factor of R+ (see Hermans and van der Auweraer)
    W2 : double
        lower triangular Cholesky factor of R- (see Hermans and van der Auweraer)
    L1 : double
        inverse of W1
    L2 : double
        inverse of W2

    References
    --------------------------------
    Hermans and van Der Auweraer :cite:`Hermans1999`

    """

    i = int(np.shape(R)[2]/2)-1
    l = np.shape(R)[0]

    Rp = np.zeros([i*l, i*l])
    Rm = np.zeros([i*l, i*l])

    for row in range(0, i):
        for col in range(0, i):
            # R+
            if row<col: #above diagonal
                Rp[(row*l):(row*l+l), (col*l):(col*l+l)] = R[:, :, np.abs(row-col)].T
            elif row>=col:  #below or at diagonal
                Rp[(row*l):(row*l+l), (col*l):(col*l+l)] = R[:, :, np.abs(row-col)]

            # R-
            if row<=col:
                Rm[(row*l):(row*l+l), (col*l):(col*l+l)] = R[:, :, np.abs(row-col)]
            elif row>col:
                Rm[(row*l):(row*l+l), (col*l):(col*l+l)] = R[:, :, np.abs(row-col)].T   

    if balancing_H0 is not None:
        __, T = matrix_balance(balancing_H0)
        Rp = np.linalg.inv(T) @ Rp @ T
        Rm = np.linalg.inv(T) @ Rm @ T
    
    if balance:
        Rp, __ = matrix_balance(Rp)
        Rm, __ = matrix_balance(Rm)

    L1 = gaxpy_chol(Rp)            
    L2 = gaxpy_chol(Rm)
    W1 = np.linalg.inv(L1)
    W2 = np.linalg.inv(L2)

    return W1, W2, L1, L2


def stack_hankel(R):
    """
    Stacks lag-shifted cross-correlation matrices (stacked along third dimension) into Hankel matrix.

    Arguments
    ---------------------------
    R : double
        n_channels-by-n_channels-by-n_lags large array, each slice in third dimension 
        corresponds to the cross-correlation matrix for a given time lag 

    Returns
    ---------------------------
    H0 : double
        block-Hankel matrix
    H1 : double
        one-lag shifted block-Hankel matrix

    References
    --------------------------------
    Hermans and van Der Auweraer :cite:`Hermans1999`

    """

    i = int(np.shape(R)[2]/2)-1
    l = np.shape(R)[0]

    H0 = np.zeros([l*i, l*i])
    H1 = np.zeros([l*i, l*i])

    for row in range(0, i):             #go through all block rows from 1 to i
            for col in range(0, i):         #go through all block columns from 1 to i
                H0[(row*l):(row*l+l), (col*l):(col*l+l)] = R[:, :, row+col+1]           #Hankel matrix, H0
                H1[(row*l):(row*l+l), (col*l):(col*l+l)] = R[:, :, row+col+2]           #Lag-shifted Hankel matrix, H1    

    return H0, H1


def stack_toeplitz(R):
    """
    Stacks lag-shifted cross-correlation matrices (stacked along third dimension) into Toeplitz matrix.

    Arguments
    ---------------------------
    R : double
        n_channels-by-n_channels-by-n_lags large array, each slice in third dimension 
        corresponds to the cross-correlation matrix for a given time lag 

    Returns
    ---------------------------
    H0 : double
        block-Toeplitz matrix
    H1 : double
        one-lag shifted block-Toeplitz matrix

    References
    --------------------------------
    Rainieri and Fabbrocino :cite:`Rainieri`

    """

    i = int(np.shape(R)[2]/2)-1
    l = np.shape(R)[0]

    H0 = np.zeros([l*i, l*i])
    H1 = np.zeros([l*i, l*i])

    for row in range(0, i):             #go through all block rows from 1 to i
        for col in range(0, i):         #go through all block columns from 1 to i
            H0[(row*l):(row*l+l), (col*l):(col*l+l)] = R[:, :, i+1-col+row]            # Toeplitz matrix, H0
            H1[(row*l):(row*l+l), (col*l):(col*l+l)] = R[:, :, i+1-col+row+1]       

    return H0, H1


def covssi(data, fs, i, orders, weighting='none', matrix_type='hankel', 
algorithm='shift', showinfo=True, balance=True, return_A=False, discard_conjugates=True, return_flat=True):
    """
    Main function for covariance-driven SSI.

    Arguments
    ---------------------------
    data : double
        data matrix, with channels stacked column-wise (n_samples-by-n_channels)
    fs : double
        sampling frequency
    i : int
        maximum number of block rows
    orders : int
        array or list of what orders to include
    weighting : 'none', optional
        what weighting type to use ('none' or 'br', or 'cva')
    matrix_type : 'matrix_type', optional
        what matrix type formulation to base computation on ('hankel' or 'toeplitz')
    algorithm : 'shift', optional
        what algorithm to use ('shift' or 'standard' - when using 'standard' the method is equivalent to NExt-ERA)
    showinfo : True, optional
        whether or not to print information during computation
    balance : True, optional
        whether or not to conduct balancing to the cross-correlation matrices prior to matrix operations (Cholesky and SVD)
    return_A : False, optional
        whether or not to output the state matrix (discrete) from the last order evaluated
    discard_conjugates : True, optional
        whether or not to discard half of the poles as conjugates
    return_flat : True, optional
        whether or not to return flattened (directly plottable in stabplot)
    
    Returns
    ---------------------------
    lambd : double
        if `return_flat` is True, `lambd` is an array with complex-valued eigenvalues (one for each pole), 
        otherwise the variable is a list of arrays with complex-valued eigenvalues (each list entry correspond to one order)
    phi : double
        if `return_flat` is True, `phi` is a 2d array with complex-valued eigenvectors (stacked as columns), each column corresponds to a pole,
        otherwise the variable is a list of 2d arrays with complex-valued eigenvectors (each list entry correspond to one order), each column corresponds to a mode
    orders : int
        corresponding order for each stable mode - only output if `return_flat` is True  


    References
    --------------------------------
    * Hermans and van Der Auweraer :cite:`Hermans1999` (1 in code)
    * Van Overschee and de Moor :cite:`VanOverschee1996` (2 in code)
    * Rainieri and Fabbrocino :cite:`Rainieri` (3 in code)
    * Pridham and Wilson :cite:`Pridham2003` (4 in code)
    """ 

    if not (weighting in ['cva', 'br', 'none']):
        raise ValueError('Invalid weighting algorithm requested ("cva", "br", "none" accepted)')

    if showinfo:
        print('*** Covariance-driven SSI algorithm for OMA ***')

    dt = 1.0/fs
    l = np.shape(data)[1]

    # Establish Toeplitz/Hangel and lag-shifted Toeplitz/Hankel matrix [1]
    if showinfo:
        print('> Establishing Hankel/Toeplitz matrices')
        print('  >> Correlation estimation')

    R = xcorr_lag(data, 2*i+1)  #including no time lag entry, R_0 = R[:,:,0] - such that R_1-->R_2i+1 & R_0 makes 2i+2 in total

    # Matrix stacking
    if showinfo:
        print('  >> Matrix stacking')

    if matrix_type is 'hankel':
        H0, H1 = stack_hankel(R)
    elif matrix_type is 'toeplitz':
        H0, H1 = stack_toeplitz(R)
    else:
        raise ValueError('Unknown matrix string.')

    # Matrix stacking
    # Weighting
    W1, W2, L1, L2 = [np.eye(i*l)]*4 # Standard values

    if showinfo:
        print('> Establishing weighting matrices')
        print('  >> Weighting requested: %s' % (weighting.upper()))

    if weighting is 'cva':  # [2] used, [1] results in ill-conditioning issues
        try:
            if showinfo:
                print('  >> Establishing R+ and R-')
            
            W1, W2, L1, L2 = cva_weights(R, balance=balance)
            
            if algorithm is 'standard':
                print('Solution algorithm STANDARD is not supported with CVA. Changing to SHIFT algorithm.')
                algorithm='shift'
        except:
           print('  >> CVA failed! Continuing with no weighting (BR). If you want to run with CVA, try adjusting number of block rows.')

    if showinfo:
        print('> Computing SVD')
 
    U, d, Q = np.linalg.svd(W1 @ H0 @ W2.T)
    V = Q.T     # comes out transposed from svd function
    D = np.diag(d)  # make matrix

    # Eigenvalue solution
    if showinfo:
        print('> Computing state matrix for each order to establish modes')

    lambd = [None]*len(orders)
    phi = [None]*len(orders)

    for j in range(0, len(orders)):
        o = orders[j]

        # Select subset from H0, via the SVD decomposition, based on order
        D1 = D[0:o, :][:, 0:o]
        sqrtD1 = np.sqrt(D1)
        U1 = U[:, 0:o]
        V1 = V[:, 0:o]
        
        # State matrix estimate
        O = L1 @ U1 @ sqrtD1   #O, observability (T=identity matrix)
        C = O[0:l, :]      #Pick the l first values
        
        Oup = O[0:-l, :]   #Remove the last l rows from O
        Odown = O[l:, :]   #Remove the first l rows from O

        if algorithm is 'standard': #this means NExT-ERA is chosen
            A = np.linalg.inv(sqrtD1) @ U1.T @ H1 @ V1 @ np.linalg.inv(sqrtD1)        # [4] (and [3])
        elif algorithm is 'shift':
            A = np.linalg.pinv(Oup) @ Odown  # [1]

        # Eigenvalue decomposition [1] and convertion from discrete to continuous poles
        lambd_d, psi_d = np.linalg.eig(A)        #system eigenvectors and eigenvalues
        
        lambd_j = np.log(lambd_d)/dt    #make vector from diagonal matrix and transform from discrete to continuous
        phi_j = C @ psi_d                          #observed part of system eigenvectors, referring to input channels
        
        # Sort and order modes 
        sortix = np.argsort(np.abs(lambd_j))      #find index for unique absolute values
        if discard_conjugates:
            lambd[j] = lambd_j[sortix[::2]]                   #keep the corresponding eigenvalues
            phi[j] = phi_j[:, sortix[::2]]  
        else:
            lambd[j] = lambd_j[sortix]
            phi[j] = phi_j[:, sortix]                    #and the corresponding eigenvectors

    if showinfo:
        print('> Computation completed')
    
    if return_A:
        return A
    else:
        if return_flat:
            lambd_arr, phi_arr, orders_arr = flatten_stab_results(lambd, phi, orders)
            return lambd_arr, phi_arr, orders_arr
        else:
            return lambd, phi


def flatten_stab_results(lambd, phi, orders):
    """
    Flats out nested lambd and phi results.

    Arguments
    -----------------------
    lambd : double
        list of arrays with complex-valued eigenvalues (each list entry correspond to one order)
    phi : double
        list of 2d arrays with complex-valued eigenvectors (each list entry correspond to one order), each column corresponds to a mode
    orders : int
        array or list of what orders to include (each entry refers to each list entry in `lambd` and `phi`) 

    Returns
    -----------------------
    lambd : double
        array with complex-valued eigenvalues (one for each pole)
    phi : double
        2d array with complex-valued eigenvectors (stacked as columns), each column corresponds to a pole,
        otherwise the variable is a 
    orders : int
        corresponding order for each stable mode

    """

    orders_flat = np.hstack([[orders[ix]]*len(lambdi) for ix,lambdi in enumerate(lambd)])
    lambd_flat = np.hstack(lambd)
    phi_flat = np.hstack(phi)
    
    return lambd_flat, phi_flat, orders_flat

def find_stable_poles(lambd, phi, orders, s, stabcrit={'freq': 0.05, 'damping': 0.1, 'mac': 0.1}, 
                      valid_range={'freq': [0, np.inf], 'damping':[0, np.inf], 'mpc': [0,1]}, indicator='freq', 
                      return_both_conjugates=False, use_legacy=False):
     """
     Post-processing of Cov-SSI results, to establish modes (stable poles) from all poles.
    
     Arguments
     ---------------------------
     lambd : double
         array with complex-valued eigenvalues (one for each pole)
     phi : double
         2d array with complex-valued eigenvectors (stacked as columns), each column corresponds to a pole
     orders : int
         corresponding order for each stable mode   
     s : int
         stability level, see :cite:`Kvale2017_OMA`
     stabcrit : {'freq': 0.05, 'damping':0.1, 'mac': 0.1}, optional
         criteria to be fulfilled for pole to be deemed stable
     valid_range : {'freq': [0, np.inf], 'damping': [0, np.inf], 'mpc': [0,1]}, optional
         valid ranges of frequencies (rad/s) and damping for pole to be deemed stable
     indicator : 'freq', optional
         what modal indicator to use ('freq' or 'mac')
     use_legacy : False, optional
         False currently, legacy option will be removed altogether later

    
     Returns
     ---------------------------
     lambd_stab : double
         array with complex-valued eigenvalues deemed stable
     phi_stab : double
         2d array with complex-valued eigenvectors (stacked as columns), each column corresponds to a mode
     orders_stab : int
         corresponding order for each stable mode
     idx_stab : int
         indices of poles (within its order from the input) given in lambd_stab deemed stable (for indexing post)
         
     References
     --------------------------------
     Kvåle et al. :cite:`Kvale2017_OMA`
    
     """   
     
     if use_legacy:
         return find_stable_poles_legacy(lambd, phi, orders, s, stabcrit={'freq': 0.05, 'damping': 0.1, 'mac': 0.1}, 
                               valid_range={'freq': [0, np.inf], 'damping':[0, np.inf]}, indicator='freq', 
                               return_both_conjugates=False)
     
     def get_order(orders, order0, n=0):
         ix0 = np.where(np.unique(orders)==order0)[0]
         order = np.unique(orders)[ix0+n]
         return order
     
     wtol = stabcrit['freq']
     xitol = stabcrit['damping']
     mactol = stabcrit['mac']

     if 'freq' not in valid_range.keys():
         valid_range['freq'] = [0, np.inf]
     
     if 'damping' not in valid_range.keys():
         valid_range['damping'] = [0, np.inf]
         
     if 'mpc' not in valid_range.keys():
         valid_range['mpc'] = [0, 1.0]
         
     mpc_all = mpc(phi)
     

     lambd_stab = []
     phi_stab = []
     orders_stab = []
     idx_stab = []
     
     orders_unique = np.unique(orders)
     order_ix = orders == orders[0]
     
     # Establish for all orders above stablevel
     for order in orders_unique[s:]:
         order_ix = orders == order

         omega = np.abs(lambd[order_ix])
         xi = -np.real(lambd[order_ix])/np.abs(lambd[order_ix])
         phi_this = phi[:, order_ix]    # modal transformation matrix for order with index i
     
         # Stable poles
         for pole_ix in range(0, len(omega)):
             stab = 0
             for level in range(1, s+1):
                 level_order_ix = orders==get_order(orders, order, n=-level)
                 phi_last = phi[:, level_order_ix]
     
                 if indicator is 'mac':
                     macs = xmacmat_alt(phi_last, phi_this[:, pole_ix], conjugates=True)
                     pole_ix_last = np.argmax(macs[:,0])   #find largest mac in first column (all in matrix compared with vector)
                 elif indicator is 'freq':
                     omega_last = abs(lambd[level_order_ix])
                     pole_ix_last = np.argmin(abs(omega[pole_ix]-omega_last))
     
                 lambd_last = lambd[level_order_ix][pole_ix_last]
                 mpc_last = mpc_all[level_order_ix][pole_ix_last]
                 xi_last = -np.real(lambd_last)/abs(lambd_last)
                 
                 omega_last = abs(lambd_last)
                 dxi = abs(xi[pole_ix] - xi_last)
                 dw = abs(omega[pole_ix] - omega_last)
                 mac = xmacmat_alt(phi_last[:, pole_ix_last], phi_this[:, pole_ix], conjugates=True)
                 
                 if (((dw/omega_last)<=wtol) and ((dxi/xi_last)<=xitol) and 
                     (mac>=(1-mactol)) and (valid_range['freq'][0]<omega_last<valid_range['freq'][1]) 
                     and (valid_range['damping'][0]<xi_last<valid_range['damping'][1]) and
                     (valid_range['mpc'][0]<mpc_last<valid_range['mpc'][1])):
                     stab += 1
                 else:
                     stab = 0
                     break

             if stab>=s:
                 lambd_stab.append(lambd[order_ix][pole_ix])
                 phi_stab.append(phi[:,order_ix][:, pole_ix])
                 orders_stab.append(orders[order_ix][pole_ix])
                 idx_stab.append(pole_ix)
                 
     phi_stab = np.array(phi_stab).T
     
     lambd_stab = np.array(lambd_stab)
     orders_stab = np.array(orders_stab)
     idx_stab = np.array(idx_stab)
     phi_stab = np.array(phi_stab)
         
     return lambd_stab, phi_stab, orders_stab, idx_stab    

def find_stable_poles_legacy(lambd, phi, orders, s, stabcrit={'freq': 0.05, 'damping': 0.1, 'mac': 0.1}, 
                      valid_range={'freq': [0, np.inf], 'damping':[0, np.inf]}, indicator='freq', 
                      return_both_conjugates=False):
    """
    Post-processing of Cov-SSI results, to establish modes (stable poles) from all poles.

    Arguments
    ---------------------------
    lambd : double
        list of arrays with complex-valued eigenvalues, one list entry per order
    phi : double
        list of 2d arrays with complex-valued eigevectors (stacked column-wise), one list entry per order
    orders : int
        orders corresponding to the list entries in lambd and phi
    s : int
        stability level, see :cite:`Kvale2017_OMA`
    stabcrit : {'freq': 0.05, 'damping':0.1, 'mac': 0.1}, optional
        criteria to be fulfilled for pole to be deemed stable
    valid_range : {'freq': [0, np.inf], 'damping':[0, np.inf]}, optional
        valid ranges of frequencies (rad/s) and damping for pole to be deemed stable
    indicator : 'freq', optional
        what modal indicator to use ('freq' or 'mac')
    return_both_conjugates : boolean
        whether or not to return both conjugates of each pole

    Returns
    ---------------------------
    lambd_stab : double
        array with complex-valued eigenvalues deemed stable
    phi_stab : double
        2d array with complex-valued eigenvectors (stacked as columns), each column corresponds to a mode
    orders_stab : int
        corresponding order for each stable mode
    idx_stab : int
        indices of poles (within its order from the input) given in lambd_stab deemed stable (for indexing post)
        
    References
    --------------------------------
    Kvåle et al. :cite:`Kvale2017_OMA`

    """ 

    wtol = stabcrit['freq']
    xitol = stabcrit['damping']
    mactol = stabcrit['mac']

    if 'freq' not in valid_range.keys():
        valid_range['freq'] = [0, np.inf]
    
    if 'damping' not in valid_range.keys():
        valid_range['damping'] = [0, np.inf]

    lambd_stab = []
    phi_stab = []
    orders_stab = []
    idx_stab = []
    
    # Establish for all orders above stablevel
    for order_ix in range(s-1, len(orders)):   # i is order
        omega = abs(lambd[order_ix])
        xi = -np.real(lambd[order_ix])/abs(lambd[order_ix])
        phi_this = phi[order_ix]    # modal transformation matrix for order with index i
    
        # Stable poles
        for pole_ix in range(0, len(omega)):
            stab = 0
            for level in range(1, s+1):
                phi_last = phi[order_ix-level]
    
                if indicator is 'mac':
                    macs = xmacmat_alt(phi_last, phi_this[:, pole_ix], conjugates=True)
                    pole_ix_last = np.argmax(macs[:,0])   #find largest mac in first column (all in matrix compared with vector)
                elif indicator is 'freq':
                    omega_last = abs(lambd[order_ix-level])
                    pole_ix_last = np.argmin(abs(omega[pole_ix]-omega_last))
    
                lambd_last = lambd[order_ix-level][pole_ix_last]
                xi_last = -np.real(lambd_last)/abs(lambd_last)
                omega_last = abs(lambd_last)
                dxi = abs(xi[pole_ix] - xi_last)
                dw = abs(omega[pole_ix] - omega_last)
                mac = xmacmat_alt(phi_last[:, pole_ix_last], phi_this[:, pole_ix], conjugates=True)
    
                if ((dw/omega_last)<=wtol) and ((dxi/xi_last)<=xitol) and (mac>=(1-mactol)) and (valid_range['freq'][0]<omega_last<valid_range['freq'][1]) and (valid_range['damping'][0]<xi_last<valid_range['damping'][1]):
                    stab += 1
                else:
                    stab = 0
                    break

            if stab>=s:
                lambd_stab.append(lambd[order_ix][pole_ix])
                phi_stab.append(phi[order_ix][:, pole_ix])
                orders_stab.append(orders[order_ix])
                idx_stab.append(pole_ix)
                
    phi_stab = np.array(phi_stab).T
    
    lambd_stab = np.array(lambd_stab)
    orders_stab = np.array(orders_stab)
    idx_stab = np.array(idx_stab)
    phi_stab = np.array(phi_stab)

    if not return_both_conjugates:
        lambd_stab = lambd_stab[::2]
        phi_stab = phi_stab[:, ::2]
        orders_stab = orders_stab[::2]
        idx_stab = idx_stab[::2]
        
    return lambd_stab, phi_stab, orders_stab, idx_stab