"""
##############################################
Clustering analysis module
##############################################
All functions related to the clustering of poles for automatic OMA.
"""

import numpy as np
# import hdbscan
from sklearn.cluster import HDBSCAN
from . import modal


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


def establish_tot_diff(lambd, phi, order, boolean_stops='default', scaling=None, normalize_distances=False, ensure_symmetric=True):
    """
    Establish total difference matrix based on input modes (from find_stable_modes).
        
    Arguments
    ---------------------------
    lambd : double
        array with complex-valued eigenvalues deemed stable
    phi : double
        2d array with complex-valued eigenvectors (stacked as columns), each column corresponds to a mode
    orders : int
        corresponding order for each stable mode
    boolean_stops : None, optional
        boolean stops to remove problematic poles (refers to difference matrices), e.g., to avoid same-order poles to appear
        in the same cluster the standard value {'order': [1, np.inf]} could be used (enforced by setting to 'avoid_same_order')
    scaling : {'mac': 1.0, 'lambda_real': 1.0, 'lambda_imag': 1.0}, optional
        scaling of predefined available variables used in total difference (available
        variables: 'mac', 'lambda_real', 'lambda_imag', 'omega_d', 'omega_n', 'order', 'xi')
    normalize_distances : True
        whether or not to normalize the distances (ensures compatibility between mac and other parameters)
    ensure_symmetric : True
        ensure that output is numerically symmetric

    Returns
    ---------------------------
    tot_diff : double
        cross-difference matrix (n_points-by-n_points)

    References
    ---------------------------
    Kvåle and Øiseth :cite:`Kvale2020`
    """
    
    # Establish dictionary with available difference variables
    diff_vars = dict()
    
    if type(boolean_stops)==str and boolean_stops == 'avoid_same_order':
        boolean_stops = {'order': [1, np.inf]}
        diff_vars['order'] = np.abs(crossdiff(order, relative=False))   #generates != integers?
                    
    elif boolean_stops is None:
         boolean_stops = {}  
    
    if scaling is None:
        scaling = {'mac': 1.0, 'lambda_real': 1.0, 'lambda_imag': 1.0}

    omega_n = np.abs(lambd)
    omega_d = np.abs(np.imag(lambd))
    xi = -np.real(lambd)/np.abs(lambd)


    if 'mac' in scaling:
        diff_vars['mac'] = np.abs(1.0 - modal.xmacmat(phi))
    
    if ('lambda_real' in scaling) or ('lambda_imag' in scaling):
        xlambda_diff = crossdiff(lambd, relative=True, allow_negatives=True)
        diff_vars['lambda_real']  = np.abs(np.real(xlambda_diff))
        diff_vars['lambda_imag']  = np.abs(np.imag(xlambda_diff))
        
    if 'omega_n' in scaling:
        diff_vars['omega_n'] = np.abs(crossdiff(omega_n, relative=True))
    
    if 'omega_d' in scaling:
        diff_vars['omega_d'] = np.abs(crossdiff(omega_d, relative=True))
        
    if 'order' in scaling:
        diff_vars['order'] = np.abs(crossdiff(order, relative=False))   #generates != integers?
        
    if 'xi' in scaling:
        diff_vars['xi']  = np.abs(crossdiff(xi, relative=True))

    # Normalize distances
    if normalize_distances:
        for key in scaling:
            diff_vars[key] = diff_vars[key]/np.max(np.abs(diff_vars[key])[:])

    # Establish boolean hard stop differences
    boolean_stop_diff = np.zeros([lambd.shape[0], lambd.shape[0]])
    
    for key in boolean_stops.keys():
        stops = boolean_stops[key]
        invalid_ix = np.logical_or((diff_vars[key]<stops[0]), (diff_vars[key]>stops[1]))
        boolean_stop_diff[invalid_ix] = np.inf
     
    # Establish total difference
    tot_diff = np.zeros([lambd.shape[0], lambd.shape[0]])

    for var in scaling:
        tot_diff += boolean_stop_diff + (diff_vars[var]*scaling[var])**2

    tot_diff = np.sqrt(tot_diff) + boolean_stop_diff

    if ensure_symmetric:
        tot_diff = (tot_diff + tot_diff.T)/2

    return tot_diff


class PoleClusterer:
    """
    Object to create pole clusters.
    
    Arguments
    ---------------------------
    lambd : double
        array with complex-valued eigenvalues deemed stable
    phi : double
        2d array with complex-valued eigenvectors (stacked as columns), each column corresponds to a mode
    orders : int
        corresponding order for each stable mode
    min_samples : 20, optional
        number of points in neighbourhood for point to be 
        considered core point (larger value => more conservative clustering)
    min_cluster_size : 20, optional
        when min_cluster_size points fall out of cluster it is not a split,
        it is merely points falling out of the cluster
    alpha : 1.0, optional
        distance scaling parameter, implies conservativism of clustering (higher => fewer points)
    boolean_stops : None, optional
        boolean stops to remove problematic poles (refers to difference matrices), e.g., to avoid same-order poles to appear
        in the same cluster the standard value {'order': [1, np.inf]} could be used (enforced by setting to 'avoid_same_order')
    scaling : {'mac': 1.0, 'lambda_real': 1.0, 'lambda_imag': 1.0}, optional
        scaling of predefined available variables used in total difference (available
        variables: 'mac', 'lambda_real', 'lambda_imag', 'omega_d', 'omega_n', 'order', 'xi')
    normalize_distances : False
        whether or not to normalize the distances (ensures compatibility between mac and other parameters)

    References
    ---------------------------
    Kvåle and Øiseth :cite:`Kvale2020`
    """


    def __init__(self, lambd, phi, order, min_samples=20, min_cluster_size=20, 
                 alpha=1.0, boolean_stops=None, scaling=None, normalize_distances=False,
                 run_clustering=True):
        
        self.boolean_stops = boolean_stops
        
        if scaling is None:
            self.scaling = {'mac': 1.0, 'lambda_real': 1.0, 'lambda_imag': 1.0}
        else:
            self.scaling = scaling

        self.normalize_distances = normalize_distances
        self.hdbscan_clusterer = HDBSCAN(metric='precomputed', 
                                                 min_samples=min_samples, 
                                                 min_cluster_size=min_cluster_size, 
                                                 alpha=alpha)
        self.lambd = lambd
        self.phi = phi
        self.order = order
        
        if run_clustering:
            self.cluster()


    def cluster(self):
        """
        Create tot_diff matrix and HDBSCAN cluster object from input data.
        """

        self.tot_diff = establish_tot_diff(self.lambd, self.phi, self.order, 
                                           boolean_stops=self.boolean_stops, scaling=self.scaling,
                                           normalize_distances=self.normalize_distances,
                                           ensure_symmetric=True)
        self.hdbscan_clusterer.fit(self.tot_diff)


    def postprocess(self, prob_threshold=0.0, normalize_and_maxreal=True):
        """
        Postprocess cluster object (sort and restrict).

        Arguments
        ---------------------------
        prob_threshold : 0.0, optional
            threshold value for probability of point belonging 
            to its determined cluster
        normalize_and_maxreal : True, optional
            whether or not to normalize each mode shape and maximize its real value 
            (rotate all components equally much in complex plane)

        Returns
        ---------------------------
        lambd_used : double
            sorted/remaining eigenvalues after restrictions/sort
        phi_used : double
            sorted/remaining eigenvectors after restrictions/sort
        order_stab_used : double
            corresponding orders
        group_ix : int
            indices (sorted based on damped natural freq.) of modes
        all_single_ix : double 
            index corresponding to input data
        probs : double
            probabilities of all points in all clusters
    
        """

        omega_d = np.abs(np.imag(self.lambd))

        if normalize_and_maxreal:
            phi0,__ = modal.normalize_phi(modal.maxreal(self.phi))
        else:
            phi0 = self.phi*1.0   
    

        # Establish all labels
        labels_all = self.hdbscan_clusterer.labels_

        # Align modes
        for label in np.unique(labels_all):
            phi0[:, labels_all==label] = modal.align_modes(phi0[:, labels_all==label])  #also align modes
            
        # Remove noise (label <0)
        keep_ix_temp = np.where(labels_all>=0)[0]

        # Apply probability threshold
        probs_temp = self.hdbscan_clusterer.probabilities_[keep_ix_temp]
        keep_ix = keep_ix_temp[probs_temp>=prob_threshold]
        probs_unsorted = self.hdbscan_clusterer.probabilities_[keep_ix]

        # Retain only "kept" indices from arrays
        labels_unsorted = labels_all[keep_ix]
        if len(labels_unsorted) == 0:
            return [], [], [], [], [], []
        
        n_labels = max(labels_unsorted)+1
            
        # Sort of cluster groups based on mean frequency
        wmean = [np.mean(omega_d[keep_ix][labels_unsorted==label]) for label in range(0, max(labels_unsorted)+1)]
        sort_ix = np.argsort(wmean)

        # Rearrange labels and probs (sorted based on frequency)
        labels = np.array([np.where(sort_ix==label)[0][0] for label in labels_unsorted]).flatten()
        probs = np.hstack([probs_unsorted[labels==label] for label in range(0, n_labels)])

        # Remove double results (more poles at same order within same cluster)
        keep_single_ix = [None]*n_labels
        for label in range(0, n_labels):           
            relevant_orders = self.order[keep_ix][labels==label]
            unique_relevant_orders = np.unique(relevant_orders)
            groups = [np.where(o==relevant_orders)[0] for o in unique_relevant_orders]
            keep_best_ix = []
            these_ix = np.where(labels==label)[0]
            
            for group in groups:
                keep_best_ix.append(these_ix[group[np.argmax(probs[labels==label][group])]])
            
            keep_single_ix[label] = np.array(keep_best_ix)

        all_single_ix = np.hstack(keep_single_ix)
        group_ixs = labels[all_single_ix]
        probs = probs[all_single_ix]
        order_stab_used = self.order[keep_ix][all_single_ix]
        lambd_used = self.lambd[keep_ix][all_single_ix]
        phi_used = phi0[:, keep_ix][:, all_single_ix]

        return lambd_used, phi_used, order_stab_used, group_ixs, all_single_ix, probs


def group_clusters(lambd_used, phi_used, order_stab_used, group_ixs, all_single_ixs, probs, return_lambda=False):
    '''
    Group the output of PoleClusterer.postprocess()

    Arguments
    ---------------------------
    lambd_used : double
        sorted/remaining eigenvalues after restrictions/sort
        to its determined cluster
    phi_used : double
        sorted/remaining eigenvectors after restrictions/sort
    order_stab_used : double
        corresponding orders
    group_ixs : int
        indices (sorted based on damped natural freq.) of modes
    all_single_ixs : double 
        index corresponding to input data
    probs : double
        probabilities of all points in all clusters
    return_lambda : boolean
        whether or not to return as lambda instead of xi, omega

    Returns
    ---------------------------
    xi_cluster : double
        list of arrays with xi grouped
    omega_n_cluster : double
        list of arrays with omega_n grouped
    phi_cluster : double
        list of arrays with phi grouped
    order_cluster : double
        list of arrays with orders grouped
    probs_cluster : double
        list of arrays with probs grouped
    ixs_cluster : double
        list of arrays with ixs corresponding to each cluster
    (lambd_cluster) : double, complex
        grouped clusters of lambd, only returned if return_lambda is True (returned in stead of xi_cluster and omega_n_cluster)

    '''  

    n_groups = len(np.unique(group_ixs))
    xi_cluster = [None]*n_groups
    omega_n_cluster = [None]*n_groups
    lambd_cluster = [None]*n_groups
    phi_cluster = [None]*n_groups
    order_cluster = [None]*n_groups
    probs_cluster = [None]*n_groups
    ixs_cluster = [None]*n_groups
    
    for group_ix in range(n_groups):
        this_ix = group_ixs==group_ix
        lambd_cluster[group_ix] = lambd_used[this_ix]
        xi_cluster[group_ix] = -np.real(lambd_used[this_ix])/np.abs(lambd_used[this_ix])
        omega_n_cluster[group_ix] = np.abs(lambd_used[this_ix])
        phi_cluster[group_ix] = phi_used[:, this_ix]
        order_cluster[group_ix] = order_stab_used[this_ix]
        probs_cluster[group_ix] = probs[this_ix]
        ixs_cluster[group_ix] = all_single_ixs[this_ix]
    
    if return_lambda:
        return lambd_cluster, phi_cluster, order_cluster, probs_cluster, ixs_cluster
    else:
        return xi_cluster, omega_n_cluster, phi_cluster, order_cluster, probs_cluster, ixs_cluster


def group_array(arr, group_ixs, axis=0):
    '''
    Group a single output array of PoleClusterer.postprocess() based on group indices.

    Arguments
    ---------------------------
    arr : double
        array corresponding
    group_ixs : int
        indices (sorted based on damped natural freq.) of modes

    Returns
    ---------------------------
    arr_grouped : double
        grouped array
    '''  

    n_groups = len(np.unique(group_ixs))
    arr_grouped = [None]*n_groups

    for group_ix in range(n_groups):
        this_ix = np.where(group_ixs==group_ix)[0]
        arr_grouped[group_ix] = np.take(arr, this_ix, axis=axis)

    return arr_grouped