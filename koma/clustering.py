"""
##############################################
Clustering analysis module
##############################################
All functions related to the clustering of poles for automatic OMA.
"""


import numpy as np
import hdbscan
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


def establish_tot_diff(lambd, phi, order, boolean_stops='default', scaling=None):
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
    boolean_stops : 'default', optional
        boolean stops to remove problematic poles (refers to difference matrices), e.g., to avoid same-order poles to appear
        in the same cluster the standard value {'order': [1, np.inf]} could be used
    scaling : {'mac': 1.0, 'lambda_real': 1.0, 'lambda_imag': 1.0}, optional
        scaling of predefined available variables used in total difference (available
        variables: 'mac', 'lambda_real', 'lambda_imag', 'omega_d', 'omega_n', 'order', 'xi')

    Returns
    ---------------------------
    tot_diff : double
        cross-difference matrix (n_points-by-n_points)

    References
    ---------------------------
    Kvåle and Øiseth :cite:`Kvale2020`
    """

    if boolean_stops is 'default':
        boolean_stops = {'order': [1, np.inf]}
                    
    elif boolean_stops is None:
         boolean_stops = {}  
    
    if scaling is None:
        scaling = {'mac': 1.0, 'lambda_real': 1.0, 'lambda_imag': 1.0}

    omega_n = np.abs(lambd)
    omega_d = np.abs(np.imag(lambd))
    xi = -np.real(lambd)/np.abs(lambd)

    # Establish dictionary with available difference variables
    diff_vars = dict()

    diff_vars['mac'] = np.abs(1.0 - modal.xmacmat(phi))
    xlambda_diff = crossdiff(lambd, relative=False, allow_negatives=True)
    diff_vars['lambda_real']  = np.abs(np.real(xlambda_diff))
    diff_vars['lambda_imag']  = np.abs(np.imag(xlambda_diff))
    diff_vars['omega_n'] = np.abs(crossdiff(omega_n, relative=True))
    diff_vars['omega_d'] = np.abs(crossdiff(omega_d, relative=True))
    diff_vars['order'] = np.abs(crossdiff(order, relative=False))   #generates != integers?
    diff_vars['xi']  = np.abs(crossdiff(xi, relative=True))

    # Establish boolean hard stop differences
    boolean_stop_diff = np.zeros(np.shape(diff_vars['xi']))
    
    for key in boolean_stops.keys():
        stops = boolean_stops[key]
        invalid_ix = np.logical_or((diff_vars[key]<stops[0]), (diff_vars[key]>stops[1]))
        boolean_stop_diff[invalid_ix] = np.inf
     
    # Establish total difference
    tot_diff = np.zeros(np.shape(diff_vars['xi']))

    for var in scaling:
        tot_diff += boolean_stop_diff + (diff_vars[var]*scaling[var])**2

    tot_diff = np.sqrt(tot_diff) + boolean_stop_diff
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
    boolean_stops : 'default', optional
        boolean stops to remove problematic poles (refers to difference matrices), e.g., to avoid same-order poles to appear
        in the same cluster the standard value {'order': [1, np.inf]} could be used
    scaling : {'mac': 1.0, 'lambda_real': 1.0, 'lambda_imag': 1.0}, optional
        scaling of predefined available variables used in total difference (available
        variables: 'mac', 'lambda_real', 'lambda_imag', 'omega_d', 'omega_n', 'order', 'xi')

    Returns
    ---------------------------
    pole_clusterer : object
        clusterer object

    References
    ---------------------------
    Kvåle and Øiseth :cite:`Kvale2020`
    """

    def __init__(self, lambd, phi, order, min_samples=20, min_cluster_size=20, alpha=1.0, boolean_stops='default', scaling=None):
        self.boolean_stops = boolean_stops
        self.scaling = None
        self.clusterer = hdbscan.HDBSCAN(metric='precomputed', min_samples=min_samples, min_cluster_size=min_cluster_size, alpha=alpha, gen_min_span_tree=False)
        self.lambd = lambd
        self.phi = phi
        self.order = order

    def cluster(self):
        """
        Create tot_diff matrix and HDBSCAN cluster object from input data.
        """

        self.tot_diff = establish_tot_diff(self.lambd, self.phi, self.order, boolean_stops=self.boolean_stops, scaling=self.scaling)
        self.clusterer.fit(self.tot_diff)

    def postprocess(self, prob_threshold=0.0):
        """
        Postprocess cluster object (sort and restrict).

        Arguments
        ---------------------------
        prob_threshold : 0.0, optional
            threshold value for probability of point belonging 
            to its determined cluster

        Returns
        ---------------------------
        lambd_used : double
            sorted/remaining eigenvalues after restrictions/sort
        phi_used : double
            sorted/remaining eigenvectors after restrictions/sort
    
        """

        omega_d = np.abs(np.imag(self.lambd))
        phi0,__ = modal.normalize_phi(np.real(modal.maxreal(self.phi)))

        # Establish all labels
        labels_all = self.clusterer.labels_

        # Align modes
        for label in np.unique(labels_all):
            phi0[:, labels_all==label] = modal.align_modes(phi0[:, labels_all==label])  #also align modes
            
        # Remove noise (label <0)
        keep_ix_temp = np.where(labels_all>=0)[0]

        # Apply probability threshold
        probs_temp = self.clusterer.probabilities_[keep_ix_temp]
        keep_ix = keep_ix_temp[probs_temp>=prob_threshold]
        probs_unsorted = self.clusterer.probabilities_[keep_ix]

        # Retain only "kept" indices from arrays
        labels_unsorted = labels_all[keep_ix]
        n_labels = max(labels_unsorted)+1
            
        # Sort of cluster groups based on mean frequency
        wmean = [np.mean(omega_d[keep_ix][labels_unsorted==label]) for label in range(0, max(labels_unsorted)+1)]
        sort_ix = np.argsort(wmean)
        wmean_sort = np.sort(wmean)

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
        labels = labels[all_single_ix]
        probs = probs[all_single_ix]

        order_stab_used = self.order[keep_ix][all_single_ix]
        lambd_used = self.lambd[keep_ix][all_single_ix]

        phi_used = phi0[:, keep_ix][:, all_single_ix]

        return lambd_used, phi_used