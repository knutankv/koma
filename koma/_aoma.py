from .oma import covssi
from .oma import find_stable_poles

class StabCrit:
    def __init__(self, freq=None, damping=None, mac=None, s=1, valid_ranges=None):
        if valid_ranges is None:
            valid_range = dict(freq=[0, np.inf], damping=[0, np.inf])

        self.valid_range = valid_range
        self.mac = mac
        self.damping = damping
        self.freq = freq
        self.s = s     

class CovSSI:
    def __init__(self, settings=None, data=None, fs=None, add_data_to_obj=False):
        if (data is not None) and (settings is not None):
            if fs is None:
                print('No analysis conducted as sampling frequency (fs) or settings are not provided')
            else:
                self.run(data, fs)
                            
            if add_data_to_obj is True:
                self.data = data
                self.fs = fs

    def run(self, data, fs):
        lambd, phi = covssi(data,fs,self.blockrows,self.orders,weighting)
        self.results = CovSSIResults(lambd=lambd, phi=phi, orders=self.orders)
        if self.settings.stab_crit is not None:
            lambd_stab, phi_stab, order_stab, idx_stab = establish_stable_poles
        
    def establish_stable_poles(self, stab_crit):
        lambd_stab, phi_stab, order_stab, idx_stab = find_stable_poles(lambd_ssi, phi_ssi, orders, s, stabcrit=stabcrit, valid_range=valid_range, indicator='mac') # find stable poles
        order_stab = np.array(order_stab)
        self.results.stable_results = dict(lambd=lambd_stab, phi=phi_stab, orders=order_stab)

class CovSSIResults:
    def __init__(self, lambd=None, phi=None, orders=None):
        self.phi = phi
        self.lambd = lambd
        self.orders = orders
        self.stable_results = dict()

class CovSSISettings:
    def __init__(self, blockrows, orders, weighting='br', matrix_type='hankel', algorithm='shift', stab_crit=None):
        self.blockrows = blockrows
        self.orders = orders
        self.weighting = weighting
        self.matrix_type = matrix_type
        self.algorithm = algorithm
        self.stab_crit = stab_crit

class ClusterPoles:
    def __init__(self, ):

