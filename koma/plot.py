"""
##############################################
Visualization module
##############################################
All functions related to the plotting of results from OMA methods.
"""
import plotly.graph_objects as go 
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.backend_bases import MouseButton
from matplotlib import cm
import warnings
from .modal import mpc

from matplotlib.offsetbox import AnnotationBbox, TextArea

class Selector:
    def __init__(self, ax, active_settings={}, deactive_settings={}, templine_settings={}, picked_settings={}):
        self.ax = ax
        self.fig = ax.get_figure()
        self.lines = list(ax.lines)
        self.picked = []
        self.picked_lines = []
        self.index = 0
        self.templine_settings = dict(ls='-', lw=0.5) | templine_settings

        if len(self.lines)>0:
            self.x = self.lines[0].get_xdata()

        self.active_settings = {'linewidth': 2.0, 'alpha':1.0} | active_settings
        self.deactive_settings = {'linewidth': 1.0, 'alpha': 0.5}  | deactive_settings
        self.picked_settings = dict(ls='--', lw=1) | picked_settings

    @property
    def n_lines(self):
        return len(self.lines)
    
    @property
    def current_line(self):
        return self.lines[self.index]
    
    def on_click(self, event):
            if event.inaxes and event.button is MouseButton.LEFT and self.fig.canvas.manager.toolbar.mode.value == '':
                color = self.current_line.get_color()
                ix_sel = np.argmin(np.abs(event.xdata - self.x))

                self.picked.append([self.index, ix_sel, self.x[ix_sel]])
                self.picked_lines.append(self.ax.axvline(self.x[ix_sel], color=color, **self.picked_settings))

            elif event.inaxes and event.button is MouseButton.RIGHT:
                ix_sel = np.argmin(np.abs(event.xdata - self.x))
                picked_mat = np.vstack(self.picked)
                valid_x_ix = np.where((picked_mat[:,0] == self.index))[0]

                if len(valid_x_ix)>0:
                    row_ix = np.argmin(np.abs(ix_sel-picked_mat[valid_x_ix, 1]))
                    ix_remove = valid_x_ix[row_ix]

                    self.picked.pop(ix_remove)
                    self.picked_lines[ix_remove].remove()
                    self.picked_lines.pop(ix_remove)

            self.fig.canvas.draw()

    def on_scroll(self, event):
        increment = 1 if event.button == 'up' else -1
        self.index = np.clip(self.index + increment, 0, self.n_lines-1)
        self.x = self.lines[self.index].get_xdata()
        self.templine.set_color(self.current_line.get_color())
        self.update()
    
    def update(self):
        for ix,line in enumerate(self.lines):
            if ix==self.index:
                line.set(**self.active_settings)
            else:
                line.set(**self.deactive_settings) 

        self.fig.canvas.draw()

    def on_close(self, event):
        self.fig.canvas.stop_event_loop()

    def on_keyboard(self, event):
        if event.key == 'enter':
            self.title_format = 'ix = {index}'
            self.update()
            self.fig.canvas.stop_event_loop()
            
            # Prepare for printing
            self.templine.set_xdata([np.nan])
            for ix,line in enumerate(self.lines):
                line.set(**self.deactive_settings) 
                line.set(alpha=1.0)

    def on_move(self, event):
        if event.inaxes:
            ix_sel = np.argmin(np.abs(event.xdata-self.x))
            self.templine.set_xdata([self.x[ix_sel]])
            self.fig.canvas.draw()

    def get_fig(self, show=True, block=True):
        '''
        Get an interactive figure.
        
        Parameters
        -----------
        show : bool, default=True
            whether or not to automatically show figure
        block : bool, default=True
            whether or not to block before returning to retain interactivity
            
        Returns
        -----------
        fig : `matplotlib.Figure` object
            figure object
        '''
        
        self.ax.legend(frameon=False)
        self.templine = self.ax.axvline(np.nan, **self.templine_settings)
        self.templine.set_color(self.current_line.get_color())

        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)   
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.fig.canvas.mpl_connect('key_press_event', self.on_keyboard)
        self.fig.canvas.mpl_connect('close_event', self.on_close)
        
        if block:
            self.title_format = 'ix = {index} (Press enter when selection is done)'

        self.update()   

        if show:
            plt.show()
            if block:
                self.fig.canvas.start_event_loop()

        return self.fig
    
    @property
    def all_picked_x(self):
        return np.array([p[2] for p in self.picked])
    
    @property
    def picked_x(self):
        pmat = np.vstack(self.picked)
        x = [None]*self.n_lines
        for ix in range(self.n_lines):
            okix = pmat[:,0] == ix
            x[ix] = pmat[okix, 2]

        return x
    
    @property
    def picked_ix(self):
        pmat = np.vstack(self.picked)
        x_ix = [None]*self.n_lines
        for ix in range(self.n_lines):
            okix = pmat[:,0] == ix
            x_ix[ix] = pmat[okix, 1].astype(int)

        return x_ix
    
    @property
    def picked_line_ix(self):
        return np.array([p[0] for p in self.picked])
    
    @property
    def all_picked_ix(self):
        return np.array([p[1] for p in self.picked])

class StabPlotter:
    '''
    Class to initialize interactive stabilization plot.
    
    Example
    -----------
    
    Assuming we have a data matrix `data` sampled at `fs=128.0`. The data is first processed
    using the covssi (can also be filtered using the `find_stable_poles` function):
        
        fs = 128.0
        i = 30
        orders_input = np.arange(2, 100, 2)
        lambd, phi, orders = koma.oma.covssi(data, fs, i, orders_input)

    
    The stabilization plot is generated like this:
        
        stab_plotter = koma.plot.StabPlotter(lambd, orders, phi=phi, freq_unit='hz', damped_freq=False, annotate_hover=True)
        fig = stab_plotter.get_fig()
    
    The following code extracts the data from the plotter object as variables:
    
        phi_sel, fn_sel, xi_sel = stab_plotter.phi, stab_plotter.fn, stab_plotter.xi
        '
    '''
    
    def __init__(self, lambd, orders, phi=None, freq_unit='rad/s', 
                damped_freq=False, psd_freq=None, psd_y=None, log_psd_scale=True, 
                pole_settings=None, selected_pole_settings=None, hover_pole_settings=None, ax=None, num=None, 
                sort_by='undamped', annotate_hover=False, psd_color='gray'):
           
        '''
        Parameters
        ------------
  
        lambd : double
            array with complex-valued eigenvalues
        orders : int
            corresponding order for each pole in `lambd`
        phi : double, optional
            matrix where each column is complex-valued eigenvector corresponding to lambd
        freq_unit : str, default='rad/s'
            what frequency unit to use; 'Hz' or 'rad/s'
        damped_freq : False, optional
            whether or not to use damped frequency (or period) values in plot (False enforces undamped freqs)
        psd_freq : double, optional
            frequency values of plot to overlay, typically spectrum of data
        psd_y : double, optional
            function values of plot to overlay, typically spectrum of data
        log_psd_scale: boolean, default=True
            whether or not to plot the overlaid PSD using a logarithmic y-scale
        pole_settings : dict
            dictionary with settings to pass to the plot settings of the poles
        selected_pole_settings : dict
            dictionary with settings to pass to the plot settings of the selected poles
        hover_pole_settings : dict
            dictionary with settings to pass to the plot settings of the pole currently
            being hovered
        ax : `matplotlib.Axis` object
            axis to place plot in; if not given, a new axis in the figure specified
            will be created
        num : int, optional
            figure number used; only used if `ax` = None
        sort_by : str, default='undamped'
            what quantity to sort output by; either 'undamped', 'damped' or None
        psd_color : str, default='gray'
            color to use for PSD plot overlayed

        Returns
        ---------------------------
        fig : obj
            plotly figure object
        
        '''
        if pole_settings is None: pole_settings = {}
        if selected_pole_settings is None: selected_pole_settings = {}
        if hover_pole_settings is None: hover_pole_settings = {}
    
        self._lambd = lambd
        self._orders = orders
        self._phi = phi
        self.sort_by = sort_by

        self.damped_freq = damped_freq
        
        # Make list
        if psd_y is not None and not isinstance(psd_y, list):
            psd_y = [psd_y]
            psd_freq = [psd_freq]
        
        self.psd_color = psd_color
        self.psd_freq = psd_freq
        self.psd_y = psd_y
        self.log_psd_scale = log_psd_scale
        
        self.pole_settings = dict(linestyle='none', marker='.', color='k') | pole_settings
        self.selected_pole_settings = dict(linestyle='none', marker='o', color='r') | selected_pole_settings
        self.hover_pole_settings = dict(linestyle='none', marker='o', color='r', alpha=0.3) | hover_pole_settings
        self.annotate_hover = annotate_hover

        self._picked = []
        self._hoverpos = [np.nan, np.nan]
        self._hoverix = None
        
        self._hoverdot = None
        self._annotation = None
        
        self.picked_dots = []

        if damped_freq:
            dampedornot = 'd'
            self._f = np.abs(np.imag(lambd))
        else:
            dampedornot = 'n'
            self._f = np.abs(lambd)

        if freq_unit.lower() == 'hz':
            self._f = self._f/2/np.pi
            self.freq_name = fr'$f_{dampedornot}$'
            self.freq_unit = 'Hz'
        else:
            self.freq_name = fr'$\omega_{dampedornot}$'
            self.freq_unit = 'rad/s'

        if ax is None:
            plt.figure(num).clf()
            self.fig, self.ax = plt.subplots(num=num, figsize=(18,7))
        else:
            plt.sca(ax)
            self.ax = ax
            self.fig = plt.gcf()
 
    @property
    def hoverpos(self):
        return self._hoverpos
    
    @hoverpos.setter
    def hoverpos(self, ix):
        self._hoverix = ix
        
        if self._hoverix is None:
            x = np.nan
            y = np.nan
            text = ''
            self._annotation.set_visible(False)
        else:
            x = self._f[ix]
            y = self._orders[ix]
            text = (f'ix = {self._hoverix}\n' +
                    f'n = {self._orders[ix]}\n' + 
                    fr'{self.freq_name} = {self._hoverpos[0]:.2f} {self.freq_unit}' + '\n' +
                    fr'$\xi$ = {self.get_xi(ix)*100:.2f}%')
            
            if self._phi is not None:
                this_mpc = mpc(self._phi[:,ix:ix+1])[0]
                text = text + f'\nMPC={this_mpc*100:.1f}%'
            
        self._hoverpos = x,y
        self._hoverdot.set_xdata([x])
        self._hoverdot.set_ydata([y])
        
        if self.annotate_hover:
            self._annotation.offsetbox.set(text=text)
            self._annotation.xy = self._hoverpos
    
    @property
    def picked(self):   #sorted picked
        if self.sort_by is None:
            ix = None
        elif self.sort_by == 'undamped':
            ix = np.argsort(np.abs(self._lambd[self._picked]))
        elif self.sort_by == 'damped':
            ix = np.argsort(np.abs(np.imag(self._lambd[self._picked])))

        return np.array(self._picked)[ix]
    
    @property
    def ix(self):   #alias
        return self.picked

    # Picked orders
    @property
    def n_picked(self):
        if len(self.picked)>0:
            return self._orders[self.picked]
        else:
            return np.empty([0])

    # Eigenvalues and eigenvectors
    @property
    def lambd(self):
        if len(self.picked)>0:
            return self._lambd[self.picked]
        else:
            return np.empty([0])

    @property
    def phi(self):
        if len(self.picked)>0:
            return self._phi[:, self.picked]
        else:
            return np.empty([0])
    
    # Damped natural freqs
    @property
    def wd(self):
        return np.abs(np.imag(self.lambd))
    
    @property
    def omegad(self):
        return self.wd
    
    @property
    def fd(self):
        return self.wd/2/np.pi

    # Undamped natural freqs
    @property
    def wn(self):
        return np.abs(self.lambd)
    
    @property
    def omegan(self):
        return self.wn
    
    @property
    def fn(self):
        return self.wn/2/np.pi
    
    # Damping
    @property
    def xi(self):
        return -np.real(self.lambd)/np.abs(self.lambd)
    
    def get_xi(self, ix=None):
        xi = -np.real(self._lambd[ix])/np.abs(self._lambd[ix])
        return xi

    def get_df(self, pars=['ix', 'n_picked', 'wn', 'xi']):
        '''
        Get pandas dataframe with results.
        '''
        
        df = pd.DataFrame(data=np.vstack([getattr(self, par) for par in pars]).T, 
        columns=pars)

        if 'n_picked' in df:
            df['n_picked'] = df['n_picked'].astype(int)
        
        if 'picked' in df:
            df['picked'] = df['picked'].astype(int)

        if 'ix' in df:
            df['ix'] = df['ix'].astype(int)

        return df

    
    def get_ix(self, event):
        dist = (event.xdata - self._f)**2 + (event.ydata-self._orders)**2
        ix_sel = np.argmin(dist)
        return ix_sel
    

    def get_fig(self, show=True, block=True):
        if self.psd_y is not None:
            ax2 = self.ax.twinx()
            for ix, (psd, fi) in enumerate(zip(self.psd_y, self.psd_freq)):
                ax2.plot(fi, psd, self.psd_color)

            if self.log_psd_scale:
                ax2.set_yscale('log')

        self.ax.plot(self._f, self._orders, **self.pole_settings)
        self.ax.set_xlabel(f'{self.freq_name} [{self.freq_unit}]')
        self.ax.set_ylabel('Order, $n$')
        self._hoverdot = self.ax.plot([np.nan], [np.nan], **self.hover_pole_settings)[0]

        if self.annotate_hover:
            offsetbox = TextArea('')
            self._annotation = AnnotationBbox(offsetbox, [np.nan, np.nan],
                    xybox=(80, 0),
                    xycoords='data',
                    boxcoords="offset points", 
                    arrowprops=dict(arrowstyle="->"),
                    pad=0.3, bboxprops=dict(alpha=0.85))
            
            self.ax.add_artist(self._annotation)
            
        if self.psd_y is not None:
            self.ax.set_zorder(ax2.get_zorder() + 1)

        self.ax.patch.set_visible(False)

        self.fig.canvas.mpl_connect('button_press_event', self.on_click)   
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.fig.canvas.mpl_connect('key_press_event', self.on_keyboard)
        self.fig.canvas.mpl_connect('close_event', self.on_close)
        
        if block:
            self.title_format = '(Press enter when selection is done)'

        self.update()   

        if show:
            plt.show()
            if block:
                self.fig.canvas.start_event_loop()

        return self.fig
    
    # Interaction methods
    def on_click(self, event):
        if event.inaxes:
            ix_sel = self.get_ix(event)
           
            if (event.button is MouseButton.LEFT and self.fig.canvas.manager.toolbar.mode.value == ''):
                self._picked.append(ix_sel)
                self.picked_dots.append(
                    self.ax.plot(self._f[ix_sel], self._orders[ix_sel], 
                                 **self.selected_pole_settings)[0]
                                       )

            elif event.button is MouseButton.RIGHT:
                if ix_sel in self._picked:
                    ix_remove = self._picked.index(ix_sel)

                    self._picked.pop(ix_remove)
                    self.picked_dots[ix_remove].remove()
                    self.picked_dots.pop(ix_remove)

        self.fig.canvas.draw()

    def on_close(self, event):
        self.fig.canvas.stop_event_loop()

    def on_keyboard(self, event):
        if event.key == 'enter':
           
            # Prepare for printing
            self.hoverpos = None      
            
            self.update()
            self.fig.canvas.stop_event_loop()
                
        

    def on_move(self, event):
        if event.inaxes:
            ix_hover = self.get_ix(event)
            self.hoverpos = ix_hover
            self.fig.canvas.draw()

    def update(self):
        self.fig.canvas.draw()


class FDDPlotter:
    '''
    Class to initiate interactive peak picking plot for FDD (or directly on CPSD).    
    
    Example
    ----------
    Assuming we have a data matrix `data` sampled at `fs=128.0`.
    
    First, we import the necessary functions and classes:
            
        from koma.signal import xwelch
        from koma.oma import freq_svd
        from koma.plot import FDDPlotter
        
    Then, we need to create a CPSD matrix:
        
        f, cpsd = xwelch(data, fs=fs, nfft=2048, nperseg=1024)
        
    This is thereafter decomposed and used as input to create the FDDPlotter object:
        
        U, D = koma.oma.freq_svd(cpsd)
        fdd_plotter = koma.plot.FDDPlotter(D, U=U, f=f, lines=[0, 1], num=1)
        
    Finally, we get the figure:
        
        fig = fdd_plotter.get_fig()
        
    Scroll with the mouse to change the active line, click right mouse button to
    remove point and left mouse button to add point.
    
    
    '''
    def __init__(self, D, U=None, f=None, lines=None, colors=None, ax=None, 
                 num=None, active_settings={}, deactive_settings={}, picked_settings={},
                 vline_settings= {}, normalize=False, logaritmic=False, label_str='Singular line'):
            
        '''
        Parameters
        -------------
        D : float
            numpy array of dimensions N_dofs x N_dofs x N_freq; can either be diagonal
            (3d) array from FDD (i.e., SVD) or the CPSD matrix (frequencies along last
            axis) 
        U : float, optional
            if FDD is used as a prior step, this is the U matrix defining the orthogonal
            vectors for varying frequencies
        f : float, optional
            frequency axis corresponding to D (and U); if not given indices are used
        lines : int, optional
            list of integers corresponding to the indices of D; these define which
            diagonal terms to plot; if not defined, all components/lines are plotted
        colors : list, optional
            list of strings or rgb tuples to use for the different line plots;
            if not given, default color order is applied
        ax : `matplotlib.Axis` object
            axis to place plot in; if not given, a new axis in the figure specified
            will be created
        num : int, optional
            figure number used; only used if `ax` = None
        active_settings : dict
            dictionary with settings to use for lines that are active
        deactive_settings : dict
            dictionary with settings to use for lines that are inactive
        picked_settings : dict
            dictionary with settings to use for picked dot
        vline_settings : dict
            dictionary with settings to use for vertical line
        normalize : bool, default=False
            whether or not to normalize each line to its max value
        logaritmic : bool, default=False
            whether or not to use logaritmic y-axis            
        label_str : str, default='Singular line'
            label used in legend
            
        Returns
        ----------
        fdd_plotter : `FDDPlotter` object
        
        '''
        
        self.logaritmic = logaritmic
        self.normalize = normalize
        
        if U is None:
            U = np.stack([np.eye(D.shape[0])]*D.shape[2],axis=2)
            
        self._U = U
        self._D = D
        self.line_ixs = np.array(lines)
        self.index = 0
        self.lines = [None]*len(self.line_ixs)
        self.n_lines = len(self.line_ixs)
        self.active_settings = {'linewidth': 2.0, 'alpha':1.0} | active_settings
        self.deactive_settings = {'linewidth': 1.0, 'alpha': 0.5}  | deactive_settings
        self._picked = []
        self.picked_dots = []
        self.vline_settings = dict(ls='-', lw=0.5) | vline_settings
        self.title_format = 'ix = {index}'
        self.label_str = label_str

        if f is None:
            self.f = np.arange(self.D.shape[2])
        else:
            self.f = f

        if ax is None:
            plt.figure(num).clf()
            self.fig, self.ax = plt.subplots(num=num)
        else:
            plt.sca(ax)
            self.ax = ax
            self.fig = plt.gcf()

        if colors is None:
            self.colors = list(mcolors.TABLEAU_COLORS.values())
        else:
            self.colors = colors


    def get_fig(self, show=True, block=True):
        '''
        Get an interactive figure.
        
        Parameters
        -----------
        show : bool, default=True
            whether or not to automatically show figure
        block : bool, default=True
            whether or not to block before returning to retain interactivity
            
        Returns
        -----------
        fig : `matplotlib.Figure` object
            figure object
        '''
        
        for l in range(self.n_lines):
            if self.normalize:
                plotted = [self.D[l, l, :]/np.max(self.D[l, l, :]) for l in range(self.D.shape[0])]
            else:
                plotted = [self.D[l, l, :] for l in range(self.D.shape[0])]

            if self.logaritmic:
                plotted = [np.log10(plotted_l) for plotted_l in plotted]

            self.lines[l], = self.ax.plot(self.f, plotted[l], color=self.colors[l], linewidth=1.0, label=f'{self.label_str} {self.line_ixs[l]}')
        
        self.ax.legend(frameon=False)
        self.vline = self.ax.axvline(np.nan, **self.vline_settings)
        self.vline.set_color(self.current_line.get_color())

        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)   
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.fig.canvas.mpl_connect('key_press_event', self.on_keyboard)
        self.fig.canvas.mpl_connect('close_event', self.on_close)
        
        if block:
            self.title_format = 'ix = {index} (Press enter when selection is done)'

        self.update()   

        if show:
            plt.show()
            if block:
                self.fig.canvas.start_event_loop()

        return self.fig
    
    
    @property
    def picked(self):
        if len(self._picked) == 0:
            return self._picked
        else:
            _picked = np.array(self._picked)
            sort_ix = np.argsort(_picked[:,1])
            return _picked[sort_ix, :]
    
    @property
    def current_line(self):
        return self.lines[self.index]



    def on_click(self, event):
        if event.inaxes and event.button is MouseButton.LEFT and self.fig.canvas.manager.toolbar.mode.value == '':
            current_line = self.current_line
            color = current_line.get_color()
            ix_sel = np.argmin(np.abs(event.xdata - self.f))

            self._picked.append([self.index, ix_sel])
            self.picked_dots.append(self.ax.plot(self.f[ix_sel], current_line.get_ydata()[ix_sel], marker='o',
                            markersize=7, markerfacecolor=color, markeredgecolor='black')[0])

        elif event.inaxes and event.button is MouseButton.RIGHT:
            ix_sel = np.argmin(np.abs(event.xdata - self.f))
            picked_mat = np.array(self._picked)
            valid_ix = np.where((picked_mat[:,0] == self.index))[0]
            if len(valid_ix)>0:
                row_ix = np.argmin(np.abs(ix_sel-picked_mat[valid_ix, 1]))
                ix_remove = valid_ix[row_ix]

                self._picked.pop(ix_remove)
                self.picked_dots[ix_remove].remove()
                self.picked_dots.pop(ix_remove)

        self.fig.canvas.draw()

    def on_close(self, event):
        self.fig.canvas.stop_event_loop()

    def on_keyboard(self, event):
        if event.key == 'enter':
            self.title_format = 'ix = {index}'
            self.update()
            self.fig.canvas.stop_event_loop()
            
            # Prepare for printing
            self.vline.set_xdata([np.nan])
            for ix,line in enumerate(self.lines):
                line.set(**self.deactive_settings) 
                line.set(alpha=1.0)

    def on_move(self, event):
        if event.inaxes:
            ix_sel = np.argmin(np.abs(event.xdata-self.f))
            self.vline.set_xdata([self.f[ix_sel]])
            self.fig.canvas.draw()

    def on_scroll(self, event):
        increment = 1 if event.button == 'up' else -1
        self.index = np.clip(self.index + increment, 0, self.n_lines-1)
        self.vline.set_color(self.current_line.get_color())
        self.update()

    def update(self):
        self.ax.set_title(self.title_format.format(index = self.line_ixs[self.index]))
        for ix,line in enumerate(self.lines):
            if ix==self.index:
                line.set(**self.active_settings)
            else:
                line.set(**self.deactive_settings) 

        self.fig.canvas.draw()

    @property
    def D(self):
        return self._D[np.ix_(self.line_ixs, self.line_ixs, np.arange(self._D.shape[2]))]

    @property
    def U(self):
        return self._U[np.ix_(np.arange(self._U.shape[0]), self.line_ixs, np.arange(self._U.shape[2]))]

    @property
    def freq(self):
        if len(self.picked)>0:
            ixs = self.picked[:, 1]
            return self.f[ixs]

    @property
    def phi(self):
        if len(self.picked)>0:
            
            phi = np.zeros([self.U.shape[0], len(self._picked)]).astype(complex)
            for ix,pick in enumerate(self.picked):
                phi[:, ix] = self.U[:, pick[0], pick[1]]

            return phi



def stabplot(lambd, orders, phi=None, freq_range=None, frequency_unit='rad/s', damped_freq=False, psd_freq=None, psd_y=None, psd_plot_scale='log', 
    renderer=None, pole_settings=None, selected_pole_settings=None, to_clipboard='none', return_ix=False):
    """
    (DEPRECATED) Generate plotly-based stabilization plot from output from find_stable_poles.

    Arguments
    ---------------------------
    lambd : double
        array with complex-valued eigenvalues
    orders : int
        corresponding order for each pole in `lambd`
    phi : optional, double
        matrix where each column is complex-valued eigenvector corresponding to lambd
    freq_range : double, optional
        list of min and max values used for frequency axis
    frequency_unit : {'rad/s', 'Hz', 's'}, optional
        what frequency unit to use ('s' or 'period' enforces period rather than frequency) 
    damped_freq : False, optional
        whether or not to use damped frequency (or period) values in plot (False enforces undamped freqs)
    psd_freq : double, optional
        [not yet implemented] frequency values of plot to overlay, typically spectrum of data
    psd_y : double, optional
        [not yet implemented] function values of plot to overlay, typically spectrum of data
    psd_plot_scale: {'log', 'linear'}, optional
        how to plot the overlaid PSD (linear or logarithmic y-scale)
    renderer : None (render no plot - manually render output object), optional
        how to plot figure, refer plotly documentation for details 
        ('svg', 'browser', 'notebook', 'notebook_connected', are examples - 
        use 'default' to give default and None to avoid plot)
    to_clipboard : {'df', 'ix', 'none'}, optional
        update clipboard every time a pole is added, keeping selected indices or table
        'df' is not operational yet
    return_ix : False, optional
        whether or not to return second variable with indices - this is updated as more poles are selected
        

    Returns
    ---------------------------
    fig : obj
        plotly figure object

        
    Notes
    ----------------------------
    By hovering a point, the following data about the point will be given in tooltip:
        
         * Natural frequency / period in specified unit (damped or undamped)
         * Order 
         * Critical damping ratio in % (xi)
         * Index of pole (corresponding to inputs lambda_stab and order_stab)
    """

    warnings.warn("deprecated", DeprecationWarning)

    # Treat input settings
    if pole_settings is None: pole_settings = {}
    if selected_pole_settings is None: selected_pole_settings = {}

    unsel_settings = {'color':'#a9a9a9', 'size':6, 'opacity':0.6}
    unsel_settings.update(**pole_settings)

    current_settings = listify_each_dict_entry(unsel_settings, len(lambd))

    sel_settings = {'color':'#cd5c5c', 'size':10, 'opacity':1.0, 'line': {'color': '#000000', 'width': 0}}
    sel_settings.update(**selected_pole_settings)

    select_status = np.zeros(len(lambd), dtype=bool)
    
    # Create suffix and frequency value depending on whether damped freq. is requested or not
    if damped_freq:
        dampedornot = 'd'
        omega = np.abs(np.imag(lambd))
    else:
        dampedornot = 'n'
        omega = np.abs(lambd)

    # Create frequency/period axis and corresponding labels
    if frequency_unit == 'rad/s':
        x = omega
        xlabel = fr'$\omega_{dampedornot} [{frequency_unit}]$'
        tooltip_name = fr'\omega_{dampedornot}'
        frequency_unit = 'rad/s'
    elif frequency_unit.lower() == 'hz':
        x = omega/(2*np.pi)
        xlabel = fr'$f_{dampedornot} [{frequency_unit}]$'
        tooltip_name = fr'f_{dampedornot}'
        frequency_unit = 'Hz'
    elif (frequency_unit.lower() == 's') or (frequency_unit.lower() == 'period'):
        x = (2*np.pi)/omega
        xlabel = fr'Period, $T_{dampedornot} [{frequency_unit}]$'
        tooltip_name = fr'T_{dampedornot}'
        frequency_unit = 's'
    
    # Damping ratio and index to hover
    xi = -np.real(lambd)/np.abs(lambd)
    text = [f'xi = {xi_i*100:.2f}% <br> ix = {ix}' for ix, xi_i in enumerate(xi)]     # rewrite xi as %, and make string
    htemplate = f'{tooltip_name}' + ' = %{x:.3f} ' + f'{frequency_unit}<br>n =' + ' %{y}' +'<br> %{text}'

    # Construct dataframe and create scatter trace     
    poles = pd.DataFrame({'freq': x, 'order':orders})
    scatter_trace =  go.Scattergl(
                x=poles['freq'], y=poles['order'], mode='markers', name='',
                hovertemplate = htemplate, text=text,
                marker=current_settings)
    
    scatter_trace['name'] = 'Poles'

    # PSD overlay trace
    overlay_trace = go.Scatter(x=psd_freq, y=psd_y, mode='lines', name='PSD', 
                               hoverinfo='skip', line={'color':'#ffdab9'})
    
    # Create figure object, add traces and adjust labels and axes
    # fig = go.FigureWidget(scatter_trace)
    fig = make_subplots(rows=2, cols=1, specs=[[{"type": "xy", "secondary_y": True}],
           [{"type": "table"}]])
    fig.layout.hovermode = 'closest'
    fig.add_trace(scatter_trace, secondary_y=False, row=1, col=1)
    
    if psd_freq is not None:
        fig.add_trace(overlay_trace, secondary_y=True, row=1, col=1)
        fig.update_yaxes(title_text="PSD", secondary_y=True, type=psd_plot_scale)
        fig['layout']['yaxis2']['showgrid'] = False
    
    fig.layout['xaxis']['title'] = xlabel
    fig.layout['yaxis']['title'] = '$n$'

    if freq_range is not None:
        fig.layout['xaxis']['range'] = freq_range
        
    fig['layout']['yaxis']['range'] = [0.05, np.max(orders)*1.1]

    # Renderer (refer to plotly documentation for details)
    if renderer is 'default':
        import plotly.io as pio
        renderer = pio.renderers.default

    df = pd.DataFrame(columns=['ix','x','xi'])
    pd.options.display.float_format = '{:,.2f}'.format
    
    fig.add_trace(
        go.Table(header=
                     dict(values=['Pole index', xlabel, r'$\xi [\%]$'],
                          fill_color='paleturquoise',
                          align='left'),
                 cells=
                     dict(values=[],fill_color='lavender',
               align='left', format=["",".4",".4"])), row=2, col=1)
    
    fig.update_layout(
          height=1000,
          showlegend=False
          )

    fig = go.FigureWidget(fig)  #convert to widget
    
    # Callback function for selection poles
    ix = np.arange(0, len(lambd))
    ix_sel = []
            
    def update_table():
        df = pd.DataFrame(data={'ix': ix[select_status], 
                                'freq': x[select_status], 
                                'xi':100*xi[select_status]})
        
        df = df.sort_values(by=['freq'])
        
        if len(fig.data)==2:
            sel_ix = 1
        else:
            sel_ix = 2
        
        fig.data[sel_ix].cells.values=[df.ix, df.freq, df.xi]
    
    
    def toggle_pole_selection(trace, clicked_point, selector):
        def export_df():
            df.to_clipboard(index=False)
        
        def export_ix_list():
            import pyperclip #requires pyperclip
            ix_str = '[' + ', '.join(str(i) for i in ix_sel) + ']'
            pyperclip.copy(ix_str)

        for i in clicked_point.point_inds:
            if select_status[i]:
                for key in current_settings:
                    current_settings[key][i] = unsel_settings[key]
            else:
                 for key in current_settings:
                    current_settings[key][i] = sel_settings[key]       
                    
            select_status[i] = not select_status[i]     #swap status

            with fig.batch_update():
                trace.marker = current_settings
        
        update_table()
        ix_sel = ix[select_status]

        if to_clipboard == 'ix':
            export_ix_list()
        elif to_clipboard == 'df':
            export_df()

    fig.data[0].on_click(toggle_pole_selection)    

    if renderer == 'browser_legacy':
        from plotly.offline import plot
        plot(fig, include_mathjax='cdn')
    elif renderer is not None:
        fig.show(renderer=renderer, include_mathjax='cdn')
    
    if return_ix:
        return fig, ix_sel
    else:
        return fig


def listify_each_dict_entry(dict_in, n):
    dict_out = dict()
    for key in dict_in:
        dict_out[key] = [dict_in[key]]*n
        
    return dict_out


def plot_argand(phi, ax=None, colors=None, labels=None, **plot_settings):
    if ax is None:
        ax = plt.gca()
        
    if labels is None:
        labels = np.arange(len(phi))
    
    plot_settings = {'width': 0.0, 'linestyle': '-'} | plot_settings
    
    if type(colors) == str:
        colors = cm.get_cmap(colors, len(phi))
        if hasattr(colors, 'colors'):
            colors = colors.colors

    for ix,val in enumerate(phi):
        if colors is not None:
            color = colors[ix]
        else:
            color = None
        
        ax.arrow(0, 0, np.real(val), np.imag(val), color=color, label=labels[ix], **plot_settings)
    
    ymax = np.max(np.abs(ax.get_ylim()))
    xmax = np.max(np.abs(ax.get_xlim()))
    
    xymax = np.max([ymax,xmax])
    
    ax.set_ylim([-xymax, xymax])
    ax.set_xlim([-xymax, xymax])
    return ax

    
