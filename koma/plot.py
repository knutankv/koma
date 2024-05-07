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

class FDDPlotter:
    def __init__(self, U, D, f=None, lines=None, colors=None, ax=None, 
                 num=None, active_settings={}, deactive_settings={}, picked_settings={},
                 vline_settings= {}):
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

    def get_fig(self, show=True, block=True):
        for l in range(self.n_lines):
            self.lines[l], = self.ax.plot(self.f, self.D[l, l, :]/np.max(self.D[l, l, :]), color=self.colors[l], linewidth=1.0, label=f'Singular line {self.line_ixs[l]}')
        
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

    def on_click(self, event):
        if event.inaxes and event.button is MouseButton.LEFT and self.fig.canvas.manager.toolbar.mode.value == '':
            current_line = self.current_line
            color = current_line.get_color()
            ix_sel = np.argmin(np.abs(event.xdata - self.f))

            self._picked.append([self.index, ix_sel])
            self.picked_dots.append(self.ax.plot(self.f[ix_sel], current_line.get_ydata()[ix_sel], marker='o',
                            markersize=7, markerfacecolor=color, markeredgecolor='red')[0])

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

    def on_move(self, event):
        if event.inaxes:
            ix_sel = np.argmin(np.abs(event.xdata-self.f))
            self.vline.set_xdata(self.f[ix_sel])
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


def stabplot(lambd, orders, phi=None, model=None, freq_range=None, frequency_unit='rad/s', damped_freq=False, psd_freq=None, psd_y=None, psd_plot_scale='log', 
    renderer=None, pole_settings=None, selected_pole_settings=None, to_clipboard='none', return_ix=False):
    """
    Generate plotly-based stabilization plot from output from find_stable_poles. This is still beta!

    Arguments
    ---------------------------
    lambd : double
        array with complex-valued eigenvalues
    orders : int
        corresponding order for each pole in `lambd`
    phi : optional, double
        matrix where each column is complex-valued eigenvector corresponding to lambd
    model : optional, double
        model object which is required input for plotting phi (based on geometry definition of system, constructed by `Model` class)
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
        xlabel = f'$\omega_{dampedornot} \; [{frequency_unit}]$'
        tooltip_name = f'\omega_{dampedornot}'
        frequency_unit = 'rad/s'
    elif frequency_unit.lower() == 'hz':
        x = omega/(2*np.pi)
        xlabel = f'$f_{dampedornot} \; [{frequency_unit}]$'
        tooltip_name = f'f_{dampedornot}'
        frequency_unit = 'Hz'
    elif (frequency_unit.lower() == 's') or (frequency_unit.lower() == 'period'):
        x = (2*np.pi)/omega
        xlabel = f'Period, $T_{dampedornot} \; [{frequency_unit}]$'
        tooltip_name = f'T_{dampedornot}'
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
               align='left')), row=2, col=1)
    
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
        
        df = df.sort('freq')
        
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


def plot_argand(phi, ax=None, colors=None, **plot_settings):
    if ax is None:
        ax = plt.gca()
    
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
            
        ax.arrow(0, 0, np.real(val), np.imag(val), color=color, label=ix, **plot_settings)
    
    ymax = np.max(np.abs(ax.get_ylim()))
    xmax = np.max(np.abs(ax.get_xlim()))
    
    ax.set_ylim([-ymax, ymax])
    ax.set_xlim([-xmax, xmax])
    return ax

    
