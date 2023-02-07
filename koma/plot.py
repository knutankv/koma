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
    overlay_trace = go.Scatter(x=psd_freq, y=psd_y, mode='lines', name='PSD', hoverinfo='skip', line={'color':'#ffdab9'})
    
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
        df = pd.DataFrame(data={'ix': ix[select_status], 'freq': x[select_status], 'xi':100*xi[select_status]})
        
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