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

def stabplot(lambda_stab, orders_stab, freq_range=None, frequency_unit='rad/s', damped_freq=False, psd_freq=None, psd_y=None, psd_plot_scale='log', renderer='browser_legacy'):
    """
    Generate plotly-based stabilization plot from output from find_stable_poles.

    Arguments
    ---------------------------
    lambd_stab : double
        array with complex-valued eigenvalues deemed stable
    orders_stab : int
        corresponding order for each stable mode
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
    renderer : 'browser_legacy', optional
        how to plot figure, refer plotly documentation for details ('svg', 'browser', 'notebook', 'notebook_connected', are examples)
        

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
    
    if damped_freq:
        dampedornot = 'd'
        omega = np.abs(np.imag(lambda_stab))
    else:
        dampedornot = 'n'
        omega = np.abs(lambda_stab)
        
    ix = np.arange(0, len(lambda_stab))    # index corresponding to 

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

    # lambda_stab[np.abs(lambda_stab) == 0] = np.nan
    xi_stab = -np.real(lambda_stab)/np.abs(lambda_stab)

    # Rewrite xi as %, and make string
    text = [f'xi = {xi_i*100:.2f}% <br> ix = {ix}' for ix, xi_i in enumerate(xi_stab)]

    htemplate = f'{tooltip_name}' + ' = %{x:.2f} ' + f'{frequency_unit}<br>n =' + ' %{y}' +'<br> %{text}'
    
    stable_poles = pd.DataFrame({'freq': x, 'order':orders_stab})
    scatter_trace =  go.Scatter(
                x=stable_poles['freq'], y=stable_poles['order'], mode='markers', name='',
                hovertemplate = htemplate, text=text,
                marker={'color':'#4682b4'}
                )

    overlay_trace =  go.Scatter(x=psd_freq, y=psd_y, mode='lines', name='PSD', hoverinfo='skip', line={'color':'#cd5c5c'})
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    scatter_trace['name'] = 'Poles'
    fig.add_trace(scatter_trace, secondary_y=False)
    
    if psd_freq is not None:
        fig.add_trace(overlay_trace, secondary_y=True)
        fig.update_yaxes(title_text="PSD", secondary_y=True, type=psd_plot_scale)
        fig['layout']['yaxis2']['showgrid'] = False
        
    fig.layout['xaxis']['title'] = xlabel
    fig.layout['yaxis']['title'] = '$n$'

    if freq_range is not None:
        fig.layout['xaxis']['range'] = freq_range
        
    fig['layout']['yaxis']['range'] = [0.05, np.max(orders_stab)*1.1]

    if renderer is None:
        import plotly.io as pio
        renderer = pio.renderers.default
        
    if renderer == 'browser_legacy':
        from plotly.offline import plot
        plot(fig, include_mathjax='cdn')
    else:
        fig.show(renderer=renderer, include_mathjax='cdn')

    return fig
