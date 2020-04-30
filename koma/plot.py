"""
##############################################
Visualization module
##############################################
All functions related to the plotting of results from OMA methods.
"""

import plotly.express as px
from plotly.offline import plot as plot_browser   
import pandas as pd
import numpy as np

def stabplot(lambda_stab, order_stab, frequency_unit='rad/s', damped_freq=False, overlay_x=None, overlay_y=None):
    """
    Plot stabilization plot from lambda and phi arrays (output from covssi). Based on plotly => web browser is used.

    Arguments
    ---------------------------
    lambd_stab : double
        array with complex-valued eigenvalues deemed stable
    phi_stab : double
        2d array with complex-valued eigenvectors (stacked as columns), each column corresponds to a mode
    frequency_unit : {'rad/s', 'Hz', 's'}, optional
        what frequency unit to use ('s' or 'period' enforces period rather than frequency) 
    damped_freq : False, optional
        whether or not to use damped frequency (or period) values in plot (False enforces undamped freqs)
    overlay_x : double, optional
        [not yet implemented] frequency values of plot to overlay, typically spectrum of data
    overlay_y : double, optional
        [not yet implemented] function values of plot to overlay, typically spectrum of data
        
        
    Notes
    ----------------------------
    By hovering a point, the following data about the point will be given in tooltip:
        
         * Natural frequency / period in specified unit (damped or undamped)
         * Order 
         * Critical damping ratio in % (xi)
         * Index of pole (corresponding to inputs lambda_stab and order_stab)
    """
    
    if damped_freq:
        add = ' (damped)'
        omega = np.abs(np.imag(lambda_stab))
    else:
        add = ' (undamped)'
        omega = np.abs(lambda_stab)
        
    ix = np.arange(0, len(lambda_stab))    # index corresponding to 

    if frequency_unit == 'rad/s':
        x = omega
        xlabel = 'Frequency [rad/s]' + add
    elif frequency_unit.lower() == 'hz':
        x = omega/(2*np.pi)
        xlabel = 'Frequency [Hz]' + add
    elif (frequency_unit.lower() == 's') or (frequency_unit.lower() == 'period'):
        x = (2*np.pi)/omega
        xlabel = 'Period [s]' + add

    xi_stab = -np.real(lambda_stab)/np.abs(lambda_stab)

    # Rewrite xi as %, and make string
    xi_str = [f'{xi_i*100:.2f}%' for xi_i in xi_stab]
    
    stable_poles = pd.DataFrame({xlabel: x, 'xi':xi_str, 'Order':order_stab, 'ix':ix})
    fig = px.scatter(stable_poles, x=xlabel, y='Order', hover_data=['xi', 'ix'])
    fig.update_layout(title='Interactive stability plot')
    plot_browser(fig)

