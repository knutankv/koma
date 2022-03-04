MATLAB example: shear frame
--------------------------------------
See the file `./examples/matlab/covssi_example.m`. The content and results of the file are reproduced in the following subsection.

.. code :: matlab

    import koma.*

    %% Definitions
    K = [8000, -8000,0; -8000, 16000, -8000; 0, -8000, 16000];
    C = [10, 0,0; 0, 10, 0; 0, 0, 10];
    M = [500,0,0; 0, 500,0; 0, 0, 500];

    lambda_ref = polyeig(M,C,K);
    data = csvread('response_data.csv');
    levels = 3;
    fs = 3.0;

    t = 0:1/fs:(1/fs)*(size(data, 1)-1);

    xlabel('Time [s]')

.. image:: matlab_example_response.png

.. code :: matlab

    %% Response plot, mid-level
    figure(100),clf

    plot(t, data(:,2))
    hold on
    xlim([0,100])
    legend({'+ 100% noise' 'Clean signal'})

.. image:: matlab_example_noise.png

.. code :: matlab

%% Geometry data for mode interpretation
    path='geometry\';
    slave=dlmread([path 'slave.asc']);
    elements=dlmread([path 'elements.asc']);
    griddef=dlmread([path 'grid.asc']);
    active_nodes=[4,6,8];   %corresponding to nodes in grid-matrix; the nodes that are represented in the mode shape vector
  
.. image:: matlab_example_psd.png

.. code :: matlab

    %% SSI run
    i = 20;    %blockrows
    order = 2:2:50;
    s = 2; % stabilization level
    stabcrit = [0.05, 0.1, 0.1];
    indicator = 'freq';
    figure(5)

    [lambda,phi,order] = koma.oma.covssi(data, fs, i, 'order',order);
    koma.vis.stabplot(lambda,phi,order,'plot','stable','indicator','freq','stablevel',s,'selection',true,'grid',griddef,'slave',slave,'elements',elements,'active_nodes',active_nodes, 'convert_to_hz', false)

Resulting output:

.. code ::

    *** COVARIANCE DRIVEN SSI ALGORITHM FOR OPERATIONAL MODAL ANALYSIS ****
    * ESTABLISHING HANKEL AND TOEPLITZ MATRICES
    ** Correlation estimation
    ** Matrix stacking
    * CALCULATING SVD OF MATRICES AND CONTROLLING MAXIMUM ORDER
    ** Rank(D) = 60
    ** Maximum system order used is n_max = 50
    * ESTABLISHING WEIGHTING MATRIX
    * COMPUTING STATE MATRIX FROM DECOMPOSITION FOR EACH ORDER, AND ESTABLISH EIGENVALUES AND EIGENVECTORS
    * COMPUTATION COMPLETE!

.. image:: matlab_example_stab.png



By selecting a pole in the stabilization plot, information about the mode is provided in the tooltip. The keyboard can be used for certain commands:

 * 'P': plot the mode shape in a new figure (figure 999)
 * 'A': append mode to results
 * 'C': clear all appended modes
 * 'S': save the appended modes to specified mat-file
 * 'T': save the appended modes to 'temp_modes.mat', and clear appended data

This is how it looks when the selected pole is directly plotted: 

.. image:: matlab_stabplot.jpg


Finally, the poles can be automatically selected using pick_stable_modes, as follows:

.. code ::

    %% Pick stable modes
    slack = [0.1, 0.1, 0.1];

    [lambda_stab, phi_stab, order_stab, idx_stab] = koma.modal.find_stable_poles(lambda, phi, order, s, stabcrit, 'freq');
    [lambda_picked,phi_picked,stats] = koma.modal.pick_stable_modes(lambda_stab, phi_stab, slack);
    disp(abs(lambda_picked))

This prints the following natural frequencies:

.. code ::

    1.7780
    4.9833
    7.1971

