MATLAB exaple: visualization of beam modes
--------------------------------------
See the file `./examples/matlab/beam_modes.m`. The content and results of the file are reproduced in the following subsection.

.. code :: matlab

    import ../../koma.*

    %% Geometry data for mode interpretation
    path='geometry_beam\';
    slave=dlmread([path 'slave.asc']);
    elements=dlmread([path 'elements.asc']);
    griddef=dlmread([path 'grid.asc']);
    active_nodes=2:8;   %corresponding to nodes with sensors, and thus has to match the dimensions of phi (/3)

    %% Plot mode
    mode_n = 1;

    x = griddef(:,2);
    phi_1d_all = cos(mode_n*(1-x/5) * pi/2);
    phi_1d_sensors = phi_1d_all(4:end);
    phi_3d = [phi_1d_sensors*0, phi_1d_sensors, phi_1d_sensors*0]';
    phi = phi_3d(:);

    koma.vis.plot_mode(phi, griddef, elements, slave, active_nodes, 'labels', 'deformed')

.. image:: matlab_beam_mode.png