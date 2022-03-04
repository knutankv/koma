import ../../koma.*
close all
clear all

%% Definitions
K = [8000, -8000,0; -8000, 16000, -8000; 0, -8000, 16000];
C = [10, 0,0; 0, 10, 0; 0, 0, 10];
M = [500,0,0; 0, 500,0; 0, 0, 500];

lambda_ref = polyeig(M,C,K);
data = csvread('response_data.csv');
levels = 3;
fs = 3.0;

t = 0:1/fs:(1/fs)*(size(data, 1)-1);

%% Response plot, mid-level
figure(100),clf

plot(t, data(:,2))
hold on
xlim([0,100])
legend({'+ 100% noise' 'Clean signal'})

%% Geometry data for mode interpretation
path='geometry\';
slave=dlmread([path 'slave.asc']);
elements=dlmread([path 'elements.asc']);
griddef=dlmread([path 'grid.asc']);
active_nodes=[4,6,8];   %corresponding to nodes in grid-matrix; the nodes that are represented in the mode shape vector

%% Plot spectral density
figure(3), clf

% semilogy(f*2*pi,S(:,2)/2/pi)
ylabel('log(S)')
xlabel('\omega [rad/s]')
hold on
%% SSI run
i = 20;    %blockrows
order = 2:2:50;
s = 2; % stabilization level
stabcrit = [0.05, 0.1, 0.1];
indicator = 'freq';
figure(5)

[lambda,phi,order] = koma.oma.covssi(data, fs, i, 'order',order);

koma.vis.stabplot(lambda,phi,order,'plot','stable','indicator','freq','stablevel',s,'selection',true,'grid',griddef,'slave',slave,'elements',elements,'active_nodes',active_nodes, 'convert_to_hz', false)

%% Pick stable modes
slack = [0.1, 0.1, 0.1];

[lambda_stab, phi_stab, order_stab, idx_stab] = koma.modal.find_stable_poles(lambda, phi, order, s, stabcrit, 'freq');
[lambda_picked,phi_picked,stats] = koma.modal.pick_stable_modes(lambda_stab, phi_stab, slack);
disp(abs(lambda_picked))
