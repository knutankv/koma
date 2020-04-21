function animate_mode(phi, w, griddef, elements, slave, active_nodes, varargin)
%% Plot translational part of mode, by inputting translational mode shape.
% 
% Arguments
% ---------------------------
% phi : double
%   complex-valued mode shape to plot (translational part only)
% w : double
%   scalar frequency value (rad/s) in which to animate the mode vibration
% griddef : double
%    matrix defining the nodes of the model (see Notes for details)
% elements : int
%   matrix defining the elements of the model, connecting the node numbers
%   together as patches or beams (see Notes for details)
% slave : int
%   matrix defining linear relations between DOFs (see Notes for details)
% active_nodes : int
%   Index of nodes that are active (highlighted). If empty, all nodes are assumed active!
% movietitle : str, optional
%   title (name) for movie file (if not given a dialog will appear)
% fps : 30, optional    
%   frames per second of generated video file
% laps : 1, optional
%   laps to animate the mode
% rotations : 0, optional
%    rotations (pan-style) about the z-axis during animation
% xi : 0, optional
%   damping estimate for mode, used for animation
% cam : 0, optional
%   as for plot_mode, camera defintions
% complexmode : true, optional   
%   wheter or not to plot the complex mode
%
% Notes
% ------------
% All of the optional inputs to plot_mode are valid as input as well.
% Please refer to doc of `plot_mode`.


%% INPUT HANDLING
p=inputParser;
p.KeepUnmatched=true;
addParameter(p,'movietitle',[]);
addParameter(p,'fps',30)
addParameter(p,'laps',1)
addParameter(p,'rotations',0)
addParameter(p,'xi',0)
addParameter(p,'cam',[]);
addParameter(p,'complexmode',true)
parse(p,varargin{:})

laps=p.Results.laps;
fps=p.Results.fps;
rotations=p.Results.rotations;
xi=p.Results.xi;
movietitle=p.Results.movietitle;
cam=p.Results.cam;
complexmode=p.Results.complexmode;

%% FILE NAME
if isempty(movietitle)
    [f,p]=uiputfile('*.mp4','Save movie file as...');
    movietitle=[p,f];
end

%% REAL MODE OR COMPLEX MODE?
if complexmode==false
        angleVector=0:0.01:pi/2;
        rot_mode=phi*exp(angleVector*1i);
        [~,iReMax] = max(sum(real(rot_mode).^2,1));  %find index with the angle giving the largest real part
        phi=phi*exp(angleVector(iReMax)*1i);
    phi=real(phi);
end

%% PLOT MODES
%Find colormap
[~,scaling,CM]=koma.vis.plot_mode(abs(phi),griddef,elements,slave,active_nodes,varargin{:},'colormap',[]);    %Plot absolute value of mode, to establish the maximum values for colormap
pause
%Calculations...
f=w/(2*pi);     %frequency in Hz
dt=1/fps;       %seconds
tmax=laps/f;    %seconds
t=0:dt:tmax;    %time axis
rot=linspace(0,360*rotations,length(t));        %rotation vector
wd=(1-xi^2).^0.5*w; %damped natural frequency

%If no figure/camera information is given - construct one!
if isempty(cam)
    cam=getCam;
end

%Video
M=VideoWriter(movietitle,'MPEG-4');   
M.FrameRate = 30;
open(M);

%Animation, run plotMode for each discrete time instance
for n=1:length(t)
    tnow=t(n);
    u=real(exp(-xi.*w.*tnow) .* exp(wd*1i*tnow).*phi);
    clf
    koma.vis.plot_mode(u,griddef,elements,slave,ap,varargin{:},'cam',cam,'colormap',CM,'scaling',scaling);
    
    if rotations~=0  %ROTATIONS
        camorbit(rot(n),0,'data',[0 0 1])
    end    
    
    frame=getframe(gcf);
    writeVideo(M,frame); %write to video file
end

close(M);
disp('***** FILE SAVED *****')
end