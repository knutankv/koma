function [axis_handle, scaling, cm] = plot_mode(phi, griddef, elements, slave, active_nodes, varargin)
%% Plot translational part of mode, by inputting translational mode shape.
% 
% Arguments
% ---------------------------
% phi : double
%   complex-valued mode shape to plot (translational part only)
% griddef : double
%    matrix defining the nodes of the model (see Notes for details)
% elements : int
%   matrix defining the elements of the model, connecting the node numbers
%   together as patches or beams (see Notes for details)
% slave : int
%   matrix defining linear relations between DOFs (see Notes for details)
% active_nodes : int
%   Index of nodes that are active (highlighted). If empty, all nodes are assumed active!
% scaling : [], optional
%   numerical value defining linear scaling of displacement vector
%   ([] results in automatic scaling).
% scalefactor : 1, optional
%   second factor used to adjust automatic scaling
% axishandle : [], optional
%   axis handle to make plot in
% cam : [], optional
%   camera properties 
% cmdir : 'abs', optional
%   component to base the colormap on ('abs', 'x', 'y' or 'z')
% colormap : [], optional
%    predefined colormap vector (max and min val.), [] leads to auto.
% plots : 'both', optional
%   which deformation state to plot ('both, 'undeformed', or 'deformed')
% labels : 'none', optional
%   put labels on the figure (which state)? ('both', 'undeformed',
%   'deformed', or 'none')
% title : '', optional
%   title of figure, at position as defined in `cam`
% axis : false, optional          
%   whether or not to show axis
% edge : 'both', optional         
%   show edge color for surfaces ('both', 'undeformed', 'deformed', 'none')
% nodes : 'both', optional        
%   switch for nodes ('both','active', or'inactive')
% lines : true, optional         
%   switch for lines showing displacements 
% renderer : 'zbuffer', optional     
%   figure renderer ('zbuffer', 'painters', or 'opengl')
% maxreal : true, optional         
%   make the real part as big as possible or not          
% linewidth : 1, optional
%   line width on edges and lines 
%
% Returns
% -----------------------
% axis_handle : obj
%   handle of axis object
% scaling : double
%   used scaling
% cm : obj
%   used colormap
%
% Notes 
% ---------------------------
% `griddef`, `elements` and `slave` follow the exact same definition as for 
% the MACEC toolbox. See details below.
% 
% **griddef** is exemplified as follows for arbitrary three nodes (exclude header):  
%
%  ============   ===== ===== ======
%                   Coordinates          
%  ------------   ------------------
%   Node no.       X     Y      Z    
%  ============   ===== ===== ======
%   1              0      0      1    
%   2              0.1   0.2    0.3    
%   3              -2    0.3    10    
%  ============   ===== ===== ======
%  
% 3D input is required even if the analysis is 2D. Other than that, its definition
% is assumed self-explanatory.
%
% **elements** is exemplified as follows (exclude header):
%
% =======  ========  =========
% Node 1    Node 2    Node 3
% =======  ========  =========
%   1         2          3
%   2         3          0
% =======  ========  =========
%
% This example creates two elements: (1) a patch between nodes 1,2, and 3
% and (2) a beam element between nodes 2 and 3 (node 3 = 0 => beam element).
% If no patch elements are needed, the table's first two columns are
% sufficient.
%
% **slave** is exemplified as follows (exclude header):  
%
% ============   === ==== ==== ============ ======== ========= ========
%                  Master DOF                Slave DOFs' amplitudes          
% ------------   ------------- ------------ ---------------------------
%  Master         X    Y    Z    Slave        X         Y         Z   
% ============   === ==== ==== ============ ======== ========= ========
%      1          1    0    0        2         0.3      0.7       0  
%      1          0    1    0        4         0         1        0
%      3          1    0    0        5         1         0        0  
% ============   === ==== ==== ============ ======== ========= ========
%
% The displacement of the slave node is defined based on the displacement of
% the master node in either X, Y or Z direction. A unity value indicates that 
% this is the local DOF to use as master (keep the others equal 0). The 
% amplitudes of the different DOFs of the slave nodes are defined by the 
% values in the latter three columns.
%
% * The example table above describes the following slave-master definition:
%
%       - Node  1, DOF X, dictates the motion of node 2 by a 
%         factor 0.3 on the X-component and a factor of 0.7 on the Y-component.
%       - Node 1, DOF Y, dictates the motion of node 4 by a factor of 1 on the
%         Y-component.
%       - Node 3, DOF X, dictates the motion of node 5 by a factor of 1 on the
%         X-component.
%

%% INPUT HANDLING
p=inputParser;
p.KeepUnmatched=true;

addParameter(p, 'scaling', [],@isnumeric);
addParameter(p, 'scalefactor', 1,@isnumeric)
addParameter(p, 'axishandle', []);
addParameter(p, 'cam', [])
addParameter(p, 'cmdir', 'abs', @ischar);
addParameter(p, 'colormap', []);
addParameter(p, 'plots', 'both');
addParameter(p, 'labels', 'none', @ischar);
addParameter(p, 'title', '');
addParameter(p, 'axis', false)
addParameter(p, 'edge', 'both')
addParameter(p, 'lines', 'active')
addParameter(p, 'nodes', 'both', @ischar)
addParameter(p, 'renderer', 'zbuffer', @ischar)
addParameter(p, 'maxreal', true)
addParameter(p, 'dims', 3)
addParameter(p, 'linewidth', 1)

parse(p,varargin{:})
scaling = p.Results.scaling;
axishandle = p.Results.axishandle;
cmdir = p.Results.cmdir;
labels = p.Results.labels;
cm = p.Results.colormap;
ax = p.Results.axis;
plots = p.Results.plots;
edg = p.Results.edge;
nodes = p.Results.nodes;
cam = p.Results.cam;
renderer = p.Results.renderer;
scalefactor = p.Results.scalefactor;
tit = p.Results.title;
lines = p.Results.lines;
realPart = p.Results.maxreal;
dims = p.Results.dims;
lw = p.Results.linewidth;

if dims ~=3
    phitemp = phi;
    phi = [];
    for i = 1:length(phi)
       phi = [phi;phitemp(i);zeros(dims,1)];
    end
end

if strcmp(edg,'both')
    edge{1} = 'black';
    edge{2} = 'black';
elseif strcmp(edg,'deformed')
    edge{1} = 'none';
    edge{2}  ='black';
elseif strcmp(edg,'undeformed')
    edge{1} = 'black';
    edge{2} = 'none';
else
    edge{1} = 'none'; 
    edge{2} = 'none';
end

if ~isempty(axishandle)
    axes(axishandle);
end

hold on;
oldlines=findobj(gcf,'Type', 'line'); oldpatches = findobj(gcf,'Type', 'patch');

if isempty(phi)
    active_nodes = zeros(size(griddef,1),1);
    phi = zeros(3*size(griddef,1),1);
end

if isempty(active_nodes)
   active_nodes = griddef(:,1);
end

% Check if input are OK

if 3*length(active_nodes)~=size(phi,1)
        error('Dimensions of phi and active_nodes have to match!');
end

if isempty(elements)
    for i = 1:size(griddef,1)-2
        elements(i,:) = i:i+2;
    end
end

%Adjust element list so that the numbers refer to the index in the node
%list, and not the node labels
elementsLabels=elements;
clear elements

for i = 1:size(elementsLabels,1)
    elements{i} = find(ismember(griddef(:,1)', elementsLabels(i,:)));
end

%% INITIAL
%ESTABLISH XYZ-GRID
XYZ0 = griddef(:,2:end);
active_nodes = find(ismember(griddef(:,1),active_nodes));   %find indices in griddef corresponding to active_nodes

%PLOT GRID (NODES)
Npoints = size(griddef,1);

%% MAXIMIZE THE REAL PART
if realPart==true
    angleVector = 0:0.01:pi/2;
    rot_mode = phi*exp(angleVector*1i);
    [~,iReMax] = max(sum(real(rot_mode).^2,1));  %find index with the angle giving the largest real part
    phi = phi*exp(angleVector(iReMax)*1i);
    
    phi = real(phi);
end

%% CALCULATE DISPLACEMENTS AT EACH POINT
% ENFORCE INPUT DISPLACEMENTS
Nnodes = size(XYZ0,1);    %number of points/nodes
dXYZ = zeros(3,Nnodes);
for p = 1:length(active_nodes)
    p0 = (p-1)*3;
    dXYZ(:,active_nodes(p)) = phi(p0+1:p0+3);
end

% ENFORCE LINEAR SLAVE DISPLACEMENTS (nothing happens if slave is empty!)
for s = 1:size(slave,1)
    masternode = slave(s,1);
    slavenode = slave(s,5);
    masterdof = slave(s,2:4)~=0;
    amp = slave(s,6:8)';
    dXYZ(1:3,slavenode) = dXYZ(1:3,slavenode)+amp.*dXYZ(masterdof, masternode);
end

dX = dXYZ(1,:);
dY = dXYZ(2,:);
dZ = dXYZ(3,:);

dXYZtot = sqrt(dX.^2+dY.^2+dZ.^2);

%CREATE COLOR MAP FOR AMPLITUDES
if strfind(cmdir, 'x')
    d = abs(dX);
elseif strfind(cmdir, 'y')
    d = abs(dY);
elseif strfind(cmdir, 'z')
    d = abs(dZ);
elseif strfind(cmdir, 'abs')
    d = dXYZtot;
end

%% PLOT DEFORMED FIGURE
if strcmp(plots,'deformed') || strcmp(plots,'both')
    %SCALING
    if isempty(scaling)
        %Auto scaling algorithm
        k = 0.07;
        scaling = scalefactor.*k.*norm([max(XYZ0(:,1)) - min(XYZ0(:,1)), max(XYZ0(:,2))-min(XYZ0(:,2)), max(XYZ0(:,3))-min(XYZ0(:,3))])/norm([max(abs(dX)), max(abs(dY)), max(abs(dZ))]);
    end
    
    dX = dX*scaling;
    dY = dY*scaling;
    dZ = dZ*scaling;
    d = d*scaling;
    
    XYZ = zeros(length(active_nodes),3);
    
    %PLOT DEFORMED GRID
    
    for p=1:Npoints
        if ismember(griddef(p,1),active_nodes)
            col='.r';
        else
            col='.k';
        end
        XYZ(p,:) = XYZ0(p,:) + [dX(p),dY(p),dZ(p)];
        if strcmp(nodes,'both')
            plot3(XYZ(p, 1), XYZ(p, 2), XYZ(p, 3), col)
        elseif strcmp(nodes,'active') && strcmp(col,'.r')
            plot3(XYZ(p, 1), XYZ(p, 2), XYZ(p, 3), col)
        end
        if strcmp(labels,'deformed') || strcmp(labels,'both')
            text(XYZ(p,1), XYZ(p,2), XYZ(p,3), num2str(griddef(p,1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize',8)
        end
        
        if lines==true
            h = quiver3(XYZ0(p,1),XYZ0(p,2),XYZ0(p,3),dX(p),dY(p),dZ(p),0);
            set(h, 'LineStyle', '-', 'Color', 'k');
        end
    end
    
    % PLOT BEAMS AND PATCHES

    for el=1:length(elements)
        if length(elements{el}) == 2 %BEAMS
            line(XYZ(elements{el},1),XYZ(elements{el},2),XYZ(elements{el},3),'color', 'red', 'LineWidth', lw);
        elseif length(elements{el}) ==3 %3-NODED PATCHES
            fill3(XYZ(elements{el},1)', XYZ(elements{el},2)', XYZ(elements{el},3)', d(elements{el})', 'EdgeColor', edge{2},'LineWidth', lw);
        else     
            warning('Elements have to be made out of 2 or 3 three nodes. Element avoided!')
        end
    end
end

if isempty(cm)
    cm = [min(d),max(d)];
end

%% PLOT UNDEFORMED FIGURE
if strcmp(plots,'undeformed') || strcmp(plots,'both')
    for p=1:Npoints
        if ismember(griddef(p,1),active_nodes)
            col='.r';
        else
            col='.k';
        end
        if strcmp(nodes,'both')
            plot3(XYZ0(p,1),XYZ0(p,2),XYZ0(p,3),col)
        elseif strcmp(nodes,'active') && strcmp(col,'.r')
            plot3(XYZ0(p,1),XYZ0(p,2),XYZ0(p,3),col)
        end
        if strcmp(labels,'undeformed') || strcmp(labels,'both')
            text(XYZ0(p,1),XYZ0(p,2),XYZ0(p,3), num2str(griddef(p,1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8)
        end
    end
    
    
    for el=1:length(elements)
        if length(elements{el}) == 3  %3-NODED PATCHES
            patch(XYZ0(elements{el},1),XYZ0(elements{el},2),XYZ0(elements{el},3),0*XYZ0(elements{el},3)' ,'FaceColor', 'green','EdgeColor',edge{1},'LineWidth',lw);
        elseif length(elements{el}) == 2 %BEAMS
            line(XYZ0(elements{el},1),XYZ0(elements{el},2),XYZ0(elements{el},3),'color', 'blue','LineWidth',lw);
        else
            warning('Elements have to be made out of 2 or 3 three nodes. Element avoided!')
        end
    end
end



%% FINAL OPERATIONS
axis equal
delete(oldlines); delete(oldpatches);

%TITLE
th=text(1,1,1,tit);

if ~isempty(cam)
    %SET AXIS PROPERTIES
    set(gca,'CameraPositionMode','manual')
    set(gca,'CameraTargetMode','manual')
    set(gca,'CameraUpVectorMode','manual')
    set(gca,'CameraViewAngleMode','manual')
    set(gca,'PlotBoxAspectRatioMode','manual')
% Moved inside if ^^ 


    set(gca,'CameraPosition',cam.CameraPosition)
    set(gca,'CameraTarget',cam.CameraTarget)
    set(gca,'CameraUpVector',cam.CameraUpVector)
    set(gca,'CameraViewAngle',cam.CameraViewAngle)
    set(gca,'Position',cam.Position)
    set(gca,'Projection',cam.Projection)
    set(gcf,'Position',cam.FigurePosition)
    set(th,'Position',cam.TitlePosition);
end

if ~strcmp(plots,'undeformed')  %unless only undeformed is wanted
    if isnan(cm)
       cm = [-1 1]; 
    end
    set(gca,'CLim',cm)
end

%SET RENDERER
set(gcf,'renderer',renderer)

%SET AXIS SETTINGS
if ax==false
    set(gca,'visible','off');
    set(gcf,'Color','white');
elseif ax==true
    set(gca,'visible','on');
    set(gcf,'Color',[0.8000    0.8000    0.8000]);
end

hold off
axis_handle = gca;
colormap('parula')
