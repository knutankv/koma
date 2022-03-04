function [f] = stabplot(lambda, phi, orders, varargin)
%% Plot stabilization plot from lambda and phi arrays (output from covssi).
% 
% Arguments
% ---------------------------
% lambda : double
%     cell array of arrays with complex-valued eigenvalues, one cell entry per order
% phi : double
%     cell array of matrices with complex-valued eigevectors (stacked column-wise), one cell entry per order
%   matrix)
% orders : int
%     1d array with orders corresponding to the cell array entries in lambd and phi
% s : int, optional
%     stability level, see :cite:`Kvale2017_OMA`
% stabcrit : {'freq': 0.05, 'damping':0.1, 'mac': 0.1}, optional
%     criteria to be fulfilled for pole to be deemed stable
% indicator : 'freq', optional
%     what modal indicator to use ('freq' or 'mac')
% plots : 'both', optional
%     what poles to plot ('both' or 'stable')
% markers : '.', optional     
%     first value applies to all poles, second to stable poles
% markersize : [36], optional  
%     first value applies to all poles, second to stable poles
% colors : {'r' 'k'}, optional
%     first value applies to all poles, second to stable poles
% picked_lambda : [], optional
%     plot vertical lines of these lambda values
% picked_error : [], optional  
%     plot corresponding error bar (frequency, abs(lambda), and not lambda)
%     relying on thirdparty contribution by Jos van der Geest, `herrorbar`
% figure : [], optional
%     specified figure handle
% selection : false, optional
%     whether or not to enable selection of poles for plotting and output
% grid : [], optional
%     if selection is true, this is used to plot mode shape
% elements : [], optional
%     if selection is true, this is used to plot mode shape
% slave : [], optional
%     if selection is true, this is used to plot mode shape
% active_nodes : [], optional
%     if selection is true, this is used to plot mode shape
% convert_to_hz: false, optional
%     whether the plot should show frequencies in Hz (true) or rad/s
%     (false, default)
%
%
% References
% -----------------
% Jos (10584) (2020). HERRORBAR (https://www.mathworks.com/matlabcentral/fileexchange/3963-herrorbar), MATLAB Central File Exchange. Retrieved April 21, 2020.


global osel
global msel
global lambda_sel
global phi_sel

osel = 0;
msel = 0;
lambda_sel = [];
phi_sel = [];

p=inputParser;
addParameter(p,'stabcrit',[0.01,0.05,0.05]);
addParameter(p,'plots','both')
addParameter(p,'stablevel',1)
addParameter(p,'markers','.')      %first value is "ALL", second is "STABLE"
addParameter(p,'markersize',[36])  %first value is "ALL", second is "STABLE"
addParameter(p,'colors',{'r' 'k'})
addParameter(p,'picked_lambda',[])
addParameter(p,'picked_error',[])   %note, this is frequency, abs(lambda), and not lambda
addParameter(p,'figure',[])
addParameter(p,'indicator','mac')   %modal indicator (freq or mac)
addParameter(p,'selection',false)
addParameter(p,'grid',[])
addParameter(p,'elements',[])
addParameter(p,'slave',[])
addParameter(p,'active_nodes',[])
addParameter(p,'convert_to_hz',false)
parse(p,varargin{:})

stabcrit=p.Results.stabcrit;
plots=p.Results.plots;
s=p.Results.stablevel;
markers=p.Results.markers;
markersize=p.Results.markersize;
colors=p.Results.colors;
picked_lambda = p.Results.picked_lambda;
picked_error = p.Results.picked_error;
f = p.Results.figure;
indicator=p.Results.indicator;
selection = p.Results.selection;
griddef = p.Results.grid;
elements = p.Results.elements;
slave = p.Results.slave;
active_nodes = p.Results.active_nodes;
convert_to_hz = p.Results.convert_to_hz;

if length(markersize)==1
    markersize(2)=markersize(1);
end

if length(markers)==1
    markers(2)=markers(1);
end

%% FIND STABLE MODES
[lambda_stab, ~, order_stab, ~] = koma.modal.find_stable_poles(lambda, phi, orders, s, stabcrit, indicator);
lambda_vec = cell2mat(lambda');               % all modes
modecount = cellfun(@(x) length(x),lambda);   % corresponding orders
order_all = repelem(orders,modecount);

%% PLOT
leg={};
if isempty(f)
    f=figure;
end
hold on

if convert_to_hz
   freq_scaling = 1/(2*pi);
   freq_unit = 'Hz';
else
   freq_scaling = 1.0;
   freq_unit = 'rad/s';
end

if strcmp(plots,'both') || strcmp(plots,'all')
    a=scatter(abs(lambda_vec)*freq_scaling,order_all,markers(1));
    set(a,'sizedata',markersize(1),'MarkerEdgeColor',colors{1})
    leg{length(leg)+1} = 'All poles';
end

if strcmp(plots,'both') || strcmp(plots,'stable')
    b=scatter(abs(lambda_stab)*freq_scaling,order_stab,markers(2));
    set(b,'sizedata',markersize(2),'MarkerEdgeColor',colors{2})
    leg{length(leg)+1} = 'Stable poles';
    dothandle=b;    % used for mouse click function
else
    dothandle=a;
end


legend(leg)
title({[num2str(s) ' x stability '], ['Criteria: \Delta\omega/\omega < ' num2str(100*stabcrit(1)) '% and \Delta\xi/\xi < ' num2str(100*stabcrit(2)) '% and MAC > ' num2str(100*(1-stabcrit(3))) '%']});
xlabel(['Mode frequency [', freq_unit, ']'])
ylabel('Order');
box

% Plot lines on picked lambdas
if isempty(picked_error)
    picked_error = abs(picked_lambda)*nan;
end

absdist = diff(get(gca,'xlim'));
mindist = absdist/30;

for line = 1:length(picked_lambda)
    lmd=picked_lambda(line);
    err=picked_error(line);
    prevpos(line) = abs(lmd)*freq_scaling;
    distances=abs((abs(lmd)*freq_scaling-prevpos(1:line-1)));
    
    if any(distances<mindist)
        add = sum(distances<mindist)*0.1*max(order);
    else
        add=0;
    end
    
    ln=plot([abs(lmd)*freq_scaling abs(lmd)*freq_scaling],[0,max(order)]);
    t=text(abs(lmd)*freq_scaling,max(order)/2+add,num2str(line));
    e=thirdparty.herrorbar(abs(lmd)*freq_scaling,max(order)/4+add/2,err);
    set(e,'color',ln.Color)
    set(t,'color',ln.Color,'fontsize',10,'backgroundcolor',[1,1,1],'margin',0.5,'horizontalalignment','center')
end

dcm=datacursormode(f);
datacursormode on
set(dcm,'updatefcn',{@infoCursor lambda orders})

if selection == true && ~isempty(griddef) && ~isempty(slave) && ~isempty(elements) && ~isempty(active_nodes)
    geometry.griddef=griddef;
    geometry.elements=elements;
    geometry.slave=slave;
    geometry.active_nodes=active_nodes;
    
    % Enable input while datacursor mode is active
    hManager = uigetmodemanager(f);
    try
        set(hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
    catch
        [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
    end
    set(f,'KeyPressFcn',{@keypress phi geometry})
end

    function [output_txt] = infoCursor(~,event_obj,lambda,order)
        % Display the position of the data cursor
        % obj          Currently not used (empty)
        % event_obj    Handle to event object
        % output_txt   Data cursor text string (string or cell array of strings).
        
        pos = get(event_obj,'Position');
        
        freq_sel = pos(1); %frequency
        nsel = pos(2); %order
        [~,osel] = min(abs(order-nsel)); %find index of order
        [~,msel] = min(abs(abs(lambda{osel})*freq_scaling-freq_sel)); %find number of mode for this order
        xisel = -real(lambda{osel}(msel))/abs(lambda{osel}(msel));
                
        output_txt = {...
            ['order = ',num2str(nsel,3)],...
            ' ',...
            ['freq = ',num2str(freq_sel,4), ' ',freq_unit],...
            ['xi = ', num2str(xisel*100,4) '%'],...
            ' ',...
            ['mode = ', num2str(msel,3)],...
            ['order = ', num2str(order(osel),3)]};
    end

    function keypress(~,event_obj,phi,geometry)       
        switch event_obj.Character
            case 'p'    %plot mode shape
                figure(999)
                koma.vis.plot_mode(phi{osel}(:,msel),geometry.griddef,geometry.elements,geometry.slave,geometry.active_nodes);
                figure(f)
            case 'a'    %append mode
                lambda_sel = [lambda_sel,lambda{osel}(msel)];
                phi_sel = [phi_sel,phi{osel}(:,msel)];
                disp(['Mode is appended. Number of selected modes are now ' num2str(length(lambda_sel)) '.'])
            case 'c'    %clear all selected modes
                disp(['All ' num2str(length(lambda_sel)) ' selected modes are cleared.'])
                lambda_sel = [];
                phi_sel = [];
            case 's'
                nummodes = length(lambda_sel);
                disp([num2str(nummodes) ' modes are selected, and can now be stored using the save dialog.'])
                [file,folder]=uiputfile('*.mat','Save the selected modes as mat-file...');
                if file~=0
                    save(fullfile(folder,file),'lambda_sel','phi_sel')
                    lambda_sel = [];
                    phi_sel = [];
                    disp(['Modes stored in mat-file ' fullfile(folder,file) '. Modes are cleared from RAM.'])
                else
                    disp('Saving cancelled. You may keep appending modes.')
                end
            case 't'
                nummodes = length(lambda_sel);
                disp([num2str(nummodes) ' modes are selected, and are now stored in temporary file temp_modes.mat. Modes are also cleared from RAM.'])
                save('temp_modes.mat','lambda_sel','phi_sel')
                lambda_sel = [];
                phi_sel = [];
        end
        
    end
end