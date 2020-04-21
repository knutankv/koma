function [wp, phi, xi] = fdd(data, fs, varargin)
%% Frequency-domain decomposition. Old implementation (2014).
%
%
% Arguments
% ---------------------------
% data : double
%     either (i) data matrix, with channels stacked column-wise (n_samples-by-n_channels)
%     or (ii) cross-spectral density matrix with dimensions
%     n_channels-by-n_channels-by-n_freqs
% fs : double
%     either (i) sampling frequency in Hz or (ii) array of frequencies (in
%     rad/s) corresponding to the cross-spectral density input as `data`
%     with dimensions n_freqs
% windows : 2, optional
%   welch divisions
% zp : 1, optional
%   zero-padding factor for spectral density estimate
% wlim : Inf, optional
%   freq. area of interst in rad/s: either input as [wmin, wmax] or wmax (=> wmin = 0)
% plot : false, optional
%   plot results
% lowpass : [], optional
%   lowpass frequency
% lowpassorder : 5, optional
%   butterworth lowpass filter order
% mac : 0.9, optional
%   MAC criterion for mode to be deemed okay
% lines : 4, optional
%   number of frequency lines to check from SVD
% damppoints : 3, optional
%   points used in linear regression for damping estimate (EFDD)
% dampZP : 1, optional
%   zero-padding for correlation estimate used to calculate damping
%   (log. decrement)
% peakfactor : 0.8, optional
%   requirement that the `peakfactor` selected peak of a mode is at least
%   peakfactor * the maximum singular value found within the region 
%   surrounding the peak that has mac > mac_crit.
% minlength : 10, optional
%   minimum number of points in spectrum that is within MAC criterion
%   to accept as mode
%
% Returns
% ---------------------------
% wp : double
%     1d array of peak frequencies in rad/s
% phi : double
%     modal transformation matrix (mode shapes stacked column-wise)
% xi : double
%     1d array of critical damping ratios
%
% References
% --------------------------------
% * Brincker et al. :cite:`Brincker2000`
% * Rainieri and Fabbrocino :cite:`Rainieri`
%


p=inputParser;
addParameter(p,'windows',2,@isnumeric);  %Welch divisions
addParameter(p,'zp',1,@isnumeric);  %zero-padding factor
addParameter(p,'wlim',Inf,@isnumeric);  %frequency area of interest
addParameter(p,'lowpass',[],@isnumeric);    %frequency of lowpass
addParameter(p,'lowpassorder',5,@isnumeric);%order of lowpass (Butterworth filter)
addParameter(p,'plot',false)
addParameter(p,'mac',0.95)
addParameter(p,'lines',4)
addParameter(p,'damppoints',3)
addParameter(p,'dampZP',1)
addParameter(p,'peakfactor',0)
addParameter(p,'minlength',10)

parse(p,varargin{:})
Nwelch=p.Results.windows;
zp=p.Results.zp;
wlim=p.Results.wlim;
plotit=p.Results.plot;
flowpass=p.Results.lowpass;
lporder=p.Results.lowpassorder;
mac_criterion=p.Results.mac;
lines=p.Results.lines;
damppoints=p.Results.damppoints;
dampZP=p.Results.dampZP;
peakfactor=p.Results.peakfactor;
minLength=p.Results.minlength;

if size(data,3) ~= 1 && length(size(data))==3 && size(data,1)==size(data,2)
    %spectra is input
    [Ndofs,~,~] = size(data);
    S=data;
    w=fs;
    wmax0=max(w);
elseif size(data,3) == 1
    %time series matrix is input
    [Ntime,Ndofs]=size(data);
    w=[];
    wmax0=fs*Ntime;
else
    error('Invalid format of data array! Either input time series data matrix with Ntime x Ndofs dimensions, or spectral density matrix with dimensions Ndofs x Ndofs x Nfreq. If spectra is input, fs = frequency axis in rad/s')
end


lines=min(Ndofs,lines);     %adjust lines in case number of dof is small


if length(wlim) == 1
    wmax=min(wlim,wmax0);
    wmin=0;
else
    wmax=max(wlim);
    wmin=min(wlim);
end

maxMAC=0.5;              %If mac>maxMAC: the modes are considered equal

%% DETREND, LOWPASS AND RESAMPLE
if isempty(w)   %time series has been input - and not spectra
    data=detrend(data,'constant');
    dt=1/fs;
    t=0:dt:dt*(size(data,1)-1);
    
    if ~isempty(flowpass)
        RF = floor(0.5*(fs/flowpass));        %resample at 2x the lowpass frequency (so that f_Nyquist = f_lowpass!)
        data = lowpass(data,dt,flowpass,lporder);   %lowpass
        data = data(1:RF:end,:);              %enforce resample factor
        t = t(1:RF:end);
    end
    dt=t(2)-t(1);
    data=detrend(data,'constant');
    
    
    %% ESTABLISH SPECTRAL DENSITY USING WELCH
    h=waitbar(0,'Initializing calculation...','Name','Calculating covariance matrix');
    
    %Initialize
    [~,w] = welchCPSD(data(:,1),data(:,1),dt,Nwelch,zp);
    S=zeros(Ndofs,Ndofs,length(w));
    
    for i = 1:Ndofs
        progress=i/Ndofs;
        waitbar(progress,h,sprintf('%1.0f%%',progress*100),'Name','Calculating covariance matrix')
        for j = 1:Ndofs
            S0 = welchCPSD(data(:,i),data(:,j),dt,Nwelch,zp);
            S(i,j,:)=S0;
        end
    end
    close(h);
end

%% FIND INDEX CORRESPONDING TO wmax AND wmin
[~,imax]=min(abs(w-wmax));
[~,imin]=min(abs(w-wmin));

S=S(:,:,imin:imax);
w=w(imin:imax);
kmax=length(w);
%% SINGULAR VALUE DECOMPOSITION
for k=1:kmax
    [U(:,:,k),D(:,:,k),~] = svd(S(:,:,k));
end
clrs=hsv(lines);

if plotit==true
    figure;
    
    Ssum=zeros(kmax,1);
    for k=1:kmax
        Ssum(k) = Ssum(k) + sum(diag(squeeze((S(:,:,k)))));
    end
    h_sp1=subplot(2,1,1);
    plot(w,Ssum/Ndofs);
    ylabel('\Sigma S_{ii} (\omega)')
    xlabel('\omega [rad/s]')
    
    h_sp2=subplot(2,1,2);
    
    for i = 1:lines
        plot(w,squeeze(D(i,i,:)),'color',clrs(i,:));   %plot singular values
        hold on;
        legendtxt{i} = ['Singular value line ' num2str(i)];
    end
    
    set(h_sp2,'YScale','log');
    set(h_sp1,'YScale','log');
    
    ylabel('Singular value lines, D_{ii}(\omega)');
    xlabel('\omega [rad/s]');
    legend(legendtxt);
end

%% PICK NATURAL FREQUENCIES
D1=squeeze(D(1,1,:));
nrange=startAndStopGUI(w,D1);
hold on
Ns=size(nrange,1); %number of selected ranges

%% RUN THROUGH ALL CHOSEN RANGES TO ESTABLISH MODES
phi=[];
wp=[];
mode=0;

for n = 1:Ns
    for l = 1:lines
        Dl=squeeze(D(l,l,:));
        [pks,np0] = findpeaks(Dl(nrange(n,1):nrange(n,2))); %find peaks within range
        [pmax,pp0]=max(pks);                                   %find the largest of the peaks
        np0=np0(pp0);                                       %establish the index corresponding to this peak, relative to the chosen range
        
        np = nrange(n,1)+np0;                               %establish index for peak, using global indexing
        
        %GO DOWNWARDS TO ESTABLISH LIMITS FOR "SDOF"-APPROXIMATION
        for ndown = 1:np-1
            macdown=modass(U(:,l,np),U(:,l,np-ndown));
            if macdown<mac_criterion
                break
            end
        end
        if isempty(ndown)
            ndown=0;
        end
        
        %GO UPWARDS TO ESTABLISH LIMITS FOR "SDOF"-APPROXIMATION
        for nup = 1:length(w)-np
            macup=modass(U(:,l,np),U(:,l,np+nup));
            if macup<mac_criterion
                break
            end
        end
        
        if isempty(nup)
            nup=0;
        end
        
        nrangeThisLine=(np-ndown):(np+nup);
        validMode=pmax>=peakfactor*max(Dl(nrangeThisLine)); %check if the peak value is minimum "peakfactor" times as big as the maximum value within the valid range
        
        macAll=0;
        if size(phi,2)>0
            for i=1:size(phi,2)
                macAll(i) = modass(U(:,l,np),phi(:,i));
            end
        end
        
        if length(nrangeThisLine)>minLength && ~isempty(macAll(macAll<maxMAC)) && validMode    %if 2 values outside peak are "same mode" as given by mac_criterion, check also that no similar modes exist!
            if size(phi,2)==mode
                mode=mode+1;
                phi(:,mode) = zeros(1,Ndofs);
            end
            phi(:,mode) = phi(:,mode) + U(:,l,np) * D(l,l,np);  %add to current mode vector
            wp(mode)=w(np);
            modeLines(mode)=l;
            modeRanges{mode}=nrangeThisLine;
            modeIndices(mode)=np;
%             plot(h_sp2,w(nrangeThisLine),Dl(nrangeThisLine),'-','color',[0.7 0.7 0.7],'linewidth',1.5);
            [R,t] = S2R(squeeze(D(l,l,nrangeThisLine)),w(nrangeThisLine),'ZP',dampZP,'window',true);  %estimate correlation function
            R=R-mean(R(round(0.1*end):end));
            [xi(mode),wn(mode)]=logdec(R/R(1),t,'plot',true,'origin',true,'points',damppoints);
            macWidth(mode)=length(nrangeThisLine);
        end
    end
end

Nmodes=length(wp);
if Nmodes==0
    error('No modes satisfying the MAC-criterion were found! Aborted...')
end

%% RUN THROUGH PICKED MODES SATISFYING THE MAC CRITERION, TO REMOVE EQUAL MODES
macAll=modass(phi,phi);  %all MAC numbers phi vs phi
keepModes=[];

for m = 1:Nmodes
    equalModes = find(macAll(m,:)>maxMAC); %modes matching (by means of MAC) more than specified MAC-criterion
    [~,widestMode] = max(macWidth(equalModes));
    keepModes(m) = equalModes(widestMode);
end

keepModes=unique(keepModes);
phi=phi(:,keepModes);
wp=wp(keepModes)';
xi=xi(keepModes)';
modeLines=modeLines(keepModes);
modeRanges=modeRanges(keepModes);
modeIndices=modeIndices(keepModes);
%% PLOT MODES

if plotit==1
    for m=1:length(wp)
        range=modeRanges{m};
        Dl=squeeze(D(modeLines(m),modeLines(m),:));
        plot(h_sp2,w(range),Dl(range),'-','color',[0.7 0.7 0.7],'linewidth',1.5);
        text(w(modeIndices(m)),Dl(modeIndices(m))*1.5, num2str(m),'Parent', h_sp2,'FontSize',9,'HorizontalAlignment','center');
    end
end
