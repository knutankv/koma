function [lambda_picked, phi_picked, statistics] = pick_stable_modes(lambda_stab, phi_stab, slack, varargin)
%% Pick stable modes from arrays of lambdas.
% Use like this: [lambda_picked, phi_picked, statistics] = pick_stable_modes(lambda_stab, phi_stab, slack, ...)
%
%
% Arguments
% ------------------------------------
% lambd_stab : double
%    array with complex-valued eigenvalues deemed stable by
%    find_stable_poles
% phi_stab : double
%    2d array with complex-valued eigenvectors (stacked as columns) deemed 
%    stable by `find_stable_poles`, each column corresponds to a mode
% slack : double
%   allowed relative deviance between different modes input for them to be 
%   clustered together as the same mode (3-by-1, [freq, xi, mac])
% phi_ref : double, optional
%   reference modal matrix, MAC slack value corresponds to slack from this 
%   if input. If not wanted, use []. Note slack(1) and slack(2) not used in this mode!
% maxmacdeviance : 0.2, optional
%   restricts the deviance from the ref. mode shape (value of 0.2 => MAC>0.8)
% onlyone : true, optional
%   if true,  only one mode is identified per reference - otherwise, 
%   multiple modes within slack limits are selected
% mean : true, optional
%   if true, the statistical mean of all within slack is used instead of 
%   the single pole with largest MAC. The standard deviations are reported 
%   in lambda_std and phi_std. Only applicable when phi_ref is input.
% slack_mode : 'all', optional
%    'all' (all have to be within) or 'one' (okay if only one is)
% 
% Returns
% ---------------------
% lambda_picked : double
%   array with selected lambdas
% phi_picked : double
%   phi with stable mode shapes (stacked column-wise)
% statistics : struct
%   structure with statistics for each mode cluster

p=inputParser;
addParameter(p,'onlyone',true,@islogical);
addParameter(p,'phi_ref',[])
addParameter(p,'maxmacdeviance',0.3)
addParameter(p,'mean',true,@islogical)
addParameter(p,'slack_mode','all')

parse(p,varargin{:})
onlyone = p.Results.onlyone;
phi_ref = p.Results.phi_ref;
maxmacdeviance = p.Results.maxmacdeviance;
meanmode = p.Results.mean;
slack_mode = p.Results.slack_mode;

statistics=[];

if isempty(phi_ref)
    [sorted_freqs,six]=sort(abs(lambda_stab));
    sorted_phi = phi_stab(:,six);
    sorted_lambda = lambda_stab(six);
    sorted_damp = -real(lambda_stab(six))./abs(lambda_stab(six));
    MAC = koma.modal.xmacmat(sorted_phi);
    modecount=1;
    mode=1;
    lambda_picked=[];
    phi_picked=[];
    
    while mode<length(lambda_stab)
        okmacind = find(MAC(:,mode)>=(1-slack(3)));
        freqs = sorted_freqs(okmacind);
        damps = sorted_damp(okmacind);
        okfreqind = abs((mean(freqs)-freqs)./mean(freqs))<slack(1);
        okdampind = abs((mean(damps)-damps)./mean(freqs))<slack(2);
        
        okind = okmacind(bsxfun(@and,okfreqind,okdampind)); %all that are okay
        if ~isempty(okind)
            lambda_picked(modecount,1) = mean(sorted_lambda(okind));
            phi_picked(:,modecount) = mean(sorted_phi(:,okind),2);
            mode = okmacind(end)+1;
            modecount=modecount+1;
        else
            mode = mode+1;
        end
    end
else    %phi_ref is input (reference mode)
     
    minmac=(1-maxmacdeviance);
    macs = modal.xmacmat(phi_ref,phi_stab);
    [macmax,maxidx]=max(macs,[],2);          %find indices amongst the modes for all reference modes matching most (number corresponds to index in modes), first element = mode matching reference mode 1 most, etc.
    okidx_logic = macmax>=minmac;            %referring to phi_ref
    modeidx = maxidx(okidx_logic);

    if onlyone==true
        modeidx_temp = maxidx(okidx_logic);
        modeidx=nan(size(modeidx_temp));

        for m = 1:length(modeidx_temp)
            group=find(modal.slack_check(lambda_stab(modeidx_temp(m)),phi_stab(:,modeidx_temp(m)),lambda_stab,phi_stab,slack,slack_mode));
            [common_modeidx,~,idx_common_from_modeidx]=intersect(group,modeidx_temp);
            [~,sel]=max(diag(macs(idx_common_from_modeidx,common_modeidx)));
            if idx_common_from_modeidx(sel)==m
                modeidx(m) = common_modeidx(sel);
            end
        end
    end
    okidx = find(okidx_logic);
    okidx(isnan(modeidx)) = []; %remove nan-entries
    modeidx(isnan(modeidx)) = [];
    
    % Start with nans for all modal values
    phi_picked = nan(size(phi_ref));
    lambda_picked = nan(size(phi_ref,2),1);
    
    if meanmode==true
        % Initialize stats with nan values
        statistics.phi.std = phi_picked;
        statistics.w.std = lambda_picked;
        statistics.w.mean = lambda_picked;
        statistics.xi.std = lambda_picked;
        statistics.xi.mean = lambda_picked;
        statistics.N = lambda_picked;
        statistics.mac.std = lambda_picked;
        statistics.mac.mean = lambda_picked;
        
        temp = cell(size(lambda_picked));
        temp(:) = {nan};
        statistics.lambda = temp;
        damp = -real(lambda_stab)./abs(lambda_stab);
        freq = abs(lambda_stab);  
        
        % Loop through all ok indices
        for i=1:length(modeidx)
            store_idx = okidx(i);   %store index
            ret_idx = modeidx(i);   %retrieve index

            freq_func = abs((freq-freq(ret_idx)))./freq(ret_idx); %frequency deviance functional
            damp_func = abs((damp-damp(ret_idx)))./damp(ret_idx); %damping deviance functional
            mac_func = 1-modal.xmacmat(phi_stab,phi_stab(:,ret_idx)); %MAC deviance functional
            
            within_slack=find(modal.slack_check(lambda_stab(ret_idx),phi_stab(:,ret_idx),lambda_stab,phi_stab,slack,slack_mode));
            these_macs=macs(store_idx,within_slack);    
            ok_macs = (1-these_macs)<=maxmacdeviance;
            within_slack = within_slack(ok_macs); %ensure to throw away all macs deviating more than maxmacdeviance from the reference
            
            lambda_picked(store_idx) = mean(lambda_stab(within_slack));
            
            statistics.w.std(store_idx) = nanstd(abs(lambda_stab(within_slack)));
            statistics.xi.std(store_idx) = nanstd(-real(lambda_stab(within_slack))./abs(lambda_stab(within_slack)));

            statistics.w.mean(store_idx) = nanmean(abs(lambda_stab(within_slack)));
            statistics.xi.mean(store_idx) = nanmean(-real(lambda_stab(within_slack))./abs(lambda_stab(within_slack))); 
            
            if any(within_slack)
                statistics.lambda{store_idx} = lambda(within_slack);
            end         
            
            statistics.mac.mean(store_idx) = mean(macs(store_idx,within_slack));
            statistics.mac.std(store_idx) = std(macs(store_idx,within_slack));
            statistics.N(store_idx) = sum(within_slack);
        end
        phi_picked(:,okidx(~isnan(modeidx))) = phi_stab(:,modeidx(~isnan(modeidx)));
        
    else
        phi_picked = phi_ref*nan;
        lambda_picked = nan(size(phi_ref,2),1);

        lambda_picked(okidx(~isnan(modeidx))) = lambda(modeidx(~isnan(modeidx)));
    end
end
