function [ R ] = autocorr_lag(data, maxlag, varargin)
%% Computes cross-correlation matrices for multiple time lags based on input data. Relying on FFT for efficient calculation.
%
% Arguments
% -------------------
% data : double
%   matrix with data, n_samples-by-n_channels (channels are column-wise)  
% maxlag : double
%   maximum number of sample lags    
% mode : 'unbiased', optional
%   mode used to compute ('unbiased' or 'small')
% 
% Returns
% -------------------
% R : double
%   n_channels-by-n_channels-by-n_lags large array, each slice in third dimension 
%   corresponds to the cross-correlation matrix for a given time lag 

% Knut Andreas Kvaale, 2016


p=inputParser;
addParameter(p,'mode','unbiased');   %system orders (for stabilization diagram etc.)

parse(p,varargin{:})
mode=p.Results.mode;

[N, Ndofs]=size(data);
if maxlag>N-1
    warning('Maximum sample lag is larger than total sample length! Reducing to correspond to one below sample length.')
    maxlag=N-1;
end

if strcmp(mode,'small')     %this is faster if the lag is short!
    R=zeros(size(data,2),size(data,2),maxlag);
    lags=0:maxlag;
    
    for k = 1:length(lags)
        lag=lags(k);
        R(:,:,k)= 1/(N-lag)*(data(1+lag:N,:)'*data(1:N-lag,:))';
    end   
else
    R0 = xcorr(data,maxlag,'unbiased');
    R = reshape(R0(maxlag+1:end,:)',[Ndofs, Ndofs, maxlag+1]);
end

if ~issymmetric(R(:,:,1))   %symmetrize first element if not already symmetric
    R(:,:,1) = (R(:,:,1)+R(:,:,1)')/2;
    
end

