function [S, f] = cpsd_matrix(data,fs,nfft,varargin)
%% CPSD_MATRIX Cross-power spectral density matrix estimation using Welch's method.
%
% INPUT
% ---------------------
% data          Data matrix, ordered with channels, nch, as columns
% fs            Sampling frequency [Hz]
% nfft          Samples in each window. If ndiv is given, this is used
%                   instead and the value of nfft is irrelevant.
% [ndiv]        Divisions of time series. If given, input nfft is not used,
%                   and is instead calculated from ndiv.
% [zp]          Zero-padding factor. The FFT segment length will be nfft*zp. 
% [siding]      1 or 2. Defines 1- or 2-sided spectral estimate.
% [nextpow2]    Use nextpow2 for finding appropriate nfft.
% [overlap]     Overlap factor. 
% [window]      Window type. Supports the same as cpsd.
% [show_progress] true/false
%
% OUTPUT
% ----------------------
% S         Matrix (nch x nch x nfreq)
% f         Array (1 x nfreq)
%
%
% Knut Andreas Kvaale, 2017
%

p=inputParser;
addParameter(p,'ndiv',[],@(x) mod(x,1)==0);
addParameter(p,'zp',1,@(x) mod(x,1)==0);
addParameter(p,'siding',1);     %trouble with two-sided :S
addParameter(p,'nextpow2',false);
addParameter(p,'window','hanning');
addParameter(p,'overlap',[])
addParameter(p,'show_progress',[])

parse(p,varargin{:})
ndiv = p.Results.ndiv;
zp = p.Results.zp;
siding = p.Results.siding;
np2 = p.Results.nextpow2;
window_type = lower(p.Results.window);
overlap = p.Results.overlap;
show_progress = p.Results.show_progress;

sidingstring = {'onesided' 'twosided'};
[nsamples,nch] = size(data); 

if ~isempty(ndiv)
   nfft = nsamples/ndiv;
end

if np2 == true
   nfft = 2^nextpow2(nfft);
end

nfft = ceil(nfft);      % length
nfft = nfft+mod(nfft,2); %make even
window_envelope = feval(window_type,nfft);
nsegment = nfft*zp;

totlength = nsegment - (2-siding)*(nfft*zp/2-1);

S = zeros(nch,nch,totlength);

if show_progress
    h = waitbar(0, 'CPSD estimate: 0%');
end

for m = 1:nch
    for n = 1:nch
        [Stmp,f] = cpsd(data(:,m), data(:,n), window_envelope, overlap, nsegment, fs, sidingstring{siding});
        S(m,n,:) = reshape(Stmp,[1,1,totlength]);
        if show_progress
            progress = ((m-1)*nch+n)/(nch^2);
            waitbar(progress, h, sprintf('CPSD estimate: %.1f%%...', progress*100))
        end
    end
end

if show_progress
    close(h)
end

S = S(:,:,1:nsegment/2+1);
f = f(1:nsegment/2+1);

end

