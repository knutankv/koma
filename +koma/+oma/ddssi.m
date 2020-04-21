function [lambda, phi, orders, weighting] = ddssi(data, fs, blockrows, varargin)
%% Data-driven SSI.
%
% Arguments
% ---------------------------
% data : double
%     data matrix, with channels stacked column-wise (n_samples-by-n_channels)
% fs : double
%     sampling frequency
% blockrows : int
%     maximum number of block rows
% orders : [], optional
%     array or list of what orders to include
% weighting : 'upc', optional
%     what weighting type to use ('pc', 'upc', 'cva')
% algorithm : 2, optional
%     what algorithm to use with reference to Van Overschee and de Moor :cite:`VanOverschee1996` (1, 2 [fastest], or 3)
% showinfo : true, optional
%     whether or not to print information during computation
%
% Returns
% ---------------------------
% lambda : double
%     cell array with vectors for all modes corresponding to all orders
% phi : double
%     cell array with matrices containing all mode shapes corresponding to all orders
% orders : int
%     vector with all orders 
% weighting : str
%     the applied weighting algorithm (not always equal input)      
%
% References
% --------------------------------
% * Van Overschee and de Moor :cite:`VanOverschee1996` (1 in code)
% * Rainieri and Fabbrocino :cite:`Rainieri` (2 in code)


p=inputParser;
addParameter(p,'orders',[], @isnumeric);   %system orders (for stabilization diagram etc.)
addParameter(p,'weighting', 'upc');
addParameter(p,'algorithm', 2);
addParameter(p,'showinfo', true)

parse(p,varargin{:})
n=p.Results.orders;
i=blockrows;
weighting = lower(p.Results.weighting);
algorithm = p.Results.algorithm;
showinfo=p.Results.showinfo;

if ~isnumeric(algorithm) || ~(algorithm>=1 && algorithm<=3)
    error('Algorithm input is invalid. 1, 2 and 3 is valid inputs.')
end

if showinfo==true
    disp('*** DATA-DRIVEN SSI ALGORITHM FOR OPERATIONAL MODAL ANALYSIS ****')
end

dt=1/fs;
[N,l] = size(data);

%% ESTABLISH HANKEL MATRICES
if showinfo==true
    disp('* ESTABLISHING HANKEL MATRICES')
end

j = N-2*i+1;    %use all data
Y0 = zeros(2*l*i,j);

for row = 1:2*i
   Y0((row-1)*l+(1:l),:) = data(row:row+j-1,:)'; 
end
Y0 = Y0*1/sqrt(j);  %scale such that consistent with def. of correlation
Yf = Y0(1:l*i,:);
Yp = Y0(l*i+1:end,:);

Yp_plus = Y0(1:l*i+1,:);
Yf_minus = Y0(l*i+1+1:end,:);

%% LQ FACTORIZATION AND PROJECTION MATRICES
if showinfo==true
    disp('* PERFORMING QR/LQ FACTORIZATION AND ESTIMATING PROJECTIONS')
end

[Q,R] = qr(Y0');
Q=Q';
L=R';

L11 = L(1:l*i,1:l*i);
L21 = L(l*i+1:l*i+l,1:l*i);
L22 = L(l*i+1:l*i+l,l*i+1:l*i+l);
L31 = L(l*i+l+1:end,1:l*i);
L32 = L(l*i+l+1:end,l*i+1:l*i+l);

Q1 = Q(:,1:l*i);
Q2 = Q(:,l*i+1:l*i+l);

%%%%-----         STEP 1 in [1]       -----%%%%
Pi = [L21;L31]*Q1';      
Pim1 = [L31,L32] *[Q1';Q2'];
Yii = [L21, L22] * [Q1';Q2'];

%% WEIGHTING MATRICES
if showinfo==true
    disp('* ESTABLISHING WEIGHING MATRICES')
end

switch weighting
    case 'pc'
        phi_pp =  1/j * Yp*Yp';
        W1 = 1;
        W2 = Yp'*phi_pp^(-1/2)*Yp; %cov(Yp) is covariance matrix of past part Yp
    case 'upc'
        W1 = 1;
        W2 = 1;
    case 'cva'
        phi_ff = 1/j * Yf*Yf';
        W1 = phi_ff^(-1/2);    %cov(Yf) covariance matrix of past part Yf
        W2 = 1;
    otherwise   %upc
        weighting='upc';
        W1 = 1;
        W2 = 1;
end
 
%% CALCULATE SVD OF HANKEL/TOEPLITZ MATRICES AND CONTROL/SELECT MAXIMUM ORDER (based on [3])
if showinfo==true
    disp('* CALCULATING SVD OF MATRICES AND CONTROLLING MAXIMUM ORDER')
    disp(['  ** Rank of projection Pi = ' num2str(rank(Pi))])
end

%%%%-----         STEP 2 in [1]       -----%%%%
[U,D,V]=svd(W1*Pi*W2);

if isempty(n)
    figure;
    stem(diag(D));
    
    if showinfo==true
        disp('  ** Please input the system orders in normal MATLAB format, e.g. 2:2:10 or [2,10,40].');
    end
    
    n = input('  ** System orders, n = ');
    if max(n)<l
        warning(['Maximum system order too small. Should be minimum equal to the number of channels. Adjusted to ' num2str(l) '.'])
        n=[n, l];
    end
end
if rank(D)<max(n)
    warning('System order reduced to rank of D matrix.')
    n=n(n<rank(D));
end

if showinfo==true
    disp(['  ** Maximum system order used is n_max = ' num2str(max(n))]);
end

%% COMPUTE STATE MATRIX A FROM DECOMPOSITION FOR EACH ORDER, AND ESTABLISH EIGENVALUES AND EIGENVECTORS (based on [3])
if showinfo==true
    disp('* COMPUTING STATE MATRIX FROM DECOMPOSITION FOR EACH ORDER, AND ESTABLISH EIGENVALUES AND EIGENVECTORS')
end

for j = 1:length(n)  %order indices, j
    o=n(j);          %the current order (jth value of i-array)
    
    % Select subset from H0, via the SVD decomposition, based on order
    D1 = D(1:o,1:o);
    U1 = U(:,1:o);
    V1 = V(:,1:o);
    
    %%%%-----         STEP 3 in [1]       -----%%%%
    O = W1*U1*sqrt(D1);      
    Oup = O(1:end-l,:);     %remove the last l rows from O
    Odown = O(l+1:end,:);   %remove the first l rows from O  
    
    %%%%-----         STEP 4 in [1]       -----%%%%
    Si = pinv(O)*Pi;        %Si, 4.242 in [2]       Kalman filter sequence 
    Sip1 = pinv(Oup)*Pim1;  %S_i+1, 4.244 in [2]    Shifted Kalman filter sequence
    
    %%%%-----         STEP 5 in [1]       -----%%%%
    switch algorithm
        case 1
            AC = [Sip1;Yii]*pinv(Si); 
            A = AC(1:o,:);
            C = AC(o+1:end,:);
        case 2
            A = pinv(Oup)*Odown;
            C = O(1:l,:); 
        case 3
            warning('Algorithm 3 not implemented! 2 used instead');
            A = pinv(Oup)*Odown;
            C = O(1:l,:);             
    end

    % Eigenvalue decomposition [1] and convertion from discrete to
    % continuous poles
    
    [Psi,Lambda0] = eig(A);                 % system eigenvectors and eigenvalues
    lambda0=log(diag(Lambda0))./dt;         % make vector from diagonal matrix and transform from discrete to continuous
    Phi0=C*Psi;                             % observed part of system eigenvectors, referring to input channels
    
    [~,ix]=unique(abs(lambda0));            % find index for unique absolute values
    lambda{j}=lambda0(ix);                  % keep the corresponding eigenvalues
    phi{j}=Phi0(:,ix);                      % and the corresponding eigenvectors
end

orders=n;    %store n as orders

if showinfo==true
    disp('* COMPUTATION COMPLETE!')
end