function [lambda, Phi, order, weighting] = covssi(data, fs, blockrows, varargin)
%% Covariance-driven SSI.
%
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
% weighting : 'none', optional
%     what weighting type to use ('none' or 'br', or 'cva')
% matrix_type : 'matrix_type', optional
%     what matrix type formulation to base computation on ('hankel' or 'toeplitz')
% algorithm : 'shift', optional
%     what algorithm to use ('shift' or 'standard' - when using 'standard' the method is equivalent to NExt-ERA)
% showinfo : True, optional
%     whether or not to print information during computation
% balance : True, optional
%     whether or not to conduct balancing to the cross-correlation matrices prior to matrix operations (Cholesky and SVD)
% 
% Returns
% ---------------------------
% lambd : double
%     cell array with vectors for all modes corresponding to all orders
% phi : double
%     cell array with matrices containing all mode shapes corresponding to all orders
% order : int
%     vector with all orders 
% weighting : str
%     the applied weighting algorithm (not always equal input)      
%
% References
% --------------------------------
% * Hermans and van Der Auweraer :cite:`Hermans1999` (1 in code)
% * Van Overschee and de Moor :cite:`VanOverschee1996` (2 in code)
% * Rainieri and Fabbrocino :cite:`Rainieri` (3 in code)
% * Pridham and Wilson :cite:`Pridham2003` (4 in code)

p=inputParser;
addParameter(p,'order',[],@isnumeric);   %system orders (for stabilization diagram etc.)
addParameter(p,'weighting','none');
addParameter(p,'showinfo',true,@islogical)
addParameter(p,'matrix_type','hankel')
addParameter(p,'algorithm','standard')

parse(p,varargin{:})
n=p.Results.order;
i=blockrows;
weighting = p.Results.weighting;
showinfo=p.Results.showinfo;
matrix_type = p.Results.matrix_type;
algorithm = p.Results.algorithm;

if showinfo==true
    disp('*** COVARIANCE DRIVEN SSI ALGORITHM FOR OPERATIONAL MODAL ANALYSIS ****')
end

dt=1/fs;
[~, l] = size(data);

%% ESTABLISH TOEPLITZ/HANKEL AND LAG-SHIFTED TOEPLITZ/HANKEL MATRIX (notation from [1])
% H0=zeros(i*l,i*l);
% H1=zeros(i*l,i*l);
Rp=zeros(i*l,i*l);

if showinfo==true
    disp('* ESTABLISHING HANKEL AND TOEPLITZ MATRICES')
    disp('  ** Correlation estimation')
end

% Establish fastest correlation function estimator
if i>30
    corrmode='large';
else
    corrmode='small';
end
R = autocorr_lag(data,2*i,'mode',corrmode);

if showinfo==true
    disp('  ** Matrix stacking')
end

if strcmpi(matrix_type,'hankel')
    for row=1:i            %go through all block rows from 1 to i
        for col=1:i        %go through all block columns from 1 to i
            H0((row-1)*l+(1:l),(col-1)*l+(1:l)) = R(:,:,row+col);             %Hankel/Toeplitz matrix, H0
            H1((row-1)*l+(1:l),(col-1)*l+(1:l)) = R(:,:,row+col+1);             %Lag-shifted Hankel/Toeplitz matrix, H1
        end
    end
elseif strcmpi(matrix_type,'toeplitz')
    for row=1:i
        for col=1:i
             H0((row-1)*l+(1:l),(col-1)*l+(1:l)) = R(:,:,i-col+row+1);       %Toeplitz
             H1((row-1)*l+(1:l),(col-1)*l+(1:l)) = R(:,:,i-col+row+2);       %Lag-shifted Toeplitz        
        end
    end
else
    error('Unknown matrix string :(')
end

% Establish Rp if CVA weighting is asked for
if strcmp(weighting,'cva')
    if showinfo==true
        disp('  ** Establishing R+ and R- for weighting matrices')
    end
    for row=1:i
        for col=1:i
            if row<col
                Rp((row-1)*l+(1:l),(col-1)*l+(1:l)) = R(:,:,abs(row-col)+1)';   %R+ matrix
                Rm((row-1)*l+(1:l),(col-1)*l+(1:l)) = R(:,:,abs(row-col)+1);
            elseif row>=col
                Rp((row-1)*l+(1:l),(col-1)*l+(1:l)) = R(:,:,abs(row-col)+1);    %R- matrix
                Rm((row-1)*l+(1:l),(col-1)*l+(1:l)) = R(:,:,abs(row-col)+1)';
            end
        end
    end
    
    % Improve conditioning
    [T,H0] = balance(H0);
    Rp=T\Rp*T;
    Rm=T\Rm*T;
end


%% CALCULATE SVD OF HANKEL/TOEPLITZ MATRICES AND CONTROL/SELECT MAXIMUM ORDER (based on [3])
if showinfo==true
    disp('* CALCULATING SVD OF MATRICES AND CONTROLLING MAXIMUM ORDER')
    disp(['  ** Rank(D) = ' num2str(rank(H0))])
end

[U,D,V]=svd(H0);

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

%% ESTABLISH WEIGHTING MATRIX
if showinfo==true
    disp('* ESTABLISHING WEIGHTING MATRIX')
end

switch lower(weighting)
    case 'cva'   % [2] used, [1] results in ill-conditioning issues
        try
            L1 = gaxpy_chol(balance(Rp));            
            L2 = gaxpy_chol(balance(Rm));
            [U,D,V] = svd(L1\H0/L2);
            if strcmp(algorithm,'standard')
                warning('Solution algorithm STANDARD is not supported with CVA. Changing to SHIFT.')
                algorithm='shift';
            end
            
        catch err
            warning(['CVA failed! ' err.message ' Continuing with no weighting (BR). If you want to run with CVA, try adjusting number of block rows.'])
            L1=eye(size(D,1));
            weighting = 'br';
        end
    case 'br' % [1]
        L1 = eye(size(D,1));
    case 'none' % [4]
        L1 = 1;
    otherwise
        error('Invalid weighting algorithm requested')
end

%% COMPUTE STATE MATRIX A FROM DECOMPOSITION FOR EACH ORDER, AND ESTABLISH EIGENVALUES AND EIGENVECTORS (based on [3])
if showinfo==true
    disp('* COMPUTING STATE MATRIX FROM DECOMPOSITION FOR EACH ORDER, AND ESTABLISH EIGENVALUES AND EIGENVECTORS')
end

for j = 1:length(n)  %order indices, j
    o=n(j);          %the current order (jth value of i-array)
    
    % Select subset from H0, via the SVD decomposition, based on order
    D1=D(1:o,1:o);
    U1=U(:,1:o);
    V1=V(:,1:o);
    
    % State matrix estimate
    O = L1*U1*sqrt(D1);     %O, observability (T=identity matrix)
    C = O(1:l,:);           %Pick the l first values
    
    Oup = O(1:end-l,:);     %remove the last l rows from O
    Odown = O(l+1:end,:);   %remove the first l rows from O
    
    if strcmp(algorithm,'standard') %this means NExT-ERA is chosen
        A = sqrt(D1)\U1'*H1*V1/sqrt(D1);        % [4] (and [3])
    elseif strcmp(algorithm,'shift')
        A = pinv(Oup)*Odown;                    % [1]
    end
    
    % Eigenvalue decomposition [1] and convertion from discrete to
    % continuous poles
    
    [Psi,Lambda0] = eig(A);                 %system eigenvectors and eigenvalues
    lambda0=log(diag(Lambda0))./dt;         %make vector from diagonal matrix and transform from discrete to continuous
    Phi0=C*Psi;                              %observed part of system eigenvectors, referring to input channels
    
    [~,ix]=unique(abs(lambda0));            %find index for unique absolute values
    lambda{j}=lambda0(ix);                  %keep the corresponding eigenvalues
    Phi{j}=Phi0(:,ix);                      %and the corresponding eigenvectors
end
order=n;    %store n as order

if showinfo==true
    disp('* COMPUTATION COMPLETE!')
end