function [ model ] = BTRTF( TY, varargin )
% BTRTF: Bayesian Low-Tubal-Rank Robust Tensor Factorization
%
% %[Syntax]%: 
%   model = BTRTF( TY )
%   model = BTRTF( ___, Name, Value )
%
% %[Inputs]%:
%   TY:            the 3rd-order input tensor of size I1 x I2 x I3
%     
% %[Name-Value Pairs]
%   'init':        the initialization type
%                  'ml':   Initialization with tSVD (default)
%                  'rand': Random initialization
%
%   'Rs':          the initialized multi-rank
%                  min(I1, I2)*ones(I3, 1) (Default)
%
%   'initVar':     the initialized variance of outliers
%                  1 (Default)
%
%   'gamma':       the relaxation parameter \gamma
%                  I3 (Default)
%
%   'maxIters':    the maximum number of iterations
%                  500 (Default)
%
%   'tol':         the tolerance of the relative change of log-likelihood
%                  1e-5 (Default)
%
% %[Outputs]%:
%   model.TX:      the low-rank component
%   model.S:       the sparse component
%   model.Rs:      the estimated multi-rank
%   model.Fit:     the model fitness at each iteration
%                       
% %[Reference]%:            
%   Yang Zhou, Yiu-ming Cheung. 
%   Bayesian Low-Tubal-Rank Robust Tensor Factorization with Multi-Rank Determination. 
%   IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%   DOI:10.1109/TPAMI.2019.2923240, 2019.
%                             
% %[Author Notes]%   
%   Author:        Yang ZHOU
%   Email :        youngzhou12@gmail.com
%   Affiliation:   Department of Computer Science
%                  Hong Kong Baptist University
%   Release date:  June 22, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if nargin < 1, error('Not enough input arguments.'); end
N = ndims(TY); % Order of the input tensor
if N ~= 2 && N ~= 3, error('Invalid inputs.'); end
Is = size(TY); % Dimensions of the input tensor
if N == 2, Is(3) = 1; end
I = prod(Is);

ip = inputParser;
ip.addParameter('init', 'rand', @(x) (ismember(x,{'ml','rand'})));
ip.addParameter('Rs', min(Is(1:end-1))*ones(Is(3),1), @(x) (isvector(x) && max(x) <= min(Is(1:end-1))));
ip.addParameter('initVar', 1, @isscalar);
ip.addParameter('gamma', Is(3), @isscalar);
ip.addParameter('maxIters', 100, @isscalar);
ip.addParameter('tol', 1e-5, @isscalar);
ip.parse(varargin{:});

init    = ip.Results.init;
Rs      = ip.Results.Rs;        % Initialized multi-rank 
initVar = ip.Results.initVar;   % Initialized variance of outliers
gamma   = ip.Results.gamma;     % Relaxation parameter
maxI    = ip.Results.maxIters;  % Default max number of iterations
tol     = ip.Results.tol;       % Threshold of convergence

% Data Preprocessing
NormY = norm(TY(:));
Y2sum = sum(TY(:).^2);
scale2 = Y2sum / I;
scale = sqrt(scale2);
fftY = fft(TY,[],3);

% Initialization
a_lmbd0 = 1e-6; b_lmbd0 = 1e-6;
a_beta0 = 1e-6; b_beta0 = 1e-6;
a_tau0 = 1e-6; b_tau0 = 1e-6;

switch init
    case 'rand'	% Random initialization        
        initR = max(Rs);
        U = randn(Is(1), initR, Is(3)) * sqrt(scale);
        V = randn(Is(2), initR, Is(3)) * sqrt(scale);
        
        fftU = fft(U, [], 3);
        fftV = fft(V, [], 3);
        fftX = zeros(Is);
        
        % T-product in the frequency domain
        for k = 1:Is(3)
            fftX(:,:,k) = fftU(:,:,k)*fftV(:,:,k)';
        end
        
        fftSigU = scale*repmat( eye(initR), [1 1 Is(3)]);
        fftSigV = scale*repmat( eye(initR), [1 1 Is(3)]);
        
        fftU = squeeze(mat2cell(fftU,Is(1),initR,ones(Is(3),1)));
        fftV = squeeze(mat2cell(fftV,Is(2),initR,ones(Is(3),1)));
        
        fftSigU = squeeze(mat2cell(fftSigU,initR,initR,ones(Is(3),1)));
        fftSigV = squeeze(mat2cell(fftSigV,initR,initR,ones(Is(3),1)));
        
        lambda = scale*ones(initR,Is(3));
        lambda = mat2cell(lambda,initR,ones(Is(3),1));        
        tau = 1./scale2;
        
    case 'ml'   % Initialization with tSVD
        fftU = cell(Is(3),1); 
        fftV = cell(Is(3),1); 
        fftSigU = cell(Is(3),1); 
        fftSigV = cell(Is(3),1); 
        lambda = cell(1,Is(3)); 
        fftX = zeros(Is);
        
        for k = 1:Is(3)
            [Uk,Sk,Vk]=svds(fftY(:,:,k),Rs(k));
            fftU{k} = Uk*(Sk.^(0.5));
            fftV{k} = Vk*(Sk.^(0.5));
            fftX(:,:,k) = fftU{k}*fftV{k}';
            fftSigU{k} = scale*eye(Rs(k));
            fftSigV{k} = scale*eye(Rs(k));
            lambda{k} = scale*ones(Rs(k),1);
        end
        TX = ifft(fftX, [], 3);
        tau = 1./scale2;
end

% Initialize the Sparse Component S
betas = initVar.^(-1)*ones(Is).*((a_beta0+eps)/(b_beta0+eps));
S = betas.^(-0.5).*rand(Is);
fftS = fft(S,[],3);
 
expUUc = cell(Is(3),1);    
expVVc = cell(Is(3),1); 

SigUSigV = cell(Is(3),1);
UUcSigV = cell(Is(3),1);
VVcSigU = cell(Is(3),1);

for k = 1:Is(3)
    expVVc{k} = Is(2)*Is(3)*fftSigV{k} + fftV{k}'*fftV{k};
end

Fit = zeros(maxI,1); % Model fitness at each iteration
FitOld = 1e-8;
for iI = 1:maxI  
    old_fftX = fftX;
    old_fftE = fftS;

    fftX = fftY - fftS;
    if Is(3) ~= 1
        endIdx = int16(Is(3)/2 + 1);     
    else
        endIdx = 1;
    end    
    for k = 1:endIdx
        % Update U
        fftSigU{k} = eye(Rs(k))/(tau*expVVc{k} + FitOld/gamma*diag(lambda{k}) + 1e-5*eye(Rs(k)));
        fftU{k} = tau*fftX(:,:,k)*fftV{k}*fftSigU{k}; 
        fftUU = fftU{k}'*fftU{k};
        expUUc{k} = Is(1)*Is(3)*fftSigU{k} + fftUU;   

        % Update V
        fftSigV{k} = eye(Rs(k))/(tau*expUUc{k} + FitOld/gamma*diag(lambda{k}) + 1e-5*eye(Rs(k)));
        fftV{k} = tau*fftX(:,:,k)'*fftU{k}*fftSigV{k};         
        fftVV = fftV{k}'*fftV{k};
        expVVc{k} = Is(2)*Is(3)*fftSigV{k} + fftVV;
        
        % Update Lambda
        a_lmbd = a_lmbd0 + 0.5*(Is(1)+Is(2)); 
        b_lmbd = (b_lmbd0 + 0.5/(Is(3))*diag(expUUc{k}+expVVc{k}));
        lambda{k} = a_lmbd./b_lmbd;
        
        SigUSigV{k} = fftSigU{k}*fftSigV{k};
        UUcSigV{k} = fftUU*fftSigV{k};
        VVcSigU{k} = fftVV*fftSigU{k};
        
        fftX(:,:,k) = fftU{k}*fftV{k}';
    end
    
    for k = Is(3):-1:endIdx+1        
        expUUc{k} = conj(expUUc{Is(3)-k+2});
        fftSigU{k} = conj(fftSigU{Is(3)-k+2});
        fftU{k} = conj(fftU{Is(3)-k+2});
        
        expVVc{k} = conj(expVVc{Is(3)-k+2});
        fftSigV{k} = conj(fftSigV{Is(3)-k+2});
        fftV{k} = conj(fftV{Is(3)-k+2});
        
        lambda{k} = conj(lambda{Is(3)-k+2});
        
        SigUSigV{k} = conj(SigUSigV{Is(3)-k+2});
        UUcSigV{k} = conj(UUcSigV{Is(3)-k+2});
        VVcSigU{k} = conj(VVcSigU{Is(3)-k+2});
        
        fftX(:,:,k) = conj(fftX(:,:,Is(3)-k+2));
    end

    % Update the Sparse Component S
    TX = ifft(fftX,[],3);
    SigmaS = 1./(betas+tau);
    S = double(tau*(TY-TX).*SigmaS);
    fftS = fft(S, [], 3);
    
    % Update Beta
    a_beta = a_beta0 + 0.5;
    b_beta = b_beta0 + 0.5*(S.^2 + SigmaS);
    betas = a_beta./b_beta;
    
    % Update tau    
    err = norm(fftY(:)-fftX(:)-fftS(:))^2/Is(3) ...
        + Is(2)*sum(arrayfun(@(x) trace(x{:}), UUcSigV)) ...
        + Is(1)*sum(arrayfun(@(x) trace(x{:}), VVcSigU)) ...
        + I*sum(arrayfun(@(x) trace(x{:}), SigUSigV)) + sum(SigmaS(:));

    a_tau = a_tau0 + 0.5*I;
    b_tau = b_tau0 + 0.5*err;
    tau = a_tau / b_tau;         

    % Calc Fittness
    Fit(iI) = 1 - sqrt(err)/NormY;
    RelChan = abs(FitOld - Fit(iI));
    if Fit(iI) > 0
        FitOld = Fit(iI);    
    else
        FitOld = 1e-8;
    end
     
    % Prune Unnecessary Components
    for k = 1:Is(3)
        MAX_Lmbd = min(lambda{k}) * 1e4;
        dimIdx = find(lambda{k} <= MAX_Lmbd);
        if length(dimIdx) ~= Rs(k)
            Rs(k) = length(dimIdx);
            lambda{k} = lambda{k}(dimIdx);
            fftU{k} = fftU{k}(:,dimIdx);
            fftV{k} = fftV{k}(:,dimIdx); 
            expVVc{k} = expVVc{k}(dimIdx,dimIdx);
            expUUc{k} = expUUc{k}(dimIdx,dimIdx);
            fftSigU{k} = fftSigU{k}(dimIdx,dimIdx);
            fftSigV{k} = fftSigV{k}(dimIdx,dimIdx);
        end          
        cmpNorm = norm(fftU{k},'fro') + norm(fftV{k},'fro'); 
        if Rs(k) > 1 && cmpNorm < 1e-12 && FitOld > 0.98
            Rs(k) = 1;
            lambda{k} = lambda{k}(1);
            fftU{k} = fftU{k}(:,1);
            fftV{k} = fftV{k}(:,1);
            expVVc{k} = expVVc{k}(1,1);
            expUUc{k} = expUUc{k}(1,1);
            fftSigU{k} = fftSigU{k}(1,1);
            fftSigV{k} = fftSigV{k}(1,1);
        end
    end

    % Check Convergence:           
    Xconv = norm(old_fftX(:) - fftX(:))/norm(old_fftX(:));
    Econv = norm(old_fftE(:) - fftS(:))/norm(old_fftE(:));    
    
    if iI > 1
        if RelChan < tol && Xconv < tol && FitOld > 0.5
            fprintf('Iter %u: Fit = %f, FitTol = %f, Xconv = %g, Econv = %g, R = %s.\n', ...
                iI, Fit(iI), RelChan, Xconv, Econv, mat2str(Rs));
            disp('Evidence Lower Bound Converges.'); 
            break; 
        end
        if mod(iI,10) == 0
            fprintf('Iter %u: Fit = %f, FitTol = %f, Xconv = %g, Econv = %g, R = %s.\n', ...
                iI, Fit(iI), RelChan, Xconv, Econv, mat2str(Rs));
        end        
    else 
        fprintf('Iter %u: Fit = %f, Xconv = %g, Econv = %g, R = %s.\n', ...
            iI, Fit(iI), Xconv, Econv, mat2str(Rs));
    end    
end

model.TX = TX;
model.S = S;
model.Rs = Rs;
model.Fit = Fit(1:iI);
end



