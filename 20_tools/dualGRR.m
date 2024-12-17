function [uhat, betas_grr,yhat,betas_grr_ci,lambdas,sigmasq] = ...
                      dualGRR(Zprime,y, dimX, lambda1,lambda2,sweigths,olsprior,ridgeprior,CI,eweigths,GCV,calcul_beta,nf); 

%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This function estimates the parameters of a dual Ridge regression model,
% including the coefficients and variance components. It handles various 
% configurations such as co-integration, generalized cross-validation, and 
% different weighting schemes.
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% Zprime:   Expanded matrix of explanatory variables, formatted as 
%           [(n-1) x (m x p + 1)   n x (m x p + 1)].
%
% Y:        Matrix of response variables, with dimensions n x m.
%
% dimX:     The dimension of the matrix X, representing the original ZZ matrix
%           before expansion.
%
% lambda1:  Vector representing the ratio of sigma_u to sigma_epsilon for each 
%           variable.
%
% lambda2:  Regularization parameter for non-u parameters (default = 0.001).
%
% sweigths: Weights used for calculating sigma_u, specified as a 1 x dimX vector
%           (default = 1).
%
% olsprior: A flag indicating whether to shrink beta_0 towards OLS estimates 
%           instead of zero (default = 0).
%
% CI:       A flag indicating whether co-integration is present (default = 0).
%
% eweigths: Weights used for calculating sigma_e, specified as a 1 x n vector 
%           (default = 1).
%
% GCV:      A flag indicating whether to perform Generalized Cross-Validation 
%           (default = 0).
%
% calcul_beta: Flag indicating whether to calculate the beta coefficients 
%              (default = true).
%
% nf:       Number of lambda values, typically one for each explanatory variable 
%           (default = dimX).
%
%--------------------------------------------------------------------------
%
% OUTPUTS:
%
% uhat:     Estimated sigma_u values.
%
% betas_grr: Matrix of estimated beta coefficients for the Ridge regression.
%
% yhat:     Estimated Y matrix.
%
% betas_grr_ci: Matrix of estimated beta coefficients for the Ridge regression 
%               accounting for co-integration.
%
% lambdas:  Vector of the ratio of sigma_u to sigma_epsilon for each variable.
%
% sigmasq:  Vector of the estimated variances, weighted by sweigths.
%
%--------------------------------------------------------------------------
%
% Author:  Christophe Barrette, September 2024
% Based on the function developed by Goulet Coulombe (2024).
%--------------------------------------------------------------------------

% Set random seed for reproducibility
rng(1071);

%--------------------------------------------------------------------------
% DEFAULT VALUES
%--------------------------------------------------------------------------
if nargin < 13; nf = 0; end          % Number of lambdas
if nargin < 12; calcul_beta = 1; end % Calculate beta coefficients (true)
if nargin < 11; GCV = 0; end         % Generalized Cross-Validation (false)
if nargin < 10; eweigths = 1; end    % Weights for sigma_e
if nargin < 9; CI = 0; end           % Co-integration (false)
if nargin < 8; ridgeprior = 0; end   % Ridge prior (false)
if nargin < 7; olsprior = 0; end     % OLS prior (false)
if nargin < 6; sweigths = 1; end    % Weights for sigma_u
if nargin < 5; lambda2 = 0.001; end  % Regularization parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% PREPARE MATRICES AND SET DEFAULT VALUES
%--------------------------------------------------------------------------
% Initialize lambdas and set default number of lambdas if not specified
lambdas = [lambda1 lambda2];
if nf == 0
    nf = dimX;
end

% Ensure lambda2 is not too small
if lambda2 < 1e-18
    fprintf('Lambda_2 imposed to be at least 1e-18\n');
    lambdas(2) = 1e-18;
end

% Compute size parameters
ncolZ = size(Zprime, 2);
T = (ncolZ - dimX) / nf + 1;

% Set default weights if not provided
if length(sweigths) == 1
    sweigths = ones(1, nf);
end

% Adjust Kmat_half matrix
Kmat_half = Zprime';
for m = 1:nf
    beg = (m - 1) * (T - 1) + 1;
    en = m * (T - 1);
    Kmat_half(beg:en, :) = (1 / lambdas(1) * sweigths(m)) * Kmat_half(beg:en, :);
    if nf > dimX
        Kmat_half(beg, :) = 1000 * Kmat_half(beg, :);
    end
end
Kmat_half((ncolZ - dimX + 1):ncolZ, :) = (1 / lambdas(2)) * Kmat_half((ncolZ - dimX + 1):ncolZ, :);

%--------------------------------------------------------------------------
% DUAL OR PRIMAL GRR ESTIMATION
%--------------------------------------------------------------------------
% Initialize Lambda_T matrix
Lambda_T = eye(size(Zprime, 1));
if length(eweigths) > 1
    Lambda_T = diag(eweigths);
end

% Number of parameters and observations
param = nf * (T - 1) + dimX;
obs = size(Zprime, 1);

if param > obs % Dual GRR estimation
    Kmat_du = Zprime * Kmat_half;

    if ridgeprior == 0
        if olsprior == 0
            alpha = (Kmat_du + Lambda_T) \ y;
            uhat = Kmat_half * alpha;
        elseif olsprior == 1
            X = Zprime(:, (ncolZ - dimX + 1):ncolZ);
            beta_ols = (X' * X) \ (X' * y);
            alpha = (Kmat_du + Lambda_T) \ (y - X * beta_ols);
            uhat = Kmat_half * alpha;
            uhat((ncolZ - dimX + 1):ncolZ) = uhat((ncolZ - dimX + 1):ncolZ) + beta_ols;
        end
    elseif ridgeprior == 1
        X = Zprime(:, (ncolZ - dimX + 1):ncolZ);
        beta_ols = (X' * X + lambda2 * eye(size(X, 2))) \ (X' * y); 
        alpha = (Kmat_du + Lambda_T) \ (y - X * beta_ols);
        uhat = Kmat_half * alpha;
        uhat((ncolZ - dimX + 1):ncolZ) = uhat((ncolZ - dimX + 1):ncolZ) + beta_ols;
    end
    yhat = Kmat_du * alpha;

else % Primal GRR estimation
    % Adjust weights for Kmat_half
    for tt = 1:obs
        Kmat_half(:, tt) = Kmat_half(:, tt) * (eweigths(tt) ^ -1);
    end
    
    Kmat_pri = Kmat_half * Zprime;

    if ridgeprior == 0
        if olsprior == 0
            uhat = (Kmat_pri + diag(param)) \ (Kmat_half * y);
        elseif olsprior == 1
            X = Zprime(:, (ncolZ - dimX + 1):ncolZ);
            beta_ols = (X' * X) \ (X' * y);
            uhat = (Kmat_pri + eye(param)) \ (Kmat_pri * (y - X * beta_ols));
            uhat((ncolZ - dimX + 1):ncolZ) = uhat((ncolZ - dimX + 1):ncolZ) + beta_ols;
        end
    elseif ridgeprior == 1
        X = Zprime(:, (ncolZ - dimX + 1):ncolZ);
        beta_ols = (X' * X + lambda2 * eye(size(X, 2))) \ (X' * y); 
        uhat = (Kmat_pri + eye(param)) \ (Kmat_pri * (y - X * beta_ols));
        uhat((ncolZ - dimX + 1):ncolZ) = uhat((ncolZ - dimX + 1):ncolZ) + beta_ols;
    end
    yhat = Zprime * uhat;
end

%--------------------------------------------------------------------------
% RECOVER BETAS
%--------------------------------------------------------------------------
if calcul_beta
    betas_grr = zeros(size(y, 2), nf, T);
    for eq = 1:size(y, 2)
        for k = 1:nf
            betas_grr(eq, k, 1) = uhat((size(uhat, 1) - dimX + k), eq);
        end
        for t = 2:T
            positions = [];
            for k = 1:nf
                positions = [positions, (k - 1) * (T - 1) + (t - 1)];
            end
            betas_grr(eq, :, t) = betas_grr(eq, :, t - 1) + uhat(positions, eq)';
        end
    end
end

%--------------------------------------------------------------------------
% CO-INTEGRATION (CI)
%--------------------------------------------------------------------------
if CI > 0
    Kmh = Kmat_half;
    if param > obs % Dual GRR for CI
        for tt = 1:obs
            Kmh(:, tt) = Kmat_half(:, tt) * (eweigths(tt) ^ -1);
        end
    end

    if dimX > nf
        Kmh = Kmh(1:(size(Kmat_half, 1) - dimX), :);
    else
        % Reorder Kmat_half
        for k = 1:dimX
            Kmh((1 + T * (k - 1)), :) = Kmh((size(Kmat_half, 1) - dimX + k), :);
            beg = 2 + T * (k - 1);
            en = T + T * (k - 1);
            Kmh(beg:en, :) = Kmh((beg - k):(en - k), :);
        end
    end

    % Compute variance-covariance matrix
    inCov = inv(chol(Kmh * Kmh' + eye(size(Kmh, 1)))' * chol(Kmh * Kmh' + eye(size(Kmh, 1))));
    CT = tril(ones(T - 1, T - 1));
    C = kron(eye(nf), CT);
    sandwich = inCov * (std(y - yhat) ^ 2) * lambdas(1:end-1) ^ -1;
    Vb = sqrt(diag(C * sandwich * C'));

    betas_grr_ci = zeros(size(y, 2), nf, T, 2);
    for c = 1:2
        for eq = 1:size(y, 2)
            for k = 1:nf
                for t = 1:T
                    betas_grr_ci(eq, k, t, c) = betas_grr(eq, k, t) + ((-1) ^ c * CI * Vb((k - 1) * (T - 1) + 1));
                end
            end
        end
    end
else
    betas_grr_ci = nan(1, 1);
end

% Return sigma^2 weights
sigmasq = sweigths;
