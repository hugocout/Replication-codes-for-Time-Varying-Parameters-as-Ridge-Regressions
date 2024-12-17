
function [varargout] = IRF_tvp_2SRR(Bcomp, RES, IRhoriz,shock,cumulative,phi);


%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% Computes the impulse response functions (IRFs) based on a
% given reduced-form VAR model. The impulse responses are calculated for
% a specified number of periods, showing how shocks to one variable
% propagate through the system over time.
%
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% Computes and plots the impulse response functions (IRFs) for a specified
% VAR model. This function generates the impulse responses based on a given
% shock to a specific variable and calculates cumulative impulse responses 
% if specified. The IRFs are computed over a defined number of periods.
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% Bcomp:         Matrix containing the estimated reduced-form coefficients 
%                in companion form, excluding deterministic terms. This 
%                matrix has dimensions (n x vlag) x (n x vlag) x T, where 
%                n is the number of variables and vlag is the number of lags. 
%                It includes both the estimated parameters and a diagonal 
%                identity matrix that represents the lead-lag relationships 
%                between variables.
%
% RES:           Matrix of residuals from the VAR model, with dimensions T x n. 
%                This matrix captures the deviations from the model’s fitted values.
%
% IRhoriz:       Number of periods for which the impulse responses are 
%                computed and plotted. Defines the horizon of the IRF analysis.
%
% shock:         Scalar specifying the position of the variable responsible 
%                for the shock. Determines which variable in the system will 
%                be subjected to the shock.
%
% cumulative:    Scalar specifying the variable for which the cumulative 
%                impulse response should be computed. If set, the IRF for this 
%                variable will be accumulated over the periods.
%
% phi:           Parameter controlling time variation in the impulse responses.
%                This parameter can be used to account for changes in the impact 
%                of shocks over time.
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% IR:            Matrix of impulse responses with dimensions (IRhoriz x n x n). 
%                This matrix contains the responses of all variables to a shock 
%                for each period up to the specified horizon. It provides insights 
%                into how shocks propagate through the system over time.
%
%--------------------------------------------------------------------------
%
% Author:        Christophe Barrette, September 2024
% Based on:      Goulet Coulombe (2024)
%--------------------------------------------------------------------------

% Set random seed for reproducibility
rng(1071);

%--------------------------------------------------------------------------
% DEFAULT VALUES
%--------------------------------------------------------------------------
if nargin<6; phi=1000; end % Default value for 'phi'

%--------------------------------------------------------------------------
% TURN OFF WARNINGS
%--------------------------------------------------------------------------
warning('off', 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_ini=phi;

while true
    try
        % Initialize variables
        numObs = size(RES, 1);
        numVars = size(RES, 2);
        
        % Compute the lower triangular covariance matrices
        VCRESlow1 = zeros(numVars, numVars, numObs);
        for k = 1:numObs
            VCRESlow1(:,:,k) = tril(RES(k:k, :)' * RES(k:k, :));
        end
        VCRESlow = VCRESlow1;

        % Flatten the lower triangular matrices into vectors
        VCRESveclow1 = [];
        for k = 1:size(VCRESlow1, 3)
            VCmid = VCRESlow(:,:,k);
            VCRESveclow1 = [VCRESveclow1; VCmid(tril(true(size(VCmid))))'];
        end

        % Center the vectors by subtracting the mean
        meanVCRESveclow1 = mean(VCRESveclow1);
        VCRESveclow = VCRESveclow1 - meanVCRESveclow1;

        % Extend the data by duplicating the first and last observations
        VCRESveclow = [repmat(VCRESveclow(1, :), 30, 1); VCRESveclow; repmat(VCRESveclow(end, :), 30, 1)];

        % Create the difference matrix for regularization
        D = full(spdiags([ones(size(VCRESveclow, 1), 1), -ones(size(VCRESveclow, 1), 1)], [0, 1], size(VCRESveclow, 1), size(VCRESveclow, 1)))';

        % Regularize the covariance matrix
        VCRESvec1 = (eye(size(VCRESveclow, 1)) + (phi * D' * D)) \ VCRESveclow;

        % Remove the duplicated observations
        VCRESvec2 = VCRESvec1(31:end-30, :);

        % Restore the mean to the covariance matrix
        VCRESvec = VCRESvec2 + meanVCRESveclow1;

        % Initialize the matrix to store the covariance matrices
        VC_eps1 = zeros(size(VCRESlow));
        for k = 1:size(VCRESlow1, 3)
            L = 1;
            T = size(RES, 2);
            for i = 1:size(RES, 2)
                VC_eps1(i:end, i, k) = VCRESvec(k, L:T)';
                VC_eps1(i, i:end, k) = VCRESvec(k, L:T);
                L = L + size(RES, 2) - i + 1;
                T = T + size(RES, 2) - i;
            end
        end
        VC_eps = VC_eps1;

        % Initialize variables for Cholesky decomposition
        Nbig = size(VC_eps, 2);
        Nbigcomp = max(size(Bcomp, 1), size(Bcomp, 2));
        periodV = size(VC_eps, 3);
        VC_epschol = zeros(Nbig, Nbig, periodV);
        stddepsvec = zeros(Nbig, 1, periodV);

        % Compute Cholesky decompositions
        for k = 1:periodV
            VC_epschol(:,:,k) = chol(VC_eps(:,:,k))';
            stddepsvec(:,:,k) = diag(VC_eps(:,:,k)).^0.5;
        end

        % If Cholesky decomposition succeeds, exit the loop
        break;
    catch
        % If Cholesky decomposition fails, increment phi and retry
        phi = phi + 100;
        if phi > 10000
            fprintf('The Cholesky decomposition is not possible with a value of %d.\n', phi);
            break;
        end
    end
end

% Report the minimum or exact value of phi for which Cholesky decomposition is possible
if phi > phi_ini
    if phi < 10000
        fprintf('The minimum value of Phi for which Cholesky decomposition is possible is %d.\n', phi);
    end
end

if phi ~= phi_ini
fprintf('The value of Phi for which Cholesky decomposition is possible is %d.\n', phi);
end

% Prepare output matrices for impulse response calculations
shockvecmat = VC_epschol;

[nz, nofIR] = size(shockvecmat);
aux = zeros(Nbigcomp, Nbigcomp, periodV);
aux(1:Nbig, 1:Nbig, 1:periodV) = shockvecmat;
shockvecmat = aux;
period = size(Bcomp, 3);

% Initialize matrix for impulse responses
IR = zeros(IRhoriz + 1, Nbig, 1, period);
Impmat = eye(Nbigcomp);

% Calculate impulse responses
if period == 1
    for hh = 1:IRhoriz + 1
        IRbig = Impmat * shockvecmat;
        IR(hh, 1:Nbig, 1:Nbig) = IRbig(1:Nbig, 1:Nbig);
        Impmat = Impmat * Bcomp;
    end
else
    for t = 1:period
        Impmat = eye(Nbigcomp);
        for hh = 1:IRhoriz + 1
            IRbig = Impmat * shockvecmat(:,:,min(t + hh, size(shockvecmat, 3)));
            IR(hh, 1:Nbig, 1, t) = IRbig(1:Nbig, shock);
            Impmat = Impmat * Bcomp(:,:,min(t + hh, size(Bcomp, 3)));
        end
    end
end

% Calculate cumulative impulse responses if specified
for i = 1:size(RES, 2)
    eval(sprintf('IR%d = nan(size(RES, 1), size(IR, 1));', i));
end

kk = 1:size(RES, 2);

if ~isempty(cumulative)
    i = cumulative;
    indices = ismember(kk, i);
    j = kk(~indices);

    for n = 1:size(RES, 1)
        for k = 1:size(RES, 2)
            eval(sprintf('IR%d(n, 1) = IR(1, k, 1, n);', k));
        end
        for k = 2:IRhoriz + 1
            for ii = 1:length(i)
                iii = i(ii);
                eval(sprintf('IR%d(n, k) = IR%d(n, k - 1) + IR(k, iii, 1, n);', iii, iii));
            end
            for jj = 1:length(j)
                jjj = j(jj);
                eval(sprintf('IR%d(n, k) = IR(k, jjj, 1, n);', jjj));
            end
        end
    end
else
    for n = 1:size(RES, 1)
        for j = 1:IRhoriz + 1
            for i = 1:size(RES, 2)
                eval(sprintf('IR%d(n, j) = IR(j, k(i), 1, n);', kk(i)));
            end
        end
    end
end

% Prepare output arguments
for i = 1:size(RES, 2)
    eval(sprintf('varargout{i} = IR%d;', i));
end
