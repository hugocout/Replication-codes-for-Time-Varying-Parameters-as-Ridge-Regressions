function [BETAS_GRR, BETAS_VARF, Bcomp, LAMBDAS, fcast, YHAT_VAR, YHAT_VARF,EWvec] = ...
                      tvp_2SRR(X,Y,block_size, lambdacandidates, oosX, lambda2, kfold, CVplot, CV2SRR, siguparam, sigepsparam, olsprior,ridgeprior)

                                
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This routine estimates the residuals and companion matrix 
% of a reduced-form VAR using Ridge regression and the two-step Ridge regression 
% proposed in Goulet Coulombe (2024). 
% 
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% X:    Matrix of explanatory data series, n x (m x p).
%
% Y:    Matrix of data series, n x m.    
%
% lambdacandidates: Candidates for Ridge's lambda driving the amount of
% time-variation; default included.
%
% oosX:  A vector of X for out-of-sample prediction (default = empty).
%
% lambda2:  The penalty for non-‘u’ parameters (default = 0.1).
%
% kfold:  Number of folds for cross-validation (default = 5).
%
% CVplot:  If 1, some graphs are displayed; if 0, they are not (default = FALSE).
%
% CV2SRR:  When doing multivariate analysis, should the second cross-validation step 
% be performed for each equation separately? (This can take more time; default = TRUE).
%
% siguparam:  How much should the constant variance be shrunk through time?
%   (Between 0 and 1; 0 means constant variance, 1 means plain two-step; default = 0.75).
%
% sigepsparam:  How much should the constant variance be shrunk across ‘u’?
%   (Between 0 and 1; 0 means constant variance, 1 means plain two-step; default = 0.75).
%
% olsprior:  If 1, then shrink beta_0's to OLS rather than to 0 (default = 0).
%
% Block_size:  Lags for cross-validation (default = 8).
%--------------------------------------------------------------------------
%
% OUTPUTS:
%
% BETAS_GRR:  Matrix with estimated betas of the Ridge regression.
%
% BETAS_VARF:  Matrix with estimated betas of the two-step Ridge regression.
%
% Bcomp:  Matrix with the structure of estimated reduced-form coefficients in 
%         companion form, excluding all deterministic terms. This includes the 
%         estimated parameters in the partition n x (n x (n x vlag)). The remaining 
%         partition is a diagonal identity matrix that defines the lead-lag relationship 
%         between variables. The size is (n x vlag) x (n x vlag).  
%
% LAMBDAS:  The ratio of sigma_u to sigma_epsilon for each variable.
%
% fcast:  Forecast value from the out-of-sample prediction.
%
% YHAT_VAR:  Predicted Y values (YHAT) from the Ridge regression.
%
% YHAT_VARF:  Predicted Y values (YHAT) from the two-step Ridge regression.
%
% EWvec:  Normalized evolving volatility weights.
%
%--------------------------------------------------------------------------
%
% Author: Christophe Barrette, September 2024
% Based on the function in Goulet Coulombe (2024).
%--------------------------------------------------------------------------


% Set random seed for reproducibility
rng(1071);

%--------------------------------------------------------------------------
% DEFAULT VALUES
%--------------------------------------------------------------------------
if nargin < 13; ridgeprior = 1; end  % Default: ridgeprior similar to olsprior
if nargin < 12; olsprior = 0; end    % Default: OLS prior is disabled
if nargin < 11; sigepsparam = 0.75; end  % Default: 0.75 for sigma_epsilon shrinkage
if nargin < 10; siguparam = 0.75; end    % Default: 0.75 for sigma_u shrinkage
if nargin < 9; CV2SRR = true; end    % Default: Perform second CV step for each equation
if nargin < 8; CVplot = false; end   % Default: Disable cross-validation plots
if nargin < 7; kfold = 10; end        % Default: 5 folds for cross-validation
if nargin < 6; lambda2 = 1; end      % Default: Choose lambda2 by ridge regression
if nargin < 5; oosX = []; end        % Default: No out-of-sample X data
if nargin < 4; lambdacandidates = exp(linspace(-5, 10, 20)); end % Default lambda candidates
if nargin < 3; block_size = 8; end  % Default block size for cross-validation

%--------------------------------------------------------------------------
% TURN OFF WARNINGS
%--------------------------------------------------------------------------
warning('off', 'all');

%--------------------------------------------------------------------------
% Variables re-naming
%--------------------------------------------------------------------------
lambdavec= lambdacandidates;
CVagain=CV2SRR;
sweights=1;
svparam=sigepsparam;
homoparam=siguparam;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% INITIALIZE VARIABLES AND PARAMETERS
%--------------------------------------------------------------------------
Yoriginal = Y; % Store the original Y matrix
[bigT, M] = size(Y); % Get the dimensions of Y

scalingfactor = ones(M, size(X, 2)); % Initialize scaling factors
sdy = 1:M; % Placeholder for standard deviations of Y

%--------------------------------------------------------------------------
% SCALING AND NORMALIZING DATA
%--------------------------------------------------------------------------
% Normalize each column of X and calculate scaling factors
for j = 1:size(X, 2)
    for m = 1:M
        scalingfactor(m, j) = std(Y(:, m)) / std(X(:, j));
    end
    X(:, j) = X(:, j) / std(X(:, j));
end

% Normalize each column of Y and store their standard deviations
for m = 1:M
    sdy(m) = std(Y(:, m));
    Y(:, m) = Y(:, m) / sdy(m);
end

ZZ = Zfun(X); % Basis expansion
YY = Y; % Univariate model
yy = reshape(YY, [], 1); % Flatten YY into a vector
dimX = size(X, 2) + 1; % Include a constant term

%--------------------------------------------------------------------------
% INITIALIZE OUTPUTS AND SETTINGS
%--------------------------------------------------------------------------
BETAS_GRR = zeros(M, dimX, size(Y, 1)); % Coefficients for Ridge regression
BETAS_GRRATS = BETAS_GRR; % Placeholder for 2-step Ridge regression coefficients
LAMBDAS = zeros(1, M); % Lambda values for each equation
EWmat = zeros(size(Y, 1), M); % Evolving volatility weights
YHAT_VAR = Y; % Placeholder for Ridge regression predictions
YHAT_VARF = Y; % Placeholder for 2-step Ridge regression predictions

%--------------------------------------------------------------------------
% CROSS-VALIDATION TO FIND BEST LAMBDA
%--------------------------------------------------------------------------
if lambda2 == 1
    data = [Y X]; % Combine Y and X for training
    dfTrain = data;

    % Create indices for k-fold cross-validation
    if block_size == 1
        index = randsample(1:kfold, size(data, 1), true);
    else
        index = [];
        for jjj = 1:kfold
            index = [index (ones(1, block_size) * jjj)];
        end

        www = ceil(size(data, 1) / length(index));
        index1 = index;

        for jjj = 2:www
            index = [index index1];
        end

        index = index(1:size(data, 1));
    end

    % Initialize MSE matrix for cross-validation
    Mse = nan(kfold, length(lambdavec));

    % Perform k-fold cross-validation
    for i = 1:kfold
        testIndexes = find(index == i);
        testData = dfTrain(testIndexes, :);
        trainData = dfTrain(find(index ~= i), :);

        for ii = 1:length(lambdavec)
            lambda123 = lambdavec(ii);

            % Ridge regression coefficients
            B = (trainData(:, size(Y, 2) + 1:end)' * trainData(:, size(Y, 2) + 1:end) + ...
                 (eye(size(X, 2)) * lambda123)) \ (trainData(:, size(Y, 2) + 1:end)' * trainData(:, 1:size(Y, 2)));
            y_pred = testData(:, size(Y, 2) + 1:end) * B;
            Mse(i, ii) = mean((testData(:, 1:size(Y, 2)) - y_pred).^2, "all");
        end
    end

    % Calculate mean MSE across folds for each lambda
    meanMSE = mean(Mse, 1);

    % Find the best lambda with the minimum mean MSE
    [bestMSE, bestLambdaIdx] = min(meanMSE);
    lambda2 = lambdavec(bestLambdaIdx); % Set optimal lambda2
else
    lambda2 = lambda2;
end

%--------------------------------------------------------------------------
% SECOND CROSS-VALIDATION STEP (OPTIONAL)
%--------------------------------------------------------------------------
if CV2SRR == true
    if length(lambdavec) > 1
        lsubset = 1:length(lambdavec);
        [cvlist_lambdastar, cvlist_finalmse, cvlist_lambda1se, cvlist_lambda2se, cvlist_lambdas_het] = ...
            CVmultivariate2(Y, ZZ, kfold, block_size, lambdavec(lsubset), lambda2, dimX, CVplot, sweights);
        lambdas_list = cvlist_lambdas_het;
    end
else
    lambdas_list = lambdavec * ones(1, M);
end

%--------------------------------------------------------------------------
% ESTIMATION LOOP FOR EACH EQUATION
%--------------------------------------------------------------------------
for m = 1:M
    if M > 1
        disp(['Computing/tuning equation: ', num2str(m)]);
    end

    lambda1 = lambdas_list(m);

    % Final estimation using dualGRR
    [grr_uhat, grr_betas_grr, grr_yhat, grr_betas_grr_ci, grr_lambdas, grr_sigmasq] = ...
        dualGRR(ZZ, Y(:, m), dimX, lambda1, lambda2, sweights, olsprior, ridgeprior);

    BETAS_GRR(m, :, :) = grr_betas_grr(1, :, :);
    LAMBDAS(m) = lambda1;
    YHAT_VAR(:, m) = grr_yhat;

    % Handle outliers
    e = Y(:, m) - grr_yhat;
    e(abs(e) > 50 * mad(e)) = 0;

    % Estimate evolving volatility
    mdl = garch('GARCHLags', 1, 'ARCHLags', 1, 'Offset', NaN);
    estmdl = estimate(mdl, e, 'Display', 'off');
    sestmdl1 = infer(estmdl, e);
    sestmdl = sqrt(sestmdl1);
    EW = sestmdl.^svparam;
    EW = (EW / mean(EW));
    EWmat(:, m) = EW;

    % Rescale coefficients
    betas_grr = reshape(BETAS_GRR(m, :, :), [size(BETAS_GRR, 2), size(BETAS_GRR, 3)]);
    umat = betas_grr(:, 2:end) - betas_grr(:, 1:end-1);
    sigmasq = diag(umat * umat').^homoparam;
    sigmasq = (sigmasq / mean(sigmasq));

    % Optional second cross-validation
    if CVagain
        [cvlist_lambdastar, cvlist_finalmse, cvlist_lambda1se, cvlist_lambda2se] = ...
            CVunivariate2(Y(:, m), ZZ, kfold, block_size, lambdavec, lambda2, dimX, CVplot, sigmasq, EW);
        usethislambda = cvlist_lambdastar;
    else
        usethislambda = lambda1;
    end

    % Final estimation
    if homoparam > 0 || svparam > 0
        [grrats_uhat, grrats_betas_grr, grrats_yhat, vbetas_grr_ci, grrats_lambdas, grrats_sigmasq] = ...
            dualGRR(ZZ, Y(:, m), dimX, usethislambda, lambda2, sigmasq, olsprior, ridgeprior, 0, EW);
        betas_grrats = grrats_betas_grr;
        BETAS_GRRATS(m, :, :) = grrats_betas_grr(1, :, :);
        YHAT_VARF(:, m) = grrats_yhat;
    else
        betas_grrats = grr_betas_grr;
        BETAS_GRRATS(m, :, :) = grr_betas_grr(1, :, :);
        YHAT_VARF(:, m) = YHAT_VAR(:, m);
    end
end

EWvec = reshape(EWmat, [], 1) / mean(EWmat, [1 2]); % Normalize evolving volatility weights
lambda1_step2 = lambda1; % Store lambda1 for second step
fcast = []; % Placeholder for forecast output
BETAS_VARF = BETAS_GRRATS; % Set final coefficients for the 2-step Ridge regression
BETAS_VARF_STD = BETAS_VARF; % Placeholder for standard errors (if needed)

%--------------------------------------------------------------------------
% RESCALE BETAS AND PREDICTIONS TO ORIGINAL SCALE
%--------------------------------------------------------------------------
for m = 1:M
    for j = 1:(dimX - 1)
        BETAS_VARF(m, j + 1, :) = BETAS_VARF(m, j + 1, :) .* scalingfactor(m, j);
        BETAS_GRR(m, j + 1, :) = BETAS_GRR(m, j + 1, :) .* scalingfactor(m, j);
    end

    BETAS_VARF(m, 1, :) = BETAS_VARF(m, 1, :) .* sdy(m);
    BETAS_GRR(m, 1, :) = BETAS_GRR(m, 1, :) .* sdy(m);
end

% Rescale predictions
YHAT_VARF = YHAT_VARF .* sdy;
YHAT_VAR = YHAT_VAR .* sdy;

%--------------------------------------------------------------------------
% OUT-OF-SAMPLE FORECASTING
%--------------------------------------------------------------------------
if length(oosX) > 1
    fcast = BETAS_VARF(:, :, size(BETAS_VARF, 3)) * double([ones(size(oosX, 1), 1), oosX])';
end

%--------------------------------------------------------------------------
% CONSTRUCT Bcomp MATRIX
%--------------------------------------------------------------------------
lag = size(X, 2) / M;
Bcomp = nan(size(X, 2), size(X, 2), size(X, 1));
for n = 1:size(BETAS_VARF, 3)
    Bcomp(:, :, n) = [BETAS_VARF(:, 2:end, n); kron(eye(lag - 1), eye(size(Y, 2))), zeros((lag - 1) * size(Y, 2), size(Y, 2))];
end

%--------------------------------------------------------------------------
% FINAL SETTINGS
%--------------------------------------------------------------------------

if CVplot == false
    close all;
end




