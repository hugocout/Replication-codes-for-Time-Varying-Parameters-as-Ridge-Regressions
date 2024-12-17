function [lambdastar,finalmse,lambda1se,lambda2se] = ...
                                    CVunivariate2(Y,ZZ,k,block_size,lambdavec,lambda2,dimX,plota,sweigths,eweigths,nf)

                                
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% Performs cross-validation to determine the optimal lambda for a
% univariate Ridge regression. The optimal lambda is selected by
% minimizing the mean squared error (MSE) across folds.
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% Y:          Matrix of response data, with dimensions n x m.
%
% ZZ:         Expanded matrix of explanatory variables, with dimensions 
%             [(n-1) x (m x p + 1)   n x (m x p + 1)].
%
% k:          Number of folds for cross-validation.
%
% lambdavec:  Vector of candidate lambda values for Ridge regression.
%
% lambda2:    Penalty for non-u parameters (default: 0.1).
%
% dimX:       Dimension of the original matrix X before expansion.
%
% plota:      If 1, plots are generated; if 0, no plots are generated 
%             (default: 0).
%
% sweigths:   Weights for calculating sigma_u, with default value of 
%             1 for each explanatory variable (default: 1).
%
% eweigths:   Weights for calculating sigma_e, with default value of 
%             1 for each observation (default: 1).
%
% nf:         Number of lambda values, typically equal to dimX (default: 
%             dimX).
%
% Block_size: Number of lags used for cross-validation (default: 8).
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% lambdastar:  Optimal lambda value that minimizes the MSE.
%
% finalmse:    MSE obtained using lambdastar.
%
% lambda1se:   Lambda value at 1 standard error above the optimal lambda.
%
% lambda2se:   Lambda value at 2 standard errors above the optimal lambda.
%
%--------------------------------------------------------------------------
%
% Author:  Christophe Barrette, September 2024
% Based on the function by Goulet Coulombe (2024).
%--------------------------------------------------------------------------

% Set random seed for reproducibility
rng(1071);

%--------------------------------------------------------------------------
% DEFAULT VALUES
%--------------------------------------------------------------------------
% Set default values for input arguments if not provided
if nargin < 11; nf = dimX; end        % Default number of lambdas
if nargin < 10; eweigths = 1; end     % Default weights for sigma_e
if nargin < 9; sweigths = 1; end      % Default weights for sigma_u
if nargin < 8; plot = 0; end          % Default to no plots


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables
T = size(ZZ,1);
data = [Y ZZ];
dfTrain = data;

% Generate cross-validation folds
if block_size == 1
    % Random sampling for each data point
    index = randsample(1:k, size(data,1), true);
else
    % Create a block structure for cross-validation folds
    index = [];
    for jjj = 1:k
        index = [index (ones(1,block_size) * jjj)];
    end
    
    % Extend index to match the number of data points
    www = ceil(size(data,1) / length(index));
    index1 = index;
    for jjj = 2:www
        index = [index index1];
    end
    index = index(1:size(data,1));
end

% Initialize variables for cross-validation
seqY = 1:size(Y,2);
PMSE = zeros(k, length(lambdavec));

% Cross-validation process
for i = 1:k
    % Split data into training and test sets
    testIndexes = find(index == i);
    testData = dfTrain(testIndexes,:);
    trainData = dfTrain(find(index ~= i),:);

    % Set up lambda matrix for regularization
    lambda_t = eye(size(trainData,1));
    if length(eweigths) > 1
        lambda_t = diag(double(eweigths(find(index ~= i))));
    end

    % Handle missing values (dropouts)
    bindex = index;
    bindex(index == i) = 0;
    bindex(bindex ~= 0) = 1;
    Dosigfactor = (1 + cumulzeros(bindex));

    % Compute K matrix for training data
    Z = trainData(:, (size(Y,2) + 1):end);
    ncolZ = size(Z,2);
    X = Z(:, (ncolZ - dimX + 1):ncolZ);
    MX = eye(size(X,1)) - X * inv((X' * X) + lambda2 * eye(size(X,2))) * X';
    MXZ = MX * Z;

    % Adjust weights and apply dropout correction
    for m = 1:nf
        beg = (m - 1) * (T - 1) + 1;
        en = m * (T - 1);
        if length(sweigths) > 1
            MXZ(:, beg:en) = sweigths(m) * MXZ(:, beg:en);
        end
        for tt = 1:(T - 1)
            MXZ(:, (m - 1) * (T - 1) + tt) = MXZ(:, (m - 1) * (T - 1) + tt) .* Dosigfactor(tt + 1);
        end
    end
    Kmat = MXZ * MXZ';

    % Compute K matrix for test data
    z = testData(:, (size(Y,2) + 1):end);
    ncolz = size(z,2);
    x = z(:, (ncolz - dimX + 1):ncolz);
    Mx = eye(size(x,1)) - x * inv((x' * x) + lambda2 * eye(size(x,2))) * x';
    mxz = Mx * z;
    if length(sweigths) == 1
        kmat = mxz * MXZ';
    else
        for m = 1:length(sweigths)
            beg = (m - 1) * (T - 1) + 1;
            en = m * (T - 1);
            mxz(:, beg:en) = sweigths(m) * mxz(:, beg:en);
        end
        kmat = mxz * MXZ';
    end

    % Out-of-sample prediction
    for j = 1:length(lambdavec)
        pred = kmat * ((Kmat + lambdavec(j) * lambda_t) \ (MX * trainData(:, seqY)));
        PMSE(i, j) = mean((pred - Mx * testData(:, 1:size(Y,2))).^2, [1 2]);
    end
end

% Compute the optimal lambda
score = mean(PMSE);
lambdastarpos = find(score == min(score));

if length(lambdastarpos) > 1
    lambdastarpos = lambdastarpos(1);
end

finalmse = score(lambdastarpos);
lambdastar = lambdavec(lambdastarpos);

% Compute standard errors and lambdas at 1 and 2 SE
if length(lambdavec) > 1
    SE = nan(1, length(lambdavec));
    for j = 1:length(lambdavec)
        SE(j) = std(PMSE(:, j)) / sqrt(k);
    end
    se = SE(lambdastarpos); % Standard error at lambdastar
    scoreub = score + SE;
    scorelb = score - SE;

    % Lambda at 1 SE above the minimum
    lambda1sepos = lambdastarpos;
    while true
        if lambda1sepos >= length(score) || score(lambda1sepos) > (finalmse + se)
            break;
        end
        lambda1sepos = lambda1sepos + 1;
    end
    lambda1se = lambdavec(lambda1sepos);

    % Lambda at 2 SE above the minimum
    lambda2sepos = lambdastarpos;
    while true
        if lambda2sepos >= length(score) || score(lambda2sepos) > (finalmse + 2 * se)
            break;
        end
        lambda2sepos = lambda2sepos + 1;
    end
    lambda2se = lambdavec(lambda2sepos);

    % Plot results
    if plota == 1
        limit2 = lambda2sepos;
        while true
            if limit2 >= length(score) || score(limit2) > (finalmse + 20 * se)
                break;
            end
            limit2 = limit2 + 1;
        end
        limit1 = lambdastarpos;
        while true
            if limit1 <= 1 || score(limit1) > (finalmse + 20 * se)
                break;
            end
            limit1 = limit1 - 1;
        end
    
        if finalmse + 2 * se > max(score)
            figure
            h4 = plot(limit1:limit2, score(limit1:limit2));
            hold on
            h6 = xline(9);
            h5 = scatter(lambdastarpos, score(lambdastarpos), 'red');
            title('Cross-Validation Results Variable')
            legend([h4], 'MSE', 'Location', 'northwest')
            hold off
        else
            figure
            h4 = plot(limit1:limit2, score(limit1:limit2));
            hold on
            h6 = xline(9);
            h5 = scatter(lambdastarpos, score(lambdastarpos), 'red');
            h2 = yline(finalmse + se, '-b');
            h3 = yline(finalmse + 2 * se, 'Color', [0.4940 0.1840 0.5560]);
            title('Cross-Validation Results Variable')
            legend([h4, h2, h3], 'MSE', 'Minimum + 1 SE', 'Minimum + 2 SE', 'Location', 'northwest')
            hold off
        end
    end
else
    lambda1se = nan;
    lambda2se = nan;
end
