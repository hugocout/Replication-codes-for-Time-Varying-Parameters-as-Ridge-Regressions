function [lambdastar,finalmse,lambda1se,lambda2se,lambdas_het] = ...
                                    CVmultivariate2(Y,ZZ,k,block_size,lambdavec,lambda2,dimX,plota,sweigths,eweigths,nf,beta0_given);

                                
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This function performs cross-validation to find the optimal lambda for a 
% multivariate Ridge regression model. The lambda parameter controls the 
% amount of regularization applied to the model, which in turn influences 
% the time-variation in the regression coefficients. By testing multiple 
% lambda candidates, the function identifies the best one that minimizes 
% the mean squared error (MSE) of the model.
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% Y:            An n x m matrix containing the response (dependent) variables.
%               Each column represents a different variable in the multivariate
%               regression, and each row corresponds to a time point or observation.
%
% ZZ:           A matrix of explanatory (independent) variables after basis expansion.
%               This matrix has dimensions [(n-1) x (m x p + 1)] or [n x (m x p + 1)],
%               depending on the model setup. It includes lagged values and interactions.
%
% k:            The number of folds for cross-validation. This controls how the data 
%               is split into training and testing sets for model evaluation.
%
% block_size:   The size of the blocks used in cross-validation. This defines the 
%               number of lags to be considered when constructing the training and 
%               testing sets. The default value is 8.
%
% lambdavec:    A vector of candidate lambda values for the Ridge regression. These 
%               values determine the strength of regularization, affecting the amount 
%               of time-variation allowed in the model.
%
% lambda2:      The penalty parameter for non-u parameters (default = 0.1). This value 
%               is used to regularize other aspects of the model beyond the main Ridge 
%               penalty.
%
% dimX:         The dimension of the original matrix of explanatory variables (X) before
%               any basis expansion. This is used to correctly interpret and process the
%               expanded matrix ZZ.
%
% plota:        A flag to control plotting. If set to 1, the function generates graphs 
%               for visualization of the cross-validation process; if set to 0, no graphs 
%               are displayed (default = 0).
%
% sweigths:     Weights used to calculate sigma_u, a vector of size 1 x dimX (default = 1).
%               These weights influence the estimation of the variance of the coefficients.
%
% eweigths:     Weights used to calculate sigma_e, a vector of size 1 x n (default = 1).
%               These weights influence the estimation of the variance of the residuals.
%
% nf:           The number of lambdas to be considered, usually one for each explanatory 
%               variable (default = dimX).
%
% beta0_given:  A logical flag indicating if the initial beta values are known or fixed 
%               (default = FALSE). When set to TRUE, the function uses the provided beta 
%               values during cross-validation.
%
%--------------------------------------------------------------------------
%
% OUTPUTS:
%
% lambdastar:   The optimal lambda value that minimizes the mean squared error (MSE) 
%               during cross-validation.
%
% finalmse:     The MSE obtained using the optimal lambda (lambdastar).
%
% lambda1se:    The lambda value at one standard error from the optimal lambda. This 
%               is often used as a more conservative choice in model selection.
%
% lambda2se:    The lambda value at two standard errors from the optimal lambda.
%
% lambdas_het:  The ratio of sigma_u (coefficient variance) to sigma_e (residual variance) 
%               for each variable. This provides insight into the relative variability of 
%               the coefficients compared to the model errors.
%
%--------------------------------------------------------------------------
%
% Author: Christophe Barrette, September 2024
% Based on the function from Goulet Coulombe (2024)
%--------------------------------------------------------------------------

% Set random seed for reproducibility
rng(1071);

%--------------------------------------------------------------------------
% DEFAULT VALUES
%--------------------------------------------------------------------------
if nargin < 12; beta0_given = false; end  % Default: FALSE
if nargin < 11; nf = dimX; end            % Default: dimX
if nargin < 10; eweigths = 1; end         % Default: 1
if nargin < 9; sweigths = 1; end          % Default: 1
if nargin < 8; plota = 0; end             % Default: 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the number of time periods
T = (size(ZZ,2) - dimX) / nf + 1;

% Adjust ZZ if beta0 is given
if beta0_given
    ZZ = ZZ(:, 1:(size(ZZ,2) - dimX));
end

% Combine Y and ZZ into a single dataset
data = [Y ZZ];
dfTrain = data;

% Create cross-validation indices
if block_size == 1
    index = randsample(1:k, size(data,1), true);
else
    index = [];
    for jjj = 1:k
        index = [index; ones(block_size, 1) * jjj];
    end
    
    % Repeat indices to cover all data points
    num_repeats = ceil(size(data,1) / length(index));
    index = repmat(index, num_repeats, 1);
    index = index(1:size(data,1));
end

% Initialize matrices to store results
seqY = 1:size(Y,2);
PMSE = zeros(k, length(lambdavec));
PMSE_het = zeros(k, length(lambdavec), length(seqY));

% Perform k-fold cross-validation
for i = 1:k
    % Split data into training and testing sets
    testIndexes = find(index == i);
    testData = dfTrain(testIndexes, :);
    trainData = dfTrain(index ~= i, :);

    % Initialize lambda matrix for training
    lambda_t = eye(size(trainData,1));
    if length(eweigths) > 1
        lambda_t = diag(double(eweigths(index ~= i)));
    end

    % Correction for dropouts
    bindex = index;
    bindex(index == i) = 0;
    bindex(bindex ~= 0) = 1;
    Dosigfactor = (1 + cumulzeros(bindex));

    % Prepare training data
    Z = trainData(:, (size(Y,2)+1):end);
    ncolZ = size(Z,2);
    if ~beta0_given
        X = Z(:, (ncolZ - dimX + 1):end);
        MX = eye(size(X,1)) - X * inv((X' * X) + lambda2 * eye(size(X,2))) * X';
        MXZ = MX * Z;
    else
        MXZ = Z;
    end

    % Apply weights and adjust for time-variation
    for m = 1:nf
        beg = (m - 1) * (T - 1) + 1;
        en = m * (T - 1);
        if length(sweigths) > 1
            MXZ(:, beg:en) = sweigths(m) * MXZ(:, beg:en);
        end
        if -nf > -dimX
            MXZ(:, beg) = 1000 * MXZ(:, beg);
        end
        for tt = 1:(T - 1)
            MXZ(:, (m - 1) * (T - 1) + tt) = MXZ(:, (m - 1) * (T - 1) + tt) * Dosigfactor(tt + 1);
        end
    end

    % Determine matrix Kmat
    param = nf * (T - 1);
    obs = size(ZZ,1);
    
    if param > obs
        % Dual problem
        Kmat = MXZ * MXZ';
        
        % Prepare test data
        z = testData(:, (size(Y,2) + 1):end);
        ncolz = size(z,2);
        if ~beta0_given
            x = z(:, (ncolz - dimX + 1):end);
            Mx = eye(size(x,1)) - x * inv((x' * x) + lambda2 * eye(size(x,2))) * x';
            mxz = Mx * z;
        else
            mxz = z;
        end

        % Adjust for time-variation
        for m = 1:nf
            beg = (m - 1) * (T - 1) + 1;
            if -nf > -dimX
                mxz(:, beg) = 1000 * mxz(:, beg);
            end
        end

        % Calculate Kmat for test data
        if length(sweigths) == 1
            kmat = mxz * MXZ';
        else
            for m = 1:nf
                beg = (m - 1) * (T - 1) + 1;
                en = m * (T - 1);
                mxz(:, beg:en) = sweigths(m) * mxz(:, beg:en);
            end
            kmat = mxz * MXZ';
        end

        % Out-of-sample prediction
        if ~beta0_given
            MXY = MX * trainData(:, seqY);
            real = Mx * testData(:, 1:size(Y,2));
        else
            MXY = trainData(:, seqY);
            real = testData(:, 1:size(Y,2));
        end

        for j = 1:length(lambdavec)
            pred = kmat * ((Kmat + lambdavec(j) * lambda_t) \ MXY);
            PMSE(i, j) = mean((pred - real).^2, [1 2]);
            for mm = 1:size(Y,2)
                PMSE_het(i, j, mm) = mean((pred(:, mm) - real(:, mm)).^2);
            end
        end
    
    else
        % Primal problem
        MXZ2 = MXZ;
        if length(eweigths) > 1
            for tt = 1:size(lambda_t,1)
                MXZ2(tt,:) = MXZ(tt,:) * (lambda_t(tt,tt) ^ -1);
            end
        end
        Kmat = MXZ2' * MXZ;
        
        % Prepare test data
        z = testData(:, (size(Y,2) + 1):end);
        ncolz = size(z,2);
        if ~beta0_given
            x = z(:, (ncolz - dimX + 1):end);
            Mx = eye(size(x,1)) - x * inv((x' * x) + lambda2 * eye(size(x,2))) * x';
            mxz = Mx * z;
        else
            mxz = z;
        end

        % Adjust for time-variation
        for m = 1:nf
            beg = (m - 1) * (T - 1) + 1;
            if -nf > -dimX
                mxz(:, beg) = 1000 * mxz(:, beg);
            end
        end

        % Apply weights if needed
        if length(sweigths) > 1
            for m = 1:nf
                beg = (m - 1) * (T - 1) + 1;
                en = m * (T - 1);
                mxz(:, beg:en) = sweigths(m) * mxz(:, beg:en);
            end
        end

        % Out-of-sample prediction
        if ~beta0_given
            MXY = MXZ2' * (MX * trainData(:, seqY));
            real = MX * testData(:, 1:size(Y,2));
        else
            MXY = MXZ2' * trainData(:, seqY);
            real = testData(:, 1:size(Y,2));
        end

        for j = 1:length(lambdavec)
            pred = mxz * ((Kmat + lambdavec(j) * eye(size(Kmat,1))) \ MXY);
            PMSE(i, j) = mean((pred - real).^2, [1 2]);
            for mm = 1:size(Y,2)
                PMSE_het(i, j, mm) = mean((pred(:, mm) - real(:, mm)).^2);
            end
        end
    end
end

% Compute mean MSE for each lambda
score_het = nan(size(PMSE_het,2), size(PMSE_het,3));
lambdas_het = seqY;
for mm = 1:size(PMSE_het,3)                
    score_het(:,mm) = mean(PMSE_het(:,:,mm))';    
end

% Find optimal lambda for each variable
for mm = 1:size(Y,2)
    a = lambdavec(score_het(:,mm) == min(score_het(:,mm)));
    if length(a) == 1
        lambdas_het(mm) = a;
    else
        lambdas_het(mm) = a(1,1);
    end
end

% Compute overall score and find optimal lambda
score = nan(size(PMSE,2), size(PMSE,3));
for mm = 1:size(PMSE,3)            
    for mmm = 1:size(PMSE,2)
        score(mmm,mm) = mean(PMSE(:,mmm,mm));
    end
end

% Find position of the optimal lambda
lambdastarpos = find(score == min(score));
if length(lambdastarpos) > 1
    lambdastarpos = lambdastarpos(1);
end

finalmse = score(lambdastarpos);
lambdastar = lambdavec(lambdastarpos);

% Apply one standard error rule if more than one lambda
if length(lambdavec) > 1
    SE = nan(1,length(lambdavec));
    for j = 1:length(lambdavec)
        SE(j) = std(PMSE(:,j)) / sqrt(k);
    end
    se = SE(lambdastarpos);
    scoreub = score' + double(SE);
    scorelb = score' - double(SE);

    % Find lambda values at one and two standard errors
    lambda1sepos = lambdastarpos;
    while lambda1sepos < length(score) && score(lambda1sepos) <= (finalmse + se)
        lambda1sepos = lambda1sepos + 1;
    end
    lambda1se = lambdavec(lambda1sepos);

    lambda2sepos = lambdastarpos;
    while lambda2sepos < length(score) && score(lambda2sepos) <= (finalmse + 2 * se)
        lambda2sepos = lambda2sepos + 1;
    end
    lambda2se = lambdavec(lambda2sepos);

    % Plot results if requested
    if plota == 1
        limit2 = lambda2sepos;
        while limit2 < length(score) && score(limit2) <= (finalmse + 20 * se)
            limit2 = limit2 + 1;
        end
        limit1 = lambdastarpos;
        while limit1 > 1 && score(limit1) <= (finalmse + 20 * se)
            limit1 = limit1 - 1;
        end

        % Create plot
        figure;
        if finalmse + 2 * se > max(score)
            h4 = plot(limit1:limit2, score(limit1:limit2));
            hold on;
            h5 = scatter(lambdastarpos, score(lambdastarpos), 'red');
            title('Cross-Validation Results Multivariate');
            legend([h4], 'MSE', 'Location', 'northwest');
            hold off;
        else
            h4 = plot(limit1:limit2, score(limit1:limit2));
            hold on;
            h5 = scatter(lambdastarpos, score(lambdastarpos), 'red');
            h2 = yline(finalmse + se, '-b'); % finalmse + se
            h3 = yline(finalmse + 2 * se, 'Color', [0.4940 0.1840 0.5560]); % finalmse + 2 * se
            title('Cross-Validation Results Multivariate');
            legend([h4, h2, h3], 'MSE', 'Minimum + 1 SE', 'Minimum + 2 SE', 'Location', 'northwest');
            hold off;
        end
    end
else
    lambda1se = nan;
    lambda2se = nan;
end

