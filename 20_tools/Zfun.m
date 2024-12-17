function [Zprime] = Zfun(data);

                                
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This function generates the expanded matrix Z' (Z prime) required for the
% Two-Step Ridge Regression (2SRR) procedure. The matrix Z' includes the
% original data matrix X along with additional columns of lagged variables.
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% data:   Matrix of data series with lags, of size n x (m x p), where:
%         - n is the number of observations
%         - m is the number of variables
%         - p is the number of lags
%         Example: For n = 220 observations, m = 5 variables, and p = 5 lags,
%         the matrix size would be 220 x 30.
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% Zprime: Matrix of size n x (m x p + 1), where:
%         - The first column is a constant term (ones).
%         - The subsequent columns include the lagged variables.
%         - The final columns in Zprime include the original data matrix X.
%
%         Specifically, Zprime combines the matrix of lagged variables with
%         the original data matrix to form the expanded matrix needed for 2SRR.
%
%--------------------------------------------------------------------------
%
% Author: Christophe Barrette, September 2024
% Based on the function by Goulet Coulombe (2020)
%--------------------------------------------------------------------------

% Get the number of observations (n) and the number of variables with lags (k)
[n, k] = size(data);

% Initialize the expanded matrix with a constant term (ones) and the original data
X = [ones(n, 1), data];
Z = zeros(n, n-1, k+1);

% Fill the Z matrix with the lagged variables
for tt = 2:n
    Z(tt, 1:(tt-1), :) = repmat(X(tt, :), tt-1, 1);
end

% Combine the lagged variables to form Zprime
Zprime = Z(:, :, 1);
for kk = 2:(k+1)
    Zprime = [Zprime, Z(:, :, kk)];
end

% Append the original data matrix X to the expanded Zprime matrix
Zprime = [Zprime, X];






