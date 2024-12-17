function [Ymat, Xmat] = fXMAT(Y,vlag)

%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This function prepares the matrices required for performing a Vector
% Autoregression (VAR) with a specified lag length.
%
% It returns the matrices Ymat and Xmat, where Ymat contains the observations
% with the necessary lags removed, and Xmat includes the variables along with
% their lags arranged in columns.
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% Y:     Matrix of data series, of size m x n, where m is the number of variables
%        and n is the number of time periods.
%
% vlag:  Scalar representing the lag length for the VAR model.
%
%--------------------------------------------------------------------------
%
% OUTPUTS:
%
% Ymat:  Matrix of size (T-vlag) x m, where each row corresponds to an observation
%        from the data series Y with lags applied. Specifically, Ymat contains the 
%        observations with the necessary lags removed to align with the lag structure.
%
% Xmat:  Matrix of size (T-vlag) x (m*vlag), where each row consists of the variables
%        and their corresponding lagged values. The columns represent the lagged 
%        variables from 1 up to the specified lag length.
%
%--------------------------------------------------------------------------
%
% Author:  Christophe Barrette, September 2024
%
%--------------------------------------------------------------------------

% Determine the number of observations (T) and the number of variables (M)
[T, M] = size(Y);

% Calculate the number of observations after applying the lag length
Tbig = T - vlag;

% Initialize the matrix to store both the observations and their lags
allY = zeros(Tbig, (vlag + 1) * M);

% Initialize the column index for storing lagged variables
uplcc = 1;

% Populate the matrix with the observations and their lags
for ii = 0:vlag
    % Extract the relevant portion of Y for the current lag
    allY(:, uplcc:uplcc + M - 1) = Y(1 + vlag - ii:T - ii, 1:M);
    % Update the column index for the next set of lagged variables
    uplcc = uplcc + M;
end

% Separate the lagged variables from the original observations
allYlag = allY(:, M + 1:end);

% Assign the matrices Ymat and Xmat
Ymat = allY(:, 1:M);   % Observations matrix with lags applied
Xmat = allYlag;        % Lagged variables matrix
