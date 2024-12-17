function [ZX] = cumulzeros(X);

                                
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This code processes a binary matrix X by performing logical negation,
% then identifies and transforms blocks of consecutive values. The transformed 
% matrix Z is generated where the end of each zero-block is replaced with 
% the negative length of the preceding one-block. Finally, the code computes 
% the matrix ZX by element-wise multiplying X with the cumulative sum of Z.
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% X:     A binary matrix where each element is either 0 or
%        (n×m, where n is the number of observations, and m 
%        is the number of variables)
%
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% ZX:  A matrix resulting from the element-wise multiplication 
%      of the original matrix X with the cumulative sum of matrix Z.
%
%--------------------------------------------------------------------------
%
% Author:  Christophe Barrette, September 2024
% Based on the function in Goulet Coulombe (2024)
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Invert X (logical negation)
X = ~X;

% Number of observations
n = size(X, 1);

% Initialize arrays to store values and their lengths
values = [];
lengths = [];

% Identify the distinct values and their lengths in X
i = 1;
while i <= n
    currentValue = X(i);
    j = i + 1;
    while j <= n && X(j) == currentValue
        j = j + 1;
    end
    values = [values, currentValue];
    lengths = [lengths, j - i];
    i = j;
end

% Convert values and lengths to arrays
v = values;
len = lengths;

% Calculate cumulative lengths for block ends
cumlen = cumsum(len);

% Initialize Z matrix as a double precision copy of X
Z = double(X);

% Identify blocks of zeros and replace the end of each zero-block
% with the negative length of the preceding block of ones
iDrops = [0 diff(v)] < 0; % Logical array for zero-block ends
Z(cumlen(iDrops)) = -len(logical([iDrops(2:end), 0]));

% Compute ZX by element-wise multiplying X with the cumulative sum of Z
ZX = X .* cumsum(Z);
