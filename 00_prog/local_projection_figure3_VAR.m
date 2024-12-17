% Author:  Christophe Barrette, September 2024
% Based on the function in Goulet Coulombe (2024)
% Code generating the figure 3 in Goulet Coulombe (2024)

clear; close all; clc;

wd ='INSERT YOUR PATH HERE/Empirical/';
cd(wd);

% set the path to yours
paths.dat = [wd '10_data/'];      
paths.too = [wd '20_tools/'];     
addpath(paths.dat);
addpath(genpath(paths.too));

rng(1234);

%% Data processing

raw.data1 = readtable( [paths.dat,'/cs18_fig4_ppuzzle.xls'] );

gdp = table2array(raw.data1(:, 1)); 
cpi = table2array(raw.data1(:, 2)); 
pcom = table2array(raw.data1(:, 3)); 
ur = table2array(raw.data1(:, 4)); 
usdcad = table2array(raw.data1(:, 5)); 
expo = table2array(raw.data1(:, 6)); 
imp = table2array(raw.data1(:, 7)); 
mpshock = table2array(raw.data1(:, 8));
gdp_yoy = table2array(raw.data1(:, 9));
ur_lvl = table2array(raw.data1(:, 10));
cpi_yoy = table2array(raw.data1(:, 11));

%% Plot data (Figure 5)
subplot(2,2,1);
plot(ur_lvl);
title('Unemployment Rate');
subplot(2,2,2);
plot(gdp_yoy);
title('Year over Year monthly GDP growth');
subplot(2,2,3);
plot(cpi_yoy);
title('Year over Year monthly Inflation Rate');
subplot(2,2,4);
plot(mpshock);
title('Monetary Policy Shocks');
saveas(gcf,'40_results/figure5.png')


%% Canada 4 variables

% Define the number of lags and block size for the analysis
lags = 24;        % Number of lagged observations to include in the model
block_size = 24;  % Size of the blocks of data to process

% Concatenate and prepare the data matrix Y
Y = [gdp, cpi, pcom, mpshock];  % Combine GDP, CPI, PCOM, and MPS shock data into a single matrix
Y = Y(24:end, :);               % Exclude the first 24 observations to start the estimation at the right period

% Prepare matrices for model fitting
[Ymat, Xmat] = fXMAT(Y, lags);  % Generate matrices for the model using the fXMAT function

% Estimate the model parameters and forecasts
[BETAS_GRR, BETAS_VARF, Bcomp, LAMBDAS, fcast, YHAT_VAR, YHAT_VARF, EWvec] = ...
    tvp_2SRR(Xmat, Ymat, block_size);  % Estimate model parameters, including time-varying parameters (TVP) and forecasts 

% Calculate residuals
RES = Ymat - YHAT_VARF;  % Compute the residuals by subtracting the forecasted values from the actual values

% Compute impulse response functions (IRFs)
[IRY, IRI, IRC, IRMP] = IRF_tvp_2SRR(Bcomp, RES, 48, 4, 1:4);  % Calculate IRFs based on the estimated model parameters

% Normalize the impulse response functions
% Normalize IRFs to ensure a constant shock across time
IRY = (IRY * 100) ./ IRMP(:,1);  % Normalize IRY (GDP response) to the first column of IRMP
IRI = (IRI * 100) ./ IRMP(:,1);  % Normalize IRI (CPI response) to the first column of IRMP
IRC = (IRC) ./ IRMP(:,1);        % Normalize IRC (PCOM response) to the first column of IRMP
IRMP = (IRMP) ./ IRMP(:,1);      % Normalize IRMP (MPS shock response) to the first column of IRMP


% Print the CPI plot
figure;
surfc(IRI(end:-1:1, :));
shading interp
% Add labels and title
ylabel('Years');
xlabel('Horizon');
% Adjust Y-axis limits to ensure all tick labels are visible
ylim([0, 466]);  % Set the Y-axis limit to include all tick positions
% Adjust Y-axis labels to show years from 1975 to 2012
yticks([0, 116, 233, 350, 466]);  % Positions corresponding to the given years
yticklabels({'2012','2005','1990','1983','1975'});  % Labels for each tick
%yticklabels({'1975','1983','1990','2005','2012'});  % Labels for each tick
% Adjust the view angle for better visualization
view(50, 25);
saveas(gcf, '40_results/figure3_inflation_VAR4.png');


% Print the GDP plot
figure;
surfc(IRY(end:-1:1, :));
shading interp
% Add labels and title
ylabel('Years');
xlabel('Horizon');
% Adjust Y-axis limits to ensure all tick labels are visible
ylim([0, 466]);  % Set the Y-axis limit to include all tick positions
% Adjust Y-axis labels to show years from 1975 to 2012
yticks([0, 116, 233, 350, 466]);  % Positions corresponding to the given years
yticklabels({'2012','2005','1990','1983','1975'});  % Labels for each tick
%yticklabels({'1975','1983','1990','2005','2012'});  % Labels for each tick

% Adjust the view angle for better visualization
view(140, 25);
saveas(gcf, '40_results/figure3_GDP_VAR4.png');



%% Canada 8 variables

% Define the number of lags and block size for the analysis
lags = 24;        % Number of lagged observations to include in the model
block_size = 24;  % Size of the blocks of data to process

% Concatenate and prepare the data matrix Y with additional variables
Y = [gdp, cpi, pcom, ur, usdcad, expo, imp, mpshock];  % Combine GDP, CPI, PCOM, unemployment rate (UR), USD/CAD exchange rate, export (EXPO), import (IMP), and MPS shock data into a single matrix
Y = Y(24:end, :);  % Exclude the first 24 observations to start the estimation at the right period

% Prepare matrices for model fitting
[Ymat, Xmat] = fXMAT(Y, lags);  % Generate matrices for the model using the fXMAT function, which processes the data with the specified lag length

% Estimate the model parameters and forecasts
[BETAS_GRR, BETAS_VARF, Bcomp, LAMBDAS, fcast, YHAT_VAR, YHAT_VARF, EWvec] = ...
    tvp_2SRR(Xmat, Ymat, block_size);  % Estimate model parameters, including time-varying parameters (TVP) and forecasts

% Calculate residuals
RES = Ymat - YHAT_VARF;  % Compute the residuals by subtracting the forecasted values from the actual values

% Compute impulse response functions (IRFs)
[IRY, IRI, IRC, IRU, IRUC, IRE, IRIM, IRMP] = IRF_tvp_2SRR(Bcomp, RES, 48, 8, 1:8);  % Calculate IRFs for each variable based on the estimated model parameters

% Normalize the impulse response functions
% Normalize IRFs to ensure a constant shock across time
IRY = (IRY * 100) ./ IRMP(:,1);  % Normalize IRY (GDP response) to the first column of IRMP, scaled by 100
IRI = (IRI * 100) ./ IRMP(:,1);  % Normalize IRI (CPI response) to the first column of IRMP, scaled by 100
IRU = (IRU) ./ IRMP(:,1);        % Normalize IRU (UR response) to the first column of IRMP
IRMP = (IRMP) ./ IRMP(:,1);      % Normalize IRMP (MPS shock response) to the first column of IRMP


% Print the CPI plot
figure;
surfc(IRI(end:-1:1, :));
shading interp
% Add labels and title
ylabel('Years');
xlabel('Horizon');
% Adjust Y-axis limits to ensure all tick labels are visible
ylim([0, 466]);  % Set the Y-axis limit to include all tick positions
% Adjust Y-axis labels to show years from 1975 to 2012
yticks([0, 116, 233, 350, 466]);  % Positions corresponding to the given years
yticklabels({'2012','2005','1990','1983','1975'});  % Labels for each tick
% Adjust the view angle for better visualization
view(50, 25);
saveas(gcf, '40_results/figure3_inflation_VAR8.png');


% Print the GDP plot
figure;
surfc(IRY(end:-1:1, :));
shading interp
% Add labels and title
ylabel('Years');
xlabel('Horizon');
% Adjust Y-axis limits to ensure all tick labels are visible
ylim([0, 466]);  % Set the Y-axis limit to include all tick positions
% Adjust Y-axis labels to show years from 1975 to 2012
yticks([0, 116, 233, 350, 466]);  % Positions corresponding to the given years
yticklabels({'2012','2005','1990','1983','1975'});  % Labels for each tick
% Adjust the view angle for better visualization
view(140, 25);
saveas(gcf, '40_results/figure3_GDP_VAR8.png');


% Print the UR plot
figure;
surfc(IRU(end:-1:1, :));
shading interp
% Add labels and title
ylabel('Years');
xlabel('Horizon');
% Adjust Y-axis limits to ensure all tick labels are visible
ylim([0, 466]);  % Set the Y-axis limit to include all tick positions
% Adjust Y-axis labels to show years from 1975 to 2012
yticks([0, 116, 233, 350, 466]);  % Positions corresponding to the given years
yticklabels({'2012','2005','1990','1983','1975'});  % Labels for each tick
% Adjust the view angle for better visualization
view(340, 20);
saveas(gcf, '40_results/figure3_UR_VAR8.png');