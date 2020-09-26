clear; clc; close all;
%% Load factors and variables
load('DerivedFactors.mat', 'DerivedFactors');
load('SensitivityAnalysis.mat', 'SensitivityAnalysis');

[~, index_mu, ~]        = find(DerivedFactors.symbolicvariables == 'mu');
[~, index_phi_1, ~]     = find(DerivedFactors.symbolicvariables == 'phi_1');
[~, index_phi_2, ~]     = find(DerivedFactors.symbolicvariables == 'phi_2');
[~, index_r_frac, ~]    = find(DerivedFactors.symbolicvariables == 'r_frac');

mu      = DerivedFactors.symbolicvariables(1, index_mu);
phi_1   = DerivedFactors.symbolicvariables(1, index_phi_1);
phi_2   = DerivedFactors.symbolicvariables(1, index_phi_2);
r_frac  = DerivedFactors.symbolicvariables(1, index_r_frac);

debug_plot  = true

%% Settings
nc  = 10000;     %number of simulations per geometry

% Manual addition of configuration of interest. Set boolean to true if
% configuration should be examined
mu_manual       = SensitivityAnalysis.ManualConfig.mu;
phi_1_manual    = SensitivityAnalysis.ManualConfig.phi_1;
phi_2_manual    = SensitivityAnalysis.ManualConfig.phi_2;
r_frac_manual   = SensitivityAnalysis.ManualConfig.r_fraction;

% Assumed variances
Var_mu          = SensitivityAnalysis.Variance.mu;
Var_phi_1       = SensitivityAnalysis.Variance.phi_1;
Var_phi_2       = SensitivityAnalysis.Variance.phi_2;
Var_r_fraction  = SensitivityAnalysis.Variance.r_fraction;

%% Take the viable configuration and add noise
% Define the function to calculate the amplification factor and mu * q_1
AmpFacFunc  = matlabFunction(DerivedFactors.AmpFactor_nondim);
MuQ1Func    = matlabFunction(mu*DerivedFactors.q_1_nondim);

dialogueWindow  = waitbar(0, ...
    'Monte Carlo Simulation - Amplification Factor and mu q_1 (with covariance)');

% For the manually entered configuration
local_config.r_frac     = r_frac_manual;
local_config.phi_1      = phi_1_manual;
local_config.phi_2      = phi_2_manual;
local_config.mu         = mu_manual;

for i = 1:nc
    waitbar((i / nc), dialogueWindow, ...
    'Monte Carlo Simulation - Amplification Factor and mu q_1 (with covariance)')
    % For the amount of simulations per geometry
    % Noise with the specified variance is added to all four factors
    noisyConfig.r_frac  = Nnoise(local_config.r_frac, sqrt(Var_r_fraction));
    noisyConfig.phi_1   = Nnoise(local_config.phi_1, sqrt(Var_phi_1));
    noisyConfig.phi_2   = Nnoise(local_config.phi_2, sqrt(Var_phi_2));
    noisyConfig.mu      = Nnoise(local_config.mu, sqrt(Var_mu));

    %These 'noisy' values are used to calculate the 'noisy'
    %amplification factor and value for mu * q_1
    local_config.results(:,i) = [AmpFacFunc(noisyConfig.mu, noisyConfig.phi_1, ...
        noisyConfig.phi_2, noisyConfig.r_frac); MuQ1Func(noisyConfig.mu, ...
        noisyConfig.phi_1, noisyConfig.phi_2, noisyConfig.r_frac)];
end

% Close waitbar
close(dialogueWindow)

% Calculate covariance matrix manually. First calculate average
average     = mean(local_config.results, 2)

% Calculate the covariance, by definition
configs.covariance_matrix(1,1) = 1 / (nc - 1) * ...
    dot(local_config.results(1,:) - average(1,1), ...
    transpose(local_config.results(1,:) - average(1,1)));
configs.covariance_matrix(1,2) = 1 / (nc - 1) * ...
    dot(local_config.results(1,:) - average(1,1), ...
    transpose(local_config.results(2,:) - average(2,1)));
configs.covariance_matrix(2,1) = 1 / (nc - 1) * ...
    dot(local_config.results(2,:) - average(2,1), ...
    transpose(local_config.results(1,:) - average(1,1)));
configs.covariance_matrix(2,2) = 1 / (nc - 1) * ...
    dot(local_config.results(2,:) - average(2,1), ...
    transpose(local_config.results(2,:) - average(2,1)));

disp(join(['Covariance Matrix after', string(nc), 'simulations']))
configs.covariance_matrix

if debug_plot
    scatter(local_config.results(2,:), local_config.results(1,:), 10, 'black', '.')    
    title(join(['Monte Carlo Distribution (N = ', string(nc), ...
    ')'],''))
    grid('on')
    xlabel('$\mu \, q_2$', 'Interpreter', 'Latex')
    ylabel('$\xi$', 'Interpreter', 'Latex')
end

%% Noise function
function [delta] = Nnoise(mean, stddev)
delta = mean + stddev*randn; 
end