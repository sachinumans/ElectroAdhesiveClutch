clc
clear
close all

%% Script Information and Simulation Parameters
% This script aims to plot
%   - A random distribution of points, following from the Monte Carlo
%   method
%   - The covariance ellipses from:
%       - 1st order Taylor approximation of the function and corresponding
%       variance
%       - 2nd order Taylor approximation of the function and corresponding
%       variance
%       - The Monte-Carlo approximation
load("DerivedFactors.mat", "DerivedFactors");
load("SensitivityAnalysis.mat", "SensitivityAnalysis");

% Simulation Parameters
N               = 200;              % number of empirical points to simulate
n_ellipse       = 500;              % number of points to plot the ellipse with
conf_interval   = SensitivityAnalysis.ConfidenceInterval;      
                                    % The confidence interval for the size of the ellipse

% The eucledian distance from the average, after transformation by the
% inverse covariance matrix is given by a 2 DOF chi-squared distribution.
% Thus, this distribution should be used to determine the confidence
% interval

% Input variances and the corresponding covariance matrix
Var_mu          = SensitivityAnalysis.Variance.mu;
Var_phi_1        = SensitivityAnalysis.Variance.phi_1;
Var_phi_2        = SensitivityAnalysis.Variance.phi_2;
Var_r_fraction  = SensitivityAnalysis.Variance.r_fraction;
input_covariance_matrix     = SensitivityAnalysis.Variance.covariancematrix;

% Manual addition of configuration of interest
mu_manual       = SensitivityAnalysis.ManualConfig.mu;
phi_1_manual    = SensitivityAnalysis.ManualConfig.phi_1;
phi_2_manual    = SensitivityAnalysis.ManualConfig.phi_2;
r_frac_manual   = SensitivityAnalysis.ManualConfig.r_fraction;

% Manual addition of Monte Carlo results if boolean set to true
MC_bool         = true;
MC_average      = [2.6909; 0.6302];
MC_covariance   = [0.0232, 0.0031; 0.0031, 4.265518011552103e-04];

%% Random distribution of points
% Load the variables stored in DerivedFactors
[~, index_mu, ~]        = find(DerivedFactors.symbolicvariables == 'mu');
[~, index_phi_1, ~]     = find(DerivedFactors.symbolicvariables == 'phi_1');
[~, index_phi_2, ~]     = find(DerivedFactors.symbolicvariables == 'phi_2');
[~, index_r_frac, ~]    = find(DerivedFactors.symbolicvariables == 'r_frac');

mu      = DerivedFactors.symbolicvariables(1, index_mu);
phi_1   = DerivedFactors.symbolicvariables(1, index_phi_1);
phi_2   = DerivedFactors.symbolicvariables(1, index_phi_2);
r_frac  = DerivedFactors.symbolicvariables(1, index_r_frac);

% Randomly generate N points, with average located at the manual
% configuration
AmpFacFunc  = matlabFunction(DerivedFactors.AmpFactor_nondim);
MuQ1Func    = matlabFunction(mu*DerivedFactors.q_1_nondim);

% For the manually entered configuration
local_config.r_frac     = r_frac_manual;
local_config.phi_1      = phi_1_manual;
local_config.phi_2      = phi_2_manual;
local_config.mu         = mu_manual;

for i = 1:N
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
% The data is thus stored in local_config.results
 
%% Calculate Averages and Variances via First and Second Order Taylor Approximations
[Average_1Order, Variance_1Order]   = TaylorApproximation_1Order( ...
    DerivedFactors, input_covariance_matrix, mu_manual, phi_1_manual, ...
    phi_2_manual, r_frac_manual);
[Average_2Order, Variance_2Order]   = TaylorApproximation_2Order( ...
    DerivedFactors, input_covariance_matrix, mu_manual, phi_1_manual, ...
    phi_2_manual, r_frac_manual);

% Compute first and second order ellipses
t   = linspace(0, 2*pi, 500);

[V_1, D_1]  = eig(Variance_1Order);
[V_2, D_2]  = eig(Variance_2Order);

theta_1   = atan2(D_1(1,1) - Variance_1Order(1,1), Variance_1Order(1,2));
theta_2   = atan2(D_2(1,1) - Variance_2Order(1,1), Variance_2Order(1,2));

x_ellipse_1order  = conf_interval * sqrt(D_1(1,1)) * cos(theta_1) * cos(t) - ...
    conf_interval * sqrt(D_1(2,2)) * sin(theta_1) * sin(t) + Average_1Order(1);
y_ellipse_1order  = conf_interval * sqrt(D_1(1,1)) * sin(theta_1) * cos(t) + ...
    conf_interval * sqrt(D_1(2,2)) * cos(theta_2) * sin(t) + Average_1Order(2);
x_ellipse_2order  = conf_interval * sqrt(D_2(1,1)) * cos(theta_2) * cos(t) - ...
    conf_interval * sqrt(D_2(2,2)) * sin(theta_2) * sin(t) + Average_2Order(1);
y_ellipse_2order  = conf_interval * sqrt(D_2(1,1)) * sin(theta_2) * cos(t) + ...
    conf_interval * sqrt(D_2(2,2)) * cos(theta_2) * sin(t) + Average_2Order(2);

if MC_bool
    [V_MC, D_MC]  = eig(MC_covariance);
    
    theta_MC  = atan2(D_MC(1,1) - MC_covariance(1,1), MC_covariance(1,2));
    
    x_MC  = conf_interval * sqrt(D_MC(1,1)) * cos(theta_MC) * cos(t) - ...
        conf_interval * sqrt(D_MC(2,2)) * sin(theta_MC) * sin(t) + MC_average(1);
    y_MC  = conf_interval * sqrt(D_MC(1,1)) * sin(theta_MC) * cos(t) + ...
        conf_interval * sqrt(D_MC(2,2)) * cos(theta_MC) * sin(t) + MC_average(2);
end

%% Plot results
figure()
hold('on')
s0      = scatter(local_config.results(2,:), local_config.results(1,:), 20, 'black', '.');
s1      = scatter(Average_1Order(2), Average_1Order(1), 30, 'red', '+');
s2      = scatter(Average_2Order(2), Average_2Order(1), 30, 'green', '+');
p1      = plot(y_ellipse_1order, x_ellipse_1order, 'red', ...
    'DisplayName', '1\textsuperscript{st} Order Taylor Approximation');
p2      = plot(y_ellipse_2order, x_ellipse_2order, 'green', ...
    'DisplayName', '2\textsuperscript{nd} Order Taylor Approximation');
legendplots     = [p1, p2];
if MC_bool
    s3  = scatter(MC_average(2), MC_average(1), 30, 'blue', '+');
    p3  = plot(y_MC, x_MC, 'blue', 'DisplayName', join(["Monte Carlo simulation (N = ", ...
        string(N), ")"], ""));
    legendplots     = [p1, p2, p3];
end
legend(legendplots, 'Location', 'best', 'Interpreter', 'Latex', 'Fontsize', 11);
hold('off')
grid('on')
% title(join(["Monte Carlo Distribution (N = ", string(N), ...
%     ") with predicted averages and variances"],""))
xlabel("$\mu \, q_1$", 'Interpreter', 'Latex')
ylabel("$\xi$", 'Interpreter', 'Latex')

%% Noise Function
function [delta] = Nnoise(mean, stddev)
delta = mean + stddev*randn; 
end