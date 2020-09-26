clc
clear
close all

%% Script Information
% This script stores all of variables related to the sensitivity analysis
% in one struct

%% Step 1 - Enter \Sigma_x, \vec{c} and additional constraints on parameter space
% Input Variances and covariance matrix
SensitivityAnalysis.Variance.mu              = 0.02^2;
SensitivityAnalysis.Variance.phi_1           = (deg2rad(2.5)/3)^2;
SensitivityAnalysis.Variance.phi_2           = (deg2rad(2.5)/3)^2;
SensitivityAnalysis.Variance.r_fraction      = (0.01/3)^2;

SensitivityAnalysis.Variance.covariancematrix    = ...
    diag([SensitivityAnalysis.Variance.mu, SensitivityAnalysis.Variance.phi_1, ...
    SensitivityAnalysis.Variance.phi_2, SensitivityAnalysis.Variance.r_fraction], 0);

% Weights to examine Var(c_1 * xi + c_2 * mu * q1)
SensitivityAnalysis.Weights.c1  = 1;
SensitivityAnalysis.Weights.c2  = 5;

% Constraints to take into account
SensitivityAnalysis.Constraint.phi_mindiff      = (2 * pi) / 6;
SensitivityAnalysis.Constraint.r_fraction_max   = 0.9;
SensitivityAnalysis.Constraint.normXi           = 2.7;
SensitivityAnalysis.Constraint.normMu           = 0.63;
% note - normXi = normAmpFactor

% Starting configuration for Optimization
SensitivityAnalysis.Optimization.startstate     = ...
    [0.6, pi/4, 3*pi/4, 0.9];            %mu, phi1, phi2, r_frac
%% Step 2 - After running Constrained_Optimization_and_Covariances.m enter \vec{\mu}_x

% Manual configuration to examine
SensitivityAnalysis.ManualConfig.mu          = 0.6300;
SensitivityAnalysis.ManualConfig.phi_1       = 0.4588;
SensitivityAnalysis.ManualConfig.phi_2       = 2.3427;
SensitivityAnalysis.ManualConfig.r_fraction  = 0.900;

% Confidence interval for covariance ellipse
SensitivityAnalysis.ConfidenceInterval      = sqrt(9.21); 
% The eucledian distance from the average, after transformation by the
% inverse covariance matrix is given by a 2 DOF chi-squared distribution.
% Thus, this distribution should be used to determine the confidence
% interval
% sqrt(5.991) = 95% confidence ellipse
% sqrt(9.21)  = 99% confidence ellipse

% Save the configured struct
save('SensitivityAnalysis.mat', 'SensitivityAnalysis')
disp('Sensitivity Analysis Structure Successfully Updated')