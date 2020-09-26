clc
clear
close all

%% Script Information
% This script attempts to apply constrained optimization, to find the
% combination of parameters that is 'most robust'. This is done by
% minimizing the weighted 'spread' in both variables, where the radii a and
% b are understood to be the square roots of the eigenvalue matrix

% Make use of the first and second order taylor approximations

%% Apply constrained optimization
load('DerivedFactors.mat', 'DerivedFactors');
load('SensitivityAnalysis.mat', 'SensitivityAnalysis');

% Define symbolic variables
vars = sym('vars', [1 4]);              %mu, phi_1, phi_2, r_frac

% Assumed variances and covariance matrix
Var_mu              = SensitivityAnalysis.Variance.mu;
Var_phi_1           = SensitivityAnalysis.Variance.phi_1;
Var_phi_2           = SensitivityAnalysis.Variance.phi_2;
Var_r_fraction      = SensitivityAnalysis.Variance.r_fraction;
input_covariance_matrix     = SensitivityAnalysis.Variance.covariancematrix;

% Entered weights Var(c_1 * xi + c_2 * mu * q1) to be minimized)
c_1     = SensitivityAnalysis.Weights.c1;
c_2     = SensitivityAnalysis.Weights.c2;

% Enforced conditions
global normAmpFactor AmpFactorfunc AmpFactorfunc_2Order
phi_mindiff     = SensitivityAnalysis.Constraint.phi_mindiff;
r_fraction_max  = SensitivityAnalysis.Constraint.r_fraction_max;
normAmpFactor   = SensitivityAnalysis.Constraint.normXi;

% Starting point for constrained optimization. May not be zero
StartPoint      = SensitivityAnalysis.Optimization.startstate;     %mu, phi1, phi2, r_frac

%% Determine the function for the first and second order variances
% First, load the derived factors. mat
load('DerivedFactors.mat', 'DerivedFactors')

% Next determine the first order estimate for the variance as a function of
% the variances
[~, index_mu, ~]        = find(DerivedFactors.symbolicvariables == 'mu');
[~, index_phi_1, ~]     = find(DerivedFactors.symbolicvariables == 'phi_1');
[~, index_phi_2, ~]     = find(DerivedFactors.symbolicvariables == 'phi_2');
[~, index_r_frac, ~]    = find(DerivedFactors.symbolicvariables == 'r_frac');

mu      = DerivedFactors.symbolicvariables(1, index_mu);
phi_1   = DerivedFactors.symbolicvariables(1, index_phi_1);
phi_2   = DerivedFactors.symbolicvariables(1, index_phi_2);
r_frac  = DerivedFactors.symbolicvariables(1, index_r_frac);

fmin_amp    = subs(DerivedFactors.AmpFactor_nondim, [mu, phi_1, phi_2, r_frac], vars);
fmin_muQ    = subs(mu * DerivedFactors.q_1_nondim, [mu, phi_1, phi_2, r_frac], vars); 

AmpFactorfunc       = matlabFunction(fmin_amp, 'Vars', {vars});
q_1func             = matlabFunction(fmin_muQ, 'Vars', {vars});

% Assuming the 4 i.i.d. input variables result in a two dimensional vector
% output, calculate the 2x2 covariance matrix of a first order
% approximation
Jacobian        = jacobian([fmin_amp, fmin_muQ], vars);
Variance_1Order = Jacobian * input_covariance_matrix * ...
    transpose(Jacobian);
Average_1Order  = [fmin_amp; fmin_muQ];

% Continue and calculate the second order variance. First, calculate the
% required hessians
Hessian_1       = hessian(fmin_amp, vars);
Hessian_2       = hessian(fmin_muQ, vars);

Variance_2Order = Variance_1Order + 1/2 .* ...
    [trace(input_covariance_matrix * Hessian_1 * input_covariance_matrix * Hessian_1), ...
    trace(input_covariance_matrix * Hessian_1 * input_covariance_matrix * Hessian_2); ...
    trace(input_covariance_matrix * Hessian_2 * input_covariance_matrix * Hessian_1), ...
    trace(input_covariance_matrix * Hessian_2 * input_covariance_matrix * Hessian_2)];

Average_2Order = Average_1Order + 1/2 .* [trace(Hessian_1 * input_covariance_matrix); ...
    trace(Hessian_2 * input_covariance_matrix)] ;

AmpFactorfunc_2Order    = matlabFunction(Average_2Order(1,1), 'Vars', {vars});
q_1func_2Order          = matlabFunction(Average_2Order(2,1), 'Vars', {vars});

% Now, formulate the function that must be minimized : Var(c_1 * xi + c_2
% * mu * q_1)

% First order with incorporation of the covariance
func_Var_combined_1Order    = matlabFunction(c_1^2 * Variance_1Order(1,1) + ...
    2 * c_1 * c_2 * Variance_1Order(1,2) + c_2^2 * Variance_1Order(2,2),'Vars', {vars});

func_Var_combined_2Order    = matlabFunction(c_1^2 * Variance_2Order(1,1) + ...
    2 * c_1 * c_2 * Variance_2Order(1,2) + c_2^2 * Variance_2Order(2,2),'Vars', {vars});

%% Set constraint matrices for optimization
A = [   [0, 1,-1, 0]
        [0, 0, 1, 0]
        [0,-1, 0, 0]
        [0, 0, 0,-1]
        [0, 0, 0, 1]];
b = [-phi_mindiff; pi; 0; 0; r_fraction_max];

Aeq = [1, 0, 0, 0];
beq = [SensitivityAnalysis.Constraint.normMu];

lb = [];
ub = [];
%% Minimise variances
nonlcon                         = @local_nonlcon_1order;
options                         = optimoptions('fmincon'); 
options.ConstraintTolerance     = 1e-8;
options.StepTolerance           = 1e-10;
% options.Algorithm              = 'sqp';
disp('First Order Taylor Approximation of Variations')
[optim, funValue]       = fmincon(func_Var_combined_1Order, StartPoint, ...
    A, b, Aeq, beq, lb, ub, nonlcon, options)

% Calculate the amplification factor, the value of q_1 and the value of 
% the linearised variances accompanying the found point
AmplificationFactor_1Order      = AmpFactorfunc(optim)
q_1_1Order                      = q_1func(optim)
Variance_xi                     = vpa(subs(Variance_1Order(1,1), vars, optim), 5)
Covariance_xi_muq_1             = vpa(subs(Variance_1Order(2,1), vars, optim), 5)
Variance_muq_1                  = vpa(subs(Variance_1Order(2,2), vars, optim), 5)

% Repeat the procedure with the second order approximation
disp('Second Order Taylor Approximation of Variations')
nonlcon                         = @local_nonlcon_2order;
[optim, funValue]       = fmincon(func_Var_combined_2Order, StartPoint, ...
    A, b, Aeq, beq, lb, ub, nonlcon, options)

AmplificationFactor_2Order      = AmpFactorfunc_2Order(optim)
q_1_2Order                      = q_1func_2Order(optim)
Variance_xi                     = vpa(subs(Variance_1Order(1,1), vars, optim), 5)
Covariance_xi_muq_1             = vpa(subs(Variance_1Order(2,1), vars, optim), 5)
Variance_muq_1                  = vpa(subs(Variance_1Order(2,2), vars, optim), 5)

%% Declaration of nonlinear constraint
function [c, ceq] = local_nonlcon_1order(x)
    global normAmpFactor AmpFactorfunc
    c = [];
    ceq = normAmpFactor - AmpFactorfunc(x);
end

function [c, ceq] = local_nonlcon_2order(x)
    global normAmpFactor AmpFactorfunc_2Order
    c = [];
    ceq = normAmpFactor - AmpFactorfunc_2Order(x);
end