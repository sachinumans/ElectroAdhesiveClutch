function [Average, Variance] = TaylorApproximation_2Order(DerivedFactors, ...
    input_covariance_matrix, mu_manual, phi_1_manual, phi_2_manual, r_frac_manual)
% This function attempts to determine the average and the variance in two
% non-linear scalar functions as a function of the second order taylor
% approximation
[~, index_mu, ~]        = find(DerivedFactors.symbolicvariables == 'mu');
[~, index_phi_1, ~]     = find(DerivedFactors.symbolicvariables == 'phi_1');
[~, index_phi_2, ~]     = find(DerivedFactors.symbolicvariables == 'phi_2');
[~, index_r_frac, ~]    = find(DerivedFactors.symbolicvariables == 'r_frac');

% Define new symbolic variables
vars = sym('vars', [1 4]);              %mu, phi_1, phi_2, r_frac

% Define the functions that will be used to calculate the variance
fmin_amp = DerivedFactors.AmpFactor_nondim; 
fmin_amp = subs(fmin_amp, DerivedFactors.symbolicvariables(1, index_mu), vars(1));
fmin_amp = subs(fmin_amp, DerivedFactors.symbolicvariables(1, index_phi_1), vars(2));
fmin_amp = subs(fmin_amp, DerivedFactors.symbolicvariables(1, index_phi_2), vars(3));
fmin_amp = subs(fmin_amp, DerivedFactors.symbolicvariables(1, index_r_frac), vars(4));
 
fmin_muQ = DerivedFactors.symbolicvariables(1, index_mu) * ...
    DerivedFactors.q_1_nondim;
fmin_muQ = subs(fmin_muQ, DerivedFactors.symbolicvariables(1, index_mu), vars(1));
fmin_muQ = subs(fmin_muQ, DerivedFactors.symbolicvariables(1, index_phi_1), vars(2));
fmin_muQ = subs(fmin_muQ, DerivedFactors.symbolicvariables(1, index_phi_2), vars(3));
fmin_muQ = subs(fmin_muQ, DerivedFactors.symbolicvariables(1, index_r_frac), vars(4));

% Define the function we will be using (f(x) = [Amp, mu * q_1])
VectorFunc          = matlabFunction([fmin_amp, fmin_muQ], 'Vars', {vars});

% First, calculate the first order average and variance
Jacobian    = jacobian([fmin_amp, fmin_muQ], vars);

Variance_1OrderFunc     = matlabFunction(Jacobian * input_covariance_matrix * ...
    transpose(Jacobian), 'Vars', {vars});

% Express the average and standard deviation in terms of the factors above
Average_1Order  = transpose(VectorFunc([mu_manual, phi_1_manual, phi_2_manual, r_frac_manual]));
Variance_1Order = Variance_1OrderFunc([mu_manual, phi_1_manual, phi_2_manual, r_frac_manual]);

% Define the hessians needed to determine the second order average and
% variance
Hessian_1_func  = matlabFunction(hessian(fmin_amp, vars), 'Vars', {vars});
Hessian_2_func  = matlabFunction(hessian(fmin_muQ, vars), 'Vars', {vars});

% Calculate the second order average and variance
Average     = Average_1Order + 1/2 .* [trace(Hessian_1_func([mu_manual, ...
    phi_1_manual, phi_2_manual, r_frac_manual]) * input_covariance_matrix); ...
    trace(Hessian_2_func([mu_manual, phi_1_manual, phi_2_manual, r_frac_manual]) * ...
    input_covariance_matrix)];

% Calculate the 2nd order variance
Variance    = Variance_1Order + 1/2 .* ...
    [trace(input_covariance_matrix * Hessian_1_func([mu_manual, phi_1_manual, ...
    phi_2_manual, r_frac_manual]) * input_covariance_matrix * ...
    Hessian_1_func([mu_manual, phi_1_manual, phi_2_manual, r_frac_manual])), ...
    trace(input_covariance_matrix * Hessian_1_func([mu_manual, phi_1_manual, ...
    phi_2_manual, r_frac_manual]) * input_covariance_matrix * ...
    Hessian_2_func([mu_manual, phi_1_manual, phi_2_manual, r_frac_manual])); ...
    trace(input_covariance_matrix * Hessian_2_func([mu_manual, phi_1_manual, ...
    phi_2_manual, r_frac_manual]) * input_covariance_matrix * ...
    Hessian_1_func([mu_manual, phi_1_manual, phi_2_manual, r_frac_manual])), ...
    trace(input_covariance_matrix * Hessian_2_func([mu_manual, phi_1_manual, ...
    phi_2_manual, r_frac_manual]) * input_covariance_matrix * ...
    Hessian_2_func([mu_manual, phi_1_manual, phi_2_manual, r_frac_manual]))];
end

