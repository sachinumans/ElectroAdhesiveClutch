function [Average, Variance] = TaylorApproximation_1Order(DerivedFactors, ...
    input_covariance_matrix, mu_manual, phi_1_manual, phi_2_manual, r_frac_manual)
% This function approximates the non-linear functions of mu * q_1 and the
% amplification factor as linearisations and attempts to determine the
% average and variance of the linearised functions

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

%Define the function we will be using (f(x) = [Amp, mu * q_1])
VectorFunc          = matlabFunction([fmin_amp, fmin_muQ], 'Vars', {vars});

%% First Order Approximation with Covariance Matrix
% Calculate jacobian matrix of the vector function and use it to calculate
% the symbolic expression for the variance
Jacobian    = jacobian([fmin_amp, fmin_muQ], vars);

Variance_1OrderFunc     = matlabFunction(Jacobian * input_covariance_matrix * ...
    transpose(Jacobian), 'Vars', {vars});

% Express the average and standard deviation in terms of the factors above
Average     = transpose(VectorFunc([mu_manual, phi_1_manual, phi_2_manual, r_frac_manual]));
Variance    = Variance_1OrderFunc([mu_manual, phi_1_manual, phi_2_manual, r_frac_manual]);
end

