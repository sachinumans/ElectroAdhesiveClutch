clc
clear all
close all

%% Background Information
% This script is to be used to supplement the manual calculations performed
% in Alvaro's notebook

%% Moment Balance - Calculating FN1
syms r r_R theta_R theta_bar phi_1 phi_2 phi_bar ...
     normF_0 normF_N0 normF_Br0 normF_N1 mu real positive
            
assume(0 < phi_1 < phi_2 < pi)

%Define angle definitions
% 1. Relative angle to middle of shoe
theta_bar_def   = (phi_1 + phi_2) / 2;
% 2. Relative angle to CoP of sinusoidal pressure distribution
phi_bar_def     = (phi_1 * cos(phi_1) - sin(phi_1) - phi_2 * cos(phi_2) + sin(phi_2)) / ...
                (cos(phi_1) - cos(phi_2));

% Define the two position vectors
% 1. Arm from shoe hinge to electrostatic force
r_Sm_R      = [-r*sin(theta_R)*sin(theta_bar) + cos(theta_R)*(r*cos(theta_bar) - r_R); ...
                r*cos(theta_R)*sin(theta_bar) + sin(theta_R)*(r*cos(theta_bar) - r_R);
                0];   
% 2. Arm from shoe hinge to center of sinusoidal center of pressure
r_Scop_R    = [-r*sin(theta_R)*sin(phi_bar) + cos(theta_R)*(r*cos(phi_bar) - r_R); ...
                r*cos(theta_R)*sin(phi_bar) + sin(theta_R)*(r*cos(phi_bar) - r_R);
                0];        

% Define the four force vectors
% 1. Electrostatic force
F_0         = normF_0 .* [cos(theta_R + theta_bar); ...
                         sin(theta_R + theta_bar); ...
                         0];
% 2. Normal Force, iteration 0
normF_N0     = normF_0;
F_N0        = -1 * normF_N0 .* [cos(theta_R + phi_bar); ...
                               sin(theta_R + phi_bar); ...
                               0];
% 3. Braking force, iteration 0 (thus resulting from F_N0)
F_Br0       = normF_N0 * mu .* [sin(theta_R + phi_bar); ...
                               -cos(theta_R + phi_bar); ...
                               0];
% 4. Normal force, iteration 1 (unknown in equation)
F_N1        = -1 * normF_N1 .* [cos(theta_R + phi_bar); ...
                               sin(theta_R + phi_bar); ...
                               0];

% Calculate moment balance
Moment_Balance_Vector_1     = cross(r_Sm_R, F_0) + cross(r_Scop_R, (F_N0 + F_N1 + F_Br0));
Moment_Balance_Z_1          = Moment_Balance_Vector_1(3,1);

%Solve for normF_N1
normF_N1result      = solve(Moment_Balance_Z_1 == 0, normF_N1, 'ReturnConditions', true);
disp(join(['normF_N1 = ', '', string(subs(simplify(normF_N1result.normF_N1, 'Steps', 100), ...
    [theta_bar, phi_bar], [theta_bar_def, phi_bar_def]))], newline))
%% Calculate q_0
% q_0 is defined as: normF_N1 / normF_Br0
% Define new variable, r_frac = r / r_R
syms r_frac real positive
q_0         = subs(simplify(normF_N1result.normF_N1 / (normF_N0 * mu), 'Steps', 100), ...
    [theta_bar, phi_bar], [theta_bar_def, phi_bar_def])
q_0_nondim  = subs(q_0, [r, r_R], [1, r_frac])

%% Moment balance - Calculating FN2
% All angles and position vectors remain the same. Only new forces must be
% defined
syms normF_N2 real positive
% 1. Braking force, iteration 1
F_Br1       = normF_N1 * mu .* [sin(theta_R + phi_bar); ...
                               -cos(theta_R + phi_bar); ...
                               0];
% 2. Normal force, iteration 2
F_N2        = normF_N2 * -1 .* [cos(theta_R + phi_bar); ...
                               sin(theta_R + phi_bar); ...
                               0];
                           
%Calculat moment balance
Moment_Balance_Vector_2     = cross(r_Scop_R, (F_Br1 + F_N2));
Moment_Balance_Z_2          = Moment_Balance_Vector_2(3,1);

%Solution for F_N2
normF_N2result      = solve(Moment_Balance_Z_2 == 0, normF_N2, 'ReturnConditions', true);
disp(join(['normF_N2 = ', '', string(subs(simplify(normF_N2result.normF_N2, 'Steps', 100), ...
    [theta_bar, phi_bar], [theta_bar_def, phi_bar_def]))], newline))

%% Calculate q_1
% q_1 is defined as: normF_N2 / normF_Br1
q_1         = subs(simplify(normF_N2result.normF_N2 / (normF_N1 * mu), 'Steps', 100), ...
    [theta_bar, phi_bar], [theta_bar_def, phi_bar_def])
q_1_nondim  = subs(q_1, [r, r_R], [1, r_frac])

%% Calculate total ampification factor
AmpFactor           = simplify((1 - mu * (q_1 - q_0)) / (1 - mu * q_1), 'Steps', 100)
AmpFactor_nondim    = simplify((1 - mu * (q_1_nondim - q_0_nondim)) / (1 - mu * q_1_nondim), 'Steps', 100)

%% Store derived factors in structure and save
% Define all of the symbolic variables used, so these can easily be used in
% other scripts
DerivedFactors.symbolicvariables    = [r, r_R, theta_R, theta_bar, phi_1, ...
    phi_2, phi_bar, normF_0, normF_N0, normF_Br0, normF_N1, normF_N2,  mu, r_frac];

% Save the found solutions for the size of the normal forces
DerivedFactors.normF_N1             = subs(simplify(normF_N1result.normF_N1, ...
    'Steps', 100), [theta_bar, phi_bar], [theta_bar_def, phi_bar_def]);
DerivedFactors.normF_N2             = subs(simplify(normF_N2result.normF_N2, ...
    'Steps', 100), [theta_bar, phi_bar], [theta_bar_def, phi_bar_def]);

% Save the found values for q_0, q_1 and the amplification factor
DerivedFactors.q_0                  = q_0;
DerivedFactors.q_0_nondim           = q_0_nondim;
DerivedFactors.q_1                  = q_1;
DerivedFactors.q_1_nondim           = q_1_nondim;
DerivedFactors.AmpFactor            = AmpFactor;
DerivedFactors.AmpFactor_nondim     = AmpFactor_nondim;

save('DerivedFactors.mat', 'DerivedFactors')