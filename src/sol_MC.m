% =========================================================================
% MAIN_MC
% Monte Carlo simulation of the stochastic agent-based model
% "Derivation and quasi-invariant asymptotics of
%  phenotype-structured integro-differential models"
% E. Bernardi, T. Lorenzi, A. Tosin
% =========================================================================

clear all; close all; clc;

%% Load parameters
load('../data/parameters.mat');

%% Spatial grid
h  = (b - a) / (n + 1);
x  = a:h:b;
 
%% Initial condition
Int = sum(exp(-s0/2 * (x - v_0).^2) * h);
u0  = (rho0 / Int) * exp(-s0/2 * (x - v_0).^2);

%% Time grid
Deltat_mc = 1e-3;
t_mc      = 0:Deltat_mc:Tf;
iter      = length(t_mc);

%% Initial population
mass1 = floor(rho0 * N);
mass0 = N - mass1;
C1    = normrnd(v_0, 1/sqrt(s0), [1, mass1]);
C0    = zeros(1, mass0);
V     = [C1, C0; ones(1, mass1), zeros(1, mass0)]';

%% Epsilon values
epsilons = [sqrt(1e-0), sqrt(1e-1), sqrt(1e-2)];

%% Loop over epsilon
for ee = 1:3

    epsilon = epsilons(ee);
    mu      = 1 / epsilon^2;
    Vstar   = V;

    % Transition probabilities
    Th = (Deltat_mc / (1 + Deltat_mc))         * ones(N, 1);
    MU = (mu * Deltat_mc / (1 + mu * Deltat_mc)) * ones(N, 1);

    % Initialise moments
    Mass = zeros(1, iter);  Mass(1) = rho0;
    Mean = zeros(1, iter);  Mean(1) = v_0;
    E    = zeros(1, iter);  E(1)    = 1/s0 + v_0^2;

    % Initialise distributions
    MC_sol_mid = zeros(length(x), 1);
    MC_sol_end = zeros(length(x), 1);

    % Index for mid-time snapshot
    mid_idx = round(iter / 7);

    % Time loop
    for i = 1:iter - 1

        Vstar = binary(Vstar, Th, MU, Deltat_mc, R, vm, ...
                          alpha, beta, epsilon, delta, 0);

        Active = find(Vstar(:, 2) == 1);
        mc1    = histc(Vstar(Active, 1), x);
        Normal = 1 / sum(mc1 * h);

        Mass(i+1) = length(Active) / N;
        Mean(i+1) = Normal * sum(x' .* mc1 * h);
        E(i+1)    = Normal * sum((x.^2)' .* mc1 * h);

        if i + 1 == mid_idx
            MC_sol_mid = mc1 / sum(mc1 * h) * Mass(i+1);
        elseif i + 1 == iter
            MC_sol_end = mc1 / sum(mc1 * h) * Mass(i+1);
        end

    end

    % Save
    save(sprintf('../data/MC_eps_%d_alpha_%.1f.mat', ee, alpha), ...
    'epsilon', 'x', 't_mc', 'Deltat_mc', ...
    'Mass', 'Mean', 'E', ...
    'MC_sol_mid', 'MC_sol_end');

end