% =========================================================================
% MAIN_IDE
% Numerical solution of the IDE model
% "Derivation and quasi-invariant asymptotics of
%  phenotype-structured integro-differential models"
% E. Bernardi, T. Lorenzi, A. Tosin
% =========================================================================

clear all; close all; clc;

%% Load parameters
load('../data/parameters.mat');

%% Spatial grid
h = (b - a) / (n + 1);
x = a:h:b;

%% Time grid
Deltat = 1e-3;
t      = 0:Deltat:Tf;

%% Initial condition
Int = trapz(x, exp(-s0/2 * (x - v_0).^2));
u0  = (rho0 / Int) * exp(-s0/2 * (x - v_0).^2);

%% Mollifier psi and net proliferation rate r(v)
f_zeta = @(x) 0.5 * (1 + tanh(2 * x ./ (1 - x.^2))) .* (abs(x) < 1);
psi    = @(v) 1 - f_zeta((2*(v - R)/delta) - 1);

Psi = zeros(size(x));
for i = 1:length(x)
    if x(i) > -R - delta && x(i) < -R
        Psi(i) = psi(-x(i));
    elseif x(i) >= -R && x(i) <= R
        Psi(i) = 1;
    elseif x(i) > R && x(i) < R + delta
        Psi(i) = psi(x(i));
    else
        Psi(i) = 0;
    end
end

r_vec = (1 - (x - vm).^2)' .* Psi';

%% Epsilon values
epsilons = [sqrt(1e-0), sqrt(1e-1), sqrt(1e-2)];

%% Loop over epsilon
for ee = 1:3

    epsilon = epsilons(ee);
    mu      = 1 / epsilon^2;

    % Gaussian mutation kernel M_epsilon(v|w)
    K = @(v, w) (1 / sqrt(2*pi*epsilon^2*beta)) * ...
                exp(-((v - w - alpha*epsilon^2).^2) / (2*epsilon^2*beta));

    % Precompute kernel matrix and normalise
    K_mat = zeros(length(x), length(x));
    for i = 1:length(x)
        K_mat(:, i) = K(x(i), x);
    end
    for j = 1:length(x)
        K_mat(j, :) = K_mat(j, :) / trapz(x, K_mat(j, :));
    end

    % Initialise
    U_current  = u0';
    RHOG       = zeros(1, length(t));  RHOG(1)  = rho0;
    MEANG      = zeros(1, length(t));  MEANG(1) = v_0;
    EG         = zeros(1, length(t));  EG(1)    = 1/s0 + v_0^2;
    UG_mid     = [];
    UG_end     = [];
    mid_index  = round(length(t) / 7);

    % Time loop
    for k = 1:length(t) - 1

        % Nonlocal integral term
        I = zeros(length(x), 1);
        for i = 1:length(x)
            I(i) = trapz(x, K_mat(:, i) .* U_current);
        end

        % Selection term
        R_vec = r_vec - RHOG(k);

        % Explicit Euler update
        U_next = U_current + Deltat * (R_vec .* U_current + mu * (I - U_current));

        % Update moments
        RHOG(k+1)  = trapz(x, U_next);
        MEANG(k+1) = trapz(x, x' .* U_next) / RHOG(k+1);
        EG(k+1)    = trapz(x, (x.^2)' .* U_next) / RHOG(k+1);

        % Snapshots
        if k == mid_index
            UG_mid = U_next;
        elseif k == length(t) - 1
            UG_end = U_next;
        end

        U_current = U_next;

    end 

    % Save
    save(sprintf('../data/IDE_eps_%d_alpha_%.1f.mat', ee, alpha), ...
    'epsilon', 'x', 't', ...
    'RHOG', 'MEANG', 'EG', ...
    'UG_mid', 'UG_end');

end