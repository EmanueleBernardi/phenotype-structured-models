% =========================================================================
% MAIN_PDE
% Numerical solution of the non-local PDE model
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

%% Initialise
U_current = u0';
RHO       = zeros(1, length(t));  RHO(1)  = rho0;
MEAN      = zeros(1, length(t));  MEAN(1) = v_0;
EP        = zeros(1, length(t));  EP(1)   = 1/s0 + v_0^2;
U_mid     = [];
U_end     = [];
mid_index = round(length(t) / 7);

%% Time loop
for k = 1:length(t) - 1

    % Advection-diffusion term
    dfdv = -alpha * central_diff(U_current, h) + ...
            beta/2 * second_diff(U_current, h);

    % Selection term
    R_vec = r_vec - RHO(k);

    % Explicit Euler update
    U_next = U_current + Deltat * (R_vec .* U_current + dfdv);

    % Update moments
    RHO(k+1)  = trapz(x, U_next);
    MEAN(k+1) = trapz(x, x' .* U_next) / RHO(k+1);
    EP(k+1)   = trapz(x, (x.^2)' .* U_next) / RHO(k+1);

    % Snapshots
    if k == mid_index
        U_mid = U_next;
    elseif k == length(t) - 1
        U_end = U_next;
    end

    U_current = U_next;

end

%% Save
save(sprintf('../data/PDE_alpha_%.1f.mat', alpha), ...
    'x', 't', ...
    'RHO', 'MEAN', 'EP', ...
    'U_mid', 'U_end');


%% Auxiliary functions
function dfdv = central_diff(f, dx)
    dfdv      = (circshift(f, [-1, 0]) - circshift(f, [1, 0])) / (2*dx);
    dfdv(1)   =  f(2) / (2*dx);
    dfdv(end) = -f(end-1) / (2*dx);
end

function d2fdv2 = second_diff(f, dx)
    d2fdv2      = (circshift(f, [-1, 0]) - 2*f + circshift(f, [1, 0])) / dx^2;
    d2fdv2(1)   = (-2*f(1) + f(2)) / dx^2;
    d2fdv2(end) = (f(end-1) - 2*f(end)) / dx^2;
end