% =========================================================================
% PLOTTING
% Plots the results of the numerical simulations
% "Derivation and quasi-invariant asymptotics of
%  phenotype-structured integro-differential models"
% E. Bernardi, T. Lorenzi, A. Tosin
% =========================================================================

clear all; close all; clc;


%% Load parameters
load('../data/parameters.mat');

%% Load PDE solution
load(sprintf('../data/PDE_alpha_%.1f.mat', alpha));

%% Figure settings
figure;
set(gcf, 'Position', [100 100 1400 900])
tiledlayout(3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

timestep = 250;
spacestep = 5;
 
%% Loop over epsilon
for ee = 1:3
    
    % Load IDE & MC solutions
    load(sprintf('../data/IDE_eps_%d_alpha_%.1f.mat', ee, alpha));
    load(sprintf('../data/MC_eps_%d_alpha_%.1f.mat', ee, alpha));
   
    % Row offset
    row = (ee - 1) * 4;

    % --- f(v,T) ---
    nexttile(row + 1)
    plot(x, UG_end, 'b-', 'LineWidth', 2)
    hold on
    plot(x, U_end, 'r:', 'LineWidth', 2.5)
    plot(x(1:spacestep:end), MC_sol_end(1:spacestep:end), 'mo', 'LineWidth', 1)
    plot([vm, vm], [0, max(UG_end) + 0.1], '--k', 'LineWidth', 1)
    xlim([vm-3, vm+3])
    ylim([0, 0.6])
    ylabel(['$\varepsilon^2=', num2str(epsilon^2), '$'], ...
        'Interpreter', 'latex', 'FontSize', 16)
    if ee == 1
        title('$f(v,T)$', 'Interpreter', 'latex', 'FontSize', 16)
        legend('IDE', 'PDE', 'MC', 'Interpreter', 'latex', ...
            'FontSize', 12, 'Location', 'Northeast')
        legend box off
    end
    grid off
    set(gca, 'TickLabelInterpreter', 'latex')
    if ee == 3
        xlabel('$v$', 'Interpreter', 'latex', 'FontSize', 16)
    end

    % --- rho(t) ---
    nexttile(row + 2)
    plot(t, RHOG, 'b-', 'LineWidth', 1.5)
    hold on
    plot(t, RHO, 'r:', 'LineWidth', 1.5)
    plot(t_mc(1:timestep:end), Mass(1:timestep:end), 'mo', 'LineWidth', 1)
    ylim([0.1, 0.6])
    if ee == 1
        title('$\rho(t)$', 'Interpreter', 'latex', 'FontSize', 16)
        legend('IDE', 'PDE', 'MC', 'Interpreter', 'latex', ...
            'FontSize', 12, 'Location', 'Southeast')
        legend box off
    end
    grid off
    set(gca, 'TickLabelInterpreter', 'latex')
    if ee == 3
        xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
    end

    % --- p(t) ---
    nexttile(row + 3)
    plot(t, MEANG, 'b-', 'LineWidth', 1.5)
    hold on
    plot(t, MEAN, 'r:', 'LineWidth', 1.5)
    plot(t_mc(1:timestep:end), Mean(1:timestep:end), 'mo', 'LineWidth', 1)
    ylim([0, max(MEANG) + 0.2])
    if ee == 1
        title('$p(t)$', 'Interpreter', 'latex', 'FontSize', 16)
        legend('IDE', 'PDE', 'MC', 'Interpreter', 'latex', ...
            'FontSize', 12, 'Location', 'Southeast')
        legend box off
    end
    grid off
    set(gca, 'TickLabelInterpreter', 'latex')
    if ee == 3
        xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
    end

    % --- E(t) ---
    nexttile(row + 4)
    plot(t, EG, 'b-', 'LineWidth', 1.5)
    hold on
    plot(t, EP, 'r:', 'LineWidth', 1.5)
    plot(t_mc(1:timestep:end), E(1:timestep:end), 'mo', 'LineWidth', 1)
    if ee == 1
        title('$E(t)$', 'Interpreter', 'latex', 'FontSize', 16)
        legend('IDE', 'PDE', 'MC', 'Interpreter', 'latex', ...
            'FontSize', 12, 'Location', 'Southeast')
        legend box off
    end
    grid off
    set(gca, 'TickLabelInterpreter', 'latex')
    if ee == 3
        xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
    end

end
set(gcf, 'Position', [0 0 1400 900], 'Units', 'pixels');
exportgraphics(gcf, sprintf('../data/figure_alpha_%.1f.pdf', alpha), ...
    'ContentType', 'vector');