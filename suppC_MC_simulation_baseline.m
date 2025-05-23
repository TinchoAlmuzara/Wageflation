%%% MONTE CARLO SIMULATION ON WAGE INFLATION MODEL %%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the simulation study in the paper
% "A Measure of Wage Trend Inflation" (Almuzara, Audoly, Melcangi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear memory
clear
close all
clc
rng(123)

% Set directories
addpath('functions')
fname      = 'baseline';
res_path   = ['testing' filesep fname filesep];
if ~exist(res_path, 'dir'), mkdir(res_path); end

% Decide tasks
estimation = true;


%% MONTE CARLO SIMULATION

% Set numerical parameters and dimensions
n_MC = 200;
T    = 300;
n    = 7;
q    = 1;

% Set true parameters (no time-variation in parameters+no outliers)
param            = struct();
param.lam_tau    = zeros(n, 1);
param.gam_dtau_c = 0;
param.gam_dtau_i = zeros(n, 1);
param.lam_eps    = zeros(n, 1);
param.gam_eps_c  = 0;
param.gam_eps_i  = zeros(n, 1);
param.ps_c       = 1;
param.ps_i       = ones(n, 1);
param.theta      = repmat(0.1./(1:q), n, 1);

% Set true initial time-varying parameters (calibrated from industry data)
param.alpha_tau0    = 0.30*ones(1, n);
param.sigma_dtau_i0 = 0.20*ones(1, n);
param.alpha_eps0    = 0.02*ones(1, n);
param.sigma_eps_i0  = 0.85*ones(1, n);

% Set priors
prior         = struct();
prior.prec_MA = 1/10;
prior.nu_lam  = 60;
prior.s2_lam  = 0.0001;
prior.nu_gam  = 60;
prior.s2_gam  = 0.001;
prior.a_ps    = (1-1/48)*(120);
prior.b_ps    = (1/48)*(120);

% Set number of MA lags and time-aggregated sectors
n_lags             = repmat(3, [n, 1]);
is_timeag          = true(n, 1);
settings           = struct();
settings.n_draw    = 5000;
settings.n_burn    = 5000;
settings.n_thin    = 2;
settings.n_lags    = n_lags;
settings.is_timeag = is_timeag;

% Set summary settings
quants  = [0.16, 0.50, 0.84];
n_quant = length(quants);

if (estimation == true)
% Perform estimation in parallel
parallel_pool = parpool(n_MC);
parfor i_MC = 1:n_MC

    % Simulate data and normalize common components
    [y, latents]         = simulate_data(param, T, n, is_timeag);
    latents.alpha_tau    = latents.alpha_tau .* repmat(latents.sigma_dtau_c, 1, n);
    latents.tau_c        = latents.tau_c ./ latents.sigma_dtau_c;
    latents.sigma_dtau_c = latents.sigma_dtau_c ./ latents.sigma_dtau_c;
    latents.alpha_eps    = latents.alpha_eps .* repmat(latents.sigma_eps_c, 1, n);
    latents.eps_c        = latents.eps_c ./ latents.sigma_eps_c;
    latents.sigma_eps_c  = latents.sigma_eps_c ./ latents.sigma_eps_c;

    % Save average trend components (assume same weights on all sectors)
    latents.trend_c = mean(latents.alpha_tau, 2) .* latents.tau_c;
    latents.trend_i = mean(latents.tau_i, 2); 
    latents.trend   = latents.trend_c + latents.trend_i;
    latents.trend_c = latents.trend_c - mean(latents.trend_c, 1); % Location is identified up to a constant shift
    latents.trend_i = latents.trend_i - mean(latents.trend_i, 1); % Location is identified up to a constant shift
    
    % Estimate model and normalize common components
    output              = estimate(y, prior, settings);
    output.alpha_tau    = output.alpha_tau .* permute(repmat(output.sigma_dtau_c, 1, 1, n), [1 3 2]);
    output.tau_c        = output.tau_c ./ output.sigma_dtau_c;
    output.sigma_dtau_c = output.sigma_dtau_c ./ output.sigma_dtau_c;
    output.alpha_eps    = output.alpha_eps .* permute(repmat(output.sigma_eps_c, 1, 1, n), [1 3 2]);
    output.eps_c        = output.eps_c ./ output.sigma_eps_c;
    output.sigma_eps_c  = output.sigma_eps_c ./ output.sigma_eps_c;

    % Save average trend components (assume same weights on all sectors)
    output.trend_c = squeeze(mean(output.alpha_tau, 2)) .* output.tau_c;
    output.trend_i = squeeze(mean(output.tau_i, 2));
    output.trend   = output.trend_c + output.trend_i;    
    output.trend_c = output.trend_c - mean(output.trend_c, 1); % Location is identified up to a constant shift
    output.trend_i = output.trend_i - mean(output.trend_i, 1); % Location is identified up to a constant shift

    % Save results
    results         = struct();
    results.trend_c = quantile(output.trend_c, quants, 2);
    results.trend_i = quantile(output.trend_i, quants, 2);
    results.trend   = quantile(output.trend, quants, 2);
    results.latents = latents;
    results.y       = y;
    parsave([res_path sprintf('output_%d.mat', i_MC)], results)
    
end
delete(parallel_pool)
end    


%% PROCESSING

% Collect estimation results into struct
latent_names   = {'trend'; 'trend_c'; 'trend_i'};
trend_true     = struct();
trend_estimate = struct();
for obj_name = latent_names'        
    trend_true.(obj_name{:})     = NaN(T, n_MC);
    trend_estimate.(obj_name{:}) = NaN(T, n_quant, n_MC);
    for i_MC = 1:n_MC 
        res_tmp = load([res_path sprintf('output_%d.mat', i_MC)]);
        trend_true.(obj_name{:})(:, i_MC)        = res_tmp.latents.(obj_name{:}); 
        trend_estimate.(obj_name{:})(:, :, i_MC) = res_tmp.(obj_name{:}); 
    end    
end

% Create array with bias
latent_array = (1:T)';
for obj_name = latent_names'
    bias_tmp     = mean(squeeze(trend_estimate.(obj_name{:})(:, 2, :)) - trend_true.(obj_name{:}), 2);
    MSE_tmp      = mean((squeeze(trend_estimate.(obj_name{:})(:, 2, :)) - trend_true.(obj_name{:})).^2, 2);
    coverage_tmp = mean((squeeze(trend_estimate.(obj_name{:})(:, 1, :)) <= trend_true.(obj_name{:})) ...
                   & (trend_true.(obj_name{:}) <= squeeze(trend_estimate.(obj_name{:})(:, 3, :))), 2);
    latent_array = [latent_array, bias_tmp, MSE_tmp, coverage_tmp]; %#ok
end

% Tabulate results for bias
col_names  = {'Period', 'TotalBias', 'TotalMSE', 'TotalCoverage', ...
                        'CommonBias', 'CommonMSE', 'CommonCoverage', ...
                        'IdiosyncraticBias', 'IdiosyncraticMSE', 'IdiosyncraticCoverage'};
latent_tab = array2table(latent_array, 'variablenames', col_names);
disp(latent_tab)
writetable(latent_tab, [res_path 'TabC1_MC_tab.csv'], 'writevariablenames', true, 'writerownames', true)

% Plot estimation error of posterior median
fig_name = {'total trend', 'common trend', 'idiosyncratic trend'};
black    = [  0,   0,   0];
blue     = [ 51,  51, 153]/255;
red      = [153,  51,  51]/255;
fig_fmt  = 'epsc';
figsize  = {'units', 'inches', 'position', [0 0 14 14]};
font     = {'fontname', 'times', 'fontsize', 30};

% Plot estimation error of posterior median
i_fig     = 0;
for obj_name = latent_names'
    bias  = quantile(squeeze(trend_estimate.(obj_name{:})(:, 2, :)) - trend_true.(obj_name{:}), quants, 2);
    i_fig = i_fig + 1;
    fig0  = figure();
    ax0   = axes();
    fill(ax0, [(1:T), fliplr(1:T)], [bias(:, 1)', fliplr(bias(:, 3)')], ...
        blue, 'facealpha', 0.2, 'linestyle', ':', 'linewidth', 0.1, 'edgecolor', blue);
    hold('on')
    plot0    = plot(ax0, (1:T)', bias(:, 2));
    ylabel(ax0, ['Estimation error (' fig_name{i_fig} ')'])
    hold('off')

    % Tune plot
    set(plot0, 'color', blue)
    set(plot0, 'linewidth', 2)
    set(plot0, 'linestyle', '--')

    % Add bias
    hold(ax0, 'on')
    plot1 = plot(ax0, (1:T)', mean(squeeze(trend_estimate.(obj_name{:})(:, 2, :)) - trend_true.(obj_name{:}), 2));
    set(plot1, 'color', black)
    set(plot1, 'linewidth', 3)
    set(plot1, 'linestyle', ':')
    hold(ax0, 'off')

    % Tune ax handle
    set(ax0, font{:})
    set(ax0, 'box', 'on')
    grid(ax0, 'off')
    line(ax0, [1, T], [0, 0], 'color', black, 'linewidth', 2, 'linestyle', '-')

    % Set legend
    legend(ax0, {'P16-P84', 'Median estimation error', 'Mean estimation error'}, 'location', 'best')
    
    % Tune and save figure
    set(fig0, figsize{:});
    print(fig0, [res_path 'FigC1_' num2str(i_fig)], ['-d' fig_fmt])
end


%% LOCAL FUNCTIONS

% Parallel save
function parsave(fname, out)
    save(fname, '-struct', 'out')
end


% Simulate data
function [y, latents] = simulate_data(param, T, n, is_timeag)

% Recover parameters
lam_tau    = param.lam_tau';
gam_dtau_c = param.gam_dtau_c;
gam_dtau_i = param.gam_dtau_i';
lam_eps    = param.lam_eps';
gam_eps_c  = param.gam_eps_c;
ps_c       = param.ps_c;
gam_eps_i  = param.gam_eps_i';
ps_i       = param.ps_i';
theta      = param.theta';

% Recover initial values
tau_c0        = 0;
alpha_tau0    = param.alpha_tau0;
sigma_dtau_c0 = 1;
tau_i0        = zeros(1, n);
sigma_dtau_i0 = param.sigma_dtau_i0;
alpha_eps0    = param.alpha_eps0;
sigma_eps_c0  = 1;
sigma_eps_i0  = param.sigma_eps_i0;

% Set support for s_eps and pi prior
n_s_vals   = 10;
s_eps_vals = [1; linspace(2, 10, n_s_vals-1)'];

% Preallocate data and latent variables
y             = NaN(T, n);
tau_c         = NaN(T, 1);
alpha_tau     = NaN(T, n);
sigma_dtau_c  = NaN(T, 1);
tau_i         = NaN(T, n);
sigma_dtau_i  = NaN(T, n);
eps_c         = NaN(T, 1);
alpha_eps     = NaN(T, n);
sigma_eps_c   = NaN(T, 1);
s_eps_c       = NaN(T, 1);
sigma_eps_i   = NaN(T, n);
s_eps_i       = NaN(T, n);

% Initialize latent variables
q          = size(theta, 1);
eps_i0_vec = zeros(q, n);
tau0_vec   = zeros(12, n);
H_tau      = [ones(1, n); zeros(11, n)]; H_tau(:, is_timeag) = (1/12);
probs_c    = [ps_c, (1-ps_c)/(n_s_vals-1)*ones(1, n_s_vals-1)];
probs_i    = [ps_i', (1-ps_i')./(n_s_vals-1)*ones(1, n_s_vals-1)];

% Simulate data and latent variables
for t = (-18):T
    alpha_tau0    = alpha_tau0 + lam_tau .* randn(1, n);
    sigma_dtau_c0 = sigma_dtau_c0 * exp(gam_dtau_c * randn(1)/2);
    tau_c0        = tau_c0 + sigma_dtau_c0 * randn(1);
    sigma_dtau_i0 = sigma_dtau_i0 .* exp(gam_dtau_i .* randn(1, n)/2);
    tau_i0        = tau_i0 + sigma_dtau_i0 .* randn(1, n);
    alpha_eps0    = alpha_eps0 + lam_eps .* randn(1, n);
    sigma_eps_c0  = sigma_eps_c0 * exp(gam_eps_c * randn(1)/2);
    sigma_eps_i0  = sigma_eps_i0 .* exp(gam_eps_i .* randn(1, n)/2);
    s_eps_c0      = mnrnd(1, probs_c)*s_eps_vals;    
    s_eps_i0      = (mnrnd(1, probs_i)*s_eps_vals)';
    eps_c0        = sigma_eps_c0 * s_eps_c0 * randn(1);
    eps_i0        = sigma_eps_i0 .* s_eps_i0 .* randn(1, n);
    if (q > 0)
        ups_i0     = eps_i0 + sum(theta .* eps_i0_vec, 1);
        eps_i0_vec = [eps_i0; eps_i0_vec(1:(end-1), :)];
    else
        ups_i0 = eps_i0;
    end
    tau0_vec = [alpha_tau0 * tau_c0 + tau_i0; tau0_vec(1:(end-1), :)];

    % Store data and latent variables
    if (t > 0)
        y(t, :)            = sum(H_tau .* tau0_vec, 1) + alpha_eps0 * eps_c0 + ups_i0;
        tau_c(t)           = tau_c0;
        alpha_tau(t, :)    = alpha_tau0;
        sigma_dtau_c(t)    = sigma_dtau_c0;
        tau_i(t, :)        = tau_i0;
        sigma_dtau_i(t, :) = sigma_dtau_i0;
        eps_c(t, :)        = eps_c0;
        alpha_eps(t, :)    = alpha_eps0;
        sigma_eps_c(t)     = sigma_eps_c0;
        s_eps_c(t)         = s_eps_c0;
        sigma_eps_i(t, :)  = sigma_eps_i0;
        s_eps_i(t, :)      = s_eps_i0;
    end
end

% Collect latent variables into struct
latents              = struct();
latents.tau_c        = tau_c;
latents.alpha_tau    = alpha_tau;
latents.sigma_dtau_c = sigma_dtau_c;
latents.tau_i        = tau_i;
latents.sigma_dtau_i = sigma_dtau_i;
latents.eps_c        = eps_c;
latents.alpha_eps    = alpha_eps;
latents.sigma_eps_c  = sigma_eps_c;
latents.s_eps_c      = s_eps_c;
latents.sigma_eps_i  = sigma_eps_i;
latents.s_eps_i      = s_eps_i;

end
