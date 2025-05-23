%%% VALIDATION EXERCISE: ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script and the next replicate the validation study of the paper
% "A Measure of Trend Wage Inflation" (Almuzara, Audoly, Melcangi)
%
% This script estimates our wage inflation model in pseudo-real time in
% parallel across samples (it requires the parallel computing toolbox).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear memory
clear
clc
close all

% Set cut of wage inflation data
data_cut = 'industries';

% Set directories
addpath('functions');
data_path = [pwd filesep 'data' filesep data_cut];
res_path  = [pwd filesep 'results' filesep 'realtime' filesep];
if ~exist(res_path, 'dir'), mkdir(res_path); end


%% DATA

% Load data
wage_inflation   = readtable([data_path filesep 'wageinflation.csv']);
weights_data_cut = readtable([data_path filesep 'weights.csv']);
names_data_cut   = readtable([data_path filesep 'names.csv']);

% Extract data
dates        = wage_inflation{:, 1};       % Tx1 vector with datetimes
labels_short = names_data_cut{:, 1};       % Nx1 cell of strings with sector names
infla_disagg = wage_inflation{:, 2:end};   % TxN array with sectoral monthly annualized inflation rates
share        = weights_data_cut{:, 2:end}; % TxN array with sectoral shares for each month (they add up to 1 row-wise)
share        = share./sum(share, 2);

% Set dimensions and tail probabilities for interval estimates
T      = length(dates);
n      = size(infla_disagg, 2);
signif = 1/6;

% Compute aggregate inflation rates
infla_agg = sum(share .* infla_disagg, 2);

% Set indexes for aggregation
agg_names    = labels_short';
n_agg        = length(labels_short);
agg_list_par = num2cell(1:n_agg);
infla_aggreg = NaN(length(dates), n_agg);
for i_agg = 1:n_agg
    infla_aggreg(:, i_agg) = sum(share(:, agg_list_par{i_agg}) .* infla_disagg(:, agg_list_par{i_agg}), 2)./sum(share(:, agg_list_par{i_agg}), 2);
end


%% ESTIMATION

% Set estimation dimensions
settings               = struct();
settings.show_progress = true;
settings.n_draw        = 3000;
settings.n_burn        = 3000;
settings.n_thin        = 2;

% Set theta/lambda/gamma/ps priors
prior_par         = struct();
prior_par.prec_MA = 1/10;
prior_par.nu_lam  = 60;
prior_par.s2_lam  = 0.0001;
prior_par.nu_gam  = 60;
prior_par.s2_gam  = 0.001;
prior_par.a_ps    = (1-1/48)*(120);
prior_par.b_ps    = (1/48)*(120);

% Set number of MA lags and time-aggregated sectors
n_lags             = repmat(3, [n, 1]);
is_timeag          = true(n, 1);
settings.n_lags    = n_lags;
settings.is_timeag = is_timeag;

% Perform estimation in parallel
n_work   = 300;
par_pool = parpool(n_work);
parfor i_work = 1:n_work

    rng(2022)
    
    % Set broadcast variables
    prior    = prior_par;
    agg_list = agg_list_par;

    % Define dataset for parallel estimation
    infla_disagg_par = infla_disagg;
    infla_disagg_par((end-(i_work-1)):end, :) = NaN;

    % Perform estimation
    output_TWIn = estimate(infla_disagg_par, prior, settings);

    % Recover number of draws
    n_draw = settings.n_draw;

    % Compute trend
    trend_sector_c_draws = squeeze(output_TWIn.alpha_tau.*permute(repmat(output_TWIn.tau_c, 1, 1, n), [1 3 2]));
    trend_sector_i_draws = output_TWIn.tau_i;
    trend_sector_draws   = trend_sector_c_draws + trend_sector_i_draws;
    TWIn_c_draws         = squeeze(sum(repmat(share, 1, 1, n_draw).*trend_sector_c_draws, 2));
    TWIn_i_draws         = squeeze(sum(repmat(share, 1, 1, n_draw).*trend_sector_i_draws, 2));
    TWIn_draws           = TWIn_c_draws + TWIn_i_draws;
    TWIn_c               = quantile(TWIn_c_draws, [signif, 0.5, 1-signif], 2);
    TWIn_i               = quantile(TWIn_i_draws, [signif, 0.5, 1-signif], 2);
    TWIn                 = quantile(TWIn_draws, [signif, 0.5, 1-signif], 2);

    % Compute decomposition of sectoral contributions
    TWIn_sector_part  = median(repmat(share, 1, 1, n_draw).*trend_sector_draws, 3);
    TWIn_agg_part     = NaN(T, n_agg);
    TWIn_agg_share    = NaN(T, n_agg);
    for i_agg = 1:n_agg
        TWIn_agg_part(:, i_agg)  = sum(TWIn_sector_part(:, agg_list{i_agg}), 2);
        TWIn_agg_share(:, i_agg) = sum(share(:, agg_list{i_agg}), 2);
    end
    trend_sector      = median(trend_sector_draws, 3);
    trend_sector_part = mean(repmat(share, 1, 1, n_draw).*trend_sector_draws, 3);

    % Compute decomposition of sectoral contributions to common/idiosyncratic
    TWIn_sector_c_part = median(repmat(share, 1, 1, n_draw).*trend_sector_c_draws, 3);
    TWIn_agg_c_part    = NaN(T, n_agg);
    for i_agg = 1:n_agg
        TWIn_agg_c_part(:, i_agg) = sum(TWIn_sector_c_part(:, agg_list{i_agg}), 2);
    end
    TWIn_sector_i_part = median(repmat(share, 1, 1, n_draw).*trend_sector_i_draws, 3);
    TWIn_agg_i_part    = NaN(T, n_agg);
    for i_agg = 1:n_agg
        TWIn_agg_i_part(:, i_agg) = sum(TWIn_sector_i_part(:, agg_list{i_agg}), 2);
    end

    % Compute outlier probabilities
    outlier_ind = output_TWIn.s_eps_i;
    outlier_ind(outlier_ind == 1) = 0;
    outlier_ind(outlier_ind > 1)  = 1;
    outlier_ind = mean(outlier_ind, 3);

    % Save results
    output = struct();
    output.param_TWIn              = struct();
    output.param_TWIn.theta        = quantile(output_TWIn.theta, [signif, 0.5, 1-signif], 3);
    output.param_TWIn.lam_tau      = quantile(output_TWIn.lam_tau, [signif, 0.5, 1-signif], 2);
    output.param_TWIn.gam_dtau_c   = quantile(output_TWIn.gam_dtau_c, [signif, 0.5, 1-signif], 2);
    output.param_TWIn.gam_dtau_i   = quantile(output_TWIn.gam_dtau_i, [signif, 0.5, 1-signif], 2);
    output.param_TWIn.lam_eps      = quantile(output_TWIn.lam_eps, [signif, 0.5, 1-signif], 2);
    output.param_TWIn.gam_eps_c    = quantile(output_TWIn.gam_eps_c, [signif, 0.5, 1-signif], 2);
    output.param_TWIn.gam_eps_i    = quantile(output_TWIn.gam_eps_i, [signif, 0.5, 1-signif], 2);
    output.param_TWIn.alpha_tau    = quantile(output_TWIn.alpha_tau, [signif, 0.5, 1-signif], 3);
    output.param_TWIn.alpha_eps    = quantile(output_TWIn.alpha_eps, [signif, 0.5, 1-signif], 3);
    output.param_TWIn.sigma_dtau_c = quantile(output_TWIn.sigma_dtau_c, [signif, 0.5, 1-signif], 2);
    output.param_TWIn.sigma_dtau_i = quantile(output_TWIn.sigma_dtau_i, [signif, 0.5, 1-signif], 3);
    output.param_TWIn.sigma_eps_c  = quantile(output_TWIn.sigma_eps_c, [signif, 0.5, 1-signif], 2);
    output.param_TWIn.sigma_eps_i  = quantile(output_TWIn.sigma_eps_i, [signif, 0.5, 1-signif], 3);
    output.trend_sector_c_med      = median(trend_sector_c_draws, 3);
    output.trend_sector_i_med      = median(trend_sector_i_draws, 3);
    output.trend_sector            = trend_sector;
    output.TWIn_c                  = TWIn_c;
    output.TWIn_i                  = TWIn_i;
    output.TWIn                    = TWIn;
    output.trend_sector_part       = trend_sector_part;
    output.TWIn_agg_part           = TWIn_agg_part;
    output.TWIn_agg_share          = TWIn_agg_share;
    output.TWIn_agg_c_part         = TWIn_agg_c_part;
    output.TWIn_agg_i_part         = TWIn_agg_i_part;
    output.outlier_ind             = outlier_ind;

    % Compute last month of data
    dates_par = dates;
    last_date = dates_par(end) - calmonths((1 + (i_work-1)));
    date_str  = [data_cut '_' datestr(last_date, 'yyyymm')]; %#ok

    % Save results
    parsave([res_path sprintf('results_%s.mat', date_str)], output)

end 
delete(par_pool)

% Store estimation results to csv
load('results/results_industries.mat');
date = dates;
twin = TWIn(:,2);
common = TWIn_c(:, 2);
twin_series = table(date, twin, common);

writetable(twin_series, 'data/twin_series.csv');



%% LOCAL FUNCTIONS

function [] = parsave(fname, x)
    save(fname, '-struct', 'x')
end
