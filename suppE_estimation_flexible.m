%%% ESTIMATION OF WAGE INFLATION MODELS: CPS DATA %%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the robustness checks of the paper
% "A Measure of Wage Trend Inflation" (Almuzara, Audoly, Melcangi)
% 
% The script estimates a flexible model in the industry data cut and
% produces figures with trend estimates and decompositions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear memory
clear
close all
clc

% Set cut of wage inflation data
data_cut = 'industries';
    
% Set directories
addpath('functions')
data_path = [pwd filesep 'data' filesep data_cut];
res_path  = [pwd filesep 'results' filesep];
fig_path  = [pwd filesep 'figures' filesep];
tab_path  = [pwd filesep 'tables' filesep];
if ~exist(res_path, 'dir'), mkdir(res_path); end
if ~exist(fig_path, 'dir'), mkdir(fig_path); end
if ~exist(tab_path, 'dir'), mkdir(tab_path); end

% Determine tasks
estimation = true;

%%% DATA
% Load data
wage_inflation   = readtable([data_path filesep 'wageinflation.csv']);
weights_data_cut = readtable([data_path filesep 'weights.csv']);
names_data_cut   = readtable([data_path filesep 'names.csv']);

% Extract data
date_str     = [data_cut '_flexible'];     % string to mark filenames (we use the YYYYMM, e.g., "202212")
dates        = wage_inflation{:, 1};       % Tx1 vector with datetimes
labels_short = names_data_cut{:, 1};       % Nx1 cell of strings with sector names
infla_disagg = wage_inflation{:, 2:end};   % TxN array with sectoral monthly annualized inflation rates
share        = weights_data_cut{:, 2:end}; % TxN array with sectoral shares for each month (they add up to 1 row-wise)
share        = share./sum(share, 2);

% Create directory for figures and update path
if ~exist([fig_path date_str], 'dir')
    mkdir([fig_path date_str])
end
fig_archive_path = [fig_path date_str filesep];

% Set dimensions and tail probabilities for interval estimates
T      = length(dates);
n      = size(infla_disagg, 2);
signif = 1/6;

% Compute aggregate inflation rates
infla_agg = sum(share .* infla_disagg, 2);

% Set indexes for aggregation
agg_names    = labels_short';
n_agg        = length(labels_short);
agg_list     = num2cell(1:n_agg);
infla_aggreg = NaN(length(dates), n_agg);
for i_agg = 1:n_agg
    infla_aggreg(:, i_agg) = sum(share(:, agg_list{i_agg}) .* infla_disagg(:, agg_list{i_agg}), 2)./sum(share(:, agg_list{i_agg}), 2);
end

%%% ESTIMATION
% Set estimation dimensions
rng(2022)
settings               = struct();
settings.show_progress = true;
settings.n_draw        = 3000;
settings.n_burn        = 3000;
settings.n_thin        = 2;

% Set priors
prior         = struct();
prior.prec_MA = 1/10;
prior.nu_lam  = 60;
prior.s2_lam  = 0.0001;
prior.nu_gam  = 60;
prior.s2_gam  = 0.001;
prior.a_ps    = (1-1/48)*(120);
prior.b_ps    = (1/48)*(120);

if (estimation == true)
    % Perform estimation
    fprintf('Estimating model\n')
    settings.n_lags_c  = 12;
    settings.n_lags    = repmat(12, [n, 1]);
    settings.is_timeag = true(n, 1);
    output_TWIn         = estimate(infla_disagg, prior, settings);

    % Normalize common components
    output_TWIn.alpha_tau    = output_TWIn.alpha_tau .* permute(repmat(output_TWIn.sigma_dtau_c, 1, 1, n), [1 3 2]);
    output_TWIn.tau_c        = output_TWIn.tau_c ./ output_TWIn.sigma_dtau_c;
    output_TWIn.sigma_dtau_c = output_TWIn.sigma_dtau_c ./ output_TWIn.sigma_dtau_c;
    output_TWIn.alpha_eps    = output_TWIn.alpha_eps .* permute(repmat(output_TWIn.sigma_eps_c, 1, 1, n), [1 3 2]);
    output_TWIn.eps_c        = output_TWIn.eps_c ./ output_TWIn.sigma_eps_c;
    output_TWIn.sigma_eps_c  = output_TWIn.sigma_eps_c ./ output_TWIn.sigma_eps_c;

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
    TWIn_sector_part = median(repmat(share, 1, 1, n_draw).*trend_sector_draws, 3);
    TWIn_agg_part    = NaN(T, n_agg);
    TWIn_agg_share   = NaN(T, n_agg);
    for i_agg = 1:n_agg
        TWIn_agg_part(:, i_agg)  = sum(TWIn_sector_part(:, agg_list{i_agg}), 2);
        TWIn_agg_share(:, i_agg) = sum(share(:, agg_list{i_agg}), 2);
    end
    trend_sector         = median(trend_sector_draws, 3);
    trend_sector_part    = mean(repmat(share, 1, 1, n_draw).*trend_sector_draws, 3);
    
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
    save([res_path sprintf('results_%s.mat', date_str)], '-struct', 'output')
    
end

% Upload results
load([res_path 'results_' date_str '.mat'])

%%% FIGURES
% Define subsample for recent data
subsample        = (year(dates) >= 2019);
sample_line      = (year(dates) >= 1990);
line_ticks       = datetime((1998:4:2023), 1, 1);
trend_c_norm     = mean(TWIn_c(year(dates) <= 2019 & year(dates) >= 2017, 2));
trend_i_norm     = mean(TWIn_i(year(dates) <= 2019 & year(dates) >= 2017, 2));
trend_agg_c_norm = mean(TWIn_agg_c_part(year(dates) <= 2019 & year(dates) >= 2017, :), 1);
trend_agg_i_norm = mean(TWIn_agg_i_part(year(dates) <= 2019 & year(dates) >= 2017, :), 1);

% Define colors and figure format
white   = [  1,   1,   1];
black   = [  0,   0,   0];
red     = [239,  65,  53]/255;
blue    = [  0,  85, 164]/255;
dec_color = {blue; red};
agg_color = mat2cell(jet(n_agg), ones(n_agg, 1), 3);
fig_fmt = 'epsc';

% Define figsize and font details
figsize       = {'units', 'inches', 'position', [0 0 16 10]};
figsize_panel = {'units', 'inches', 'position', [0 0 14 14]};
font          = {'fontname', 'times', 'fontsize', 22};
font_panel    = {'fontname', 'times', 'fontsize', 30};

% Periods to shade
recession_dates = [datetime(2001,  3, 1), datetime(2001, 11, 1); ...
                   datetime(2007, 12, 1), datetime(2009,  6, 1); ...
                   datetime(2020,  2, 1), datetime(2020,  4, 1)];
       
% Plot trend with recession bands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylim_tmp = [0, 9];
fig0     = figure();
ax0      = axes();
fill(ax0, [recession_dates'; flipud(recession_dates')], repmat(kron(ylim_tmp, ones(1, 2))', 1, size(recession_dates, 1)), 0*white, 'facealpha', 0.1, 'linestyle', 'none')
hold('on')
fill(ax0, [dates', fliplr(dates')], [TWIn(:, 1)', fliplr(TWIn(:, 3)')], ...
    blue, 'facealpha', 0.3, 'linestyle', 'none', 'linewidth', 0.1, 'edgecolor', blue);
hold('on')
plot0    = plot(ax0, dates, [infla_agg, TWIn(:, 2)]);
ylabel(ax0, '(percent)');
hold('off')

% Add extras
xlim(ax0, [dates(1), dates(end)])
ylim(ax0, ylim_tmp)
xticks(ax0, datetime(year(dates(1)):5:year(dates(end)), 1, 1))

% Tune ax handle
set(ax0, font{:})
grid(ax0, 'off')

% Tune plot
legend(plot0, {'Wage growth (YoY)', 'Trend Wage Inflation'}, 'location', 'best', font{:})
set(plot0, {'linewidth'}, {1.7; 2.5})
set(plot0, {'color'}, {black; blue})
set(plot0, {'linestyle'}, {'-';  '-'})

% Tune and save figure
set(fig0, figsize{:}); 
print(fig0, [fig_archive_path 'FigE3_trend'], ['-d' fig_fmt])

% Plot sector-specific/common decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig0  = figure();
ax0   = axes();
fill(ax0, [dates(sample_line)', fliplr(dates(sample_line)')], [TWIn_i(sample_line, 1)'-trend_i_norm, fliplr(TWIn_i(sample_line, 3)'-trend_i_norm)], ...
    dec_color{1}, 'facealpha', 0.2, 'linestyle', 'none', 'linewidth', 0.1, 'edgecolor', dec_color{1});
hold('on')
fill(ax0, [dates(sample_line)', fliplr(dates(sample_line)')], [TWIn_c(sample_line, 1)'-trend_c_norm, fliplr(TWIn_c(sample_line, 3)'-trend_c_norm)], ...
    dec_color{2}, 'facealpha', 0.2, 'linestyle', 'none', 'linewidth', 0.1, 'edgecolor', dec_color{2});
hold('on')
plot0 = plot(ax0, dates(sample_line), [TWIn(sample_line, 2)-trend_i_norm-trend_c_norm, ...
    TWIn_i(sample_line, 2)-trend_i_norm, ...
    TWIn_c(sample_line, 2)-trend_c_norm]);
ylabel(ax0, 'Cumulative change from average over 2017-2019 (percent)');

% Add extras
xlim(ax0, [dates(find(sample_line, 1, 'first')), dates(end)])
xticks(ax0, line_ticks)

% Tune ax handle
set(ax0, font{:})
grid(ax0, 'off')

% Tune plot
legend(plot0, {'Trend Wage Inflation', 'Sector-specific', 'Common'}, 'location', 'best')
set(plot0, {'linewidth'}, {1.5; 2; 2})
set(plot0, {'linestyle'}, {'--'; '-'; '-'})
set(plot0, {'color'}, [black; dec_color])

% Tune and save figure
set(fig0, figsize{:}); 
print(fig0, [fig_archive_path 'FigE3_decomposition'], ['-d' fig_fmt])