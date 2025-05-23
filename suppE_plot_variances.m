%%% PLOT VARIANCE SCATTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the Figure E4/Appendix E of the paper
% "A Measure of Trend Wage Inflation" (Almuzara, Audoly, Melcangi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
rng(2023)

% Set cut of wage inflation data
data_cut = 'industries_unweighted';

% Set directories
rng(2022);
addpath('functions');
data_path = [pwd filesep 'data' filesep data_cut];
res_path  = [pwd filesep 'results' filesep];
fig_path  = [pwd filesep 'figures' filesep];
tab_path  = [pwd filesep 'tables' filesep];


%% DATA

% Load data
wage_inflation   = readtable([data_path filesep 'wageinflation.csv']);
weights_data_cut = readtable([data_path filesep 'weights.csv']);
names_data_cut   = readtable([data_path filesep 'names.csv']);

% Extract data
date_str     = data_cut;                   % string to mark filenames (we use the YYYYMM, e.g., "202212")
dates        = wage_inflation{:, 1};       % Tx1 vector with datetimes
labels_short = names_data_cut{:, 1};       % Nx1 cell of strings with sector names
infla_disagg = wage_inflation{:, 2:end};   % TxN array with sectoral monthly annualized inflation rates
weights      = weights_data_cut{:, 2:end}; % TxN array with sectoral shares for each month (they add up to 1 row-wise)
share        = weights./sum(weights, 2);

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


%% COMPARISON OF TWV AND SAMPLE SIZES

load([res_path 'results_' date_str '.mat'])

% Set colors and markers
black   = [  0,   0,   0];
blue    = [  0,  85, 164]/255;
fig_fmt = 'epsc';

% Define figsize and font details
figsize = {'units', 'inches', 'position', [0 0 16 10]};
font    = {'fontname', 'times', 'fontsize', 22};

% Create scatterplot of inverse sqrt of sample sizes and TWV
YY    = TWV(:, 2);
XX    = 1./sqrt(sum(weights, 2));
coefs = [ones(length(XX), 1), XX]\YY;
fig0  = figure();
ax0   = axes();
plot0 = scatter(XX, YY, 75, blue, 'filled');

% Add extras
n_t_grid = [3000, 2500, 2000, 1500, 1000];
xlim_tmp = [1/sqrt(n_t_grid(1)), 1/sqrt(n_t_grid(end))];
xlim(ax0, xlim_tmp)
xticks(ax0, 1./sqrt(n_t_grid))
xticklabels(ax0, n_t_grid)
xlabel(ax0, '$n_t$', 'interpreter', 'latex')
ylim(ax0, [0.25, 0.40])
yticks(ax0, 0.25:0.05:0.40)
ylabel(ax0, '$\tilde{\sigma}_{\varepsilon, t}$', 'interpreter', 'latex')
hold('on')
plot(ax0, linspace(xlim_tmp(1), xlim_tmp(2), 1000), linspace(coefs(1)+coefs(2)*xlim_tmp(1), coefs(1)+coefs(2)*xlim_tmp(2), 1000), 'color', black, 'linewidth', 2, 'linestyle', ':')
hold('off')

% Tune ax handle
set(ax0, font{:})
set(ax0, 'box', 'on')
grid(ax0, 'on')

% Tune and save figure
set(fig0, figsize{:}); 
print(fig0, [fig_path 'FigE4_variance_samplesize'], ['-d' fig_fmt])
