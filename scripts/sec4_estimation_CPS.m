%%% ESTIMATION OF WAGE INFLATION MODELS: CPS DATA %%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the main empirical results of the paper
% "A Measure of Wage Trend Inflation" (Almuzara, Audoly, Melcangi)
% 
% The script estimates the model in different data cuts (industry,
% occupation, education, region, wage quartile, age, gender, race) and
% produces figures with trend estimates and decompositions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear memory
clear
close all
clc

% Set directories
addpath([pwd filesep '..' filesep 'functions']);
res_path  = [pwd filesep '..' filesep 'results' filesep];
fig_path  = [pwd filesep '..' filesep 'figures' filesep];
tab_path  = [pwd filesep '..' filesep 'tables' filesep];
if ~exist(res_path, 'dir'), mkdir(res_path); end
if ~exist(fig_path, 'dir'), mkdir(fig_path); end
if ~exist(tab_path, 'dir'), mkdir(tab_path); end

% Determine tasks
estimation = true;

% Set cut of wage inflation data
data_cuts = {'industries', 'occupation', 'education', 'region', 'wage_quartile', 'age', 'gender', 'race'};

for i_data = 1:length(data_cuts)

% Select data cut and set random number seed
data_cut  = data_cuts{i_data};
data_path = [pwd filesep '..' filesep 'data' filesep data_cut];
rng(2022)
close all

%%% DATA
% Load data
wage_inflation   = readtable([data_path filesep 'wageinflation.csv']);
weights_data_cut = readtable([data_path filesep 'weights.csv']);
names_data_cut   = readtable([data_path filesep 'names.csv']);

% Extract data
date_str     = data_cut;                   % string to mark filenames (we use the YYYYMM, e.g., "202212")
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
dates  = datetime(1997, (1:T)', 1);

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

% Set number of MA lags and time-aggregated sectors
n_lags    = repmat(3, [n, 1]);
is_timeag = true(n, 1);

if (estimation == true)
    % Perform estimation
    fprintf('Estimating model\n')
    settings.n_lags    = n_lags;
    settings.is_timeag = is_timeag;
    output_TWIn        = estimate(infla_disagg, prior, settings);

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
    
    % Cyclical contributions
    quants        = [0.05, 0.16, 0.50, 0.84, 0.95];
    dates_start   = [find(dates == datetime(2001, 03, 01)), ...
                     find(dates == datetime(2007, 12, 01)), ...
                     find(dates == datetime(2020, 06, 01))]';
    dates_end     = [find(dates == datetime(2001, 11, 01)), ...
                     find(dates == datetime(2009, 06, 01)), ...
                     find(dates == datetime(2022, 02, 01))]';    
    n_cycles      = size(dates_start, 1);
    change_TWIn   = NaN(n_cycles, 5);
    change_TWIn_c = NaN(n_cycles, 5);
    change_TWIn_i = NaN(n_cycles, 5);
    fract_TWIn_c  = NaN(n_cycles, 5);
    fract_TWIn_i  = NaN(n_cycles, 5);
    for j = 1:n_cycles
        TWIn_aux            = TWIn_draws(dates_end(j), :)   - TWIn_draws(dates_start(j), :);
        TWIn_c_aux          = TWIn_c_draws(dates_end(j), :) - TWIn_c_draws(dates_start(j), :);
        TWIn_i_aux          = TWIn_i_draws(dates_end(j), :) - TWIn_i_draws(dates_start(j), :);
        change_TWIn(j, :)   = quantile(TWIn_aux, quants, 2);
        change_TWIn_c(j, :) = quantile(TWIn_c_aux, quants, 2);
        change_TWIn_i(j, :) = quantile(TWIn_i_aux, quants, 2);
        fract_TWIn_c(j, :)  = quantile(TWIn_c_aux./TWIn_aux, quants, 2);
        fract_TWIn_i(j, :)  = quantile(TWIn_i_aux./TWIn_aux, quants, 2);
    end
            
    % Save results
    output                         = struct();
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
    output.TWIn_cycles             = struct();
    output.TWIn_cycles.dates       = [dates_start, dates_end];
    output.TWIn_cycles.change      = change_TWIn;
    output.TWIn_cycles.change_c    = change_TWIn_c;
    output.TWIn_cycles.change_i    = change_TWIn_i;
    output.TWIn_cycles.fract_c     = fract_TWIn_c;
    output.TWIn_cycles.fract_i     = fract_TWIn_i;
    save([res_path sprintf('results_%s.mat', date_str)], '-struct', 'output') 
end

% Upload estimates
load([res_path 'results_' date_str '.mat'])

%%% FIGURES
% Define subsample for recent data
subsample    = (year(dates) >= 2019);
sample_line  = (year(dates) >= 1990);
line_ticks   = datetime((1998:4:2023), 1, 1);
trend_c_norm = mean(TWIn_c(year(dates) <= 2019 & year(dates) >= 2017, 2));
trend_i_norm = mean(TWIn_i(year(dates) <= 2019 & year(dates) >= 2017, 2));

% Define colors and figure format
white   = [  1,   1,   1];
black   = [  0,   0,   0];
red     = [239,  65,  53]/255;
blue    = [  0,  85, 164]/255;
golden  = [ 51, 153,  51]/255;
dec_color = {blue; red};
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
print(fig0, [fig_archive_path 'Fig2_trend'], ['-d' fig_fmt])

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
print(fig0, [fig_archive_path 'Fig2_decomposition'], ['-d' fig_fmt])

% Plot sectoral details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_aux = (year(dates) >= 2017) & (year(dates) <= 2019);
for i = 1:n

    % Create figure
    fig0 = figure();
    ax0  = axes();
    yyaxis(ax0, 'right')
    area(ax0, dates, outlier_ind(:, i), 'facecolor', blue, 'facealpha', 0.25, 'edgecolor', white)
    ylabel(ax0, 'Outlier Probability')
    ylim(ax0, [0 1])
    yyaxis(ax0, 'left')
    hold on
    plot(ax0, dates, infla_disagg(:, i), 'color', black, 'linestyle', '-', 'linewidth', 1);
    plot(ax0, dates, trend_sector(:, i), ...
        'color', blue, 'linestyle', '-', 'linewidth', 2)
    plot(ax0, dates, trend_sector_c_med(:, i) - mean(trend_sector_c_med(t_aux, i) - trend_sector(t_aux, i)), ...
        'color', red, 'linestyle', '-', 'linewidth', 2);
    ylabel(ax0, 'Nominal wage growth (percent)')
    hold off

    % Add extras
    xlim(ax0, [dates(1), dates(end)])
    xticks(ax0, datetime(year(dates(1)):5:year(dates(end)), 1, 1))
    title(ax0, labels_short{i}, font{:})

    % Tune ax handle
    set(ax0, font_panel{:})
    set(ax0, 'box', 'on')
    grid(ax0, 'off')
    ax0.YAxis(1).Color = black;
    ax0.YAxis(2).Color = black;

    % Tune plot
    legend(ax0, {'Sectoral nominal wage growth (YoY)', 'Total trend', 'Common trend', 'Outlier probability (right axis)'}, 'location', 'best', font_panel{:})

    % Tune and save figure
    set(fig0, figsize_panel{:}); 
    print(fig0, [fig_archive_path sprintf('FigD1_sector_trend_%d', i)], ['-d' fig_fmt])

end

% Plot sectoral details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n

    % Create figure
    fig0 = figure();

    % Plot alpha_tau
    sign_tau = sign( squeeze(median(param_TWIn.alpha_tau(end, i, :))) );
    ax0      = subplot(2, 2, 1);
    plot0    = plot(ax0, dates, sign_tau*squeeze(param_TWIn.alpha_tau(:, i, :)).*param_TWIn.sigma_dtau_c, ...
        'color', black, 'linewidth', 2);
    xlim(ax0, [dates(1), dates(end)])
    ylabel(ax0, '$\alpha_{\tau, i}\sigma_{\Delta\tau, c}$', 'interpreter', 'latex')
    set(ax0, font{:})
    grid(ax0, 'on')
    set(plot0, {'linestyle'}, {':'; '-'; ':'})

    % Plot alpha_eps
    sign_eps = sign( squeeze(median(param_TWIn.alpha_eps(end, i, :))) );
    ax0      = subplot(2, 2, 2);
    plot0    = plot(ax0, dates, sign_eps*squeeze(param_TWIn.alpha_eps(:, i, :)).*param_TWIn.sigma_eps_c, ...
        'color', black, 'linewidth', 2);
    xlim(ax0, [dates(1), dates(end)])
    ylabel(ax0, '$\alpha_{\varepsilon, i}\sigma_{\varepsilon, c}$', 'interpreter', 'latex')    
    set(ax0, font{:})
    grid(ax0, 'on')
    set(plot0, {'linestyle'}, {':'; '-'; ':'})

    % Plot sigma_dtau_i
    ax0   = subplot(2, 2, 3);
    plot0 = plot(ax0, dates, squeeze(param_TWIn.sigma_dtau_i(:, i, :)), ...
        'color', black, 'linewidth', 2);
    xlim(ax0, [dates(1), dates(end)])
    ylabel(ax0, '$\sigma_{\Delta\tau, i}$', 'interpreter', 'latex')
    set(ax0, font{:})
    grid(ax0, 'on')
    set(plot0, {'linestyle'}, {':'; '-'; ':'})

    % Plot sigma_eps_i
    ax0   = subplot(2, 2, 4);
    plot0 = plot(ax0, dates, squeeze(param_TWIn.sigma_eps_i(:, i, :)), ...
        'color', black, 'linewidth', 2);
    xlim(ax0, [dates(1), dates(end)])
    ylabel(ax0, '$\sigma_{\varepsilon, i}$', 'interpreter', 'latex')
    set(ax0, font{:})
    grid(ax0, 'on')
    set(plot0, {'linestyle'}, {':'; '-'; ':'})

    % Tune and save figure
    sgtitle(fig0, labels_short{i}, font_panel{:})
    set(fig0, figsize_panel{:}); 
    print(fig0, [fig_archive_path sprintf('FigD2_sector_tvp_%d', i)], ['-d' fig_fmt])

end

%%% TABLES
% Tabulate cyclical contributions
var_names  = {'Start', 'End', 'Variable', 'Q05', 'Q16', 'Q50', 'Q84', 'Q95'};
row_names  = {'Change total', 'Change common', 'Change idiosyncratic', 'Fraction common', 'Fraction idiosyncratic'}';
n_cycles   = size(TWIn_cycles.dates, 1);
n_row      = length(row_names);
mat_cell   = NaN(n_cycles*n_row, 5);
for j = 1:n_cycles, mat_cell((j-1)*n_row+(1:n_row), :) = [TWIn_cycles.change(j, :); TWIn_cycles.change_c(j, :); TWIn_cycles.change_i(j, :); TWIn_cycles.fract_c(j, :); TWIn_cycles.fract_i(j, :)]; end
tab_cell   = [cell(n_cycles*n_row, 2), repmat(row_names, n_cycles, 1), num2cell(round(mat_cell, 2))];
for j = 1:n_cycles
    tab_cell{(j-1)*n_row+1, 1} = dates(TWIn_cycles.dates(j, 1));
    tab_cell{(j-1)*n_row+1, 2} = dates(TWIn_cycles.dates(j, 2));
end
tab_cycles = cell2table(tab_cell, 'variableNames', var_names);

% Write to excel
writetable(tab_cycles, [tab_path 'Tab3_episodes_contribution.xlsx'], 'sheet', date_str, 'writerownames', false);

end