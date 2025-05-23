%%% PLOT WAGE AND LABOR MARKET TIGHTNESS MEASURES %%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the Figure 3/Section 4 of the paper
% "A Measure of Trend Wage Inflation" (Almuzara, Audoly, Melcangi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear memory
clear
close all
clc

% Set directories
rng(2022)
addpath('functions')
data_path = [pwd filesep 'data' filesep];
fig_path  = [pwd filesep 'figures' filesep];
if ~exist(fig_path, 'dir'), mkdir(fig_path); end
warning('off')

% Load data
data      = readtable([data_path 'labor_market_monthly.csv']);
for t = 3:size(data, 1), data{t, 'AWGT_MA'} = mean(data{(t-2):t, 'AWGT'}); end
dates     = datetime(1948, 1:size(data, 1), 1);
subsample = (year(dates) > 2000);
dates     = dates(subsample);
UR        = 1-data{subsample, 'u_rate'};
VU        = data{subsample, 'tightness_v_over_u'};
VUE       = data{subsample, 'tightness_v_over_u_e'};
TWIn      = data{subsample, 'TWIn'};
AWGT      = data{subsample, 'AWGT_MA'};
ECI       = data{subsample, 'ECI'};

% Interpolate ECI
keep  = ~isnan(ECI);
x     = (1:size(ECI))';
x_in  = x(keep);
x_out = x(~keep);
y_in  = ECI(keep);
y_out = interp1(x_in, y_in, x_out);
ECI(~keep) = y_out;

% Create data matrix and variable names
X = [VUE, TWIn, AWGT, ECI];
var_names = {'Labor market tightness', 'Trend wage inflation', 'Atlanta Fed wage growth tracker', 'Employment cost index'};

% Define colors and figure format
black   = [  0,   0,   0];
red     = [239,  65,  53]/255;
green   = [ 51, 153,  51]/255;
blue    = [  0,  85, 164]/255;
golden  = [218, 165,  32]/255;
fig_fmt = 'epsc';
figsize_panel = {'units', 'inches', 'position', [0 0 14 14]};
font_panel    = {'fontname', 'times', 'fontsize', 30};

% Periods to shade
ylim_tmp        = [0, 10];
recession_dates = [datetime(2001,  3, 1), datetime(2001, 11, 1); ...
                   datetime(2007, 12, 1), datetime(2009,  6, 1); ...
                   datetime(2020,  2, 1), datetime(2020,  4, 1)];

% Plot subfigures by subperiod
base_periods = [...
    find(year(dates) == 2001 & month(dates) == 03, 1, 'first'), ...
    find(year(dates) == 2007 & month(dates) == 12, 1, 'first'), ...
    find(year(dates) == 2020 & month(dates) == 06, 1, 'first')];

for j = 1:3
    if j == 1, n_first = 2; else; n_first = 12; end
    n_month = 33;

    %%% Plot normalized series
    t_tmp  = base_periods(j);
    X_norm = X./X(t_tmp, :);    
    x_plot = dates(t_tmp+(-n_first:n_month-1));
    y_plot = X_norm(t_tmp+(-n_first:n_month-1), :);
    y_MA   = y_plot;
    for jj = 2:(size(y_plot, 1)-1), y_MA(jj, :) = mean(y_plot((jj-1):(jj+1), :), 1, 'omitnan'); end
    y_plot(:, 1) = y_MA(:, 1);
    fig0   = figure();
    ax0    = axes();
    fill(ax0, [recession_dates'; flipud(recession_dates')], repmat(kron([min(y_plot(:))-0.1, max(y_plot(:))+0.1], ones(1, 2))', 1, size(recession_dates, 1)), black, 'facealpha', 0.1, 'linestyle', 'none')
    hold('on')
    plot0  = plot(ax0, x_plot, y_plot);
    hold('off')

    % Add extras
    xlim(ax0, [dates(t_tmp-n_first), dates(t_tmp+n_month-1)])
    ylim(ax0, [min(y_plot(:))-0.05, max(y_plot(:))+0.05])
    ylabel(ax0, sprintf('Labor market tightness and wage inflation (%dm%d=1)', year(dates(t_tmp)), month(dates(t_tmp))))
    
    % Tune ax handle
    set(ax0, font_panel{:})
    grid(ax0, 'on')    
    if (j == 1)
        legend(plot0, var_names(1:end-1), 'location', 'best')
    else
        legend(plot0, var_names, 'location', 'best')
    end
    set(plot0, {'color'}, {red; blue; 0.9*golden; 0.6*green+0.4*blue})
    set(plot0, {'linewidth'}, {2.5; 3; 2; 2})
    set(plot0, {'linestyle'}, {'-';' :'; '-.'; '-.'})
    
    % Tune and save figure
    set(fig0, figsize_panel{:})
    print(fig0, [fig_path 'Fig3_labor_market_' num2str(year(dates(t_tmp)))], ['-d' fig_fmt])
    
end
