%%% PLOT WAGE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the Figure 1/Section 2 of the paper
% "A Measure of Trend Wage Inflation" (Almuzara, Audoly, Melcangi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot wage inflation series
clearvars;
close all;

% Create directories
data_path = [pwd filesep 'data' filesep];
fig_path  = [pwd filesep 'figures' filesep];
if ~exist(fig_path, 'dir'), mkdir(fig_path); end

%% Data

df_m = readtable([data_path 'wage_inflation_monthly.csv']);
df_q = readtable([data_path 'wage_inflation_quarterly.csv']);

%% Figure options

fig_fmt = 'epsc';
figsize = {'units', 'inches', 'position', [0 0 16 10]};
font    = {'fontname', 'times', 'fontsize', 22};
line_width = 2.5;

% Define colors
white   = [  1,   1,   1];
black   = [  0,   0,   0];
gray    = [102, 102, 102]/255;
red     = [239,  65,  53]/255;
green   = [ 51, 153,  51]/255;
blue    = [  0,  85, 164]/255;


%% Plot aggregate wage inflation

% Periods to shade
ylim_tmp        = [-1, 10];
recession_dates = [datetime(2001,  3, 1), datetime(2001, 11, 1); ...
                   datetime(2007, 12, 1), datetime(2009,  6, 1); ...
                   datetime(2020,  2, 1), datetime(2020,  4, 1)];

f0  = figure();
ax0 = axes();

fill(ax0, [recession_dates'; flipud(recession_dates')], repmat(kron(ylim_tmp, ones(1, 2))', 1, size(recession_dates, 1)), black, 'facealpha', 0.1, 'linestyle', 'none')
hold('on')
start_date = datetime(2004,12,31);

% Monthly data
df = df_m(df_m.date>=start_date,:);
t = df.date;
y = [df.atlfed_wgt df.ahe];
p_m = plot(t,y,'LineWidth',line_width);
set(p_m, {'color'}, {blue;green});
set(p_m, {'linestyle'}, {'-';'-'});

% Quarterly data
hold('on');
df = df_q(df_q.date>=start_date,:);
t = df.date;
y = [df.eci];
p_q = plot(t,y,'LineWidth',line_width);
set(p_q, {'color'}, {red});
set(p_q, {'linestyle'}, {'-'});
hold('off');

% Format axis
xlim(ax0, [t(1) t(end)]);
ylim(ax0, [0, 9]);

% Tune ax handle
set(ax0,font{:});
grid(ax0,'off');

ylabel('Percent change from year earlier');

series_names = {'Atlanta Fed Wage Growth Tracker (within worker)', ...
    'Average Hourly Earnings', ...
    'Employment Cost Index (within job)'};

legend([p_m; p_q], series_names, 'Location', 'north');

set(f0, figsize{:}); 
print(f0, [fig_path 'Fig1_wage_inflation_aggregate'], ['-d' fig_fmt]);

