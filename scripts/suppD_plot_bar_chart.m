%%% PLOT COMMON COMPONENT CONTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the Figure D6/Appendix C of the paper
% "A Measure of Trend Wage Inflation" (Almuzara, Audoly, Melcangi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
rng(2023)

% Set paths
res_path = [pwd filesep '..' filesep 'results' filesep];
fig_path = [pwd filesep '..' filesep 'figures' filesep];

% Set cut of wage inflation data
data_cuts = {'industries', 'occupation', 'race', 'education', 'age', 'gender', 'wage_quartile', 'region'};
cut_name  = {'Industry', 'Occupation', 'Race', 'Education', 'Age', 'Gender', 'Wage quartile', 'Region'};

% Upload data and variable names
X_common = NaN(length(data_cuts), 3, 3);
X_total  = NaN(length(data_cuts), 3, 3);
for i_data = 1:length(data_cuts)
    data_cut = data_cuts{i_data};
    date_str = data_cut;
    load([res_path 'results_' date_str '.mat'])
    X_common(i_data, :, :) = TWIn_cycles.change_c(:, [2, 3, 4])';
    X_total(i_data, :, :)  = TWIn_cycles.change(:, [2, 3, 4])';
end

% Compute dimensions
n_cut       = size(X_common, 1);
n_episode   = size(X_common, 3);
n_rest      = ceil(n_cut);
n_rest_init = ceil(n_cut);
n_rest_end  = ceil(n_cut);

% Fill horizontal and vertical axes
n_x           = n_rest_init+n_episode*n_cut+(n_episode-1)*n_rest+n_rest_end;
abscisa       = (1:n_x)';
ordinates_med = NaN(n_x, n_cut);
ordinates_low = NaN(n_x, n_cut);
ordinates_upp = NaN(n_x, n_cut);
ordinates_tot = NaN(n_x, n_cut);
for i_episode = 1:n_episode
    for i_cut = 1:n_cut
        ordinates_low(n_rest_init+(i_episode-1)*(n_cut+n_rest)+i_cut, i_cut) = X_common(i_cut, 1, i_episode);
        ordinates_med(n_rest_init+(i_episode-1)*(n_cut+n_rest)+i_cut, i_cut) = X_common(i_cut, 2, i_episode);
        ordinates_upp(n_rest_init+(i_episode-1)*(n_cut+n_rest)+i_cut, i_cut) = X_common(i_cut, 3, i_episode);
        ordinates_tot(n_rest_init+(i_episode-1)*(n_cut+n_rest)+i_cut, i_cut) = X_total(i_cut, 2, i_episode);
    end
end

% Compute tick values
tick_name  = {'2001'; '2008'; '2022'};
tick_val   = abscisa(~all(isnan(ordinates_med), 2));
tick_label = repmat({''}, n_cut*n_episode, 1);
for i_episode = 1:n_episode, tick_label{(i_episode-1)*n_cut+ceil(n_cut/2)} = tick_name{i_episode}; end

% Define colors
black       = [0, 0, 0];
line_color  = mat2cell(jet(n_cut), ones(n_cut, 1), 3);
plot_color  = line_color;
markers_low = repmat({'none'}, n_cut, 1);
markers_med = repmat({'diamond'}, n_cut, 1);
markers_upp = repmat({'none'}, n_cut, 1);
markers_tot = repmat({'square'}, n_cut, 1);

% Define figsize and font details
figsize = {'units', 'inches', 'position', [0 0 16 10]};
font    = {'fontname', 'times', 'fontsize', 22};
fig_fmt = 'epsc';

% Create figure
fig0     = figure();
ax0      = axes();
plot0    = plot(ax0, abscisa, ordinates_low);
set(plot0, {'color'}, plot_color)
set(plot0, 'linestyle', 'none')
set(plot0, {'marker'}, markers_low)
hold('on')
plot0    = plot(ax0, abscisa, ordinates_upp);
set(plot0, {'color'}, plot_color)
set(plot0, 'linestyle', 'none')
set(plot0, {'marker'}, markers_upp)

% Add lines connecting lower and upper bars
for i_episode = 1:n_episode
    for i_cut = 1:n_cut
        abscisa_tmp = n_rest_init+(i_episode-1)*(n_cut+n_rest)+i_cut;
        line(ax0, abscisa_tmp*ones(1, 2), [X_common(i_cut, 1, i_episode), X_common(i_cut, 3, i_episode)], 'color', line_color{i_cut}, 'linewidth', 5)
    end    
end
line(ax0, [1, n_x], [0, 0], 'color', black, 'linewidth', 0.15, 'linestyle', '--')

hold('on')
plot0    = plot(ax0, abscisa, ordinates_tot);
set(plot0, 'color', black)
set(plot0, 'linestyle', 'none')
set(plot0, {'marker'}, markers_tot)
set(plot0, 'markerfacecolor', black)
set(plot0, 'markersize', 8)
plot0    = plot(ax0, abscisa, ordinates_med);
set(plot0, {'color'}, plot_color)
set(plot0, 'linestyle', 'none')
set(plot0, {'marker'}, markers_med)
set(plot0, {'markerfacecolor'}, plot_color)
set(plot0, 'markersize', 8)

% Tune ax handle
set(ax0, font{:})
grid(ax0, 'off')

% Add extras
xlim(ax0, [min(abscisa), max(abscisa)])
xticks(ax0, tick_val)
xticklabels(ax0, tick_label)
legend(plot0, cut_name, 'location', 'Northwest', 'numColumns', 2, 'box', 'off')
ylim(ax0, [-3, 5])
ylabel(ax0, '(percent)')

% Tune and save figure
set(fig0, figsize{:});
print(fig0, [fig_path 'FigD6_episodes_change'], ['-d' fig_fmt]);
