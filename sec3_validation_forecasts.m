%%% VALIDATION EXERCISE: COMPARISON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This and the previous script replicate the validation study of the paper
% "A Measure of Trend Wage Inflation" (Almuzara, Audoly, Melcangi)
%
% This script compares the forecasting performance of our wage inflation 
% model in pseudo-real time against the alternatives. It also implements
% tests of forecasting dominance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear memory
clear
close all
clc

% Set cut of wage inflation data
data_cut = 'industries';

% Set directories
rng(2022)
addpath('functions')
data_path = [pwd filesep 'data' filesep data_cut];
res_path  = [pwd filesep 'results' filesep 'realtime' filesep];
tab_path  = [pwd filesep 'tables' filesep];
if ~exist(tab_path, 'dir'), mkdir(tab_path); end


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

% Load wage data (AWGT/AHE/ECI)
data        = readtable([pwd filesep 'data' filesep 'wage_growth_forecasts.csv']);
sample_full = (data{:, 'DATE'} >= min(dates));
AWGT_full   = data{sample_full, 'AWGT'};

% Set numerical parameters for reading parallel output
base_year  = 1998;
base_month = 9;
n_res      = 300;
month_gap  = 1;

% Define map from year-month to string
ym2str = @(y, m) sprintf('%4d%02d', y, m);

% Allocate horizons and lags in regression forecasts
horizon = [3, 6, 9, 12];
n_h     = length(horizon);
n_lag   = 0;

% Set number of boostrap replicas
n_boot     = 10000;
block_size = 7;


%% FORECAST PERFORMANCE (MEDIAN WAGE GROWTH)

%%% Computation of forecasts
% Allocate results to matrix
dates_fc          = NaT(n_res, 1);
TWIn_fc_RW        = NaN(n_res, n_h);
data_fc_RW        = NaN(n_res, n_h);
wage_realized     = NaN(n_res, n_h);
TWIn_fc_RW_sec    = NaN(n_res, n_h, n);
data_fc_RW_sec    = NaN(n_res, n_h, n);
wage_realized_sec = NaN(n_res, n_h, n);
i_vintage         = NaN(n_res, 1);
TWIn_vintage      = NaN(T, n_res);
for i_res = 1:n_res
    
    % Determine date index
    date_tmp = datetime(base_year, base_month+month_gap*(i_res-1), 1);
    y_tmp    = year(date_tmp);
    m_tmp    = month(date_tmp);
    i_date   = find(year(dates) == y_tmp & month(dates) == m_tmp, 1, 'first');
    dates_fc(i_res) = date_tmp;

    % Upload predictors
    res          = load([res_path 'results_' data_cut '_' ym2str(y_tmp, m_tmp) '.mat']);
    TWIn_tmp     = res.TWIn(1:i_date, 2);
    data_tmp     = infla_agg(1:i_date);
    TWIn_sec_tmp = res.trend_sector(1:i_date, :);
    data_sec_tmp = infla_disagg(1:i_date, :);

    % Save vintages
    i_vintage(i_res)       = i_date;
    TWIn_vintage(:, i_res) = res.TWIn(:, 2);
    
    % Compute forecasts and realized inflation
    for i_h = 1:n_h
        % Set horizon
        h_tmp = horizon(i_h);
        
        % Compute random-walk forecasts
        TWIn_fc_RW(i_res, i_h) = mean([TWIn_tmp((i_date-12+1+h_tmp):i_date); repmat(TWIn_tmp(i_date), [h_tmp, 1])]);
        data_fc_RW(i_res, i_h) = data_tmp(i_date);
                
        % Compute realized wage inflation
        if (i_date+h_tmp<T), wage_realized(i_res, i_h) = infla_agg(i_date+h_tmp); end
        
        for i = 1:n
            % Compute random-walk forecasts
            TWIn_fc_RW_sec(i_res, i_h, i) = mean([TWIn_sec_tmp((i_date-12+1+h_tmp):i_date, i); repmat(TWIn_sec_tmp(i_date, i), [h_tmp, 1])]);
            data_fc_RW_sec(i_res, i_h, i) = data_sec_tmp(i_date, i);
                        
            % Compute realized wage inflation
            if (i_date+h_tmp<T), wage_realized_sec(i_res, i_h, i) = infla_disagg(i_date+h_tmp, i); end
        end
    end
    
end

%%% RMSE calculations
% Preallocate RMSE
RMSE_agg      = NaN(n_h, 2);
RMSE_sector   = NaN(n_h, 2, n);
stderr_agg    = NaN(n_h, 3);
stderr_sector = NaN(n_h, 3, n);
pval_agg      = NaN(n_h, 1);
pval_sector   = NaN(n_h, n);

% Compute RMSE and bootstrap standard errors for aggregate
for i_h = 1:n_h
    
    % Compute forecast error samples
    error_TWIn = wage_realized(:, i_h)-TWIn_fc_RW(:, i_h);
    error_data = wage_realized(:, i_h)-data_fc_RW(:, i_h);
    nonmiss    = ~isnan(error_TWIn) & ~isnan(error_data);
    n_RMSE     = nnz(nonmiss);
    error_TWIn = error_TWIn(nonmiss);
    error_data = error_data(nonmiss);
    T_boot     = nnz(nonmiss);
    T_boot     = block_size*floor(T_boot/block_size);    
    
    % Compute estimates RMSE
    RMSE_TWIn = sqrt( mean(error_TWIn.^2, 'omitnan'));
    RMSE_data = sqrt( mean(error_data.^2, 'omitnan'));
    RMSE_agg(i_h, :) = [RMSE_TWIn, RMSE_data];
    
    % Compute bootstrap standard errors
    RMSE_boot = NaN(n_boot, 3);
    for i_boot = 1:n_boot
        idx_tmp         = datasample(((1:(T_boot-block_size+1))')+(0:(block_size-1)), T_boot/block_size);
        error_TWIn_boot = error_TWIn(idx_tmp(:));
        error_data_boot = error_data(idx_tmp(:));
        RMSE_boot(i_boot, :) = [sqrt( mean(error_TWIn_boot.^2, 'omitnan')), sqrt( mean(error_data_boot.^2, 'omitnan')), ...
            sqrt( mean(error_TWIn_boot.^2, 'omitnan'))-sqrt( mean(error_data_boot.^2, 'omitnan'))];
    end
    stderr_agg(i_h, :) = std(RMSE_boot, [], 1);
    
    % Compute p-values
    pval_agg(i_h, :) = 1-mean(RMSE_boot(:, 1) <= RMSE_boot(:, 2));
    
end

% Compute RMSE and bootstrap standard errors for each sector
for i = 1:n
    for i_h = 1:n_h
        
        % Compute forecast error samples
        error_TWIn = wage_realized_sec(:, i_h, i)-TWIn_fc_RW_sec(:, i_h, i);
        error_data = wage_realized_sec(:, i_h, i)-data_fc_RW_sec(:, i_h, i);
        nonmiss    = ~isnan(error_TWIn) & ~isnan(error_data);
        error_TWIn = error_TWIn(nonmiss);
        error_data = error_data(nonmiss);
        T_boot     = nnz(nonmiss);
        T_boot     = block_size*floor(T_boot/block_size);
        
        % Compute estimates RMSE
        RMSE_TWIn = sqrt( mean(error_TWIn.^2, 'omitnan'));
        RMSE_data = sqrt( mean(error_data.^2, 'omitnan'));
        RMSE_sector(i_h, :, i) = [RMSE_TWIn, RMSE_data];
        
        % Compute bootstrap standard errors
        RMSE_boot = NaN(n_boot, 3);
        for i_boot = 1:n_boot
            idx_tmp         = datasample(((1:(T_boot-block_size+1))')+(0:(block_size-1)), T_boot/block_size);
            error_TWIn_boot = error_TWIn(idx_tmp(:));
            error_data_boot = error_data(idx_tmp(:));
            RMSE_boot(i_boot, :) = [sqrt( mean(error_TWIn_boot.^2, 'omitnan')), sqrt( mean(error_data_boot.^2, 'omitnan')), ...
                sqrt( mean(error_TWIn_boot.^2, 'omitnan'))-sqrt( mean(error_data_boot.^2, 'omitnan'))];
        end
        stderr_sector(i_h, :, i) = std(RMSE_boot, [], 1);
        
        % Compute p-values
        pval_sector(i_h, i) = 1-mean(RMSE_boot(:, 1) <= RMSE_boot(:, 2));

    end
end

% Store results in table
col_names     = {'TWInRMSE', 'TWInSE', 'MWGRMSE', 'MWGSE', 'Diff', 'DiffSE', 'pval'};
row_names     = cell(1, n_h); for i_h = 1:n_h, row_names{i_h} = [num2str(horizon(i_h)) '-month']; end
tab_aggregate = array2table(round([RMSE_agg(:, 1), stderr_agg(:, 1), RMSE_agg(:, 2), stderr_agg(:, 2), ...
    RMSE_agg(:, 2)-RMSE_agg(:, 1), stderr_agg(:, 3), pval_agg], 4), ...
    'variableNames', col_names, 'rowNames', row_names);
disp('MEDIAN WAGE GROWTH (RANDOM-WALK) FORECAST')
disp(['No. of periods: ' num2str(n_RMSE)])
disp('Aggregate')
disp(tab_aggregate)
writetable(tab_aggregate, [tab_path 'Tab1_RMSEs.xlsx'], 'sheet', 'MWG RW Aggregate', 'writevariablenames', true, 'writerownames', true);
for i = 1:n
    tab_sectors = array2table(round([RMSE_sector(:, 1, i), stderr_sector(:, 1, i), RMSE_sector(:, 2, i), stderr_sector(:, 2, i), ...
        RMSE_sector(:, 2, i)-RMSE_sector(:, 1, i), stderr_sector(:, 3, i), pval_sector(:, i)], 4), ...
        'variableNames', col_names, 'rowNames', row_names);
    writetable(tab_sectors, [tab_path 'Tab1_RMSEs.xlsx'], 'sheet', ['MWG RW -' labels_short{i}(1:min(length(labels_short{i}), 20))], 'writevariablenames', true, 'writerownames', true);
    disp(labels_short{i})
    disp(tab_sectors)
end


%% FORECAST PERFORMANCE (12M MEDIAN WAGE GROWTH)

%%% Computation of forecasts
% Allocate results to matrix
dates_fc      = NaT(n_res, 1);
TWIn_fc_RW    = NaN(n_res, n_h);
data_fc_RW    = NaN(n_res, n_h);
wage_realized = NaN(n_res, n_h);
for i_res = 1:n_res
    
    % Determine date index
    date_tmp = datetime(base_year, base_month+month_gap*(i_res-1), 1);
    y_tmp    = year(date_tmp);
    m_tmp    = month(date_tmp);
    i_date   = find(year(dates) == y_tmp & month(dates) == m_tmp, 1, 'first');
    dates_fc(i_res) = date_tmp;

    % Upload predictors
    res      = load([res_path 'results_' data_cut '_' ym2str(y_tmp, m_tmp) '.mat']);
    TWIn_tmp = res.TWIn(1:i_date, 2);
    data_aux = infla_agg(1:i_date);
    data_tmp = NaN(length(data_aux), 1);
    for t = 1:i_date
        data_tmp(t) = mean(data_aux(max(1, (t-11)):t));
    end
    
    % Compute forecasts and realized inflation
    for i_h = 1:n_h
        % Set horizon
        h_tmp = horizon(i_h);
        
        % Compute random-walk forecasts
        TWIn_fc_RW(i_res, i_h) = mean([TWIn_tmp((i_date-12+1+h_tmp):i_date); repmat(TWIn_tmp(i_date), [h_tmp, 1])]);
        data_fc_RW(i_res, i_h) = data_tmp(i_date);
                
        % Compute realized wage inflation
        if (i_date+h_tmp<T), wage_realized(i_res, i_h) = infla_agg(i_date+h_tmp); end
    end
    
end

%%% RMSE calculations
% Preallocate RMSE
RMSE_agg   = NaN(n_h, 2);
stderr_agg = NaN(n_h, 3);
pval_agg   = NaN(n_h, 1);

% Compute RMSE and bootstrap standard errors for aggregate
for i_h = 1:n_h
    
    % Compute forecast error samples
    error_TWIn = wage_realized(:, i_h)-TWIn_fc_RW(:, i_h);
    error_data = wage_realized(:, i_h)-data_fc_RW(:, i_h);
    nonmiss    = ~isnan(error_TWIn) & ~isnan(error_data);
    n_RMSE     = nnz(nonmiss);
    error_TWIn = error_TWIn(nonmiss);
    error_data = error_data(nonmiss);
    T_boot     = nnz(nonmiss);
    T_boot     = block_size*floor(T_boot/block_size);    
    
    % Compute estimates RMSE
    RMSE_TWIn = sqrt( mean(error_TWIn.^2, 'omitnan'));
    RMSE_data = sqrt( mean(error_data.^2, 'omitnan'));
    RMSE_agg(i_h, :) = [RMSE_TWIn, RMSE_data];
    
    % Compute bootstrap standard errors
    RMSE_boot = NaN(n_boot, 3);
    for i_boot = 1:n_boot
        idx_tmp         = datasample(((1:(T_boot-block_size+1))')+(0:(block_size-1)), T_boot/block_size);
        error_TWIn_boot = error_TWIn(idx_tmp(:));
        error_data_boot = error_data(idx_tmp(:));
        RMSE_boot(i_boot, :) = [sqrt( mean(error_TWIn_boot.^2, 'omitnan')), sqrt( mean(error_data_boot.^2, 'omitnan')), ...
            sqrt( mean(error_TWIn_boot.^2, 'omitnan'))-sqrt( mean(error_data_boot.^2, 'omitnan'))];
    end
    stderr_agg(i_h, :) = std(RMSE_boot, [], 1);
    
    % Compute p-values
    pval_agg(i_h, :) = 1-mean(RMSE_boot(:, 1) <= RMSE_boot(:, 2));
    
end

% Store results in table
col_names     = {'TWInRMSE', 'TWInSE', 'MGW12mRMSE', 'MGW12mSE', 'Diff', 'DiffSE', 'pval'};
row_names     = cell(1, n_h); for i_h = 1:n_h, row_names{i_h} = [num2str(horizon(i_h)) '-month']; end
tab_aggregate = array2table(round([RMSE_agg(:, 1), stderr_agg(:, 1), RMSE_agg(:, 2), stderr_agg(:, 2), ...
    RMSE_agg(:, 2)-RMSE_agg(:, 1), stderr_agg(:, 3), pval_agg], 4), ...
    'variableNames', col_names, 'rowNames', row_names);
disp('12-MONTH MEDIAN WAGE GROWTH')
disp(['No. of periods: ' num2str(n_RMSE)])
disp('Aggregate')
disp(tab_aggregate)
writetable(tab_aggregate, [tab_path 'Tab1_RMSEs.xlsx'], 'sheet', 'MWG12M - Aggregate', 'writevariablenames', true, 'writerownames', true);


%% FORECAST PERFORMANCE (AWGT)

%%% Computation of forecasts
% Allocate results to matrix
dates_fc      = NaT(n_res, 1);
TWIn_fc_RW    = NaN(n_res, n_h);
data_fc_RW    = NaN(n_res, n_h);
wage_realized = NaN(n_res, n_h);
for i_res = 1:n_res
    
    % Determine date index
    date_tmp = datetime(base_year, base_month+month_gap*(i_res-1), 1);
    y_tmp    = year(date_tmp);
    m_tmp    = month(date_tmp);
    i_date   = find(year(dates) == y_tmp & month(dates) == m_tmp, 1, 'first');
    dates_fc(i_res) = date_tmp;

    % Upload predictors
    res      = load([res_path 'results_' data_cut '_' ym2str(y_tmp, m_tmp) '.mat']);
    TWIn_tmp = res.TWIn(1:i_date, 2);
    data_tmp = AWGT_full(1:i_date);
    
    % Compute forecasts and realized inflation
    for i_h = 1:n_h
        % Set horizon
        h_tmp = horizon(i_h);
        
        % Compute random-walk forecasts
        TWIn_fc_RW(i_res, i_h) = mean([TWIn_tmp((i_date-12+1+h_tmp):i_date); repmat(TWIn_tmp(i_date), [h_tmp, 1])]);
        data_fc_RW(i_res, i_h) = data_tmp(i_date);
                
        % Compute realized wage inflation
        if (i_date+h_tmp<T), wage_realized(i_res, i_h) = AWGT_full(i_date+h_tmp); end
    end
    
end


%%% RMSE calculations
% Preallocate RMSE
RMSE_agg   = NaN(n_h, 2);
stderr_agg = NaN(n_h, 3);
pval_agg   = NaN(n_h, 1);

% Compute RMSE and bootstrap standard errors for aggregate
for i_h = 1:n_h
    
    % Compute forecast error samples
    error_TWIn = wage_realized(:, i_h)-TWIn_fc_RW(:, i_h);
    error_data = wage_realized(:, i_h)-data_fc_RW(:, i_h);
    nonmiss    = ~isnan(error_TWIn) & ~isnan(error_data);
    n_RMSE     = nnz(nonmiss);
    error_TWIn = error_TWIn(nonmiss);
    error_data = error_data(nonmiss);
    T_boot     = nnz(nonmiss);
    T_boot     = block_size*floor(T_boot/block_size);    
    
    % Compute estimates RMSE
    RMSE_TWIn = sqrt( mean(error_TWIn.^2, 'omitnan'));
    RMSE_data = sqrt( mean(error_data.^2, 'omitnan'));
    RMSE_agg(i_h, :) = [RMSE_TWIn, RMSE_data];
    
    % Compute bootstrap standard errors
    RMSE_boot = NaN(n_boot, 3);
    for i_boot = 1:n_boot
        idx_tmp         = datasample(((1:(T_boot-block_size+1))')+(0:(block_size-1)), T_boot/block_size);
        error_TWIn_boot = error_TWIn(idx_tmp(:));
        error_data_boot = error_data(idx_tmp(:));
        RMSE_boot(i_boot, :) = [sqrt( mean(error_TWIn_boot.^2, 'omitnan')), sqrt( mean(error_data_boot.^2, 'omitnan')), ...
            sqrt( mean(error_TWIn_boot.^2, 'omitnan'))-sqrt( mean(error_data_boot.^2, 'omitnan'))];
    end
    stderr_agg(i_h, :) = std(RMSE_boot, [], 1);
    
    % Compute p-values
    pval_agg(i_h, :) = 1-mean(RMSE_boot(:, 1) <= RMSE_boot(:, 2));
    
end

% Store results in table
col_names     = {'TWInRMSE', 'TWInSE', 'AWGTRMSE', 'AWGTSE', 'Diff', 'DiffSE', 'pval'};
row_names     = cell(1, n_h); for i_h = 1:n_h, row_names{i_h} = [num2str(horizon(i_h)) '-month']; end
tab_aggregate = array2table(round([RMSE_agg(:, 1), stderr_agg(:, 1), RMSE_agg(:, 2), stderr_agg(:, 2), ...
    RMSE_agg(:, 2)-RMSE_agg(:, 1), stderr_agg(:, 3), pval_agg], 4), ...
    'variableNames', col_names, 'rowNames', row_names);
disp('ATLANTA FED WAGE GROWTH TRACKER')
disp(['No. of periods: ' num2str(n_RMSE)])
disp('Aggregate')
disp(tab_aggregate)
writetable(tab_aggregate, [tab_path 'Tab1_RMSEs.xlsx'], 'sheet', 'AWGT - Aggregate', 'writevariablenames', true, 'writerownames', true);
