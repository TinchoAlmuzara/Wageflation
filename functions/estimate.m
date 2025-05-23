function output = estimate(y, prior, settings, initial)
% ESTIMATE    Estimate Trend Wage Inflation Model using MCMC simulation
%
%   OUTPUT = ESTIMATE(Y, PRIOR, SETTINGS) estimates a Trend Wage Inflation
%   model given the observed data Y, a structure of prior hyperparameters,
%   and a structure of simulation settings. The function returns an OUTPUT
%   structure containing the MCMC draws for latent states and model parameters.
%
%   OUTPUT = ESTIMATE(Y, PRIOR, SETTINGS, INITIAL) allows the user to supply
%   an INITIAL structure with starting values for the latent variables and
%   parameters. If INITIAL is not provided, default initializations based on
%   the standard deviation of Y are used.
%
%   INPUTS:
%       y         - (T x n double) Matrix of observed data where T is the 
%                   number of time periods and n is the number of series.
%
%       prior     - (structure) Structure with prior hyperparameters:
%                     .prec_MA  - Precision for the MA coefficients.
%                     .nu_lam   - Degrees of freedom for the lambda prior.
%                     .s2_lam   - Scale parameter for the lambda prior.
%                     .nu_gam   - Degrees of freedom for the gamma prior.
%                     .s2_gam   - Scale parameter for the gamma prior.
%                     .a_ps     - Shape parameter a for the scale probability (ps)
%                     .b_ps     - Shape parameter b for the scale probability (ps)
%
%       settings  - (structure) Structure with simulation settings:
%                     .n_draw       - Number of MCMC draws (post burn-in).
%                     .n_burn       - Number of burn-in iterations.
%                     .n_thin       - Thinning factor for the MCMC draws.
%                     .show_progress- (Optional, logical) If true, the function
%                                    displays progress and plots during simulation.
%                     .n_lags       - (Optional) Number of lags for individual series.
%                                    Default is 0 if not provided.
%                     .n_lags_c     - (Optional) Number of lags for the common
%                                    component. Default is 0 if not provided.
%                     .is_timeag    - (Optional, logical vector of length n) Indicator
%                                    for series that are time-aggregated.
%
%       initial   - (structure, optional) Structure with initial values for the
%                   latent variables and parameters. If provided, the following
%                   fields are used to initialize the simulation:
%                     .alpha_tau, .sigma_dtau_c, .sigma_dtau_i, .alpha_eps,
%                     .sigma_eps_c, .s_eps_c, .sigma_eps_i, .s_eps_i, .theta_c,
%                     .theta, .lam_tau, .gam_dtau_c, .gam_dtau_i, .lam_eps,
%                     .gam_eps_c, .gam_eps_i, .ps_c, .ps_i
%
%   OUTPUT:
%       output - (structure) Contains the MCMC draws of the latent states and
%                model parameters. Key fields include:
%
%                   .tau_c         - (T x n_draw) Draws of the common trend.
%                   .alpha_tau     - (T x n x n_draw) Draws of scaling factors for
%                                    the common trend.
%                   .sigma_dtau_c  - (T x n_draw) Draws of the common trend volatility.
%                   .tau_i         - (T x n x n_draw) Draws of the individual trends.
%                   .sigma_dtau_i  - (T x n x n_draw) Draws of individual trend volatilities.
%                   .eps_c         - (T x n_draw) Draws of the common error component.
%                   .alpha_eps     - (T x n x n_draw) Draws of scaling factors for
%                                    the common error.
%                   .sigma_eps_c   - (T x n_draw) Draws of the common error volatility.
%                   .s_eps_c       - (T x n_draw) Draws of the scale parameters for the common error.
%                   .eps_i         - (T x n x n_draw) Draws of the individual error components.
%                   .sigma_eps_i   - (T x n x n_draw) Draws of the individual error volatilities.
%                   .s_eps_i       - (T x n x n_draw) Draws of the scale parameters for the individual errors.
%                   .theta_c       - (q_c x n_draw) Draws of the common moving-average (MA) parameters.
%                   .theta         - (n x max(q) x n_draw) Draws of the individual MA parameters.
%                   .lam_tau       - (n x n_draw) Draws of the lambda parameters for the trends.
%                   .gam_dtau_c    - (1 x n_draw) Draws of the gamma parameter for the common trend volatility.
%                   .gam_dtau_i    - (n x n_draw) Draws of the gamma parameters for the individual trend volatilities.
%                   .lam_eps       - (n x n_draw) Draws of the lambda parameters for the errors.
%                   .gam_eps_c     - (1 x n_draw) Draws of the gamma parameter for the common error volatility.
%                   .ps_c          - (1 x n_draw) Draws of the probability parameter for the common error scale.
%                   .gam_eps_i     - (n x n_draw) Draws of the gamma parameters for the individual error volatilities.
%                   .ps_i          - (n x n_draw) Draws of the probability parameter for the individual error scales.
%                   .y_draw        - (T x n x n_draw) Simulated draws of the observed data 
%                                    based on the estimated model.
%
%   ALGORITHM OVERVIEW:
%       1. **Initialization:** Recover simulation settings and dimensions from Y.
%          Set up latent variables and parameters either from INITIAL or via defaults.
%
%       2. **State-Space Setup:** Construct state-space models (SSM and SSM_TVC)
%          to represent the dynamics of trend components, errors, and their scales.
%
%       3. **Simulation Smoothing:** Use the simulation_smoother to initialize
%          the latent states (e.g., tau_c, tau_i, eps_c, eps_i) from the observed Y.
%
%       4. **MCMC Loop:** For each draw (including burn-in and thinning):
%            - Optionally display progress and update plots.
%            - Update auxiliary variables and state-space matrices.
%            - Update latent states using the simulation smoother.
%            - Update parameters (e.g., theta, lambda, gamma, volatilities, scales)
%              via functions such as update_theta, update_gam, update_vol,
%              update_scl, and update_ps.
%            - Compute simulated draws of Y using the updated state-space model.
%            - Save the draws if the iteration is post burn-in.
%
%       5. **Output:** Return the OUTPUT structure with all MCMC draws.
%
%   DEPENDENCIES:
%       This function requires auxiliary functions that perform the simulation
%       smoothing and parameter updates (e.g., simulation_smoother, update_theta,
%       update_gam, update_vol, update_scl, update_ps). Ensure these functions
%       are on the MATLAB path.
%
%   USAGE EXAMPLE:
%       % Define example data and settings:
%       T = 100; n = 3;
%       y = randn(T, n);
%       prior.prec_MA = 1;
%       prior.nu_lam = 3;
%       prior.s2_lam = 0.1;
%       prior.nu_gam = 3;
%       prior.s2_gam = 0.1;
%       prior.a_ps = 2;
%       prior.b_ps = 2;
%
%       settings.n_draw = 1000;
%       settings.n_burn = 500;
%       settings.n_thin = 10;
%       settings.n_lags = 2;
%       settings.n_lags_c = 1;
%       settings.show_progress = true;
%       settings.is_timeag = [false, true, false];
%
%       % Call the estimator:
%       output = estimate(y, prior, settings);
%
%   REFERENCES:
%       "A Measure of Trend Wage Inflation" (Almuzara, Audoly, Melcangi)
%
%   Author: Martin Almuzara

% Determine progress report
if isfield(settings, 'show_progress') 
    show_progress = settings.show_progress; 
    tic
    fig0          = figure();
    set(fig0, 'units', 'normalized', 'position', [0 0 1 1])
else
    show_progress = false; 
end

% Recover settings and dimensions
n_draw = settings.n_draw;
n_burn = settings.n_burn;
n_thin = settings.n_thin;
[T, n] = size(y);
if isfield(settings, 'n_lags_c'),  q_c = settings.n_lags_c; else, q_c = 0; end
if isfield(settings, 'n_lags'),    q   = settings.n_lags;   else, q   = 0; end
if isfield(settings, 'is_timeag'), is_timeag = settings.is_timeag; else, is_timeag = false(n, 1); end % which series are time-aggregated
t_skip = max(q+12*is_timeag)+1;

% Set theta restrictions and prior
if (length(q) == 1), q = repmat(q, [n, 1]); end
theta_r = ((1:max(q)) <= q);
prec_MA = prior.prec_MA;

% Set lambda/gamma priors
nu_lam = prior.nu_lam;
s2_lam = prior.s2_lam;
nu_gam = prior.nu_gam;
s2_gam = prior.s2_gam;

% Set support for s_eps and ps prior
a_ps     = prior.a_ps;
b_ps     = prior.b_ps;
n_s_vals = 10;
s_vals   = [1; linspace(2, 10, n_s_vals-1)'];

% Preallocate output
output              = struct();
output.tau_c        = NaN(T, n_draw);
output.alpha_tau    = NaN(T, n, n_draw);
output.sigma_dtau_c = NaN(T, n_draw);
output.tau_i        = NaN(T, n, n_draw);
output.sigma_dtau_i = NaN(T, n, n_draw);
output.eps_c        = NaN(T, n_draw);
output.alpha_eps    = NaN(T, n, n_draw);
output.sigma_eps_c  = NaN(T, n_draw);
output.s_eps_c      = NaN(T, n_draw);
output.eps_i        = NaN(T, n, n_draw);
output.sigma_eps_i  = NaN(T, n, n_draw);
output.s_eps_i      = NaN(T, n, n_draw);
output.theta_c      = NaN(q_c, n_draw);
output.theta        = NaN(n, max(q), n_draw);
output.lam_tau      = NaN(n, n_draw);
output.gam_dtau_c   = NaN(1, n_draw);
output.gam_dtau_i   = NaN(n, n_draw);
output.lam_eps      = NaN(n, n_draw);
output.gam_eps_c    = NaN(1, n_draw);
output.ps_c         = NaN(1, n_draw);
output.gam_eps_i    = NaN(n, n_draw);
output.ps_i         = NaN(n, n_draw);
output.y_draw       = NaN(T, n, n_draw);

% Initialize latent variables and parameters
if (nargin < 4)
    y_scale      = std(y, 'omitnan');
    alpha_tau    = (y_scale/16).*ones(T, n);
    sigma_dtau_c = ones(T, 1);
    sigma_dtau_i = (y_scale/16).*ones(T, n);
    alpha_eps    = (y_scale/16).*ones(T, n);
    sigma_eps_c  = ones(T, 1);
    s_eps_c      = ones(T, 1);
    sigma_eps_i  = (y_scale/16).*ones(T, n);
    s_eps_i      = ones(T, n);
    theta_c      = zeros(1, q_c);
    theta        = zeros(n, max(q));
    lam_tau      = sqrt(s2_lam)*ones(n, 1);
    gam_dtau_c   = sqrt(s2_gam);
    gam_dtau_i   = sqrt(s2_gam)*ones(n, 1);
    lam_eps      = sqrt(s2_lam)*ones(n, 1);
    gam_eps_c    = sqrt(s2_gam);
    gam_eps_i    = sqrt(s2_gam)*ones(n, 1);
    ps           = a_ps/(a_ps+b_ps);
    s_probs      = [ps, (1-ps)/(n_s_vals-1)*ones(1, n_s_vals-1)];
    s_probs_c    = s_probs;
    s_probs_i    = repmat(s_probs, n, 1);    
else
    alpha_tau     = initial.alpha_tau;
    sigma_dtau_c  = initial.sigma_dtau_c;
    sigma_dtau_i  = initial.sigma_dtau_i;
    alpha_eps     = initial.alpha_eps;
    sigma_eps_c   = initial.sigma_eps_c;
    s_eps_c       = initial.s_eps_c;
    sigma_eps_i   = initial.sigma_eps_i;
    s_eps_i       = initial.s_eps_i;
    theta_c       = initial.theta_c;
    theta         = initial.theta;
    lam_tau       = initial.lam_tau;
    gam_dtau_c    = initial.gam_dtau_c;
    gam_dtau_i    = initial.gam_dtau_i;
    lam_eps       = initial.lam_eps;
    gam_eps_c     = initial.gam_eps_c;
    gam_eps_i     = initial.gam_eps_i;
    s_probs_c     = [initial.ps_c, (1-initial.ps_c)/(n_s_vals-1)*ones(1, n_s_vals-1)];
    s_probs_i     = [initial.ps_i, (1-initial.ps_i)/(n_s_vals-1)*ones(1, n_s_vals-1)];
end

% Define indexing for state variables
if any(is_timeag), id_tau_c = 1:12; else, id_tau_c = 1; end
id_tmp   = length(id_tau_c);
id_tau_i = cell(n, 1);
for i = 1:n
    if is_timeag(i), id_tau_i{i} = id_tmp+(1:12); else, id_tau_i{i} = id_tmp+1; end
    id_tmp = id_tmp + length(id_tau_i{i});
end
id_eps_c = id_tmp + (1:(1+q_c));
id_tmp   = id_tmp + length(id_eps_c);
id_eps_i = cell(n, 1);
for i = 1:n
    id_eps_i{i} = id_tmp + (1:(1+q(i)));
    id_tmp      = id_tmp + length(id_eps_i{i});
end
n_state = id_tmp;

% Define auxiliary variables
eye_cell      = cell(n, 1); for i = 1:n, if is_timeag(i), eye_cell{i} = repmat(1/12, [1, 12]); else, eye_cell{i} = 1; end; end
eye_blk       = blkdiag(eye_cell{:});
theta_cell    = cell(n, 1); for i = 1:n, theta_cell{i} = [1, theta(i, theta_r(i, :))]; end
theta_blk     = blkdiag(theta_cell{:});
theta_blk_c   = [1, theta_c];
sigmaXs_eps_c = sigma_eps_c.*s_eps_c;
sigmaXs_eps_i = sigma_eps_i.*s_eps_i;

% Initialize state-space model
SSM = struct();
%%% H
SSM.H = zeros(n, n_state, T+1);
for t = 1:T
    SSM.H(~is_timeag, id_tau_c(1), t+1) = alpha_tau(t, ~is_timeag)';
    SSM.H(is_timeag, id_tau_c, t+1)     = (1/12)*[alpha_tau(t:(-1):max(1, t-11), is_timeag)', repmat(alpha_tau(1, is_timeag)', [1, max(12-t, 0)])];
    SSM.H(:, [id_tau_i{:}], t+1)        = eye_blk;
    SSM.H(:, id_eps_c, t+1)             = (alpha_eps(t, :)')*theta_blk_c;
    SSM.H(:, [id_eps_i{:}], t+1)        = theta_blk;
end
%%% Sigma_eps
SSM.Sigma_eps = (1e-6)*eye(n);
%%% F
SSM.F                     = zeros(n_state);
n_tmp                     = length(id_tau_c);
SSM.F(id_tau_c, id_tau_c) = [1, zeros(1, n_tmp-1); eye(n_tmp-1), zeros(n_tmp-1, 1)];
for i = 1:n
    n_tmp                           = length(id_tau_i{i});
    SSM.F(id_tau_i{i}, id_tau_i{i}) = [1, zeros(1, n_tmp-1); eye(n_tmp-1), zeros(n_tmp-1, 1)];
end
n_tmp                     = length(id_eps_c);
SSM.F(id_eps_c, id_eps_c) = [zeros(1, n_tmp); eye(n_tmp-1), zeros(n_tmp-1, 1)];
for i = 1:n
    n_tmp                           = length(id_eps_i{i});
    SSM.F(id_eps_i{i}, id_eps_i{i}) = [zeros(1, n_tmp); eye(n_tmp-1), zeros(n_tmp-1, 1)];
end
%%% G
SSM.G                = zeros(n_state, 2*(n+1));
n_tmp                = length(id_tau_c);
SSM.G(id_tau_c, 1)   = [1; zeros(n_tmp-1, 1)];
for i = 1:n
    n_tmp                   = length(id_tau_i{i});
    SSM.G(id_tau_i{i}, 1+i) = [1; zeros(n_tmp-1, 1)];
end
n_tmp                = length(id_eps_c);
SSM.G(id_eps_c, 2+n) = [1; zeros(n_tmp-1, 1)];
for i = 1:n
    n_tmp                     = length(id_eps_i{i});
    SSM.G(id_eps_i{i}, 2+n+i) = [1; zeros(n_tmp-1, 1)];
end
%%% Sigma_eta
SSM.Sigma_eta = zeros(2*(n+1), 2*(n+1), T);
for t = 1:T
    SSM.Sigma_eta(:, :, t) = diag([sigma_dtau_c(t), sigma_dtau_i(t, :), ...
                                  sigmaXs_eps_c(t), sigmaXs_eps_i(t, :)].^2);    
end
%%% mu_1 and Sigma_1
SSM.mu_1                        = zeros(n_state, 1);
SSM.Sigma_1                     = zeros(n_state);
for i = 1:n
    SSM.Sigma_1(id_tau_i{i}, id_tau_i{i}) = eye(length(id_tau_i{i})) + 1e1*ones(length(id_tau_i{i}));
    SSM.Sigma_1(id_eps_i{i}, id_eps_i{i}) = eye(length(id_eps_i{i}));
end
SSM.Sigma_1(id_tau_c, id_tau_c) = 0;
SSM.Sigma_1(id_eps_c, id_eps_c) = 0;

% Define indexing for TVCs
id_TVC = cell(n, 1);
id_tmp = 0;
for i = 1:n
    if is_timeag(i), id_TVC{i} = id_tmp + (1:12); else, id_TVC{i} = id_tmp + 1; end
    id_tmp = id_tmp + length(id_TVC{i});
end
n_TVC  = id_tmp;

% Initialize state-space model for TVCs
SSM_TVC = struct();
%%% H
SSM_TVC.H = zeros(n, n_TVC+n+size(theta_blk, 2), T);
for t = 1:T
    SSM_TVC.H(:, :, t) = [zeros(n, n_TVC+n), theta_blk];
end
%%% Sigma_eps
SSM_TVC.Sigma_eps = (1e-6)*eye(n);
%%% F
SSM_TVC.F = blkdiag(zeros(n_TVC+n), SSM.F([id_eps_i{:}], [id_eps_i{:}]));
for i = 1:n
    n_tmp                           = length(id_TVC{i});
    SSM_TVC.F(id_TVC{i}, id_TVC{i}) = [1, zeros(1, n_tmp-1); eye(n_tmp-1), zeros(n_tmp-1, 1)];
end
SSM_TVC.F(n_TVC+(1:n), n_TVC+(1:n)) = eye(n);
%%% G
SSM_TVC.G = blkdiag(zeros(n_TVC+n, 2*n), SSM.G([id_eps_i{:}], 2+n+(1:n)));
for i = 1:n
    n_tmp                   = length(id_TVC{i});
    SSM_TVC.G(id_TVC{i}, i) = [1; zeros(n_tmp-1, 1)];
end
SSM_TVC.G(n_TVC+(1:n), n+(1:n)) = eye(n);
%%% Sigma_eta
SSM_TVC.Sigma_eta = zeros(3*n, 3*n, T-1);
%%% mu_1 and Sigma_1
SSM_TVC.mu_1    = zeros(n_TVC+n+size(theta_blk, 2), 1);
Sigma_aux       = zeros(n_TVC); for i = 1:n, Sigma_aux(id_TVC{i}, id_TVC{i}) = 1; end
SSM_TVC.Sigma_1 = eye(n_TVC+n+size(theta_blk, 2)) ...
     + 1e1*blkdiag(Sigma_aux, ones(n), zeros(size(theta_blk, 2)));

% Initialize state variables
state_smooth = simulation_smoother([NaN(1, n); y]', SSM);
tau_c        = state_smooth(id_tau_c(1), 2:(T+1))';
tau_i        = NaN(T, n);
for i = 1:n
    tau_i(:, i) = state_smooth(id_tau_i{i}(1), 2:(T+1))';
end
eps_c        = state_smooth(id_eps_c(1), 2:(T+1))';
ups_c        = (theta_blk_c*state_smooth(id_eps_c, 2:(T+1)))';
eps_i        = NaN(T, n);
for i = 1:n
    eps_i(:, i) = state_smooth(id_eps_i{i}(1), 2:(T+1))';
end

for i_draw = (-n_burn):n_draw
    
    if show_progress
        % Show progress
        message = sprintf('MCT model - Draw %d/%d\n', i_draw, n_draw);
        fprintf(message)

        % Plot progress
        if (mod(i_draw, 10) == 0)
            y_total  = mean(y, 2, 'omitnan');
            t_common = mean(alpha_tau, 2) .* tau_c;
            t_total  = mean(alpha_tau, 2) .* tau_c + mean(tau_i, 2);
            data_all = [y_total, t_common, t_total];
            ax0 = subplot(2, 3, 1); plot0 = plot(ax0, data_all);
            xlim(ax0, [1, T]), xlabel('t')
            ylim(ax0, [min(data_all(:))-0.5, max(data_all(:))+0.5])
            ylabel(ax0, 'Common (green) and total trend (violet)', 'Interpreter', 'latex')
            set(plot0, {'linewidth'}, {1; 1.5; 1.5})
            set(plot0, {'color'}, {[0, 0, 0]; [0, 2/3, 0]; [2/3, 0, 1/3]})
            ax0 = subplot(2, 3, 2); plot0 = plot(ax0, [sigma_dtau_i, sigma_eps_i]);
            xlim(ax0, [1, T]), xlabel('t')
            ylabel(ax0, '$\sigma_{\tau, i}$ (red) and $\sigma_{\varepsilon, i}$ (blue)', 'Interpreter', 'latex')
            set(plot0, 'linewidth', 1)
            set(plot0, {'color'}, [repmat({[2/3, 0, 0]}, [n, 1]); repmat({[0, 0, 2/3]}, [n, 1])])
            ax0 = subplot(2, 3, 3); plot0 = plot(ax0, [alpha_tau, alpha_eps]);
            xlim(ax0, [1, T]), xlabel('t')
            ylabel(ax0, '$\alpha_{\tau, i}$ (red) and $\alpha_{\varepsilon, i}$ (blue)', 'Interpreter', 'latex')
            set(plot0, 'linewidth', 1)
            set(plot0, {'color'}, [repmat({[2/3, 0, 0]}, [n, 1]); repmat({[0, 0, 2/3]}, [n, 1])])
            ax0 = subplot(2, 3, 4); plot0 = plot([sigma_dtau_c, sigma_eps_c]);
            xlim(ax0, [1, T]), xlabel('t')
            ylabel(ax0, '$\sigma_{\tau, c}$ (red) and $\sigma_{\varepsilon, c}$ (blue)', 'Interpreter', 'latex')
            set(plot0, 'linewidth', 1.5)
            set(plot0, {'color'}, {[2/3, 0, 0]; [0, 0, 2/3]})
            ax0 = subplot(2, 3, 5); plot0 = plot(ax0, tau_c);
            xlim(ax0, [1, T]), xlabel('t')
            ylabel(ax0, '$\tau_{c}$', 'Interpreter', 'latex')
            set(plot0, 'linewidth', 1.5)
            set(plot0, 'color', [2/3, 0, 0])
            ax0 = subplot(2, 3, 6); plot0 = plot(ax0, eps_c);
            xlim(ax0, [1, T]), xlabel('t')
            ylabel(ax0, '$\varepsilon_{c}$', 'Interpreter', 'latex')
            set(plot0, 'linewidth', 1.5)
            set(plot0, 'color', [0, 0, 2/3])            
            drawnow
        end
    end

    for i_thin = 1:n_thin

        % Define auxiliary variables
        theta_cell    = cell(n, 1); for i = 1:n, theta_cell{i} = [1, theta(i, theta_r(i, :))]; end
        theta_blk     = blkdiag(theta_cell{:});
        theta_blk_c   = [1, theta_c(:)'];
        sigmaXs_eps_c = sigma_eps_c.*s_eps_c;
        sigmaXs_eps_i = sigma_eps_i.*s_eps_i;
        
        % Update alpha
        y_TVC = (y') - eye_blk*state_smooth([id_tau_i{:}], 2:(T+1));
        for t = 1:T
            for i = 1:n
                if is_timeag(i)
                    SSM_TVC.H(i, id_TVC{i}, t) = (1/12)*[tau_c(t:(-1):max(1, t-11))', repmat(tau_c(1), [1, max(12-t, 0)])];
                else
                    SSM_TVC.H(i, id_TVC{i}, t) = tau_c(t);
                end
            end
            SSM_TVC.H(:, n_TVC+(1:n), t) = ups_c(t)*eye(n);
            SSM_TVC.H(:, n_TVC+n+(1:size(theta_blk, 2)), t) = theta_blk;
        end
        for t = 2:T
            SSM_TVC.Sigma_eta(:, :, t-1) = diag([lam_tau', lam_eps', sigmaXs_eps_i(t, :)].^2); 
        end
        alpha_smooth = simulation_smoother(y_TVC, SSM_TVC);
        for i = 1:n
            alpha_tau(:, i) = alpha_smooth(id_TVC{i}(1), :)';
            alpha_eps(:, i) = alpha_smooth(n_TVC+i, :)';
        end
        
        % Draw lambda
        lam_tau = update_gam(diff(alpha_tau((t_skip+1):T, :), 1, 1), nu_lam, s2_lam);
        lam_eps = update_gam(diff(alpha_eps((t_skip+1):T, :), 1, 1), nu_lam, s2_lam);
                
        % Update state-space representation
        for t = 1:T
            SSM.H(~is_timeag, id_tau_c(1), t+1) = alpha_tau(t, ~is_timeag)';
            SSM.H(is_timeag, id_tau_c, t+1)     = (1/12)*[alpha_tau(t:(-1):max(1, t-11), is_timeag)', repmat(alpha_tau(1, is_timeag)', [1, max(12-t, 0)])];
            SSM.H(:, id_eps_c, t+1)             = (alpha_eps(t, :)')*theta_blk_c;
            SSM.H(:, [id_eps_i{:}], t+1)        = theta_blk;
            SSM.Sigma_eta(:, :, t)              = diag([sigma_dtau_c(t), sigma_dtau_i(t, :), ...
                                                       sigmaXs_eps_c(t), sigmaXs_eps_i(t, :)].^2);
        end
        SSM.H(:, :, 1) = SSM.H(:, :, 2);

        % Draw states
        state_smooth = simulation_smoother([NaN(1, n); y]', SSM);
        tau_c        = state_smooth(id_tau_c(1), 2:(T+1))';
        dtau_c       = (state_smooth(id_tau_c(1), 2:(T+1)) - state_smooth(id_tau_c(1), 1:T))';
        tau_i        = NaN(T, n);
        dtau_i       = NaN(T, n);
        for i = 1:n
            tau_i(:, i)  = state_smooth(id_tau_i{i}(1), 2:(T+1))';
            dtau_i(:, i) = (state_smooth(id_tau_i{i}(1), 2:(T+1)) - state_smooth(id_tau_i{i}(1), 1:T))';
        end
        eps_c        = state_smooth(id_eps_c(1), 2:(T+1))';
        ups_c        = (theta_blk_c*state_smooth(id_eps_c, 2:(T+1)))';
        eps_i        = NaN(T, n);
        for i = 1:n
            eps_i(:, i) = state_smooth(id_eps_i{i}(1), 2:(T+1))';            
        end
        ups_i        = (theta_blk*state_smooth([id_eps_i{:}], 2:(T+1)))';

        % Draw theta
        if (q_c > 0)
            y_theta = ups_c((q_c+1):T) ./ sigmaXs_eps_c((q_c+1):T);
            x_theta = NaN(T-q_c, q_c);
            for lag = 1:q_c, x_theta(:, lag) = eps_c((q(i)-lag+1):(T-lag)) ./ sigmaXs_eps_c((q_c+1):T); end
            theta_c(1:q_c) = update_theta(y_theta, x_theta, prec_MA);
        end
        
        % Draw theta_c
        for i = 1:n
            if (q(i) > 0)
                y_theta = ups_i((q(i)+1):T, i) ./ sigmaXs_eps_i((q(i)+1):T, i);
                x_theta = NaN(T-q(i), q(i));
                for lag = 1:q(i), x_theta(:, lag) = eps_i((q(i)-lag+1):(T-lag), i) ./ sigmaXs_eps_i((q(i)+1):T, i); end
                theta(i, 1:q(i)) = update_theta(y_theta, x_theta, prec_MA);
            end
        end
        
        % Draw volatilities
        sigma_dtau_c = update_vol(dtau_c, sigma_dtau_c, gam_dtau_c, 0);
        for i = 1:n
            sigma_dtau_i(:, i) = update_vol(dtau_i(:, i), sigma_dtau_i(:, i), gam_dtau_i(i));
        end
        sigma_eps_c  = update_vol(eps_c./s_eps_c, sigma_eps_c, gam_eps_c, 0);
        for i = 1:n
            sigma_eps_i(:, i)  = update_vol(eps_i(:, i)./s_eps_i(:, i), sigma_eps_i(:, i), gam_eps_i(i));
        end
                   
        % Draw gamma
        gam_dtau_c = update_gam(2*diff(log(sigma_dtau_c((t_skip+1):T)), 1, 1), nu_gam/12, s2_gam);
        gam_dtau_i = update_gam(2*diff(log(sigma_dtau_i((t_skip+1):T, :)), 1, 1), nu_gam, s2_gam);
        gam_eps_c  = update_gam(2*diff(log(sigma_eps_c((t_skip+1):T)), 1, 1), nu_gam, s2_gam);
        gam_eps_i  = update_gam(2*diff(log(sigma_eps_i((t_skip+1):T, :)), 1, 1), nu_gam, s2_gam);

        % Draw scales (outliers)
        s_eps_c = update_scl([NaN(t_skip, 1); eps_c((t_skip+1):T)./sigma_eps_c((t_skip+1):T)], s_vals, s_probs_c);
        for i = 1:n
            s_eps_i(:, i) = update_scl([NaN(t_skip, 1); eps_i((t_skip+1):T, i)./sigma_eps_i((t_skip+1):T, i)], s_vals, s_probs_i(i, :));                
        end
        
        % Draw ps 
        ps_c      = update_ps(s_eps_c((t_skip+1):T), a_ps, b_ps);
        ps_i      = update_ps(s_eps_i((t_skip+1):T, :), a_ps, b_ps);
        s_probs_c = [ps_c, (1-ps_c)/(n_s_vals-1)*ones(1, n_s_vals-1)];
        s_probs_i = [ps_i, (1-ps_i)/(n_s_vals-1)*ones(1, n_s_vals-1)];
                
    end
    
    % Compute data draws
    y_draw = NaN(n, T);
    for t = 1:T
        y_draw(:, t) = SSM.H(:, :, t+1)*state_smooth(:, t+1);
    end

    % Save draws
    if (i_draw > 0)
        output.tau_c(:, i_draw)           = tau_c;
        output.alpha_tau(:, :, i_draw)    = alpha_tau;
        output.sigma_dtau_c(:, i_draw)    = sigma_dtau_c;
        output.tau_i(:, :, i_draw)        = tau_i;
        output.sigma_dtau_i(:, :, i_draw) = sigma_dtau_i;
        output.eps_c(:, i_draw)           = eps_c;
        output.alpha_eps(:, :, i_draw)    = alpha_eps;
        output.sigma_eps_c(:, i_draw)     = sigma_eps_c;
        output.s_eps_c(:, i_draw)         = s_eps_c;
        output.eps_i(:, :, i_draw)        = eps_i;
        output.sigma_eps_i(:, :, i_draw)  = sigma_eps_i;
        output.s_eps_i(:, :, i_draw)      = s_eps_i;
        output.theta_c(:, i_draw)         = theta_c;
        output.theta(:, :, i_draw)        = theta;
        output.lam_tau(:, i_draw)         = lam_tau;
        output.gam_dtau_c(i_draw)         = gam_dtau_c;
        output.gam_dtau_i(:, i_draw)      = gam_dtau_i;
        output.lam_eps(:, i_draw)         = lam_eps;
        output.gam_eps_c(i_draw)          = gam_eps_c;
        output.ps_c(i_draw)               = ps_c;
        output.gam_eps_i(:, i_draw)       = gam_eps_i;
        output.ps_i(:, i_draw)            = ps_i;        
        output.y_draw(:, :, i_draw)       = y_draw';
    end
    
end

end
