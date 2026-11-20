function Omega_h2 = freeze_loop(mx, mphi, gx, g, interaction, T_RH, doplot, production)
% freeze_loop  Self-consistent relic abundance from the full Boltzmann ODE.
%
% Solves
%   dY/dx = (1/x^2) * corr * (g_s/sqrt(g_rho)) * (lambda_prod*Y_eq^2 - lambda_full*Y^2)
% with Y = 0 at x_init = mx/T_RH (no DM above reheating).
%
% Production / depletion are split:
%   lambda_prod uses <sigma v> into SM final states (gamma gamma + meson
%     pairs; svavg channel 'sm').  The SM bath is in equilibrium, so by
%     detailed balance this is production from gamma gamma, pi pi, K K, ...
%     The mesons are Boltzmann-suppressed, but so are the photon pairs with
%     sqrt(s) > 2 m_pi, and above that the tree-level meson channels dominate.
%     Dark final states (phi phi, A'A') are not produced from the bath.
%   lambda_full uses the full <sigma v> over every open annihilation
%     channel.  Depletion has no thermal-bath assumption: the cross-section
%     is what it is, and chi pairs annihilate to whatever final state
%     regardless of whether those products are populated in the bath.
%
% When all bath species are in equilibrium (high T), the two coincide and
% the form degenerates to the standard (Y_eq^2 - Y^2).
%
% Caveat:  leaving phi out of the bath breaks down when mphi <~ T_RH (real
% on-shell phi becomes a bath species in its own right).  A runtime warning
% fires in that regime.
%
% Limits:
%   * gx small (FIMP):   depletion negligible, recovers freeze-in.
%   * gx large (WIMP):   Y rapidly tracks Y_eq, then freezes out -> matches
%                        freeze_out.m at the same gx.
%   * crossing region:   gives the correct self-consistent Omega h^2 that
%                        neither linearized FI nor equilibrium-IC FO captures.
%
% Notes:
%   * Uses svavg.m / xxggxsec.m for <sigma v>.  For interactions where
%     xxggxsec.m is not implemented (lightScalar*, *PP) this will return 0.
%   * Y_eq uses the K2-based form (valid for any x), not the NR formula
%     in freeze_out.m, so high-T (x << 1) production is captured correctly.
%
% Usage:
%   Oh2 = freeze_loop(0.1, 10, 1, [1 1 1]/3, 'heavyScalar')
%   Oh2 = freeze_loop(0.1, 10, 1, [1 1 1]/3, 'heavyScalarPS', true)

if nargin < 6 || isempty(T_RH), T_RH = 0.01; end   % default reheat T [GeV]
if nargin < 7, doplot = false; end
% production: 'sm' (default; gamma gamma + meson pairs) or 'gg' (photons only)
if nargin < 8 || isempty(production), production = 'sm'; end

g_chi = 4;   % Dirac fermion: 2 spins x 2 (particle+antiparticle)

% -----------------------------------------------------------------------
% x range:  start integration at the reheat temperature (no DM, no
% photon-bath production above T_RH) and run through freeze-out.
% -----------------------------------------------------------------------
x_init  = mx / T_RH;
x_final = max(1000, 10*x_init);    % ensure we always integrate well past freeze-out
T_min   = mx / x_final;
T_max   = T_RH;

% Warn if leaving phi out of the bath is likely poor.
if mphi <= 50*T_RH
    warning('freeze_loop:phiBathRegime', ...
        ['mphi = %.3g GeV is within 50x T_RH = %.3g GeV; real phi may be ', ...
         'a non-negligible bath species and the SM-bath-only production ', ...
         'approximation is suspect.  Consider a coupled (Y_chi, Y_phi) ODE.'], ...
         mphi, T_RH);
end

% -----------------------------------------------------------------------
% Pre-compute <sigma v>(T) on the integration grid -- both the full
% cross-section (depletion) and the SM-final-state piece (production).
% -----------------------------------------------------------------------
Tgrid    = logspace(log10(T_min), log10(T_max), 200);
logTgrid = log(Tgrid);
svvals_full = svavg_grid(mx, mphi, gx, g, Tgrid, interaction);
svvals_sm   = svavg_grid(mx, mphi, gx, g, Tgrid, interaction, production);
% Interpolated in log <sv> (as the Julia port): the resonance and meson
% thresholds make <sv> fall exponentially, which pchip in value misses by several %.
svFull = griddedInterpolant(logTgrid(:), log(max(svvals_full(:), realmin)), 'pchip');
svProd = griddedInterpolant(logTgrid(:), log(max(svvals_sm(:),   realmin)), 'pchip');

% Precompute g_s, g_rho, d(ln g_s)/d(ln T) interpolants over the integration
% T range.  Avoids 4 gstarofT_vec calls per ode15s RHS evaluation.
S = gstar_interp(T_min, T_max);

% -----------------------------------------------------------------------
% Initial condition: Y = 0 at T_RH, first step taken analytically.
% Over a step dtau short against the variation of Y_eq (~1/x per e-fold)
% the coefficients are constant and dY/dtau = a - b Y^2 solves exactly:
%   Y(dtau) = Y* tanh(sqrt(a b) dtau),  Y* = sqrt(a/b),  a = cp Y_eq^2.
% This is ~a dtau at weak coupling and exactly Y* when the rate is huge, so
% the start, and hence Omega h^2, is continuous in the coupling.  (A rate
% test choosing Y = 0 or Y = Y_eq made Omega h^2(eta) jump by up to ~500x
% near mx ~ 20 T_RH, which fzero then reported as a root.)  The seed also
% spares ode15s the fast initial transient at large coupling.
% -----------------------------------------------------------------------
tau0 = log(x_init) + 1e-3 / max(1, x_init);
[cp0, cf0, Yeq0] = boltz_coef(exp(tau0));
a0 = cp0 * Yeq0^2;
dtau0 = tau0 - log(x_init);
if cf0 > 0
    Y0 = sqrt(a0/cf0) * tanh(sqrt(a0*cf0) * dtau0);
else
    Y0 = a0 * dtau0;
end

% -----------------------------------------------------------------------
% Integrate in tau = log(x) to reduce stiffness at small x.  Provide an
% explicit Jacobian; loosen tolerances slightly for robustness across
% the full FIMP-WIMP range.
% -----------------------------------------------------------------------
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-18, ...
              'Jacobian', @jac_tau);
warnState = warning('off', 'MATLAB:ode15s:IntegrationTolNotMet');
lastwarn('');
[tau_sol, Y_sol] = ode15s(@rhs_tau, [tau0, log(x_final)], Y0, opts);
[~, warnId] = lastwarn;
warning(warnState);

x_sol = exp(tau_sol);

% If ode15s bailed before reaching x_final, the returned Y_sol(end) is the
% state at the failure point — meaningless for relic abundance.  Mark NaN
% so fzero / the scan sees a clear signal rather than noise.
reached_end = tau_sol(end) >= log(x_final) - 1e-6;
if ~reached_end || strcmp(warnId, 'MATLAB:ode15s:IntegrationTolNotMet')
    Omega_h2 = NaN;
    return;
end

Y_now = Y_sol(end);

s0          = 2891.2;    % cm^-3   (CMB entropy density today)
rho_crit_h2 = 1.054e-5;  % GeV cm^-3
Omega_h2    = mx * s0 * Y_now / rho_crit_h2;

% -----------------------------------------------------------------------
% Optional plot
% -----------------------------------------------------------------------
if doplot
    x_plot    = logspace(log10(x_init), log10(x_final), 600);
    Y_eq_plot = arrayfun(@(xx) Y_eq_of(xx), x_plot);

    figure('Name', 'Self-consistent loop integration');
    loglog(x_sol, max(Y_sol, 1e-40), 'b-', 'LineWidth', 2); hold on;
    loglog(x_plot, Y_eq_plot, 'r--', 'LineWidth', 1.5);
    xlabel('x = m_\chi / T'); ylabel('Y = n / s');
    legend('Y (loop)', 'Y_{eq}', 'Location', 'best');
    title(sprintf('Loop:  m_\\chi = %.4g GeV,  g_\\chi = %.4g,  \\Omega h^2 = %.3g', ...
                  mx, gx, Omega_h2));
    xlim([x_init x_final]); grid on;
end

% -----------------------------------------------------------------------
% Y_eq(x) using the K2-based formula (valid at any x; reduces to
% 0.145 g_chi/g_s * x^{3/2} e^{-x} for x >> 1).
% -----------------------------------------------------------------------
    function Yeq = Y_eq_of(x)
        T_   = mx / x;
        g_s_ = S.gs(log(T_));   % same g_s source as boltz_coef (consistency)
        if x > 30
            % asymptotic expansion (error < 1e-5 already at x=30)
            K2_ = exp(-x + 0.5*log(pi/(2*x))) * (1 + 15/(8*x));
        else
            K2_ = besselk(2, x);
        end
        Yeq = (45/(4*pi^4)) * (g_chi/g_s_) * x^2 * K2_;
    end

% -----------------------------------------------------------------------
% Boltzmann RHS in tau = log(x):  dY/dtau = x * dY/dx.
% This removes one factor of 1/x from the prefactor, easing stiffness
% at small x.
% -----------------------------------------------------------------------
    function [coef_prod, coef_full, Yeq] = boltz_coef(x)
        T     = mx / x;
        lT    = log(T);
        g_rho = S.grho(lT);
        g_s   = S.gs(lT);

        correction = 1 + S.dlngs(lT)/3;

        sv_p = exp(svProd(log(T)));
        sv_f = exp(svFull(log(T)));
        % Y counts chi + chibar (g_chi = 4) and <sigma v> is the chi chibar cross
        % section, so n-dot = -(<sigma v>/2)(n^2 - n_eq^2): the 1/2 is the Dirac factor.
        lambda_p = 0.5 * sqrt(pi/45) * 1.22e19 * mx * sv_p;
        lambda_f = 0.5 * sqrt(pi/45) * 1.22e19 * mx * sv_f;

        coef_prod = (lambda_p / x) * correction * (g_s / sqrt(g_rho));
        coef_full = (lambda_f / x) * correction * (g_s / sqrt(g_rho));
        Yeq  = Y_eq_of(x);
    end

    function dYdtau = rhs_tau(tau, Y)
        x = exp(tau);
        [cp, cf, Yeq] = boltz_coef(x);
        dYdtau = cp * Yeq^2 - cf * Y^2;
    end

    function J = jac_tau(tau, Y)
        x = exp(tau);
        [~, cf, ~] = boltz_coef(x);
        J = -2 * cf * Y;
    end

end
