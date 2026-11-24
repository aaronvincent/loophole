function [Omega_h2, info] = freeze_coupled(mx, mphi, gx, g, interaction, T_RH, doplot, production)
% freeze_coupled  Coupled (Y_chi, Y_phi) Boltzmann ODE.
%
% Generalization of freeze_loop.m that explicitly tracks the mediator
% abundance Y_phi alongside Y_chi.  Intended for the regime where mphi
% is comparable to or below T_RH, in which phi becomes a non-negligible
% bath species and freeze_loop.m, which has no Y_phi, breaks down.
%
% Y_chi starts at 0 at T_RH with an analytic first step (see the IC block).
% Two choices of variable for Y_chi, picked by that seed:
%   * seed within a factor 2 of its quasi-steady value: state =
%       (W_chi, Y_phi) with W_chi = log(Y_chi).  Log transform on Y_chi
%       follows freeze_out.m and is needed to integrate cleanly through
%       the deep-WIMP equilibrium-tracking regime; linear-Y is too stiff
%       for ode15s with reasonable tolerances when Y tracks Y_eq over
%       many decades.
%   * otherwise: state = (Y_chi, Y_phi) linear.
%
% Model.  rho = Y_chi/Y_chi_eq.  Per e-fold in x (tau = log x):
%
%   dY_chi/dtau = cp_s Y_chi_eq^2 - cf_s Y_chi^2                 continuum
%               + 2 c_xx (Y_phi - Y_phi_eq rho^2)                 phi <-> chi chibar
%               + c_rev Y_phi^2                                   phi phi -> chi chibar
%   dY_phi/dtau = -c_SM (Y_phi - Y_phi_eq)                        phi <-> SM bath
%               - c_xx (Y_phi - Y_phi_eq rho^2)                   phi <-> chi chibar
%               - c_mes Y_phi                                     phi -> mesons ('gg')
%               + c_d Y_chi^2 - c_rev Y_phi^2                     chi chibar <-> phi phi
%               - c_AA (Y_phi^2 - Y_phi_eq^2)                     A'A' <-> gamma gamma, MM
%                                                                 (axialVector only)
%
% (chi chibar -> phi phi depletion sits in cf_s, which is built from 'all'.)
%   c_xx  = Gamma_xx <K1/K2> / H,
%   c_SM  = [(Gamma_gg + Gamma_meson [if 'sm']) <K1/K2> + Gamma_sc] / H,
%   c_mes = Gamma_meson <K1/K2> / H with 'gg' (mesons not in the bath), else 0,
%   c_d, c_rev = lambda <sv>_fwd,rev / x for chi chibar <-> phi phi (A'A'; SP/PP
%         phi'phi'), tabulated directly in log by sv_darkpair_grid, each normalised
%         by its own initial state so the exothermic direction stays O(1) at any
%         m/T; c_rev Y_phi_eq^2 = c_d Y_chi_eq^2 exactly.
% All carry the g_s correction (1 + dln g_s/3 dln T); every lambda carries the
% Dirac 1/2.
%
% cp_s, cf_s are the 'sm' (or 'gg') and 'all' coefficients with the on-shell phi
% removed, since phi is tracked explicitly.  By detailed balance that part is
%   2 c_xx (Y_phi_eq / Y_chi_eq^2) BR,    BR = Gamma(phi -> final states)/Gamma_phi,
% (production: the production final states; depletion: every SM one).  Where the
% resonance dominates, the subtraction leaves a small difference of large numbers
% and the NWA matches the solver's <sv> to ~1%, so the remainder is clamped at 0.
%
% Inverse decay chi chibar -> phi carries rho^2 (it needs chi), so phi is not
% held at equilibrium by a chi population that is not there: in freeze-in with
% fast phi -> chi chibar, Y_phi sits at Y_phi_eq c_SM/(c_SM + c_xx), far below
% equilibrium, and on-shell production is 2 c_xx Y_phi = the old resonant rate.
% Before 2026-11-16 phi <-> chi chibar used rho = 1 and phi phi -> chi chibar was
% absent.
%
% Fast phi (adiabatic_phi: total phi rate / H > 1e3 over the range): Y_phi is
% eliminated at its exact quasi-steady value (a quadratic, chi_exch_adia), and
% the system is the scalar chi equation.  This avoids the O(1e6) steps ode15s
% takes chasing Y_phi's tiny deviations when the rates reach 1e19 H.
% LOOP_FC_FULL2D=1 in the environment disables the shortcut (tests).
%
% Y_phi(T_RH) = Y_phi_eq (1 - exp(-Gamma_SM dtau0/H)): the same first step as chi.
% <sv> tables are interpolated in log.
%
% Gamma_sc(T) is phi's 2 -> 2 conversion rate on the bath (in c_SM only: these
% processes make no chi):
%   * e+- gamma <-> e+- phi and e+ e- <-> gamma phi through the phi gamma gamma
%     vertex (primakoff_phi_rate), every non-axial model;
%   * pi+- gamma <-> pi+- phi (pigamma_phi_rate), S⊗S / P⊗S with mphi < 2 m_pi
%     and production = 'sm' (pions in the bath);
%   * A' pi <-> pi pi (apipi_phi_rate), axialVector with production = 'sm':
%     the contact term of the axial 3pi current; LOOP_APIPI=0 removes it.
% All are tabulated once on a fixed log-T grid (5 keV .. T_RH, independent of
% mx) and interpolated.  See docs/phi_thermalisation_check.md.
%
% Supported interactions: heavyScalar and heavyScalarPS (anything that
% phi_partial_widths.m handles).
%
% Usage:
%   Oh2 = freeze_coupled(0.1, 0.05, 1, [2 -1 -1]/3, 'heavyScalar')
%   Oh2 = freeze_coupled(..., 0.01, true)   % with diagnostic plot

if nargin < 6 || isempty(T_RH), T_RH = 0.01; end
if nargin < 7, doplot = false; end
% production: 'sm' (default; gamma gamma + meson pairs) or 'gg' (photons only)
if nargin < 8 || isempty(production), production = 'sm'; end

g_chi = 4;   % Dirac fermion
if strcmp(interaction, 'axialVector'), g_phi = 3; else, g_phi = 1; end   % A' is massive spin 1

Mpl = 1.22e19;

% Warn if outside the intended regime.  freeze_loop is faster / more
% accurate when mphi >> T_RH and there is no real-phi bath population.
if mphi > 50 * T_RH
    warning('freeze_coupled:heavyMediator', ...
        ['mphi = %.3g GeV >> T_RH = %.3g GeV.  freeze_loop is the ', ...
         'better tool here; freeze_coupled adds little and is slower.'], ...
         mphi, T_RH);
end

% -----------------------------------------------------------------------
% x range: integrate from T_RH to well past freeze-out
% -----------------------------------------------------------------------
info = struct('Y_phi_end', NaN, 'Omega_phi_h2', NaN, 'Gamma_phi', NaN, ...
              'phi_stable', false, 'overclosed', false, 'tau_phi', NaN);   % filled on success
x_init  = mx / T_RH;
x_final = max(1000, 10*x_init);
T_min   = mx / x_final;
T_max   = T_RH;

% -----------------------------------------------------------------------
% phi partial widths (T-independent; only K1/K2 thermal-average factor
% depends on T).  Compute once outside the ODE.
% -----------------------------------------------------------------------
Gphi_struct = phi_partial_widths(mphi, mx, gx, g, interaction);
Gamma_gg    = Gphi_struct.gg;       % phi -> gamma gamma (bath in eq -> inverse decay alive)
Gamma_xx    = Gphi_struct.xx;       % phi -> chi chibar  (sources Y_chi)
Gamma_meson = Gphi_struct.pipi + Gphi_struct.KK + Gphi_struct.eta;  % phi <-> meson pairs
mes_eq      = strcmp(production, 'sm');          % mesons in the bath: inverse decay alive

% phi <-> bath 2 -> 2 conversions (Primakoff, pi gamma -> pi phi, A' pi -> pi pi) on a fixed grid.
% LOOP_PHI_SC=0 in the environment turns them off (the pre-2026-11-12 solver).
Tsc = logspace(log10(5e-6), log10(T_RH), 120);
Gsc = zeros(size(Tsc));
use_sc = ~strcmp(getenv('LOOP_PHI_SC'), '0');
if use_sc && ~strcmp(interaction, 'axialVector')
    mref = max(mphi, 1e-12);
    if mphi < mref, Gwr = phi_partial_widths(mref, mx, gx, g, interaction); else, Gwr = Gphi_struct; end
    Gsc = Gsc + primakoff_phi_rate(Tsc, mphi, 64*pi*Gwr.gg/mref^3);
end
if use_sc && mes_eq && mphi < 2*0.13957 && ...
        ismember(interaction, {'heavyScalar','lightScalar','heavyScalarPS','lightScalarPS'})
    Gsc = Gsc + pigamma_phi_rate(Tsc, mphi, g);
end
% axialVector: A' pi -> pi pi (contact term of the 3pi current, crossed), the
% A''s tree-level link to the pion bath below 3 m_pi (apipi_phi_rate.m).
if use_sc && mes_eq && strcmp(interaction, 'axialVector') && ~strcmp(getenv('LOOP_APIPI'), '0')
    Gsc = Gsc + apipi_phi_rate(Tsc, mphi, g);
end
has_sc = any(Gsc > 0);  Gsc_interp = [];
if has_sc
    Gsc_interp = griddedInterpolant(log(Tsc(:)), log(max(Gsc(:), realmin)), 'pchip', 'nearest');
end

% -----------------------------------------------------------------------
% Pre-compute <sigma v>(T) on the integration grid -- both full and the
% SM-final-state piece (production), exactly as in freeze_loop.
% -----------------------------------------------------------------------
Tgrid    = logspace(log10(T_min), log10(T_max), 200);
logTgrid = log(Tgrid);
% Depletion into SM final states ('sm') and production (production channel);
% the dark pair's depletion comes from sv_darkpair_grid below, the same forward
% rate that feeds Y_phi, so chi chibar <-> phi phi conserves chi + phi exactly.
svvals_sm = svavg_grid(mx, mphi, gx, g, Tgrid, interaction, production);
if strcmp(production, 'sm')
    svvals_full = svvals_sm;
else
    svvals_full = svavg_grid(mx, mphi, gx, g, Tgrid, interaction, 'sm');
end
% Interpolated in log <sv>: the resonance and meson thresholds make <sv> fall as
% e^{-(mphi - 2mx)/T} or e^{-2m_pi/T}, which pchip in value misses by several %.
lfloor = log(realmin);
svFull = griddedInterpolant(logTgrid(:), log(max(svvals_full(:), realmin)), 'pchip');
svProd = griddedInterpolant(logTgrid(:), log(max(svvals_sm(:),   realmin)), 'pchip');
% chi chibar <-> phi phi (phi phi / A'A' / phi'phi'), both directions tabulated
% directly in log (sv_darkpair_grid): forward depletion feeds Y_phi, the reverse
% phi phi -> chi chibar is lambda <sv>_rev Y_phi^2, O(1) even when mphi > mx.
% Only the reverse table is kept: the forward rate is derived from it in logs,
% cd = crev (Y_phi_eq/Y_chi_eq)^2, so detailed balance holds to round-off (two
% separately interpolated tables leave ~1e-6 that the 1e17-1e19 rates turn into
% a spurious source the stiff solver cannot start from).
[~, lsvR] = sv_darkpair_grid(mx, mphi, gx, g, Tgrid, interaction, g_phi);
has_dark = any(isfinite(lsvR));
lsvR(~isfinite(lsvR)) = lfloor - 50;
svDarkR = griddedInterpolant(logTgrid(:), lsvR(:), 'pchip');
% axialVector: A'A' -> gamma gamma (charged-meson loop) and meson pairs, the
% A' number-changing link to the bath.  Below 3 m_pi
% it is the only one: Gamma(A' -> gamma gamma) = 0 (Landau-Yang).
% LOOP_AA_SM=0 in the environment turns it off (the pre-2026-11-17 solver).
has_AA = strcmp(interaction, 'axialVector') && ~strcmp(getenv('LOOP_AA_SM'), '0');
if has_AA
    lsvA = sv_AAsm_grid(mphi, g, Tgrid);
    has_AA = any(isfinite(lsvA));
    lsvA(~isfinite(lsvA)) = lfloor - 50;
    svAA = griddedInterpolant(logTgrid(:), lsvA(:), 'pchip');
end

% Precompute g_s, g_rho, d(ln g_s)/d(ln T) interpolants (C2 spline) over the
% integration T range.  Y_eq_mass and boltz_coef BOTH read from S so the
% WIMP quasi-steady IC cancellation stays consistent with the RHS.
S = gstar_interp(T_min, T_max);

% Adiabatic-phi flag: when Gamma_phi_total / H >> 1 throughout the integration
% range, Y_phi tracks Y_phi_eq instantly compared to Hubble.  The
% (Y_phi - Y_phi_eq) source term for Y_chi then vanishes and the system
% reduces to a scalar ODE for Y_chi alone.  Detecting this lets us drop the
% second component and avoid the O(1e6) step counts ode15s hits chasing
% Y_phi's tiny deviations when Gamma_total/H reaches 1e19+ at gx ~ 1.
H_RH_check    = T_RH^2/Mpl * sqrt(8*pi^3 * gstarofT_vec(T_RH,0) / 90);
Gamma_total   = Gamma_gg + Gamma_xx + Gamma_meson;
K1K2_RH_check = K1_over_K2(mphi/T_RH);
adiabatic_phi = (Gamma_total * K1K2_RH_check / H_RH_check > 1e3);
force_2d = strcmp(getenv('LOOP_FC_FULL2D'), '1');      % tests: always integrate (Y_chi, Y_phi)
if has_sc && ~adiabatic_phi
    % Gamma_total K1/K2 / H only grows as T falls, so without Gamma_sc the
    % T_RH test is the minimum; with it, take the minimum over the range.
    Tchk = logspace(log10(T_min), log10(T_max), 100);
    rchk = arrayfun(@(T) (Gamma_total*K1_over_K2(mphi/T) + Gamma_sc(T)) / hubble(T), Tchk);
    adiabatic_phi = min(rchk) > 1e3;
end
if force_2d, adiabatic_phi = false; end

% -----------------------------------------------------------------------
% Initial conditions
% -----------------------------------------------------------------------
% Y_chi starts at 0 at T_RH (no DM above reheating) and its first step is
% taken analytically, as in freeze_loop.m.  Over a step dtau0 short against
% the variation of Y_eq (~1/x per e-fold) dY/dtau = a - cf Y^2 with constant
% coefficients solves exactly:
%   Y(dtau0) = Y* tanh(sqrt(a cf) dtau0),   Y* = Y_eq sqrt(cp/cf),  a = cp Y_eq^2.
% The start is continuous in the coupling, so Omega h^2 is too.  (A rate test
% choosing Y = 0 or Y = Y* made Omega h^2(g_q) jump by up to ~500x near
% mx ~ 20 T_RH, which fzero reported as a root.)  The phi -> chi chibar
% source is left out of the seed; over dtau0 it is negligible.
%
% use_logY is now only a choice of variable: log-Y when the seed is already
% within a factor 2 of Y* (needed to track Y_eq over many decades in the
% deep-WIMP regime), linear-Y otherwise.  Both integrate the same ODE from
% the same state.
%
% Every initial condition below is built at x0 = exp(tau0), the x the RHS
% reconstructs at t0.  Built at any other x, the 1-ulp difference is a 1e-16
% in W_eq, which the annihilation term multiplies by cf*exp(W) ~ 1e31, so
% rhs(t0) lands ~1e15 per e-fold off an exact steady state and ode15s aborts
% on step 1 -> NaN.
tau0  = log(x_init) + 1e-3 / max(1, x_init);
dtau0 = tau0 - log(x_init);
x0    = exp(tau0);
% Y_phi IC: like chi, nothing before T_RH -- phi gets the same analytic first
% step over dtau0, Y_phi_eq (1 - exp(-Gamma_SM dtau0 / H)): zero unless its bath
% link is fast enough to equilibrate it within dtau0, and continuous in the
% couplings (a 0-or-equilibrium switch put a step into Omega h^2(g_q)).  A phi
% population present at T_RH would instead depend on reheating: on the P(x)S
% mphi = 0.1 GeV freeze-in branch above mphi/2, one Hubble time of gamma gamma
% -> phi before T_RH sets ~98% of the relic.
H_RH = hubble(T_RH);
K1K2_RH = K1_over_K2(mphi / T_RH);
Gsm_eff_RH = (Gamma_gg + mes_eq*Gamma_meson) * K1K2_RH;
if has_sc, Gsm_eff_RH = Gsm_eff_RH + Gamma_sc(T_RH); end
C0 = boltz_coef(x0);
% A'A' -> SM run backwards fills the A' at c_AA Y_phi_eq per e-fold (linearised)
Y_phi_0 = Y_eq_mass(x0, mphi, g_phi) * -expm1(-(Gsm_eff_RH / H_RH + C0.c_AA*C0.Yeq_phi) * dtau0);
Yeq_chi0 = C0.Yeq_chi;
if adiabatic_phi
    % Saturated seed: at rho = 1 the phi exchange vanishes, leaving cp_s, cf_s.
    % Linear seed: the production at Y_chi = 0, phi exchange included.
    cp0 = C0.cp_s;  cf0 = C0.cf_s;                % phi exchange vanishes at rho = 1
    T0  = chi_exch_adia(C0, 0);
    a_lin = max(cp0*Yeq_chi0^2 + T0, 0);
    % Fast dark exchange: with mphi > 2 mx, chi chibar -> phi phi -> 4 chi
    % multiplies chi at ~cd Y per e-fold, so any seed reaches equilibrium within
    % dtau0 -- start saturated (the linear seed leaves the solver to follow that
    % runaway and it gives up).
    [~, dTs] = chi_exch_adia(C0, a_lin*dtau0);
    dark_sat = a_lin > 0 && abs(dTs)*dtau0 > 1;
else
    dark_sat = false;
    [~, cf0, ~, lP0] = chi_coefs(C0, Y_phi_0);
    a_lin = exp(lP0);
    cp0 = exp(min(lP0 - 2*C0.lYeq_chi, 700));       % seed only (ratio, tanh argument)
end
if cf0 > 0
    ratio   = snap_unit_ratio(cp0 / cf0);
    seed_f  = tanh(sqrt(cp0*cf0) * Yeq_chi0 * dtau0);   % Y(dtau0)/Y*
    % A seed within 1e-10 of Y* is the steady state: snap it, so rhs(t0) is
    % exactly 0 there rather than ~1e-9 (which Julia's QNDF cannot start from;
    % kept identical here so the two languages share one IC).
    if 1 - seed_f < 1e-10, seed_f = 1; end
    use_logY = seed_f > 0.5;
    W_chi_0  = C0.lYeq_chi + 0.5*log(min(1, ratio)) + log(seed_f);
    % linear seed: Y* tanh with the production at Y_chi = 0 (= cp0 Y_eq^2 off the
    % fast-phi path), continuous with the saturated one where the remainder is small
    Y_chi_0  = sqrt(a_lin/cf0) * tanh(sqrt(a_lin*cf0) * dtau0);
else
    use_logY = false;
    Y_chi_0  = a_lin * dtau0;
end
% Fast phi whose exchange alone would carry chi past half its equilibrium within
% dtau0 (or the dark 2 -> 4 runaway above): start at the saturated seed.  With
% cf0 = 0 (the on-shell phi is all of cp, cf) the linear seed above is then
% production x dtau0, far beyond Y_eq.
if adiabatic_phi && (dark_sat || a_lin*dtau0 > 0.5*Yeq_chi0)
    use_logY = true;
    W_chi_0  = C0.lYeq_chi;
    if cf0 > 0, W_chi_0 = W_chi_0 + 0.5*log(min(1, snap_unit_ratio(cp0/cf0))); end
end

if use_logY
    Y0 = [W_chi_0; Y_phi_0];
else
    Y0 = [Y_chi_0; Y_phi_0];
end

% -----------------------------------------------------------------------
% Dark-exchange shortcut (closed or nearly closed dark sector)
% -----------------------------------------------------------------------
% When phi's own rates are slow (no fast-phi path) but chi chibar <-> phi phi
% runs far faster than H -- e.g. axialVector with mx > mphi/2, where Gamma(A')
% = 0 and the exchange reaches ~1e20 H -- chi and phi sit in relative chemical
% equilibrium, Y_phi = r Y_chi with r = Y_phi_eq/Y_chi_eq, and only the total
% N = Y_chi + Y_phi evolves slowly:
%   dN/dtau = cp_s Y_chi_eq^2 - cf_s Y_chi^2 + c_xx (Y_phi - Y_phi_eq rho^2)
%             - c_SM (Y_phi - Y_phi_eq) - c_mes Y_phi,   Y_chi = N/(1+r).
% Stage 1 integrates log N (the chi terms in ann_logY form with Y_eq ->
% Y_chi_eq (1+r)); an event stops it when the exchange falls below 1e3 per
% e-fold (dark freeze-out), and stage 2 continues the full 2x2 system from the
% split state.  The 2x2 system alone cannot start in this regime.
% LOOP_FC_NOEXCH=1 disables the shortcut (tests).
exch_fast = has_dark && ~adiabatic_phi && ~strcmp(getenv('LOOP_FC_NOEXCH'), '1') ...
            && C0.cd*Yeq_chi0 > 1e9 && C0.crev*C0.Yeq_phi > 1e9;
if exch_fast
    % N's own analytic first step: chi's SM terms around N's equilibrium
    % Y_chi_eq (1+r) (chi's seed treats the exchange as pure depletion)
    [~, ~, l1r0] = split_exch(C0, 0);
    if C0.cf_s > 0
        zN = sqrt(C0.cp_s*C0.cf_s) * Yeq_chi0 * dtau0;
        fN = tanh(zN);  if 1 - fN < 1e-10, fN = 1; end
        WN0 = C0.lYeq_chi + l1r0 + 0.5*log(min(1, snap_unit_ratio(C0.cp_s/C0.cf_s))) + log(max(fN, realmin));
    else
        WN0 = log(max(C0.cp_s*Yeq_chi0^2*dtau0 + Y_phi_0, realmin));
    end
    optsN = odeset('RelTol', 1e-6, 'AbsTol', 1e-10, 'Jacobian', @jac_exch, 'Events', @ev_exch);
    warnState = warning('off', 'MATLAB:ode15s:IntegrationTolNotMet');
    lastwarn('');
    [tN, WN, te] = ode15s(@rhs_exch, [tau0, log(x_final)], WN0, optsN);
    [~, warnId] = lastwarn;
    warning(warnState);
    if strcmp(warnId, 'MATLAB:ode15s:IntegrationTolNotMet') || ...
       (isempty(te) && tN(end) < log(x_final) - 1e-6)
        Omega_h2 = NaN;  return;
    end
    [Yc, Yp_] = split_exch(boltz_coef(exp(tN(end))), WN(end));
    if isempty(te)                                  % fast to the end
        Y_chi_now = Yc;  Y_phi_now = Yp_;
    else                                            % dark freeze-out: hand off
        opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-10, 'Jacobian', @jac_logY);
        warnState = warning('off', 'MATLAB:ode15s:IntegrationTolNotMet');
        lastwarn('');
        [tau2, Y2] = ode15s(@rhs_logY, [tN(end), log(x_final)], [log(Yc); Yp_], opts);
        [~, warnId] = lastwarn;
        warning(warnState);
        if tau2(end) < log(x_final) - 1e-6 || strcmp(warnId, 'MATLAB:ode15s:IntegrationTolNotMet')
            Omega_h2 = NaN;  return;
        end
        Y_chi_now = exp(Y2(end, 1));  Y_phi_now = Y2(end, 2);
    end
    Omega_h2 = mx * 2891.2 * Y_chi_now / 1.054e-5;
    info = phi_info(Y_phi_now);
    return;
end

% -----------------------------------------------------------------------
% Integrate in tau = log(x) with explicit 2x2 Jacobian
% -----------------------------------------------------------------------
% Tolerances: WIMP branch uses W = log Y (O(1)-O(100) range -> AbsTol=1e-10
% reasonable).  FIMP branch keeps tight AbsTol since Y can be very small.
if adiabatic_phi
    % Y_phi locked to Y_phi_eq; (Y_phi - Y_phi_eq) source vanishes.
    % Scalar state (W_chi if WIMP IC, Y_chi if FIMP).
    if use_logY
        opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-10, 'Jacobian', @jac_adia_logY);
        rhs_fn = @rhs_adia_logY;
        Y0s = Y0(1);
    else
        opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-18, 'Jacobian', @jac_adia_linY);
        rhs_fn = @rhs_adia_linY;
        Y0s = Y0(1);
    end
elseif use_logY
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-10, 'Jacobian', @jac_logY);
    rhs_fn = @rhs_logY;
    Y0s = Y0;
else
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-18, 'Jacobian', @jac_linY);
    rhs_fn = @rhs_linY;
    Y0s = Y0;
end
warnState = warning('off', 'MATLAB:ode15s:IntegrationTolNotMet');
lastwarn('');
[tau_sol, Y_sol] = ode15s(rhs_fn, [tau0, log(x_final)], Y0s, opts);
[~, warnId] = lastwarn;
warning(warnState);

reached_end = tau_sol(end) >= log(x_final) - 1e-6;
if ~reached_end || strcmp(warnId, 'MATLAB:ode15s:IntegrationTolNotMet')
    Omega_h2 = NaN;
    return;
end

if use_logY
    Y_chi_now = exp(Y_sol(end, 1));   % first column is W_chi (or scalar in adiabatic)
else
    Y_chi_now = Y_sol(end, 1);
end

if adiabatic_phi
    Y_phi_now = phi_qss(boltz_coef(x_final), Y_chi_now);
else
    Y_phi_now = Y_sol(end, 2);
end
info = phi_info(Y_phi_now);

s0          = 2891.2;     % cm^-3
rho_crit_h2 = 1.054e-5;   % GeV cm^-3
Omega_h2    = mx * s0 * Y_chi_now / rho_crit_h2;

% -----------------------------------------------------------------------
% Optional diagnostic plot
% -----------------------------------------------------------------------
if doplot
    x_sol = exp(tau_sol);
    if use_logY
        Y_chi_sol = exp(Y_sol(:,1));
    else
        Y_chi_sol = Y_sol(:,1);
    end
    if adiabatic_phi
        % Y_phi was its quasi-steady value throughout.
        Y_phi_sol = arrayfun(@(tt, yy) phi_qss(boltz_coef(exp(tt)), yy), tau_sol, Y_chi_sol);
    else
        Y_phi_sol = Y_sol(:,2);
    end
    x_plot = logspace(log10(x_init), log10(x_final), 600);
    Yeq_chi_plot = arrayfun(@(xx) Y_eq_mass(xx, mx,   g_chi), x_plot);
    Yeq_phi_plot = arrayfun(@(xx) Y_eq_mass(xx, mphi, g_phi), x_plot);

    figure('Name', 'Coupled (Y_\chi, Y_\phi) integration');
    loglog(x_sol, max(Y_chi_sol, 1e-40), 'b-', 'LineWidth', 2); hold on;
    loglog(x_sol, max(Y_phi_sol, 1e-40), 'g-', 'LineWidth', 2);
    loglog(x_plot, Yeq_chi_plot, 'b--', 'LineWidth', 1.2);
    loglog(x_plot, Yeq_phi_plot, 'g--', 'LineWidth', 1.2);
    xlabel('x = m_\chi / T'); ylabel('Y = n / s');
    legend('Y_\chi', 'Y_\phi', 'Y_{\chi,eq}', 'Y_{\phi,eq}', 'Location','best');
    title(sprintf('Coupled:  m_\\chi=%.4g, m_\\phi=%.4g, g_\\chi=%.4g, \\Omega h^2=%.3g', ...
                  mx, mphi, gx, Omega_h2));
    xlim([x_init x_final]); grid on;
end

% -----------------------------------------------------------------------
% Helpers
% -----------------------------------------------------------------------
    function I = phi_info(Yp)
        % Mediator left at x_final.  Omega_phi_h2 is what it would contribute
        % today if it never decays; phi_stable: lifetime beyond the age of the
        % universe (4.35e17 s = 6.61e41 GeV^-1).  overclosed: stable and
        % Omega_phi_h2 > 0.12 (the point is then excluded whatever chi does).
        I.Y_phi_end    = max(Yp, 0);
        I.Omega_phi_h2 = mphi * 2891.2 * I.Y_phi_end / 1.054e-5;
        I.Gamma_phi    = Gamma_total;
        I.phi_stable   = Gamma_total * 6.61e41 < 1;
        I.overclosed   = I.phi_stable && I.Omega_phi_h2 > 0.12;
        I.tau_phi = 6.582e-25 / Gamma_total;
    end
    function H = hubble(T)
        gR = gstarofT_vec(T, 0);
        % H = T^2/Mpl * sqrt(8 pi^3 g_rho / 90), radiation domination
        H = T^2 / Mpl * sqrt(8 * pi^3 * gR / 90);
    end

    function Yeq = Y_eq_mass(x, m, g_internal)
        % Y_eq for a species of mass m and internal dof g_internal, with the
        % time variable still x = mx/T (so the actual species argument is
        % m/T = (m/mx) * x).
        T_ = mx / x;
        g_s_ = S.gs(log(T_));
        z = m / T_;
        if z > 30
            K2_ = exp(-z + 0.5*log(pi/(2*z))) * (1 + 15/(8*z));
            z2K2 = z^2 * K2_;
        elseif z < 1e-3
            % Relativistic limit: z^2 K2(z) -> 2 as z->0 (massless species).
            % Evaluating besselk(2,0)=Inf directly gives 0*Inf = NaN.
            z2K2 = 2;
        else
            z2K2 = z^2 * besselk(2, z);
        end
        Yeq = (45/(4*pi^4)) * (g_internal / g_s_) * z2K2;
    end

    function G = Gamma_sc(T)
        % phi's 2 -> 2 conversion rate; 0 below the table (both rates are
        % e^{-m/T}-dead there) and with no conversion channel.
        if ~has_sc || T < Tsc(1), G = 0; return; end
        G = exp(Gsc_interp(log(min(T, Tsc(end)))));
        if G <= 2*realmin, G = 0; end
    end

    function r = K1_over_K2(z)
        % Thermal-average decay-rate suppression factor: <Gamma> = Gamma * K1/K2.
        % Asymptotic K_n(z) = sqrt(pi/(2z)) e^{-z} * series in 1/z; ratio K1/K2
        % is regular at all z (-> 1 as z->inf, -> z/2 as z->0).
        if z > 30
            r = (1 - 3/(8*z)) / (1 + 15/(8*z));
        elseif z < 1e-3
            r = z / 2;
        else
            r = besselk(1, z) / besselk(2, z);
        end
    end

    function C = boltz_coef(x)
        % Every rate coefficient at x, per e-fold in x.  Logs of the equilibrium
        % yields (lYc, lYp) and lp = log(Y_phi_eq/Y_chi_eq^2) keep the cross terms
        % finite where the yields themselves underflow.
        T     = mx / x;
        lT    = log(T);
        g_rho = S.grho(lT);
        g_s   = S.gs(lT);
        correction = 1 + S.dlngs(lT)/3;
        % Y counts chi + chibar (g_chi = 4) and <sigma v> is the chi chibar cross
        % section, so n-dot = -(<sigma v>/2)(n^2 - n_eq^2): the 1/2 is the Dirac factor.
        lam = 0.5 * sqrt(pi/45) * Mpl * mx / x * correction * (g_s / sqrt(g_rho));
        cp = lam * exp(svProd(lT));        % production channel
        cf = lam * exp(svFull(lT));        % depletion into SM final states

        H_val = T^2/Mpl * sqrt(8 * pi^3 * g_rho / 90);
        K1K2  = K1_over_K2(mphi / T);
        c_xx  = correction * Gamma_xx * K1K2 / H_val;
        c_SM  = correction * (Gamma_gg + mes_eq*Gamma_meson) * K1K2 / H_val;
        if has_sc, c_SM = c_SM + correction * Gamma_sc(T) / H_val; end
        c_mes = correction * (~mes_eq) * Gamma_meson * K1K2 / H_val;

        lYc = logY_eq_mass(x, mx,   g_chi);
        lYp = logY_eq_mass(x, mphi, g_phi);
        lp  = lYp - 2*lYc;
        if has_dark
            lsr  = svDarkR(lT);
            crev = lam * exp(min(lsr, 700));                     % phi phi -> chi chibar, on Y_phi^2
            cd   = lam * exp(min(lsr + 2*(lYp - lYc), 700));     % chi chibar -> phi phi, on Y_chi^2
        else
            cd = 0;  crev = 0;
        end
        % A'A' <-> SM, on (Y_phi^2 - Y_phi_eq^2): 2 lam = s/H, no Dirac 1/2
        c_AA = 0;
        if has_AA, c_AA = 2*lam * exp(min(svAA(lT), 700)); end
        % on-shell phi in cp / cf (detailed balance), removed: phi is explicit
        Gtot = Gamma_gg + Gamma_xx + Gamma_meson;
        if c_xx > 0 && Gtot > 0
            Ron = exp(min(log(2*c_xx) + lp, 700));
            R_p = Ron * (Gamma_gg + mes_eq*Gamma_meson) / Gtot;
            R_f = Ron * (Gamma_gg + Gamma_meson) / Gtot;
        else
            R_p = 0;  R_f = 0;
        end
        cp_s = max(cp - R_p, 0);
        cf_s = cp_s + max((cf - cp) - (R_f - R_p), 0);

        C = struct('cp_s', cp_s, 'cf_s', cf_s, 'cd', cd, 'crev', crev, 'c_xx', c_xx, ...
                   'c_SM', c_SM, 'c_mes', c_mes, 'c_AA', c_AA, 'lp', lp, 'lYeq_chi', lYc, 'lYeq_phi', lYp, ...
                   'Yeq_chi', exp(lYc), 'Yeq_phi', exp(lYp));
    end

    function lY = logY_eq_mass(x, m, g_internal)
        % log of Y_eq_mass, finite where Y_eq_mass underflows.
        T_ = mx / x;
        z = m / T_;
        if z > 30
            lz2K2 = 2*log(z) - z + 0.5*log(pi/(2*z)) + log1p(15/(8*z));
        elseif z < 1e-3
            lz2K2 = log(2);
        else
            lz2K2 = log(z^2 * besselk(2, z));
        end
        lY = log(45/(4*pi^4)) + log(g_internal / S.gs(log(T_))) + lz2K2;
    end

    function v = ex(l, Y, n)
        % exp(l) * Y^n for Y >= 0, in logs, capped; 0 for Y <= 0.
        if Y > 0
            v = exp(min(l + n*log(Y), 700));
        else
            v = 0;
        end
    end

    function Yp = phi_qss(C, Y_chi)
        % quasi-steady Y_phi for fast phi (diagnostic plot): Y_phi_eq + delta
        [~, ~, d] = chi_exch_adia(C, Y_chi);
        Yp = C.Yeq_phi + d;
    end

    function [cpE, cfE, P_phi, lP] = chi_coefs(C, Y_phi)
        % The chi equation as cpE Y_chi_eq^2 - cfE Y_chi^2 at fixed Y_phi:
        %   cpE = cp_s + [2 c_xx Y_phi + crev Y_phi^2] / Y_chi_eq^2
        %   cfE = cf_s + cd + 2 c_xx Y_phi_eq / Y_chi_eq^2.
        % P_phi = d(production)/dY_phi (Y units), for the Jacobians.
        lYc2 = -2*C.lYeq_chi;
        cpE = C.cp_s;
        if C.c_xx > 0, cpE = cpE + ex(log(2*C.c_xx) + lYc2, Y_phi, 1); end
        if C.crev > 0, cpE = cpE + ex(log(C.crev) + lYc2, Y_phi, 2); end
        cfE = C.cf_s + C.cd;
        if C.c_xx > 0, cfE = cfE + exp(min(log(2*C.c_xx) + C.lp, 700)); end
        P_phi = 2*C.c_xx + 2*C.crev*max(Y_phi, 0);
        % log of the production in Y units, cp_s Y_chi_eq^2 + 2 c_xx Y_phi
        % + crev Y_phi^2, as a log-sum: never divided by Y_chi_eq^2, which
        % underflows late while phi lingers (cpE above is capped there).
        l = -Inf(1, 3);
        if C.cp_s > 0, l(1) = log(C.cp_s) + 2*C.lYeq_chi; end
        if Y_phi > 0
            if C.c_xx > 0, l(2) = log(2*C.c_xx) + log(Y_phi); end
            if C.crev > 0, l(3) = log(C.crev) + 2*log(Y_phi); end
        end
        lm = max(l);
        if isfinite(lm), lP = lm + log(sum(exp(l - lm))); else, lP = -Inf; end
    end

    function dW = ann_logP(lP, cf, W, W_eq)
        % ann_logY with the production given as lP = log(cp Y_chi_eq^2): the
        % ratio cp/cf is exp(lP - 2 W_eq - log cf), formed in logs and snapped
        % like snap_unit_ratio.  Far from balance the two terms are evaluated
        % directly (nothing cancels).
        if cf > 0 && isfinite(lP)
            lr  = lP - 2*W_eq - log(cf);
            if abs(lr) <= 1e-12, lr = 0; end
            arg = lr - 2*(W - W_eq);
            if arg > 30
                dW = exp(lP - W) - cf*exp(W);
            else
                dW = cf*exp(W) * expm1(arg);
            end
        elseif isfinite(lP)
            dW = exp(lP - W);
        else
            dW = -cf*exp(W);
        end
    end

    function [Tphi, dT, delta] = chi_exch_adia(C, Y_chi)
        % Fast phi: Y_phi = Y_phi_eq + delta at its exact quasi-steady value,
        %   0 = -c_SM delta - c_xx (Y_phi - Y_phi_eq rho^2) - c_mes Y_phi
        %       + cd Y_chi^2 - crev Y_phi^2,
        % i.e. crev delta^2 + D' delta - B' = 0 with
        %   D' = c_SM + c_xx + c_mes + 2 crev Y_phi_eq,
        %   B' = (c_xx Y_phi_eq + crev Y_phi_eq^2)(rho^2 - 1) - c_mes Y_phi_eq
        %   (cd Y_chi^2 = crev Y_phi_eq^2 rho^2, exact by construction),
        %   delta = 2 B' / (D' + sqrt(D'^2 + 4 crev B')).
        % Everything chi exchanges with phi -- 2 c_xx (Y_phi - Y_phi_eq rho^2)
        % + crev Y_phi^2 - cd Y_chi^2 -- then follows from phi's balance as
        %   Tphi = -2 (c_SM + c_mes + crev Y_phi_eq) delta - 2 c_mes Y_phi_eq
        %          + crev Y_phi_eq^2 (rho^2 - 1) - crev delta^2.
        % The elastic loops (chi chibar -> phi -> chi chibar, -> phi phi -> chi
        % chibar) cancel analytically, every term vanishes at equilibrium through
        % rho^2 - 1 = expm1(.) (exactly 0 at W = W_eq), and nothing is linearised.
        % (Splitting this into a linear A Y_eq^2 - B Y^2 plus a remainder left two
        % ~1e19 pieces cancelling, and the solvers crawled: 1.6e6 steps.)
        % dT = dTphi/dY_chi, for the Jacobians.
        % A'A' <-> SM (c_AA, axialVector) adds c_AA (2 Y_phi_eq delta + delta^2):
        % quadratic coefficient q = crev + c_AA, D' += 2 c_AA Y_phi_eq, and in Tphi
        % G += 2 c_AA Y_phi_eq, crev delta^2 -> (crev + 2 c_AA) delta^2.
        Yp  = C.Yeq_phi;
        q   = C.crev + C.c_AA;
        D1  = C.c_SM + C.c_xx + C.c_mes + 2*q*Yp;
        A   = C.c_xx*Yp + C.crev*Yp^2;
        % A (rho^2 - 1) and crev Y_phi_eq^2 (rho^2 - 1): through expm1 near
        % equilibrium (exact 0 at rho = 1); far from it as A rho^2 - A with
        % A rho^2 = (c_xx Y_phi_eq/Y_chi_eq^2 + cd) Y_chi^2 built in logs --
        % rho^2 itself overflows once chi has frozen out (x ~ 400).
        Arho2 = 0;  Crho2 = 0;                               % A rho^2, crev Y_phi_eq^2 rho^2
        if Y_chi > 0
            if C.c_xx > 0, Arho2 = ex(log(C.c_xx) + C.lp, Y_chi, 2); end
            Crho2 = C.cd*Y_chi^2;
            Arho2 = Arho2 + Crho2;
        end
        lr = -Inf;  if Y_chi > 0, lr = 2*(log(Y_chi) - C.lYeq_chi); end
        if abs(lr) <= 1
            Ar = A*expm1(lr);  Cr = C.crev*Yp^2*expm1(lr);
        else
            Ar = Arho2 - A;    Cr = Crho2 - C.crev*Yp^2;
        end
        % delta = Y_phi - Y_phi_eq from crev delta^2 + D' delta - B' = 0
        Bp = Ar - C.c_mes*Yp;
        if D1 > 0
            if q > 0
                sq = sqrt(max(D1^2 + 4*q*Bp, 0));
                delta = 2*Bp / (D1 + sq);
                fp = 1/D1;
                if sq > 0, fp = (2*(D1 + sq) - 4*q*Bp/sq) / (D1 + sq)^2; end
            else
                delta = Bp / D1;  fp = 1/D1;
            end
        else
            delta = 0;  fp = 0;
        end
        if abs(delta) <= 0.5*Yp
            % phi near its equilibrium (whatever chi does): delta form, every term
            % 0 at rho = 1; the large c_SM multiplies the small delta directly
            G    = C.c_SM + C.c_mes + C.crev*Yp + 2*C.c_AA*Yp;
            q2   = C.crev + 2*C.c_AA;
            Tphi = -2*G*delta - 2*C.c_mes*Yp + Cr - q2*delta^2;
            dT = 0;
            if Y_chi > 0, dT = (-2*G - 2*q2*delta)*fp*(2*Arho2/Y_chi) + 2*Crho2/Y_chi; end
        else
            % phi largely empty: Y_phi itself, root of crev Y^2 + D Y - Q = 0 with
            % Q = c_SM Y_phi_eq + A rho^2 -- no subtraction (with delta ~ -Y_phi_eq
            % the delta form cancels two ~1e16 terms and the production goes negative)
            D  = C.c_SM + C.c_xx + C.c_mes;
            Q  = C.c_SM*Yp + Arho2 + C.c_AA*Yp^2;
            sq = sqrt(D^2 + 4*q*Q);
            if D + sq > 0, Yph = 2*Q / (D + sq); else, Yph = Yp; end
            delta = Yph - Yp;
            Tphi = -2*C.c_SM*delta - 2*C.c_mes*Yph + Crho2 - C.crev*Yph^2 ...
                   - 2*C.c_AA*(Yph^2 - Yp^2);
            dT = 0;
            if Y_chi > 0 && sq > 0                       % dY_phi/dQ = 1/sq
                dT = (-2*C.c_SM - 2*C.c_mes - 2*(C.crev + 2*C.c_AA)*Yph)*(2*Arho2/Y_chi)/sq ...
                     + 2*Crho2/Y_chi;
            end
        end
    end

    function f = phi_rhs(C, Y_chi, Y_phi)
        f = -C.c_SM*(Y_phi - C.Yeq_phi) - C.c_xx*(Y_phi - ex(C.lp, Y_chi, 2)) - C.c_mes*Y_phi ...
            + C.cd*Y_chi^2 - C.crev*max(Y_phi, 0)^2 - C.c_AA*(max(Y_phi, 0)^2 - C.Yeq_phi^2);
    end

    function r = snap_unit_ratio(r)
        % <sv>_prod / <sv>_full, snapped to exactly 1 when the two agree to
        % round-off.  Whenever no dark channel is open, every chi chi channel
        % is an SM one, so the two are the same quantity -- but they are built
        % by different code paths (svavg_grid with and without the 'sm' filter,
        % then two separate pchip fits), so they differ in the last ulp.
        %
        % That 1e-16 is not negligible here.  The QSS initial condition carries
        % it as 0.5*log(ratio) ~ 5e-16 in W, and the annihilation term
        % multiplies it by cf*exp(W) ~ 1e31 (Gamma/H ~ 1e29 with a
        % Breit-Wigner-saturated <sv>), giving rhs(t0) ~ 1e15 per e-fold out of
        % what should be an exact steady state.  ode15s cannot resolve that
        % layer and aborts on step 1 -> NaN.  Whether a given mass lands on an
        % exact 0 or on 1e15 is decided in the last bit, which is why the
        % failures were isolated and ragged in m_chi rather than forming a
        % region.  1e-12 is ~1e4 ulp yet ~1e11 below any real channel
        % difference (S(x)S runs at ratio ~ 0.85).
        if abs(r - 1) <= 1e-12, r = 1; end
    end

    function dW = ann_logY(cp, cf, W, W_eq)
        % Annihilation part of dW/dtau in log-Y form, written so that the
        % near-equilibrium limit does not lose precision to cancellation.
        %
        % The direct form  cp*exp(2 W_eq - W) - cf*exp(W)  is a difference of
        % two terms that are individually enormous: cf ~ 1e32 when Gamma/H
        % reaches 1e29 and <sv> is Breit-Wigner saturated, so each term is
        % ~1e31 while the physical result near equilibrium is ~0.  For S(x)P
        % and P(x)P, gamma gamma is the ONLY chi chi -> SM channel, hence
        % cp == cf identically and at W = W_eq the whole RHS is the roundoff
        % residual of that subtraction -- ~1e15 per e-fold, i.e. a few ulp of
        % 1e31.  ode15s cannot integrate a derivative that IS noise: it aborts
        % on step 1 with IntegrationTolNotMet and the caller sees NaN.  Whether
        % a given mass lands on an exact 0 or on 1e15 is decided in the last
        % bit, which is why the failures were isolated and ragged in m_chi.
        %
        % Factoring out cf*exp(W) leaves expm1 of the log-space imbalance:
        %   dW = cf e^W [ (cp/cf) e^{-2(W - W_eq)} - 1 ]
        %      = cf e^W expm1( log(cp/cf) - 2(W - W_eq) ),
        % which is exactly 0 when cp == cf and W == W_eq, and accurate to full
        % precision for small imbalances instead of losing all of it.
        % cp == 0 is handled by log(0) = -Inf -> expm1 = -1 -> dW = -cf e^W.
        % The eliminated fast-phi coefficients can have either sign (chi chibar
        % -> phi phi -> 4 chi makes the net Y^2 term negative): then there is no
        % cancellation to protect and the direct form is used.
        if cf > 0 && cp >= 0
            dW = cf*exp(W) * expm1(log(snap_unit_ratio(cp/cf)) - 2*(W - W_eq));
        else
            dW = cp*exp(2*W_eq - W) - cf*exp(W);
        end
    end

    % ---- dark-exchange shortcut: state W_N = log(Y_chi + Y_phi) --------------
    function [Yc, Yp_, l1r, lrr] = split_exch(C, WN)
        % relative chemical equilibrium: Y_phi = r Y_chi, r = Y_phi_eq / Y_chi_eq
        lrr = C.lYeq_phi - C.lYeq_chi;                      % log r
        if lrr > 30, l1r = lrr + log1p(exp(-lrr)); else, l1r = log1p(exp(lrr)); end
        Yc  = exp(WN - l1r);
        Yp_ = exp(WN + lrr - l1r);
    end

    function dW = rhs_exch(tau, WN)
        C = boltz_coef(exp(tau));
        [Yc, Yp_, l1r] = split_exch(C, WN);
        % chi's SM terms with Y_eq -> Y_chi_eq (1+r): cp_s Y_eq^2 - cf_s (N/(1+r))^2
        dW = ann_logY(C.cp_s*exp(-2*l1r), C.cf_s*exp(-2*l1r), WN, C.lYeq_chi + l1r);
        % phi's (slow) own links change N too
        X = C.c_xx*(Yp_ - ex(C.lp, Yc, 2)) - C.c_SM*(Yp_ - C.Yeq_phi) - C.c_mes*Yp_ ...
            - C.c_AA*(Yp_^2 - C.Yeq_phi^2);
        dW = dW + X*exp(-WN);
    end

    function J = jac_exch(tau, WN)
        C = boltz_coef(exp(tau));
        [Yc, Yp_, l1r] = split_exch(C, WN);
        J = -C.cp_s*exp(2*C.lYeq_chi - WN) - C.cf_s*exp(WN - 2*l1r);
        % phi's own links, d(X e^-WN)/dWN = (dX/dWN - X) e^-WN: the SM link can be
        % fast (Gamma_sc ~ 1e2 H) without phi being eliminated, and is stiff then
        Xc = ex(C.lp, Yc, 2);
        X  = C.c_xx*(Yp_ - Xc) - C.c_SM*(Yp_ - C.Yeq_phi) - C.c_mes*Yp_ ...
             - C.c_AA*(Yp_^2 - C.Yeq_phi^2);
        dX = C.c_xx*(Yp_ - 2*Xc) - (C.c_SM + C.c_mes)*Yp_ - 2*C.c_AA*Yp_^2;
        J = J + (dX - X)*exp(-WN);
    end

    function [v, term, dirn] = ev_exch(tau, WN)
        % exchange rate per e-fold, min(cd Y_chi, crev Y_phi), falls below 1e3
        C = boltz_coef(exp(tau));
        [Yc, Yp_] = split_exch(C, WN);
        v = log(max(min(C.cd*Yc, C.crev*Yp_), realmin)) - log(1e3);
        term = true;  dirn = -1;
    end

    % ---- full 2x2 system --------------------------------------------------
    % log-Y branch: state (W_chi, Y_phi)
    function dY = rhs_logY(tau, Y)
        C = boltz_coef(exp(tau));
        [~, cfE, ~, lP] = chi_coefs(C, Y(2));
        dY = [ann_logP(lP, cfE, Y(1), C.lYeq_chi);
              phi_rhs(C, exp(Y(1)), Y(2))];
    end

    function J = jac_logY(tau, Y)
        C = boltz_coef(exp(tau));
        [~, cfE, P_phi, lP] = chi_coefs(C, Y(2));
        W = Y(1);
        J11 = -exp(lP - W) - cfE*exp(W);
        J12 =  P_phi*exp(-W);
        J21 =  2*C.c_xx*ex(C.lp, exp(W), 2) + 2*C.cd*exp(2*W);
        J22 = -(C.c_SM + C.c_xx + C.c_mes) - 2*(C.crev + C.c_AA)*max(Y(2), 0);
        J = [J11, J12; J21, J22];
    end

    % linear branch: state (Y_chi, Y_phi)
    function dY = rhs_linY(tau, Y)
        C = boltz_coef(exp(tau));
        [~, cfE, ~, lP] = chi_coefs(C, Y(2));
        dY = [exp(lP) - cfE*Y(1)^2;
              phi_rhs(C, Y(1), Y(2))];
    end

    function J = jac_linY(tau, Y)
        C = boltz_coef(exp(tau));
        [~, cfE, P_phi] = chi_coefs(C, Y(2));
        J21 = 2*C.c_xx*ex(C.lp, Y(1), 1) + 2*C.cd*Y(1);
        J22 = -(C.c_SM + C.c_xx + C.c_mes) - 2*(C.crev + C.c_AA)*max(Y(2), 0);
        J = [-2*cfE*Y(1), P_phi; J21, J22];
    end

    % ---- fast phi: Y_phi eliminated (chi_exch_adia), scalar chi equation ----
    function dW = rhs_adia_logY(tau, W)
        C = boltz_coef(exp(tau));
        Tphi = chi_exch_adia(C, exp(W));
        dW = ann_logY(C.cp_s, C.cf_s, W, C.lYeq_chi) + Tphi*exp(-W);
    end

    function J = jac_adia_logY(tau, W)
        C = boltz_coef(exp(tau));
        [Tphi, dT] = chi_exch_adia(C, exp(W));
        J = -C.cp_s*exp(2*C.lYeq_chi - W) - C.cf_s*exp(W) + dT - Tphi*exp(-W);
    end

    function dY = rhs_adia_linY(tau, Y_chi)
        C = boltz_coef(exp(tau));
        dY = C.cp_s*C.Yeq_chi^2 - C.cf_s*Y_chi^2 + chi_exch_adia(C, Y_chi);
    end

    function J = jac_adia_linY(tau, Y_chi)
        C = boltz_coef(exp(tau));
        [~, dT] = chi_exch_adia(C, Y_chi);
        J = -2*C.cf_s*Y_chi + dT;
    end

end
