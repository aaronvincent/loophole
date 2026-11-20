function sv = svavg(mx,mphi,gx,g,T,interaction,channel,moment)
% svavg(mx,mphi,gx,g,T,interaction,channel,moment)
%
% Thermally averaged <sigma v> for chi-chibar annihilation.
% channel (optional, default 'all'):
%   'all' -- full cross-section from xxggxsec.m (every open annihilation channel)
%   'gg'  -- s-channel chi-chibar -> gamma gamma only
%   'sm'  -- every SM final state (gamma gamma + meson pairs).  This is the
%            PRODUCTION side of the relic solvers: the SM bath is in
%            equilibrium, so n_eq^2 sv_{chi chi -> X} is the X X -> chi chi
%            rate.  Mesons are Boltzmann-suppressed, but for sqrt(s) above
%            2 m_pi so are the photon pairs, and the tree-level meson channels
%            then dominate.  Dark final states (phi phi, A'A') are excluded.
%   'had' -- meson pairs only

if nargin < 7 || isempty(channel)
    channel = 'all';
end
% moment (optional, default 'number'):
%   'number'  the usual <sigma v>            [GeV^-2]
%   'energy'  <sigma v E_tot>                [GeV^-1]
% The energy moment weights each annihilation by the lab-frame total energy of
% the annihilating pair, whose mean at fixed s is
%       <E_tot>_s = sqrt(s) K0(sqrt(s)/T)/K1(sqrt(s)/T) + 2T,
% so the kernel K1 is replaced by sqrt(s) K0 + 2 T K1.  The measure that
% produces this is the same Gondolo-Gelmini one that gives the K1 kernel in the
% first place, and it is unchanged by the initial-state mass -- the m-dependence
% factorises into the s-dependent part -- so the substitution is identical here
% and in ggxxrate.m.  n^2 <sigma v E> is then an ENERGY density per unit time,
% (Kept for symmetry with svavg_phiphi.m; the relic solvers use 'number'.)
if nargin < 8 || isempty(moment)
    moment = 'number';
end
if ~ismember(moment, {'number','energy'})
    error('svavg:moment', 'moment must be ''number'' or ''energy''');
end

mpi0  = 0.13498;  mpi_pm = 0.13957;
mK    = 0.493677;
m_eta = 0.547862;
B0      = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV (draft convention; was 2.6)
alphaEM = 1/137;
Dscript = 0.5*(B0*(g(1)+g(2))/mpi_pm^2 + B0*(g(1)+g(3))/mK^2)^2;   % charged-meson loop (pi+, K+), low-energy limit A0 -> 1/3

% Cross-section evaluator + phi total width.  Single source of truth:
% make_xxsigma.m (shared with the vectorized fast path svavg_grid.m), so the
% gg-only branch stays detailed-balance consistent with ggxxxsec/ggxxrate.
[sig, Gphi] = make_xxsigma(mx, mphi, gx, g, interaction, channel);
% Interactions with no gg-channel implementation return sig == 0 everywhere.
if strcmp(channel,'gg') && ~ismember(interaction, ...
        {'heavyScalar','lightScalar','heavyScalarPS','lightScalarPS', ...
         'heavyScalarSP','lightScalarSP','scalarPP','heavyScalarPP','lightScalarPP', ...
         'axialVector'})
    sv = 0;
    return
end

% Gondolo-Gelmini thermal-average kernel for Dirac chi-chibar annihilation.
% Use scaled Bessels K_n(z)*exp(z) (besselk(...,1)) and absorb exponentials
% into the integrand so neither the K2(mx/T)^2 denominator nor the K1(sqrt(s)/T)
% integrand underflow at low T (mx/T > ~700 otherwise hits double precision
% underflow and returns NaN).  exp((2 mx - sqrt(s))/T) <= 1 across the
% integration range (sqrt(s) >= 2 mx at threshold), so this is well-behaved.
K2scaled = besselk(2, mx/T, 1);
integrand    = @(s)  sig(s).*(s-4*mx^2).*sqrt(s) ...
                     .* svkernel(sqrt(s), T, moment) .* exp((2*mx - sqrt(s))/T);
logintegrand = @(ls) integrand(exp(ls)).*exp(ls);

s_lo  = 4*mx^2;
% Cover thermally relevant s plus the resonance neighbourhoods.  Matches the
% structure of ggxxrate.m so the two integrators agree by detailed balance.
s_max = max([100*T^2, 100*mphi^2, 100*mx^2, 1e2]);

% s-channel resonance poles (s_p = m_p^2, width w_p = m_p*Gamma_p) for this
% interaction.  thermal_pole_integral notches the narrow ones (NWA) and
% straddles the resolvable ones with waypoints, so the meson and mediator
% resonances are never silently missed -- in particular the narrow single-
% meson poles (pi0/eta/eta') that the old waypoint-only path skipped for
% sub-meson mx in the SP/PP/axialVector cross sections.
[pole_s, pole_w] = sv_poles(mx, mphi, gx, g, interaction, Gphi, channel);

% Meson-pair / phiphi / A'A' thresholds kink the 'all' cross section; pass
% them as continuum waypoints so the adaptive rule resolves the cusps.
m_etap   = 0.95778;
extra_wp = [4*mpi_pm^2, 4*mpi0^2, 4*m_eta^2, 4*mK^2, 4*0.497611^2, (mpi0+m_eta)^2, 4*mphi^2, 4*m_etap^2];
if ismember(interaction, {'scalarPP','heavyScalarPP','lightScalarPP','heavyScalarSP','lightScalarSP'})
    if threepi_enabled()
        extra_wp = [extra_wp, 9*mpi0^2, (2*mpi_pm + mpi0)^2];          % 3pi thresholds
    end
    if pipiphi_enabled()
        extra_wp = [extra_wp, (2*mpi0 + mphi)^2, (2*mpi_pm + mphi)^2]; % pi pi phi thresholds
    end
elseif strcmp(interaction, 'axialVector')
    if threepi_enabled()
        extra_wp = [extra_wp, 9*mpi0^2, (2*mpi_pm + mpi0)^2];          % 3pi thresholds
    end
    if pipiphi_enabled()
        extra_wp = [extra_wp, (2*mpi_pm + mphi)^2];                     % pi+pi-A' threshold
    end
end

% RelTol 1e-5 / AbsTol 1e-20: sv feeds a pchip interpolant consumed by ode15s
% at its own RelTol; tighter just burns wall time the downstream can't see.
opts = struct('RelTol', 1e-5, 'AbsTol', 1e-20, 'narrow_tol', 1e-5, ...
              'extra_wp', extra_wp);
kernel_at = @(sp) (sp - 4*mx^2).*sqrt(sp) ...
                  .* svkernel(sqrt(sp), T, moment) .* exp((2*mx - sqrt(sp))/T);
I = thermal_pole_integral(logintegrand, sig, kernel_at, s_lo, s_max, ...
                          pole_s, pole_w, opts);
sv = I / (8*mx^4*T*K2scaled.^2);
end

% --- s-channel resonance poles per interaction ----------------------------
function [pole_s, pole_w] = sv_poles(mx, mphi, gx, g, interaction, Gphi, channel)
mpi0 = 0.13498; m_eta = 0.547862;
Gpi = 7.7e-9; G_eta = 1.31e-6;
switch interaction
    case {'heavyScalar','lightScalar','heavyScalarPS','lightScalarPS'}
        pole_s = mphi^2;            pole_w = mphi*Gphi;
    case {'heavyScalarSP','lightScalarSP','scalarPP','heavyScalarPP','lightScalarPP'}
        % pi0, phi', eta mixing poles -- or eta', phi' for g ~ [1 1 1].
        [mj, ~, ~, Gam] = sppp_couplings(mphi, mx, gx, g, interaction);
        pole_s = mj.^2;     pole_w = mj.*Gam;
    case 'axialVector'
        % pi0, eta, eta' SM mixing poles (widths include Gamma(j->chi chi))
        [mj, ~, ~, Gam_j] = axial_couplings(mphi, mx, gx, g);
        pole_s = mj.^2;             pole_w = mj.*Gam_j;
        % plus the A' pole of the J = 1 chi chi -> 3pi channel (axial3pi_sigma),
        % in the channels that carry it
        if threepi_enabled() && mphi > 3*mpi0 && ismember(channel, {'all','sm','had'})
            pole_s = [pole_s, mphi^2];  pole_w = [pole_w, mphi*Gphi];
        end
    otherwise
        pole_s = [];                pole_w = [];
end
end

% =========================================================================
function K = svkernel(sq, T, moment)
% Bessel kernel of the Gondolo-Gelmini s-integral, for either moment.
% sq = sqrt(s).  SCALED Bessels (besselk(...,1)), matching the rest of svavg --
% the exp((2mx - sqrt(s))/T) factor outside compensates, and the K0/K1 scalings
% are the same exp(+x), so the combination below is consistent.
switch moment
    case 'number'
        K = besselk(1, sq./T, 1);
    case 'energy'
        K = sq .* besselk(0, sq./T, 1) + 2*T .* besselk(1, sq./T, 1);
    otherwise
        error('svavg:moment', 'moment must be ''number'' or ''energy''');
end
end
