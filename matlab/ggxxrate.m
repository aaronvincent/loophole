function G = ggxxrate(TT,mx,mphi,gx,g,interaction,moment)
% Thermally averaged gamma gamma -> chi chi rate (per the FI normalization).
% Single-pole interactions take the fast vectorized path; the multi-pole
% interactions (SP/PP/axialVector) use the per-T thermal_pole_integral below.
%
%   G = ggxxrate(TT,mx,mphi,gx,g,interaction)            % 'number' (default)
%   G = ggxxrate(TT,mx,mphi,gx,g,interaction,'energy')
%
% MOMENT.  'number' returns the reaction density R (reactions per volume per
% time), the zeroth moment, which is what the relic solvers consume.  'energy'
% returns the ENERGY injection rate into the chi sector -- the first moment of
% the same collision term.  The two differ only by the Bessel kernel:
%
%   number:  int ds sigma(s) s^(3/2) K1(sqrt(s)/T)
%   energy:  int ds sigma(s) s^(3/2) [ sqrt(s) K0(sqrt(s)/T) + 2 T K1(sqrt(s)/T) ]
%
% because every reaction injects the lab-frame energy of the incoming photon
% pair, whose mean at fixed s is <E_tot>_s = sqrt(s) K0/K1 + 2T (the total-energy
% distribution at fixed s is proportional to sqrt(E^2-s) exp(-E/T) for
% E >= sqrt(s)).  Same integrand, same poles, same NWA residues -- only the
% kernel changes, which is why this is a parameter here rather than a second
% copy of the pole bookkeeping.  The relic solvers use 'number'.
if nargin < 7 || isempty(moment), moment = 'number'; end
if ~ismember(moment, {'number','energy'})
    error('ggxxrate:moment', 'moment must be ''number'' or ''energy''');
end
G = zeros(1, length(TT));

% Fast path: single-Breit-Wigner-pole interactions (gamma gamma -> chi chi
% via one s-channel phi).  sigma(s) is T-independent -> evaluate once on a
% fixed Gauss-Legendre grid and reduce against the Bessel kernel for all T.
if ismember(interaction, {'heavyScalar','lightScalar','heavyScalarPS','lightScalarPS'})
    G = ggxxrate_fast(TT, mx, mphi, gx, g, interaction, moment);
    return
end

% All remaining interactions are multi-pole (SP / PP / axialVector): the
% gamma gamma -> chi chi amplitude runs through several s-channel resonances
% (SP/PP: pi0, phi', eta; axialVector: pi0, eta, eta').  Build the pole list
% once, then integrate each T with the shared thermal_pole_integral, which
% notches narrow poles (NWA) and straddles resolvable ones with waypoints.
% This captures every narrow meson pole; the NWA residue is derived
% consistently from sig(s_p)*prop2_p.
[pole_s, pole_w] = gg_poles(mx, mphi, gx, g, interaction);

% Plain `for`: nested parfor inside compare_limits_g's outer parfor degrades to
% serial anyway, and standalone calls (e.g. diagnostics) trip on stale worker
% source caches.  No throughput lost since the outer parfor already parallelises.
for k = 1:length(TT)
    T = TT(k);
    sig          = @(s)  ggxxxsec(mx, mphi, gx, g, s, interaction);
    sigintegrand = @(s)  sig(s).*s.^(3/2).*bkernel(sqrt(s), T, moment);
    logintegrand = @(ls) sigintegrand(exp(ls)).*exp(ls);
    % The kernel cuts off the integrand for s >> T^2 in both moments (the
    % 'energy' one only trades K1 for a combination with K0, same exponential);
    % cover a few decades above max(T, mphi, mx)^2 -- the mx term keeps a heavy
    % chi at low T from putting the threshold 4*mx^2 above s_max.
    s_max     = max([100*T^2, 100*mphi^2, 100*mx^2, 1e6]);
    kernel_at = @(sp) sp.^(3/2).*bkernel(sqrt(sp), T, moment);
    opts = struct('RelTol', 1e-6, 'AbsTol', 1e-20, 'narrow_tol', 1e-5);
    I = thermal_pole_integral(logintegrand, sig, kernel_at, 4*mx^2, s_max, ...
                              pole_s, pole_w, opts);
    G(k) = T./16./pi^4.*I;
end
end

% --- s-channel resonance poles for gamma gamma -> chi chi -----------------
function [pole_s, pole_w] = gg_poles(mx, mphi, gx, g, interaction)
mpi0 = 0.13498; m_eta = 0.547862;
Gpi = 7.7e-9; G_eta = 1.31e-6;
switch interaction
    case {'heavyScalarSP','lightScalarSP','scalarPP','heavyScalarPP','lightScalarPP'}
        % {pi0, phi', eta} generically, or {eta', phi'} for g ~ [1 1 1].
        [mj, ~, ~, Gam] = sppp_couplings(mphi, mx, gx, g, interaction);
        pole_s = mj.^2;     pole_w = mj.*Gam;
    case 'axialVector'
        [mj, ~, ~, Gam_j] = axial_couplings(mphi, mx, gx, g);
        pole_s = mj.^2;     pole_w = mj.*Gam_j;
    otherwise
        pole_s = [];        pole_w = [];
end
end

% =========================================================================
function G = ggxxrate_fast(TT, mx, mphi, gx, g, interaction, moment)
% Vectorized gamma gamma -> chi chi rate for single-BW-pole interactions.
% Mirrors the scalar path's continuum + narrow-width-approximation logic
% (see the main function body) but evaluates sigma(s) once on a fixed
% composite Gauss-Legendre grid in log(s) and reduces over all T at once.
TT   = TT(:).';
Gphi = phi_partial_widths(mphi, mx, gx, g, interaction).total;

s_lo  = 4*mx^2;
s_max = max([100*max(TT)^2, 100*mphi^2, 100*mx^2, 1e6]);
L_lo  = log(s_lo);
L_hi  = log(s_max);
L_mphi2 = log(mphi^2);

rel_w   = Gphi / max(mphi, eps);
pole_in = (mphi^2 > s_lo) && (mphi^2 < s_max) && (mphi^2 > 4*mx^2);
use_nwa = pole_in && (Gphi > 0) && (rel_w < 1e-6);

logdelta = 1e-3;
if use_nwa
    notch_lo = L_mphi2 - logdelta;
    notch_hi = L_mphi2 + logdelta;
    L_extra  = [notch_lo, notch_hi];
elseif pole_in
    fwhm = mphi * Gphi;
    cand = mphi^2 + [-100 -30 -5 -1 1 5 30 100]*fwhm;
    cand = cand(cand > s_lo & cand < s_max);
    L_extra = log(cand);
else
    L_extra = [];
end

edges = loggrid_edges(L_lo, L_hi, L_extra);
[xb, wb] = gl_nodes_weights(24);
ls = []; wl = [];
for j = 1:numel(edges)-1
    a = edges(j); b = edges(j+1);
    if b <= a, continue; end
    if use_nwa && abs(a-notch_lo) < 1e-12 && abs(b-notch_hi) < 1e-12
        continue;
    end
    mid = (a+b)/2; half = (b-a)/2;
    ls = [ls; mid + half*xb];
    wl = [wl; half*wb];
end

s_nodes = exp(ls);
sq      = sqrt(s_nodes);
sigvals = real(ggxxxsec(mx, mphi, gx, g, s_nodes, interaction));
sigvals = sigvals(:);

% integrand(s) = sig * s^(3/2) * kernel(sqrt(s),T); in log s pick up a factor s.
base = wl .* sigvals .* s_nodes.^(2.5);          % N x 1
Trow = TT;                                       % 1 x M
KK   = bkernel(sq, Trow, moment);                % N x M (unscaled, matches scalar path)
I    = sum(base .* KK, 1);                       % 1 x M

if use_nwa
    % Narrow-width residue: the same delta-function limit, with the kernel
    % evaluated at the pole.  s^(3/2) -> mphi^3 is already factored out below,
    % so only the Bessel combination differs between the two moments.
    prop2_pole = (mphi*Gphi)^2;
    sig_num    = real(ggxxxsec(mx, mphi, gx, g, mphi^2, interaction)) * prop2_pole;
    I = I + pi * sig_num * mphi^3 .* bkernel(mphi, Trow, moment) / (mphi*Gphi);
end

G = TT ./ 16 ./ pi^4 .* I;
end

% =========================================================================
function K = bkernel(sq, T, moment)
% Bessel kernel of the thermal s-integral, for either moment.  sq = sqrt(s);
% sq and T broadcast (sq column, T row) in the vectorized path.
%
%   'number'  K1(sqrt(s)/T)
%   'energy'  sqrt(s) K0(sqrt(s)/T) + 2 T K1(sqrt(s)/T)  =  <E_tot>_s * K1
%
% Unscaled Bessels in both cases, matching the rest of ggxxrate.
switch moment
    case 'number'
        K = besselk(1, sq ./ T);
    case 'energy'
        K = sq .* besselk(0, sq ./ T) + 2 .* T .* besselk(1, sq ./ T);
    otherwise
        error('ggxxrate:moment', 'moment must be ''number'' or ''energy''');
end
end
