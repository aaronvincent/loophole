function sv = svavg_grid(mx, mphi, gx, g, Tgrid, interaction, channel)
% svavg_grid  Vectorized <sigma v>(T) over a whole temperature grid.
%
% Same physics as svavg.m, but computes the entire Tgrid in one shot using a
% fixed composite Gauss-Legendre rule in log(s).  Because sigma(s) is
% T-independent it is evaluated ONCE on the node set and reused for every T;
% the per-T cost collapses to a Bessel evaluation + a matrix reduction.
%
% Returns a row vector sv the same length as Tgrid (units GeV^-2).
%
% Fast path covers the single-Breit-Wigner-pole interactions
% (heavyScalar/lightScalar/heavyScalarPS/lightScalarPS).  Multi-pole SP/PP
% interactions fall back to a per-T loop over svavg.m (correct, just slower).

if nargin < 7 || isempty(channel), channel = 'all'; end
Tgrid = Tgrid(:).';                 % row
sv    = zeros(1, numel(Tgrid));

[sig, Gphi, single_pole] = make_xxsigma(mx, mphi, gx, g, interaction, channel);

% ---- Fallback for multi-pole interactions (pi0/eta/phi resonances) -------
if ~single_pole
    for k = 1:numel(Tgrid)
        sv(k) = svavg(mx, mphi, gx, g, Tgrid(k), interaction, channel);
    end
    return
end

mpi0 = 0.13498;  mpi_pm = 0.13957;  mK = 0.493677;  mK0 = 0.497611;   % meson thresholds
m_eta = 0.547862;  m_etap = 0.95778;                                 % (kinks in 'all'/'sm' xsec)

% ---- Integration range (global, covers every T) -------------------------
s_lo  = 4*mx^2;
s_max = max([100*max(Tgrid)^2, 100*mphi^2, 100*mx^2, 1e2]);
L_lo  = log(s_lo);
L_hi  = log(s_max);

% ---- Breakpoints at cross-section kinks (meson + phiphi thresholds) ------
bp = [4*mpi_pm^2, 4*mpi0^2, 4*m_eta^2, 4*mK^2, 4*mK0^2, (mpi0+m_eta)^2, 4*mphi^2, 4*m_etap^2];   % = svavg.m extra_wp
bp = bp(bp > s_lo & bp < s_max);
L_bp = log(bp);

% ---- Breit-Wigner pole at s = mphi^2 -------------------------------------
L_mphi2  = log(mphi^2);
rel_w    = (mphi > 0) * (Gphi / max(mphi, eps));
pole_in  = (mphi^2 > s_lo) && (mphi^2 < s_max) && (mphi^2 > 4*mx^2);
use_nwa  = pole_in && (Gphi > 0) && (rel_w < 1e-6);

logdelta = 1e-3;        % notch half-width in log(s) for the NWA case
if use_nwa
    % Notch [L_mphi2-logdelta, L_mphi2+logdelta] out of the continuum; the
    % pole is added analytically below.
    notch_lo = L_mphi2 - logdelta;
    notch_hi = L_mphi2 + logdelta;
    L_extra  = [notch_lo, notch_hi];
elseif pole_in
    % Resolvable resonance: straddle the peak with nested panels at
    % +/- {1,5,30,100} FWHMs so Gauss-Legendre samples the Lorentzian top.
    fwhm = mphi * Gphi;
    cand = mphi^2 + [-100 -30 -5 -1 1 5 30 100]*fwhm;
    cand = cand(cand > s_lo & cand < s_max);
    L_extra = log(cand);
else
    L_extra = [];
end

% ---- Assemble sorted panel edges (graded near threshold) -----------------
edges = loggrid_edges(L_lo, L_hi, [L_bp, L_extra]);

% ---- Composite Gauss-Legendre nodes/weights in log(s) --------------------
npts = 24;
[xb, wb] = gl_nodes_weights(npts);
ls = []; wl = [];
for j = 1:numel(edges)-1
    a = edges(j); b = edges(j+1);
    if b <= a, continue; end
    if use_nwa && abs(a-notch_lo) < 1e-12 && abs(b-notch_hi) < 1e-12
        continue;          % skip the notched pole window
    end
    mid = (a+b)/2; half = (b-a)/2;
    ls = [ls; mid + half*xb];
    wl = [wl; half*wb];
end

s_nodes = exp(ls);                  % column N x 1
sq      = sqrt(s_nodes);
sigvals = sig(s_nodes);             % sigma evaluated ONCE
sigvals = real(sigvals(:));

% base column (s_node factor = dlog(s) Jacobian; matches svavg's
% logintegrand = integrand(s).*s integrated in log s)
base = wl .* s_nodes .* sigvals .* (s_nodes - 4*mx^2) .* sq;   % N x 1

% ---- Reduce against the Bessel kernel for every T ------------------------
% integrand kernel = besselk(1, sqrt(s)/T, 1) .* exp((2 mx - sqrt(s))/T)
Trow = Tgrid;                                   % 1 x M
K1   = besselk(1, sq ./ Trow, 1);               % N x M (scaled)
expf = exp((2*mx - sq) ./ Trow);                % N x M
I    = sum(base .* K1 .* expf, 1);              % 1 x M

% ---- Analytic NWA pole contribution (if notched) -------------------------
if use_nwa
    s_pole      = mphi^2;
    sig_pole    = real(sig(s_pole));
    prop2_pole  = (mphi*Gphi)^2;
    sqp         = sqrt(s_pole);
    smooth_pole = (s_pole - 4*mx^2) * sqp ...
                  .* besselk(1, sqp ./ Trow, 1) .* exp((2*mx - sqp) ./ Trow);
    I = I + sig_pole * prop2_pole * pi/(mphi*Gphi) .* smooth_pole;
end

% ---- Gondolo-Gelmini normalization ---------------------------------------
K2scaled = besselk(2, mx ./ Trow, 1);
sv = I ./ (8*mx^4 .* Trow .* K2scaled.^2);
sv = max(0, sv);
end
