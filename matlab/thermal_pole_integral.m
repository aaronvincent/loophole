function I = thermal_pole_integral(logintegrand, sig, kernel_at, s_lo, s_max, pole_s, pole_w, opts)
% thermal_pole_integral  Continuum + narrow-resonance integral in log(s).
%
%   I = thermal_pole_integral(logintegrand, sig, kernel_at, s_lo, s_max, ...
%                             pole_s, pole_w, opts)
%
% Shared resonance-aware quadrature for the thermal s-integrals in svavg.m
% (chi chibar -> gamma gamma / all) and ggxxrate.m (gamma gamma -> chi chi).
% Both integrate sigma(s) times a smooth Bessel kernel over s; the cross
% sections have one or more s-channel Breit-Wigner poles at s = m_j^2 that the
% adaptive integrator silently misses when the width is narrow.
%
% Inputs:
%   logintegrand(logs) : full integrand expressed in log s (= integrand(s).*s),
%                        used for the smooth continuum away from notched poles.
%   sig(s)             : the BW-regulated cross section sigma(s); finite at the
%                        pole, used for the narrow-width residue sig(s_p)*prop2.
%   kernel_at(s_p)     : the smooth (non-sigma) integrand factors at s = s_p,
%                        i.e. integrand(s_p)/sig(s_p).
%   s_lo, s_max        : integration limits in s.
%   pole_s, pole_w     : equal-length arrays of pole positions s_p = m_p^2 and
%                        widths w_p = m_p*Gamma_p.
%   opts (optional)    : struct with fields (all optional)
%                        .RelTol (1e-6) .AbsTol (1e-20) .logdelta (1e-3)
%                        .narrow_tol (1e-5, on Gamma_p/m_p = w_p/s_p)
%                        .extra_wp (extra continuum waypoints in s, e.g. the
%                                   meson-pair / A'A' thresholds that kink the
%                                   'all' cross section).
%
% Per pole (in range, above threshold): narrow (w_p/s_p < narrow_tol) -> notch
% a +/-logdelta window in log s and add the narrow-width contribution
%   sig(s_p) * w_p^2 * kernel_at(s_p) * pi / w_p   (= numerator * pi/w_p);
% resolvable -> leave it in the continuum and straddle it with a geometric
% waypoint stack so the adaptive rule samples the Lorentzian top.

if nargin < 8 || isempty(opts), opts = struct(); end
RelTol     = getdef(opts, 'RelTol',     1e-6);
AbsTol     = getdef(opts, 'AbsTol',     1e-20);
logdelta   = getdef(opts, 'logdelta',   1e-3);
narrow_tol = getdef(opts, 'narrow_tol', 1e-5);
extra_wp   = getdef(opts, 'extra_wp',   []);

L_lo = log(s_lo);
L_hi = log(s_max);
I = 0;
if L_hi <= L_lo, return; end

% --- classify poles -------------------------------------------------------
notch   = zeros(0,2);    % [lo hi] log-s windows around narrow poles
wp_log  = log(extra_wp(extra_wp > s_lo & extra_wp < s_max));   % continuum waypoints
nwa_idx = [];            % narrow poles to add analytically
for p = 1:numel(pole_s)
    sp = pole_s(p);  wp = pole_w(p);
    Lp = log(sp);
    if ~(sp > s_lo && Lp > L_lo && Lp < L_hi), continue; end
    if wp > 0 && (wp/sp) < narrow_tol
        notch   = [notch; Lp - logdelta, Lp + logdelta];
        nwa_idx = [nwa_idx, p];
    else
        % resolvable resonance: +/- {1,5,50} FWHMs straddle the peak
        fw   = max(wp, 1e-12*sp);
        cand = sp + [-50 -5 -1 0 1 5 50]*fw;
        Lc   = log(cand(cand > s_lo));
        wp_log = [wp_log, Lc(Lc > L_lo & Lc < L_hi)];
    end
end

% --- continuum on [L_lo, L_hi] minus the notched windows ------------------
bounds = [L_lo; L_hi];
for n = 1:size(notch,1)
    bounds = [bounds; max(L_lo, notch(n,1)); min(L_hi, notch(n,2))];
end
bounds = unique(sort(bounds));
for ii = 1:numel(bounds)-1
    a = bounds(ii);  b = bounds(ii+1);
    if b <= a, continue; end
    mid = 0.5*(a + b);
    in_notch = false;
    for n = 1:size(notch,1)
        if mid > notch(n,1) && mid < notch(n,2), in_notch = true; break; end
    end
    if in_notch, continue; end
    wseg = wp_log(wp_log > a & wp_log < b);
    if isempty(wseg)
        I = I + integral(logintegrand, a, b, 'RelTol', RelTol, 'AbsTol', AbsTol);
    else
        I = I + integral(logintegrand, a, b, 'RelTol', RelTol, 'AbsTol', AbsTol, 'Waypoints', wseg);
    end
end

% --- narrow-width pole contributions --------------------------------------
% int ds 1/((s-s_p)^2 + w_p^2) = pi/w_p; sig(s_p)*w_p^2 strips the BW
% denominator, leaving the (finite) resonant numerator of pole p (other poles'
% off-resonant tails are suppressed by w_p^2 -> 0).
for q = nwa_idx
    sp = pole_s(q);  wp = pole_w(q);
    I = I + real(sig(sp)) * wp^2 * kernel_at(sp) * pi / wp;
end
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
