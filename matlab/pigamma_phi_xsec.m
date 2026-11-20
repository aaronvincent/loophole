function sig = pigamma_phi_xsec(s, mphi, g)
% pigamma_phi_xsec  sigma(pi+ gamma -> pi+ phi) for a SCALAR phi-quark coupling.
%
%   sig = pigamma_phi_xsec(s, mphi, g)     [GeV^-2], vectorised in s [GeV^2]
%
% Leading-order chiPT with M_q -> M_q + G phi: L ⊃ -mu phi pi+ pi-,
% mu = B0 (g_u + g_d) (same vertex as meson_pairs.m).  No phi phi pi pi
% contact and no derivative coupling, so no seagull: only the s- and u-channel
% pion poles, which combine into the soft-photon form
%     M = e mu [ p.eps/(p.k) - p'.eps/(p'.k) ].
% With a = s - m^2, b = u - m^2 (< 0):
%     sum_pol |M|^2 = -4 e^2 mu^2 [ m^2 (a+b)^2/(a b)^2 + (a+b)/(a b) - mphi^2/(a b) ]
%     sigma = -(e^2 mu^2 / (8 pi a^2)) * I,
%     I = (b+ - b-)(m^2/a^2 + 1/a) + m^2 (1/b- - 1/b+)
%         + ((2m^2 - mphi^2)/a + 1) ln(b+/b-),
%     b± = mphi^2 - [ (s+m^2)(s+mphi^2-m^2) ∓ a lambda^(1/2)(s,m^2,mphi^2) ] / (2s).
% Averaged over the two photon polarisations.  The pi- gamma -> pi- phi cross
% section is identical; the pi0 does not couple to the photon at this order.
%
% Applies to S⊗S / P⊗S (heavyScalar*, lightScalar*, *ScalarPS) only: the chi
% vertex never enters.  Valid for mphi < 2 m_pi (above it pi pi -> phi inverse
% decay leads).  The final pi-phi state is p-wave (E1).
%
% Numerics: the closed form cancels catastrophically near threshold for light
% phi (relative error ~1 at s/s_th - 1 ~ 1e-5 for mphi = 1e-9 GeV).  Points
% with s/s_th - 1 < THR_SWITCH are done by direct quadrature of the grouped
% |M|^2 over b, which has no cancellation.  Formulae and checks:
% docs/phi_thermalisation_check.md.

THR_SWITCH = 1e-2;
alphaEM = 1/137.036;
m   = 0.13957;                                % m_pi+ [GeV]
B0  = m^2/(0.0022 + 0.00469);                 % ~2.827 GeV (as meson_pairs.m)
mu  = B0*(g(1) + g(2));
pref = 4*pi*alphaEM*mu^2/(8*pi);              % e^2 mu^2 / (8 pi)

sth = (m + mphi)^2;
sig = zeros(size(s));
ok  = s > sth;
rel = (s - sth)/sth;

closed = ok & rel >= THR_SWITCH;
if any(closed(:))
    sc = s(closed);
    [a, bp, bm] = bounds(sc, m, mphi);
    I = (bp - bm).*(m^2./a.^2 + 1./a) + m^2*(1./bm - 1./bp) ...
        + ((2*m^2 - mphi^2)./a + 1).*log(bp./bm);
    sig(closed) = -pref*I./a.^2;
end

near = find(ok & rel < THR_SWITCH);
% For mphi ~ 1e-9 the b-interval is ~1e-10 wide and integral() reports
% "minimum step size" / "max interval count" at RelTol 1e-10; the value is
% unaffected (continuity across the switch holds to 1e-9 in test_pigamma_phi,
% and the rate is mphi-independent there), so silence both locally.
wst = warning;  warning('off', 'MATLAB:integral:MinStepSize');
warning('off', 'MATLAB:integral:MaxIntervalCountReached');
cleanup = onCleanup(@() warning(wst));
for k = near(:)'
    [a, bp, bm] = bounds(s(k), m, mphi);
    if ~(bm < bp), continue; end
    f = @(b) -( m^2*((a + b)./(a*b)).^2 + (a + b)./(a*b) - mphi^2./(a*b) );
    sig(k) = pref*integral(f, bm, bp, 'RelTol', 1e-10, 'AbsTol', 0)/a^2;
end
end

function [a, bp, bm] = bounds(s, m, mphi)
a = s - m^2;
L = sqrt(max((s - (m + mphi)^2).*(s - (m - mphi)^2), 0));
c = (s + m^2).*(s + mphi^2 - m^2);
bp = mphi^2 - (c - a.*L)./(2*s);
bm = mphi^2 - (c + a.*L)./(2*s);
end
