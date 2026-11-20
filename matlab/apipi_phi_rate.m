function G = apipi_phi_rate(T, mA, g)
% apipi_phi_rate  Per-A' absorption rate from A' pi -> pi pi in the bath (axialVector).
%
%   G = apipi_phi_rate(T, mA, g)     [GeV], vectorised in T [GeV]
%
% The chi chi -> 3pi current crossed to 2 -> 2.  On
% shell the pi0-pole piece is proportional to k.eps = 0, so only the contact
% term survives:
%     A'(k) pi0(p)  -> pi+ pi- :  V = (2Da/3F)(p+ + p- + 2p) = (2Da/3F)(k + 3p), so
%                                 M = (2Da/F) eps.p and, summed over pol.,
%                                 |M|^2 = 4 Da^2/F^2 [(k.p)^2/mA^2 - m0^2]   (isotropic)
%     A'(k) pi+-(p) -> pi+- pi0(p0) :  same with p -> p0 (outgoing pi0), so
%                                 <(k.p0)^2> = (E_k E_0)^2 + (p_i p_f)^2/3
% Da = m_{pi0 A'}/F from axial_mixing.m, the A'-pi0 normalisation of the 3pi current; F = 92.4 MeV.  sigma = (1/3) sum|M|^2 p_f/(16 pi s p_i),
% no identical particles in either final state.
%
% G = R / n_A'^eq with the Gondolo-Gelmini unequal-mass rate density
%     R = g_A g_pi T/(32 pi^4) ∫ ds sigma(s) lambda(s, mA^2, m_pi^2) K1(√s/T)/√s,
% g_A = 3, g_pi = 1 per charge state (pi0, pi+, pi-), and
% n_A'^eq = 3 mA^2 T K2(mA/T)/(2 pi^2).  The pions are in equilibrium, so by
% detailed balance freeze_coupled uses it like the other Gamma_sc terms,
% on (Y_A' - Y_A'eq).
%
% Scales exactly as m_{pi0 A'}^2, so the shape is cached per (mA, T) at Da = 1.
% T < m_pi/60 returns 0.  Kaon analogues (A' K -> K pi, ...) are dropped:
% they pay e^{-m_K/T} at the temperatures where this matters.

m0 = 0.13498;  mc = 0.13957;  F = 0.0924;
G  = zeros(size(T));
if isempty(T) || mA <= 0, return; end
mjA = axial_mixing(g);
c2  = (mjA(1)/F)^2;
if c2 == 0, return; end

persistent cache
if isempty(cache), cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
key = sprintf('%.17g|', mA, T(:));
if isKey(cache, key)
    G = c2 * cache(key);
    return
end

lam = @(s, a, b) (s - (sqrt(a) + sqrt(b)).^2) .* (s - (sqrt(a) - sqrt(b)).^2);
% channel 1: A' pi0 -> pi+ pi-
th1 = max(mA + m0, 2*mc);
sig1 = @(s) (4/F^2) * (lam(s, mA^2, m0^2)/(4*mA^2)) / 3 ...
            .* sqrt(max(lam(s, mc^2, mc^2), 0)) ./ (16*pi*s .* sqrt(lam(s, mA^2, m0^2)));
% channel 2: A' pi+- -> pi+- pi0 (each charge)
th2 = max(mA + mc, mc + m0);
sig2 = @(s) sig2_fun(s, mA, mc, m0, F, lam);

Gu = zeros(size(T));
for i = 1:numel(T)
    Ti = T(i);
    if Ti < mc/60, continue; end
    I1 = gg_int(sig1, th1, mA, m0, Ti, lam);
    I2 = gg_int(sig2, th2, mA, mc, Ti, lam);
    xth = min(th1, th2);
    % R e^{+xth/T}, with g_A = 3 and the two charged pions
    R_scaled = 3*Ti/(32*pi^4) * (I1*exp(-(th1 - xth)/Ti) + 2*I2*exp(-(th2 - xth)/Ti));
    z = mA/Ti;
    n_scaled = 3*mA^2*Ti*besselk(2, z, 1)/(2*pi^2);    % n e^{+mA/T}
    Gu(i) = R_scaled*exp(-(xth - mA)/Ti)/n_scaled;
end
if cache.Count > 2000, remove(cache, keys(cache)); end
cache(key) = Gu;
G = c2 * Gu;
end

function I = gg_int(sig, xth, m1, m2, Ti, lam)
% ∫ ds sigma lambda(s,m1^2,m2^2) K1(√s/T)/√s, times e^{+xth/T}; sqrt(s) = x, ds = 2x dx
f = @(x) sig(x.^2) .* lam(x.^2, m1^2, m2^2) .* besselk(1, x/Ti, 1) ...
         .* exp(-(x - xth)/Ti) .* 2;
I = integral(f, xth*(1 + 1e-12), xth + 60*Ti, 'RelTol', 1e-7, 'AbsTol', 0, ...
             'Waypoints', xth + Ti*[0.01 0.1 1 5]);
end

function s2 = sig2_fun(s, mA, mc, m0, F, lam)
rs = sqrt(s);
Ek = (s + mA^2 - mc^2)./(2*rs);   E0 = (s + m0^2 - mc^2)./(2*rs);
pi2 = lam(s, mA^2, mc^2)./(4*s); pf2 = max(lam(s, mc^2, m0^2), 0)./(4*s);
M2 = (4/F^2) * (((Ek.*E0).^2 + pi2.*pf2/3)/mA^2 - m0^2);
s2 = M2/3 .* sqrt(pf2) ./ (16*pi*s .* sqrt(pi2));
end
