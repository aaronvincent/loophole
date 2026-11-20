function [G, out] = primakoff_phi_rate(T, mphi, g2)
% primakoff_phi_rate  Per-phi rate of phi <-> photon conversion on the e+- plasma.
%
%   G = primakoff_phi_rate(T, mphi, g2)     [GeV], vectorised in T [GeV]
%   [G, out] = ...   out.compton, out.annihilation: the two pieces separately
%
% g2 = g_phigamgam^2 [GeV^-2], the phi F F (or phi F F~) coupling with
% Gamma(phi -> gamma gamma) = g_phigamgam^2 mphi^3 / (64 pi); callers take it
% from phi_partial_widths as 64 pi Gamma_gg / mphi^3.
%
% Two processes, both through the phi gamma gamma vertex:
%   e+- gamma -> e+- phi   (Primakoff, t-channel photon, Debye-screened)
%   e+ e-     -> gamma phi (s-channel photon)
% Massless electrons in the kinematics; X = s - mphi^2, x = -t in [0, X]:
%   sum|M|^2 (e gamma -> e phi) = e^2 g2 x (s^2 + (X - x)^2) / (x + mD^2)^2
% (unscreened: e^2 g2 (s^2 + u^2)/(-t); the same for a scalar and a
% pseudoscalar vertex).  Spin/polarisation averaged, dsigma/dt = |M|^2/(16 pi s^2):
%   sigma_C = alpha g2 I / (16 s^2),
%   I = [ y^2/2 - (2B + D) y + (C + 2BD) ln y + D C / y ]_{y=D}^{D+X},
%   D = mD^2,  B = X + D,  C = s^2 + B^2,
% -> (alpha g2 / 8) [ln(s/mD^2) - 7/4] at s >> mD^2;
%   sigma_A = alpha g2 (1 - mphi^2/s)^3 / 24.
% Rate density (Gondolo-Gelmini, massless initial states, lambda = s^2):
%   R = T/(32 pi^4) ∫_{mphi^2} ds [8 f_e sigma_C + 4 f_e^2 sigma_A] s^(3/2) K1(√s/T),
% with g_e g_gamma = (e- + e+, 2 spins each) x 2 = 8 for Compton and
% 2 x 2 = 4 for annihilation.  G = R / n_phi^eq (MB, one dof), the phi
% absorption rate by detailed balance, so G/H compares directly with
% Gamma_gg K1/K2.
%
% The electron mass enters only through
%   f_e = n_e(m_e, T)/n_e(0, T) = z^2 K2(z)/2,  z = m_e/T,
% which multiplies the number of scatterers (once for Compton, twice for
% annihilation) and the Debye mass mD^2 = (4 pi alpha/3) T^2 f_e.  It is a
% good approximation for T >~ m_e/3 and only O(1) below it, where both pieces
% are e^{-m_e/T}-suppressed.  T < m_e/60 returns 0.  T << m_mu is assumed
% (no mu+- or quark targets).
%
% Exactly ∝ g2, so the shape is computed once per (mphi, T) at g2 = 1 and
% cached (persistent, keyed on mphi and the T values).
% See docs/phi_thermalisation_check.md.

me = 0.51099895e-3;
G  = zeros(size(T));
out = struct('compton', G, 'annihilation', G);
if ~(g2 > 0) || isempty(T), return; end

persistent cache
if isempty(cache), cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
key = sprintf('%.17g|', mphi, T(:));
if isKey(cache, key)
    U = cache(key);
else
    U = zeros(2, numel(T));
    for i = 1:numel(T)
        [U(1,i), U(2,i)] = unit_rate(T(i), mphi, me);
    end
    if cache.Count > 2000, remove(cache, keys(cache)); end
    cache(key) = U;
end
out.compton      = g2 * reshape(U(1,:), size(T));
out.annihilation = g2 * reshape(U(2,:), size(T));
G = out.compton + out.annihilation;
end

function [Gc, Ga] = unit_rate(T, mphi, me)
Gc = 0;  Ga = 0;
if T < me/60, return; end
alpha = 1/137.036;
ze = me/T;
if ze < 1e-3, fe = 1; else, fe = ze^2*besselk(2, ze)/2; end
D  = (4*pi*alpha/3) * T^2 * fe;
xth = mphi;                                  % sqrt(s) threshold
% s = x^2: ds s^(3/2) = 2 x^4 dx;  K1(x/T) = besselk(1,x/T,1) e^{-xth/T} e^{-(x-xth)/T}
w  = @(x) 2*x.^4 .* besselk(1, x/T, 1) .* exp(-(x - xth)/T);
wp = xth + T*[0.01 0.1 1 5];
Ic = integral(@(x) w(x) .* primakoff_phi_xsec(x.^2, mphi, D), xth, xth + 80*T, ...
              'RelTol', 1e-8, 'AbsTol', 0, 'Waypoints', wp);
Ia = integral(@(x) w(x) .* sig_ann(x.^2, mphi, D), xth, xth + 80*T, ...
              'RelTol', 1e-8, 'AbsTol', 0, 'Waypoints', wp);
z = mphi/T;
if z < 1e-4
    n_scaled = T^3/pi^2;                         % massless MB, e^{-z} ~ 1
else
    n_scaled = mphi^2*T*besselk(2, z, 1)/(2*pi^2);   % n e^{+mphi/T}
end
% R = T/(32 pi^4) e^{-xth/T} I, n = n_scaled e^{-mphi/T}, xth = mphi.
Gc = 8*fe    * T/(32*pi^4) * Ic / n_scaled;
Ga = 4*fe^2  * T/(32*pi^4) * Ia / n_scaled;
end

function sa = sig_ann(s, mphi, D)
[~, sa] = primakoff_phi_xsec(s, mphi, D);
end
