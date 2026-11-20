function G = pigamma_phi_rate(T, mphi, g)
% pigamma_phi_rate  Per-phi production rate from pi gamma -> pi phi in the bath.
%
%   G = pigamma_phi_rate(T, mphi, g)     [GeV], vectorised in T [GeV]
%
% G = R / n_phi^eq, with R the rate density of pi+gamma -> pi+phi plus
% pi-gamma -> pi-phi, Maxwell-Boltzmann statistics (Gondolo-Gelmini form for
% unequal masses):
%     R = g_pi g_gamma T/(32 pi^4) ∫ ds sigma(s) lambda(s, m_pi^2, 0) K1(√s/T)/√s,
% g_pi = 2 (pi+ and pi-, each counted once), g_gamma = 2 (sigma is already
% polarisation-averaged), lambda(s,m^2,0) = (s - m^2)^2.
% n_phi^eq = mphi^2 T K2(mphi/T)/(2 pi^2)  (real scalar, MB; -> T^3/pi^2).
% Equal to the phi -> (pi gamma) absorption rate by detailed balance, so
% G/H compares directly with the inverse-decay rate Gamma_gg K1/K2.
%
% Scales exactly as (g_u + g_d)^2, so the shape is computed once per
% (mphi, T) at g_u + g_d = 1 and cached (persistent, keyed on mphi and the T
% values).  T < m_pi/60 returns 0: the rate is down by e^{-60} there and
% Gamma/H < 1e-10 for any coupling the code reaches.
%
% S⊗S / P⊗S only; mphi < 2 m_pi.  See pigamma_phi_xsec.m and
% docs/phi_thermalisation_check.md.

m  = 0.13957;
G  = zeros(size(T));
c2 = (g(1) + g(2))^2;
if c2 == 0 || isempty(T), return; end

persistent cache
if isempty(cache), cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
key = sprintf('%.17g|', mphi, T(:));
if isKey(cache, key)
    G = c2 * cache(key);
    return
end

Gu = zeros(size(T));
xth = m + mphi;                               % sqrt(s) at threshold
for i = 1:numel(T)
    Ti = T(i);
    if Ti < m/60, continue; end
    % sqrt(s) = x; ds = 2x dx.  Exponentials pulled out:
    % K1(x/T) = besselk(1,x/T,1) e^{-x/T},  e^{-x/T} = e^{-xth/T} e^{-(x-xth)/T}.
    f = @(x) pigamma_phi_xsec(x.^2, mphi, [1 0 0]) .* (x.^2 - m^2).^2 ...
             .* besselk(1, x/Ti, 1) .* exp(-(x - xth)/Ti) .* 2;   % (2x dx)/x = 2 dx
    I = integral(f, xth, xth + 60*Ti, 'RelTol', 1e-7, 'AbsTol', 0, ...
                 'Waypoints', xth + Ti*[0.01 0.1 1 5]);
    R_scaled = 4*Ti/(32*pi^4)*I;              % R e^{+xth/T}
    z = mphi/Ti;
    if z < 1e-4
        n_scaled = Ti^3/pi^2;                 % massless limit, e^{-z} ~ 1
        Gu(i) = R_scaled*exp(-xth/Ti)/n_scaled;
    else
        n_scaled = mphi^2*Ti*besselk(2, z, 1)/(2*pi^2);   % n e^{+mphi/T}
        Gu(i) = R_scaled*exp(-(xth - mphi)/Ti)/n_scaled;  % = e^{-m_pi/T} ...
    end
end
if cache.Count > 2000, remove(cache, keys(cache)); end
cache(key) = Gu;
G = c2 * Gu;
end
