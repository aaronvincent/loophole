function lsv = sv_AAsm_grid(mA, g, Tgrid)
% sv_AAsm_grid  log <sigma v> of A'A' -> SM (gamma gamma + meson pairs), axialVector.
%
%   lsv = sv_AAsm_grid(mA, g, Tgrid)     natural log, -Inf where zero
%
% Gondolo-Gelmini (Maxwell-Boltzmann, as freeze_coupled's Y_eq):
%   <sv> = int sigma(s) (s - 4mA^2) sqrt(s) K1(sqrt(s)/T) ds / (8 mA^4 T K2(mA/T)^2),
% sigma = AA_sm_xsec (gamma gamma + meson pairs).  freeze_coupled uses it as
%   dY_A'/dtau = -(s <sv>/H) (Y_A'^2 - Y_A'eq^2)        (A' self-conjugate, no 1/2).
% Quadrature as sv_darkpair_grid (sqrt(s) = b + u T, u on 400 log nodes in
% [1e-10, 80], trapezoid), piecewise between the thresholds 2mA, 2m_pi+-, 2m_K+-,
% 2m_K0, so each cusp sits at a panel's first node.

Tgrid = Tgrid(:).';
nT  = numel(Tgrid);
lsv = -Inf(1, nT);
if mA <= 0, return, end
u  = [1e-10, logspace(-9, log10(80), 399)].';
mth = [2*0.13957, 2*0.493677, 2*0.497611];
b  = [2*mA, sort(mth(mth > 2*mA))];
I  = zeros(1, nT);                                % int ... ds, times e^{+2mA/T}
for k = 1:numel(b)
    if k < numel(b), b_hi = b(k+1); else, b_hi = Inf; end
    RS = b(k) + u .* Tgrid;  S = RS.^2;
    [sgg, smm] = AA_sm_xsec(S(:).', mA, g);
    sig = reshape(sgg + smm, size(S));
    sig(RS >= b_hi) = 0;                          % the next panel takes over
    F = sig .* (S - 4*mA^2) .* RS .* besselk(1, RS ./ Tgrid, 1) .* exp(-u) .* (2*RS .* Tgrid);
    I = I + trapz(u, F, 1) .* exp(-(b(k) - 2*mA) ./ Tgrid);
end
ok = I > 0 & isfinite(I);
T  = Tgrid(ok);
% 8 mA^4 T K2^2 = 8 T^5 (z^2 K2(z))^2;  besselk(1,.,1) carried e^{+sqrt(s)/T}
lsv(ok) = log(I(ok)) - 2*mA ./ T - log(8) - 5*log(T) - 2*log_z2K2(mA ./ T);
end

function l = log_z2K2(z)
l = zeros(size(z));
hi = z > 30;  lo = z < 1e-3;  mid = ~hi & ~lo;
l(hi)  = 2*log(z(hi)) - z(hi) + 0.5*log(pi./(2*z(hi))) + log1p(15./(8*z(hi)));
l(lo)  = log(2);
l(mid) = log(z(mid).^2 .* besselk(2, z(mid)));
end
