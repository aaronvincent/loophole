function [lsv_fwd, lsv_rev] = sv_darkpair_grid(mx, mphi, gx, g, Tgrid, interaction, g_phi)
% sv_darkpair_grid  log <sigma v> of chi chibar <-> phi phi in both directions.
%
%   [lsv_fwd, lsv_rev] = sv_darkpair_grid(mx, mphi, gx, g, Tgrid, interaction, g_phi)
%
% The dark mediator pair is xxggxsec(...,'dark'): phi phi (S⊗S, P⊗S), A'A'
% (axialVector) or phi'phi' (SP/PP, the mediator term of the jj sum) -- all
% t-channel, no s-channel pole.  With
%   I(T) = int sigma(s) (s - 4 mx^2) sqrt(s) K1(sqrt(s)/T) ds
% over s >= max(4 mx^2, 4 mphi^2):
%   <sv>_fwd = I / (8 mx^4 T K2(mx/T)^2)                       (Gondolo-Gelmini, = svavg)
%   <sv>_rev = (g_chi/g_phi)^2 I / (8 mphi^4 T K2(mphi/T)^2)
% so that <sv>_rev Y_phi_eq^2 = <sv>_fwd Y_chi_eq^2 exactly (Y counts chi + chibar,
% g_chi = 4).  freeze_coupled's phi phi -> chi chibar source is then
% lambda <sv>_rev Y_phi^2.  Both are returned as natural logs (-Inf where zero):
% each is normalised by its own initial state's K2, so the reverse rate stays
% O(1) when mphi > mx (phi phi -> chi chibar is exothermic) and the forward one
% when mx > mphi, at any m/T -- which a difference of 'all' and 'sm' averages,
% or the plain values, cannot give.
%
% Quadrature: sqrt(s) = sqrt(s_lo) + u T, u on 400 log-spaced nodes in
% [1e-10, 80] (resolves the threshold layer, width ~T, at any m/T), trapezoid.
% The kernel exp(-u) carries the Boltzmann factor relative to threshold; the
% threshold exponent and the K2 normalisations are combined in logs.

g_chi = 4;
Tgrid = Tgrid(:).';
nT = numel(Tgrid);
lsv_fwd = -Inf(1, nT);
lsv_rev = -Inf(1, nT);

rs_lo = 2*max(mx, mphi);
u  = [1e-10, logspace(-9, log10(80), 399)].';           % column
RS = rs_lo + u .* Tgrid;                                % nu x nT, sqrt(s)
S  = RS.^2;
sig = reshape(xxggxsec(mx, mphi, gx, g, S(:).', interaction, 'dark'), size(S));
if ~any(sig(:) > 0), return; end
% integrand of I e^{+rs_lo/T} in u:  sigma (s - 4mx^2) sqrt(s) K1s(sqrt(s)/T) e^{-u} * ds/du,
% ds/du = 2 sqrt(s) T
F  = sig .* (S - 4*mx^2) .* RS .* besselk(1, RS ./ Tgrid, 1) .* exp(-u) .* (2*RS .* Tgrid);
It = trapz(u, F, 1);                                    % 1 x nT
ok = It > 0 & isfinite(It);
lI = log(It(ok)) - rs_lo ./ Tgrid(ok);                  % log I
T  = Tgrid(ok);
% 8 m^4 T K2(m/T)^2 = 8 T^5 (z^2 K2(z))^2,  z = m/T
lsv_fwd(ok) = lI - log(8) - 5*log(T) - 2*log_z2K2(mx ./ T);
lsv_rev(ok) = lI - log(8) - 5*log(T) - 2*log_z2K2(mphi ./ T) + 2*log(g_chi/g_phi);
end

function l = log_z2K2(z)
% log(z^2 K2(z)), finite at any z (same guards as freeze_coupled's Y_eq).
l = zeros(size(z));
hi = z > 30;  lo = z < 1e-3;  mid = ~hi & ~lo;
l(hi)  = 2*log(z(hi)) - z(hi) + 0.5*log(pi./(2*z(hi))) + log1p(15./(8*z(hi)));
l(lo)  = log(2);
l(mid) = log(z(mid).^2 .* besselk(2, z(mid)));
end
