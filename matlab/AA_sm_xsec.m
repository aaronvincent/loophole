function [sig_gg, sig_mm] = AA_sm_xsec(s, mA, g)
% AA_sm_xsec  A'A' -> SM cross sections for the axialVector mediator [GeV^-2].
%
%   [sig_gg, sig_mm] = AA_sm_xsec(s, mA, g)    s [GeV^2], mA = m_A', g = [a_u a_d a_s]
%
% Both come from the one O(p^2) chiPT vertex quadratic in the source
% (docs/notes/AA_to_gammagamma.tex, docs/notes/threepi_annihilation.tex sec. pipiA):
%   L_aa = -sum_{i<j} (a_i - a_j)^2 A'.A' M_ij M_ji,
% i.e. -(a_u-a_d)^2 A'.A' pi+pi- - (a_u-a_s)^2 A'.A' K+K- - (a_d-a_s)^2 A'.A' K0 K0bar.
% Diagonal mesons commute with the diagonal source: no pi0/eta/eta' term.
%
%   sig_gg  A'A' -> gamma gamma, one charged-meson loop (pi+-, K+-).  The contact
%           is a scalar source sigma = A'.A' for |M+|^2, so the loop is h -> gamma gamma
%           with m_h^2 -> s:  L_eff = c(s) (A'.A') F_{mu nu} F^{mu nu},
%             c(s) = alpha/(16 pi) sum_j Da_j^2 A0(s/4m_j^2)/m_j^2,   A0(0) = 1/3,
%           sigma = |c|^2 s P(s) / (9 pi beta_A).
%   sig_mm  A'A' -> pi+pi-, K+K-, K0K0bar at tree level (contact):
%           sigma = Da^4 P(s) beta_M / (36 pi s beta_A).
% P(s) = sum_pol |eps1.eps2|^2 = 2 + (s - 2mA^2)^2/(4 mA^4): the second term is
% the two longitudinal A' (Goldstone-equivalent to the pi0 admixture, angle
% F Da/mA), so both grow as s^2/mA^4 for a light A'.  1/9 averages over the
% initial polarizations; 1/2 for identical photons.
%
% One loop carries corr4^(V-I) = corr4^0: LO Da (no L4 rescaling), unlike the
% tree-level 3pi / pi pi A' currents (docs/notes/threepi_annihilation.tex).

alpha = 1/137.036;
mpi = 0.13957;  mK = 0.493677;  mK0 = 0.497611;
s  = reshape(s, 1, []);
sig_gg = zeros(size(s));  sig_mm = zeros(size(s));
ok = s > 4*mA^2;
if ~any(ok) || mA <= 0, return, end
sp = s(ok);
bA = sqrt(1 - 4*mA^2 ./ sp);
P  = 2 + (sp - 2*mA^2).^2 / (4*mA^4);

Dpi = g(1) - g(2);  DK = g(1) - g(3);  DK0 = g(2) - g(3);

% gamma gamma: charged-meson loop
c = alpha/(16*pi) * (Dpi^2 * A0(sp/(4*mpi^2)) / mpi^2 + DK^2 * A0(sp/(4*mK^2)) / mK^2);
sig_gg(ok) = abs(c).^2 .* sp .* P ./ (9*pi*bA);

% meson pairs, tree level
mm = zeros(size(sp));
for ch = [Dpi, mpi; DK, mK; DK0, mK0].'
    D = ch(1);  m = ch(2);
    op = sp > 4*m^2;
    if D ~= 0 && any(op)
        bM = sqrt(1 - 4*m^2 ./ sp(op));
        mm(op) = mm(op) + D^4 * P(op) .* bM ./ (36*pi*sp(op) .* bA(op));
    end
end
sig_mm(ok) = mm;
end

function a = A0(tau)
% Charged-scalar loop function, A0 = (f(tau) - tau)/tau^2, A0(0) = 1/3
% (Djouadi's A_0^H up to sign); complex above threshold tau > 1.
a = zeros(size(tau));
sm = tau < 1e-3;  lo = ~sm & tau <= 1;  hi = tau > 1;
a(sm) = 1/3 + 8*tau(sm)/45 + 4*tau(sm).^2/35;
f = asin(sqrt(tau(lo))).^2;
a(lo) = (f - tau(lo)) ./ tau(lo).^2;
if any(hi)
    t = tau(hi);  r = sqrt(1 - 1./t);
    f = -0.25*(log((1 + r)./(1 - r)) - 1i*pi).^2;
    a(hi) = (f - t) ./ t.^2;
end
end
