function W = threepi_weight(s, g)
% threepi_weight  Three-pion final states of a pseudoscalar phi-SM vertex.
%
%   W = threepi_weight(s, g)     s vector [GeV^2], g = [g_u g_d g_s]
%
% Returns W(s) = sum_ch S int dE1 dE2 |A|^2 over pi+pi-pi0 (S = 1) and
% pi0pi0pi0 (S = 1/6), where A is the phi* -> 3pi amplitude with the mediator
% propagator stripped off.  The consumers are
%   sigma(chi chi -> 3pi) = gx^2 W / (128 pi^3 beta |s - mphi^2 + i mphi G|^2)
%       (PP; SP has an extra (s - 4 mx^2)/s)              threepi_sigma.m
%   Gamma(phi -> 3pi)     = W(mphi^2) / (64 pi^3 mphi)    phi_partial_widths.m
%
% Amplitude: LO chiPT.  In the (phi, mesons) basis,
%   A(+-0) = (B0 D/F) X / (s - mpi^2)  +  sum_j a_j b_j(X) / (s - mj^2)
%   A(000) = (B0 D/F) s / (s - mpi^2)  +  sum_j a_j c_j    / (s - mj^2)
% with D = g_u - g_d, X = m_{+-}^2 - mpi^2 and a_j the LO phi-meson mass mixing.
% The first terms already hold the direct phi 3pi vertex, the pi0 pole of the
% contact and kinetic 4-pion vertices, and the phi-pole share of the mixing.
% Use the LO mixing here, not sppp_mixing's NLO bracket: that is 2.58x larger
% and would spoil the cancellation against the LO direct vertex.  j = eta
% generically, eta' for g ~ [1 1 1], the same pole set as sppp_couplings.m, so
% the integrators' pole lists stay right.  b_eta carries the LO energy
% dependence of eta -> pi+pi-pi0; it sums away in the 3pi0 channel.
%
% The channels are off when threepi_enabled() is false.

W = zeros(size(s));
if ~threepi_enabled(), return, end

mpi0 = 0.13498;  mpi_pm = 0.13957;
m_eta = 0.547862;  m_etap = 0.95778;
G_eta = 1.31e-6;   G_etap = 1.88e-4;
B0   = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV
f    = 0.0924;
dm   = 0.00469 - 0.0022;               % m_d - m_u

D = g(1) - g(2);
if g(1) == g(2) && g(2) == g(3)        % degenerate: eta' singlet
    mj = m_etap;  Gj = G_etap;
    aj = B0*f*sqrt(2/3)*sum(g);
    cj = -sqrt(6)/3*B0*dm/f^2;         % 3! x the coefficient in L
    bj = -sqrt(6)/9*B0*dm/f^2;
else
    mj = m_eta;   Gj = G_eta;
    aj = B0*f*(g(1) + g(2) - 2*g(3))/sqrt(3);
    cj = -B0*dm/(sqrt(3)*f^2);         % 3! x the coefficient in L
    bj = -B0*dm/(3*sqrt(3)*f^2);
end
d2 = mj^2 - mpi0^2;
pole = @(x) aj ./ (x - mj^2 + 1i*mj*Gj);   % meson pole, common to both channels

% pi0pi0pi0.  Constant over the Dalitz plot.
mask = s > 9*mpi0^2;
if any(mask)
    sp = s(mask);  sp = sp(:);
    A  = B0*D/f .* sp ./ (sp - mpi0^2) + pole(sp)*cj;
    W(mask) = abs(A).^2 .* dalitz_moments(sp, mpi0, mpi0) / 6;
end

% pi+pi-pi0.  A = P + Q*X, so only moments of X are needed.
mask = s > (2*mpi_pm + mpi0)^2;
if any(mask)
    sp = s(mask);  sp = sp(:);
    P  = pole(sp)*bj .* (1 - sp/d2);
    Q  = B0*D/f ./ (sp - mpi0^2) + pole(sp)*bj*3/d2;
    [I0, I1, I2] = dalitz_moments(sp, mpi0, mpi_pm);
    W(mask) = W(mask) + reshape(abs(P).^2.*I0 + 2*real(P.*conj(Q)).*I1 + abs(Q).^2.*I2, size(W(mask)));
end
end
