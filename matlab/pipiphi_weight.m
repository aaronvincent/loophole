function [W, W00] = pipiphi_weight(s, mphi, g)
% pipiphi_weight  pi pi phi final states of a pseudoscalar phi-SM vertex.
%
%   [W, W00] = pipiphi_weight(s, mphi, g)     s vector [GeV^2], g = [g_u g_d g_s]
%
% Returns W(s) = sum_ch S int dE1 dE2 |A|^2 over pi+pi-phi (S = 1) and pi0pi0phi
% (S = 1/2), W00 being the second.  Unlike threepi_weight, A keeps every
% propagator, so
%   sigma(chi chi -> pi pi phi) = gx^2 W / (128 pi^3 beta)
% with the extra (s - 4 mx^2)/s for SP (threepi_sigma.m).
%
% Amplitude: LO chiPT at O(g^2), mass basis, isospin limit in the vertices.  No
% cubic vertices exist, so A = sum_j U_{phi j} V_j / (s - mj^2) with V_j the
% contact term, including the outgoing phi's own meson content.  A depends on
% m_{pipi}^2 = x alone.  Poles are pi0, eta, phi, or eta', phi for g ~ [1 1 1]:
% the pole set of sppp_couplings.m.  Mixing is LO for the reason given in
% threepi_weight.m.  The phi pole sits below threshold, so it needs no width.
%
% The final state contains a phi, so using this cross section for production
% through detailed balance assumes phi is in the bath -- while the 'sm' rate
% drops chi chi -> phi phi and freeze_coupled's Y_phi has no scattering source.
% See docs/cross_sections.md, "Freeze-in consistency".
%
% The channels are off when pipiphi_enabled() is false.

W = zeros(size(s));  W00 = zeros(size(s));
if ~pipiphi_enabled(), return, end

mpi0 = 0.13498;  mpi_pm = 0.13957;
m_eta = 0.547862;  m_etap = 0.95778;
G_eta = 1.31e-6;   G_etap = 1.88e-4;
B0 = 0.13957^2/(0.0022 + 0.00469);
f  = 0.0924;
b  = B0*f;
p2 = mpi0^2;
D  = g(1) - g(2);
Sg = g(1) + g(2);

degenerate = g(1) == g(2) && g(2) == g(3);
if degenerate
    t0 = -b*sqrt(2/3)*sum(g)/(mphi^2 - m_etap^2);
    Rm = @(x) 1 ./ (x - m_etap^2 + 1i*m_etap*G_etap);
    Nphi = -2*t0*(sqrt(6)*b*Sg - p2*t0)/(3*f^2);
    Nj   =    t0*(sqrt(6)*b*Sg - 2*p2*t0)/(3*f^2);
else
    tp = -b*D/(mphi^2 - p2);
    te = -b*(g(1) + g(2) - 2*g(3))/sqrt(3)/(mphi^2 - m_eta^2);
    Re = @(x) 1 ./ (x - m_eta^2 + 1i*m_eta*G_eta);
end

% pi+pi-phi: A = P + Q*X, X = m_{+-}^2 - mpi^2.
mask = s > (2*mpi_pm + mphi)^2;
if any(mask)
    sp = s(mask);  sp = sp(:);
    Rf = 1 ./ (sp - mphi^2);
    if degenerate
        P = Nphi*Rf + Nj*Rm(sp);
        Q = zeros(size(sp));
    else
        Rp = 1 ./ (sp - p2);
        a  = mphi^2 + p2 + sp - 3*p2;       % x -> X + mpi^2 in the 3 thpi^2 x terms
        P  = Rf .* (-2*sqrt(3)*b*Sg*te - 2*b*D*tp + p2*te^2 - a*tp^2)/(3*f^2) ...
             + Rp .* tp .* (b*D + a*tp)/(3*f^2) ...
             + Re(sp) .* te*(sqrt(3)*b*Sg - p2*te)/(3*f^2);
        Q  = tp^2/f^2 * (Rf - Rp);
    end
    [I0, I1, I2] = dalitz_moments(sp, mphi, mpi_pm);
    W(mask) = reshape(abs(P).^2.*I0 + 2*real(P.*conj(Q)).*I1 + abs(Q).^2.*I2, size(W(mask)));
end

% pi0pi0phi.  Constant over the Dalitz plot.
mask = s > (2*mpi0 + mphi)^2;
if any(mask)
    sp = s(mask);  sp = sp(:);
    Rf = 1 ./ (sp - mphi^2);
    if degenerate
        A = Nphi*Rf + Nj*Rm(sp);
    else
        A = -Rf * (2*sqrt(3)*b*Sg*te + 6*b*D*tp - p2*te^2 - 3*p2*tp^2)/(3*f^2) ...
            + 1./(sp - p2) * tp*(b*D - p2*tp)/f^2 ...
            + Re(sp) * te*(sqrt(3)*b*Sg - p2*te)/(3*f^2);
    end
    W00(mask) = abs(A).^2 .* dalitz_moments(sp, mphi, mpi0) / 2;
end
W = W + W00;
end
