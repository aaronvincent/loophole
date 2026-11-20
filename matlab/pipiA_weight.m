function [WT, WL] = pipiA_weight(s, mA, g)
% pipiA_weight  pi+pi-A' final state of the axial-vector A'-SM vertex.
%
%   [WT, WL] = pipiA_weight(s, mA, g)    s vector [GeV^2], mA = m_A', g = [a_u a_d a_s]
%
% Dalitz weights of the hadronic current V_sigma = <pi+ pi- A'|A_sigma|0>, summed
% over the outgoing A' polarizations, in the form axial3pi_sigma.m consumes:
%   WT  int dE_A dE_+ |V_spatial|^2   (CM frame of q)
%   WL  int dE_A dE_+ |q.V|^2
%
% Current, LO chiPT:
%   V_sigma = -2i Da Da0 [eps*_sigma - q_sigma (q.eps*)/(s - mpi0^2)],
% the A'A'pi+pi- contact plus the pi0 pole (A'-pi0 mixing, then the cubic
% A'pi+pi-pi0 vertex).  Da0 = a_u - a_d.  q.V is O(mpi^2), as current
% conservation in the chiral limit requires.  Summed over polarizations,
%   |V_spatial|^2 = 4 Da^2 Da0^2 (3 + |k_A|^2/mA^2)
%   |q.V|^2       = 4 Da^2 Da0^2 mpi0^4 s |k_A|^2 / (mA^2 (s - mpi0^2)^2).
% Both depend on E_A alone.  The 1/mA^2 is the longitudinal A', not suppressed
% by the pion mass.  Diagrams in which the outgoing A' comes off a pi0 line
% vanish, because the mixing vertex is k.eps = 0.  There is no pi0pi0A' final
% state, because the diagonal pi0 commutes with the diagonal source.  Nor is
% there an eta or eta' path.
%
% Normalization: L4 rescales the whole O(p^2) Lagrangian by corr4
% (axial_mixing.m), and so every tree amplitude by corr4^(V-I) = corr4 once.  The
% amplitude therefore carries Da = m_{pi0 A'}/F a single time, as in the 3pi
% current, and the other insertion stays LO Da0.  Rate: corr4^2 ~ 4 x LO.
%
% The final state contains an A', so using this cross section for production
% through detailed balance assumes A' is in the bath (as for pi pi phi).
%
% The channel is off when pipiphi_enabled() is false.

WT = zeros(size(s));  WL = zeros(size(s));
if ~pipiphi_enabled(), return, end

mpi0 = 0.13498;  mpi_pm = 0.13957;
f    = 0.0924;
mjA  = axial_mixing(g);
c4   = 4*(mjA(1)/f)^2*(g(1) - g(2))^2;

mask = s > (2*mpi_pm + mA)^2;
if any(mask)
    sp = s(mask);  sp = sp(:);
    % E_A = (c - X)/(2 sqrt s), c = s + mA^2 - mpi0^2, so |k_A|^2 is quadratic in X.
    [I0, I1, I2] = dalitz_moments(sp, mA, mpi_pm);
    c  = sp + mA^2 - mpi0^2;
    K  = (c.^2.*I0 - 2*c.*I1 + I2)./(4*sp) - mA^2*I0;
    WT(mask) = reshape(c4*(3*I0 + K/mA^2), size(WT(mask)));
    WL(mask) = reshape(c4*mpi0^4*sp.*K ./ (mA^2*(sp - mpi0^2).^2), size(WL(mask)));
end
end
