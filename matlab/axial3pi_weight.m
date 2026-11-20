function [WT, WL] = axial3pi_weight(s, g)
% axial3pi_weight  Three-pion final states of the axial-vector A'-SM vertex.
%
%   [WT, WL] = axial3pi_weight(s, g)     s vector [GeV^2], g = [a_u a_d a_s]
%
% Returns the Dalitz weights W = sum_ch S int dE1 dE2 |.|^2 of the two parts of
% the current V^mu = <3pi|A^mu|0>, times Da^2 with Da = m_{pi0 A'}/F:
%   WT  transverse, |V_T|^2 = 4 |p_0|^2/F^2 (CM frame), pi+pi-pi0 only
%   WL  longitudinal, |q.V|^2, over pi+pi-pi0 (S = 1) and 3pi0 (S = 1/6)
% The consumers are
%   sigma(chi chi -> 3pi) = lambda_a^2/(32 pi^3) [beta WT/(6 |D_A|^2)
%                           + mx^2 WL/(beta s^2 m_A'^4)]     axial3pi_sigma.m
%   Gamma(A' -> 3pi)      = WT(m_A'^2)/(192 pi^3 m_A')      phi_partial_widths.m
%
% Current: LO chiPT.  With X = m_{+-}^2 - mpi^2,
%   q.V(+-0) = (mpi^2/F) X/(s - mpi^2),   q.V(000) = (mpi^2/F) s/(s - mpi^2),
% the pi0-pole piece of the pseudoscalar amplitude times mpi^2/(B0 (g_u - g_d)),
% as PCAC requires.  The transverse part is the contact term of the cubic
% current alone; it is J = 1 (a1-like) and vanishes for 3pi0.  Only the
% isovector a_u - a_d reaches 3pi: the isoscalar part drops out of the current,
% and the isospin-breaking eta -> 3pi path that threepi_weight.m keeps is
% omitted here.  The current is normalized by m_{pi0 A'}/F from axial_mixing.m
% (LO: a_u - a_d), the same A'-pi0 coupling the gamma gamma channel uses.  The
% L4 term multiplies the whole LO current, pole and contact alike, so unlike
% the pseudoscalar case (threepi_weight.m) the NLO value keeps the PCAC
% cancellation; the L5 term is a ~1% non-uniform piece taken along with it.
%
% The channels are off when threepi_enabled() is false.

WT = zeros(size(s));  WL = zeros(size(s));
if ~threepi_enabled(), return, end

mpi0 = 0.13498;  mpi_pm = 0.13957;
f    = 0.0924;
mjA  = axial_mixing(g);
Da2  = (mjA(1)/f)^2;
cL   = Da2 * mpi0^4 / f^2;

% pi0pi0pi0: longitudinal only, constant over the Dalitz plot.
mask = s > 9*mpi0^2;
if any(mask)
    sp = s(mask);  sp = sp(:);
    WL(mask) = cL * sp.^2 ./ (sp - mpi0^2).^2 .* dalitz_moments(sp, mpi0, mpi0) / 6;
end

% pi+pi-pi0.  E_0 = (s - X)/(2 sqrt s), so |p_0|^2 = (s - X)^2/(4 s) - mpi0^2
% is a quadratic in X and the moments of X are all that is needed.
mask = s > (2*mpi_pm + mpi0)^2;
if any(mask)
    sp = s(mask);  sp = sp(:);
    [I0, I1, I2] = dalitz_moments(sp, mpi0, mpi_pm);
    p02 = (sp.^2.*I0 - 2*sp.*I1 + I2)./(4*sp) - mpi0^2*I0;
    WT(mask) = reshape(4*Da2/f^2 * p02, size(WT(mask)));
    WL(mask) = WL(mask) + reshape(cL * I2 ./ (sp - mpi0^2).^2, size(WL(mask)));
end
end
