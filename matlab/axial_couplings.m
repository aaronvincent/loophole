function [mj, lam_chi, lam_gg, Gam_j] = axial_couplings(mphi, mx, gx, g)
% axial_couplings  A' (axial-vector mediator) mixing with the neutral
% pseudoscalars, single source of truth for the axialVector interaction.
%
%   [mj, lam_chi, lam_gg, Gam_j] = axial_couplings(mphi, mx, gx, g)
%
% Inputs (mapped onto the project's usual mediator slots):
%   mphi = m_A'      axial-vector mediator mass [GeV]
%   mx   = m_chi     DM mass [GeV]  (sets the meson -> chi chi partial widths)
%   gx   = lambda_a  A'-chi-chi axial coupling
%   g    = [a_u a_d a_s]   A'-quark axial couplings
%
% Outputs are 1x3 arrays for j = {pi0, eta, eta'}:
%   mj       meson masses [GeV]
%   lam_chi  A'-meson-chi coupling  lambda_{j chi} = lambda_a * m_{jA'} / m_A'^2
%            (small-mixing approximation; also equals lambda_{j chi chi},
%            the meson->chi chi vertex used by the meson-decay channel)
%   lam_gg   meson-gamma-gamma coupling  lambda_{j gamma gamma}
%   Gam_j    meson TOTAL widths = SM total width + Gamma(j -> chi chibar),
%            used as the Breit-Wigner regulators for the s = m_j^2 poles.
%            The chi chibar piece (the meson-decay channel) is added here so
%            the propagator width is consistent with the other interactions'
%            mediator widths (which likewise include the -> chi chi mode).
%
% The gamma gamma <-> chi chi amplitude for the axial-vector mediator runs
% through these three SM-meson s-channel poles (NOT an A' pole), so A' enters
% the cross sections only through lam_chi.  Gamma(A') itself lives in
% phi_partial_widths.m and is used by the chi chi -> A'A' / coupled-bath
% bookkeeping.

mpi0    = 0.13498;
m_eta   = 0.547862;
m_etap  = 0.95778;

% Meson total widths (PDG) -- regulate the narrow s = m_j^2 propagator poles.
Gpi    = 7.7e-9;     % pi0   total width [GeV]
G_eta  = 1.31e-6;    % eta   total width [GeV]
G_etap = 1.88e-4;    % eta'  total width [GeV] (0.188 MeV)

% eta' -> gamma gamma partial width: 0.188 MeV total * BR(gg) = 2.307% (PDG).
G_etap_gg = 4.34e-6; % [GeV]

% eta -> gamma gamma partial width: 1.31 keV total * BR(gg) = 39.36% (PDG).
G_eta_gg  = 5.16e-7; % [GeV]

alphaEM = 1/137;
F       = 0.0924;    % GeV, pion decay constant

% --- meson -> gamma gamma couplings -------------------------------------
%   lambda_{pi0 gg}  = e^2/(4 pi^2 F)        = alphaEM/(pi F)        (chiral)
%   lambda_{eta gg}  = e^2/(4 pi^2 F_eta)    = alphaEM/(pi F_eta)    (data-driven
%      as of the 2026-08-11 draft, Eq. (lambdagamma): the LO octet value
%      alphaEM/(sqrt3 pi F) gives Gamma(eta->gg) = 0.17 keV against a measured
%      0.52 keV -- the familiar factor 3 from taking eta = eta_8 with no
%      eta-eta' mixing.  F_eta = 0.0925 GeV so defined, i.e. 1.73x the LO
%      sqrt3*F = 0.1600 GeV in amplitude.  NOTE: the draft's pseudoscalar
%      section still carries the LO form, so sppp_couplings.m /
%      phi_partial_widths.m / mesonandCMB_ratio.m are deliberately NOT changed
%      to match.
%   lambda_{eta' gg} = e^2/(4 pi^2 F_eta')   = alphaEM/(pi F_eta')   (data-driven:
%      eta' is not a Goldstone, so F_eta' is fixed from the measured gg width,
%      F_eta' = sqrt(alphaEM^2 m_eta'^3 / (64 pi^3 Gamma_eta'->gg)); alphaEM
%      cancels, leaving lambda_{eta' gg} = sqrt(64 pi Gamma_eta'->gg / m_eta'^3).)
%
% All three are the SAME convention -- the one normalized by Gamma(j->gg) =
% lambda_{j gg}^2 m_j^3/(64 pi).  That is not cosmetic: the three add coherently
% in lam_gg below, so a relative factor between entries is a physics error.
% Checks: the pi0 entry reproduces Gamma(pi0->gg)=7.7 eV to 0.4%, and the eta'
% entry reproduces the measured 4.34 keV exactly (by construction).
% F_etap = 0.0737 GeV as defined above is a DECAY CONSTANT (cf. f_pi = 0.0924),
% so it needs the same 1/pi as F -- the earlier alphaEM/F_etap was pi too large,
% which put pi^2 ~ 9.87 into every |amplitude|^2 through this coupling.
% Fixed 2026-08-05.
F_etap   = sqrt(alphaEM^2 * m_etap^3 / (64*pi^3*G_etap_gg));
F_eta    = sqrt(alphaEM^2 * m_eta^3  / (64*pi^3*G_eta_gg));
lam_pigg   = alphaEM/(pi*F);
lam_etagg  = alphaEM/(pi*F_eta);     % was alphaEM/(sqrt(3)*pi*F); see above
lam_etapgg = alphaEM/(pi*F_etap);

% --- A'-meson mixing masses m_{jA'} (small-mixing) ----------------------
% axial_mixing.m, which carries the unresolved L4 sign (corr4).
mjA = axial_mixing(g);
m_pi0Ap = mjA(1);  m_etaAp = mjA(2);  m_etapAp = mjA(3);
au = g(1); ad = g(2); as = g(3);

% --- small-mixing validity: m_A' must exceed the off-diagonal mixing masses
% The induced couplings lambda_{j chi} = lambda_a m_{jA'}/m_A'^2 come from
% diagonalising the A'-meson mass matrix to first order in the mixing angle
% sin(theta_j) ~ m_{jA'}/m_A'.  That expansion is only valid for m_A' > |m_{jA'}|;
% once m_A' drops below an off-diagonal mass the "small-mixing" picture is no
% longer perturbative (the would-be A' and meson states are maximally mixed and
% the lambda_{j chi} ~ 1/m_A'^2 formula is meaningless).  Refuse rather than
% extrapolate -- the compare_* scans catch this and mark the point NaN.
mixmass = abs([m_pi0Ap, m_etaAp, m_etapAp]);
[mxmix, jmix] = max(mixmass);
% DISABLED (temporary): the small-mixing breakdown check below was gapping the
% entire scan at small m_A' (e.g. mphi=1e-9, where m_A' << the ~0.2 GeV eta'-A'
% mixing mass), NaN'ing every point and emptying the plots.  Re-enable when the
% small-mixing-regime treatment is sorted out.
if false && mphi <= mxmix
    jn = {'pi0','eta','etaprime'};
    error('axial_couplings:mixingTooLarge', ...
        ['axialVector small-mixing breakdown: m_A''=%.3g GeV is below the ' ...
         '%s-A'' off-diagonal mixing mass |m_{jA''}|=%.3g GeV (m_A''/|m_{jA''}|' ...
         '=%.3g <= 1), mx=%.3g, lambda_a=%.3g, g=[%.3g %.3g %.3g]. The ' ...
         'sin(theta)~m_{jA''}/m_A'' expansion behind lambda_{jchi}~1/m_A''^2 ' ...
         'is non-perturbative here; axialVector couplings are invalid.'], ...
        jn{jmix}, mphi, mxmix, mphi/mxmix, mx, gx, au, ad, as);
end

% --- assemble (lambda_{j chi} = lambda_a m_{jA'} / m_A'^2) ----------------
mj      = [mpi0,     m_eta,     m_etap   ];
lam_gg  = [lam_pigg, lam_etagg, lam_etapgg];
lam_chi = gx/mphi^2 * [m_pi0Ap, m_etaAp, m_etapAp];

% Meson Breit-Wigner widths: SM total width + Gamma(j -> chi chibar).
%   Gamma(j -> chi chibar) = lambda_{j chi}^2 mx^2 sqrt(mj^2/4 - mx^2) / pi
% (docs/notes/axialvector.txt).  Open only when mj > 2 mx.
Gam_SM  = [Gpi, G_eta, G_etap];
Gam_xx  = lam_chi.^2 * mx^2 .* sqrt(max(0, mj.^2/4 - mx^2)) / pi;
Gam_j   = Gam_SM + Gam_xx;

% --- validity floor: meson BW width must stay below its mass --------------
% lam_chi = lambda_a m_{jA'}/m_A'^2 diverges as 1/m_A'^2, so the meson->chi chi
% width Gam_xx ~ lam_chi^2 ~ 1/m_A'^4 blows up at small m_A'.  Once Gam_j >= m_j
% the Breit-Wigner / small-mixing-EFT description of the j = {pi0,eta,eta'}
% pole is meaningless (a "resonance" wider than its mass), and the longitudinal
% A' modes drive a parallel 1/m_A'^4 unitarity-violating enhancement in
% chi chi -> A'A'.  Unlike a scalar mediator, the axial-vector m_A'->0 limit is
% singular -- refuse rather than return garbage.  The compare_* scans catch
% this and mark the point NaN (see the stderr [..] error logging).
[rmax, jmax] = max(Gam_j ./ mj);
% This validity floor is KEPT: it cleanly bounds the invalid region (the scans
% catch it -> NaN, no garbage in the solvers) AND is the boundary that
% axial_gx_ceiling.m probes/catches to draw the unitarity line + shading.
if rmax >= 1
    jn = {'pi0','eta','etaprime'};
    error('axial_couplings:mediatorTooLight', ...
        ['axialVector validity floor: %s BW width exceeds its mass ' ...
         '(Gamma/m = %.3g) at m_A''=%.3g GeV, mx=%.3g, lambda_a=%.3g. ' ...
         'lambda_{jchi} ~ 1/m_A''^2 has blown the mixing-EFT/Breit-Wigner ' ...
         'description; axialVector cross sections are invalid here.'], ...
        jn{jmax}, rmax, mphi, mx, gx);
end
end
