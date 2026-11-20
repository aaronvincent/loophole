function [mj, gcx, ggx, Gam] = sppp_couplings(mphi, mx, gx, g, interaction)
% sppp_couplings  s-channel pole couplings/widths for the SP and PP
% interactions, single source of truth for the gamma gamma <-> chi chi
% (and chi chi -> meson pair) amplitudes.
%
%   [mj, gcx, ggx, Gam] = sppp_couplings(mphi, mx, gx, g, interaction)
%
% Returns 1xN arrays (N = 3 generic, 2 in the degenerate case below) for the
% intermediate pseudoscalar resonances j:
%   mj   meson/mediator masses [GeV]
%   gcx  g_{j chi chi}    (coupling to the DM bilinear)
%   ggx  g_{j gamma gamma}
%   Gam  total BW width of pole j   (regulates 1/(mj^2 - s - i mj Gam))
%
% The cross sections are sigma ~ (s-kinematics)/(64 pi) * |sum_j gcx*ggx /
% (mj^2 - s - i mj Gam)|^2; consumers (ggxxxsec.m, xxggxsec.m, make_xxsigma.m,
% ggxxrate.m, svavg.m) loop over the returned arrays.
%
% Generic couplings (g not all equal): the {pi0, phi' (mediator), eta} set,
% with the chiral mixing couplings used previously inline.
%
% DEGENERATE case g_u = g_d = g_s (e.g. g = [1 1 1]/3): the pi0 and eta mixing
% couplings (~ g_u-g_d and g_u+g_d-2g_s) AND the induced phi'gamma coupling all
% vanish, so the generic amplitude is identically zero and chi decouples from
% photons -- unphysical.  The surviving channel is the SU(3) SINGLET eta'
% (~ g_u+g_d+g_s), normally neglected as the heaviest meson.  We switch to the
% {eta', phi'} set from docs/notes/SPPPetaprime.tex (small-mixing limit):
%   m_{eta' phi'} = B0 f/sqrt(6) (g_u+g_d+g_s)                  (mixing mass)
%   g_{eta' chi}  = -gx m_{eta' phi'} / (2(mphi^2 - m_eta'^2))
%   g_{phi' chi}  =  gx
%   g_{eta' gg}   =  alpha / (pi F_eta')   (data-driven; m_eta'^3 form, reused
%                   from axial_couplings.m -- the note's m_eta'^2 is a typo).
%                   The 1/pi matches the g_{pi0 gg} = alpha/(pi f) convention of
%                   the generic branch below; both are normalized by
%                   Gamma(j->gg) = g^2 m_j^3/(64 pi).  Was alpha/F_eta' until
%                   2026-08-05, i.e. pi too large -- and since BOTH ggx entries
%                   here descend from it, that was pi^2 ~ 9.87 on every SP/PP
%                   cross section in the degenerate direction (the default
%                   g = [1 1 1]/3).
%   g_{phi' gg}   =  g_{eta' gg} m_{eta' phi'} / (2(mphi^2 - m_eta'^2))   (induced)
% eta' BW width = SM total + Gamma(eta'->chichi): beta^1 for PP (pseudoscalar
% chi vertex), beta^3 for SP (scalar chi vertex).

mpi0    = 0.13498;   Gpi    = 7.7e-9;
m_eta   = 0.547862;  G_eta  = 1.31e-6;
m_etap  = 0.95778;   G_etap = 1.88e-4;   % eta' SM total width
G_etap_gg = 4.34e-6;                     % eta' -> gamma gamma partial width
alphaEM = 1/137;
B0      = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV (draft convention; was 2.6)
f       = 0.0924;
L6      = 1.46e-3;   % NLO chiral low-energy constants (mixing-mass NLO terms)
L7      = -0.39e-3;
m_up    = 0.0022;    % GeV, light-quark masses (match axial_couplings.m)
m_dn    = 0.00469;
m_s     = 0.095;

isPP = ismember(interaction, {'scalarPP','heavyScalarPP','lightScalarPP'});

Gphi = phi_partial_widths(mphi, mx, gx, g, interaction).total;

if g(1) == g(2) && g(2) == g(3)
    % --- degenerate g ~ [1 1 1]: eta' singlet + mediator -----------------
    F_etap      = sqrt(alphaEM^2 * m_etap^3 / (64*pi^3 * G_etap_gg));
    g_etapgamma = alphaEM / (pi*F_etap);   % 1/pi: same convention as gpigamma below
    % m_{eta' phi}^2 from sppp_mixing.m (single source of truth).  Was
    % B0*f/sqrt(6)*sum(g) -- the "all equal coupling" appendix's LO form, 4.3x
    % low -- until 2026-08-12, when that appendix was commented out of the draft
    % leaving only the body's sqrt(6)/3 * sum(g) * (B0 F + NLO).
    m2_all      = sppp_mixing(mphi, g);
    m_mix       = m2_all(3);
    denom       = (mphi^2 - m_etap^2);

    g_etapx        = -gx * m_mix / denom;
    gphiprimex     =  gx;
    gphiprimegamma =  g_etapgamma * m_mix / denom;   % induced phi'->gamma gamma

    if isPP   % pseudoscalar chi vertex, beta^1
        Gam_etapxx = g_etapx^2 * sqrt(max(0, m_etap^2/4 - mx^2)) / (4*pi);
    else      % SP: scalar chi vertex, beta^3
        Gam_etapxx = g_etapx^2 * max(0, m_etap^2 - 4*mx^2)^(3/2) / (8*pi*m_etap^2);
    end
    Gam_etap = G_etap + Gam_etapxx;

    mj  = [m_etap,      mphi          ];
    gcx = [g_etapx,     gphiprimex    ];
    ggx = [g_etapgamma, gphiprimegamma];
    Gam = [Gam_etap,    Gphi          ];
else
    % --- generic couplings: pi0, mediator phi', eta ----------------------
    % Pseudoscalar meson-mediator mass mixings m_{j phi}^2 (NLO, L6/L7 terms).
    % The eta expression is 1/sqrt(3) times the pi0 bracket, term by term: the
    % chiral trace Tr[M_g Phi] = (g_u-g_d) pi0/sqrt(2) + (g_u+g_d-2g_s) eta8/sqrt(6)
    % gives the eta8 coefficient a relative sqrt(2)/sqrt(6) = 1/sqrt(3).  This is
    % the same algebra behind the sqrt(3)F/3 in m_{eta A'} (axial_couplings.m:73).
    m2_all    = sppp_mixing(mphi, g);
    m2_pi0phi = m2_all(1);
    m2_etaphi = m2_all(2);

    % Signs: the mixing angle is theta_j = m2_{j phi}/(mphi^2-mj^2) for BOTH
    % poles, so g_{j chi} = +gx theta_j with the SAME sign for pi0 and eta, and
    % the induced g_{phi' gg} = -sum_j g_{j gg} theta_j.  Each meson pole and its
    % share of the mediator pole then enter with opposite sign, which is the
    % mixing-insertion structure A ~ [1/(s-mj^2) - 1/(s-mphi^2)].  Until
    % 2026-08-10 gpix carried a spurious minus and gphiprimegamma's pi0 term a
    % spurious plus, so the pi0 and eta blocks interfered with the WRONG relative
    % sign.  This is now identical to
    % mesonandCMB_ratio.m's sppp_poles_indirect, which always had it right.
    gpigamma       = alphaEM/pi/f;
    getagamma      = alphaEM/sqrt(3)/pi/f;
    gpix           = gx*m2_pi0phi / (mphi^2-mpi0^2);
    gphiprimex     = gx;
    gphiprimegamma = -gpigamma *m2_pi0phi/(mphi^2-mpi0^2) ...
                     -getagamma*m2_etaphi/(mphi^2-m_eta^2);
    getax          = gx*m2_etaphi / (mphi^2-m_eta^2);

    mj  = [mpi0,     mphi,           m_eta   ];
    gcx = [gpix,     gphiprimex,     getax   ];
    ggx = [gpigamma, gphiprimegamma, getagamma];
    Gam = [Gpi,      Gphi,           G_eta   ];
end
end
