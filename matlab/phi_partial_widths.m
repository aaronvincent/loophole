function G = phi_partial_widths(mphi, mx, gx, g, interaction)
% phi_partial_widths  Per-channel partial widths of the scalar mediator phi.
%
%   G = phi_partial_widths(mphi, mx, gx, g, interaction)
%
% Returns a struct G with fields
%   gg     Gamma(phi -> gamma gamma)
%   xx     Gamma(phi -> chi chibar)
%   pipi   Gamma(phi -> pi+pi-, pi0pi0); for SP/PP and axialVector, Gamma(-> 3pi)
%   KK     Gamma(phi -> K+K-, K0K0bar)
%   eta    Gamma(phi -> eta eta, pi0 eta)
%   total  sum of the above
% The meson pairs come from meson_pairs.m (scalar phi-SM vertex: S(x)S, P(x)S).
%
% Used by xxggxsec.m / ggxxxsec.m to assemble the Breit-Wigner total
% width regulating the s-channel propagator, and by freeze_coupled.m
% for the phi decay source / sink in the coupled Boltzmann ODE.
%
% Channels per interaction (each gated on its threshold):
%   heavy/lightScalar     chi chibar (beta^3), gamma gamma (D, charged pi/K loop),
%                         meson pairs pi pi, K Kbar, eta eta, pi0 eta (meson_pairs.m)
%   heavy/lightScalarPS   chi chibar (beta^1), gamma gamma, meson pairs as above
%   heavy/lightScalarSP   chi chibar (beta^3), gamma gamma via pi0/eta(/eta') mixing,
%                         3pi (pipi field, LOOP_3PI=1); no two-meson channel (parity)
%   scalarPP, heavy/light chi chibar (beta^1), gamma gamma via mixing, 3pi (LOOP_3PI=1)
%   axialVector           chi chibar (beta^3, polarisation-averaged /12pi), 3pi
%                         (LOOP_3PI=1).  No gamma gamma (Landau-Yang), so below
%                         min(2 mx, 3 mpi) the A' is stable.
% Anything else returns all-zero.
mpi0    = 0.13498;
mpi_pm  = 0.13957;
mK      = 0.493677;
m_eta   = 0.547862;
alphaEM = 1/137;
B0      = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV (draft convention; was 2.6)
f       = 0.0924;
mN      = 0.938;
m_up    = 0.0022;
m_dn    = 0.00469;
m_s     = 0.095;
fpu     = 0.02;
fpd     = 0.026;
fps     = 0.118;

Dscript = 0.5*(B0*(g(1)+g(2))/mpi_pm^2 + B0*(g(1)+g(3))/mK^2)^2;   % charged-meson loop (pi+, K+), low-energy limit A0 -> 1/3

G = struct('gg', 0, 'xx', 0, 'pipi', 0, 'KK', 0, 'eta', 0, 'total', 0);

switch interaction
    case {'heavyScalar', 'lightScalar'}
        % Scalar phi-chi-chi, scalar phi-SM coupling.  Channels:
        %   gg   : D alpha^2 mphi^3 / (4608 pi^3)
        %   xx   : gx^2 (mphi^2 - 4 mx^2)^(3/2) / (8 pi mphi^2)     (beta^3 scalar vertex)
        %   meson pairs: below the switch
        G.gg   = alphaEM^2 * Dscript * mphi^3 / (4608*pi^3);
        % beta^3 scalar-vertex chi chi width, same form as the SP branch
        % (identical phi-chi-chi vertex); gate the 1/mphi^2 prefactor on
        % threshold so a closed channel gives 0, not Inf*0 = NaN.
        if mphi > 2*mx
            G.xx = gx^2 * (mphi^2 - 4*mx^2)^(3/2) / (8*pi*mphi^2);
        end
    case {'heavyScalarPS', 'lightScalarPS'}
        G.gg   = alphaEM^2 * Dscript * mphi^3 / (4608*pi^3);
        G.xx   = gx^2/(4*pi)               * sqrt(max(0, mphi^2/4 - mx^2));
    case {'scalarPP', 'heavyScalarPP', 'lightScalarPP'}
        % Pseudoscalar phi-chi-chi vertex, pseudoscalar phi-SM coupling.
        % Channels:
        %   gg : g_phigamgam^2 * mphi^3 / (64 pi)   (via pi0/eta mixing)
        %   xx : gx^2 * sqrt(mphi^2/4 - mx^2) / (4 pi)   (beta^1 pseudoscalar vertex)
        %   3pi: below the switch.  phi -> pi pi is parity forbidden.
        gphigg = phi_gg_coupling(mphi, gx, g, alphaEM, B0, f, mpi0, m_eta);

        G.gg = gphigg^2 * mphi^3 / (64*pi);
        G.xx = gx^2/(4*pi) * sqrt(max(0, mphi^2/4 - mx^2));
    case {'heavyScalarSP', 'lightScalarSP'}
        % Scalar phi-chi-chi vertex, pseudoscalar phi-SM coupling routed via
        % pi0/eta mixing.  Width terms (user-provided):
        %   gg  : g_phigamgam^2 * mphi^3 / (64 pi)
        %   xx  : gx^2 (mphi^2 - 4 mx^2)^(3/2) / (8 pi mphi^2)   (beta^3 scalar vertex)
        %   PP  : h_P^2 sqrt(mphi^2/4 - m_P^2) / (4 pi)          (pseudoscalar meson pairs)
        % g_phigamgam mirrors gphiprimegamma in ggxxxsec.m / xxggxsec.m SP branch.
        % h_P uses the nucleon-mass-weighted-quark form (same hp as PP branch),
        % shared between pi pi and K K channels.
        gphigg = phi_gg_coupling(mphi, gx, g, alphaEM, B0, f, mpi0, m_eta);
        hp = mN*(g(1)*fpu/m_up + g(2)*fpd/m_dn + g(3)*fps/m_s);

        G.gg   = gphigg^2 * mphi^3 / (64*pi);
        % beta^3 scalar-vertex width carries a 1/mphi^2 prefactor; gate on the
        % chi chi threshold so the massless / closed-channel limit is 0 rather
        % than 0/0 = NaN (cf. the SS pipi/KK guards; see docs/mediator_width.md).
        if mphi > 2*mx
            G.xx = gx^2 * (mphi^2 - 4*mx^2)^(3/2) / (8*pi*mphi^2);
        end
        % G.pipi = hp^2/(4*pi) * sqrt(max(0, mphi^2/4 - mpi0^2));
        % G.KK   = hp^2/(4*pi) * sqrt(max(0, mphi^2/4 - mK^2));
    case 'axialVector'
        % Axial-vector mediator A' (mphi = m_A', gx = lambda_a).  The only
        % decay channel is A' -> chi chibar (beta^3, axial-vector vertex):
        %   Gamma = lambda_a^2 (m_A'^2 - 4 mx^2)^(3/2) / (12 pi m_A'^2)
        % i.e. lambda_a^2 m_A' beta^3 / 12pi, the standard polarization-
        % AVERAGED axial-vector width (1/3 from averaging over the three
        % polarizations of a massive vector).
        % A' has no gamma gamma width (Landau-Yang) and no two-meson width:
        % the gamma gamma <-> chi chi amplitude runs through pi0/eta/eta'
        % poles, not an A' pole.  With LOOP_3PI=1 it gains Gamma(A' -> 3pi)
        % (below the switch), and the total then regulates the A' pole of the
        % J = 1 chi chi -> 3pi channel (axial3pi_sigma.m).
        % beta^3 1/mphi^2 prefactor: gate on threshold (massless A' is anyway
        % hard-errored upstream in axial_couplings.m, but keep this finite).
        if mphi > 2*mx
            G.xx = gx^2 * (mphi^2 - 4*mx^2)^(3/2) / (12*pi*mphi^2);
        end
end

% Two-meson channels of the scalar phi-SM vertex: Gamma = S mu^2 sqrt(lambda)
% / (16 pi mphi^3).  Gated on threshold so mphi = 0 gives 0, not Inf*0 = NaN.
if ismember(interaction, {'heavyScalar','lightScalar','heavyScalarPS','lightScalarPS'})
    P = meson_pairs(g);
    for j = 1:numel(P)
        mA = P(j).mA; mB = P(j).mB;
        if mphi > mA + mB
            lam = (mphi^2 - (mA+mB)^2) * (mphi^2 - (mA-mB)^2);
            G.(P(j).family) = G.(P(j).family) + P(j).S * P(j).mu^2 * sqrt(lam) / (16*pi*mphi^3);
        end
    end
end

% Three-pion channel of the pseudoscalar phi-SM vertex (threepi_weight.m).
% Gated on threshold so mphi = 0 gives 0, not 0/0 = NaN.
if mphi > 3*mpi0 && ismember(interaction, {'scalarPP','heavyScalarPP','lightScalarPP','heavyScalarSP','lightScalarSP'})
    G.pipi = threepi_weight(mphi^2, g) / (64*pi^3*mphi);
end
% Axial vector: only the transverse (J = 1) current decays an on-shell A'.
if mphi > 3*mpi0 && strcmp(interaction, 'axialVector')
    G.pipi = axial3pi_weight(mphi^2, g) / (192*pi^3*mphi);
end

G.total = G.gg + G.xx + G.pipi + G.KK + G.eta;
end

% --- effective phi'->gamma gamma coupling for SP/PP ----------------------
% Generic: phi' mixes with pi0 and eta (~ g_u-g_d, g_u+g_d-2g_s).  For
% g ~ [1 1 1] both vanish and the coupling is instead INDUCED via the eta'
% singlet (~ g_u+g_d+g_s); kept consistent with sppp_couplings.m's g_{phi' gg}.
% (The induced value is alpha^2-small, so this only matters for
% Gamma(phi->gg) itself, not as a BW regulator.)
function gphigg = phi_gg_coupling(mphi, gx, g, alphaEM, B0, f, mpi0, m_eta)
if g(1) == g(2) && g(2) == g(3)
    m_etap = 0.95778; G_etap_gg = 4.34e-6;
    F_etap = sqrt(alphaEM^2 * m_etap^3 / (64*pi^3*G_etap_gg));
    m_mix  = sppp_mixing(mphi, g);  m_mix = m_mix(3);   % draft's sqrt(6)/3 + NLO
    % alphaEM/(pi*F_etap): F_etap is a decay constant, so it carries the same
    % 1/pi as the gpigamma = alphaEM/(pi*f) convention in the else-branch below.
    % Kept in lockstep with sppp_couplings.m:57.  Fixed 2026-08-05.
    gphigg = (alphaEM/(pi*F_etap)) * m_mix / (mphi^2 - m_etap^2);
else
    gpigamma  = alphaEM/pi/f;
    getagamma = alphaEM/sqrt(3)/pi/f;
    % NLO meson-mediator mass mixings, from sppp_mixing.m (single source of truth).
    m2_all    = sppp_mixing(mphi, g);
    m2_pi0phi = m2_all(1);
    m2_etaphi = m2_all(2);
    % Both terms enter with the SAME sign: g_{phi' gg} = -sum_j g_{j gg} theta_j
    % with theta_j = m2_{j phi}/(mphi^2-mj^2).  The pi0 term used to carry the
    % opposite sign, which changed |gphigg| and hence Gamma(phi->gg).  Kept in
    % lockstep with sppp_couplings.m.
    gphigg = -gpigamma  * m2_pi0phi / (mphi^2 - mpi0^2) ...
             -getagamma * m2_etaphi / (mphi^2 - m_eta^2);
end
end
