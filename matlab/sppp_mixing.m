function [m2, theta, theta_max, mj] = sppp_mixing(mphi, g)
%SPPP_MIXING  phi-meson mass mixings and small-mixing angles for the SP/PP sector.
%
%   [m2, theta, theta_max, mj] = SPPP_MIXING(MPHI, G)
%
%   m2         1x3  [m^2_{pi0 phi}, m^2_{eta phi}, m^2_{etap phi}]   [GeV^2]
%   theta      1x3  small-mixing angles theta_j = m2_j / (mphi^2 - mj^2)
%   theta_max  max(abs(theta))  -- the expansion parameter of the whole sector
%   mj         1x3  [mpi0, m_eta, m_etap]
%
%   Single source of truth for the pseudoscalar mediator's mass mixing with the
%   neutral mesons.  Used by sppp_couplings.m (which builds the pole couplings
%   from these) and by model_observables.m (which refuses points where the
%   expansion the whole sector rests on has broken down).
%
%   CONVENTION (2026-08-12 draft, Eq. eq:pseudocouplingslim and the appendix
%   mixing-angle block).  The mixing angle carries NO factor 2:
%
%       theta_j = m^2_{j phi} / (m_j^2 - m_phi^2)
%
%   and the couplings that follow are g_{j chi} = +gx*theta_j (same sign for
%   every j) with the induced g_{phi gg} = -sum_j g_{j gg} theta_j.  This code
%   uses the opposite overall sign, theta_j = m2_j/(mphi^2 - mj^2), which is an
%   overall sign on the whole amplitude and so drops out of |A|^2.  Earlier
%   drafts (and this code until 2026-08-12) carried a 1/2 in every denominator;
%   removing it doubles every SP/PP coupling and so multiplies the SP/PP cross
%   sections and the K+ -> pi+ chi chibar rate by 4.
%
%   VALIDITY.  The draft states the small-mixing limit as
%       |m^2_{pi0 phi}| << |m_phi^2 - m_pi0^2|,  and likewise for eta and eta',
%   i.e. |theta_j| << 1 for all three poles simultaneously.  Everything in the
%   SP/PP sector -- the pole couplings, the induced phi-gamma-gamma coupling,
%   the kaon amplitude -- is first order in theta, so theta_max is the honest
%   expansion parameter and model_observables.m cuts on it.
%
%   The eta' is included even though the generic pole set is {pi0, phi', eta}:
%   the condition is on the validity of the DIAGONALIZATION, not on which poles
%   are kept in the amplitude.  A large eta'-phi mixing invalidates the rotation
%   regardless.  It is also the largest of the three, since m^2_{etap phi} goes
%   as (g_u+g_d+g_s) rather than a difference of couplings.

mpi0   = 0.13498;
m_eta  = 0.547862;
m_etap = 0.95778;

B0   = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV
f    = 0.0924;
L6   = 1.46e-3;
L7   = -0.39e-3;
m_up = 0.0022;
m_dn = 0.00469;
m_s  = 0.095;

summ = m_up + m_dn + m_s;
% NLO bracket shared by all three (O(p^2) B0 F plus the O(p^4) L6 piece).
% NOTE: this bracket is 2.58x its own leading term with the draft's constants,
% which is not a convergent chiral expansion, but is used as written here.
K    = B0*f + 32*B0^2*L6/f*summ;
sg   = g(1) + g(2) + g(3);

m2_pi0phi  = (g(1)-g(2))*K + 32*B0^2*L7/f * sg * (m_up - m_dn);
m2_etaphi  = ( (g(1)+g(2)-2*g(3))*K ...
             + 32*B0^2*L7/f * sg * (m_up + m_dn - 2*m_s) ) / sqrt(3);
% eta' (SU(3) singlet).  sqrt(6)/3 = sqrt(2/3) is the eta_0 coefficient of
% Tr[g_q Pi]; the L6+L7 bracket is the draft's.  This replaced B0*f/sqrt(6)
% (LO only, and sqrt(6)/3 -> 1/sqrt(6), together 4.3x low) on 2026-08-12, when
% the appendix carrying that form was commented out of the draft.
m2_etapphi = sqrt(6)/3 * sg * (B0*f + (L6 + L7)*32*B0^2/f*summ);

mj    = [mpi0, m_eta, m_etap];
m2    = [m2_pi0phi, m2_etaphi, m2_etapphi];
theta = m2 ./ (mphi^2 - mj.^2);
theta_max = max(abs(theta));
end
