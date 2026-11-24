function [hp, htildeP, a_N, htildeN, hn] = nucleon_couplings(g)
% nucleon_couplings  Effective DM/mediator-NUCLEON couplings from the quark
% couplings g = [g_u g_d g_s].  Single source of truth shared by DMnucxsec.m
% (sigma_chiN) and external_loglikes.m (SN1987A / neutron-star / axial bounds),
% so the nucleon-coupling maps never drift between them.
%
%   [hp, htildeP, a_N, htildeN, hn] = nucleon_couplings(g)
%
%   hp       scalar PROTON coupling    = mN * sum_i g_i f_{p,i} / m_i
%            (SS / PS interactions; sigma-term form factors).
%   htildeP  pseudoscalar nucleon coupling = mN * sum_i (Dn_i/m_i)(g_i - g_proj)
%            (SP / PP interactions; proton spin content with the SU(3)-singlet
%            (eta') piece projected out -- g ~ [1 1 1] gives htildeP = 0).
%   a_N      axial nucleon coupling    = sum_i g_i Delta_i^p  (axialVector).
%   htildeN  pseudoscalar NEUTRON coupling -- htildeP with the u and d spin
%            contents swapped (isospin).  Needed by the pseudoscalar 5th-force
%            bound, which constrains g_P^2 and g_P*g_N separately.
%   hn       scalar NEUTRON coupling -- hp with the u and d sigma-term form
%            factors swapped (isospin); each 1/m_i stays with its own quark.
%            Needed by the Hardy & Lasenby stellar-cooling bound (TRGB + HB),
%            which assumes equal phi-proton / phi-neutron couplings, so we feed
%            it (hp + hn)/2.
%
% hp, htildeP and a_N are PROTON values; the neutron swaps u<->d.  References:
% sigma terms arXiv:1202.1292 (table II); spin contents FLAG 2024 (below).

% Proton quark spin contents Delta_q^p = g_A^q, FLAG 2024 (arXiv:2411.04268),
% N_f = 2+1 (eqs. 440-442).  ONE set for the axial (a_N) and pseudoscalar
% (htildeP, htildeN) couplings: both are the same matrix element
% <N|qbar gamma^mu gamma5 q|N>.  Neutron: swap u <-> d.
Du = 0.847;  Dd = -0.407;  Ds = -0.035;
mN  = 0.938;
% Scalar sigma-term form factors
fpu = 0.02;  fpd = 0.026;  fps = 0.118;
% Light-quark masses [GeV]
m_up = 0.0022;  m_dn = 0.00469;  m_s = 0.095;

% Scalar nucleon coupling (+ heavy-quark term neglected)
hp = mN*(g(1)*fpu/m_up + g(2)*fpd/m_dn + g(3)*fps/m_s);
% Neutron: isospin swaps the u<->d sigma-term form factors, each 1/m_i kept with
% its own quark (same rule as htildeN below).  So hn is NOT hp([g2 g1 g3]).
hn = mN*(g(1)*fpd/m_up + g(2)*fpu/m_dn + g(3)*fps/m_s);

% Pseudoscalar nucleon coupling: proton spin content Dn_i weighted by 1/m_i,
% with the SU(3)-singlet (eta') direction projected out -- (g_i - <g/m>/<1/m>).
Dnu = Du;  Dnd = Dd;  Dns = Ds;
invm_sum = 1/m_up + 1/m_dn + 1/m_s;
g_proj   = (g(1)/m_up + g(2)/m_dn + g(3)/m_s) / invm_sum;
htildeP  = mN*( Dnu/m_up*(g(1)-g_proj) ...
              + Dnd/m_dn*(g(2)-g_proj) ...
              + Dns/m_s *(g(3)-g_proj) );
% Neutron: isospin swaps the u and d SPIN CONTENTS only -- each 1/m_i stays with
% its own quark, and g_proj (a property of the coupling direction, not of the
% nucleon) is unchanged.  So htildeN is NOT nucleon_couplings([g2 g1 g3]).
htildeN  = mN*( Dnd/m_up*(g(1)-g_proj) ...
              + Dnu/m_dn*(g(2)-g_proj) ...
              + Dns/m_s *(g(3)-g_proj) );

% Axial nucleon coupling
a_N = g(1)*Du + g(2)*Dd + g(3)*Ds;
end
