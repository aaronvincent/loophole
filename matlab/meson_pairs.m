function P = meson_pairs(g)
% meson_pairs  Two-meson final states of a scalar phi coupled to quark masses.
%
%   P = meson_pairs(g)    g = [g_u g_d g_s]
%
% Lowest-order chiPT with M_q -> M_q + g*phi gives L ⊃ -B0 phi tr(G Phi^2),
% G = diag(g), Phi the octet (eta = eta_8, no eta-eta' mixing).  Expanding:
%   pi+ pi-   B0 (g_u+g_d)          K+ K-    B0 (g_u+g_s)
%   pi0 pi0   B0 (g_u+g_d)/2        K0 K0bar B0 (g_d+g_s)
%   eta eta   B0 (g_u+g_d+4g_s)/6   pi0 eta  B0 (g_u-g_d)/sqrt(3)
% as coefficients of -phi A B.  The vertex amplitude mu doubles the identical
% pairs' coefficient; S = 1/2 is their phase-space factor.  So for any pair
%   Gamma(phi -> A B)       = S mu^2 sqrt(lambda(mphi^2,mA^2,mB^2)) / (16 pi mphi^3)
% and the chi chi -> A B cross sections in xxggxsec.m carry S mu^2 the same way.
%
% Returns a struct array with fields name, mA, mB, mu, S, family
% (family = 'pipi', 'KK' or 'eta', the phi_partial_widths field it sums into).

B0 = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV
mpi_pm = 0.13957; mpi0 = 0.13498; mK = 0.493677; mK0 = 0.497611; m_eta = 0.547862;
gu = g(1); gd = g(2); gs = g(3);

P = struct( ...
    'name',   {'pi+pi-',      'pi0pi0',      'K+K-',        'K0K0bar',     'etaeta',              'pi0eta'}, ...
    'mA',     {mpi_pm,        mpi0,          mK,            mK0,           m_eta,                 mpi0}, ...
    'mB',     {mpi_pm,        mpi0,          mK,            mK0,           m_eta,                 m_eta}, ...
    'mu',     {B0*(gu+gd),    B0*(gu+gd),    B0*(gu+gs),    B0*(gd+gs),    B0*(gu+gd+4*gs)/3,     B0*(gu-gd)/sqrt(3)}, ...
    'S',      {1,             1/2,           1,             1,             1/2,                   1}, ...
    'family', {'pipi',        'pipi',        'KK',          'KK',          'eta',                 'eta'});
end
