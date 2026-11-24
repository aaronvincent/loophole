function [gq_mix, mxmix, jbind] = axial_mixing_floor(mphi, mx, g, smax)
% axial_mixing_floor  Small-mixing-breakdown boundary for the axialVector
% interaction (companion of axial_gq_ceiling / axial_gx_ceiling).
%
%   [gq_mix, mxmix, jbind] = axial_mixing_floor(mphi, mx, g, smax)
%
% The axial-vector meson->chi chi coupling lambda_{j chi} = lambda_a m_{jA'}/m_A'^2
% comes from diagonalising the A'-meson mass matrix to first order in the mixing
% angle sin(theta_j) ~ m_{jA'}/m_A'.  That leading-order (small-mixing) expansion
% degrades as the mixing angle grows and is meaningless once m_A' <~ m_{jA'}.
% Here we return the g_q scan boundary at which the largest mixing angle reaches
% a threshold smax (default 0.3):
%
%   sin(theta_j) = g_q*|m_{jA'}(g)|/m_A' = smax   ->   gq_mix = smax*mphi/mxmix,
%
% with mxmix = max_j |m_{jA'}(g)| (the eta'-A' mixing mass binds).  In the
% g_q scan the quark couplings scale as g_eff = g_q*g and the mixing masses are
% linear in g, so gq_mix is mx-independent.  Outputs:
%   gq_mix  boundary g_q (sin(theta)=smax); NaN if it cannot be evaluated
%   mxmix    max_j |m_{jA'}(g)|  [GeV]  (so the gx-scan can form sin(theta)=mxmix/mphi)
%   jbind    index of the binding meson in {pi0, eta, eta'}
%
% Unlike axial_couplings.m's mixingTooLarge gate (disabled), this NEVER errors
% and NEVER gates a scan -- it is a drawn-overlay diagnostic only.
%
% axial_couplings.m is the single source of truth for the mixing masses, so
% (as in axial_gx_ceiling.m) we probe it at a tiny lambda_a -- which clears the
% mediatorTooLight floor without changing the lambda_a-independent mixing masses
% -- and recover |m_{jA'}| = |lambda_{j chi}| m_A'^2 / lambda_a.

if nargin < 4 || isempty(smax)
    smax = 0.3;
end

% Probe at a lambda_a small enough to clear axial_couplings' validity floor.
gx_probe = 1;
lam_chi = [];
ok = false;
for it = 1:200
    try
        [~, lam_chi, ~, ~] = axial_couplings(mphi, mx, gx_probe, g);
        ok = true;
        break;
    catch ME
        if strcmp(ME.identifier, 'axial_couplings:mediatorTooLight')
            gx_probe = gx_probe * 1e-6;   % step well below the floor
            continue;
        else
            rethrow(ME);
        end
    end
end
if ~ok
    gq_mix = NaN; mxmix = NaN; jbind = NaN;
    return;
end

% |m_{jA'}| = |lambda_{j chi}| m_A'^2 / lambda_a  (lambda_a-independent).
m_jAp = abs(lam_chi) * mphi^2 / gx_probe;
[mxmix, jbind] = max(m_jAp);

if mxmix > 0
    gq_mix = smax * mphi / mxmix;
else
    gq_mix = Inf;   % no mixing (e.g. couplings cancel) -> no breakdown
end
end
