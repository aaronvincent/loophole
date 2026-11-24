function gq_max = axial_gq_ceiling(mphi, mx, gx, g)
% axial_gq_ceiling  EFT-validity ceiling on the g_q scan normalization for the
% axialVector interaction (companion of axial_gx_ceiling, for compare_limits_g).
%
%   gq_max = axial_gq_ceiling(mphi, mx, gx, g)
%
% In the g_q scan gx is FIXED and g_q scales the quark couplings (g_eff = g_q*g),
% so BOTH axial_couplings.m validity gates move with g_q and the gx-solving
% axial_gx_ceiling cannot be reused (at g_q=1 the small-mixing gate is already
% tripped for an ultralight m_A', so its sub-floor probe never clears and it
% returns NaN).  The gate is mediatorTooLight: the meson->chichi width reaches
% its mass; the width ~ lambda_{jchi}^2 ~ (gx*g_q)^2, so it grows with g_q.
% It is monotonic in g_q, so the valid set is g_q < gq_max and we bisect on
% log10(g_q) for the boundary.  Returns Inf if g_q=1e12 is still valid (no
% ceiling, e.g. heavy m_A'); NaN if even g_q=1e-30 is invalid.

valid = @(le) local_valid(mphi, mx, gx, 10^le, g);

lo = -30; hi = 12;
if valid(hi),  gq_max = Inf;  return; end   % no ceiling in range
if ~valid(lo), gq_max = NaN;  return; end    % invalid even at vanishing gq

for it = 1:60
    mid = 0.5 * (lo + hi);
    if valid(mid), lo = mid; else, hi = mid; end
end
gq_max = 10^lo;
end

function tf = local_valid(mphi, mx, gx, gq, g)
tf = true;
try
    axial_couplings(mphi, mx, gx, gq * g);
catch ME
    if strcmp(ME.identifier, 'axial_couplings:mediatorTooLight')
        tf = false;
    else
        rethrow(ME);
    end
end
end
