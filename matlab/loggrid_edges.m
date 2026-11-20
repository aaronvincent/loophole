function edges = loggrid_edges(L_lo, L_hi, L_extra)
% loggrid_edges  Sorted panel edges in log(s) for the fixed-grid thermal
% integrators (svavg_grid / ggxxrate_fast).
%
% At low T the thermal kernel concentrates the integrand in a thin boundary
% layer just above the annihilation threshold s = 4 mx^2 (= exp(L_lo)): the
% relevant width in log(s) shrinks like T/mx.  A single coarse panel there
% grossly under-resolves the peak, so we lay down a geometrically graded set
% of breakpoints rising from L_lo (spacing 1e-4 .. span in log s).  Extra
% breakpoints (meson thresholds, BW-pole sub-panels, NWA notch edges) are
% merged in.
if nargin < 3, L_extra = []; end
span = L_hi - L_lo;
d    = logspace(-4, log10(span), 14);     % distances above L_lo (log s)
thr  = L_lo + d(d > 0 & d < span);
edges = unique([L_lo, L_hi, L_extra(:).', thr]);
end
