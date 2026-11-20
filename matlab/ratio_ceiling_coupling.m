function gc = ratio_ceiling_coupling(rfun, g_lo, g_max)
%RATIO_CEILING_COUPLING  Coupling above which an exclusion ratio falls back below 1.
%
%   GC = RATIO_CEILING_COUPLING(RFUN, G_LO, G_MAX) returns the smallest coupling
%   GC > G_LO with RFUN(GC) = 1 and RFUN < 1 just above it -- the upper edge of
%   an excluded BAND whose lower edge G_LO comes from ratio_crit_coupling.
%   NaN if RFUN stays >= 1 up to G_MAX (excluded all the way up).
%
%   Only K+ -> pi+ phi has one: at large coupling the phi decays visibly inside
%   the detector, the invisible fraction drops to 0, and the ratio with it.
%   The scan steps in quarter decades from G_LO; a ratio that is 0 or not
%   finite counts as below 1.

if ~(isfinite(g_lo) && g_lo > 0)
    gc = NaN;  return;
end
above = @(L) local_above(rfun, L);

L_prev = log10(g_lo);  L = L_prev + 0.25;
while L <= log10(g_max)
    if ~above(L), break; end
    L_prev = L;  L = L + 0.25;
end
if L > log10(g_max), gc = NaN;  return; end

% Bracket [L_prev, L]: ratio >= 1 at L_prev, < 1 at L; r = 0 maps to a large
% negative log so fzero still sees the sign change.  If the very first step
% already fell below 1, L_prev is the lower root itself, where log r ~ 0 with
% either sign: look inside that quarter decade for a point clearly above 1.
f = @(x) log(max(rfun(10^x), 1e-300));
if L_prev == log10(g_lo)
    Ls = L_prev + (0.01:0.01:0.24);
    k  = find(arrayfun(@(x) f(x) > 0, Ls), 1, 'last');
    if isempty(k), gc = g_lo;  return; end    % band thinner than the grid
    L_prev = Ls(k);
end
gc = 10^fzero(f, [L_prev, L]);
end

% -------------------------------------------------------------------------
function tf = local_above(rfun, L)
r = rfun(10^L);
tf = isfinite(r) && r >= 1;
end
