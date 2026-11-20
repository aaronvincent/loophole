function E = neff_band_edges(nd_fn, le_lo, le_hi, n, mu, sigma, nsig, want)
% neff_band_edges  Coupling-band edges where a Delta-Neff component is
% excluded at the nsig-sigma level of the two-sided measured-Delta-Neff
% Gaussian (centre mu, width sigma) -- the SAME constraint shape as the
% GAMBIT-light neff_loglike.  A quantity is "excluded" where it lies more
% than thr = nsig*sigma from the measured central value mu.
%
%   E = neff_band_edges(nd_fn, le_lo, le_hi, n, mu, sigma, nsig [, want])
%
% nd_fn(le) returns the row [dNeff_chi, dNeff_phi, dNeff_total] (the total
% is SIGNED) evaluated at log10(coupling) = le.  A single grid of n points
% is swept and reused for every band, so nd_fn is called n times per mass
% point (plus a few bisection calls per edge).
%
% Bands (excluded where the listed signed distance exceeds thr):
%   E.chi     : dNeff_chi - mu        chi alone reaches the +nsig-sigma edge
%   E.phi     : dNeff_phi - mu        phi alone reaches it
%   E.tot     : |dNeff_total - mu|    total outside [mu-thr, mu+thr] (magnitude)
%   E.tot_pos : dNeff_total - mu      dark radiation pushes ABOVE mu+thr
%   E.tot_neg : mu - dNeff_total      photon heating pushes BELOW mu-thr
% chi and phi are positive dark-radiation contributions, so only the +side
% can bite; the two-sided window matters for the SIGNED total (a sub-MeV phi
% heats the photons -> dNeff_total<0 -> the tot_neg branch).
%
% Returns struct E with band edges [lo hi] in *linear* coupling (10^le),
% NaN for an edge not bracketed in the scanned range (open band -> caller
% caps it), and E.allon.<f> = true if excluded across the whole scan.
% [lo hi] is the ENVELOPE: lo = first up-crossing, hi = last down-crossing.
% It is only right when the excluded set is one interval.
%
% E.iv.<f> is the excluded set itself, one row [lo hi] per disjoint interval
% in ascending coupling (K x 2, linear coupling; NaN lo = open below the scan,
% NaN hi = open above it; 0 x 2 if nothing is excluded).  Needed since
% dNeff can change sign along the coupling (e.g. Primakoff: positive while phi
% is released as late free radiation, negative once phi stays coupled), which
% leaves an allowed strip between two excluded intervals.  That strip can be
% narrower than the grid step, so for 'tot' two neighbouring excluded grid
% points whose SIGNED excess (dNeff_total - mu) has opposite signs are searched
% (bisection on the signed value) for the allowed window the total must cross.
%
% want : optional cellstr selecting which of the five bands to compute
%        (any of {'chi','phi','tot','tot_pos','tot_neg'}).  Default: all.
%        Bands not requested are returned as [NaN NaN] / allon=false WITHOUT
%        bisection refinement -- the refinement makes off-grid nd_fn calls
%        that miss the grid cache, so skipping unused bands (e.g. chi/phi at
%        the CMB epoch, which no caller reads) avoids real Delta-Neff calls.
%
% Edges are bisection-refined on nd_fn: the excess is (near-)discontinuous
% at the thermalization onset, where linear interpolation would quantize
% the edge to the scan grid (the cause of choppy limits).

if nargin < 8 || isempty(want)
    want = {'chi','phi','tot','tot_pos','tot_neg'};
end

thr = nsig * sigma;                 % exclusion half-width about mu

le = linspace(le_lo, le_hi, n);
D  = nan(n, 3);
for k = 1:n
    try, D(k,:) = nd_fn(le(k)); catch, end
end

% Band specs: name, scalar selector on a single nd_fn row (signed distance
% from mu; > thr == excluded), and the same evaluated on the precomputed grid.
specs = { 'chi',     @(y)  y(1)-mu,       (D(:,1).'-mu);
          'phi',     @(y)  y(2)-mu,       (D(:,2).'-mu);
          'tot',     @(y)  abs(y(3)-mu),  abs(D(:,3).'-mu);
          'tot_pos', @(y)  y(3)-mu,       (D(:,3).'-mu);
          'tot_neg', @(y)  mu-y(3),       (mu-D(:,3).') };

E = struct();
E.allon = struct();
E.iv = struct();
for s = 1:size(specs,1)
    nm = specs{s,1};
    if ismember(nm, want)
        if strcmp(nm, 'tot')            % signed total: allowed strips inside a sign flip
            sg = {@(y) y(3)-mu, D(:,3).'-mu};
        else
            sg = {};
        end
        [E.(nm), E.allon.(nm), E.iv.(nm)] = find_band(nd_fn, le, specs{s,3}, specs{s,2}, thr, sg{:});
    else
        E.(nm) = [NaN, NaN];        % not requested -> no refinement
        E.allon.(nm) = false;
        E.iv.(nm) = zeros(0, 2);
    end
end
end

% =======================================================================
function [edges, allon, iv] = find_band(nd_fn, le, gridv, sval, thr, ssel, sgrid)
% Band where sval(.) > thr, located from sign changes of the excess
% (sval - thr) on the grid, each refined by bisection.  edges is the envelope
% (first up-crossing, last down-crossing); iv lists the disjoint intervals.
n    = numel(le);
e    = gridv - thr;                        % >0 == excluded
good = ~isnan(gridv);
iv   = zeros(0, 2);
if ~any(good)
    edges = [NaN, NaN]; allon = false; return;
end
pos   = e > 0;
allon = all(pos(good));

lo = NaN; hi = NaN;
kg = find(good);
inside = pos(kg(1));                       % excluded at the bottom of the scan
cur = NaN;                                 % lower edge of the open interval
for k = 1:n-1
    if ~good(k) || ~good(k+1), continue; end
    if e(k) <= 0 && e(k+1) > 0                 % up-crossing -> lower edge
        lc = refine_edge(nd_fn, sval, thr, le(k), le(k+1), e(k));
        if isnan(lo), lo = lc; end
        cur = lc;  inside = true;
    elseif e(k) > 0 && e(k+1) <= 0             % down-crossing -> upper edge
        hi = refine_edge(nd_fn, sval, thr, le(k), le(k+1), e(k));
        iv(end+1, :) = [edge_lin(cur), edge_lin(hi)]; %#ok<AGROW>
        cur = NaN;  inside = false;
    elseif nargin > 5 && e(k) > 0 && e(k+1) > 0 && sign(sgrid(k)) ~= sign(sgrid(k+1))
        % Excluded on both sides of a sign flip of the signed total: it passes
        % through the allowed window in between unless it jumps.  Find a point
        % in the window, then refine the down- and up-crossing around it.
        [lg, eg] = find_gap(nd_fn, ssel, thr, le(k), le(k+1), sgrid(k));
        if ~isnan(lg)
            hi = refine_edge(nd_fn, sval, thr, le(k), lg, e(k));
            iv(end+1, :) = [edge_lin(cur), edge_lin(hi)]; %#ok<AGROW>
            cur = refine_edge(nd_fn, sval, thr, lg, le(k+1), eg);
        end
    end
end
if inside                                  % still excluded at the top: open above
    iv(end+1, :) = [edge_lin(cur), NaN];
end
edges = [edge_lin(lo), edge_lin(hi)];
end

% =======================================================================
function lc = refine_edge(nd_fn, sval, thr, a, b, fa)
% Bisection on the excess (sval - thr) inside [a,b].
for it = 1:12
    m  = 0.5*(a + b);
    fm = sval(nd_fn(m)) - thr;
    if ~isfinite(fm), break; end
    if (fa <= 0) == (fm <= 0)
        a = m; fa = fm;
    else
        b = m;
    end
end
lc = 0.5*(a + b);
end

% =======================================================================
function [lg, eg] = find_gap(nd_fn, ssel, thr, a, b, sa)
% Bisection on the SIGNED excess ssel (opposite signs at a and b, |.| > thr at
% both) until a point with |ssel| <= thr is hit.  lg = that point (NaN if the
% sign flip is a jump with no allowed window), eg = |ssel(lg)| - thr <= 0.
lg = NaN; eg = NaN;
for it = 1:20
    m  = 0.5*(a + b);
    sm = ssel(nd_fn(m));
    if ~isfinite(sm), return; end
    if abs(sm) <= thr, lg = m; eg = abs(sm) - thr; return; end
    if sign(sm) == sign(sa), a = m; sa = sm; else, b = m; end
end
end

% =======================================================================
function v = edge_lin(le)
if isnan(le), v = NaN; else, v = 10^le; end
end
