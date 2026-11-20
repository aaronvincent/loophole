function [h, n_unplaced] = inline_labels(ax, varargin)
%INLINE_LABELS  Replace the legend with text labels drawn on the objects.
%
%   INLINE_LABELS(AX) collects every object in AX that would have claimed a
%   legend row (non-empty DisplayName) and labels it in place:
%
%     * Line objects            -> rotated text riding the curve (clabel-style)
%     * legend_swatch proxies   -> horizontal text at the deepest interior
%                                  point of the region the proxy stands for
%
%   The proxies are vertex-less by construction, so they cannot be labelled
%   directly; legend_swatch records the real region patch in the proxy's
%   UserData.geom and this function labels that instead.
%   A region patch with UserData.label_float = true is not labelled on its
%   edge: its label floats, horizontal, in the part of its fill farthest from
%   every other region, and is placed after everything else.
%
%   ALL geometry is done in NORMALIZED axes coordinates.  The axes here are
%   log-log and strongly anisotropic, so a slope in data units is not the slope
%   on the page -- rotating text by a data-space angle produces labels that
%   visibly do not follow their line.
%
%   Name/value options:
%     'FontSize'     label size in points                      (default 11)
%     'Grid'         occupancy raster [nx ny]                  (default [320 240])
%     'MinClear'     minimum clearance, in raster cells, for a
%                    placement to count as "on the object"     (default 6)
%     'MaxRot'       rotation clamp in degrees                 (default 80)
%     'Offsite'      allow free-space fallback beside the
%                    object when nothing fits on it            (default true)
%     'OverlapPen'   ranking penalty per OTHER region whose
%                    interior covers the candidate             (default 25)
%     'CrowdPen'     ranking penalty per OTHER line with ink
%                    within CrowdRad cells of the candidate    (default 20)
%     'CrowdRad'     radius, in cells, for the crowding count  (default 8)
%
%   The two penalties are what keeps a label out of the places where the plot
%   is ambiguous.  Clearance alone (distance to the NEAREST ink) cannot see
%   either problem: a point buried under three overlapping exclusion fills is
%   nowhere near a boundary, and a point in a narrow alley between two curves
%   is equidistant from both, so both score well while reading badly.  Penalties
%   are SOFT -- a crowded spot is still taken when the object offers nothing
%   better, which is routine for a band that lies entirely inside another.
%
%   Returns the handles of the text objects created, and the number of labels
%   that could NOT be placed.  Nothing is ever dropped silently: every item
%   that is placed has its legend entry switched off, so a caller that calls
%   legend('show') when n_unplaced > 0 gets a compact legend holding exactly
%   the leftovers.
%
%   Call AFTER the figure has been given its final size: text extent is in
%   points, so the fraction of the axes a label covers depends on figure size.

p = inputParser;
p.addParameter('FontSize', 11);
p.addParameter('Grid', [320 240]);
p.addParameter('MinClear', 6);
p.addParameter('MaxRot', 80);
p.addParameter('Offsite', true);
p.addParameter('OverlapPen', 25);
p.addParameter('CrowdPen', 20);
p.addParameter('CrowdRad', 8);
p.addParameter('CandScan', 1500);   % candidates examined per label; must exceed
                                    % the boundary list or the interior fallback
                                    % rows at its tail are never reached
p.addParameter('FloatAlign', 'horizontal');  % floating region labels:
                                    % 'horizontal', or 'bottom' = parallel to
                                    % the fill's lower edge, just inside it
p.addParameter('Debug', false);
p.parse(varargin{:});
o = p.Results;

if nargin < 1 || isempty(ax), ax = gca; end
CAND_SCAN = o.CandScan;
nx = o.Grid(1); ny = o.Grid(2);

% Normalized coordinates are [0,1] on BOTH axes, but the axes rectangle is not
% square on the page, so an angle taken as atan2(dy,dx) in normalized units is
% not the angle the text must be rotated by.  Convert through the axes' own
% width:height in points.
u_ = ax.Units; ax.Units = 'points'; apos_ = ax.Position; ax.Units = u_;
aspect = apos_(3) / max(eps, apos_(4));

% ---------------------------------------------------------------- collect
% allchild, NOT ax.Children: every region patch, dashed continuation and the
% corner caption is drawn with HandleVisibility 'off', and .Children omits
% exactly those.  Using .Children here leaves the ink map almost empty and the
% swatch proxies with no geometry attached.
kids = flipud(allchild(ax));
items = struct('h', {}, 'lbl', {}, 'kind', {}, 'geom', {}, 'col', {});
for k = 1:numel(kids)
    ch = kids(k);
    if isempty(get_disp(ch)), continue; end
    if isa(ch, 'matlab.graphics.chart.primitive.Line')
        items(end+1) = struct('h', ch, 'lbl', ch.DisplayName, 'kind', 'line', ...
                              'geom', ch, 'col', ch.Color);                     %#ok<AGROW>
    elseif isa(ch, 'matlab.graphics.primitive.Patch')
        g = [];
        if isstruct(ch.UserData) && isfield(ch.UserData, 'geom') && ...
                ~isempty(ch.UserData.geom) && isgraphics(ch.UserData.geom)
            g = ch.UserData.geom;
        end
        if isempty(g), continue; end   % vertex-less proxy with no region: skip
        items(end+1) = struct('h', ch, 'lbl', ch.DisplayName, 'kind', 'patch', ...
                              'geom', g, 'col', ch.FaceColor);                  %#ok<AGROW>
    end
end
h = gobjects(0);
n_unplaced = 0;
if isempty(items), return; end
was_placed = false(1, numel(items));

% ------------------------------------------------------- normalized geometry
% Raster cell centres in [0,1]^2.
xc = (0.5:nx-0.5) / nx;
yc = (0.5:ny-0.5) / ny;
[XG, YG] = meshgrid(xc, yc);

% Per-item ink mask + the polyline (normalized) that a line label rides.
masks = false(ny, nx, numel(items));
polys = cell(1, numel(items));
for k = 1:numel(items)
    it = items(k);
    [xn, yn] = norm_xy(it.geom.XData, it.geom.YData, ax);
    good = isfinite(xn) & isfinite(yn);
    xn = xn(good); yn = yn(good);
    if isempty(xn), continue; end
    if strcmp(it.kind, 'line')
        [xn, yn] = densify(xn, yn, 4/max(nx,ny));
        polys{k} = [xn(:), yn(:)];
        masks(:,:,k) = stamp_points(xn, yn, nx, ny);
    else
        % Region: interior mask, plus its boundary as ink.
        masks(:,:,k) = inpolygon(XG, YG, xn, yn);
        % Densified boundary, same as a line: candidate placement measures the
        % tangent over an arclength window, and a polygon given by a handful of
        % far-apart vertices has no points inside that window.
        [xb, yb] = densify(xn, yn, 4/max(nx,ny));
        polys{k} = [xb(:), yb(:)];
    end
end

% Everything already drawn that a label must avoid: all ink, plus any text
% already in the axes (the corner caption sits bottom-right and is the most
% likely collision), plus a frame margin so labels never touch the box.
ink = false(ny, nx);
for k = 1:numel(items)
    if strcmp(items(k).kind, 'line')
        ink = ink | masks(:,:,k);
    else
        ink = ink | edge_of(masks(:,:,k));
    end
end
% Unlabelled lines (dashed continuations, band edges) are ink too.
for k = 1:numel(kids)
    ch = kids(k);
    if isa(ch, 'matlab.graphics.chart.primitive.Line') && isempty(get_disp(ch))
        [xn, yn] = norm_xy(ch.XData, ch.YData, ax);
        g = isfinite(xn) & isfinite(yn);
        if any(g)
            [xd, yd] = densify(xn(g), yn(g), 4/max(nx,ny));
            ink = ink | stamp_points(xd, yd, nx, ny);
        end
    end
end
% ------------------------------------------------------- penalty rasters
% cover(j,i): how many labelled REGION INTERIORS contain the cell.  A label is
% ranked down once per region other than its own, so "overabundant" prefers the
% part of its fill that no other exclusion covers.  Interiors, not edges: the
% edges are already in `ink`.
cover = zeros(ny, nx);
for k = 1:numel(items)
    if strcmp(items(k).kind, 'patch'), cover = cover + double(masks(:,:,k)); end
end

% line_masks: one dilated mask per LINE, labelled or not.  Counting them at a
% candidate gives the local line density, which is what "don't sit where the
% curves bunch up" needs -- the distance transform only ever sees the closest
% one.  Dilation is by CrowdRad cells, roughly a label height, so the count is
% "how many curves would this label have to share its box with".
line_masks = false(ny, nx, 0);
line_of_item = zeros(1, numel(items));   % item -> its slice, 0 for regions
for k = 1:numel(items)
    if ~strcmp(items(k).kind, 'line'), continue; end
    line_masks(:,:,end+1) = masks(:,:,k);   %#ok<AGROW>
    line_of_item(k) = size(line_masks, 3);
end
for k = 1:numel(kids)
    ch = kids(k);
    if ~isa(ch, 'matlab.graphics.chart.primitive.Line') || ~isempty(get_disp(ch))
        continue;
    end
    [xn, yn] = norm_xy(ch.XData, ch.YData, ax);
    g = isfinite(xn) & isfinite(yn);
    if ~any(g), continue; end
    [xd, yd] = densify(xn(g), yn(g), 4/max(nx,ny));
    line_masks(:,:,end+1) = stamp_points(xd, yd, nx, ny);   %#ok<AGROW>
end
near_line = false(ny, nx, size(line_masks, 3));
for k = 1:size(line_masks, 3)
    near_line(:,:,k) = dt(line_masks(:,:,k)) <= o.CrowdRad;
end
crowd_all = sum(double(near_line), 3);
% Ink of the lines that carry a label of their own.  A region whose boundary IS
% such a line (the overabundant fill is bounded by the relic locus, and every
% exclusion band by its own limit curve) must not park its label on that
% stretch: the line's own label then has nowhere left to sit, and the two names
% end up on top of each other.  Regions are pushed off it hard -- deeper into
% the fill, or along to a stretch that is theirs alone.
labelled_ink = false(ny, nx);
for k = 1:numel(items)
    if line_of_item(k) > 0
        labelled_ink = labelled_ink | line_masks(:,:,line_of_item(k));
    end
end
dt_labelled = dt(labelled_ink);   % cells to the nearest such line

% Two occupancy layers.  HARD is what a label may never touch: the frame, the
% corner caption, and any label already placed.  SOFT is drawn ink -- curves and
% region boundaries.  Without an opaque background behind the text, crossing a
% faint dash-dot boundary is legible; crossing another label never is.  So a
% label tries hard+soft first and falls back to hard-only before giving up and
% taking a leader line.
% 5% gutter, not 2%: the axes are drawn with Layer='top', so the frame line and
% grid paint OVER any text, and a label parked a hair inside the frame reads as
% clipped even though its box technically fits.
hard = frame_margin(nx, ny, 0.05);
ntxt = 0;
for k = 1:numel(kids)
    if isa(kids(k), 'matlab.graphics.primitive.Text') && strcmp(kids(k).Visible, 'on')
        hard = hard | stamp_extent(kids(k), ax, nx, ny);
        ntxt = ntxt + 1;
        if o.Debug
            old = kids(k).Units; kids(k).Units = 'normalized';
            fprintf('[il] pre-existing text stamped: Extent=[%.3f %.3f %.3f %.3f]\n', kids(k).Extent);
            kids(k).Units = old;
        end
    end
end
if o.Debug
    fprintf('[il] %d labelled items, %d pre-existing text stamped, ink=%d cells\n', ...
            numel(items), ntxt, nnz(ink));
end

% ------------------------------------------------------------ ordering
% Least flexible first: an object with few viable spots must claim one before a
% flexible neighbour takes it.  Greedy with this ordering beats a global
% optimiser here and is far easier to debug.
% Is anything of this object actually drawn inside the axes?  An object with
% no visible geometry (all-NaN data, or entirely outside the limits) gets no
% label AND no legend row: a legend entry for a curve that is not on the plot
% is worse than saying nothing.  Regions are tested on their interior raster,
% not their vertices -- a band can span the whole panel with every vertex off
% the edge.
visible = false(1, numel(items));
for k = 1:numel(items)
    if isempty(polys{k}), continue; end
    if strcmp(items(k).kind, 'patch')
        visible(k) = any(any(masks(:,:,k)));
    else
        px_ = polys{k}(:,1); py_ = polys{k}(:,2);
        visible(k) = any(px_ >= 0 & px_ <= 1 & py_ >= 0 & py_ <= 1);
    end
end

flex = zeros(1, numel(items));
for k = 1:numel(items)
    % polys{k} is [] when the object's data is entirely NaN -- routine for
    % mphi = 1e-9, where whole curves are invalid.  Indexing (:,1) on a 0x0
    % errors ("Index in position 2 exceeds array bounds"), so the empty case is
    % skipped here exactly as it already is in the placement loop below.
    if isempty(polys{k}), flex(k) = 0; continue; end
    own  = own_ink(items(k), masks(:,:,k));
    cl   = dt((hard | ink) & ~own);
    cand = stamp_points(polys{k}(:,1)', polys{k}(:,2)', nx, ny);
    flex(k) = nnz(cand & cl > o.MinClear);
end
% A floating label goes FIRST: the part of its fill that no other region
% covers is usually small, and a curve label (which can take either side of
% its line, anywhere along it) would otherwise claim it.
floats = false(1, numel(items));
for k = 1:numel(items)
    ud = items(k).geom.UserData;
    floats(k) = strcmp(items(k).kind, 'patch') && isstruct(ud) && ...
                isfield(ud, 'label_float') && ud.label_float;
end
flex(floats) = -Inf;
[~, order] = sort(flex, 'ascend');

% ------------------------------------------------------------ placement
for k = order
    it   = items(k);
    poly = polys{k};
    if isempty(poly) || ~visible(k), continue; end

    % Clearance from everything except this object's OWN ink.  For a line that
    % is the whole polyline (a point on the curve is itself ink); for a region
    % it is the boundary ONLY -- subtracting the filled interior would erase
    % every competing curve and already-placed label inside it, which is what
    % drove every region label out to the frame.
    own     = own_ink(it, masks(:,:,k));
    clear_k = dt((hard | ink) & ~own);

    % Ranking penalty for THIS item: other regions covering the cell, and other
    % lines crowding it.  Its own fill and its own curve are subtracted -- a
    % label is supposed to sit on its own object.
    cov_k = cover;
    if strcmp(it.kind, 'patch'), cov_k = cov_k - double(masks(:,:,k)); end
    crowd_k = crowd_all;
    if line_of_item(k) > 0
        crowd_k = crowd_k - double(near_line(:,:,line_of_item(k)));
    end
    pen_k = o.OverlapPen * max(0, cov_k) + o.CrowdPen * max(0, crowd_k);
    % A floating label scans its WHOLE fill and never leaves it: in a thin or
    % crowded fill every top-ranked spot can be taken, and the off-object
    % fallback would then park it outside the region it names.
    scan_k = CAND_SCAN;
    if floats(k), scan_k = Inf; end


    % Measure the label before choosing a spot: the required clearance and the
    % box-fit test both depend on its size, which depends on figure size.
    % Try the label at full size, then progressively smaller, before ever
    % accepting an overlap.  A label one point down is far less damaging than
    % two labels on top of each other -- which is what "always place" degrades
    % into on a crowded panel.
    placed = false;  fs_used = o.FontSize;  lw = 0; lh = 0;
    % A floating label runs the sizes twice: first offering only spots whose
    % box lies wholly in its own fill and in no other region, then anywhere.
    % Being inside its own region matters more than its size or than clearing
    % a faint boundary.
    % No label ever shrinks: every inline label is drawn at FontSize, so the
    % text reads as one size across a panel and a row.  (It used to step down
    % to 0.88x and 0.78x to squeeze onto crowded objects, which left a single
    % panel with three label sizes; the drivers now set FONT.inline to the old
    % smallest step instead.)  A floating label still runs its one size twice:
    % first offering only spots wholly in its own fill, then anywhere.
    shrinks = 1;
    if floats(k), shrinks = [1, 1]; end
    for si_ = 1:numel(shrinks)
        shrink_ = shrinks(si_);
        fs_ = o.FontSize * shrink_;
        [lw, lh] = label_size(ax, it.lbl, fs_);
        % Recorded here, not on success: the fallbacks below reuse LW/LH from
        % the last size TRIED, so they must draw at that size or the text
        % overflows the box they were fitted against.
        fs_used = fs_;
        % The keep-off band is a label height wide, which is what it takes to
        % leave the line's own label room to sit; CrowdRad is a fraction of
        % that, so the plain crowd penalty never fires at the depth a region
        % label actually sits.
        pen_use = pen_k;
        if strcmp(it.kind, 'patch')
            pen_use = pen_use + 3*o.OverlapPen * double(dt_labelled <= 1.3*lh*ny);
        end
        if floats(k) && strcmpi(o.FloatAlign, 'bottom')
            cands = float_bottom_candidates(masks(:,:,k), cov_k > 0, lw, lh, nx, ny, ...
                                            aspect, o.MaxRot, si_ == 1);
        elseif floats(k)
            cands = float_candidates(masks(:,:,k), cov_k > 0, lw, lh, nx, ny, si_ == 1);
        else
            cands = ranked_candidates(it, poly, clear_k, pen_use, masks(:,:,k), ...
                                      nx, ny, o.MaxRot, aspect, lw, lh);
        end
        if o.Debug
            fprintf('[il] cand %-30s n=%d fs=%.1f lw=%.3f lh=%.3f\n', ...
                    it.lbl, size(cands,1), fs_, lw, lh);
        end
        for pass_ = 1:2
            % The object's OWN ink is never part of the avoid set.  The fitted
            % box is the text extent padded 15%/30% (label_size), so a box that
            % must clear its own curve cannot be seated closer than half that
            % padded height -- which is a visible gap of most of a line of text,
            % and no offset value can shrink it.  Keeping the label off its own
            % curve is the offset's job (and the bow test's, which rejects a
            % stretch where the curve wanders across the text); the box test is
            % then free to seat it snugly.  Pass 2 additionally tolerates
            % crossing ANOTHER object's faint ink.
            if pass_ == 1, avoid = hard | (ink & ~own); else, avoid = hard; end
            for c = 1:min(size(cands, 1), scan_k)
                cx = cands(c,1); cy = cands(c,2); rot = cands(c,3);
                % REJECT a candidate whose box leaves the axes; do not clamp it
                % back in.  Region patches are drawn far above the visible ylim
                % (fills run to ytop = ylim*1e3), so a top band's "inward"
                % normal points off-screen; clamping slides the label off the
                % boundary it was chosen to follow and leaves it clipped.
                if ~box_inside(cx, cy, lw, lh, rot), continue; end
                if bbox_free(avoid, cx, cy, lw, lh, rot, nx, ny)
                    h(end+1) = place(ax, cx, cy, it, rot, fs_); %#ok<AGROW>
                    hard = hard | bbox_mask(cx, cy, lw, lh, rot, nx, ny);
                    placed = true;
                    if o.Debug
                        fprintf('[il] p%d %-30s c=(%.3f,%.3f) rot=%.0f fs=%.1f\n', ...
                                pass_, it.lbl, cx, cy, rot, fs_);
                    end
                    break;
                end
            end
            if placed, break; end
        end
        if placed, break; end
    end

    if ~placed && o.Offsite && ~floats(k)
        % Nothing fits anywhere ON the object: park it in the nearest free
        % space instead.  No leader line is drawn -- an arrow trailing across
        % the panel costs more ink than it buys, and a label sitting just off
        % its curve at the curve's own angle already reads as belonging to it.
        % The angle is therefore the local tangent at the nearest point on the
        % object, NOT horizontal: a horizontal label beside a steep boundary is
        % the one placement that reads as belonging to nothing.
        % Distance to the object's own INK -- its boundary for a region, the
        % curve itself for a line.  Not to the filled interior: that is zero
        % over the whole fill, so the scan would start at whichever interior
        % cell happens to come first and park a band's label in the middle of
        % the panel, far from the edge that gives the band its meaning.
        near = dt(own);
        [~, ord] = sort(near(:), 'ascend');
        % Subsampled for the nearest-point search only: the densified polyline
        % has thousands of points and the exact nearest vertex is not needed --
        % only the stretch of curve whose angle the label should take.
        sub = unique(round(linspace(1, size(poly,1), min(size(poly,1), 200))));
        rsub = arrayfun(@(ii) local_angle(poly, ii, lw, aspect, o.MaxRot), sub);
        % Nearest subsample point for EVERY cell at once, then its angle: the
        % scan below visits thousands of cells and doing this per cell is the
        % slowest thing in the function.
        best_d = inf(ny, nx); ROTMAP = zeros(ny, nx);
        for q = 1:numel(sub)
            dq = (XG - poly(sub(q),1)).^2 + (YG - poly(sub(q),2)).^2;
            upd = dq < best_d;
            best_d(upd) = dq(upd); ROTMAP(upd) = rsub(q);
        end
        % Clear of the curves first, clear of the other LABELS at worst: text
        % over text is illegible, text over a faint dash-dot boundary is not.
        % Never over its OWN curve, though: the scan starts at the cells
        % nearest it, so without that the fallback strikes the label through.
        for avoid = {hard | ink, hard | own}
            for c = 1:min(numel(ord), 20000)
                [j, i] = ind2sub([ny nx], ord(c));
                cx = (i-0.5)/nx; cy = (j-0.5)/ny;
                rot = ROTMAP(j, i);
                [cx, cy] = clamp_inside(cx, cy, lw, lh, rot);
                if ~box_inside(cx, cy, lw, lh, rot), continue; end
                if ~bbox_free(avoid{1}, cx, cy, lw, lh, rot, nx, ny), continue; end
                h(end+1) = place(ax, cx, cy, it, rot, fs_used); %#ok<AGROW>
                if o.Debug
                    fprintf('[il] OFFSITE   %-34s c=(%.3f,%.3f) rot=%.0f\n', ...
                            it.lbl, cx, cy, rot);
                end
                hard = hard | bbox_mask(cx, cy, lw, lh, rot, nx, ny);
                placed = true;
                break;
            end
            if placed, break; end
        end
    end

    % Last resort: nothing fits cleanly ON the object and there is no clear
    % space beside it either.  Take the position with the LEAST overlap and
    % place it anyway -- a slightly crowded inline label still beats a legend
    % box, which is the thing inline labels exist to remove.  Tried AFTER the
    % off-object fallback, because a label parked in clear space a little away
    % from its curve reads better than one printed over another label.
    if ~placed && ~isempty(cands)
        best_ = Inf; bx_ = []; by_ = []; br_ = 0;
        for c = 1:min(size(cands, 1), scan_k)
            cx = cands(c,1); cy = cands(c,2); rot = cands(c,3);
            if ~box_inside(cx, cy, lw, lh, rot), continue; end
            ov = nnz(hard & bbox_mask(cx, cy, lw, lh, rot, nx, ny));
            if ov < best_, best_ = ov; bx_ = cx; by_ = cy; br_ = rot; end
        end
        if ~isempty(bx_)
            h(end+1) = place(ax, bx_, by_, it, br_, fs_used); %#ok<AGROW>
            hard = hard | bbox_mask(bx_, by_, lw, lh, br_, nx, ny);
            placed = true;
            if o.Debug
                fprintf('[il] CROWDED %-30s c=(%.3f,%.3f) overlap=%d\n', ...
                        it.lbl, bx_, by_, best_);
            end
        end
    end

    was_placed(k) = placed;
end

% Anything that found no room keeps its legend row; everything placed loses
% one.  A constraint must never end up both unlabelled AND absent from the
% legend, which is what silently skipping an unplaceable label would do.
for k = 1:numel(items)
    if was_placed(k) || ~visible(k), items(k).h.HandleVisibility = 'off'; end
end
% Only DRAWN-but-unplaceable entries justify a leftover legend.
n_unplaced = nnz(~was_placed & visible);
end

% ===================================================================== helpers

function s = get_disp(ch)
s = '';
if isprop(ch, 'DisplayName'), s = ch.DisplayName; end
end

function [xn, yn] = norm_xy(x, y, ax)
xn = norm_ax(x, xlim(ax), strcmp(ax.XScale, 'log'));
yn = norm_ax(y, ylim(ax), strcmp(ax.YScale, 'log'));
end

function n = norm_ax(v, lim, islog)
v = double(v(:))';
if islog
    v(v <= 0) = NaN;
    n = (log10(v) - log10(lim(1))) / (log10(lim(2)) - log10(lim(1)));
else
    n = (v - lim(1)) / (lim(2) - lim(1));
end
end

function [xo, yo] = densify(x, y, step)
% Resample a polyline so consecutive points are <= step apart (normalized).
xo = x(1); yo = y(1);
for i = 2:numel(x)
    d = hypot(x(i)-x(i-1), y(i)-y(i-1));
    n = max(1, ceil(d/step));
    xo = [xo, x(i-1) + (1:n)*(x(i)-x(i-1))/n]; %#ok<AGROW>
    yo = [yo, y(i-1) + (1:n)*(y(i)-y(i-1))/n]; %#ok<AGROW>
end
end

function m = stamp_points(xn, yn, nx, ny)
m = false(ny, nx);
i = round(xn*nx + 0.5); j = round(yn*ny + 0.5);
ok = i >= 1 & i <= nx & j >= 1 & j <= ny;
m(sub2ind([ny nx], j(ok), i(ok))) = true;
end

function e = edge_of(mask)
e = mask & ~(imshift(mask,1,0) & imshift(mask,-1,0) & imshift(mask,0,1) & imshift(mask,0,-1));
end

function s = imshift(m, dj, di)
s = false(size(m));
[ny, nx] = size(m);
js = max(1,1+dj):min(ny,ny+dj);
is = max(1,1+di):min(nx,nx+di);
s(js, is) = m(js-dj, is-di);
end

function m = frame_margin(nx, ny, frac)
m = false(ny, nx);
bx = max(1, round(frac*nx)); by = max(1, round(frac*ny));
m(1:by, :) = true; m(end-by+1:end, :) = true;
m(:, 1:bx) = true; m(:, end-bx+1:end) = true;
end

function m = stamp_extent(t, ~, nx, ny)
% Bounding box of a text object in normalized axes units.  For rotated text
% the axis-aligned box of the rotated rectangle is used (conservative).
old = t.Units; t.Units = 'normalized';
e = t.Extent; r = deg2rad(t.Rotation);
t.Units = old;
w = e(3); hgt = e(4);
if r ~= 0
    w2   = abs(w*cos(r)) + abs(hgt*sin(r));
    hgt2 = abs(w*sin(r)) + abs(hgt*cos(r));
    cx = e(1) + e(3)/2; cy = e(2) + e(4)/2;
    e = [cx - w2/2, cy - hgt2/2, w2, hgt2];
end
m = false(ny, nx);
% Same LaTeX under-measurement applies to text already in the axes (the corner
% caption), so its stamped box is padded too.
cx0 = e(1) + e(3)/2; cy0 = e(2) + e(4)/2;
e = [cx0 - e(3)*1.15/2, cy0 - e(4)*1.30/2, e(3)*1.15, e(4)*1.30];
i0 = max(1, floor(e(1)*nx)); i1 = min(nx, ceil((e(1)+e(3))*nx));
j0 = max(1, floor(e(2)*ny)); j1 = min(ny, ceil((e(2)+e(4))*ny));
if i1 >= i0 && j1 >= j0, m(j0:j1, i0:i1) = true; end
end

function D = dt(occ)
% Distance (in cells) from every cell to the nearest occupied cell.
if exist('bwdist', 'file') == 2 || exist('bwdist', 'builtin') == 5
    D = double(bwdist(occ));
    return;
end
[ny, nx] = size(occ);
D = inf(ny, nx); D(occ) = 0;
a = 1; b = sqrt(2);
for j = 1:ny
    for i = 1:nx
        v = D(j,i);
        if j > 1
            v = min(v, D(j-1,i)+a);
            if i > 1,  v = min(v, D(j-1,i-1)+b); end
            if i < nx, v = min(v, D(j-1,i+1)+b); end
        end
        if i > 1, v = min(v, D(j,i-1)+a); end
        D(j,i) = v;
    end
end
for j = ny:-1:1
    for i = nx:-1:1
        v = D(j,i);
        if j < ny
            v = min(v, D(j+1,i)+a);
            if i > 1,  v = min(v, D(j+1,i-1)+b); end
            if i < nx, v = min(v, D(j+1,i+1)+b); end
        end
        if i < nx, v = min(v, D(j,i+1)+a); end
        D(j,i) = v;
    end
end
end

function c = float_candidates(own, other, lw, lh, nx, ny, pure_only)
% Candidates for a region whose label FLOATS in its fill rather than riding an
% edge (patch UserData.label_float): horizontal, at every cell of the fill,
% best first.  Ranked by how many box cells land in another region or outside
% the fill, then by distance from the nearest other region -- so the label
% takes the part of the fill that is most clearly this region alone.
bw = ceil(lw*nx/2); bh = ceil(lh*ny/2);
[J, I] = find(own);
c = zeros(0,3);
if isempty(J), return; end
j0 = max(1, J-bh); j1 = min(ny, J+bh);
i0 = max(1, I-bw); i1 = min(nx, I+bw);
n_oth = box_count(other, j0, j1, i0, i1);
n_out = box_count(~own,  j0, j1, i0, i1);
D = dt(other);
score = D(sub2ind([ny nx], J, I)) - 1e3*(n_oth + n_out);
if pure_only, score(n_oth + n_out > 0) = -Inf; end
[~, ord] = sort(score, 'descend');
ord = ord(isfinite(score(ord)));
c = [(I(ord)-0.5)/nx, (J(ord)-0.5)/ny, zeros(numel(ord), 1)];
end

function c = float_bottom_candidates(own, other, lw, lh, nx, ny, aspect, maxrot, pure_only)
% Floating-label candidates parallel to the fill's LOWER edge, seated just
% inside it (FloatAlign = 'bottom').  One candidate per raster column: the
% edge is the lowest fill cell, its slope a least-squares fit over one label
% width, and the label centre sits half a label height (plus air) above the
% fitted edge along its normal.  Ranked by how much of the rotated box leaves
% the fill or enters another region, then by how straight the edge is under
% it, then by distance from other regions.
c = zeros(0,3);
jb = nan(1, nx);
for i = 1:nx
    j = find(own(:, i), 1, 'first');
    if ~isempty(j), jb(i) = j; end
end
bw = max(2, round(lw*nx/2));
D  = dt(other);
[uu, vv] = meshgrid(linspace(-0.5, 0.5, 11), linspace(-0.5, 0.5, 4));
uu = uu(:); vv = vv(:);
rows = zeros(0, 5);
for i = 1:nx
    ii = max(1, i-bw):min(nx, i+bw);
    ok = isfinite(jb(ii));
    if ~isfinite(jb(i)) || nnz(ok) < 0.8*numel(ii), continue; end
    xs = (ii(ok) - 0.5)/nx;  ys = (jb(ii(ok)) - 1)/ny;   % edge = cell bottom
    P  = polyfit(xs, ys, 1);
    res = sqrt(mean((polyval(P, xs) - ys).^2)) * ny;    % edge wiggle, in cells
    rot = atan2d(P(1), aspect);                          % page angle
    if abs(rot) > maxrot, continue; end
    x0 = (i - 0.5)/nx;  y0 = polyval(P, x0);
    d  = 0.65*lh;                                        % centre above the edge
    cx = x0 - d*sind(rot)/aspect;  cy = y0 + d*cosd(rot);
    % Rotated-box sample points, in normalized units (width lw along the edge).
    px = cx + (uu*lw*cosd(rot) - vv*lh*sind(rot)/aspect);
    py = cy + (uu*lw*sind(rot)*aspect + vv*lh*cosd(rot));
    pi_ = min(nx, max(1, ceil(px*nx)));  pj_ = min(ny, max(1, ceil(py*ny)));
    lin = sub2ind([ny nx], pj_, pi_);
    bad = nnz(~own(lin) | other(lin) | px < 0 | px > 1 | py < 0 | py > 1);
    pj0 = min(ny, max(1, ceil(cy*ny)));  pi0 = min(nx, max(1, ceil(cx*nx)));
    rows(end+1, :) = [cx, cy, rot, bad, D(pj0, pi0) - 5*res]; %#ok<AGROW>
end
if isempty(rows), return; end
if pure_only, rows = rows(rows(:,4) == 0, :); end
[~, ord] = sortrows([rows(:,4), -rows(:,5)]);
c = rows(ord, 1:3);
end

function n = box_count(m, j0, j1, i0, i1)
% Set cells of M inside each box [j0:j1, i0:i1], via a summed-area table.
[ny, nx] = size(m);
P = zeros(ny+1, nx+1);
P(2:end, 2:end) = cumsum(cumsum(double(m), 1), 2);
sz = [ny+1, nx+1];
n = P(sub2ind(sz, j1+1, i1+1)) - P(sub2ind(sz, j0, i1+1)) ...
  - P(sub2ind(sz, j1+1, i0))   + P(sub2ind(sz, j0, i0));
end

function m = own_ink(it, mask)
if strcmp(it.kind, 'line'), m = mask; else, m = edge_of(mask); end
end

function [w, hgt] = label_size(ax, str, fs)
% MATLAB's Extent for LaTeX-interpreted text under-reports the rendered width
% (the reported box is measured before the interpreter lays the glyphs out).
% Placement that trusts it produces labels that pass the fit test and still
% visibly collide, so the measured box is padded before it is used for either
% the fit test or the occupancy stamp.
PAD_W = 1.15; PAD_H = 1.30;
t = text(ax, 0.5, 0.5, str, 'Units', 'normalized', 'Interpreter', 'latex', ...
         'FontSize', fs, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'Visible', 'off', ...
         'HandleVisibility', 'off');
e = t.Extent; w = e(3)*PAD_W; hgt = e(4)*PAD_H;
delete(t);
end

function t = place(ax, cx, cy, it, rot, fs)
t = text(ax, cx, cy, it.lbl, 'Units', 'normalized', 'Interpreter', 'latex', ...
         'FontSize', fs, 'Color', clamp_col(it.col), 'Rotation', rot, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'Clipping', 'on', 'HandleVisibility', 'off');
end

function [w2, h2] = rot_box(w, hgt, rot)
r = deg2rad(rot);
w2   = abs(w*cos(r)) + abs(hgt*sin(r));
h2   = abs(w*sin(r)) + abs(hgt*cos(r));
end

function tf = box_inside(cx, cy, w, hgt, rot)
[w2, h2] = rot_box(w, hgt, rot);
tf = (cx - w2/2 >= 0.005) && (cx + w2/2 <= 0.995) && ...
     (cy - h2/2 >= 0.005) && (cy + h2/2 <= 0.995);
end

function [cx, cy] = clamp_inside(cx, cy, w, hgt, rot)
[w2, h2] = rot_box(w, hgt, rot);
cx = min(max(cx, w2/2 + 0.005), 1 - w2/2 - 0.005);
cy = min(max(cy, h2/2 + 0.005), 1 - h2/2 - 0.005);
end

function m = bbox_mask(cx, cy, w, hgt, rot, nx, ny)
m = false(ny, nx);
[i0, i1, j0, j1] = bbox_range(cx, cy, w, hgt, rot, nx, ny);
if i1 >= i0 && j1 >= j0, m(j0:j1, i0:i1) = true; end
end

function tf = bbox_free(occ, cx, cy, w, hgt, rot, nx, ny)
% Exact test: does the label's box overlap anything already there?  A single
% clearance radius approximates the box by a disc, which is far too strict for
% a wide label and too loose for a tall one.
% Indexed directly rather than through bbox_mask: the fallback scans tens of
% thousands of positions and a full-raster temporary per position is the whole
% cost of that scan.
[i0, i1, j0, j1] = bbox_range(cx, cy, w, hgt, rot, nx, ny);
tf = ~(i1 >= i0 && j1 >= j0 && any(any(occ(j0:j1, i0:i1))));
end

function [i0, i1, j0, j1] = bbox_range(cx, cy, w, hgt, rot, nx, ny)
[w2, h2] = rot_box(w, hgt, rot);
i0 = max(1, floor((cx - w2/2)*nx)); i1 = min(nx, ceil((cx + w2/2)*nx));
j0 = max(1, floor((cy - h2/2)*ny)); j1 = min(ny, ceil((cy + h2/2)*ny));
end

function c = ranked_candidates(it, poly, clearf, pen, own_mask, nx, ny, maxrot, aspect, lw, lh)
% Candidate label positions, best first, as [x y rot] rows.
%
% Lines and regions are handled by the SAME curve-following logic: a region's
% label rides its boundary exactly as a line label rides the line.  The only
% difference is that a region label is then pushed off the boundary into the
% interior, so it sits just inside the edge rather than straddling it.
%
% A region samples SEVERAL inward depths along that boundary, so a crowded edge
% can be answered by stepping deeper into the fill without the label ever
% leaving its own boundary or changing its angle.
%
% PEN is the per-cell ranking penalty (other regions covering it, other lines
% crowding it).  It is subtracted from clearance, so it trades off against room
% rather than vetoing anything.
x = poly(:,1); y = poly(:,2); n = numel(x);
c = zeros(0,3);
if n < 5, return; end
is_region = strcmp(it.kind, 'patch');

% The tangent window is ONE LABEL WIDE, measured along the curve -- not a fixed
% fraction of the point count.  "Parallel" only means anything over the stretch
% the glyphs actually cover: a secant taken over a tenth of the whole curve
% tilts the text against the piece of line it sits next to, which is exactly the
% mismatch that reads as a label belonging to nothing.
sarc = [0; cumsum(hypot(diff(x), diff(y)))];
half = max(lw/2, 1e-3);
% Window ends for every point at once: the polylines here are densified to a
% few thousand points, and searching the arclength vector per point (per relax
% pass) is the one O(n^2) step in this file.
[su, ifirst] = unique(sarc, 'first');
[~,  ilast]  = unique(sarc, 'last');
if numel(su) < 2
    IA = ones(n,1); IB = n*ones(n,1);
else
    tlo = min(max(sarc - half, su(1)), su(end));
    thi = min(max(sarc + half, su(1)), su(end));
    IA = round(interp1(su, ifirst, tlo, 'next',     'extrap'));
    IB = round(interp1(su, ilast,  thi, 'previous', 'extrap'));
    IA = min(max(IA, 1), n); IB = min(max(IB, 1), n);
end

% A ladder of relaxations, strictest first.  The rotation clamp and the bow
% test both exist to keep a label parallel and readable, but a near-vertical or
% strongly curved object can satisfy neither anywhere -- and a label with NO
% candidate falls through to the off-object fallback, which is worse than a
% steep or slightly bowed one.  So the clamp is lifted to 90 degrees, and then
% the bow tolerance widened, only when the pass before it found nothing.
ROTCAP = [min(maxrot, 90), 90, 90];
DEVCAP = [0.5, 0.5, 1.5] * lh;
for relax_ = 1:numel(ROTCAP)
rotcap = ROTCAP(relax_); devcap = DEVCAP(relax_);
rows = zeros(0, 4);
% Candidate positions are thinned to one every twelfth of a label width along
% the curve.  The polyline is densified to a few cells per step, so keeping
% every point buries the ranked list in near-duplicates of the same spot.  Not
% coarser than this: these labels are wide (a third of the panel), so a step of
% half a label width leaves a curve with under a dozen positions to offer, and
% a label whose good stretch is taken then has nowhere at all to go.
step_s = max(half/6, 1e-3);
s_last = -Inf;
for i = 1:n
    if sarc(i) - s_last < step_s, continue; end
    s_last = sarc(i);
    % Window ends: the first/last point within half a label width along the
    % curve.  Near an end of the curve the window is truncated -- reject those,
    % since a label centred there would overhang the end of its own object.
    ia = IA(i); ib = IB(i);
    if ib - ia < 2, continue; end
    if sarc(i) - sarc(ia) < 0.45*half || sarc(ib) - sarc(i) < 0.45*half, continue; end

    dx = x(ib) - x(ia); dy = y(ib) - y(ia);
    if hypot(dx, dy) < 1e-9, continue; end
    % Rotation is the angle ON THE PAGE, not in normalized units.
    ang = atan2d(dy, dx*aspect);
    if ang >  90, ang = ang - 180; end
    if ang < -90, ang = ang + 180; end
    if abs(ang) > rotcap, continue; end   % too steep to read

    % How far the curve bows away from the straight line the text will follow,
    % over the window the text covers.  A label pasted across a bend is not
    % parallel to anything, however good its clearance; measured here so it can
    % be both penalised and, past half a label height, rejected outright.
    dev = max(abs(((x(ia:ib) - x(ia))*dy - (y(ia:ib) - y(ia))*dx))) / hypot(dx, dy);
    if dev > devcap, continue; end

    nrm = [-dy, dx] / hypot(dx, dy);

    if is_region
        % Regions ride their edge at SEVERAL inward depths rather than one.
        % Depth is what a region label needs when its edge is crowded: stepping
        % further into the fill finds room without letting go of the boundary,
        % so the label stays parallel to it and stays next to the thing that
        % gives the region its meaning.  Sampling the whole interior instead
        % lets a label drift to the far side of a large polygon -- which is how
        % a band label ends up at the top of the panel, nowhere near its band.
        %
        % Depth is measured in label heights from the boundary, and the
        % shallowest row that works is what reads as belonging to that
        % boundary: 0.5 padded label heights is ~0.65 of a real line of text
        % from the edge, clear of the glyphs without the label drifting into
        % the middle of the fill.  Deeper rows are
        % handicapped heavily -- they exist only to escape a crowded edge.
        if ~inpolygon(x(i) + lh*nrm(1), y(i) + lh*nrm(2), x, y), nrm = -nrm; end
        for dmul = [0.5, 0.8, 1.2]
            px = x(i) + dmul*lh*nrm(1);  py = y(i) + dmul*lh*nrm(2);
            % Out of the band: deeper can only be further out, so stop.
            if ~inpolygon(px, py, x, y), break; end
            if px < 0 || px > 1 || py < 0 || py > 1, continue; end
            ii = round(px*nx + 0.5); jj = round(py*ny + 0.5);
            if ii < 1 || ii > nx || jj < 1 || jj > ny, continue; end
            rows(end+1,:) = [px, py, ang, clearf(jj,ii) - 60*(dev/max(lh,eps)) ...
                             - pen(jj,ii) - 40*max(0, dmul - 0.5)]; %#ok<AGROW>
        end
    else
        % Lines: BOTH sides are offered, each ranked by its own clearance.
        % Offering only the roomier one is why a curve whose good side is
        % taken -- routinely, by the label of the region it bounds -- ends up
        % with no candidate at all and falls through to the crowded fallback.
        off = 0.45 * lh;
        for sgn = [1, -1]
            px = x(i) + sgn*off*nrm(1);  py = y(i) + sgn*off*nrm(2);
            if px < 0 || px > 1 || py < 0 || py > 1, continue; end
            ii = round(px*nx + 0.5); jj = round(py*ny + 0.5);
            if ii < 1 || ii > nx || jj < 1 || jj > ny, continue; end
            rows(end+1,:) = [px, py, ang, ...
                             clearf(jj,ii) - 60*(dev/max(lh,eps)) - pen(jj,ii)]; %#ok<AGROW>
        end
    end
end
if ~isempty(rows), break; end
end

c = zeros(0,3);
if isempty(rows), return; end
[~, ord] = sort(rows(:,4), 'descend');
c = rows(ord, 1:3);
end

function ang = local_angle(poly, i, lw, aspect, maxrot)
% Page angle of the object at (or nearest) point I, over a window one label
% wide -- the same definition ranked_candidates uses, so an off-object label
% ends up at the same angle it would have had on the object.
x = poly(:,1); y = poly(:,2); n = numel(x);
ang = 0;
if n < 3, return; end
sarc = [0; cumsum(hypot(diff(x), diff(y)))];
half = max(lw/2, 1e-3);
ia = find(sarc >= sarc(i) - half, 1, 'first');
ib = find(sarc <= sarc(i) + half, 1, 'last');
if isempty(ia), ia = 1; end
if isempty(ib), ib = n; end
if ib <= ia, ia = max(1, i-1); ib = min(n, i+1); end
dx = x(ib) - x(ia); dy = y(ib) - y(ia);
if hypot(dx, dy) < 1e-9, return; end
ang = atan2d(dy, dx*aspect);
if ang >  90, ang = ang - 180; end
if ang < -90, ang = ang + 180; end
ang = max(-maxrot, min(maxrot, ang));
end

function c = clamp_col(c)
% Darken a pale region tint until the TEXT is readable on the page.
%
% The old rule was min(c, 0.75), an element-wise cap.  That neither controls
% readability nor preserves the colour: the Delta-Neff tint [0.95 0.70 0.10]
% came out [0.75 0.70 0.10], still at relative luminance 0.67 -- pale yellow on
% a pale panel -- and with its hue shifted by clipping only the red channel.
%
% Scale toward black instead, which holds the hue exactly (all three channels by
% one factor) and puts a ceiling on luminance, the quantity that actually sets
% contrast against a light background.  L_MAX = 0.40 leaves saturated mid-tones
% (SN1987A pink, L = 0.45) essentially untouched and only bites on the pale
% yellows and creams, where it is the whole point: [0.95 0.70 0.10] -> luminance
% 0.667, scaled by 0.60 to [0.57 0.42 0.06].
if ischar(c) || isstring(c), c = [0 0 0]; return; end
if numel(c) ~= 3, c = [0 0 0]; return; end
c = max(min(c(:).', 1), 0);
L = 0.2126*c(1) + 0.7152*c(2) + 0.0722*c(3);   % ITU-R BT.709 relative luminance
L_MAX = 0.40;
if L > L_MAX, c = c * (L_MAX / L); end
end
