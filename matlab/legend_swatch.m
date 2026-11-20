% legend_swatch.m
% Shared by compare_limits_g.m.
%
% Was a verbatim copy in BOTH scripts.  Two copies means every fix has to be
% applied twice.  The helpers deliberately NOT extracted -- nd_cached,
% neff_pull, neff_edges_off, svgg_at_fo, ceiling_lookup, draw_neff_band,
% decade_ticks, tighten_axes, save_batch_figure -- have already DRIFTED
% between the two files and are not interchangeable; unifying them needs a
% reading of what each difference is for, not a copy.
function legend_swatch(col, alpha, lbl, style)
% Legend-only proxy for a shaded region: a patch with no vertices, so it draws
% nothing in the axes but still claims a legend row.  The swatch uses a
% stronger tint and a coloured edge than the region itself -- at the 0.10-0.18
% FaceAlpha the regions need on the plot, ten overlapping pale fills give ten
% nearly indistinguishable legend squares.  The real patch is created with
% HandleVisibility 'off' and keeps its faint, edgeless fill.
%
% STYLE sets the swatch's border linestyle so it matches how that constraint is
% drawn on the plot: pass '-.' for the ratio exclusions, which do draw a
% dash-dot boundary, so the legend box reads as the same object.  Defaults to
% '-' for the regions that are filled with no boundary at all (Overabundant,
% Delta-Neff, non-perturbative, the whole-panel externals) -- those keep a solid
% border purely so the swatch stays legible; an edgeless one is unreadable at
% these alphas.
if nargin < 4 || isempty(style), style = '-'; end
% Capture the region this swatch stands for.  legend_swatch is always called
% immediately after the real patch is drawn, so that patch is the newest child
% of the axes right now.  inline_labels needs it: the proxy is vertex-less and
% cannot be labelled in place.
geom = gobjects(0);
kids = allchild(gca);   % NOT get(gca,'Children'): the region patch is drawn
                        % with HandleVisibility 'off' and is absent from that list
% SEARCH BACK FOR THE PATCH, do not assume it is kids(1).  draw_neff_band closes
% the ends of a band that stops inside the mass range with dash-dot verticals,
% drawn AFTER the region patch, so for a partial band the newest child is a Line
% and the capture failed -- leaving inline_labels with a vertex-less proxy it
% skips by design.  The symptom is a shaded Delta-Neff region with no label,
% only on the panels where the band does not span every mass.
for kk = 1:numel(kids)
    if isa(kids(kk), 'matlab.graphics.primitive.Patch')
        if ~all(isnan(kids(kk).XData(:))), geom = kids(kk); end
        break                     % first patch back is the region just drawn
    end
end
patch('XData', NaN, 'YData', NaN, 'FaceColor', col, ...
      'FaceAlpha', min(1, 2.5*alpha), 'EdgeColor', col, 'LineWidth', 1, ...
      'LineStyle', style, 'DisplayName', lbl, 'UserData', struct('geom', geom));
end
