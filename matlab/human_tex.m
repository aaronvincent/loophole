% human_tex.m
% Shared by compare_limits_g.m.
%
% Was a verbatim copy in BOTH scripts.  Two copies means every fix has to be
% applied twice.  The helpers deliberately NOT extracted -- nd_cached,
% neff_pull, neff_edges_off, svgg_at_fo, ceiling_lookup, draw_neff_band,
% decade_ticks, tighten_axes, save_batch_figure -- have already DRIFTED
% between the two files and are not interchangeable; unifying them needs a
% reading of what each difference is for, not a copy.
function s = human_tex(m)
% Human-readable number: 1, 0.1, 3, 0.001 stay decimal; anything smaller or
% larger than that window goes to 10^{-9} form via sci_tex.
if m == 0
    s = '0';
elseif m >= 1e-3 && m < 1e4
    s = strtrim(sprintf('%g', m));
else
    s = sci_tex(m);
end
end
