% sci_tex.m
% Shared by compare_limits_g.m.
%
% Was a verbatim copy in BOTH scripts.  Two copies means every fix has to be
% applied twice.  The helpers deliberately NOT extracted -- nd_cached,
% neff_pull, neff_edges_off, svgg_at_fo, ceiling_lookup, draw_neff_band,
% decade_ticks, tighten_axes, save_batch_figure -- have already DRIFTED
% between the two files and are not interchangeable; unifying them needs a
% reading of what each difference is for, not a copy.
function s = sci_tex(x)
% Format x in LaTeX scientific notation: 1e-3 -> 10^{-3}, 3e-4 -> 3 \times 10^{-4}.
if x == 0, s = '0'; return; end
e = floor(log10(abs(x)));
m = x / 10^e;
if abs(m - 1) < 1e-9
    s = sprintf('10^{%d}', e);
elseif abs(m - round(m)) < 1e-9
    s = sprintf('%d \\times 10^{%d}', round(m), e);
else
    s = sprintf('%.2f \\times 10^{%d}', m, e);
end
end
