% corner_caption.m
% Shared by compare_limits_g.m.
%
% Was a verbatim copy in BOTH scripts.  Two copies means every fix has to be
% applied twice.  The helpers deliberately NOT extracted -- nd_cached,
% neff_pull, neff_edges_off, svgg_at_fo, ceiling_lookup, draw_neff_band,
% decade_ticks, tighten_axes, save_batch_figure -- have already DRIFTED
% between the two files and are not interchangeable; unifying them needs a
% reading of what each difference is for, not a copy.
function corner_caption(lines, fs)
% Replaces the figure title: same information, parked inside the axes at the
% bottom-right so it does not eat vertical space above the frame.
text(0.98, 0.03, lines, 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
     'Interpreter', 'latex', 'FontSize', fs);
end
