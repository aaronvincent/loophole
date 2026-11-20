% model_tex.m
% Shared by compare_limits_g.m.
%
% Was a verbatim copy in BOTH scripts.  Two copies means every fix has to be
% applied twice.  The helpers deliberately NOT extracted -- nd_cached,
% neff_pull, neff_edges_off, svgg_at_fo, ceiling_lookup, draw_neff_band,
% decade_ticks, tighten_axes, save_batch_figure -- have already DRIFTED
% between the two files and are not interchangeable; unifying them needs a
% reading of what each difference is for, not a copy.
function s = model_tex(interaction)
% Vertex structure as (chi vertex) (x) (SM vertex), matching the SxS/PxS/PxP/
% SxP/AxA naming used in kaon_invisible_ratio.m (whose AxA is an internal key;
% the axial label carries primes here because the mediator is the A').  The
% heavy/light prefix is inert for the relic curves, so it does not appear.
if     startsWith(interaction, 'axial'),          s = 'A^{\prime} \otimes A^{\prime}';
elseif endsWith(interaction, 'PS'),               s = 'P \otimes S';
elseif endsWith(interaction, 'PP'),               s = 'P \otimes P';
elseif endsWith(interaction, 'SP'),               s = 'S \otimes P';
else,                                             s = 'S \otimes S';
end
end
