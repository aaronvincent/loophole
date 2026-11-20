% gx_tex.m
% Shared by compare_limits_g.m.
function s = gx_tex(interaction)
% Symbol for the mediator-chi coupling gx, as in the paper's Lagrangian:
% the dark vertex (first suffix letter) picks g_chi (S) or g_5chi (P), and a
% pseudoscalar SM vertex (phi-tilde) adds the tilde; lambda_chi for the A'.
if     startsWith(interaction, 'axial'), s = '\lambda_\chi';
elseif endsWith(interaction, 'PP'),      s = '\tilde{g}_{5\chi}';
elseif endsWith(interaction, 'SP'),      s = '\tilde{g}_\chi';
elseif endsWith(interaction, 'PS'),      s = 'g_{5\chi}';
else,                                    s = 'g_\chi';
end
end
