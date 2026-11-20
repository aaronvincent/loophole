% phi_tex.m
% Shared by compare_limits_g.m.
function s = phi_tex(interaction)
% Mediator symbol in the paper's Lagrangian: phi for every scalar family
% (S x S, P x S, S x P, P x P), A' for the axial vector.
if startsWith(interaction, 'axial'), s = 'A^{\prime}';
else,                                s = '\phi';
end
end
