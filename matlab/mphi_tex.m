% mphi_tex.m
% Shared by compare_limits_g.m.
function s = mphi_tex(interaction)
% Symbol for the mediator mass, as in the paper's Lagrangian: m_A for the A',
% m_phi otherwise.
if startsWith(interaction, 'axial'), s = 'm_A';
else,                                s = ['m_{' phi_tex(interaction) '}'];
end
end
