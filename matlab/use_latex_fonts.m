function use_latex_fonts()
% USE_LATEX_FONTS  Make every new figure render text in the LaTeX font.
%
% Sets the session-wide graphics defaults so that titles, axis labels,
% tick labels, legends, colorbars, and text/annotation objects all use the
% LaTeX interpreter (Computer Modern math font) instead of MATLAB's default
% TeX interpreter + Helvetica.  Call once at the top of a plotting script;
% existing per-object 'Interpreter','latex' settings are unaffected.

set(groot, 'defaultTextInterpreter',                'latex');
set(groot, 'defaultAxesTickLabelInterpreter',       'latex');
set(groot, 'defaultLegendInterpreter',              'latex');
set(groot, 'defaultColorbarTickLabelInterpreter',   'latex');
end
