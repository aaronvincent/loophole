function S = gstar_interp(Tlo, Thi, n)
% gstar_interp  Precompute g_s, g_rho and d(ln g_s)/d(ln T) on a log-T grid
% and return pchip interpolants (keyed on log T) for fast lookup inside an
% ODE RHS.  Replaces 4 gstarofT_vec calls per RHS evaluation with 3 cheap
% griddedInterpolant evaluations.
%
%   S.gs(logT)    -> g_s(T)
%   S.grho(logT)  -> g_rho(T)
%   S.dlngs(logT) -> d ln g_s / d ln T   (matches the central-difference
%                    form used previously in freeze_out.m / freeze_loop.m)
%
% The grid is padded slightly beyond [Tlo,Thi] so the ODE solver can probe
% the endpoints (and the dT finite-difference stencil) without extrapolating.

if nargin < 3 || isempty(n), n = 400; end

% pad +/- ~12% in log10 so endpoint lookups stay inside the sample range
pad = 0.05;
lT  = linspace(log10(Tlo) - pad, log10(Thi) + pad, n);
Tg  = 10.^lT;

gs   = gstarofT_vec(Tg, 1);
grho = gstarofT_vec(Tg, 0);

% d ln g_s / d ln T via the same central difference used inline before, so
% the entropy correction (1 + dlngs/3) reproduces the original to roundoff.
dT          = Tg * 1e-5;
dlngs_dlnT  = (gstarofT_vec(Tg+dT,1) - gstarofT_vec(Tg-dT,1)) ./ (2*dT) .* (Tg ./ gs);

lTln = log(Tg);   % natural-log abscissa (callers pass log(T))
S.gs    = griddedInterpolant(lTln(:), gs(:),         'spline');
S.grho  = griddedInterpolant(lTln(:), grho(:),       'spline');
S.dlngs = griddedInterpolant(lTln(:), dlngs_dlnT(:), 'spline');
end
