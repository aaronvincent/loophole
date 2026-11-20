function Omega_h2 = freeze_out(mx, mphi, gx, g, interaction, doplot)
% freeze_out  Thermal WIMP relic abundance via freeze-out
%
% Implements the numerical calculation from arXiv:1204.3622 (Steigman,
% Dasgupta & Beacom 2012).  Solves the log-transformed Boltzmann equation
%
%   dW/dx = (lambda/x^2) * [1 + (1/3) d(ln g_s)/d(ln T)]
%           * (g_s / sqrt(g_rho)) * (exp(2*W_eq - W) - exp(W))   [Eq. 333]
%
% where W = ln(Y), Y = n/s, x = m/T, and lambda = (1/2) sqrt(pi/45) Mpl mx sv(T)
% (the 1/2: Dirac chi, Y counts chi + chibar).
% sv(T) is evaluated via svavg.m (thermally averaged cross section).
%
% Requires gstarofT_vec.m, svavg.m, xxggxsec.m in the MATLAB path.
%
% Usage:
%   Oh2 = freeze_out(0.1, 10, 1, [2 -1 -1]/3, 'heavyScalar')

if nargin < 6, doplot = false; end

g_chi = 4;   % internal d.o.f.: 4 for a Dirac fermion

% -----------------------------------------------------------------------
% Precompute svavg over the ODE integration range and build interpolant
% -----------------------------------------------------------------------
Tgrid    = logspace(log10(mx/1000), log10(mx), 200);
logTgrid = log(Tgrid);
svvals   = svavg_grid(mx, mphi, gx, g, Tgrid, interaction);
svInterp = griddedInterpolant(logTgrid(:), svvals(:), 'pchip');

% Precompute g_s, g_rho, d(ln g_s)/d(ln T) interpolants over the T range
% probed by the ODE (T = mx/x, x in [1,1000]).  Avoids 4 gstarofT_vec
% calls per RHS evaluation.
S = gstar_interp(mx/1000, mx);

% -----------------------------------------------------------------------
% Initial condition: start in thermal equilibrium at x = 1
% -----------------------------------------------------------------------
x0  = 1;
T0  = mx / x0;
% Use the SAME g_s source as the RHS (S.gs), not the analytic gstarofT_vec:
% at large lambda a tiny g_s mismatch between W0 and W_eq(x=1) makes the
% stiff bracket (exp(2 W_eq - W) - exp(W)) nonzero at the start and blows up.
gs0 = S.gs(log(T0));
W0  = log(0.145 * g_chi / gs0) + 1.5 * log(x0) - x0;

% -----------------------------------------------------------------------
% Integrate from x = 1 to x = 1000
% -----------------------------------------------------------------------
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[x_sol, W_sol] = ode15s(@rhs, [1, 1000], W0, opts);

Y0_now = exp(W_sol(end));

% -----------------------------------------------------------------------
% Relic density  Omega h^2 = m * s0 * Y0 / (rho_crit / h^2)
% -----------------------------------------------------------------------
s0          = 2891.2;    % [cm^-3]
rho_crit_h2 = 1.054e-5;  % [GeV cm^-3]
Omega_h2    = mx * s0 * Y0_now / rho_crit_h2;

% -----------------------------------------------------------------------
% Optional plot
% -----------------------------------------------------------------------
if doplot
    x_plot    = logspace(0, 3, 600);
    T_plot    = mx ./ x_plot;
    gs_plot   = gstarofT_vec(T_plot, 1);
    Y_eq_plot = 0.145 * g_chi ./ gs_plot .* x_plot.^1.5 .* exp(-x_plot);

    figure('Name', 'WIMP freeze-out');
    semilogy(x_sol, exp(W_sol), 'b-', 'LineWidth', 2); hold on;
    semilogy(x_plot, Y_eq_plot, 'r--', 'LineWidth', 1.5);
    xlabel('x = m / T');
    ylabel('Y = n / s');
    title(sprintf('Freeze-out:  m = %.4g GeV,  g_x = %.4g', mx, gx));
    legend('Y (numerical)', 'Y_{eq}', 'Location', 'SouthWest');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([1 1000]); grid on;
    text(0.97, 0.97, sprintf('\\Omega h^2 = %.3f', Omega_h2), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'FontSize', 12);
end

% -----------------------------------------------------------------------
% Nested RHS  (captures mx, g_chi, svInterp from outer scope)
% -----------------------------------------------------------------------
    function dWdx = rhs(x, W)
        T     = mx / x;
        lT    = log(T);
        g_rho = S.grho(lT);
        g_s   = S.gs(lT);

        % d(ln g_s)/d(ln T) from the precomputed interpolant
        correction     = 1 + S.dlngs(lT) / 3;

        % temperature-dependent lambda from interpolated svavg
        % sv [GeV^-2]; lambda = (1/2) sqrt(pi/45)*Mpl*mx*sv (natural units throughout).
        % Y counts chi + chibar (g_chi = 4) and <sigma v> is the chi chibar cross
        % section, so n-dot = -(<sigma v>/2)(n^2 - n_eq^2): the 1/2 is the Dirac factor.
        sv           = max(0, svInterp(log(T)));
        lambda_local = 0.5 * sqrt(pi/45) * 1.22e19 * mx * sv;

        W_eq  = log(0.145 * g_chi / g_s) + 1.5 * log(x) - x;
        dWdx  = (lambda_local / x^2) * correction * (g_s / sqrt(g_rho)) ...
                * (exp(2*W_eq - W) - exp(W));
    end

end
