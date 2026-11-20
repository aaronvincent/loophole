function Oh2 = Oh2_FI(mx, mphi, gx, g, interaction, T_RH, Mpl, s0, rhoc, T0, production)
% Oh2_FI  Linearised freeze-in relic abundance at coupling gx.
%
% Production from the SM bath: gamma gamma (ggxxrate) plus meson pairs.
% production (optional): 'sm' (default) or 'gg' (photons only).
% Root-finding (not Oh2 ~ gx^2 rescaling) is needed inside the BW resonance,
% where Gamma_phi ~ gx^2 breaks Oh2 ~ gx^2.  See docs/cross_sections.md.
%
% This is freeze_loop's equation with the Y^2 term dropped.  With
% R = n_eq^2 <sigma v>/4 (ggxxrate's normalisation, and that of the meson term
% below), the Dirac-normalised production of chi + chibar is dn/dt = 2R, so
%
%   Y = int_{T_lo}^{T_RH} 2R / (s H) * (1 + dln g_s / 3 dln T) dln T,
%   s = (2 pi^2/45) g_s T^3,   H = sqrt(4 pi^3 g_rho / 45) T^2 / Mpl,
%
% and Oh2 = mx * s0 * Y / (rho_c/h^2) with freeze_loop's constants.  The rate
% is interpolated in log (it falls as e^{-2 mx/T}).  The s0, rhoc arguments are
% kept for the call signature and not used.

if nargin < 11 || isempty(production), production = 'sm'; end
T_hi = T_RH;
% Below mx/200 the rate is < e^-400 of its peak.
T_lo = max(T0, min(mx/200, T_RH/1e3));
Tgrid = logspace(log10(T_lo), log10(T_hi), 300);
R = ggxxrate(Tgrid, mx, mphi, gx, g, interaction);
% Meson pairs (pi pi, K K, ... -> chi chibar) by detailed balance with
% chi chibar -> mesons.  <sigma v> on a coarse grid above mx/60, log-interpolated.
T_lo_h = max(T_lo, min(mx/60, T_RH/100));
Th     = logspace(log10(T_lo_h), log10(T_hi), 60);
sv_h   = zeros(size(Th));
if strcmp(production, 'sm')
    sv_h = svavg_grid(mx, mphi, gx, g, Th, interaction, 'had');
end
if any(sv_h > 0)
    lsv  = griddedInterpolant(log(Th(:)), log(max(sv_h(:), realmin)), 'pchip');
    in_h = Tgrid >= T_lo_h;
    neq  = 4/(2*pi^2) * mx^2 .* Tgrid(in_h) .* besselk(2, mx ./ Tgrid(in_h));
    R(in_h) = R(in_h) + neq.^2 .* exp(lsv(log(Tgrid(in_h)))) / 4;
end
lR = griddedInterpolant(log(Tgrid(:)), log(max(R(:), realmin)), 'pchip');
S  = gstar_interp(T_lo, T_hi);
    function f = dYdlnT(lT)
        T = exp(lT);
        s = 2*pi^2/45 * S.gs(lT) .* T.^3;
        H = sqrt(4*pi^3/45 * S.grho(lT)) .* T.^2 / Mpl;
        f = 2 * exp(lR(lT)) ./ (s .* H) .* (1 + S.dlngs(lT)/3);
    end
Y = integral(@dYdlnT, log(T_lo), log(T_hi), 'RelTol', 1e-8);

s0_cm3      = 2891.2;    % cm^-3, as freeze_loop
rho_crit_h2 = 1.054e-5;  % GeV cm^-3
Oh2 = mx * s0_cm3 * Y / rho_crit_h2;
end
