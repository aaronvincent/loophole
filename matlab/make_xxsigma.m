function [sig, Gphi, single_pole] = make_xxsigma(mx, mphi, gx, g, interaction, channel)
% make_xxsigma  Build the chi-chibar annihilation cross-section evaluator
% sigma(s) used by svavg.m / svavg_grid.m.  Single source of truth so the
% scalar (adaptive) and vectorized (fixed-grid) integrators stay identical.
%
%   channel = 'all' : full cross-section from xxggxsec.m
%   channel = 'gg'  : s-channel chi-chibar -> gamma gamma piece only
%   channel = 'sm'  : SM final states (gamma gamma + meson pairs), xxggxsec.m
%   channel = 'had' : meson pairs only, xxggxsec.m
%   channel = 'phiphi' : dark mediator pair only
%
% Returns:
%   sig         : function handle, vectorized over s
%   Gphi        : phi total width (BW propagator regulator)
%   single_pole : true for the single-Breit-Wigner-pole interactions
%                 (heavyScalar/lightScalar/heavyScalarPS/lightScalarPS),
%                 which are eligible for the fast fixed-grid + NWA path.
%                 SP/PP interactions have multi-pole structure -> false.

if nargin < 6 || isempty(channel), channel = 'all'; end

mpi0  = 0.13498;  mpi_pm = 0.13957;
mK    = 0.493677;
m_eta = 0.547862;
B0      = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV (draft convention; was 2.6)
alphaEM = 1/137;
Dscript = 0.5*(B0*(g(1)+g(2))/mpi_pm^2 + B0*(g(1)+g(3))/mK^2)^2;   % charged-meson loop (pi+, K+), low-energy limit A0 -> 1/3

% Single source of truth for the mediator total width.
Gphi = phi_partial_widths(mphi, mx, gx, g, interaction).total;

single_pole = ismember(interaction, ...
    {'heavyScalar','lightScalar','heavyScalarPS','lightScalarPS'});

if strcmp(channel,'gg')
    prop2 = @(s) (s - mphi^2).^2 + mphi^2*Gphi^2;
    switch interaction
        case {'heavyScalar','lightScalar'}
            sig = @(s) (s > 4*mx^2) .* ...
                  gx^2*alphaEM^2*Dscript/(9216*pi^3) ...
                  .* s.*(s - 4*mx^2)./prop2(s) .* sqrt(s./max(s - 4*mx^2, eps));
        case {'heavyScalarPS','lightScalarPS'}
            sig = @(s) (s > 4*mx^2) .* ...
                  gx^2*alphaEM^2*Dscript/(9216*pi^3) .* s.^2 ./ prop2(s) ...
                  .* sqrt(s./max(s - 4*mx^2, eps));
        case {'heavyScalarSP','lightScalarSP'}
            % {pi0, phi', eta} generically, or {eta', phi'} for g ~ [1 1 1].
            [mj_arr, gcx_arr, ggx_arr, Gamma_arr] = sppp_couplings(mphi, mx, gx, g, interaction);
            sig = @(s) sig_SPgg(s, mx, mj_arr, gcx_arr, ggx_arr, Gamma_arr);
        case {'scalarPP','heavyScalarPP','lightScalarPP'}
            % {pi0, phi', eta} generically, or {eta', phi'} for g ~ [1 1 1].
            [mj_arr, gcx_arr, ggx_arr, Gamma_arr] = sppp_couplings(mphi, mx, gx, g, interaction);
            sig = @(s) sig_PPgg(s, mx, mj_arr, gcx_arr, ggx_arr, Gamma_arr);
        case 'axialVector'
            [mj, lam_chi, lam_gg, Gam_j] = axial_couplings(mphi, mx, gx, g);
            sig = @(s) sig_AAgg(s, mx, mj, lam_chi, lam_gg, Gam_j);
        otherwise
            sig = @(s) zeros(size(s));
    end
elseif strcmp(channel,'phiphi')
    % chi chibar -> phi phi (t-channel) piece, isolated from the full
    % xxggxsec sum.  Used only by svavg_phiphi.m, to test dark-sector
    % internal equilibration -- NOT a detailed-balance partner
    % of any production channel, so it is kept out of the 'gg' path.  The
    % closed forms are copied verbatim from xxggxsec.m (single source of
    % truth for the physics); both groups open at s > 4*mphi^2, the only
    % physical threshold.  Zero for interactions with no phi-phi channel.
    switch interaction
        case {'heavyScalar','lightScalar'}
            sig = @(s) sig_pp_scalar(s, mx, mphi, gx);
        case {'heavyScalarPS','lightScalarPS'}
            sig = @(s) sig_pp_ps(s, mx, mphi, gx);
        case 'axialVector'
            sig = @(s) sig_pp_axial(s, mx, mphi, gx);
        otherwise
            sig = @(s) zeros(size(s));
    end
else
    sig = @(s) xxggxsec(mx,mphi,gx,g,s,interaction,channel);
end
end

% --- chi chibar -> phi phi, scalar phi-chi-chi vertex --------------------
function s_out = sig_pp_scalar(s, mx, mphi, gx)
    s_out   = zeros(size(s));
    mask_pp = (s > 4*mphi^2);
    if any(mask_pp)
        sp    = s(mask_pp);
        kin   = sqrt((sp - 4*mx^2).*(sp - 4*mphi^2));
        A     = -32*mx^4 + 16*mx^2*(sp - mphi^2) + 6*mphi^4 - 4*mphi^2*sp + sp.^2;
        atarg = sqrt((sp - 4*mx^2).*(4*mphi^2 - sp)) ./ (2*mphi^2 - sp);
        t1    = 2*(mx^2*(sp - 4*mphi^2) + mphi^4) .* A .* atan(atarg) ...
                ./ ((sp - 4*mx^2).^(3/2) .* (2*mphi^2 - sp) .* sqrt(4*mphi^2 - sp));
        t2    = (16*mx^4 + 2*mx^2*(sp - 8*mphi^2) + 3*mphi^4) ./ (4*mx^2 - sp);
        s_out(mask_pp) = real(gx^4 * kin .* (t1 + t2) ...
                       ./ (32*pi * sp .* (mx^2*(sp - 4*mphi^2) + mphi^4)));
    end
end

% --- chi chibar -> phi phi, pseudoscalar phi-chi-chi vertex --------------
function s_out = sig_pp_ps(s, mx, mphi, gx)
    s_out   = zeros(size(s));
    mask_pp = (s > 4*mphi^2);
    if any(mask_pp)
        % REPLACED 2026-08-11 in step with xxggxsec.m: the old closed form was
        % negative and ~1/v_rel where the true result is positive and ~v_rel.
        % The 1/16pi already carries the identical-particle S = 1/2.
        sp    = s(mask_pp);
        kin   = sqrt((sp - 4*mx^2).*(sp - 4*mphi^2));
        t1    = -sqrt((sp - 4*mphi^2)./(sp - 4*mx^2)) ...
                .* (2*mx^2*(sp - 4*mphi^2) + 3*mphi^4) ...
                ./ (2*(mx^2*(sp - 4*mphi^2) + mphi^4));
        t2    = (6*mphi^4 - 4*mphi^2*sp + sp.^2) ...
                ./ ((sp - 4*mx^2).*(sp - 2*mphi^2)) .* atanh(kin ./ (sp - 2*mphi^2));
        s_out(mask_pp) = gx^4*(t1 + t2) ./ (16*pi*sp);
    end
end

% --- chi chibar -> A' A', axial-vector mediator (t/u-channel chi exchange) --
% Verbatim copy of the A'A' piece of xxggxsec.m's axialVector chi chibar -> all
% (the single source of truth), isolated here so svavg_phiphi can probe the
% chi <-> A' internal-equilibration rate.  Depends only on gx (= the
% axial coupling lambda_a) and the masses -- NOT on axial_couplings' meson
% mixing -- so it never trips the mediatorTooLight validity floor.  Opens at
% s = 4 mphi^2, the only physical threshold.
function s_out = sig_pp_axial(s, mx, mphi, gx)
    s_out   = zeros(size(s));
    mask_pp = (s > 4*mphi^2);
    if any(mask_pp)
        sp    = s(mask_pp);
        sqneg = sqrt((4*mphi^2 - sp).*(sp - 4*mx^2));   % imaginary for sp>4mphi^2
        atrm  = atan(sqneg ./ (2*mphi^2 - sp));
        base  = mphi^4 - 4*mphi^2*mx^2 + mx^2*sp;
        poly  = 16*mphi^8 - 16*mphi^6.*(4*mx^2 + sp) + 2*mphi^4.*sp.*(8*mx^2 + 5*sp) ...
                - 4*mphi^2.*sp.*(-16*mx^4 - 2*mx^2*sp + sp.^2) - 16*mx^4*sp.^2 + sp.^4;
        num   = -8*mphi^8 + 8*mphi^6.*(2*mx^2 + sp) ...
                + mphi^4.*(64*mx^4 - 28*mx^2*sp - 3*sp.^2) ...
                + 4*mphi^2*mx^2.*sp.*(3*sp - 16*mx^2) ...
                + 2*base.*poly.*atrm ./ ((2*mphi^2 - sp).*sqneg) ...
                - 2*mx^2*sp.^2.*(sp - 4*mx^2);
        s_out(mask_pp) = real( gx^4 * sqrt((sp - 4*mphi^2)./(sp - 4*mx^2)) ...
                       .* num ./ (32*pi*mphi^4.*sp.*base) );
    end
end

% --- gg-only s-channel cross-section for the SP interaction --------------
% BW-regulated complex propagators (matches xxggxsec.m SP gamma gamma piece);
% the widths Gamma_arr = [Gpi, Gphi, G_eta] keep the pi0/phi'/eta poles finite.
function s_out = sig_SPgg(s, mx, mj_arr, gcx_arr, ggx_arr, Gamma_arr)
    sum_prop = zeros(size(s));
    for jj = 1:numel(mj_arr)
        sum_prop = sum_prop + gcx_arr(jj)*ggx_arr(jj) ...
                   ./ (mj_arr(jj)^2 - s - 1i*mj_arr(jj)*Gamma_arr(jj));
    end
    s_out = (s > 4*mx^2) .* s.^(3/2) .* sqrt(max(s - 4*mx^2, eps)) / (128*pi) .* abs(sum_prop).^2;
end

% --- gg-only s-channel cross-section for the PP interaction --------------
function s_out = sig_PPgg(s, mx, mj_arr, gcx_arr, ggx_arr, Gamma_arr)
    sum_prop = zeros(size(s));
    for jj = 1:numel(mj_arr)
        sum_prop = sum_prop + gcx_arr(jj)*ggx_arr(jj) ...
                   ./ (mj_arr(jj)^2 - s - 1i*mj_arr(jj)*Gamma_arr(jj));
    end
    s_out = (s > 4*mx^2) .* s.^2/(128*pi) .* abs(sum_prop).^2 .* sqrt(s./max(s - 4*mx^2, eps));
end

% --- gg-only s-channel cross-section for the axialVector interaction ------
% chi chibar -> gamma gamma via the three SM-meson mixing poles (pi0/eta/eta').
% Matches the gamma gamma piece of xxggxsec.m's axialVector branch (and is the
% detailed-balance partner of ggxxxsec.m's gg->xx).
function s_out = sig_AAgg(s, mx, mj, lam_chi, lam_gg, Gam_j)
    sum_prop = zeros(size(s));
    for jj = 1:3
        sum_prop = sum_prop + lam_chi(jj)*lam_gg(jj) ...
                   ./ (s - mj(jj)^2 + 1i*mj(jj)*Gam_j(jj));
    end
    s_out = (s > 4*mx^2) .* mx^2 .* s.^2/(32*pi) .* abs(sum_prop).^2 ...
            .* sqrt(s./max(s - 4*mx^2, eps));
end
