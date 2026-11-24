function r = mesonandCMB_ratio(mx, mphi, gx, g, interaction, which)
%MESONANDCMB_RATIO  Exclusion ratios for the meson-decay and indirect-detection
%   limits beyond K+ -> pi+ chi chi (which is in kaon_invisible_ratio.m).
%
%   r = MESONANDCMB_RATIO(MX, MPHI, GX, G, INTERACTION, WHICH) returns the ratio
%   of a predicted signal to its experimental bound; r > 1 => excluded.  WHICH
%   selects the constraint:
%
%     'kphi'      Gamma(K+ -> pi+ phi) * P_inv / (K+->pi+ + invisible limit).
%                 On-shell mediator production; open only for MPHI < m_K+- -
%                 m_pi+- ~ 0.354 GeV, r = NaN above.  P_inv is the chance the
%                 phi looks invisible: it decays to chi chibar, or it reaches the
%                 end of the detector before decaying to gamma gamma / pi pi
%                 (kphi_invisible_fraction).  GX enters only through P_inv, so
%                 the bound is gx-independent wherever P_inv = 1.  r is not
%                 monotonic in the quark coupling: it falls again once the
%                 phi decays visibly inside the detector.
%     'cmb'       <sigma v>_low-v / Slatyer CMB energy-injection bound.  ~gx^2.
%                 NOT applicable for SxS/SxP (velocity-suppressed) -> r = NaN.
%     'integral'  <sigma v>_low-v / INTEGRAL telescope bound.  ~gx^2.  All models.
%     'pi0inv'    Gamma(pi0  -> chi chibar) / 3.43e-17 GeV   (pi0  -> invisible)
%     'etainv'    Gamma(eta  -> chi chibar) / 3.14e-11 GeV   (eta  -> invisible)
%     'etapinv'   Gamma(eta' -> chi chibar) / 3.95e-8  GeV   (eta' -> invisible)
%     'mesoninv'  max of the three above (NaN only if all three are NaN).
%                 The meson inherits a chi coupling by mixing with the mediator,
%                 so only the pseudoscalar-quark (SxP/PxP) and axialVector
%                 models have it; SxS/PxS (scalar phi cannot mix with a 0^-
%                 meson) -> NaN.  Open for mx < m_j/2; ~gx^2 * (mixing)^2.
%
%   Inputs mirror kaon_invisible_ratio.m / the relic-code convention:
%     MX    DM mass m_chi [GeV];  MPHI mediator mass (m_A' for axialVector) [GeV]
%     GX    phi-chi coupling (lambda_a for axialVector)
%     G     [g_u g_d g_s] phi-quark couplings (C_u, C_d, C_s)
%     INTERACTION  relic model string (heavyScalar*, lightScalar*, scalar*, axialVector)
%
%   r = NaN when the constraint is closed / not applicable to the interaction.
%
%   Covers K -> pi phi / K -> pi A', CMB and INTEGRAL.  The low-velocity <sigma v>
%   is the analytic threshold expression (s = 4 m_chi^2), NOT a v->0
%   limit of the relic svavg (that is a genuine thermal average and does not
%   reproduce the p-wave/velocity-suppressed models).  Interaction -> model map
%   matches kaon_invisible_ratio.m.

    % ---- interaction string -> model branch ---------------------------------
    switch interaction
        case {'heavyScalar',   'lightScalar',   'scalar'},   model = 'SxS';
        case {'heavyScalarPS', 'lightScalarPS', 'scalarPS'}, model = 'PxS';
        case {'heavyScalarSP', 'lightScalarSP', 'scalarSP'}, model = 'SxP';
        case {'heavyScalarPP', 'lightScalarPP', 'scalarPP'}, model = 'PxP';
        case 'axialVector',                                  model = 'AxA';
        otherwise, r = NaN; return;
    end

    % ---- constants -----------------------------------------------------------
    P.MeV = 1e-3;  P.eV = 1e-9;
    P.e = 0.3;  P.alpha = 1/137;
    P.B0 = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV (was 2.6)
    P.Fpi0 = 92.4 * P.MeV;  P.F = P.Fpi0;
    P.mpi = 134.98 * P.MeV;  P.mpi0 = P.mpi;      % neutral pion
    P.meta = 547.86 * P.MeV;                      % eta
    P.metaprime = 958 * P.MeV;  P.meta0 = P.metaprime;   % eta'
    P.mpiplus = 139.57 * P.MeV;                   % charged pion
    P.mkplus = 493.68 * P.MeV;  P.mKplus = P.mkplus;     % charged kaon
    P.mp = 0.938;
    P.GF = 1.17e-5;
    P.Vus = 0.222;  P.Vusstar = P.Vus;  P.Vstarud = 0.974;  P.Vud = P.Vstarud;
    P.g8 = 3.62;  P.g8s = -0.38 * P.g8;
    P.L4 = 1.87e-3;  P.L5 = 1.22e-3;  P.L6 = 1.46e-3;  P.L7 = -0.39e-3;
    P.mu = 2.2e-3;  P.md = 4.7e-3;  P.ms = 95e-3;
    P.mqu = P.mu;  P.mqd = P.md;  P.mqs = P.ms;
    P.hbarc = 197 * P.MeV * 1e-13;   % [GeV*cm]
    P.cc = 3e10;                     % [cm/s]
    P.hbar = P.hbarc / P.cc;         % [GeV*s]
    P.Gammapi0 = P.hbar / 8.52e-17;
    P.Gammaeta = 1.31e-6;
    P.Gammaeta0 = 0.188e-3;

    P.Cu = g(1);  P.Cd = g(2);  P.Cs = g(3);
    P.mx = mx;  P.mphi = mphi;  P.gx = gx;

    % K+ invisible-decay limit as a rate [s^-1] (same ceiling as kaon).
    lifetimekplus  = 1.238e-8;
    kplusinvlimit  = (1/lifetimekplus) * 1.14e-10;

    switch which
    % =====================================================================
    case 'kphi'    % K+ -> pi+ phi  (on-shell mediator, invisible share)
    % =====================================================================
        Gamma = kphi_rate(model, P);
        if isnan(Gamma), r = NaN; return; end
        % x NA62 acceptance at m_miss^2 = mphi^2 (R1 [0, 0.010], R2 [0.026, 0.068]
        % GeV^2, resolution 0.001; as kaon_invisible_ratio.m)
        r = Gamma * kphi_invisible_fraction(model, mx, mphi, gx, g) * na62_acceptance(mphi^2) / kplusinvlimit;

    % =====================================================================
    case 'cmb'     % <sigma v> vs Slatyer CMB bound (not for SxS/SxP)
    % =====================================================================
        if strcmp(model, 'SxS') || strcmp(model, 'SxP')
            r = NaN; return;   % velocity-suppressed: no CMB bound in the reference
        end
        % + the s-wave three-body channels, energy-weighted (threebody_cmb_sv)
        sv = sigmav_indirect(model, P) + threebody_cmb_sv(model, P);
        lim = curve_slatyer(mx);
        if ~isfinite(sv) || ~isfinite(lim), r = NaN; return; end
        r = sv / lim;

    % =====================================================================
    case 'integral'  % <sigma v> vs INTEGRAL bound (all models)
    % =====================================================================
        sv = sigmav_indirect(model, P);
        lim = curve_integral(mx);
        if ~isfinite(sv) || ~isfinite(lim), r = NaN; return; end
        r = sv / lim;

    % =====================================================================
    case {'pi0inv', 'etainv', 'etapinv', 'mesoninv'}  % neutral meson -> invisible
    % =====================================================================
        lim = [3.43e-17, 3.14e-11, 3.95e-8];   % Gamma(j -> invisible) bounds [GeV]
        rj  = meson_invisible_widths(model, mx, mphi, gx, g) ./ lim;
        switch which
            case 'pi0inv',  r = rj(1);
            case 'etainv',  r = rj(2);
            case 'etapinv', r = rj(3);
            otherwise,      r = max(rj, [], 'omitnan');   % NaN iff all NaN
        end

    otherwise
        error('mesonandCMB_ratio:which', ...
              'WHICH must be ''kphi'', ''cmb'', ''integral'', ''pi0inv'', ''etainv'', ''etapinv'' or ''mesoninv''.');
    end
end

% ======================================================================
% K+ -> pi+ phi  (on-shell mediator two-body decay rate) [s^-1]
% ======================================================================
function Gamma = kphi_rate(model, P)
    mmed = P.mphi;
    % Channel closed above threshold -> not applicable.
    if mmed >= P.mkplus - P.mpiplus
        Gamma = NaN;  return;
    end
    mk = P.mkplus;  mpip = P.mpiplus;
    kin = sqrt((mk^2 + mpip^2 - mmed^2)^2 - 4*mk^2*mpip^2);

    switch model
    case {'SxS', 'PxS'}   % scalar phi-quark coupling (identical for SS and PS)
        delta   = -P.GF * P.F^2 * P.Vusstar * P.Vud * P.g8 / sqrt(2);
        cpikphi = delta * P.B0 * ((P.Cs + P.Cu)*mpip^2 - (P.Cu + P.Cd)*mk^2) / (mk^2 - mpip^2);
        Gamma   = cpikphi^2 * kin / (16*pi*mk^3) / P.hbar;

    case {'SxP', 'PxP'}   % pseudoscalar phi-quark coupling (pi0/eta poles)
        % theta_j from sppp_mixing.m (NLO; pi0, eta AND eta'), the source the
        % relic cross sections and kaon_invisible_ratio.m use.  The eta' vertex
        % (sqrt(6)/3)(2 + 3 g8s/g8)(mK^2 - mpi^2) is from eq. Klagrangian
        % (matters off the g_u+g_d+g_s = 0 direction).
        [~, thj] = sppp_mixing(mmed, [P.Cu P.Cd P.Cs]);
        Epi = (mk^2 + mpip^2 - mmed^2) / (2*mk);
        gphipik = -(2*mk*Epi - mk^2) * thj(1) ...
                  - (mk^2 - 6*mk*Epi + 2*mpip^2)/sqrt(3) * thj(2) ...
                  - sqrt(6)/3*(2 + 3*P.g8s/P.g8)*(mk^2 - mpip^2) * thj(3);
        Gamma = (P.GF*P.F*P.g8*P.Vus*P.Vstarud)^2 / (32*pi*mk^3) * gphipik^2 ...
                * kin / P.hbar;

    case 'AxA'            % K+ -> pi+ A'
        % From (2 i A'^mu)(-lak dK+ pi- + lap dpi- K+) in eq. Klagrangian:
        % M ~ (lak pK + lap ppi).eps*, and with sum_pol eps eps* = -g + q q/mA^2
        % only the part transverse to q = pK - ppi survives, (lak + lap) ppi_perp:
        %   sum|M|^2 = 2 GF^2 F^4 Vus^2 Vud^2 (lak + lap)^2 lambda / (4 mA^2),
        %   Gamma    = GF^2 F^4 Vus^2 Vud^2 (lak + lap)^2 lambda^(3/2) / (32 pi mA^2 mK^3).
        % The NWA residue of the 3-body K -> pi chi chibar rate
        % (kaon_invisible_ratio.m) reproduces this form.
        lak = P.g8*(P.Cs + P.Cu) + P.g8s*(P.Cu + P.Cd + P.Cs);
        lap = P.g8*(P.Cd + P.Cu) + P.g8s*(P.Cu + P.Cd + P.Cs);
        Gamma = P.GF^2 * P.F^4 * P.Vus^2 * P.Vud^2 * (lak + lap)^2 ...
                * kin^3 / (32*pi*mmed^2*mk^3) / P.hbar;

    otherwise
        Gamma = NaN;
    end
end

% ======================================================================
% Chance that an on-shell phi from K+ -> pi+ phi looks invisible at NA62:
%   P_inv = 1 - BR_vis * (1 - exp(-L/ell)),   ell = (p_phi/mphi) c / Gamma_tot
% i.e. it decays to chi chibar (or nu nubar), or it is still undecayed after L.  Visible =
% gamma gamma + meson pairs (pi pi opens at 2 m_pi0 for SxS/PxS).  Widths from
% phi_partial_widths.m, the relic code's source of truth.  Geometry: L = 100 m,
% typical K decay point to the LKr; E_phi = 50 GeV, the 75 GeV K+ beam less a
% ~25 GeV pi+ (middle of the NA62 signal region).
% ======================================================================
function Pinv = kphi_invisible_fraction(model, mx, mphi, gx, g)
    L_det = 100;              % [m]
    E_phi = 50;               % [GeV]
    hbarc_m = 1.973269804e-16;   % [GeV m]
    names = struct('SxS', 'heavyScalar', 'PxS', 'heavyScalarPS', ...
                   'SxP', 'heavyScalarSP', 'PxP', 'heavyScalarPP', 'AxA', 'axialVector');
    W = phi_partial_widths(mphi, mx, gx, g, names.(model));
    Ginv = W.xx;              % chi chibar is invisible
    Gvis = W.total - Ginv;
    if W.total <= 0 || Gvis <= 0
        Pinv = 1;   % stable, or chi chibar only
        return;
    end
    ell  = sqrt(E_phi^2 - mphi^2) / mphi * hbarc_m / W.total;   % lab decay length [m]
    % (Gxx + Gvis e^{-L/ell})/Gtot, not 1 - BR_vis(1 - e^{-L/ell}): the latter
    % cancels to exactly 0 once BR(chi chibar) < 1e-16.
    Pinv = (Ginv + Gvis * exp(-L_det / ell)) / W.total;
end

% ======================================================================
% Gamma(j -> chi chibar) [GeV] for j = {pi0, eta, eta'}  (1x3; NaN = closed /
% not applicable).  Couplings come from the relic pipeline's single sources of
% truth, so the bound sees exactly the mixing the cross sections do.
% ======================================================================
function G = meson_invisible_widths(model, mx, mphi, gx, g)
    mj = [0.13498, 0.547862, 0.95778];   % as sppp_mixing.m / axial_couplings.m
    b2 = 1 - 4*mx^2 ./ mj.^2;            % beta^2; channel open where > 0
    switch model
    case {'SxP', 'PxP'}
        % g_{j chi} = gx * theta_j, the same Lorentz structure as the phi-chi
        % vertex: pseudoscalar (PxP, beta^1) or scalar (SxP, beta^3) -- the
        % forms sppp_couplings.m already uses for Gamma(eta' -> chi chibar).
        [~, theta] = sppp_mixing(mphi, g);
        gj2 = (gx * theta).^2;
        if strcmp(model, 'PxP')
            G = gj2 .* mj .* sqrt(max(b2, 0)) / (8*pi);
        else
            G = gj2 .* mj .* max(b2, 0).^1.5 / (8*pi);
        end
    case 'AxA'
        % lambda_{j chi}^2 mx^2 sqrt(mj^2/4 - mx^2)/pi, as axial_couplings.m.
        % That function refuses points whose meson width exceeds its mass.
        % lambda_{j chi} is linear in gx, so take it from a probe below that
        % floor and rescale (as the Julia _axial_probe): the ratio stays finite
        % -- an Inf would read as "not applicable" downstream -- and enormous.
        lam_chi = [];  gp = gx;
        for k = 1:50
            try
                [~, lam_chi] = axial_couplings(mphi, mx, gp, g);
                lam_chi = lam_chi * (gx / gp);
                break;
            catch ME
                if ~strcmp(ME.identifier, 'axial_couplings:mediatorTooLight'), rethrow(ME); end
                gp = gp * 1e-6;
            end
        end
        if isempty(lam_chi), G = nan(1, 3);  return; end
        G = lam_chi.^2 * mx^2 .* sqrt(max(0, mj.^2/4 - mx^2)) / pi;
    otherwise   % SxS / PxS: scalar phi does not mix with the 0^- mesons
        G = nan(1, 3);  return;
    end
    G(b2 <= 0) = NaN;
end

% ======================================================================
% Low-velocity <sigma v> [cm^3/s] used by the CMB / INTEGRAL bounds.
% Analytic threshold (s = 4 m_chi^2) expressions.
% ======================================================================
function sv = sigmav_indirect(model, P)
    e = P.e;  alpha = P.alpha;  hbarc = P.hbarc;  cc = P.cc;
    Fpi0 = P.Fpi0;  F = P.F;  B0 = P.B0;
    mpi0 = P.mpi0;  meta = P.meta;  meta0 = P.meta0;
    mpiplus = P.mpiplus;  mkplus = P.mkplus;
    Gammapi0 = P.Gammapi0;  Gammaeta = P.Gammaeta;  Gammaeta0 = P.Gammaeta0;
    Cu = P.Cu;  Cd = P.Cd;  Cs = P.Cs;
    mx = P.mx;  mmed = P.mphi;  gx = P.gx;  mp = P.mp;
    sapprox = 4*mx^2;

    switch model
    % ---- Scalar x Scalar (INTEGRAL only) -------------------------------
    case 'SxS'
        Dnew = 0.5*(B0*(Cu + Cd)/mpiplus^2 + B0*(Cu + Cs)/mkplus^2)^2;
        jfac = (2.42e17)/(1.30e23);
        % 4*pi: the source expression is dsigma/dOmega.
        % hbarc^2*cc (ONE factor of c), as in every other branch: the bracket is
        % [GeV^-2] (D ~ GeV^-2, the mass ratio and the halo <v^2> jfac
        % dimensionless), so GeV^-2 -> cm^2 needs (hbar c)^2 and cm^2 -> cm^3/s
        % needs a single c.  Cf. the SxP branch below -- same p-wave structure,
        % same jfac, one cc.
        % /589824: e^4/(589824 pi^6) = alpha^2/(36864 pi^4), with the
        % identical-photon S_gg = 1/2 carried by every chi chibar -> gamma
        % gamma cross section.
        sv = 4*pi * (gx^2*e^4*Dnew)/(589824*pi^6) ...
             * (8*mx^4*jfac)/(4*mx^2 - mmed^2)^2 ...
             * hbarc^2 * cc;

    % ---- Pseudoscalar x Scalar (CMB + INTEGRAL) ------------------------
    case 'PxS'
        Dnew = 0.5*(B0*(Cu + Cd)/mpiplus^2 + B0*(Cu + Cs)/mkplus^2)^2;
        % gamma gamma + chi chibar + meson pairs
        GammaphiPS = phi_gg_xx_width(model, P) + meson_pair_width(mmed, [Cu Cd Cs]);
        sv = 4*pi*gx^2*alpha^2*Dnew/(36864*pi^4) * (4*mx^2)^2 ...
             / ((4*mx^2 - mmed^2)^2 + (mmed*GammaphiPS)^2) * 2 * hbarc^2 * cc;

    % ---- Scalar x Pseudoscalar (INTEGRAL only) -------------------------
    case 'SxP'
        [mj, gcx, ggx, Gam] = sppp_poles_indirect(P, false);
        jfac = (2.42e17)/(1.30e23);            % halo <v^2> (p-wave rate)
        amp = sum(gcx .* ggx ./ (sapprox - mj.^2 + 1i*mj.*Gam));
        % 4*pi: the source expression is dsigma/dOmega.
        sv = 4*pi * sapprox^(3/2) * mx * jfac / (512*pi^2) ...
             * abs(amp)^2 * hbarc^2 * cc;

    % ---- Pseudoscalar x Pseudoscalar (CMB + INTEGRAL) ------------------
    case 'PxP'
        [mj, gcx, ggx, Gam] = sppp_poles_indirect(P, true);
        amp = sum(gcx .* ggx ./ (sapprox - mj.^2 + 1i*mj.*Gam));
        % 4*pi: the source expression is dsigma/dOmega.
        sv = 4*pi * sapprox^2 / (512*pi^2) * abs(amp)^2 * 2 * hbarc^2 * cc;

    % ---- Axial x Axial (CMB + INTEGRAL) --------------------------------
    case 'AxA'
        Gammaetaprimegammagamma = (2.22/100)*0.196*P.MeV;
        Fetaprime = sqrt((alpha^2*P.metaprime^3)/(64*pi^3*Gammaetaprimegammagamma));
        lambdapi0gamma  = e^2/(4*pi^2*F);
        % Data-driven F_eta (paper Eq. lambdagamma); cf. axial_couplings.m.
        % sppp_poles_indirect below uses the LO 1/sqrt(3), as in the paper's
        % pseudoscalar section.
        Feta            = sqrt((alpha^2*P.meta^3)/(64*pi^3*5.16e-7));
        lambdaetagamma  = e^2/(4*pi^2*Feta);
        lambdaeta0gamma = e^2/(4*pi^2*Fetaprime);   % 4*pi^2, as pi0/eta above
        [mpi0Ap, metaAp, meta0Ap] = axial_meson_Aprime(P);
        lambdapichi   = gx*mpi0Ap/mmed^2;
        lambdaetachi  = gx*metaAp/mmed^2;
        lambdaeta0chi = gx*meta0Ap/mmed^2;
        amp = lambdapichi*lambdapi0gamma/(sapprox - mpi0^2 + 1i*mpi0*Gammapi0) ...
            + lambdaetachi*lambdaetagamma/(sapprox - meta^2 + 1i*meta*Gammaeta) ...
            + lambdaeta0chi*lambdaeta0gamma/(sapprox - meta0^2 + 1i*meta0*Gammaeta0);
        % 4*pi: the source expression is dsigma/dOmega.
        sv = 4*pi * mx^2 * sapprox^2 / (128*pi^2) * abs(amp)^2 * 2 * hbarc^2 * cc;

    otherwise
        sv = NaN;
    end
end

% ---- shared coupling / width helpers ---------------------------------
function G = phi_gg_xx_width(model, P)
% Gamma(phi -> gamma gamma) + Gamma(phi -> chi chibar), scalar phi-SM vertex:
% SxS beta^3, PxS beta^1 chi chibar width.
    Dnew = 0.5*(P.B0*(P.Cu + P.Cd)/P.mpiplus^2 + P.B0*(P.Cu + P.Cs)/P.mkplus^2)^2;
    G = (Dnew*P.e^4*P.mphi^3)/(73728*pi^5);
    if P.mphi^2 > 4*P.mx^2
        if strcmp(model, 'SxS')
            G = G + P.gx^2*(P.mphi^2 - 4*P.mx^2)^(3/2)/(8*pi*P.mphi^2);
        else
            G = G + P.gx^2*sqrt(P.mphi^2 - 4*P.mx^2)/(8*pi);
        end
    end
end

function Gm = meson_pair_width(mphi, g)
% Gamma(phi -> meson pairs) of the scalar phi-SM vertex (SxS/PxS) [GeV], from
% phi_partial_widths.m / meson_pairs.m (the relic code's source of truth).
    W  = phi_partial_widths(mphi, 0, 0, g, 'heavyScalar');
    Gm = W.pipi + W.KK + W.eta;
end

function [mj, gcx, ggx, Gam] = sppp_poles_indirect(P, isPP)
% s-channel pole set for the SxP / PxP indirect <sigma v>, mirroring
% matlab/sppp_couplings.m (the relic pipeline's single source of truth) so the
% two pipelines describe the same physics.
%
%   generic (g not all equal):  {pi0, eta, phi'}
%   SU(3) singlet (Cu=Cd=Cs):   {eta', phi'}
%
% In the singlet direction the pi0 and eta mixing couplings (~ Cu-Cd and
% Cu+Cd-2Cs) AND the induced phi'gamma coupling all vanish, so the generic
% amplitude is identically zero and chi would decouple from photons.  The
% surviving channel is the SU(3)-singlet eta' (~ Cu+Cd+Cs), normally dropped as
% the heaviest meson -- see sppp_couplings.m.
%
% Returned in the caller's convention amp = sum gcx.*ggx./(s - mj.^2 + i mj Gam);
% |amp|^2 is invariant under the overall sign relative to sppp_couplings.m's
% 1/(mj^2 - s - i mj Gam), so the two agree term by term.
    e = P.e;  alpha = P.alpha;  F = P.F;  Fpi0 = P.Fpi0;  B0 = P.B0;
    mpi0 = P.mpi0;  meta = P.meta;  metap = P.metaprime;
    mmed = P.mphi;  mx = P.mx;  gx = P.gx;
    Cu = P.Cu;  Cd = P.Cd;  Cs = P.Cs;
    G_etap = 1.88e-4;      % eta' SM total width [GeV]   (sppp_couplings.m)
    G_etap_gg = 4.34e-6;   % eta' -> gamma gamma partial [GeV]

    if Cu == Cd && Cd == Cs
        % --- SU(3) singlet: {eta', phi'} --------------------------------
        F_etap      = sqrt(alpha^2 * metap^3 / (64*pi^3 * G_etap_gg));
        g_etapgamma = alpha / (pi*F_etap);        % data-driven eta'->gg (1/pi: see sppp_couplings.m)
        % eta' singlet mixing mass: sqrt(6)/3 * sum(g) * (B0 F + (L6+L7) NLO), as in sppp_mixing.m
        m_mix       = sqrt(6)/3 * (Cu + Cd + Cs) ...
                      * (B0*F + (P.L6 + P.L7)*32*B0^2/F*(P.mqu + P.mqd + P.mqs));
        denom       = (mmed^2 - metap^2);
        g_etapx        = -gx * m_mix / denom;
        gphiprimex     =  gx;
        gphiprimegamma =  g_etapgamma * m_mix / denom;   % induced phi'->gg
        if isPP   % pseudoscalar chi vertex, beta^1
            Gam_etapxx = g_etapx^2 * sqrt(max(0, metap^2/4 - mx^2)) / (4*pi);
        else      % SxP: scalar chi vertex, beta^3
            Gam_etapxx = g_etapx^2 * max(0, metap^2 - 4*mx^2)^(3/2) / (8*pi*metap^2);
        end
        Gammaphi = med_width(gphiprimex, gphiprimegamma, mmed, mx, isPP) ...
                   + phi_3pi_width(mmed, [Cu Cd Cs]);
        mj  = [metap,         mmed          ];
        gcx = [g_etapx,       gphiprimex    ];
        ggx = [g_etapgamma,   gphiprimegamma];
        Gam = [G_etap + Gam_etapxx, Gammaphi ];
    else
        % --- generic: {pi0, eta, phi'} -----------------------------------
        [mpi0phi2, metaphi2] = pp_meson_phi_couplings(P);
        gphichi = gx;
        gpichi  = gx * mpi0phi2 / (mmed^2 - mpi0^2);
        getachi = gx * metaphi2 / (mmed^2 - meta^2);
        ggammapi  = e^2 / (4*pi^2*Fpi0);
        ggammaeta = e^2 / (4*sqrt(3)*pi^2*Fpi0);
        ggammaphi = -ggammapi*mpi0phi2/(mmed^2 - mpi0^2) ...
                    - ggammaeta*metaphi2/(mmed^2 - meta^2);
        Gammaphi = med_width(gphichi, ggammaphi, mmed, mx, isPP) ...
                   + phi_3pi_width(mmed, [Cu Cd Cs]);
        mj  = [mpi0,        meta,        mmed     ];
        gcx = [gpichi,      getachi,     gphichi  ];
        ggx = [ggammapi,    ggammaeta,   ggammaphi];
        Gam = [P.Gammapi0,  P.Gammaeta,  Gammaphi ];
    end
end

function G = phi_3pi_width(mphi, g)
% Gamma(phi -> 3pi), as in phi_partial_widths.m; 0 below 3 mpi0 or with
% LOOP_3PI off (threepi_weight returns 0).
    G = 0;
    if mphi > 3*0.13498
        G = threepi_weight(mphi^2, g) / (64*pi^3*mphi);
    end
end

function sv = threebody_cmb_sv(model, P)
% Energy-weighted s-wave <sigma v> [cm^3/s] of the three-body channels at
% s = 4 mx^2 (v_rel = 2 beta), added to the gamma gamma <sigma v> for the CMB
% only.  P(x)P: chi chi -> 3pi and pi pi phi; AxA: the J = 0 chi chi -> 3pi.
% S(x)P and the axial J = 1 part vanish at threshold; the INTEGRAL curve is a
% gamma-line bound.  The Slatyer curve is for gamma gamma (all 2mx deposited),
% so each channel carries its deposited-energy fraction, f_eff(e) ~ f_eff(gamma):
%   3pi       0.52 = (1 + 2 x 0.275)/3, pi+- -> e+- keeps 0.787 x 7/20 = 0.275
%   pi pi phi 0.18 = 2/3 x 0.275, the pions only (phi may decay after recombination)
% Lower bounds on the deposited energy, so the bound is conservative.  Same
% numbers as analysis/mesonandCMB.py.
    kappa_3pi = 0.52;  kappa_pipiphi = 0.18;
    mx = P.mx;  mmed = P.mphi;  gx = P.gx;  g = [P.Cu P.Cd P.Cs];
    s0 = 4*mx^2;
    sv = 0;
    switch model
    case 'PxP'
        [~, ~, ~, Gam] = sppp_poles_indirect(P, true);    % Gam(end): the mediator
        D2 = (s0 - mmed^2)^2 + mmed^2*Gam(end)^2;
        sv = kappa_3pi*gx^2*threepi_weight(s0, g)/(64*pi^3*D2) ...
           + kappa_pipiphi*gx^2*pipiphi_weight(s0, mmed, g)/(64*pi^3);
    case 'AxA'
        [~, WL] = axial3pi_weight(s0, g);
        sv = kappa_3pi*gx^2*mx^2*WL/(16*pi^3*s0^2*mmed^4);
    end
    if ~isfinite(sv), sv = 0; end     % must not erase the gamma gamma bound
    sv = sv * P.hbarc^2 * P.cc;
end

function G = med_width(gphichi, ggammaphi, mphi, mx, isPP)
% Mediator total width for the SxP (scalar chi vertex) / PxP (pseudoscalar) cases.
    if isPP
        G = pseudo_chi_phi_width(gphichi, ggammaphi, mphi, mx);
    else
        G = scalar_chi_phi_width(gphichi, ggammaphi, mphi, mx);
    end
end

function [mpi0phi2, metaphi2] = pp_meson_phi_couplings(P)
% phi-meson^2 couplings for the pseudoscalar phi-quark vertex (L6/L7-corrected).
    B0 = P.B0;  F = P.F;  Cu = P.Cu;  Cd = P.Cd;  Cs = P.Cs;
    mu = P.mu;  md = P.md;  ms = P.ms;
    mpi0phi2 = (Cu - Cd)*(B0*F + 32*B0^2*P.L6/F*(mu + md + ms)) ...
               + 32*B0^2*P.L7/F*(Cu + Cd + Cs)*(mu - md);
    % eta = 1/sqrt(3) x the pi0 bracket (chiral trace; see sppp_couplings.m).
    metaphi2 = ( (Cu + Cd - 2*Cs)*(B0*F + 32*B0^2*P.L6/F*(mu + md + ms)) ...
               + 32*B0^2*P.L7/F*(Cu + Cd + Cs)*(mu + md - 2*ms) ) / sqrt(3);
end

function [mpi0Ap, metaAp, meta0Ap] = axial_meson_Aprime(P)
% A'-meson mixing amplitudes for the axial-vector vertex.
    B0 = P.B0;  F = P.F;  Cu = P.Cu;  Cd = P.Cd;  Cs = P.Cs;
    mu = P.mu;  md = P.md;  ms = P.ms;  mqu = P.mqu;  mqd = P.mqd;  mqs = P.mqs;
    mpi0Ap  = F*((Cu - Cd)*(1 + (16*B0*P.L4)/F^2*(mu + md + ms)) ...
                 + (16*B0*P.L5*(Cu*mu - Cd*md))/F^2);
    metaAp  = (sqrt(3)*F)/3*((Cu + Cd - 2*Cs)*(1 + (16*B0*P.L4)/F^2*(mqu + mqd + mqs)) ...
                 + 1/F^2 * 16*B0*P.L5*(Cu*mqu + Cd*mqd - 2*Cs*mqs));
    meta0Ap = (sqrt(6)*F)/3*((Cu + Cd + Cs)*(1 + (16*B0*P.L4)/F^2*(mqu + mqd + mqs)) ...
                 + 1/F^2 * 16*B0*P.L5*(Cu*mqu + Cd*mqd + Cs*mqs));
end

function G = scalar_chi_phi_width(gphichi, ggammaphi, mphi, mx)
% Total phi width, scalar phi-chi vertex (SxP): chi width ~ (mphi^2-4mx^2)^(3/2).
    Ggg = ggammaphi^2*mphi^3/(64*pi);
    if mphi^2 >= 4*mx^2
        G = gphichi^2*(mphi^2 - 4*mx^2)^(3/2)/(8*pi*mphi^2) + Ggg;
    else
        G = Ggg;
    end
end

function G = pseudo_chi_phi_width(gphichi, ggammaphi, mphi, mx)
% Total phi width, pseudoscalar phi-chi vertex (PxP): chi width ~ sqrt(mphi^2/4-mx^2).
    Ggg = ggammaphi^2*mphi^3/(64*pi);
    if mphi^2 >= 4*mx^2
        G = gphichi^2*sqrt(mphi^2/4 - mx^2)/(4*pi) + Ggg;
    else
        G = Ggg;
    end
end

% ---- inlined digitized limit curves -------------------------------------
function lim = curve_slatyer(mx)
% Slatyer CMB <sigma v> upper bound [cm^3/s] at DM mass mx.  Data: log10(mx/GeV)
% vs log10(<sigma v>/(cm^3 s^-1)).  NaN (not constrained) outside the data range.
    x = [-5.299492385786800, -5.0152284263959400, -4.770123277737490, ...
         -4.501812907904280, -4.171138506163890, -3.908629441624370, ...
         -3.537345902828140, -3.09934735315446, -2.973168963016680, ...
         -2.742567077592460, -2.4829586656997800, -2.265409717186370, ...
         -2.029006526468460, -1.789702683103700, -1.361856417693980, ...
         -0.9905728788977540, -0.45250181290790600, -0.014503263234229500, ...
          0.2871646120377070, 0.6932559825960830];
    y = [-32.62642947903430, -32.27573062261750, -31.99618805590850, ...
         -31.726810673443500, -31.38627700127070, -31.058449809402800, ...
         -30.49682337992380, -29.73189326556540, -29.513341804320200, ...
         -29.337992376111800, -29.363405336721700, -29.320203303684900, ...
         -29.19059720457430, -29.007623888183000, -28.606099110546400, ...
         -28.21982210927570, -27.57941550190600, -26.97712833545110, ...
         -26.601016518424400, -26.260482846251600];
    lim = 10.^interp1(x, y, log10(mx), 'linear', NaN);
end

function lim = curve_integral(mx)
% INTEGRAL telescope <sigma v> upper bound [cm^3/s] at DM mass mx.  NaN outside
% the data range (matches the reference's fill -> 0 exclusion, i.e. not excluded).
    x = [-1.5702280912364900, -1.312124849939980, -1.3085234093637500, ...
         -1.0036014405762300, -0.9987995198079230, -0.701080432172869, ...
         -0.6974789915966390, -0.22448979591836700, -0.21848739495798300, ...
          0.2557022809123650] - 3;
    y = [-32.47761194029850, -31.970149253731300, -32.471641791044800, ...
         -31.86268656716420, -31.907462686567200, -31.31641791044780, ...
         -31.01492537313430, -30.074626865671600, -30.555223880597000, ...
         -29.608955223880600];
    lim = 10.^interp1(x, y, log10(mx), 'linear', NaN);
end

% NA62 K+ -> pi+ nu nu acceptance in m_miss^2 = s [GeV^2]: R1 = [0, 0.010], R2 =
% [0.026, 0.068], edges smeared by the resolution 0.001.  Copy of the one in
% kaon_invisible_ratio.m -- keep them identical.
function A = na62_acceptance(s)
    sig = 1e-3;  Phi = @(x) 0.5*erfc(-x/sqrt(2));
    A = Phi((0.010 - s)/sig) - Phi((0 - s)/sig) + Phi((0.068 - s)/sig) - Phi((0.026 - s)/sig);
end
