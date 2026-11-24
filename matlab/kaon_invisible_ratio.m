function r = kaon_invisible_ratio(mx, mphi, gx, g, interaction)
%KAON_INVISIBLE_RATIO  Gamma(K+ -> pi+ chi chi) / K+->pi+ invisible limit.
%   r = KAON_INVISIBLE_RATIO(MX, MPHI, GX, G, INTERACTION) returns the ratio of
%   the predicted rare-kaon decay rate Gamma(K+ -> pi+ chi chi) to the
%   experimental K+ -> pi+ + (invisible) branching-ratio limit (derived from
%   the K+ -> pi+ nu nubar bound).  r > 1 => the point is excluded.
%
%   Inputs mirror the relic-code convention:
%     MX          DM mass m_chi   [GeV]
%     MPHI        mediator mass m_phi (m_A' for axialVector)   [GeV]
%     GX          phi-chi coupling g_chi (lambda_a for axialVector)
%     G           [g_u g_d g_s] phi-quark couplings (C_u, C_d, C_s)
%     INTERACTION relic model string (heavyScalar*, lightScalar*, axialVector)
%
%   The kaon channel is open only for MX < (m_K+- - m_pi+-)/2 ~ 0.177 GeV;
%   r = NaN when the channel is closed or the interaction has no meson-decay
%   implementation.
%
%   DO NOT rescale this ratio as GX^2.  That holds only while the phi
%   propagator stays off shell.  Once the pole enters the Dalitz region,
%
%       2 MX < MPHI < m_K+- - m_pi+- ~ 0.354 GeV,
%
%   the E_pi integral is saturated by the resonance, so F ~ 1/Gamma_phi; since
%   Gamma(phi -> chi chibar) ~ GX^2 that cancels the explicit GX^2 and the rate
%   plateaus, becoming GX-independent.  At MPHI = 0.1 GeV the naive estimate
%   1/sqrt(r(GX=1)) is then 60-240x too weak in the coupling.  Callers wanting
%   the excluded coupling must root-find on r(GX) = 1 -- see
%   ratio_crit_coupling.m, used by compare_limits_g.m.  Above MPHI/2, and for
%   MPHI > m_K+- - m_pi+- (pole never on shell), the GX^2 rescale is exact and
%   the root-find reproduces it to <= 0.1%.
%
%   Interaction -> model map (as analysis/kaon_invisible.py):
%     heavyScalar   / lightScalar   / scalar     -> SxS  (scalar chi, scalar quark)
%     heavyScalarPS / lightScalarPS / scalarPS   -> PxS  (pseudoscalar chi, scalar quark)
%     heavyScalarSP / lightScalarSP / scalarSP   -> SxP  (scalar chi, pseudoscalar quark)
%     heavyScalarPP / lightScalarPP / scalarPP   -> PxP  (both pseudoscalar)
%     axialVector                                -> AxA  (axial-vector A')
%
%   Resonance handling (SxP/PxP):
%     * Both carry the pi0/eta Breit-Wigner widths (PxP's were added so its
%       propagators match SxP's -- only the scalar-chi -4 mx^2 phase-space factor
%       distinguishes them).  The on-shell pi0 resonance at E_pi ~ 0.248 GeV sits
%       INSIDE the integration range and is razor-thin (Gamma_pi0 ~ 8 eV), so it
%       is handled by an NWA notch: excised from the numerical continuum and added
%       back analytically, giving an integrator-independent result.
%   AxA follows paper eq. (kpichichiA) (signed lambda_jpiK, A' width and
%   propagators, lambda_A'K^2); see the branch.
%   Mediator widths (SxP/PxP/AxA) are the relic code's totals, incl. -> chi chibar.

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
    MeV = 1e-3;  eV = 1e-9;
    e    = 0.3;
    B0   = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV (was 2.6)
    Fpi0 = 92.4  * MeV;   F = Fpi0;
    mpi  = 134.98 * MeV;  mpi0 = mpi;      % neutral pion
    meta = 547.86 * MeV;                   % eta
    mpiplus = 139.57 * MeV;                % charged pion
    mkplus  = 493.68 * MeV;  mKplus = mkplus;   % charged kaon
    meta0   = 958 * MeV;                   % eta'
    GF   = 1.17e-5;
    Vus  = 0.222;  Vusstar = Vus;  Vstarud = 0.974;  Vud = Vstarud;
    g8   = 3.62;   g8s = -0.38 * g8;
    L4   = 1.87e-3;  L5 = 1.22e-3;
    mu = 2.2e-3;  md = 4.7e-3;  ms = 95e-3;
    mqu = mu;  mqd = md;  mqs = ms;
    % hbar from hbarc/c (same hbarc as alltheconstants.m: 1.97e-14 GeV*cm)
    hbarc_cgs = 197 * MeV * 1e-13;         % [GeV*cm]
    cc        = 3e10;                      % [cm/s]
    hbar      = hbarc_cgs / cc;            % [GeV*s]

    Cu = g(1);  Cd = g(2);  Cs = g(3);
    lambdaa = gx;

    % K+ invisible-decay limit as a rate [s^-1]
    lifetimekplus  = 1.238e-8;             % [s]
    decayratekplus = 1 / lifetimekplus;
    kplustopinunubranchingratio = 1.14e-10;
    kplusinvlimit  = decayratekplus * kplustopinunubranchingratio;

    % Shared kinematics
    Epimax = @(mchi) (mkplus^2 + mpiplus^2 - 4*mchi^2) / (2*mkplus);
    Yfun   = @(Epi, mchi) sqrt( ((Epi.^2 - mpiplus^2) .* ...
                (mkplus^2 + mpiplus^2 - 2*Epi*mkplus - 4*mchi^2)) ./ ...
                (mkplus^2 + mpiplus^2 - 2*Epi*mkplus) );

    % Channel closed (E_pi^max <= m_pi+-, i.e. mx >= (mK - mpi)/2): no constraint.
    if Epimax(mx) <= mpiplus
        r = NaN;
        return;
    end
    % NA62 acceptance in m_miss^2 = s: signal regions R1 =
    % [0, 0.010] and R2 = [0.026, 0.068] GeV^2, smeared by the m_miss^2 resolution
    % sigma = 0.001 GeV^2 (JHEP 06 (2021) 093).  The gap holds the K -> pi pi0 peak
    % (s = 0.0182, 8 sigma inside: weight ~1e-16; pi0 -> invisible is mesoninv's)
    % and s > 0.068 the K -> 3 pi background.  Every rate below carries A(s); an
    % on-shell pole on a region edge counts half.  bpE: breakpoints at the edges.
    Aacc = @(s) na62_acceptance(s);
    bpE  = (mkplus^2 + mpiplus^2 - na62_edges())/(2*mkplus);
    bpE  = bpE(bpE > mpiplus & bpE < Epimax(mx));

    switch model
    % =====================================================================
    case 'SxS'   % scalar phi-chi, scalar phi-quark: s-channel phi (BW)
    % =====================================================================
        delta = -GF * F^2 * Vusstar * Vud * g8 / sqrt(2);
        Dnew  = 0.5*(B0*(Cu + Cd)/mpiplus^2 + B0*(Cu + Cs)/mkplus^2)^2;
        if mphi^2 > 4*mx^2
            Gphi = (Dnew*e^4*mphi^3)/(73728*pi^5) + (gx^2*(mphi^2 - 4*mx^2)^(3/2))/(8*pi*mphi^2);
        else
            Gphi = (Dnew*e^4*mphi^3)/(73728*pi^5);
        end
        % + phi -> pi pi, K K, eta eta, pi0 eta (open above mphi = 2 m_pi0):
        % they drain the resonance as gamma gamma and chi chibar do.
        Gphi = Gphi + meson_pair_width(mphi, g);
        % Gamma_phi/mphi can be ~1e-12 (small gx, g): notch the on-shell pole
        % (2 mx < mphi < mK - mpi+) rather than hand it to the adaptive rule.
        numf = @(Epi) (mKplus^2 + mpiplus^2 - 2*mKplus*Epi - 4*mx^2) .* Yfun(Epi, mx);
        numfA = @(Epi) numf(Epi) .* Aacc(mKplus^2 + mpiplus^2 - 2*mKplus*Epi);
        FSS = bw_notch_integral(numfA, mphi, Gphi, mpiplus, Epimax(mx), mKplus, mpiplus, bpE);
        gpiKphi = delta*B0*((Cs + Cu)*mpiplus^2 - (Cu + Cd)*mKplus^2)/(mKplus^2 - mpiplus^2);
        Gamma = 1/(32*pi^3*mKplus) * gx^2 * gpiKphi^2 * FSS / hbar;
        r = Gamma / kplusinvlimit;

    % =====================================================================
    case 'PxS'   % pseudoscalar phi-chi, scalar phi-quark: s-channel phi (BW)
    % =====================================================================
        delta = -GF * F^2 * Vusstar * Vud * g8 / sqrt(2);
        Dnew  = 0.5*(B0*(Cu + Cd)/mpiplus^2 + B0*(Cu + Cs)/mkplus^2)^2;
        if mphi^2 > 4*mx^2
            Gphi = (Dnew*e^4*mphi^3)/(73728*pi^5) + (gx^2*(mphi^2 - 4*mx^2)^(1/2))/(8*pi);
        else
            Gphi = (Dnew*e^4*mphi^3)/(73728*pi^5);
        end
        Gphi = Gphi + meson_pair_width(mphi, g);   % as SxS
        % Gamma_phi/mphi can be ~1e-12 (small gx, g): notch the on-shell pole
        % (2 mx < mphi < mK - mpi+) rather than hand it to the adaptive rule.
        numf = @(Epi) (mKplus^2 + mpiplus^2 - 2*mKplus*Epi) .* Yfun(Epi, mx);
        numfA = @(Epi) numf(Epi) .* Aacc(mKplus^2 + mpiplus^2 - 2*mKplus*Epi);
        FPS = bw_notch_integral(numfA, mphi, Gphi, mpiplus, Epimax(mx), mKplus, mpiplus, bpE);
        gpiKphi = delta*B0*((Cs + Cu)*mpiplus^2 - (Cu + Cd)*mKplus^2)/(mKplus^2 - mpiplus^2);
        Gamma = 1/(32*pi^3*mKplus) * gx^2 * gpiKphi^2 * FPS / hbar;
        r = Gamma / kplusinvlimit;

    % =====================================================================
    case {'PxP', 'SxP'}   % pseudoscalar phi-quark: pi0/eta/eta' mixing poles
    % =====================================================================
    % K+ -> pi+ j* (j = pi0, eta, eta'), j mixes into phi, phi -> chi chibar.
    % Mixing masses m^2_{j phi} and angles theta_j = m^2_{j phi}/(mphi^2 - mj^2)
    % come from sppp_mixing.m (NLO: L6, L7, and the eta'), the same source as the
    % relic cross sections (sppp_couplings.m) and kphi.  g_{j chi} = +gx theta_j
    % for every j, so each meson pole and its share of the phi pole enter as
    % theta_j g_{j pi K} [1/(s - mj^2) - 1/(s - mphi^2)].
    % K -> pi j vertices from the weak chiral Lagrangian (eq. Klagrangian), in
    % units of g8, with s = (pK - pi)^2:
    %   pi0:  mpi^2 - s = 2 mK Epi - mK^2
    %   eta:  (3 s - 2 mK^2 - mpi^2)/sqrt(3) = (mK^2 - 6 mK Epi + 2 mpi^2)/sqrt(3)
    %   eta': (sqrt(6)/3)(2 + 3 g8s/g8)(mK^2 - mpi^2)
        Gammapi  = 7.8 * eV;
        Gammaeta = 1.31 * 1000 * eV;
        Gammaetaprime = 0.196 * MeV;
        [~, thj] = sppp_mixing(mphi, g);
        gpichi   = gx * thj(1);
        getachi  = gx * thj(2);
        getapchi = gx * thj(3);
        gphichi  = gx;
        gpipik   = @(Epi) 2*mkplus*Epi - mkplus^2;
        getapik  = @(Epi) (mkplus^2 - 6*mkplus*Epi + 2*mpiplus^2)/sqrt(3);
        getappik = @(Epi) sqrt(6)/3*(2 + 3*g8s/g8)*(mkplus^2 - mpiplus^2) + 0*Epi;
        gphipik  = @(Epi) -thj(1)*gpipik(Epi) - thj(2)*getapik(Epi) - thj(3)*getappik(Epi);
        prefac = (GF*F*g8*Vus*Vstarud)^2/(64*pi^3*mkplus);
        % Both PxP and SxP now carry the pi0/eta widths (regularized propagators);
        % they differ ONLY in the scalar-chi phase-space factor (-4 mx^2, SxP).
        % The MEDIATOR propagator is BW-regulated too: s_chichi runs over
        % [4 mx^2, (mK-mpi+)^2], so the phi pole is ON SHELL inside the
        % integration range whenever 2 mx < mphi < mK-mpi+ ~ 0.354 GeV, and a
        % width-less 1/(s-mphi^2) makes the integrand singular there.  Gphi_reg
        % is the mediator's TOTAL width from the relic code (phi -> chi chibar,
        % gamma gamma, 3pi), so the on-shell rate is Gamma(K -> pi phi) *
        % BR(phi -> chi chibar), as in SxS/PxS.  r is not ~ gx^2 there:
        % callers root-find (ratio_crit_coupling).
        Wphi     = phi_partial_widths(mphi, mx, gx, g, interaction);
        Gphi_reg = Wphi.total;
        s_of = @(Epi) mkplus^2 + mpiplus^2 - 2*Epi*mkplus;             % s_chichi
        if strcmp(model, 'SxP'), ps_sub = 4*mx^2; else, ps_sub = 0; end
        PSf  = @(Epi) mkplus^2 + mpiplus^2 - 2*mkplus*Epi - ps_sub;
        amp  = @(Epi) (gpichi.*gpipik(Epi))  ./ (s_of(Epi) - mpi0^2 + 1i*mpi0*Gammapi) ...
                    + (gphichi.*gphipik(Epi)) ./ (s_of(Epi) - mphi^2 + 1i*mphi*Gphi_reg) ...
                    + (getachi.*getapik(Epi)) ./ (s_of(Epi) - meta^2 + 1i*meta*Gammaeta) ...
                    + (getapchi.*getappik(Epi)) ./ (s_of(Epi) - meta0^2 + 1i*meta0*Gammaetaprime);
        integrand = @(Epi) prefac .* PSf(Epi) .* Yfun(Epi, mx) .* abs(amp(Epi)).^2 .* Aacc(s_of(Epi));
        % NWA-notch the in-range on-shell poles (pi0, eta AND the mediator):
        % excise a window around each and add its resonance analytically,
        % freezing the smooth factors at the pole.  Over a notch [lo, hi],
        %   int dEpi 1/((s-mj^2)^2 + wj^2)
        %       = (atan(u_lo/wj) - atan(u_hi/wj))/(2 mkplus wj),  u = s(Epi)-mj^2,
        % which reduces to atan(2 mkplus dEpi/wj)/(mkplus wj) for a symmetric
        % window and stays exact when a window is clipped by the range ends.
        % Continuum(minus notches) + NWA is therefore the EXACT resolved
        % integral -- the razor-thin pi0 BW (Gamma_pi0 ~ 8 eV) and all -- and is
        % integrator-independent (mirrors thermal_pole_integral.m's convention).
        % The NWA is used when its notch is a small slice of the range (so the
        % frozen smooth factors are constant across it); broader poles are left
        % to the adaptive integrator, which gets the pole as a breakpoint.
        emax = Epimax(mx);
        mjv  = [mpi0, meta, mphi];  Gjv = [Gammapi, Gammaeta, Gphi_reg];
        njf  = {@(Epi) gpichi.*gpipik(Epi), @(Epi) getachi.*getapik(Epi), ...
                @(Epi) gphichi.*gphipik(Epi)};
        Kwid = 1e3;  max_notch_frac = 0.05;  min_notch_frac = 1e-7;
        notch = zeros(0, 2);  Inwa = 0;  brk = [];
        for jj = 1:numel(mjv)
            mj = mjv(jj);  wj = mj*Gjv(jj);
            if wj <= 0, continue; end
            Epij = (mkplus^2 + mpiplus^2 - mj^2)/(2*mkplus);
            if ~(Epij > mpiplus && Epij < emax), continue; end
            brk(end+1) = Epij; %#ok<AGROW>
            dEpi = Kwid * wj/(2*mkplus);
            if dEpi > max_notch_frac*(emax - mpiplus), continue; end
            % Gphi_reg can be tiny (small gx and g), so the phi notch can shrink
            % below eps(Epij) ~ 5.6e-17.  Epij +- dEpi then BOTH round to Epij:
            % the window excludes nothing (the continuum integrator is handed the
            % bare pole -> "minimum step size" warnings) and, worse, u_lo and u_hi
            % differenced from it come out EQUAL, so atan(u_lo/wj)-atan(u_hi/wj)
            % silently zeroes the resonance -- which is the dominant term.  Floor
            % the window at a representable width, and get the offsets from the
            % half-widths instead of by differencing: s_of is linear and
            % s_of(Epij) = mj^2 exactly, so u = -2*mkplus*(Epi - Epij) is exact.
            % The NWA stays valid because widening to a few ulp cannot move the
            % frozen smooth factors, and the analytic atan form is exact for any
            % window.  The floor is min_notch_frac of the range, not a few ulp:
            % a few-ulp notch is representable but the panels flanking it are
            % then ~8 ulp wide, which the adaptive integrator cannot subdivide
            % (that is the "minimum step size" warning).  At 1e-9 of the range
            % the notch is ~1e-10 GeV -- millions of ulp, easily resolved -- and
            % the pieces the NWA drops inside it (the other poles' continuum and
            % the interference) are ~1e-7 of their totals.  r is unchanged to 10
            % significant figures for floors from 1e-7 to 1e-11.
            dEpi = max([dEpi, 8*eps(Epij), min_notch_frac*(emax - mpiplus)]);
            d_lo = min(dEpi, Epij - mpiplus);  d_hi = min(dEpi, emax - Epij);
            lo_n = Epij - d_lo;  hi_n = Epij + d_hi;
            notch = [notch; lo_n, hi_n]; %#ok<AGROW>
            u_lo = 2*mkplus*d_lo;  u_hi = -2*mkplus*d_hi;
            nj = njf{jj}(Epij);
            Inwa = Inwa + prefac * PSf(Epij) * Yfun(Epij, mx) * abs(nj)^2 * Aacc(mj^2) ...
                          * (atan(u_lo/wj) - atan(u_hi/wj)) / (2*mkplus*wj);
        end
        % Continuum over [mpiplus, emax] minus the notch windows.  Outside a
        % notch the integrand still falls as 1/u^2 over a scale set by the notch
        % half-width, which the adaptive integrator cannot find on its own once
        % that width is tiny; grade breakpoints outward from each notch edge so
        % every decade of the falloff gets its own subinterval.
        grad = [];
        for n = 1:size(notch, 1)
            c = 0.5*(notch(n,1) + notch(n,2));  h = 0.5*(notch(n,2) - notch(n,1));
            % enough octaves to reach the interval ends: h can be ~1e-16 while
            % the range is ~0.13, so a fixed octave count leaves most of the
            % falloff unresolved and the warnings persist.
            f = 2.^(1:min(80, ceil(log2((emax - mpiplus)/max(h, realmin)))));
            grad = [grad; (c + h*f)'; (c - h*f)']; %#ok<AGROW>
        end
        if ~isempty(grad), grad = grad(grad > mpiplus & grad < emax); end
        bnds  = unique(sort([mpiplus; emax; bpE(:); notch(:); brk(:); grad(:)]));
        Icont = 0;
        for ii = 1:numel(bnds)-1
            a = bnds(ii);  b = bnds(ii+1);
            if b <= a, continue; end
            mid = 0.5*(a + b);  skip = false;
            for n = 1:size(notch, 1)
                if mid > notch(n,1) && mid < notch(n,2), skip = true; break; end
            end
            if skip, continue; end
            if na62_dead(a, b, mkplus, mpiplus), continue; end   % A < 1e-12 throughout
            Icont = Icont + integral(integrand, a, b, 'AbsTol', 1e-30, 'RelTol', 1e-9);
        end
        Gamma = (Icont + Inwa) / hbar;
        r = Gamma / kplusinvlimit;

    % =====================================================================
    case 'AxA'   % axial-vector A': pi0/eta/eta' A'-mixing, K -> pi A'* -> pi chi chi
    % =====================================================================
    % Paper eq. (kpichichiA):
    %   * lambda_jpiK = i x (real).  The common i drops out of |amp|^2 but the
    %     RELATIVE signs of the pi0/eta/eta' terms do not, so the real
    %     coefficients enter signed.
    %   * the interference term carries the A' propagator 1/((s-mA^2)^2 +
    %     mA^2 GA^2), and its meson denominators use m_pi+;
    %   * the |A'|^2 term's propagator is s - mA^2 = mK^2+mpi^2 - 2 Epi mK, with
    %     the width inside the F^4 mA^4 factor;
    %   * lambda_A'K^2 mA^4 mK^4;
    %   * Gamma_A' = the relic code's total width.
    % The E_chi dependence is quadratic, so the inner integral is analytic.  For
    % 2 mx < mA < mK - mpi+ the A' is on shell inside the range; that pole is
    % razor-thin (Gamma_A'/mA ~ gx^2/12pi) and NWA-notched like SxP's mediator.
        mmed = mphi;  mdm = mx;
        Gammapi       = 7.8 * eV;
        Gammaeta      = 1.31 * 1000 * eV;
        Gammaetaprime = 0.196 * MeV;
        Wmed     = phi_partial_widths(mphi, mx, gx, g, 'axialVector');
        Gammamed = Wmed.total;

        % A'-meson mixing amplitudes (Epi-independent)
        mpi0Aprime  = F*((Cu - Cd)*(1 + (16*B0*L4)/F^2*(mu + md + ms)) + (16*B0*L5*(Cu*mu - Cd*md))/F^2);
        metaAprime  = (sqrt(3)*F)/3*((Cu + Cd - 2*Cs)*(1 + (16*B0*L4)/F^2*(mqu + mqd + mqs)) ...
                        + 1/F^2 * 16*B0*L5*(Cu*mqu + Cd*mqd - 2*Cs*mqs));
        meta0Aprime = (sqrt(6)*F)/3*((Cu + Cd + Cs)*(1 + (16*B0*L4)/F^2*(mqu + mqd + mqs)) ...
                        + 1/F^2 * 16*B0*L5*(Cu*mqu + Cd*mqd + Cs*mqs));
        ljchi = lambdaa*[mpi0Aprime, metaAprime, meta0Aprime]/mmed^2;
        lak = g8*(Cs + Cu) + g8s*(Cu + Cd + Cs);
        lap = g8*(Cd + Cu) + g8s*(Cu + Cd + Cs);

        % lambda_jpiK / i (real, signed)
        pikbulk = @(Epi) g8*((Cd + Cu)*(Epi*mkplus - mpiplus^2) - (Cs + Cu)*(Epi*mkplus - mkplus^2)) ...
                         + g8s*(mkplus^2 - mpiplus^2)*(Cu + Cd + Cs);
        lpk = { @(Epi) g8/F^3*(mkplus^2 - 2*Epi*mkplus) + 2*mpi0Aprime/(F^2*mmed^2)*pikbulk(Epi), ...
                @(Epi) 2*sqrt(3)*g8/(3*F^3)*(mkplus^2 - mpiplus^2) + 2*metaAprime/(F^2*mmed^2)*pikbulk(Epi), ...
                @(Epi) sqrt(6)/(3*F^3)*(2*g8 + 3*g8s)*(mkplus^2 - mpiplus^2) + 2*meta0Aprime/(F^2*mmed^2)*pikbulk(Epi) };
        mjv = [mpi0, meta, meta0];  Gjv = [Gammapi, Gammaeta, Gammaetaprime];

        s_of = @(Epi) mkplus^2 + mpiplus^2 - 2*Epi*mkplus;             % s_chichi
        DA   = @(Epi) (s_of(Epi) - mmed^2).^2 + mmed^2*Gammamed^2;     % |A' propagator|^-2
        pref = (GF^2*F^8)/(128*pi^3*mkplus)*Vus^2*Vstarud^2 / hbar;
        c3   = -8*lambdaa^2/(F^4*mmed^4);
        lapk = lak + lap;  lamk = lak - lap;

        % E_chi moments over [c - d, c + d]: int 1, Echi, Echi^2
        dch = @(Epi) 0.5*sqrt(max(0, (Epi.^2 - mpiplus^2).*(s_of(Epi) - 4*mdm^2)./s_of(Epi)));
        I0  = @(Epi) 2*dch(Epi);
        I1  = @(Epi) (mkplus - Epi).*dch(Epi);
        I2  = @(Epi) 2*dch(Epi).*(3*(0.5*(mkplus - Epi)).^2 + dch(Epi).^2)/3;

        % |A'|^2 polynomial, Echi-independent part (Q0) and the Echi, Echi^2 coefficients
        Q2 = 4*mmed^4*mkplus^2*lapk^2;
        Q1 = @(Epi) 4*mmed^4*mkplus^2*(Epi - mkplus)*lapk^2;
        Q0 = @(Epi) 8*Epi.^3*mkplus^3*mdm^2*lamk^2 ...
            + 4*Epi.^2*mkplus^2*mdm^2*lamk*(2*mmed^2*lamk + mkplus^2*(lap - 5*lak) - mpiplus^2*(lak - 5*lap)) ...
            - 2*Epi*mkplus*( mmed^4*(mkplus^2*lapk^2 + 4*lak*lap*mdm^2) ...
                + 8*mmed^2*mdm^2*lamk*(lak*mkplus^2 - lap*mpiplus^2) ...
                - 4*mdm^2*(lak*mkplus^4*(2*lak - lap) + mkplus^2*mpiplus^2*(lak^2 - 4*lak*lap + lap^2) + lap*mpiplus^4*(2*lap - lak)) ) ...
            + mmed^4*(mkplus^4*lapk^2 + mkplus^2*(mpiplus^2*lapk^2 - 4*lak^2*mdm^2) - 4*lap^2*mpiplus^2*mdm^2) ...
            + 8*mmed^2*mdm^2*(lak*mkplus^2 - lap*mpiplus^2)^2 ...
            - 4*mdm^2*(mkplus^2 + mpiplus^2)*(lak*mkplus^2 - lap*mpiplus^2)^2;
        N3 = @(Epi) pref*c3*(Q0(Epi).*I0(Epi) + Q1(Epi).*I1(Epi) + Q2*I2(Epi));   % times 1/DA

        % meson sums: S1 = |sum_j num_j/(s - mj^2 + i mj Gj)|, S2 = the A'-meson
        % interference sum, num_j = lambda_jchi * lambda_jpiK/i
        amp1 = @(Epi, jj) ljchi(jj)*lpk{jj}(Epi)./(s_of(Epi) - mjv(jj)^2 + 1i*mjv(jj)*Gjv(jj));
        int2 = @(Epi, jj) ljchi(jj)*lpk{jj}(Epi).*((mmed^2 - s_of(Epi)).*(mjv(jj)^2 - s_of(Epi)) ...
                              + Gammamed*Gjv(jj)*mmed*mjv(jj)) ./ ((mjv(jj)^2 - s_of(Epi)).^2 + Gjv(jj)^2*mjv(jj)^2);
        S1 = @(Epi) abs(amp1(Epi, 1) + amp1(Epi, 2) + amp1(Epi, 3));
        S2 = @(Epi) int2(Epi, 1) + int2(Epi, 2) + int2(Epi, 3);
        T12 = @(Epi) 8*mdm^2*s_of(Epi).*S1(Epi).^2 ...
            - 32*lambdaa*mdm^2*(mmed^2 - s_of(Epi)).*(mkplus*Epi*lamk - mkplus^2*lak + mpiplus^2*lap).*S2(Epi) ...
              ./(F^2*mmed^2*DA(Epi));
        outer = @(Epi) (pref*T12(Epi).*I0(Epi) + N3(Epi)./DA(Epi)) .* Aacc(s_of(Epi));

        emax = Epimax(mx);
        % pi0 pole (razor-thin, inside the range for mx < ~0.0675): adaptive
        % waypoints, as before.
        Epipi = (mkplus^2 + mpiplus^2 - mpi0^2)/(2*mkplus);
        xw = linspace(-6, -1, 50);
        wp = [Epipi*(1 + 10.^xw), Epipi*(1 - 10.^xw)];
        wp = wp(wp > mpiplus & wp < emax);
        % on-shell A' pole: NWA notch (analytic Lorentzian over the window, the
        % rest numerically), graded breakpoints outward as in the SxP branch.
        notch = zeros(0, 2);  Inwa = 0;  grad = [];
        wA = mmed*Gammamed;
        EpiA = (mkplus^2 + mpiplus^2 - mmed^2)/(2*mkplus);
        if wA > 0 && EpiA > mpiplus && EpiA < emax
            rng  = emax - mpiplus;
            dEpi = 1e3*wA/(2*mkplus);
            if dEpi <= 0.05*rng
                dEpi = max([dEpi, 8*eps(EpiA), 1e-7*rng]);
                d_lo = min(dEpi, EpiA - mpiplus);  d_hi = min(dEpi, emax - EpiA);
                notch = [EpiA - d_lo, EpiA + d_hi];
                Inwa = N3(EpiA)*Aacc(mmed^2)*(atan(2*mkplus*d_lo/wA) - atan(-2*mkplus*d_hi/wA))/(2*mkplus*wA);
                h = 0.5*(d_lo + d_hi);
                f = 2.^(1:min(80, ceil(log2(rng/max(h, realmin)))));
                grad = [EpiA + h*f, EpiA - h*f];
                grad = grad(grad > mpiplus & grad < emax);
            end
            wp = [wp, EpiA];
        end
        bnds = unique(sort([mpiplus, emax, bpE(:)', wp, notch, grad]));
        Icont = 0;
        for ii = 1:numel(bnds)-1
            a = bnds(ii);  b = bnds(ii+1);
            if b <= a, continue; end
            if ~isempty(notch) && 0.5*(a + b) > notch(1) && 0.5*(a + b) < notch(2), continue; end
            if na62_dead(a, b, mkplus, mpiplus), continue; end   % A < 1e-12 throughout
            Icont = Icont + integral(outer, a, b, 'AbsTol', 0, 'RelTol', 1e-8);
        end
        Gamma = Icont + Inwa;
        r = Gamma / kplusinvlimit;

    otherwise
        r = NaN;
    end
end

% NA62 K+ -> pi+ nu nu acceptance in m_miss^2 = s [GeV^2]: R1 = [0, 0.010] and
% R2 = [0.026, 0.068], each edge smeared by the resolution sigma = 0.001 (JHEP 06
% (2021) 093).  A -> 1 deep inside a region, 1/2 on an edge, ~0 in the gaps.
% Shared with mesonandCMB_ratio.m's kphi (as a copy there).
function A = na62_acceptance(s)
    sig = 1e-3;  Phi = @(x) 0.5*erfc(-x/sqrt(2));
    A = Phi((0.010 - s)/sig) - Phi((0 - s)/sig) + Phi((0.068 - s)/sig) - Phi((0.026 - s)/sig);
end

% Breakpoints in s for integrating A(s): each edge and +-2, +-4 sigma around it.
function e = na62_edges()
    e = [0.010, 0.026, 0.068]' + 1e-3*[-4 -2 0 2 4];
    e = e(:)';
end

% Gamma(phi -> meson pairs) of the scalar phi-SM vertex (SxS/PxS) [GeV], from
% phi_partial_widths.m / meson_pairs.m (the relic code's source of truth).
function Gm = meson_pair_width(mphi, g)
    W  = phi_partial_widths(mphi, 0, 0, g, 'heavyScalar');
    Gm = W.pipi + W.KK + W.eta;
end

% int_a^b numf(Epi) / ((s - mphi^2)^2 + mphi^2 Gphi^2) dEpi, s = mK^2 + mpi^2 - 2 Epi mK.
% An in-range pole narrower than 5% of the range is NWA-notched: over a window
% [Epi* - d, Epi* + d] the Lorentzian is integrated analytically with numf frozen
% at the pole, (atan(u_lo/w) - atan(u_hi/w))/(2 mK w), w = mphi Gphi, and the rest
% numerically with breakpoints graded outward from the window -- the scheme of
% the SxP/PxP mediator notch (Gamma_phi/mphi can be ~1e-12, too narrow for the
% adaptive rule).
function F = bw_notch_integral(numf, mphi, Gphi, a, b, mK, mpi, extra)
    % extra: additional breakpoints (here the NA62 acceptance edges)
    extra = extra(extra > a & extra < b);
    den = @(Epi) (mK^2 + mpi^2 - 2*mK*Epi - mphi^2).^2 + mphi^2*Gphi^2;
    f   = @(Epi) numf(Epi) ./ den(Epi);
    w   = mphi*Gphi;
    Ep  = (mK^2 + mpi^2 - mphi^2)/(2*mK);
    if ~(w > 0 && Ep > a && Ep < b)
        F = panels(f, unique(sort([a, b, extra])), mK, mpi);
        return;
    end
    rng = b - a;
    d = 1e3*w/(2*mK);
    if d > 0.05*rng          % broad: the adaptive rule resolves it given a breakpoint
        F = panels(f, unique(sort([a, b, Ep, extra])), mK, mpi);
        return;
    end
    d = max([d, 8*eps(Ep), 1e-7*rng]);
    d_lo = min(d, Ep - a);  d_hi = min(d, b - Ep);
    Inwa = numf(Ep) * (atan(2*mK*d_lo/w) - atan(-2*mK*d_hi/w)) / (2*mK*w);
    h = 0.5*(d_lo + d_hi);
    fac = 2.^(1:min(80, ceil(log2(rng/max(h, realmin)))));
    bnds = unique(sort([a, b, Ep - d_lo, Ep + d_hi, Ep + h*fac, Ep - h*fac, extra]));
    bnds = bnds(bnds >= a & bnds <= b);
    F = Inwa;
    for ii = 1:numel(bnds)-1
        lo = bnds(ii);  hi = bnds(ii+1);
        if hi <= lo || (0.5*(lo + hi) > Ep - d_lo && 0.5*(lo + hi) < Ep + d_hi), continue; end
        if na62_dead(lo, hi, mK, mpi), continue; end
        F = F + integral(f, lo, hi, 'AbsTol', 0, 'RelTol', 1e-9);
    end
end

function F = panels(f, bnds, mK, mpi)
    F = 0;
    for ii = 1:numel(bnds)-1
        if na62_dead(bnds(ii), bnds(ii+1), mK, mpi), continue; end
        F = F + integral(f, bnds(ii), bnds(ii+1), 'AbsTol', 0, 'RelTol', 1e-9);
    end
end

% True when the E_pi panel [Ea, Eb] lies wholly > 7 sigma inside the NA62 gap
% (0.010, 0.026) or above 0.068 + 7 sigma: A(s) < 1e-12 there, so skip the
% integral (the pi0 waypoints and notch octaves all sit in the gap).
function d = na62_dead(Ea, Eb, mK, mpi)
    sLo = mK^2 + mpi^2 - 2*mK*max(Ea, Eb);  sHi = mK^2 + mpi^2 - 2*mK*min(Ea, Eb);
    d = (sLo > 0.017 && sHi < 0.019) || sLo > 0.075;
end
