function [dNeff, info] = neff_simple(mphi, mx, gx, g, interaction, opts)
% neff_simple  CMB Delta-Neff from dark species that are in equilibrium with
% the EM plasma at neutrino decoupling and later hand their entropy to photons.
%
%   [dNeff, info] = neff_simple(mphi, mx, gx, g, interaction [, opts])
%
% Boehm, Dolan & McCabe (2013):
%
%   F(y)  = 30/(7 pi^4) int_y^inf dxi (4 xi^2 - y^2) sqrt(xi^2 - y^2) / (e^xi +/- 1)
%   N_eff = Neff_SM [1 + (7/22) sum_i (g_i/2) F_i(y_i)]^(-4/3),  y_i = m_i/T_nu_dec
%
% (+ fermion, - boson).  F is the species' entropy in units of one massless
% Weyl-pair fermion, so F_fermion(0) = 1 and F_boson(0) = 8/7; the boson needs
% no separate prefactor.  (7/22)(g/2) = (7/8) g / (11/2): entropy relative to
% the gamma + e+- plasma it ends up heating.
%
% Each species enters with a weight for its contact with the plasma around
% T_nu_dec: equilibrated by T_nu_dec, times still in contact after it,
%   w = [1 - exp(-int_{T_nu_dec}^{T_RH} (Gamma/H) dT/T)]
%     * [1 - exp(-int_{T_nu_dec/100}^{T_nu_dec} (Gamma/H) dT/T)].
% phi's rate (Gamma/H grows as T falls) saturates the second factor; chi's
% (Gamma/H falls as T falls) saturates the first.
%   phi  inverse decay gamma gamma <-> phi, Gamma_gg K1/K2, plus pair
%        production gamma gamma <-> phi phi (photon exchange, two phi F F
%        vertices), per phi 2 R / n_phi^eq.  phi F F~
%        gives the same sigma in terms of Gamma_gg, so SP/PP use it unchanged.
%        S⊗S / P⊗S with mphi < 2 m_pi add pi+- gamma <-> pi+- phi
%        (pigamma_phi_rate; tree level, independent of mphi, dominant for light
%        phi; docs/phi_thermalisation_check.md).
%        All non-axial models add e+- gamma <-> e+- phi (Primakoff, Debye-
%        screened) and e+ e- <-> gamma phi through the same phi gamma gamma
%        vertex (primakoff_phi_rate; independent of mphi for mphi << T).
%
% A phi that equilibrates before T_nu_dec but loses contact by then (rates that
% die as T falls: pi gamma -> pi phi, gamma gamma -> phi phi) has decoupled
% alongside the neutrinos: it does not heat the photons, it is free radiation
% at T_phi = T_nu.  Its share w_free = weq (1 - wct) enters the free-streaming
% phi budget below, at Q = w_free rho_phi^eq(T_nu_dec) / T_nu_dec^4 (massless,
% fully equilibrated: Delta Neff = 4/7).
% The phi still in contact after T_nu_dec (w_phi) stays in the photon bath
% until T_late, where its optical depth counted up from recombination,
% int_{T_rec}^{T_late} (Gamma/H) dT/T, reaches 1 (T_late = T_nu_dec if it
% never does).  A link that dies with the e+- (Primakoff) releases a still-
% relativistic phi there as free radiation; one that stays fast (inverse decay
% of a massive phi) puts T_late where phi has already annihilated, and nothing
% is released.
%   chi  gamma gamma <-> chi chi (ggxxrate/n_gamma), or through phi,
%        min((n_chi^eq/2) <sv>_{chi chi -> phi phi}, Gamma_inv)   (n_chi^eq = chi + chibar).
%
% Caveats (not handled):
%   * chi, and the phi it makes, are taken at T_gamma; a dark bath at its own
%     temperature enters only through the chain bottleneck.
%   * The entropy phi carries off assumes a thermal-shape phi.
%   * phi <-> chi chi inverse decay is not among chi's links.
%
% opts (all optional):
%   T_nu_dec  neutrino decoupling T [GeV]              (2e-3)
%   Neff_SM   N_eff with no dark sector                 (3.046)
%   species   'phi' | 'phi+chi'                         ('phi+chi')
%   T_BBN, T_RH  range for the "ever coupled" check     (1e-3, 1e-2)
%   pigamma   include pi gamma <-> pi phi in phi's link  (true)
%   primakoff include e gamma <-> e phi, ee <-> gamma phi (true)
%
% info: phi_coupled, w_phi, chi_coupled, w_chi, phi_ever_coupled, Gamma_inv_over_H,
%       Gamma_pp_over_H, Gamma_pg_over_H, Gamma_prim_over_H, w_phi_free,
%       Gamma22_over_H, Gamma_chiphi_over_H, y_phi, y_chi, F_phi, F_chi, X, Neff;
%       X_phi, X_chi (the two terms of X); link_phi = integrated Gamma/H over
%       T_nu_dec/100..T_RH of [inverse decay, gg<->phiphi, pi gamma, Primakoff],
%       link_chi = same for [gg<->chichi, chichi->phiphi capped by phi's link]
%       (diagnostics, for relic_channels_g / neff_channels_g).

if nargin < 6, opts = struct(); end
T_nu    = getopt(opts, 'T_nu_dec', 2e-3);
Neff_SM = getopt(opts, 'Neff_SM', 3.046);
species = getopt(opts, 'species', 'phi+chi');
T_BBN   = getopt(opts, 'T_BBN', 1e-3);
T_RH    = getopt(opts, 'T_RH', 1e-2);
use_pg  = getopt(opts, 'pigamma', true);
use_pr  = getopt(opts, 'primakoff', true);
if ~any(strcmp(species, {'phi', 'phi+chi'}))
    error('neff_simple:species', 'species must be ''phi'' or ''phi+chi''');
end

Mpl   = 1.22e19;
zeta3 = 1.2020569;
g_chi = 4;                                          % Dirac fermion
if strcmp(interaction, 'axialVector'), g_phi = 3; else, g_phi = 1; end

Hof  = @(T) 1.66 .* sqrt(gstarofT_vec(T, 0)) .* T.^2 ./ Mpl;
ngam = @(T) (2*zeta3/pi^2) .* T.^3;

% Rates on one grid, T_nu_dec/100 .. T_RH, split at T_nu_dec.
Tb = logspace(log10(T_nu/100), log10(T_nu), 30);     % below T_nu_dec
Ta = logspace(log10(T_nu), log10(T_RH), 30);         % above
Tall = [Tb, Ta(2:end)];  Hall = Hof(Tall);  ib = 1:numel(Tb);  ia = numel(Tb):numel(Tall);
weq = @(G) 1 - exp(-trapz(log(Tall(ia)), G(ia)./Hall(ia)));   % equilibrated by T_nu_dec
wct = @(G) 1 - exp(-trapz(log(Tall(ib)), G(ib)./Hall(ib)));   % contact after it
wt  = @(G) weq(G) * wct(G);

% --- phi: inverse decay gamma gamma <-> phi, and gamma gamma <-> phi phi ---
Gw   = phi_partial_widths(mphi, mx, gx, g, interaction);
Ginv = @(T) Gw.gg .* k1k2(mphi ./ T);
Gpp  = @(T) ggphiphi_rate(T, mphi, mx, gx, g, interaction, Gw, g_phi);
has_pg = use_pg && mphi < 2*0.13957 && ...
         ismember(interaction, {'heavyScalar','lightScalar','heavyScalarPS','lightScalarPS'});
if has_pg, Gpgf = @(T) pigamma_phi_rate(T, mphi, g); else, Gpgf = @(T) zeros(size(T)); end
g2_pr = 0;                                           % g_phigamgam^2, finite as mphi -> 0
if use_pr && ~strcmp(interaction, 'axialVector')
    mref = max(mphi, 1e-12);
    if mphi < mref, Gwr = phi_partial_widths(mref, mx, gx, g, interaction); else, Gwr = Gw; end
    g2_pr = 64*pi * Gwr.gg / mref^3;
end
Gprf = @(T) primakoff_phi_rate(T, mphi, g2_pr);
Gid  = Ginv(Tall);  Gp = Gpp(Tall);  Gpg = Gpgf(Tall);  Gpr = Gprf(Tall);
Gi   = Gid + Gp + Gpg + Gpr;                         % phi's whole SM link
r_inv = Gid(numel(Tb)) / Hall(numel(Tb));
r_pp  = Gp(numel(Tb))  / Hall(numel(Tb));
r_pr  = Gpr(numel(Tb)) / Hall(numel(Tb));
r_pg  = Gpg(numel(Tb)) / Hall(numel(Tb));
phi_coupled = Gi(numel(Tb)) >= Hall(numel(Tb));
w_phi  = wt(Gi);
w_free = weq(Gi) * (1 - wct(Gi));                    % equilibrated, then decoupled before T_nu_dec

% Late decoupling of the bath-coupled phi: T_late from the optical depth below
% T_nu_dec, counted up from recombination (continuous in the coupling).
T_rec = 2.6e-10;
Tlo  = logspace(log10(T_rec), log10(T_nu), 60);
Glo  = Ginv(Tlo) + Gpp(Tlo) + Gpgf(Tlo) + Gprf(Tlo);
% Segments use the logarithmic mean of Gamma/H (exact for a power law): the
% trapezoid overshoots where the rate dies as e^{-m_e/T}, putting T_late 3-10% low.
rlo  = Glo ./ Hof(Tlo);  r1 = rlo(1:end-1);  r2 = rlo(2:end);
rm   = 0.5*(r1 + r2);
lg   = r1 > 0 & r2 > 0 & abs(r2 - r1) > 1e-6*rm;
rm(lg) = (r2(lg) - r1(lg)) ./ log(r2(lg)./r1(lg));
rm(xor(r1 > 0, r2 > 0)) = 0;
tau  = [0, cumsum(rm .* diff(log(Tlo)))];
k1   = find(tau >= 1, 1);
if isempty(k1)
    T_late = T_nu;
elseif k1 == 1
    T_late = T_rec;
else
    lt     = log(Tlo(k1-1)) + (1 - tau(k1-1))/(tau(k1) - tau(k1-1)) * log(Tlo(k1)/Tlo(k1-1));
    T_late = exp(lt);
end
% Released only if phi still carries energy there (otherwise its entropy has
% already gone to the photons and the heating formula has it).
late_rel = w_phi > 1e-10 && rho_eq(mphi, T_late, g_phi) > 1e-6 * rho_eq(0, T_late, g_phi);

% gamma gamma -> chi chi on both grids in ONE ggxxrate call.  Each T's value is
% independent of the others in the call (both grids end at T_RH, which sets the
% fast path's s_max; the multi-pole path is per-T), so this is bit-identical to
% two calls, and the fast path evaluates sigma(s) on its s-nodes once.
Tg = logspace(log10(T_BBN), log10(T_RH), 40);
if strcmp(species, 'phi+chi'), Tq = [Tg, Tall]; else, Tq = Tg; end
Rq = max(0, ggxxrate(Tq, mx, mphi, gx, g, interaction)) ./ ngam(Tq);

% "ever coupled" over [T_BBN, T_RH]: portal max(Gamma_22, Gamma_inv) vs H.
G22 = Rq(1:numel(Tg));
phi_ever_coupled = max(max(G22, Ginv(Tg) + Gpp(Tg) + Gpgf(Tg) + Gprf(Tg)) ./ Hof(Tg)) >= 1;

% Diagnostic only: each link's integrated Gamma/H over T_nu_dec/100 .. T_RH.
ilnT = @(G) trapz(log(Tall), G ./ Hall);
link_phi = [ilnT(Gid), ilnT(Gp), ilnT(Gpg), ilnT(Gpr)];   % inverse decay, gg<->phiphi, pi gamma, Primakoff
link_chi = [0, 0];                                         % gg<->chichi, chichi->phiphi (capped by phi's link)

% --- chi ---
r_22 = NaN;  r_cp = NaN;  w_chi = 0;  w_chi_eq = 0;  chi_coupled = false;
if strcmp(species, 'phi+chi')
    Gc   = Rq(numel(Tg)+1:end);
    r_22 = Gc(numel(Tb)) / Hall(numel(Tb));
    link_chi(1) = ilnT(Gc);
    if any(Gi > Gc)                      % the link through phi can matter
        % per-chi rate: a chi annihilates on chibar only, n_chibar = n_eq/2 (g_chi = 4 counts both)
        Gcp  = arrayfun(@(T) 0.5 * n_eq(mx, T, g_chi, +1) * ...
                        svavg_phiphi(mx, mphi, gx, g, T, interaction), Tall);
        Gcp(~isfinite(Gcp)) = 0;         % min() would otherwise drop a NaN
        r_cp = Gcp(numel(Tb)) / Hall(numel(Tb));
        link_chi(2) = ilnT(min(Gcp, Gi));
        Gc   = max(Gc, min(Gcp, Gi));
    end
    w_chi       = wt(Gc);
    w_chi_eq    = weq(Gc);
    chi_coupled = Gc(numel(Tb)) >= Hall(numel(Tb));
end

% --- heating formula (Boehm-Dolan-McCabe) ---
y_phi = mphi / T_nu;   F_phi = Fy(y_phi, -1);
y_chi = mx   / T_nu;   F_chi = Fy(y_chi, +1);
X = (7/22) * ( (g_phi/2)*F_phi*w_phi + (g_chi/2)*F_chi*w_chi );
dNeff_heat = Neff_SM * ((1 + X)^(-4/3) - 1);

% --- phi made by chi: free-streaming dark radiation --------------------------
% chi (at T_gamma) feeds phi through chi chibar -> phi phi and chi chibar -> phi.
% The phi energy relaxes toward equilibrium at T_gamma,
%     d(rho_phi a^4)/dt = a^4 k H rho_phi^eq (1 - rho_phi/rho_phi^eq),
%     k = min( C_E/(H rho_phi^eq),  Gamma_chi/H ),
%     C_E = (n_chi^2/4) <sv E>_{chi chi->phi phi} + Gamma_xx m_phi n_phi^eq   (n_chi = chi + chibar),
% which covers both freeze-in and thermalisation.  k is the chain bottleneck:
% the energy reaches phi from the photons through chi, so it can go no faster
% than chi's own photon link Gamma_chi.  Only phi NOT already counted
% through the photon link: the target is (1 - w_phi) rho_phi^eq.
% The phi that decoupled from the photons before T_nu_dec (w_free) is added to
% Q where the grid crosses T_nu_dec, taking the larger of it and what chi made.
% Comoving variable Q = rho_phi / T_nu^4 (a = 1/T_nu; T_nu = T_gamma above
% T_nu_dec).  Entropy the phi carries off is not photon heating, so one entropy
% balance from T_nu_dec to recombination replaces (11/4)(1+X):
%     S_i = (11/2)(1+X) + s_phi(T_nu_dec),
%     gs_bath(T_gamma) (T_gamma/T_nu)^3 + s_phi(Q) = S_i,
% in units of T_nu^3, with s_phi from a thermal-shape phi.  Q = 0 gives back
% the heating formula exactly.
Q = 0;  P_free = 1;  dNeff_phichi = 0;  Tphi_Tnu = 0;  TgTnu3 = (11/4)*(1 + X);
chi_chain = strcmp(species, 'phi+chi') && w_chi_eq*(1 - w_phi) > 1e-10;
if chi_chain || w_free > 1e-10 || late_rel
    T_end = min(T_nu/100, mx/50);
    if late_rel, T_end = min(T_end, T_late/1.5); end
    Tc    = logspace(log10(T_RH), log10(T_end), 60);            % descending
    if late_rel && T_late < T_nu                                % release exactly at T_late,
        Tc = sort([Tc(abs(log(Tc/T_late)) > 1e-9), T_late], 'descend');   % not a grid step below
    end
    rq    = arrayfun(@(T) rho_eq(mphi, T, g_phi), Tc);
    if chi_chain
        svE   = arrayfun(@(T) svavg_phiphi(mx, mphi, gx, g, T, interaction, 'energy'), Tc);
        svE(~isfinite(svE)) = 0;
        nx    = arrayfun(@(T) n_eq(mx,   T, g_chi, +1), Tc);
        np    = arrayfun(@(T) n_eq(mphi, T, g_phi, -1), Tc);
        % reaction density n_chi n_chibar <sv E> = nx^2 <sv E>/4: nx counts chi + chibar
        CE    = nx.^2/4 .* max(svE, 0) + Gw.xx * mphi * np;
        kc    = CE ./ Hof(Tc) ./ rq;                               % relaxation rate per ln T
        kc(~(rq > 0)) = 0;                                        % no phi reservoir left
        Gcl   = exp(interp1(log(Tall), log(max(Gc, realmin)), log(Tc), 'linear', 'extrap'));
        kc    = min(kc, Gcl ./ Hof(Tc));                          % chain bottleneck
    else
        kc    = zeros(size(Tc));
    end
    Q_free = w_free * rho_eq(mphi, T_nu, g_phi) / T_nu^4;
    % gs of the photon-contact bath, sampled once and splined
    gsb0  = arrayfun(@(T) 2 + 3.5*Fy(0.511e-3/T, +1) + w_chi*3.5*Fy(mx/T, +1), Tc);
    gsp   = arrayfun(@(T) (7/8)*g_phi*Fy(mphi/T, -1), Tc);      % phi's share, x w_bath
    sphi  = @(Q) g_phi * (30*Q/(pi^2*g_phi))^(3/4);
    S_i   = NaN;  Tnu_T = 1;  w_bath = w_phi;
    for i = 2:numel(Tc)
        T = Tc(i);
        if T < T_nu && isnan(S_i)
            Q   = max(Q, Q_free);
            S_i = (11/2)*(1 + X) + sphi(Q);
        end
        if ~isnan(S_i)
            Tnu_T = ((gsb0(i) + w_bath*gsp(i)) / max(S_i - sphi(Q), 1e-3))^(1/3);   % T_nu / T_gamma
        end
        if late_rel && w_bath > 0 && T <= T_late && ~isnan(S_i)
            % bath phi leaves as free radiation.  The bath gives up exactly the
            % entropy it carried for phi (w_bath gsp), not sphi(Q)'s thermal-shape
            % estimate of it, so S_i is rebased to keep S_i - sphi(Q) = bath.
            Qn     = Q + w_bath * rq(i) / (Tnu_T*T)^4;
            S_i    = S_i - sphi(Q) + sphi(Qn) - w_bath * gsp(i) / Tnu_T^3;
            Q      = Qn;
            w_bath = 0;
        end
        Qeq = (1 - w_bath) * rq(i) / (Tnu_T*T)^4;
        dlt = log(Tc(i-1)/T);
        Q   = Qeq + (Q - Qeq) * exp(-0.5*(kc(i-1) + kc(i)) * dlt);
    end
    if isnan(S_i), Q = max(Q, Q_free); S_i = (11/2)*(1 + X) + sphi(Q); end
    % phi that decays (or is re-absorbed) into photons before recombination
    % returns its entropy: drop it from the free-streaming budget.
    Tl     = logspace(log10(T_rec), log10(T_nu), 60);
    Gsm    = Gw.gg + Gw.pipi + Gw.KK + Gw.eta;
    P_free = exp(-trapz(log(Tl), Gsm .* k1k2(mphi ./ Tl) ./ Hof(Tl)));
    Qf     = P_free * Q;
    TgTnu3 = (S_i - sphi(Qf)) / 2;                                % (T_gamma/T_nu)^3 at rec
    Tphi_Tnu = (30*Qf/(pi^2*g_phi))^(1/4);
    y_rec  = mphi / (Tphi_Tnu * T_rec / TgTnu3^(1/3));            % m_phi / T_phi at rec
    rho1nu_ref = (7/8)*(pi^2/15)*(4/11)^(4/3);                    % per T_gamma^4
    if Qf > 0
        dNeff_phichi = Qf * rad_frac(y_rec) / TgTnu3^(4/3) / rho1nu_ref;
    end
end
Neff  = Neff_SM * ((4/11) * TgTnu3)^(-4/3) + dNeff_phichi;
dNeff = Neff - Neff_SM;

info = struct('phi_coupled', phi_coupled, 'w_phi', w_phi, 'chi_coupled', chi_coupled, 'w_chi', w_chi, 'w_chi_eq', w_chi_eq, ...
    'phi_ever_coupled', phi_ever_coupled, 'Gamma_inv_over_H', r_inv, 'Gamma_pp_over_H', r_pp, 'Gamma_pg_over_H', r_pg, 'Gamma_prim_over_H', r_pr, 'w_phi_free', w_free, ...
    'Gamma22_over_H', r_22, 'Gamma_chiphi_over_H', r_cp, ...
    'y_phi', y_phi, 'y_chi', y_chi, 'F_phi', F_phi, 'F_chi', F_chi, ...
    'X', X, 'dNeff_heat', dNeff_heat, 'dNeff_phi_rad', dNeff_phichi, ...
    'Tphi_over_Tnu', Tphi_Tnu, 'P_free', P_free, 'T_late', T_late, 'late_release', late_rel, 'Neff', Neff, ...
    'X_phi', (7/22)*(g_phi/2)*F_phi*w_phi, 'X_chi', (7/22)*(g_chi/2)*F_chi*w_chi, ...
    'link_phi', link_phi, 'link_chi', link_chi);
end

% ======================================================================
function F = Fy(y, sgn)
% Entropy of one species per massless Weyl-pair fermion (sgn +1 fermion, -1 boson).
if y > 100, F = 0; return; end
f = @(u) (4*u.^2 - y^2) .* sqrt(max(0, u.^2 - y^2)) ./ (exp(u) + sgn);
F = 30/(7*pi^4) * integral(f, y, y + 60, 'RelTol', 1e-8, 'AbsTol', 0);
end

function r = k1k2(z)
% K1/K2 (time dilation <1/gamma>); massless limit z/2 -> 0.
r = besselk(1, z, 1) ./ besselk(2, z, 1);
r(~isfinite(r)) = 0;
end

function G = ggphiphi_rate(TT, mphi, mx, gx, g, interaction, Gw, g_phi)
% Per-phi rate of gamma gamma <-> phi phi: 2 R / n_phi^eq, with the reaction
% density R = T/(16 pi^4) int ds sigma s^(3/2) K1(sqrt(s)/T) (ggxxrate's
% normalisation).
%
% SxS / PxS: the charged-meson loop, boxes plus seagull triangle, through the
% low-energy theorem:
%   sigma = alpha^2 Y^2 s beta / (36864 pi^3),  Y = sum_j mu_j^2 / m_j^4,
% mu_pi = B0 (g_u+g_d), mu_K = B0 (g_u+g_s).  Photon exchange is O(alpha)
% smaller (~1e-7 at g = [1 1 1]/3) and is dropped there.  Valid for
% sqrt(s) << 2 m_pi, i.e. T << 50 MeV.
% Other models keep t/u-channel photon exchange, with
% the vertex fixed by Gamma_gg = g_phigg^2 mphi^3 / (64 pi); at mphi = 0 it is
% taken from a tiny mphi (g_phigg stays finite).
G = zeros(size(TT));
if ismember(interaction, {'heavyScalar','lightScalar','heavyScalarPS','lightScalarPS'})
    alphaEM = 1/137;  B0 = 0.13957^2/(0.0022 + 0.00469);   % as phi_partial_widths
    Y = B0^2 * ((g(1)+g(2))^2/0.13957^4 + (g(1)+g(3))^2/0.493677^4);
    sig = @(s_) sig_ggphiphi_box(s_, mphi, alphaEM, Y);
    if ~(Y > 0), return; end
else
    mref = max(mphi, 1e-12);
    if mphi < mref, Gw = phi_partial_widths(mref, mx, gx, g, interaction); end
    g2 = 64*pi * Gw.gg / mref^3;                      % g_phigg^2
    sig = @(s_) sig_ggphiphi(s_, mphi, g2);
    if ~(g2 > 0), return; end
end
for k = 1:numel(TT)
    T  = TT(k);
    x0 = 2*mphi / T;
    if x0 > 200, continue; end
    % s = (x T)^2:  R = T^6/(8 pi^4) int_{x0} sigma x^4 K1(x) dx
    f = @(x) sig((x*T).^2) .* x.^4 .* besselk(1, x, 1) .* exp(-(x - x0));
    I = integral(f, x0, x0 + 80, 'RelTol', 1e-6, 'AbsTol', 0);
    n = n_eq(mphi, T, g_phi, -1);
    if n > 0
        G(k) = 2 * T^6/(8*pi^4) * exp(-x0) * I / n;
    end
end
end

function s = sig_ggphiphi_box(s_, m, alphaEM, Y)
% sigma(gamma gamma -> phi phi) from the charged-meson loop, low-energy theorem
% (includes the 1/2 for identical phi phi).
b = sqrt(max(0, 1 - 4*m^2 ./ s_));
s = alphaEM^2 * Y^2 / (36864*pi^3) .* s_ .* b;
end

function s = sig_ggphiphi(s_, m, g2)
% sigma(gamma gamma -> phi phi) by photon exchange, r = m^2/s.
r = m^2 ./ s_;
b = sqrt(max(0, 1 - 4*r));
L = log((1 - 2*r + b) ./ (2*r));
L(r == 0) = 0;                                    % r^2 L -> 0
s = g2^2/(1024*pi) .* s_ .* (b.*(5 - 4*r + 6*r.^2) - 8*r.^2.*(1 - 3*r).*(1 - r)./(1 - 2*r).*L);
s(r >= 1/4) = 0;
end

function n = n_eq(m, T, g, sgn)
y = m / T;
if y > 100, n = 0; return; end
f = @(u) u .* sqrt(max(0, u.^2 - y^2)) ./ (exp(u) + sgn);
n = g/(2*pi^2) * T^3 * integral(f, y, y + 60, 'RelTol', 1e-6, 'AbsTol', 0);
end

function r = rho_eq(m, T, g)
% Equilibrium boson energy density.
y = m / T;
if y > 100, r = 0; return; end
f = @(u) u.^2 .* sqrt(max(0, u.^2 - y^2)) ./ (exp(u) - 1);
r = g/(2*pi^2) * T^4 * integral(f, y, y + 60, 'RelTol', 1e-6, 'AbsTol', 0);
end

function R = rad_frac(y)
% Kinetic (radiation) energy of a free-streaming, relativistically decoupled
% boson, f = 1/(e^{p/T}-1), over its massless value.  y = m/T.
if y == 0, R = 1; return; end
if ~isfinite(y), R = 0; return; end          % no phi at all (Q = 0)
f = @(u) u.^4 ./ (sqrt(u.^2 + y^2) + y) ./ expm1(u);    % = (E - m) u^2 f, no cancellation
R = integral(f, 0, 60, 'RelTol', 1e-6, 'AbsTol', 0) / (pi^4/15);
end

function v = getopt(s, name, dflt)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = dflt; end
end
