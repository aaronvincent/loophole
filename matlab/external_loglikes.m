function L = external_loglikes(mx, mphi, gx, g, interaction, opts)
% external_loglikes  External (non-relic) constraint log-likelihoods.
%
%   L = external_loglikes(mx, mphi, gx, g, interaction [, opts])
%
% Single-file home for the external constraints layered on top of the relic
% scan (applied in analysis/postprocess.py for the GAMBIT scans).  Mirrors model_observables.m's arg order
% and fail-soft philosophy: takes the *physical* model point and returns a struct
% of per-constraint log-likelihoods plus a combined total.
%
% Constraints currently implemented (data in ../data, see data/README.md):
%   sn1987a      SN1987A cooling -- scalar (SS/PS, g_N=hp, Hardy-Sokolov 2410.17347)
%                or pseudoscalar (SP/PP, g_N=htildeP) nucleon coupling; closed
%                exclusion contour, excluded INSIDE, extended left in mphi.
%   neutronstar  Neutron-star cooling (Damiano 2506.19906); scalar NEUTRON
%                coupling (SS/PS, g_N=hn); one-sided upper bound (excluded ABOVE).
%   stellar      Stellar cooling -- combined TRGB + horizontal-branch bound
%                (Hardy & Lasenby); scalar coupling, quoted as
%                alpha_phiNN = g^2/(4 pi) vs mphi and assuming equal phi-p /
%                phi-n coupling, so we feed the average g_bar = (hp + hn)/2.
%                SS/PS only; one-sided upper bound (excluded ABOVE), extended
%                flat to the left in mphi below the digitized data.
%   fifthforce   Long-range / 5th-force spin-dependent bound; gx-independent
%                (a nucleon-nucleon force).  Two families, by nucleon coupling:
%                  axial  (axialVector) -- |g_A|^2 = a_N^2 vs range lambda =
%                    hbarc/m_A' (2408.15691); one-sided upper bound on g_A^2.
%                  pseudo (SP/PP) -- the pseudoscalar spin-spin bounds on
%                    g_P^2 (proton-proton) and |g_P g_N| (proton-neutron) vs
%                    mphi; two independent curves, both one-sided upper bounds,
%                    whose penalties are summed under the single 'fifthforce'
%                    entry (per-curve detail in info.fifthforce.pp / .pn).
%   bullet       DM self-interaction (bullet cluster): sigma_self/m_chi < 1 cm^2/g;
%                velocity-dependent, depends on (mx, mphi, gx) NOT g; one-sided
%                upper bound (excluded ABOVE).  The sigma_self formula lives in
%                sigma_self_over_m().
%   kaon         Rare K+ -> pi+ chi chi invisible decay vs the K+->pi+ + invisible
%                limit; depends on (mx, mphi, gx, g) and the interaction; open only
%                for mx < ~0.177 GeV; excluded when ratio > 1 (kaon_invisible_ratio.m).
%   kphi         K+ -> pi+ phi (on-shell mediator) vs the same K+->pi+ + invisible
%                limit; open only for mphi < ~0.354 GeV; gx-INDEPENDENT (mesonandCMB_ratio.m).
%   cmb          DM annihilation <sigma v> vs the Slatyer CMB energy-injection bound;
%                ~gx^2; not applicable to SxS/SxP (velocity-suppressed) (mesonandCMB_ratio.m).
%   integral     DM annihilation <sigma v> vs the INTEGRAL telescope bound; ~gx^2;
%                all models, light DM only (mesonandCMB_ratio.m).
%   beamdump     Terrestrial / beam-dump -- STUB (no data yet).
%
% Convention (matches gambit/cpp/likelihoods.cpp):
%   * each loglike is <= 0, with 0 in the allowed region and a smooth one-sided
%     Gaussian penalty into the excluded region -- so the same call doubles as a
%     soft mask for postprocessing.
%   * a constraint that does not apply to this model (wrong vertex structure,
%     mass out of the data range, disabled) contributes 0.
%   * each constraint is evaluated in its own try/catch; a failure degrades to 0
%     with the message recorded in L.info.<name>.error.
%
% Inputs:
%   mx          DM mass m_chi [GeV]
%   mphi        mediator mass m_phi (or m_A' for axialVector) [GeV]
%   gx          DM-mediator coupling (phi-chi-chi / lambda_a)
%   g           [g_u g_d g_s] quark couplings
%   interaction vertex string ('heavyScalarPS', 'scalarPP', 'axialVector', ...)
%   opts        (optional) settings struct; see default_opts().
%
% Output struct L:
%   .sn1987a .neutronstar .fifthforce .bullet .kaon .kphi .cmb .integral
%   .mesoninv .beamdump
%               per-constraint loglike
%   .total      sum over finite contributions
%   .info       per-constraint struct (derived observable, bound, applicable/
%               excluded flags, plus .error if it threw)
%
% Example:
%   L = external_loglikes(0.05, 1e-3, 1, [0.6 0.3 0.1], 'heavyScalarSP');
%   L.total, L.sn1987a, L.info.sn1987a

if nargin < 6 || isempty(opts), opts = struct(); end
opts = fill_defaults(opts, default_opts());

pt = struct('mx', mx, 'mphi', mphi, 'gx', gx, 'g', g, 'interaction', interaction);

L = struct();
L.info = struct();

% Each constraint: (field, loglike fn, enable flag).  Add new constraints here.
constraints = { ...
    'sn1987a',     @sn1987a_loglike,     opts.enable_sn1987a;     ...
    'neutronstar', @neutronstar_loglike, opts.enable_neutronstar; ...
    'stellar',     @stellar_loglike,     opts.enable_stellar;     ...
    'fifthforce',  @fifthforce_loglike,  opts.enable_fifthforce;  ...
    'bullet',      @bullet_loglike,      opts.enable_bullet;      ...
    'kaon',        @kaon_loglike,        opts.enable_kaon;        ...
    'kphi',        @kphi_loglike,        opts.enable_kphi;        ...
    'cmb',         @cmb_loglike,         opts.enable_cmb;         ...
    'integral',    @integral_loglike,    opts.enable_integral;    ...
    'mesoninv',    @mesoninv_loglike,    opts.enable_mesoninv;    ...
    'beamdump',    @beamdump_loglike,    opts.enable_beamdump};

for k = 1:size(constraints, 1)
    name = constraints{k, 1};
    fn   = constraints{k, 2};
    on   = constraints{k, 3};
    if ~on
        L.(name) = 0;
        L.info.(name) = struct('applicable', false, 'enabled', false);
        continue;
    end
    try
        [L.(name), L.info.(name)] = fn(pt, opts);
    catch ME
        L.(name) = 0;
        L.info.(name) = struct('applicable', false, 'error', ME.message);
    end
end

% Combined total over finite contributions (a -Inf would veto the point; NaN is
% treated as no contribution).
vals = [L.sn1987a, L.neutronstar, L.stellar, L.fifthforce, L.bullet, L.kaon, ...
        L.kphi, L.cmb, L.integral, L.mesoninv, L.beamdump];
L.total = sum(vals(isfinite(vals)));
if any(vals == -Inf)
    L.total = -Inf;
end
end

% ======================================================================
% Defaults
% ======================================================================
function o = default_opts()
% Per-constraint settings.  Smoothing widths sigma_* set the scale (in the
% natural log-variable of each constraint) over which the one-sided penalty ramps
% up past the bound; tune as needed.
here = fileparts(mfilename('fullpath'));
o = struct();
o.data_dir = fullfile(here, '..', 'data');

o.enable_sn1987a     = true;
o.enable_neutronstar = true;
o.enable_stellar     = true;
o.enable_fifthforce  = true;
o.enable_bullet      = true;
o.enable_kaon        = true;
o.enable_kphi        = true;
o.enable_cmb         = true;
o.enable_integral    = true;
o.enable_mesoninv    = true;
o.enable_beamdump    = true;

o.sigma_sn       = 0.3;   % log10 depth inside the SN1987A contour
o.sigma_ns       = 0.3;   % log10 depth above the neutron-star bound
o.sigma_stellar  = 0.3;   % log10(alpha) depth above the stellar-cooling bound
o.sigma_axial    = 0.3;   % log10 depth above the axial 5th-force bound
o.sigma_ff_pseudo = 0.3;  % log10 depth above the pseudoscalar 5th-force bounds
o.sigma_bullet   = 0.3;   % log10 depth above the bullet-cluster sigma/m bound
o.bound_bullet   = 1.0;   % sigma_self/m_chi bound [cm^2/g]
o.v_bullet       = 0.01;  % bullet-cluster relative velocity beta=v/c (natural units)
o.sigma_kaon     = 0.3;   % log10 depth above the K+ invisible-decay limit
o.bound_kaon     = 1.0;   % K -> pi chi chi rate / invisible limit; excluded above 1
o.sigma_kphi     = 0.3;   % log10 depth above the K+->pi+ phi (on-shell) limit
o.bound_kphi     = 1.0;   % K -> pi phi rate / invisible limit; excluded above 1
o.sigma_cmb      = 0.3;   % log10 depth above the Slatyer CMB <sigma v> bound
o.bound_cmb      = 1.0;   % <sigma v> / CMB limit; excluded above 1
o.sigma_integral = 0.3;   % log10 depth above the INTEGRAL <sigma v> bound
o.bound_integral = 1.0;   % <sigma v> / INTEGRAL limit; excluded above 1
o.sigma_mesoninv = 0.3;   % log10 depth above the pi0/eta/eta' -> invisible limits
o.bound_mesoninv = 1.0;   % max_j Gamma(j -> chi chi) / limit_j; excluded above 1
o.sigma_beamdump = 0.3;   % TODO (stub)
end

function o = fill_defaults(o, d)
% Copy any field of d missing from o (shallow).
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(o, f{i}) || isempty(o.(f{i}))
        o.(f{i}) = d.(f{i});
    end
end
end

% ======================================================================
% Shared likelihood / geometry helpers
% ======================================================================
function l = softregion(depth, sigma)
% Soft penalty for an exclusion region: 0 outside (depth <= 0), else a Gaussian
% in the depth *inside* the boundary.  sigma in the same units as depth.
if ~isfinite(depth) || sigma <= 0 || depth <= 0
    l = 0;
else
    l = -0.5 * (depth / sigma)^2;
end
end

function [l, inside, depth] = region_penalty(x, y, poly_x, poly_y, edge_x, edge_y, sigma)
% Soft exclusion for a point (x,y) against a closed polygon (poly) whose physical
% boundary is the polyline (edge).  inside = inpolygon test; depth = distance to
% the real boundary (edge), so artificial closing/cap edges never set the depth.
inside = inpolygon(x, y, poly_x, poly_y);
if inside
    depth = point_to_polyline_dist(x, y, edge_x, edge_y);
    l = softregion(depth, sigma);
else
    depth = 0;
    l = 0;
end
end

function d = point_to_polyline_dist(x, y, px, py)
% Minimum Euclidean distance from point (x,y) to the OPEN polyline (px,py).
n = numel(px);
d = inf;
for i = 1:n-1
    d = min(d, seg_dist(x, y, px(i), py(i), px(i+1), py(i+1)));
end
end

function d = seg_dist(x, y, ax, ay, bx, by)
% Distance from (x,y) to segment a-b.
abx = bx - ax;  aby = by - ay;
denom = abx^2 + aby^2;
if denom == 0
    t = 0;
else
    t = ((x - ax) * abx + (y - ay) * aby) / denom;
    t = max(0, min(1, t));
end
cx = ax + t * abx;  cy = ay + t * aby;
d = hypot(x - cx, y - cy);
end

function curve = load_curve(fpath)
% Load a 2-column digitized curve, cached per path.
persistent cache
if isempty(cache), cache = containers.Map('KeyType','char','ValueType','any'); end
if isKey(cache, fpath)
    curve = cache(fpath);  return;
end
curve = load(fpath);
cache(fpath) = curve;
end

function kind = nucleon_kind(interaction)
% Which nucleon coupling the interaction's quark vertex implies:
%   'scalar' (SS/PS -> hp), 'pseudo' (SP/PP -> htildeP), 'axial', or '' (none).
if ismember(interaction, {'heavyScalar','lightScalar','scalar', ...
                          'heavyScalarPS','lightScalarPS','scalarPS'})
    kind = 'scalar';
elseif ismember(interaction, {'heavyScalarSP','lightScalarSP','scalarSP', ...
                              'heavyScalarPP','lightScalarPP','scalarPP'})
    kind = 'pseudo';
elseif strcmp(interaction, 'axialVector')
    kind = 'axial';
else
    kind = '';
end
end

% ======================================================================
% SN1987A  --  closed exclusion contour (scalar OR pseudoscalar coupling)
% ======================================================================
function [l, info] = sn1987a_loglike(pt, opts)
% SN1987A cooling: the mediator is produced in the SN core via its nucleon
% coupling; too much emission shortens the neutrino burst.  Excluded region is a
% CLOSED band in (log10 mphi, log10 g_N) -- free-streaming below, trapping above.
% Two data sets by nucleon coupling type: scalar (SS/PS, g_N=hp, Hardy-Sokolov
% 2410.17347) and pseudoscalar (SP/PP, g_N=htildeP).  Contour extended
% indefinitely to the left in mphi (per user instruction).
info = struct('applicable', false, 'excluded', false, ...
              'log_mphi', log10(pt.mphi), 'log_gN', NaN, 'depth', 0);
l = 0;

[hp, htildeP] = nucleon_couplings(pt.g);
switch nucleon_kind(pt.interaction)
    case 'scalar', gN = hp;      fname = 'SN1987AScalarlogmphiloggN.dat';
    case 'pseudo', gN = htildeP; fname = 'SN1987AlogmphiloggN_pseudo.dat';
    otherwise, return;   % axial / unknown -> SN1987A not applicable here
end

log_gN = log10(abs(gN));
info.log_gN = log_gN;
if ~isfinite(log_gN), return; end   % g_N = 0 (e.g. g~[1 1 1] -> htildeP=0)

curve = load_curve(fullfile(opts.data_dir, fname));
x = log10(pt.mphi);  y = log_gN;

% Build the closed contour, extended to the LEFT: the closing edge (last->first
% vertex, the contour's left side) is pushed out to X_LEFT so the band continues
% flat for any mphi below the data.  X_LEFT is far below any scan floor, so the
% artificial vertical wall never sets the soft-edge depth.
cx0 = curve(:,1);  cy0 = curve(:,2);
X_LEFT = -30;
poly_x = [cx0; X_LEFT; X_LEFT];
poly_y = [cy0; cy0(end); cy0(1)];
% Real boundary for the depth = the data contour closed normally.
edge_x = [cx0; cx0(1)];  edge_y = [cy0; cy0(1)];

if x > max(cx0)          % beyond the heavy (right) end -> not constrained
    info.applicable = true;  return;
end
info.applicable = true;
[l, inside, depth] = region_penalty(x, y, poly_x, poly_y, edge_x, edge_y, opts.sigma_sn);
info.excluded = inside;  info.depth = depth;
end

% ======================================================================
% Neutron-star cooling  --  one-sided upper bound (scalar coupling)
% ======================================================================
function [l, info] = neutronstar_loglike(pt, opts)
% Neutron-star cooling bound on the SCALAR nucleon coupling (Damiano 2506.19906).
% Excluded ABOVE the curve (too-large g_N over-cools the NS).  Data columns:
% log10(mphi/eV), log10(g_N).  Applies to SS/PS.  The paper assumes equal
% phi-p / phi-n couplings and the emission is N N -> N N phi bremsstrahlung in a
% core that is ~90-95% neutrons, so we feed the NEUTRON coupling, g_N = hn.
% Conservative for hp ~= hn: unequal couplings open a dipole np channel the
% paper estimates is enhanced ~1e8 (not modelled).
info = struct('applicable', false, 'excluded', false, ...
              'log_mphi', log10(pt.mphi), 'log_gN', NaN, 'depth', 0);
l = 0;

if ~strcmp(nucleon_kind(pt.interaction), 'scalar'), return; end
[~, ~, ~, ~, hn] = nucleon_couplings(pt.g);
log_gN = log10(abs(hn));
info.log_gN = log_gN;
if ~isfinite(log_gN), return; end

curve = load_curve(fullfile(opts.data_dir, 'NeutronStarLoggNLogmphi.dat'));
% col1 = log10(mphi/eV) -> log10(mphi/GeV) = col1 - 9.
cx0 = curve(:,1) - 9;  cy0 = curve(:,2);
x = log10(pt.mphi);  y = log_gN;

info.applicable = true;
[l, inside, depth] = above_curve_penalty(x, y, cx0, cy0, opts.sigma_ns);
info.excluded = inside;  info.depth = depth;
end

% ======================================================================
% Stellar cooling (TRGB + horizontal branch, Hardy & Lasenby)  --  one-sided
% upper bound on alpha_phiNN = g_phiNN^2 / (4 pi), scalar coupling
% ======================================================================
function [l, info] = stellar_loglike(pt, opts)
% Combined stellar-cooling bound (tip of the red-giant branch + horizontal
% branch) on the SCALAR phi-nucleon coupling from Hardy & Lasenby, quoted as
% alpha_phiNN = g_phiNN^2/(4 pi) vs mphi under the assumption of EQUAL
% phi-proton and phi-neutron couplings.  Our quark-level model gives unequal
% hp, hn, so as a first pass we plug in the average g_bar = (hp + hn)/2 and
% alpha = g_bar^2/(4 pi).  Data columns: log10(mphi/eV), log10(alpha_phiNN).
% Excluded ABOVE the curve (extra energy loss from the stellar core over-cools
% the star / speeds its evolution off these branches).  Applies to SS/PS
% (scalar quark vertex).  The digitized data span mphi ~ 1e-7..4e-4 GeV; the
% bound is extended flat to the LEFT indefinitely (as SN1987A is) so a lighter
% mediator is still constrained, and drops to 0 for mphi above the data.
info = struct('applicable', false, 'excluded', false, ...
              'log_mphi', log10(pt.mphi), 'log_alpha', NaN, ...
              'hp', NaN, 'hn', NaN, 'depth', 0);
l = 0;

if ~strcmp(nucleon_kind(pt.interaction), 'scalar'), return; end
[hp, ~, ~, ~, hn] = nucleon_couplings(pt.g);
g_bar = 0.5 * (hp + hn);
log_alpha = log10(g_bar^2 / (4*pi));
info.hp = hp;  info.hn = hn;  info.log_alpha = log_alpha;
if ~isfinite(log_alpha), return; end   % g_bar = 0

curve = load_curve(fullfile(opts.data_dir, 'alphaphiNNvsmphieV.dat'));
% col1 = log10(mphi/eV) -> log10(mphi/GeV) = col1 - 9.
cx0 = curve(:,1) - 9;  cy0 = curve(:,2);
x = log10(pt.mphi);  y = log_alpha;

info.applicable = true;
[l, inside, depth] = above_curve_penalty(x, y, cx0, cy0, opts.sigma_stellar, true);
info.excluded = inside;  info.depth = depth;
end

% ======================================================================
% Fifth force / long-range  --  one-sided upper bound on a coupling product
% ======================================================================
function [l, info] = fifthforce_loglike(pt, opts)
% Long-range spin-dependent (5th-force) bound.  gx-independent: this is a force
% between two nucleons, so it depends only on (mphi, g).  Excluded ABOVE the
% envelope in every case.  Two families, selected by nucleon_kind:
%
%   'axial'  (axialVector, 2408.15691) -- one curve in
%            x = log10(lambda/m), y = log10(g_A^2), with g_A = a_N and the range
%            lambda = hbarc/m_A' (hbarc = 1.973e-16 GeV*m).
%
%   'pseudo' (SP/PP) -- two independent spin-spin curves in
%            x = log10(mphi/GeV), y = log10 of the constrained product:
%              PProtonProtonFifthForce.dat   ->  g_P^2      (proton-proton)
%              PProtonNeutronFifthForce.dat  ->  |g_P g_N|  (proton-neutron)
%            with g_P = htildeP, g_N = htildeN.  Their data files are digitized
%            in log10(mphi/eV), hence the -9 to reach GeV.  The two experiments
%            are independent, so their penalties ADD; info.pp / info.pn keep the
%            per-curve depth and excluded flag.
info = struct('applicable', false, 'excluded', false, 'kind', '', ...
              'log_lambda', NaN, 'log_gA2', NaN, 'depth', 0);
l = 0;

kind = nucleon_kind(pt.interaction);
info.kind = kind;
switch kind
    case 'axial',  [l, info] = fifthforce_axial(pt, opts, info);
    case 'pseudo', [l, info] = fifthforce_pseudo(pt, opts, info);
    otherwise,     return;   % scalar / unknown -> no 5th-force bound here
end
end

function [l, info] = fifthforce_axial(pt, opts, info)
% axialVector: |g_A|^2 = a_N^2 vs range lambda = hbarc/m_A' (2408.15691).
l = 0;
[~, ~, a_N] = nucleon_couplings(pt.g);
log_gA2 = log10(a_N^2);
hbarc_m = 1.97e-16;                  % GeV*m (= alltheconstants hbarc 1.97e-14 GeV*cm /100)
log_lambda = log10(hbarc_m / pt.mphi);
info.log_lambda = log_lambda;  info.log_gA2 = log_gA2;
if ~isfinite(log_gA2), return; end   % a_N = 0

curve = load_curve(fullfile(opts.data_dir, 'axialConstr.dat'));
cx0 = curve(:,1);  cy0 = curve(:,2);   % log10(lambda/m), log10(g_A^2)

info.applicable = true;
[l, inside, depth] = above_curve_penalty(log_lambda, log_gA2, cx0, cy0, opts.sigma_axial);
info.excluded = inside;  info.depth = depth;
end

function [l, info] = fifthforce_pseudo(pt, opts, info)
% SP/PP: pseudoscalar spin-spin bounds on g_P^2 and |g_P g_N| vs mphi.
[~, htildeP, ~, htildeN] = nucleon_couplings(pt.g);
x = log10(pt.mphi);
info.log_mphi = x;
info.log_gP   = log10(abs(htildeP));
info.log_gN   = log10(abs(htildeN));

% {sub-field, data file, log10 of the constrained product}.  htildeP = htildeN = 0
% on the singlet direction g ~ [1 1 1], so both products are then -Inf (vacuous).
specs = { 'pp', 'PProtonProtonFifthForce.dat',  log10(htildeP^2); ...
          'pn', 'PProtonNeutronFifthForce.dat', log10(abs(htildeP * htildeN)) };

l = 0;
for k = 1:size(specs, 1)
    name = specs{k, 1};
    y    = specs{k, 3};
    sub  = struct('applicable', false, 'excluded', false, 'log_prod', y, 'depth', 0);
    if isfinite(y)
        curve = load_curve(fullfile(opts.data_dir, specs{k, 2}));
        % col1 = log10(mphi/eV) -> log10(mphi/GeV) = col1 - 9.
        cx0 = curve(:,1) - 9;  cy0 = curve(:,2);
        sub.applicable = true;
        [lk, inside, depth] = above_curve_penalty(x, y, cx0, cy0, opts.sigma_ff_pseudo);
        sub.excluded = inside;  sub.depth = depth;
        l = l + lk;
        info.applicable = true;
        info.excluded   = info.excluded || inside;
        info.depth      = max(info.depth, depth);
    end
    info.(name) = sub;
end
end

function [l, inside, depth] = above_curve_penalty(x, y, cx0, cy0, sigma, extend_left)
% One-sided "excluded ABOVE the curve" within the data x-range, implemented as a
% closed region: the data curve capped at the top by Y_TOP.  This is robust to a
% non-monotonic x-trace (e.g. the near-vertical left edge of the axial envelope)
% that would break interp1, while the soft-edge depth is measured to the REAL
% curve only (not the cap edges).
%
% extend_left (optional, default false): continue the bound flat at cy0(1) for
% all x below the data, so a mediator lighter than the lightest digitized point
% is still constrained (as SN1987A does).  cx0 must be sorted ascending; the flat
% continuation IS part of the physical boundary, so it also sets the soft-edge
% depth to its left.
if nargin < 6 || isempty(extend_left), extend_left = false; end
Y_TOP = max(cy0) + 50;     % far above any physical model value
if extend_left
    X_LEFT = -30;          % far below any scan floor
    poly_x = [cx0; cx0(end); X_LEFT; X_LEFT];
    poly_y = [cy0; Y_TOP;    Y_TOP;  cy0(1)];
    edge_x = [X_LEFT; cx0];  edge_y = [cy0(1); cy0];
else
    poly_x = [cx0; cx0(end); cx0(1)];
    poly_y = [cy0; Y_TOP;     Y_TOP];
    edge_x = cx0;  edge_y = cy0;
end
[l, inside, depth] = region_penalty(x, y, poly_x, poly_y, edge_x, edge_y, sigma);
end

% ======================================================================
% Bullet cluster  --  DM self-interaction, one-sided upper bound on sigma/m
% ======================================================================
function [l, info] = bullet_loglike(pt, opts)
% Bullet-cluster DM SELF-INTERACTION bound.
%   Observable : sigma_self / m_chi [cm^2/g] for chi chi -> chi chi via mediator
%                exchange.  Depends on mx, mphi, gx (NOT g), velocity-dependent
%                (evaluated at opts.v_bullet).  Bound opts.bound_bullet ~1 cm^2/g;
%                excluded ABOVE.  The physics (sigma_self) is in sigma_self_over_m.
bound = opts.bound_bullet;
info = struct('applicable', false, 'excluded', false, ...
              'sigma_over_m', NaN, 'bound', bound);
l = 0;

sigma_over_m = sigma_self_over_m(pt, opts.v_bullet);   % [cm^2/g]
info.sigma_over_m = sigma_over_m;
if ~isfinite(sigma_over_m) || sigma_over_m <= 0
    return;   % formula not supplied / not applicable -> contributes 0
end

info.applicable = true;
depth = log10(sigma_over_m) - log10(bound);   % >0 when over the bound
l = softregion(depth, opts.sigma_bullet);
info.excluded = depth > 0;
end

function s = sigma_self_over_m(pt, v_rel)
% sigma_self/m_chi [cm^2/g] for chi chi -> chi chi via mediator exchange, Born
% transfer cross section, at relative velocity v_rel (beta = v/c).  Depends on
% (pt.mx, pt.mphi, pt.gx) and the phi-chi-chi (DM-side) vertex; independent of g.
% R = mx v / mphi; with x = R^2,
%   L(x) = log(1+x) - x/(1+x)
%   P(x) = [6(x+1)log(1+x) + (x^2-3x-6)x] / (x^2 (x+1))
%   scalar DM (SS/SP):  sigma_T = gx^4/(2 pi mx^2 v^4) * L
%   pseudo DM (PS/PP):  sigma_T = gx^4/(32 mx^2) * P      (~ R^4/2: p-suppressed)
%   axial   (AA):       sigma_T = 6 gx^4/(pi mx^2 v^4) * L
% Both brackets cancel at small x (L ~ x^2/2, P's numerator ~ x^4/2), and heavy
% mediators put x ~ 1e-10 and below, so small x uses their series.
hbarc    = 1.97e-14;    % GeV*cm   (alltheconstants.m)
gev_to_g = 1.783e-24;   % g/GeV    (alltheconstants.m)
mx = pt.mx;  mphi = pt.mphi;  gx = pt.gx;  v = v_rel;
x = (mx * v / mphi)^2;
switch dm_vertex_kind(pt.interaction)
    case 'scalar'
        sigma_T_GeV2 = gx^4 / (2*pi * mx^2 * v^4) * transfer_L(x);
    case 'pseudo'
        sigma_T_GeV2 = gx^4 / (32 * mx^2) * transfer_P(x);
    case 'axial'
        sigma_T_GeV2 = 6 * gx^4 / (pi * mx^2 * v^4) * transfer_L(x);
    otherwise
        s = NaN;  return;   % unknown vertex -> constraint inert
end
s = sigma_T_GeV2 * hbarc^2 / (mx * gev_to_g);   % GeV^-2 -> cm^2 -> cm^2/g
end

function L = transfer_L(x)
% log(1+x) - x/(1+x) = sum_{n>=2} (-1)^n (n-1)/n x^n.
if x < 0.1
    n = 2:24;
    L = sum((-1).^n .* (n-1)./n .* x.^n);
else
    L = log1p(x) - x/(1 + x);
end
end

function P = transfer_P(x)
% [6(x+1)log(1+x) + (x^2-3x-6)x] / (x^2 (x+1)); the numerator is
% 6 sum_{n>=4} (-1)^n x^n / (n(n-1)).
if x < 0.1
    n = 4:26;
    P = 6 * sum((-1).^n .* x.^(n-2) ./ (n.*(n-1))) / (x + 1);
else
    P = (6*(x + 1)*log1p(x) + (x^2 - 3*x - 6)*x) / (x^2 * (x + 1));
end
end

function kind = dm_vertex_kind(interaction)
% The phi-chi-chi (DM-side) vertex implied by the interaction name (the FIRST
% suffix letter): 'scalar' (SS/SP), 'pseudo' (PS/PP), 'axial', or '' (unknown).
% NB: distinct from nucleon_kind, which groups by the phi-SM (nucleon) vertex.
if ismember(interaction, {'heavyScalar','lightScalar','scalar', ...
                          'heavyScalarSP','lightScalarSP','scalarSP'})
    kind = 'scalar';
elseif ismember(interaction, {'heavyScalarPS','lightScalarPS','scalarPS', ...
                              'heavyScalarPP','lightScalarPP','scalarPP'})
    kind = 'pseudo';
elseif strcmp(interaction, 'axialVector')
    kind = 'axial';
else
    kind = '';
end
end

% ======================================================================
% Kaon invisible decay  --  K+ -> pi+ chi chi vs the K+->pi+ + invisible limit
% ======================================================================
function [l, info] = kaon_loglike(pt, opts)
% Rare-kaon invisible-decay bound: the mediator lets K+ -> pi+ chi chi, adding to
% the K+ -> pi+ + (invisible) width measured against the K+->pi+ nu nu limit.
%   Observable : ratio = Gamma(K+->pi+ chi chi) / (invisible-decay rate limit);
%                excluded when ratio > opts.bound_kaon (=1).  Depends on
%                (mx, mphi, gx, g) and the interaction; open only for
%                mx < (mK+- - mpi+-)/2 ~ 0.177 GeV.  Physics in kaon_invisible_ratio.m.
bound = opts.bound_kaon;
info = struct('applicable', false, 'excluded', false, 'ratio', NaN, 'bound', bound);
l = 0;

ratio = kaon_invisible_ratio(pt.mx, pt.mphi, pt.gx, pt.g, pt.interaction);
info.ratio = ratio;
if ~isfinite(ratio) || ratio <= 0
    return;   % channel closed / model has no meson-decay implementation -> 0
end

info.applicable = true;
depth = log10(ratio) - log10(bound);   % >0 when over the invisible limit
l = softregion(depth, opts.sigma_kaon);
info.excluded = depth > 0;
end

% ======================================================================
% Meson-decay & indirect-detection ratios (mesonandCMB_ratio.m)
% ======================================================================
function [l, info] = kphi_loglike(pt, opts)
% K+ -> pi+ phi (on-shell mediator) vs the K+->pi+ + invisible limit.  Open only
% for mphi < ~0.354 GeV; weighted by the chance the phi looks invisible (gx enters
% only there).  Excluded when ratio > opts.bound_kphi.
[l, info] = ratio_loglike(pt, opts, 'kphi', opts.bound_kphi, opts.sigma_kphi);
end

function [l, info] = cmb_loglike(pt, opts)
% DM annihilation <sigma v> vs the Slatyer CMB energy-injection bound.  ~gx^2;
% not applicable to SxS/SxP (velocity-suppressed).  Excluded when ratio > bound.
[l, info] = ratio_loglike(pt, opts, 'cmb', opts.bound_cmb, opts.sigma_cmb);
end

function [l, info] = integral_loglike(pt, opts)
% DM annihilation <sigma v> vs the INTEGRAL telescope bound.  ~gx^2; all models,
% light DM only.  Excluded when ratio > opts.bound_integral.
[l, info] = ratio_loglike(pt, opts, 'integral', opts.bound_integral, opts.sigma_integral);
end

function [l, info] = mesoninv_loglike(pt, opts)
% pi0 / eta / eta' -> chi chibar vs the invisible-decay limits, through the
% meson's mixing with the mediator.  SxP/PxP and axialVector only; open for
% mx < m_j/2.  The most constraining of the three is used (ratio = max_j).
[l, info] = ratio_loglike(pt, opts, 'mesoninv', opts.bound_mesoninv, opts.sigma_mesoninv);
end

function [l, info] = ratio_loglike(pt, opts, which, bound, sigma)
% Shared soft-mask loglike for a mesonandCMB_ratio constraint (mirrors kaon_loglike).
info = struct('applicable', false, 'excluded', false, 'ratio', NaN, 'bound', bound);
l = 0;
ratio = mesonandCMB_ratio(pt.mx, pt.mphi, pt.gx, pt.g, pt.interaction, which);
info.ratio = ratio;
if ~isfinite(ratio) || ratio <= 0
    return;   % constraint closed / not applicable to this interaction -> 0
end
info.applicable = true;
depth = log10(ratio) - log10(bound);   % >0 when over the bound
l = softregion(depth, sigma);
info.excluded = depth > 0;
end

% ======================================================================
% Terrestrial (beam dump)  --  STUB (no data)
% ======================================================================
function [l, info] = beamdump_loglike(pt, opts) %#ok<INUSD>
% Terrestrial beam-dump bound: accelerator production + decay of the mediator.
%   Observable : mediator-SM coupling vs mphi; exclusion in (mphi, coupling).
% TODO(user): supply the effective production coupling and the exclusion data.
l = 0;
info = struct('applicable', false, 'excluded', false, 'value', NaN, 'bound', NaN);
end
