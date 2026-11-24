% compare_limits_g.m
% Instead of solving for gx that reproduces Omega h^2 = 0.12 at fixed
% (g, mphi), fix gx and scan a scalar
% normalization `g_q` on the quark-coupling vector g (g_eff = g_q*g).
%
% Effect on observables (heavyScalar / lightScalar branch):
%   Dscript ~ g_q^2  ->  Gamma(phi->gg) ~ g_q^2, sigma(gg<->xx) ~ g_q^2
%   Gamma(phi->pipi/KK) ~ g_q^2
%   nucleon coupling fp ~ g_q  ->  sigma_chiN ~ g_q^2 * gx^2
% Non-resonant: Oh2_FI ~ g_q^4, Oh2_FO ~ 1/g_q^4 (rate dominated).
% Resonant: the BW saturates the g_q dependence; root-find handles either.
tic
alltheconstants;

% -----------------------------------------------------------------------
% Model parameters -- edit these
% -----------------------------------------------------------------------
interaction = 'heavyScalar';  % PxS: kaon+kphi+cmb+integral+SN1987A all apply.
% NB use 'heavyScalarPS'/'lightScalarPS', NOT bare 'scalarPS' -- svavg's 'gg'
% channel whitelist omits the bare name, so 'scalarPS' gives an empty <sigma v> figure.
mphi        = 0.1;         % mediator mass [GeV] -- < 0.354 so K+->pi+ phi (kphi) is open
gx          = 1;             % FIXED dark coupling
STYLE       = 'classic';     % figure style: 'compact' (decade ticks, no grid,
                             %   tight axes, doc-scaled fonts) or 'classic'
                             %   (grid, auto ticks, loose axes, 16/14/13/11 pt).
                             %   Assigned here, NOT guarded on exist(), so it is
                             %   reset every run -- a leftover 'classic' in the
                             %   workspace must not silently restyle a later run,
                             %   the same reason batch_params is consumed.
SIGMA_YCAP  = 1e-25;         % 'compact' only: hard ceiling on the sigma panel's
                             %   y axis [cm^2].  A light mediator pushes
                             %   sigma_chiN up by many decades, so a row that
                             %   shares one range with it is nearly all empty
                             %   space above the curves.  Inf disables.  Forced
                             %   to Inf under 'classic', whose ranges are the
                             %   ones the published figures were tuned against.
FLOAT_ALIGN = 'bottom';      % floating region labels (Overabundant): 'horizontal' or
                             % 'bottom' (parallel to the fill's lower edge)
LABEL_POS   = {};            % manual inline-label placement on the sigma panel:
                             %   {label, [mx sigma], rot_deg, hide; ...} in data units
                             %   (rot, hide optional; hide=true removes the label).
                             %   Drivers read these from label_pos_g.json.  Each listed label is moved there AFTER
                             %   inline_labels runs; everything else stays
                             %   automatic.  For the rare panel the placer gets
                             %   wrong.  Reset every run, like STYLE.
ghat       = [2 -1 -1];  % quark-coupling DIRECTION (gq multiplies this); asymmetric
Oh2_target  = 0.12;
T_RH        = 0.01;          % reheat temperature [GeV]
production  = 'sm';          % chi production from the bath: 'sm' (gamma gamma + meson pairs) or 'gg' (photons only)

% mx_vec = logspace(-6, log10(0.200), 200);   % uniform in log, 200 pts
% Same range, 80 points instead of 200, warped so the sampling is sparse at low
% mx and dense in the upper half of the decade span -- that is where the relic
% locus turns over, the resonance/threshold features sit, and the constraint
% boundaries cross.  1-(1-u)^p packs points toward the TOP of the range (plain
% u^p would pack them at the bottom).  Spacing runs 0.13 dex at the low end
% down to 0.001 dex at the top; 56 of the 80 points sit above 10^-3.35.
mx_dec_lo = -5;              % lightest mass, as log10(m/GeV)
mx_dec_hi = log10(0.200);    % heaviest   "    "   "
n_mx      = 80;              % grid points
mx_warp   = 2;               % 1 = plain log spacing; >1 = pack toward mx_dec_hi

u      = linspace(0, 1, n_mx);            % 0 at the light end, 1 at the heavy end
w      = 1 - (1 - u).^mx_warp;            % same endpoints, bunched near u = 1
mx_vec = 10.^(mx_dec_lo + (mx_dec_hi - mx_dec_lo) * w);

% Plot toggles -- each draws one figure; set false to skip.
plot_sigma_chiN = true;     % Fig 1: sigma_chiN vs m_chi
plot_gq        = true;     % Fig 2: gq vs m_chi
plot_svgg       = true;     % Fig 3: <sigma v>_{gamma gamma -> chi chi} vs m_chi (at T = m_chi/20)

% Small-mixing-breakdown overlay (axialVector): sin(theta) threshold above which
% the leading-order A'-meson mixing expansion is suspect.  Drawn as a boundary +
% shade only -- it does NOT gate the FI/FO/loop scans (see axial_mixing_floor.m).
s_mix = 0.3;

% -----------------------------------------------------------------------
% Constraint toggles -- which overlays are DRAWN, on every panel.  Turning
% one off only skips the drawing; the quantity behind it is still computed
% and still written to the .dat, EXCEPT the five ratio constraints
% (kaon/kphi/cmb/integral/mesoninv), whose root-find is skipped when hidden.
% Batch override: batch_params.show, a struct holding any subset of these
% fields, e.g. batch_params.show = struct('cmb', false, 'integral', false).
% -----------------------------------------------------------------------
SHOW = struct();
% -- annotation ---------------------------------------------------------
SHOW.ylabel        = true;    % false -> drop the y label and y tick labels
                              %   (inner panels of a shared-y-axis row)
SHOW.legend        = false;    % boxed legend (set false to rely on inline labels)
SHOW.inline_labels = true;   % text labels placed on the curves/regions
% -- relic curves -------------------------------------------------------
SHOW.relic_fi     = false;   % Freeze-in
SHOW.relic_fo     = false;   % Freeze-out
SHOW.relic_loop   = true;   % self-consistent loop, BOTH branches (Omega h^2 = 0.12)
SHOW.overabundant = true;   % Omega h^2 > 0.12 fill between the loop branches
% -- theory / validity --------------------------------------------------
SHOW.pert_line    = true;   % gq = 1 boundary
SHOW.nonpert      = true;   % shaded gq > 1 region
SHOW.unitarity    = true;   % axialVector: Gamma_A' >= m_A' ceiling
SHOW.large_mixing = true;   % axialVector: sin(theta_A') > s_mix
SHOW.mediator_relic = false; % stable mediator (axialVector A') overcloses: Omega_A' h^2 > 0.12.
                            % Off by default: below min(2mchi, 3mpi) the A' is stable only
                            % in this EFT; a UV completion generically makes it decay.
                            % Still computed and written to the .dat.
% -- cosmology ----------------------------------------------------------
SHOW.neff         = true;   % Delta-Neff band (see neff_color_mode below)
% -- nucleon / astrophysical (external_loglikes) ------------------------
SHOW.sn1987a      = true;
SHOW.neutronstar  = true;
SHOW.stellar      = true;   % stellar cooling (TRGB + HB, Hardy & Lasenby); SS/PS
SHOW.fifthforce   = true;   % 5th force (axial g_A^2; pseudoscalar g_P^2, g_P g_N)
SHOW.bullet       = true;   % Bullet Cluster sigma_T/m_chi (excludes m_chi columns)
% -- meson decays / indirect detection ----------------------------------
SHOW.kaon         = true;   % K+ -> pi+ chi chi (invisible)
SHOW.kphi         = true;   % K+ -> pi+ phi
SHOW.cmb          = true;   % CMB energy injection (Slatyer)
SHOW.integral     = true;   % INTEGRAL <sigma v>
SHOW.mesoninv     = true;   % pi0 / eta / eta' -> invisible (SP/PP, axialVector)

% -----------------------------------------------------------------------
% Batch override.  A driver script (run_figure_sets_g.m) can set the struct
% `batch_params` in the workspace before invoking this script; each field
% present replaces the parameter of the same name set above, and anything
% absent keeps the value above.  Interactive runs never define it, so this
% whole block is a no-op there.  Extra (non-parameter) fields recognised:
%   figdir     -- directory to save the figures into (enables batch saving)
%   figtag     -- basename stem for the saved figures
%   figfmt     -- cellstr of extensions, e.g. {'pdf','png'} (default {'pdf'})
%   figsize    -- saved page size [w h] in cm (default [18 13])
%   close_figs -- close each figure after saving it (default true)
%   dat_tag    -- suffix appended to the .dat filename
%
% batch_params is CONSUMED here: the non-parameter fields are copied into the
% BATCH struct and batch_params is cleared, so it cannot survive into a later
% interactive run.  Without that, a driver interrupted before its own
% `clear batch_params` would leave the struct in the base workspace and the
% next hand-run `compare_limits_g` would silently save and close its figures.
% -----------------------------------------------------------------------
% figsize: the saved page is the figure's own size in centimetres, so this is
% the knob that decides how big the exported figure actually is.  [18 13] cm
% ~ 7 x 5 in, the size the 16/14/13 pt fonts below are chosen for.  Left at the
% on-screen window size instead, a Retina display exports ~3.4 x 2.7 in.
BATCH = struct('save', false, 'close', true, 'figsize', [18 13], 'ylims', struct());
% Y-limits this run CHOSE, per panel.  run_figure_sets_g reads these back to
% build one shared range per model row, then feeds them in as BATCH.ylims.
PANEL_YLIM = struct();
if exist('batch_params', 'var') && isstruct(batch_params)
    if isfield(batch_params, 'interaction'),     interaction     = batch_params.interaction;     end
    if isfield(batch_params, 'mphi'),            mphi            = batch_params.mphi;            end
    if isfield(batch_params, 'gx'),              gx              = batch_params.gx;              end
    if isfield(batch_params, 'ghat'),           ghat           = batch_params.ghat;           end
    if isfield(batch_params, 'Oh2_target'),      Oh2_target      = batch_params.Oh2_target;      end
    if isfield(batch_params, 'T_RH'),            T_RH            = batch_params.T_RH;            end
    if isfield(batch_params, 'production'),      production      = batch_params.production;      end
    if isfield(batch_params, 'mx_vec'),          mx_vec          = batch_params.mx_vec;          end
    if isfield(batch_params, 'plot_sigma_chiN'), plot_sigma_chiN = batch_params.plot_sigma_chiN; end
    if isfield(batch_params, 'plot_gq'),        plot_gq        = batch_params.plot_gq;        end
    if isfield(batch_params, 'plot_svgg'),       plot_svgg       = batch_params.plot_svgg;       end
    if isfield(batch_params, 's_mix'),           s_mix           = batch_params.s_mix;           end
    % show: merge field-by-field so a partial struct keeps the defaults above.
    if isfield(batch_params, 'show') && isstruct(batch_params.show)
        bpf = fieldnames(batch_params.show);
        for kbp = 1:numel(bpf)
            if ~isfield(SHOW, bpf{kbp})
                warning('compare_limits_g:unknownShow', ...
                    'batch_params.show.%s is not a known constraint toggle -- ignored.', bpf{kbp});
            else
                SHOW.(bpf{kbp}) = logical(batch_params.show.(bpf{kbp}));
            end
        end
    end

    % Non-parameter fields -> BATCH, then drop batch_params (see above).
    if isfield(batch_params, 'figdir')
        BATCH.save   = true;
        BATCH.figdir = batch_params.figdir;
        if isfield(batch_params, 'figtag')
            BATCH.figtag = batch_params.figtag;
        else
            BATCH.figtag = sprintf('compare_limits_g_%s_mphi%.4g_gx%.4g', ...
                                   interaction, mphi, gx);
        end
        if isfield(batch_params, 'figfmt'), BATCH.figfmt = batch_params.figfmt;
        else,                               BATCH.figfmt = {'pdf'};
        end
    end
    if isfield(batch_params, 'figsize'),    BATCH.figsize = batch_params.figsize;          end
    % Stashed rather than applied here: the FONT block below needs BATCH.figsize,
    % which is not final until this whole override loop has run.
    if isfield(batch_params, 'font'),       batch_params_font_ = batch_params.font;        end
    if isfield(batch_params, 'ylims'),      BATCH.ylims = batch_params.ylims;             end
    if isfield(batch_params, 'style'),      STYLE = batch_params.style;                   end
    if isfield(batch_params, 'sigma_ycap'), SIGMA_YCAP = batch_params.sigma_ycap;      end
    if isfield(batch_params, 'close_figs'), BATCH.close = logical(batch_params.close_figs); end
    if isfield(batch_params, 'dat_tag'),    dat_tag     = batch_params.dat_tag;             end
    if isfield(batch_params, 'replot'),     REPLOT      = batch_params.replot;              end
    if isfield(batch_params, 'label_pos'),  LABEL_POS   = batch_params.label_pos;           end
    if isfield(batch_params, 'float_align'), FLOAT_ALIGN = batch_params.float_align;       end
    clear batch_params
end

% sigma_chiN uses the FULL angular-integrated form (full_sigma_name maps a heavy*
% relic model to its 'scalar*' full sibling; axialVector/light*/scalar* pass through).
interaction_sig = full_sigma_name(interaction);

logT0 = log(T0);

% -----------------------------------------------------------------------
% Unitarity / EFT-validity ceiling (axialVector only): above gq_ceiling the
% meson->chichi width exceeds the meson mass (or the A'-meson mixing mass
% reaches m_A'), so the meson-mixing EFT is invalid (axial_couplings.m
% hard-errors there).  In the g_q scan gx is FIXED and g_q scales g, so the
% ceiling lives in g_q (axial_gq_ceiling, the companion of axial_gx_ceiling).
% Mapped to sigma_chiN via DMnucxsec.  Computed FIRST (before FI/FO/loop) so
% every coupling scan can cap its bracket at it: at small m_A' the relic roots
% sit far below the fixed scan floors but still inside the EFT-valid region
% (g_q < ceiling), so an uncapped [-15,..] scan misses them entirely.
% -----------------------------------------------------------------------
gq_ceiling = nan(size(mx_vec));
sg_ceiling  = nan(size(mx_vec));
if strcmp(interaction, 'axialVector')
    for j = 1:length(mx_vec)
        try
            ec = axial_gq_ceiling(mphi, mx_vec(j), gx, ghat);
            gq_ceiling(j) = ec;
            if isfinite(ec)
                sg_ceiling(j) = DMnucxsec(mx_vec(j), mphi, gx, ec*ghat, interaction_sig) * hbarc^2;
            end
        catch
        end
    end
end

% -----------------------------------------------------------------------
% Small-mixing-breakdown boundary (axialVector only): g_q above which the
% largest A'-meson mixing angle sin(theta_j) = g_q*|m_{jA'}|/m_A' exceeds s_mix,
% so the leading-order mixing expansion behind lambda_{jchi} ~ 1/m_A'^2 is
% suspect.  Distinct from the unitarity ceiling above (gq_unit ~ mphi^2 vs
% gq_mix ~ mphi; they cross near mphi ~ 0.01 GeV).  gq_mix is mx-independent
% (mixing masses are linear in g); mapped to sigma_chiN via DMnucxsec.  Drawn as
% a drawn overlay; it does not cap the scan.  No-op for non-axialVector
% (gq_mix NaN).
% -----------------------------------------------------------------------
gq_mix = NaN;
sg_mix  = nan(size(mx_vec));
if strcmp(interaction, 'axialVector')
    try
        gq_mix = axial_mixing_floor(mphi, mx_vec(1), ghat, s_mix);
    catch
        gq_mix = NaN;
    end
    if isfinite(gq_mix)
        for j = 1:length(mx_vec)
            try
                sg_mix(j) = DMnucxsec(mx_vec(j), mphi, gx, gq_mix*ghat, interaction_sig) * hbarc^2;
            catch
            end
        end
    end
end

% -----------------------------------------------------------------------
%% Font sizing -- specified at FINAL size in the finished document
% -----------------------------------------------------------------------
% \includegraphics[width=...] scales the exported page, fonts and all, so what
% the reader measures is  font_pt * (panel width / export width).  At the old
% hardcoded 16 pt exported to 18 cm and placed in a 5.63 cm column that is
% 16*5.63/18 = 5.0 pt -- the reason the panels read as tiny.
%
% So state what each element should MEASURE in the document and scale it up by
% the shrink the \includegraphics will apply.  Set FONT.panel_cm to the width
% one panel is placed at (0.32\textwidth on a 2 cm-margin letter page = 5.63 cm;
% \the\textwidth in your document divided by the panels per row).  Export at
% FONT.panel_cm itself and the figure is WYSIWYG -- fonts come out exactly as
% typed and nothing is scaled at all.
FONT = struct('label', 10, 'tick', 9, 'legend', 8, 'inline', 6.2, 'panel_cm', 5.63);
if exist('batch_params_font_', 'var') && isstruct(batch_params_font_)
    f_ = fieldnames(batch_params_font_);
    for k_ = 1:numel(f_), FONT.(f_{k_}) = batch_params_font_.(f_{k_}); end
    clear batch_params_font_ f_ k_
end
% ---- Style preset -----------------------------------------------------
% 'compact' (default): decade ticks, no grid, axes tightened to their labels,
%           fonts scaled so they measure FONT.* points in the finished
%           document.  Packs the most plot into a multi-panel row.
% 'classic': the look the inline-label placement was originally tuned against
%           -- grid on, MATLAB's automatic log ticks, axes left at their default
%           loose inset, and ABSOLUTE font points (16/14/13/11) at the figure's
%           own size.  Everything the compact style does takes space away from
%           label placement, so keep this available for comparison.
% Set via batch_params.style / RUNSET.style.
STYLE_classic = strcmpi(STYLE, 'classic');
if STYLE_classic, SIGMA_YCAP = Inf; end   % never re-scale the classic figures

fscale_ = BATCH.figsize(1) / FONT.panel_cm;
FS = struct('label',  FONT.label  * fscale_, 'tick',   FONT.tick   * fscale_, ...
            'legend', FONT.legend * fscale_, 'inline', FONT.inline * fscale_);
if STYLE_classic
    % Absolute points at the figure's own size -- NOT scaled to a placement
    % width.  These are the sizes the label engine was tuned against.
    FS = struct('label', 16, 'tick', 14, 'legend', 13, 'inline', 11);
end

% -----------------------------------------------------------------------
%% Computation cache -- replot without recomputing
% -----------------------------------------------------------------------
% Everything from here down to the .dat write is the expensive part (root
% finds, Boltzmann solves, Delta-Neff band sweeps).  It is cached to a .mat
% beside the .dat, so a pure PRESENTATION change -- legend, colours, SHOW
% toggles -- re-renders in seconds instead of re-solving.  Set
% batch_params.replot = true (RUNSET.replot in the drivers) to load rather
% than compute.  Only numeric/logical/char variables are cached, which
% deliberately leaves out SHOW, BATCH and the plot_* toggles: those are
% re-derived on every run, so a replot picks up your new toggles while
% reusing the old numbers.
if ~exist('dat_tag', 'var'), dat_tag = ''; end
if ~exist('REPLOT',  'var'), REPLOT  = false; end
% REPLOT = 'ext' is the middle gear: reuse the cached relic / Delta-Neff
% numbers, but RE-RUN the external nucleon bands (ext_specs) and write them
% back into the .mat.  That is what you want after adding or changing an
% external constraint -- the expensive blocks are untouched, the band scan
% costs seconds.  REPLOT = true reuses everything.
% REPLOT = 'neff' reuses only the cached RELIC loci (freeze-in / freeze-out /
% loop g_q and sigma, and the perturbativity line) and recomputes everything
% downstream of them -- Delta-Neff, <sigma v>_gg, the external bands -- then
% re-saves.  For a Delta-Neff change: the relic solves are the expensive (and,
% for light-mediator SP, memory-hungry) part and do not depend on it.
REFRESH_NEFF = false;
if ischar(REPLOT) || isstring(REPLOT)
    switch lower(char(REPLOT))
        case 'ext',  REFRESH_EXT = true;  REPLOT = true;
        case 'neff', REFRESH_EXT = false; REPLOT = false; REFRESH_NEFF = true;
        otherwise
            error('compare_limits_g:replotMode', ...
                  'replot must be true/false, ''ext'' or ''neff'' (got ''%s'').', char(REPLOT));
    end
else
    REFRESH_EXT = false; REPLOT = logical(REPLOT);
end
if ~isfolder('output'), mkdir('output'); end
cache_fname = sprintf('output/%s%s.mat', sprintf('compare_limits_g_%s_mphi%.4g_gx%.4g', interaction, mphi, gx), dat_tag);
if (REPLOT || REFRESH_NEFF) && ~isfile(cache_fname)
    error('compare_limits_g:noCache', ['replot requested but %s does not exist -- run ' ...
          'once with replot=false to build it.'], cache_fname);
end
% External-band constants.  Defined ABOVE the replot branch on purpose: they
% are a cell array plus SHOW-derived flags, so the cache deliberately does not
% carry them (CACHE_SKIP, and the save loop keeps only numeric/logical/char).
% Left inside the compute block they are simply undefined on the replot path,
% and any panel that actually draws an external band -- every mphi = 1e-9 one,
% where SN1987A / neutron-star / 5th-force bands are finite -- dies there.
ext_specs = { 'sn1987a',     'SN1987A',           [0.85 0.20 0.55]; ...
              'neutronstar', 'Neutron star',      [0.30 0.30 0.75]; ...
              'stellar',     'Stellar cooling',   [0.55 0.20 0.60]; ...
              'fifthforce',  '5th force',         [0.85 0.40 0.10] };
nE = size(ext_specs, 1);
% Drawn only where toggled on -- the bands are still SCANNED either way, so
% the "excludes this point" warning below fires regardless of visibility.
ext_show = [SHOW.sn1987a, SHOW.neutronstar, SHOW.stellar, SHOW.fifthforce];

% Island walls (set in the loop block, or loaded); a leftover from the previous
% panel must not survive a load from a cache that predates them.
clear relic_tip_mx
if ~REPLOT
if REFRESH_NEFF
    % Relic loci from the cache (see REPLOT = 'neff' above); the grid must match.
    RELIC_CACHED = load(cache_fname, 'mx_vec', 'gqvec_fi', 'gqvec_fo', ...
        'gqvec_loop_fi', 'gqvec_loop_fo', 'sglim_fi', 'sglim_fo', ...
        'sglim_loop_fi', 'sglim_loop_fo', 'sglim_pert');
    if ~isequal(RELIC_CACHED.mx_vec, mx_vec)
        error('compare_limits_g:neffReplotGrid', ...
              '%s was built on a different mx grid; rerun with replot = false.', cache_fname);
    end
    gqvec_fi      = RELIC_CACHED.gqvec_fi;       gqvec_fo      = RELIC_CACHED.gqvec_fo;
    gqvec_loop_fi = RELIC_CACHED.gqvec_loop_fi;  gqvec_loop_fo = RELIC_CACHED.gqvec_loop_fo;
    sglim_fi       = RELIC_CACHED.sglim_fi;        sglim_fo       = RELIC_CACHED.sglim_fo;
    sglim_loop_fi  = RELIC_CACHED.sglim_loop_fi;   sglim_loop_fo  = RELIC_CACHED.sglim_loop_fo;
    sglim_pert     = RELIC_CACHED.sglim_pert;
    MED_CACHED = load(cache_fname);   % the mediator band rides along with the loci
    if isfield(MED_CACHED, 'relic_tip_mx'), relic_tip_mx = MED_CACHED.relic_tip_mx; end
    if isfield(MED_CACHED, 'gq_med_lo')
        gq_med_lo = MED_CACHED.gq_med_lo;  gq_med_hi = MED_CACHED.gq_med_hi;
        sg_med_lo  = MED_CACHED.sg_med_lo;   sg_med_hi  = MED_CACHED.sg_med_hi;
    else
        [gq_med_lo, gq_med_hi, sg_med_lo, sg_med_hi] = deal(nan(size(mx_vec)));
    end
    clear MED_CACHED
    fprintf('[compare_limits_g] Relic loci from %s (not re-solved)\n', cache_fname);
else
% -----------------------------------------------------------------------
% Freeze-in -- root-find on Oh2_FI(g_q) = Oh2_target.
% -----------------------------------------------------------------------
sglim_fi   = nan(size(mx_vec));
gqvec_fi  = nan(size(mx_vec));

t_fi = tic;
parfor j = 1:length(mx_vec)
    m = mx_vec(j);
    flog_FI = @(le) Oh2_FI(m, mphi, gx, 10^le * ghat, interaction, T_RH, Mpl, s0, rhoc, T0, production) - Oh2_target;
    % Cap the scan at the EFT ceiling and scan the valid decades below it.  At
    % small m_A' the FI root sits ~10 decades below the old fixed -15 floor but
    % still inside the valid region (g_q < ceiling); the fixed grid missed it
    % entirely.  No-op for non-axialVector (gq_ceiling NaN -> the [-15,5] grid).
    le_lo = -15; le_hi = 5; n_fi = 14;
    capj = gq_ceiling(j);
    if isfinite(capj) && capj > 0
        le_hi = log10(capj) - 0.05;   % stay just below the ceiling
        le_lo = le_hi - 20;           % reach well past the relic root
        n_fi  = 28;                   % denser: Oh2 ~ gq^4 crossing is steep
    end
    le_grid = linspace(le_lo, le_hi, n_fi);
    Oh_grid = nan(size(le_grid));
    for kk = 1:length(le_grid)
        try
            v = Oh2_FI(m, mphi, gx, 10^le_grid(kk) * ghat, interaction, T_RH, Mpl, s0, rhoc, T0, production);
            if isfinite(v), Oh_grid(kk) = v; end
        catch err
            fprintf(2, '[FI scan] mx=%.4g le=%.2f: %s\n', m, le_grid(kk), err.message);
        end
    end
    sgrid = Oh_grid - Oh2_target;
    gqc = NaN;
    for kk = 1:length(le_grid)-1
        if isnan(sgrid(kk)) || isnan(sgrid(kk+1)), continue; end
        if sgrid(kk)*sgrid(kk+1) < 0
            try
                lec = fzero(flog_FI, [le_grid(kk), le_grid(kk+1)], optimset('TolX',1e-3));
                gqc = 10^lec;
                break;
            catch err
                fprintf(2, '[FI fzero] mx=%.4g bracket=[%.2f,%.2f]: %s\n', ...
                        m, le_grid(kk), le_grid(kk+1), err.message);
            end
        end
    end
    gqvec_fi(j) = gqc;
    if ~isnan(gqc)
        sglim_fi(j) = DMnucxsec(m, mphi, gx, gqc * ghat, interaction_sig) * hbarc^2;
    end
end

fprintf('[compare_limits_g] FI block: %.2f s (%.2f s/mass)\n', toc(t_fi), toc(t_fi)/length(mx_vec));

% -----------------------------------------------------------------------
% Freeze-out -- root-find on freeze_out(g_q) = Oh2_target.
% -----------------------------------------------------------------------
sglim_fo  = nan(size(mx_vec));
gqvec_fo = nan(size(mx_vec));

t_fo = tic;
parfor k = 1:length(mx_vec)
    m = mx_vec(k);
    flog_FO = @(le) freeze_out(m, mphi, gx, 10^le * ghat, interaction) - Oh2_target;
    % Densified vs the gx-scan: with rate ~ g_q^4 (not g_q^2), the
    % sub/over-abundant transition is sharper in log10(g_q).  A coarse
    % grid easily misses the sign change.  Capped at the EFT ceiling like the
    % FI/loop scans (no-op for non-axialVector).
    le_lo = -12; le_hi = 4; n_fo = 20;
    capj = gq_ceiling(k);
    if isfinite(capj) && capj > 0
        le_hi = log10(capj) - 0.05;
        le_lo = le_hi - 18;
        n_fo  = 28;
    end
    % Cheap coarse pre-scan to kill hopeless masses before the dense grid.
    % Each freeze_out solves a stiff ode15s (+200-pt svavg_grid), so the dense
    % grid is ~28 ODE solves/mass.  At small m_A' the WIMP mechanism never
    % thermalizes -> Oh2 sits orders of magnitude from 0.12 across the whole
    % valid g_q range with no crossing.  Scan a handful of points; only run the
    % dense grid if the coarse scan crosses, or comes within ~2 decades of, the
    % target.  Skip (NaN) otherwise.  ~9 ODE solves instead of ~28+fzero.
    le_coarse = linspace(le_lo, le_hi, 9);
    Ohc = nan(size(le_coarse));
    for kk = 1:numel(le_coarse)
        try
            v = freeze_out(m, mphi, gx, 10^le_coarse(kk) * ghat, interaction);
            if isfinite(v), Ohc(kk) = v; end
        catch err
            fprintf(2, '[FO coarse] mx=%.4g le=%.2f: %s\n', m, le_coarse(kk), err.message);
        end
    end
    finc = ~isnan(Ohc);
    if ~any(finc)
        fprintf(2, '[FO no-root] mx=%.4g: all NaN on coarse grid\n', m);
        continue;   % gqvec_fo(k) stays NaN
    end
    sc_coarse = sign(Ohc(finc) - Oh2_target);
    has_sc  = any(diff(sc_coarse(:)) ~= 0);                       % crosses target
    near_tg = min(abs(log10(Ohc(finc) / Oh2_target))) < 2;        % within ~2 decades
    if ~has_sc && ~near_tg
        fprintf(2, '[FO no-root] mx=%.4g: coarse Oh2 in [%.3g,%.3g], far from target -> skip dense\n', ...
            m, min(Ohc(finc)), max(Ohc(finc)));
        continue;   % gqvec_fo(k) stays NaN
    end
    le_grid = linspace(le_lo, le_hi, n_fo);
    Oh_grid = nan(size(le_grid));
    for kk = 1:length(le_grid)
        try
            v = freeze_out(m, mphi, gx, 10^le_grid(kk) * ghat, interaction);
            if isfinite(v), Oh_grid(kk) = v; end
        catch err
            fprintf(2, '[FO scan] mx=%.4g le=%.2f: %s\n', m, le_grid(kk), err.message);
        end
    end
    sgrid = Oh_grid - Oh2_target;
    gqc = NaN;
    for kk = 1:length(le_grid)-1
        if isnan(sgrid(kk)) || isnan(sgrid(kk+1)), continue; end
        if sgrid(kk)*sgrid(kk+1) < 0
            try
                lec = fzero(flog_FO, [le_grid(kk), le_grid(kk+1)], optimset('TolX',1e-3));
                gqc = 10^lec;
                break;
            catch err
                fprintf(2, '[FO fzero] mx=%.4g bracket=[%.2f,%.2f]: %s\n', ...
                        m, le_grid(kk), le_grid(kk+1), err.message);
            end
        end
    end
    if isnan(gqc)
        % Diagnose why: count finite, dump min/max of Oh2 across grid.
        finite_mask = ~isnan(Oh_grid);
        if any(finite_mask)
            fprintf(2, '[FO no-root] mx=%.4g: Oh2 range over le in [%.2f,%.2f] = [%.3g, %.3g] (%d/%d finite)\n', ...
                m, le_grid(1), le_grid(end), min(Oh_grid(finite_mask)), max(Oh_grid(finite_mask)), ...
                sum(finite_mask), length(le_grid));
        else
            fprintf(2, '[FO no-root] mx=%.4g: all NaN on scan grid\n', m);
        end
    end
    gqvec_fo(k) = gqc;
    if ~isnan(gqc)
        sglim_fo(k) = DMnucxsec(m, mphi, gx, gqc * ghat, interaction_sig) * hbarc^2;
    end
end

fprintf('[compare_limits_g] FO block: %.2f s (%.2f s/mass)\n', toc(t_fo), toc(t_fo)/length(mx_vec));

% -----------------------------------------------------------------------
% Self-consistent loop -- two roots (FIMP branch / WIMP branch) in g_q.
% -----------------------------------------------------------------------
sglim_loop_fi   = nan(size(mx_vec));
sglim_loop_fo   = nan(size(mx_vec));
gqvec_loop_fi  = nan(size(mx_vec));
gqvec_loop_fo  = nan(size(mx_vec));
% Stable-mediator overclosure (freeze_coupled only; in practice the axialVector
% A' below min(2 mx, 3 m_pi)): the g_q interval where Omega_med h^2 > Oh2_target.
gq_med_lo      = nan(size(mx_vec));
gq_med_hi      = nan(size(mx_vec));

if mphi > 50 * T_RH
    solver_name = 'freeze_loop';
else
    solver_name = 'freeze_coupled';
end

t_loop = tic;
parfor j = 1:length(mx_vec)
    solver_local = @freeze_loop;
    switch solver_name
        case 'freeze_loop',    solver_local = @freeze_loop;
        case 'freeze_coupled', solver_local = @freeze_coupled;
    end
    m = mx_vec(j);
    flog = @(le) solver_local(m, mphi, gx, 10^le * ghat, interaction, T_RH, false, production) - Oh2_target;
    % Wide bracket on BOTH sides: at fixed gx, production scales as g_q^4
    % (two vertices through D ~ g_q^2) but depletion through SM final states
    % scales only as g_q^2 (one vertex on the SM side; the dark side is
    % already saturated by gx).  So the WIMP-falling branch falls half as
    % fast in log(g_q) as the FIMP-rising branch climbs, and the WIMP root
    % for Oh2=0.12 sits ~2x as many decades above the peak as the FIMP root
    % sits below.  Peak Oh2(g_q) for sub-GeV mx at gx=1 can be enormous,
    % pushing the WIMP root past g_q ~ 10^10.
    % Cap the upper end at the EFT ceiling (axialVector): above it the solver
    % only returns NaN, and an uncapped scan wastes points / risks a degenerate
    % bracket.  No-op for non-axialVector (gq_ceiling NaN); keeps 40 points so
    % the no-ceiling density is unchanged.
    % Floor tied to the FI root, not fixed: the FIMP branch drifts down with mx,
    % and a fixed floor loses it once log10(gq_FI) falls below.  With only the
    % WIMP root left in the window the proximity selection below puts it in the
    % fi slot and leaves fo empty, which reads as a 13-order jump in gq_loop_fi
    % and a vanishing lower branch.  gqvec_fi is complete before this loop.
    le_lo = -16; le_hi = 20;
    if isfinite(gqvec_fi(j)) && gqvec_fi(j) > 0
        le_lo = min(le_lo, log10(gqvec_fi(j)) - 3);
    end
    capj = gq_ceiling(j);
    if isfinite(capj) && capj > 0
        le_hi = min(le_hi, log10(capj) - 0.05);   % stay just below the ceiling
        le_lo = min(le_lo, le_hi - 14);            % always span 14 valid decades
        if le_hi <= le_lo, le_hi = le_lo + 1; end  % degenerate guard
    end
    % 48 points (was 40): the window is wider now, this holds the bracketing
    % density at ~0.75 dex.
    le_grid = linspace(le_lo, le_hi, 48);
    Oh_grid = nan(size(le_grid));
    OA_grid = nan(size(le_grid));     % freeze_coupled's mediator relic (stable A' only)
    for kk = 1:length(le_grid)
        try
            if strcmp(solver_name, 'freeze_coupled')
                [v, finf] = freeze_coupled(m, mphi, gx, 10^le_grid(kk) * ghat, interaction, T_RH, false, production);
                if finf.phi_stable, OA_grid(kk) = finf.Omega_phi_h2; end
            else
                v = solver_local(m, mphi, gx, 10^le_grid(kk) * ghat, interaction, T_RH, false, production);
            end
            if isfinite(v), Oh_grid(kk) = v; end
        catch err
            fprintf(2, '[Loop scan] mx=%.4g le=%.2f: %s\n', m, le_grid(kk), err.message);
        end
    end
    s = Oh_grid - Oh2_target;
    roots_found = nan(1, length(le_grid)-1);
    n_roots     = 0;
    for kk = 1:length(le_grid)-1
        if isnan(s(kk)) || isnan(s(kk+1)), continue; end
        if s(kk)*s(kk+1) < 0
            try
                lec = fzero(flog, [le_grid(kk), le_grid(kk+1)], optimset('TolX',1e-3));
                n_roots = n_roots + 1;
                roots_found(n_roots) = 10^lec;
            catch err
                fprintf(2, '[Loop fzero] mx=%.4g bracket=[%.2f,%.2f]: %s\n', ...
                        m, le_grid(kk), le_grid(kk+1), err.message);
            end
        end
    end
    roots_found = roots_found(1:n_roots);
    [gq_med_lo(j), gq_med_hi(j)] = overclosed_interval(le_grid, OA_grid, Oh2_target);

    % Near the FIMP/WIMP merge the two roots close to well under the coarse
    % grid step, and the wide window makes that step large: 48 points over the
    % default 36 decades is 0.77 dex, against a root separation that shrinks
    % through ~0.8 dex as the branches approach.  Both roots then land in one
    % segment (no sign change at all) or the pair is half-detected, and what
    % survives is a spurious root sitting BETWEEN the true ones -- drawn as the
    % FIMP branch continuing up and to the right past the point where the two
    % branches should have closed.
    %
    % So when the coarse pass finds fewer than two roots but both standalone
    % anchors exist, re-scan finely across the span they bracket.  It costs a
    % second sweep only on the masses that need it, and the window is a few
    % decades instead of 36.
    if n_roots < 2 && isfinite(gqvec_fi(j)) && gqvec_fi(j) > 0 ...
                   && isfinite(gqvec_fo(j)) && gqvec_fo(j) > 0
        lf = log10(gqvec_fi(j)); lo_ = log10(gqvec_fo(j));
        fine_lo = min(lf, lo_) - 1.0;
        fine_hi = max(lf, lo_) + 1.0;
        fine_grid = linspace(fine_lo, fine_hi, 40);
        Oh_fine = nan(size(fine_grid));
        for kk = 1:length(fine_grid)
            try
                v = solver_local(m, mphi, gx, 10^fine_grid(kk) * ghat, interaction, T_RH, false, production);
                if isfinite(v), Oh_fine(kk) = v; end
            catch
            end
        end
        sf = Oh_fine - Oh2_target;
        rf = nan(1, length(fine_grid)-1); nf = 0;
        for kk = 1:length(fine_grid)-1
            if isnan(sf(kk)) || isnan(sf(kk+1)), continue; end
            if sf(kk)*sf(kk+1) < 0
                try
                    lec = fzero(flog, [fine_grid(kk), fine_grid(kk+1)], optimset('TolX',1e-3));
                    nf = nf + 1; rf(nf) = 10^lec;
                catch
                end
            end
        end
        if nf > n_roots
            roots_found = rf(1:nf);
            n_roots     = nf;
        end
    end
    if n_roots < 2
        % Surface the Oh2(g_q) profile so we can see whether the WIMP
        % branch root is missing because (a) the peak doesn't reach 0.12,
        % (b) both roots fell inside one segment, or (c) the grid is too
        % narrow.  Reports peak Oh2 + the le at which it occurs.
        finite_mask = ~isnan(Oh_grid);
        if any(finite_mask)
            [oh_max, ipk] = max(Oh_grid);
            fprintf(2, '[Loop n_roots=%d] mx=%.4g: peak Oh2 = %.3g at le=%.2f; Oh2[%.2f]=%.3g  Oh2[%.2f]=%.3g\n', ...
                n_roots, m, oh_max, le_grid(ipk), ...
                le_grid(1),   Oh_grid(1), ...
                le_grid(end), Oh_grid(end));
        end
    end
    % Root selection.  Oh2(g_q) can have 4 roots: physical FIMP / WIMP pair,
    % plus a second pair on the BW-broadened "secluded" branch once Gphi/mphi
    % grows past 1.  Min/max picks the outer two and discards the physical
    % WIMP root.  Pick by proximity in log10(g_q) to the standalone FI / FO
    % solution instead; fall back to min/max if standalone has no root.
    if ~isempty(roots_found)
        le_roots = log10(roots_found);
        if ~isnan(gqvec_fi(j))
            [~, ii] = min(abs(le_roots - log10(gqvec_fi(j))));
            gqvec_loop_fi(j) = roots_found(ii);
        else
            gqvec_loop_fi(j) = min(roots_found);
        end
        sglim_loop_fi(j) = DMnucxsec(m, mphi, gx, gqvec_loop_fi(j) * ghat, interaction_sig) * hbarc^2;
        % The WIMP root is ALWAYS above the FIMP root (Oh2(g_q) rises, peaks,
        % falls), so search only the roots above the FIMP pick.  Without this
        % restriction both proximity searches can return the SAME root: the
        % standalone FO root sits BETWEEN the two loop branches (gq_FO ~ 1e2
        % against loop roots ~1e-9 and ~1e13), ~11 decades from each, so the
        % nearest-root test is a coin flip decided in the third decimal and
        % flips to the FIMP root once gq_FO drops below the geometric mean of
        % the two.  That collapsed loop_fo onto loop_fi over a wide mx band --
        % a zero-height "Overabundant" patch (invisible fill), a single drawn
        % branch, and a spurious vertical connector where the tie flipped back.
        % numel(roots_found) >= 2 only guarantees two roots EXIST, not that two
        % DIFFERENT ones were selected.
        % No gq_FO anchor -> leave the branch NaN rather than substituting a
        % root.  gq_FO goes NaN exactly where the WIMP mechanism stops having
        % a solution (min_gq Oh2 rises above 0.12), and there the only root
        % left above the FIMP one is the far secluded root ~10 decades higher
        % -- a qualitatively different branch.  Drawing it as the WIMP curve
        % put a spurious low-mass segment plus a vertical connector in the
        % panel.  (This ties the branch to the standalone FO solve, so a
        % SPURIOUS gq_FO failure would drop a real loop root too; the loop is
        % parfor, so anchoring on the neighbouring mass instead is not an
        % option here.)
        above = find(roots_found > gqvec_loop_fi(j));
        if ~isempty(above) && ~isnan(gqvec_fo(j))
            [~, ii] = min(abs(le_roots(above) - log10(gqvec_fo(j))));
            gqvec_loop_fo(j) = roots_found(above(ii));
            sglim_loop_fo(j) = DMnucxsec(m, mphi, gx, gqvec_loop_fo(j) * ghat, interaction_sig) * hbarc^2;
        end
    end
end

fprintf('[compare_limits_g] Loop block (%s): %.2f s (%.2f s/mass)\n', solver_name, toc(t_loop), toc(t_loop)/length(mx_vec));
sg_med_lo = nan(size(mx_vec));  sg_med_hi = nan(size(mx_vec));
for j = find(isfinite(gq_med_lo))
    sg_med_lo(j) = DMnucxsec(mx_vec(j), mphi, gx, gq_med_lo(j)*ghat, interaction_sig) * hbarc^2;
    sg_med_hi(j) = DMnucxsec(mx_vec(j), mphi, gx, gq_med_hi(j)*ghat, interaction_sig) * hbarc^2;
end
if any(isfinite(gq_med_lo))
    fprintf('[compare_limits_g] Stable mediator overcloses at %d of %d masses\n', ...
            sum(isfinite(gq_med_lo)), numel(mx_vec));
end

% Close the FIMP/WIMP coalescence gap.
% Only at a GENUINE merge: the two branches must already have converged to
% within COALESCE_MAX_DEC decades.  "next mass has no WIMP root" alone is not
% enough -- an isolated solver failure (one NaN from a stiff ODE abort) meets
% that test too, and snapping a still-separated pair to their geometric mean
% then throws a ~10-decade spike across the panel: the WIMP branch plunges and
% the FIMP branch leaps to meet it at a coupling neither ever takes.
%
% KNOWN LIMITATION: the merge point drawn here is cosmetic.  Within
% COALESCE_MAX_DEC the pair can still be far from merged, so the midpoint need
% not be a root, and the lobe can close at too low an mx.  PxP, gx = 0.01:
% the tip is at 181 MeV, while the roots at mphi = 1 keV are still 1.3e-2 and
% 3.4e-2 at 180 MeV and merge near 194 MeV.  A correct merge needs a fine scan
% for where peak Oh2(g_q) = 0.12.
COALESCE_MAX_DEC = 2;
both = ~isnan(gqvec_loop_fi) & ~isnan(gqvec_loop_fo);
snap_at = NaN;
for j = 1:length(mx_vec)-1
    if both(j) && ~both(j+1) && ...
       abs(log10(gqvec_loop_fo(j) / gqvec_loop_fi(j))) <= COALESCE_MAX_DEC
        gq_mid = sqrt(gqvec_loop_fi(j) * gqvec_loop_fo(j));
        gqvec_loop_fi(j) = gq_mid;
        gqvec_loop_fo(j) = gq_mid;
        sg_mid = DMnucxsec(mx_vec(j), mphi, gx, gq_mid * ghat, interaction_sig) * hbarc^2;
        sglim_loop_fi(j)  = sg_mid;
        sglim_loop_fo(j)  = sg_mid;
        if isnan(snap_at), snap_at = j; end
    end
end

% The two branches have closed at snap_at, so nothing survives beyond it.  A
% lone root at a larger mass is the tangency of a peak that has already dropped
% to Oh2 ~ 0.12, or a leftover of the coarse scan; either way the proximity test
% files it under whichever slot it happens to sit nearer, and it draws as a
% branch shooting up and to the right from a curve that should have terminated.
% Blank both branches past the merge so the locus actually closes.
if ~isnan(snap_at) && snap_at < length(mx_vec)
    tail = (snap_at+1):length(mx_vec);
    gqvec_loop_fi(tail) = NaN;  gqvec_loop_fo(tail) = NaN;
    sglim_loop_fi(tail)  = NaN;  sglim_loop_fo(tail)  = NaN;
end

% Close the LOW-mass end of each island (2026-11-18).  With the dark sector
% filled, chi chi <-> phi phi at fixed gx sets Oh2 by dark freeze-out: Oh2(g_q) is
% a plateau, flat to 4 digits over ~5 decades of g_q and ~ mx^2 (mphi = 1e-9,
% gx = 0.01).  As mx falls the
% whole plateau drops through 0.12 at once, so the Oh2 = 0.12 contour is a
% near-vertical wall between two grid masses.  The column-wise root find sees
% only the branches appear, ~6 decades apart, at the first mass above it, and
% the island was drawn open.  Where a run of two-branch masses starts with the
% branches separated, bisect (6 steps, in log mx) for the mass at which Oh2 at
% the branches' mid coupling crosses 0.12; the plot puts a vertical segment
% there (relic_tip_mx(j), NaN = no wall).  If Oh2 is still above 0.12 at the
% previous mass (the lower branch carries on, open-topped), there is no wall.
% A lone lower root just below the wall is the tip of the small bump at the
% plateau's onset (an overabundant sliver ~0.25 dex wide in g_q whose upper
% root the gqvec_fo anchor rejects); its fill would run up to g_q = 1, so it
% is dropped.
relic_tip_mx = relic_island_tips(mx_vec, gqvec_loop_fi, gqvec_loop_fo, COALESCE_MAX_DEC, ...
    relic_oh2_at(solver_name, mphi, gx, ghat, interaction, T_RH, production), Oh2_target);
for j = find(isfinite(relic_tip_mx))
    k = j - 1;
    while k >= 1 && ~isnan(gqvec_loop_fi(k)) && isnan(gqvec_loop_fo(k))
        gqvec_loop_fi(k) = NaN;  sglim_loop_fi(k) = NaN;  k = k - 1;
    end
    fprintf('[compare_limits_g] Relic island closes at mx = %.4g GeV (grid %.4g-%.4g)\n', ...
            relic_tip_mx(j), mx_vec(j-1), mx_vec(j));
end

% -----------------------------------------------------------------------
% Perturbativity boundary: g_q = 1 (g_eff = ghat).  sigma ~ g_q^2 at fixed
% gx, so points above this line need g_q > 1 -- coupling enhanced past the
% reference normalization.
% -----------------------------------------------------------------------
sglim_pert = zeros(size(mx_vec));
for j = 1:length(mx_vec)
    sglim_pert(j) = DMnucxsec(mx_vec(j), mphi, gx, ghat, interaction_sig) * hbarc^2;
end
end   % REFRESH_NEFF: relic blocks skipped

% -----------------------------------------------------------------------
% Delta-Neff exclusion band: a |dNeff - mu| > nsig*sigma cut on the SIGNED
% total dNeff (the 2-sigma region of the same two-sided Gaussian likelihood
% used in the GAMBIT-light scan), scanned in g_q at fixed gx; the edges are
% mapped through DMnucxsec to sigma_chiN.
% -----------------------------------------------------------------------
% Bands = the nsig-sigma exclusion of the SAME two-sided Gaussian likelihood
% used in the GAMBIT-light scan (neff_loglike): excluded where
% |dNeff - mu| > nsig*sigma, per epoch.  Measured: BBN -0.14 +/- 0.21
% (arXiv:2401.15054), CMB -0.054 +/- 0.17 (Planck N_eff = 2.99 +/- 0.17 vs
% SM 3.044).  dNeff is SIGNED, so the cut is two-sided.
neff_mu_bbn  = -0.14;  neff_sig_bbn = 0.21;
neff_mu_cmb  = -0.054; neff_sig_cmb = 0.17;
neff_nsig    = 2;      % exclusion CL drawn (2-sigma)
thr_pos_bbn = neff_mu_bbn + neff_nsig*neff_sig_bbn;   % ~ +0.28
thr_neg_bbn = neff_mu_bbn - neff_nsig*neff_sig_bbn;   % ~ -0.56
thr_pos_cmb = neff_mu_cmb + neff_nsig*neff_sig_cmb;   % ~ +0.29
thr_neg_cmb = neff_mu_cmb - neff_nsig*neff_sig_cmb;   % ~ -0.39
% -10, not -8: the Delta-Neff thermalisation onset sits at log10(g_q) = -8.038
% at mphi = 0.1 GeV, gx = 1 -- 0.04 decades BELOW a window floor of -8.  Every
% scanned coupling then comes back excluded, allon fires, and the band is drawn
% as a full-height column to the axis floor: a bottom edge that is the axis, not
% the physics.  Two decades of margin put the real lower edge inside the window,
% where the bisection can find it.
le_scan_lo = -10;
le_scan_hi = 8;
% ...but the Delta-Neff scan stops at PERTURBATIVITY.  le_scan_hi = 8 is a
% DRAWING choice -- it exists so the exclusion patch fills the panel rather than
% stopping mid-axis -- and g_q = 1 already means O(1) quark couplings, so eight
% decades above it is not a coupling the EFT describes.  Sending the solver
% there is not free: the band search probes its window ENDS first, so every mass
% paid a solve at g_q = 1e8 before anything else, in the regime where the dark
% sector is maximally thermalised and the integration is worst-behaved.
% MEASURED at mphi = 1e-9 GeV, mx = 1e-4 GeV, g_q = 1e8: 8.8 s and an ode15s
% "unable to meet integration tolerances" warning, for Delta-Neff = 18.4.
%
% The picture is unchanged because the drawing already handles it: a band whose
% upper edge is NaN means "open above" and is filled up to cap_fn, which is the
% perturbativity / EFT cap.  So stopping the SCAN at g_q = 1 and reporting the
% band as open-above draws the same patch without ever solving above it.
% le_scan_hi still sets the external-constraint window and the sg/g_q/svgg caps,
% which is why this is a separate knob rather than a smaller le_scan_hi.
le_neff_hi = 0;
n_neff     = 40;
% The EXTERNAL-constraint band scan (SN1987A / neutron star / 5th force, further
% down) gets its OWN window, deliberately decoupled from the Delta-Neff one
% above.  Those bands are set by a nucleon coupling times g_q, and the
% pseudoscalar h~_P carries an m_N/m_q chiral enhancement of ~500 per unit g_q,
% so e.g. the SP/PP SN1987A band at mphi = 1e-9 GeV sits at
% g_q ~ [7.5e-12, 9.1e-9] -- entirely below le_scan_lo, which is why it never
% appeared on the light-mediator panels.  Widening le_scan_lo itself would
% instead spread n_neff's fixed 40-point budget over 12 extra decades and
% coarsen the Delta-Neff bands, so the two windows are kept separate.
le_ext_lo  = -20;
le_ext_hi  = le_scan_hi;
n_ext      = 1 + round((le_ext_hi - le_ext_lo) / 0.04);   % ~0.04 dex, as before
% TOGGLE, in increasing detail:
%   'combined'  -> ONE band, one legend entry: excluded by BBN *or* CMB, either
%                  sign.  The default: the four signed/epoch bands and the two
%                  per-species edge lines all constrain the same g_q interval,
%                  so splitting them out mostly clutters the legend.
%   'signed'    -> per epoch, separate colours for dNeff_total>+max and <-max,
%                  plus the chi-only / phi-only edge lines.
%   'magnitude' -> per epoch, one band each for |dNeff_total-mu|>bound, plus
%                  the same edge lines.
neff_color_mode = 'combined';
if ~SHOW.neff, neff_color_mode = 'off'; end   % SHOW.neff = false -> draw nothing

% Delta-Neff is neff_simple.m: a CMB-epoch number only, so the BBN-epoch
% bands below stay NaN and only the CMB / combined bands are filled.
use_bbn = false;
use_cmb = true;
neff_want_bbn = {'tot','tot_pos','tot_neg'};

[gq_chi_lo, gq_chi_hi, gq_phi_lo, gq_phi_hi, gq_tot_lo, gq_tot_hi] = deal(nan(size(mx_vec)));
[sg_chi_lo,  sg_chi_hi,  sg_phi_lo,  sg_phi_hi,  sg_tot_lo,  sg_tot_hi ] = deal(nan(size(mx_vec)));
[gq_pos_lo, gq_pos_hi, gq_neg_lo, gq_neg_hi] = deal(nan(size(mx_vec)));
[sg_pos_lo,  sg_pos_hi,  sg_neg_lo,  sg_neg_hi ] = deal(nan(size(mx_vec)));
[allon_chi,  allon_phi,  allon_tot, allon_pos, allon_neg] = deal(false(size(mx_vec)));
% CMB-epoch counterparts (signed total + magnitude total).
[gq_tot_lo_c, gq_tot_hi_c, gq_pos_lo_c, gq_pos_hi_c, gq_neg_lo_c, gq_neg_hi_c] = deal(nan(size(mx_vec)));
[sg_tot_lo_c, sg_tot_hi_c, sg_pos_lo_c, sg_pos_hi_c, sg_neg_lo_c, sg_neg_hi_c] = deal(nan(size(mx_vec)));
[allon_tot_c, allon_pos_c, allon_neg_c] = deal(false(size(mx_vec)));
% Combined single band ('combined' mode): excluded by EITHER epoch, either sign.
[gq_all_lo, gq_all_hi, sg_all_lo, sg_all_hi] = deal(nan(size(mx_vec)));
allon_all = false(size(mx_vec));
% Second layer of the combined band.  dNeff can change sign along g_q (Primakoff:
% positive while phi leaves as late free radiation, negative once it stays coupled
% through recombination, or chi heating vs phi release), so the excluded g_q set
% at one mass can be two intervals with an allowed strip between.  Layer 1 is the
% lowest interval, layer 2 the next (any further ones merged into it, warned).
[gq_all_lo2, gq_all_hi2, sg_all_lo2, sg_all_hi2] = deal(nan(size(mx_vec)));
[gq_all_lo3, gq_all_hi3, sg_all_lo3, sg_all_hi3] = deal(nan(size(mx_vec)));

t_phi = tic;
% parfor for the same reason the other blocks use it: the cache is
% iteration-local (no cross-mass reuse), every temporary is written before it is
% read within the iteration, and every output is sliced by j.  This was the only
% serial block left.
parfor j = 1:length(mx_vec)
    m  = mx_vec(j);
    % One Delta-Neff call per g_q, memoized, both epochs extracted.
    cache = containers.Map('KeyType', 'double', 'ValueType', 'any');
    ndinfo = @(le) nd_cached(cache, mphi, m, gx, 10^le * ghat, interaction, T_RH);
    % Cap the g_q scan at the EFT ceiling: above it the rates trip
    % axial_couplings (NaN), so the fixed [le_scan_lo, le_scan_hi] window can
    % sit entirely in the invalid region (axialVector at small mphi) and return
    % empty bands.  Scan the EFT-valid decades below the ceiling instead.  No-op
    % for non-axialVector (gq_ceiling NaN).
    lo_j = le_scan_lo; hi_j = min(le_scan_hi, le_neff_hi);
    capj = gq_ceiling(j);
    if isfinite(capj) && capj > 0
        hi_j = min(hi_j, log10(capj) - 0.05);   % stay just below the ceiling
        lo_j = min([lo_j, hi_j - 12, -25]);      % reach down past the relic curves
    end
    % An epoch that is switched off is not evaluated at all.
    if use_bbn
        E = neff_band_edges(@(le) ndinfo(le).bbn, lo_j, hi_j, n_neff, ...
                            neff_mu_bbn, neff_sig_bbn, neff_nsig, neff_want_bbn);
    else
        E = neff_edges_off();
    end
    % CMB epoch: only the total bands are read below (chi/phi edge lines are
    % BBN-only), so skip refining the unused chi/phi CMB bands.
    if use_cmb
        Ec = neff_band_edges(@(le) ndinfo(le).cmb, lo_j, hi_j, n_neff, ...
                             neff_mu_cmb, neff_sig_cmb, neff_nsig, ...
                             {'tot','tot_pos','tot_neg'});
    else
        Ec = neff_edges_off();
    end
    % Combined band: one exclusion covering both epochs and both signs.  The
    % test statistic is the WORST pull, z = max_epoch |dNeff-mu|/sigma, so
    % "excluded" is z > nsig -- exactly the union of the four signed bands.
    % Reuses neff_band_edges (and its bisection refinement) by handing it z as
    % the "total" column against a standard normal (mu=0, sigma=1); the chi/phi
    % columns are unused, hence NaN.
    Ea = neff_band_edges(@(le) [NaN, NaN, neff_pull(ndinfo(le), ...
                            neff_mu_bbn, neff_sig_bbn, neff_mu_cmb, neff_sig_cmb, ...
                            use_bbn, use_cmb, true)], ...
                         lo_j, hi_j, n_neff, 0, 1, neff_nsig, {'tot'});

    sgmap = @(gq) DMnucxsec(m, mphi, gx, gq*ghat, interaction_sig) * hbarc^2;
    sg1 = @(e) ternan(e, sgmap);

    gq_chi_lo(j)=E.chi(1); gq_chi_hi(j)=E.chi(2);
    gq_phi_lo(j)=E.phi(1); gq_phi_hi(j)=E.phi(2);
    gq_tot_lo(j)=E.tot(1); gq_tot_hi(j)=E.tot(2);
    gq_pos_lo(j)=E.tot_pos(1); gq_pos_hi(j)=E.tot_pos(2);
    gq_neg_lo(j)=E.tot_neg(1); gq_neg_hi(j)=E.tot_neg(2);
    allon_chi(j)=E.allon.chi; allon_phi(j)=E.allon.phi; allon_tot(j)=E.allon.tot;
    allon_pos(j)=E.allon.tot_pos; allon_neg(j)=E.allon.tot_neg;

    gq_tot_lo_c(j)=Ec.tot(1); gq_tot_hi_c(j)=Ec.tot(2);
    gq_pos_lo_c(j)=Ec.tot_pos(1); gq_pos_hi_c(j)=Ec.tot_pos(2);
    gq_neg_lo_c(j)=Ec.tot_neg(1); gq_neg_hi_c(j)=Ec.tot_neg(2);
    allon_tot_c(j)=Ec.allon.tot; allon_pos_c(j)=Ec.allon.tot_pos; allon_neg_c(j)=Ec.allon.tot_neg;

    % A band reaching the top of the (perturbativity-capped) scan is OPEN ABOVE,
    % not closed at g_q = 1: everything above is excluded too, it just is not a
    % coupling we solve at.  NaN is how draw_neff_band says "fill to cap_fn".
    % Split into disjoint intervals.
    if isfield(Ea, 'iv'), ivj = Ea.iv.tot; else, ivj = Ea.tot; end
    if Ea.allon.tot || isempty(ivj), ivj = nan(1, 2); end
    top_j = 10^hi_j*(1 - 1e-9);
    ivj(isfinite(ivj(:,2)) & ivj(:,2) >= top_j, 2) = NaN;
    allon_all(j) = Ea.allon.tot;
    % Sort the intervals into three FIXED roles, so a layer means the same
    % thing at every mass and draw_neff_band never joins a negative interval at
    % one mass to a positive one at the next (bottom-up filling did that, and
    % merging intervals 2..N painted over the allowed strip between a positive
    % interval and the negative one above it):
    %   layer 1  closed dNeff < 0 interval below any positive one (heating)
    %   layer 2  closed dNeff > 0 interval (late free-streaming phi)
    %   layer 3  everything above: the open-above interval, or a closed
    %            negative one that sits above a positive one
    % The sign is read off one (memoized) dNeff call at the interval midpoint.
    % Two intervals landing in one role are merged, with a warning.
    lay = nan(3, 2);
    seen_pos = false;
    for q = 1:size(ivj, 1)
        if ~isfinite(ivj(q,1)), continue; end
        open_q = ~isfinite(ivj(q,2));
        if open_q
            k = 3;
        else
            zq = neff_pull(ndinfo(0.5*(log10(ivj(q,1)) + log10(ivj(q,2)))), ...
                           neff_mu_bbn, neff_sig_bbn, neff_mu_cmb, neff_sig_cmb, ...
                           use_bbn, use_cmb, true);
            if zq > 0,         k = 2; seen_pos = true;
            elseif ~seen_pos,  k = 1;
            else,              k = 3;
            end
        end
        if isfinite(lay(k,1))
            fprintf(['[compare_limits_g] m_chi=%.4g: two dNeff intervals in layer %d, ' ...
                     'merged\n'], m, k);
            lay(k,2) = ivj(q,2);
        else
            lay(k,:) = ivj(q,:);
        end
    end
    gq_all_lo(j)=lay(1,1);  gq_all_hi(j)=lay(1,2);
    gq_all_lo2(j)=lay(2,1); gq_all_hi2(j)=lay(2,2);
    gq_all_lo3(j)=lay(3,1); gq_all_hi3(j)=lay(3,2);

    sg_chi_lo(j)=sg1(E.chi(1)); sg_chi_hi(j)=sg1(E.chi(2));
    sg_phi_lo(j)=sg1(E.phi(1)); sg_phi_hi(j)=sg1(E.phi(2));
    sg_tot_lo(j)=sg1(E.tot(1)); sg_tot_hi(j)=sg1(E.tot(2));
    sg_pos_lo(j)=sg1(E.tot_pos(1)); sg_pos_hi(j)=sg1(E.tot_pos(2));
    sg_neg_lo(j)=sg1(E.tot_neg(1)); sg_neg_hi(j)=sg1(E.tot_neg(2));
    sg_tot_lo_c(j)=sg1(Ec.tot(1)); sg_tot_hi_c(j)=sg1(Ec.tot(2));
    sg_pos_lo_c(j)=sg1(Ec.tot_pos(1)); sg_pos_hi_c(j)=sg1(Ec.tot_pos(2));
    sg_neg_lo_c(j)=sg1(Ec.tot_neg(1)); sg_neg_hi_c(j)=sg1(Ec.tot_neg(2));
    sg_all_lo(j)=sg1(gq_all_lo(j)); sg_all_hi(j)=sg1(gq_all_hi(j));
    sg_all_lo2(j)=sg1(gq_all_lo2(j)); sg_all_hi2(j)=sg1(gq_all_hi2(j));
    sg_all_lo3(j)=sg1(gq_all_lo3(j)); sg_all_hi3(j)=sg1(gq_all_hi3(j));
end
fprintf('[compare_limits_g] Delta-Neff block (BBN+CMB): %.2f s\n', toc(t_phi));

% -----------------------------------------------------------------------
% External constraints (external_loglikes.m) as g_q-BANDS.  Here g = g_q*ghat varies
% and the nucleon couplings are linear in g_q, so SN1987A / neutron-star /
% 5th-force each exclude a fixed g_q-INTERVAL (mphi & ghat fixed ->
% mx-independent).  The min/max below assume ONE contiguous interval per
% constraint, which holds because every nucleon coupling (and hence every
% product the 5th-force curves bound) is monotonic in g_q.
% Derive each interval from external_loglikes by scanning g_q once over
% [le_ext_lo, le_ext_hi] (its OWN window -- see the note there), then remap the
% edges onto the sigma_chiN / <sigma v>_gg axes per-mx, exactly like the
% Delta-Neff bands.  {field, label, colour}:
% -----------------------------------------------------------------------
% Only the ext_specs fields are read here; the rest of external_loglikes
% (kaon / kphi / CMB / INTEGRAL, drawn separately below, kaon being a phase-space
% integral) is switched off so the wider, denser grid costs no more than the old
% narrow one -- it is in fact ~7x cheaper.
[ext_gq_lo, ext_gq_hi, ext_sg_lo, ext_sg_hi, ext_vg_lo, ext_vg_hi] = ...
    ext_bands(ext_specs, mx_vec, mphi, gx, ghat, interaction, interaction_sig, ...
              hbarc, le_ext_lo, le_ext_hi, n_ext, plot_svgg);

% -----------------------------------------------------------------------
% <sigma v>_{gamma gamma -> chi chi} at the freeze-out scale T = m_chi/20,
% via detailed balance from svavg('gg').  Mirrors the sigma_chiN figure: one
% value per mechanism at its Omega h^2 = 0.12 coupling, the perturbativity
% (g_q=1) line, and the Delta-Neff band edges remapped from g_q.  Units: cm^3/s.
% Gated by plot_svgg (each non-NaN g_q costs one svavg('gg') call).
% -----------------------------------------------------------------------
svgg_fi      = nan(size(mx_vec)); svgg_fo      = nan(size(mx_vec));
svgg_loop_fi = nan(size(mx_vec)); svgg_loop_fo = nan(size(mx_vec));
svgg_pert    = nan(size(mx_vec));
[vg_chi_lo, vg_chi_hi, vg_phi_lo, vg_phi_hi, vg_tot_lo, vg_tot_hi] = deal(nan(size(mx_vec)));
[vg_pos_lo, vg_pos_hi, vg_neg_lo, vg_neg_hi] = deal(nan(size(mx_vec)));
[vg_tot_lo_c, vg_tot_hi_c, vg_pos_lo_c, vg_pos_hi_c, vg_neg_lo_c, vg_neg_hi_c] = deal(nan(size(mx_vec)));
[vg_all_lo, vg_all_hi, vg_all_lo2, vg_all_hi2, vg_all_lo3, vg_all_hi3] = deal(nan(size(mx_vec)));
[vg_med_lo, vg_med_hi] = deal(nan(size(mx_vec)));
if plot_svgg
    t_svgg = tic;
    parfor j = 1:length(mx_vec)
        m = mx_vec(j);
        f = @(gqv) svgg_gq(m, mphi, gx, gqv, ghat, interaction, hbarc);   %#ok<PFBNS>
        svgg_fi(j)      = f(gqvec_fi(j));
        svgg_fo(j)      = f(gqvec_fo(j));
        svgg_loop_fi(j) = f(gqvec_loop_fi(j));
        svgg_loop_fo(j) = f(gqvec_loop_fo(j));
        svgg_pert(j)    = f(1);
        vg_chi_lo(j)=f(gq_chi_lo(j)); vg_chi_hi(j)=f(gq_chi_hi(j));
        vg_phi_lo(j)=f(gq_phi_lo(j)); vg_phi_hi(j)=f(gq_phi_hi(j));
        vg_tot_lo(j)=f(gq_tot_lo(j)); vg_tot_hi(j)=f(gq_tot_hi(j));
        vg_pos_lo(j)=f(gq_pos_lo(j)); vg_pos_hi(j)=f(gq_pos_hi(j));
        vg_neg_lo(j)=f(gq_neg_lo(j)); vg_neg_hi(j)=f(gq_neg_hi(j));
        vg_tot_lo_c(j)=f(gq_tot_lo_c(j)); vg_tot_hi_c(j)=f(gq_tot_hi_c(j));
        vg_pos_lo_c(j)=f(gq_pos_lo_c(j)); vg_pos_hi_c(j)=f(gq_pos_hi_c(j));
        vg_neg_lo_c(j)=f(gq_neg_lo_c(j)); vg_neg_hi_c(j)=f(gq_neg_hi_c(j));
        vg_all_lo(j)=f(gq_all_lo(j)); vg_all_hi(j)=f(gq_all_hi(j));
        vg_all_lo2(j)=f(gq_all_lo2(j)); vg_all_hi2(j)=f(gq_all_hi2(j));
        vg_all_lo3(j)=f(gq_all_lo3(j)); vg_all_hi3(j)=f(gq_all_hi3(j));
        vg_med_lo(j)=f(gq_med_lo(j)); vg_med_hi(j)=f(gq_med_hi(j));
    end
    fprintf('[compare_limits_g] <sigma v>_gg block: %.2f s\n', toc(t_svgg));
end

% -----------------------------------------------------------------------
% Meson/indirect-detection exclusions in g_q (the nucleon bounds SN1987A/NS/
% 5th-force are already handled above as ext_specs g_q-BANDS).  These five are
% ratio constraints; gq_crit(mx) is the g_q at which the ratio reaches 1 and
% the point is excluded ABOVE it.  (Naively ratio ~ g_q^2, but the mediator
% width spoils that on resonance -- see docs/kaon_constraint.md --
% so gq_crit comes from a root-find seeded with the g_q^2 rescale.)  kaon
% uses kaon_invisible_ratio; kphi/cmb/integral/mesoninv use mesonandCMB_ratio.  Each maps
% to a sigma_chiN / <sigma v> curve at gq_crit.  (kphi is g_q-dependent here --
% unlike the gx-scan, where it is gx-independent and has no curve.)  kphi can
% also have a CEILING (gq_rat_hi): above it the phi decays visibly inside
% NA62 and the bound switches off, so its excluded region is a band.
% -----------------------------------------------------------------------
t_ext = tic;
[gq_rat, sg_rat, vg_rat, gq_rat_hi, sg_rat_hi, vg_rat_hi] = ratio_curves(mx_vec, mphi, ...
    gx, ghat, interaction, interaction_sig, hbarc, rat_show_of(SHOW), plot_svgg);
fprintf('[compare_limits_g] External ratio block: %.2f s\n', toc(t_ext));

% -----------------------------------------------------------------------
% Write data to .dat file
% Columns: mx_GeV  gq_fi  gq_fo  gq_loop_fi  gq_loop_fo
%          sg_fi_cm2  sg_fo_cm2  sg_loop_fi_cm2  sg_loop_fo_cm2  sg_pert_cm2
%          gq_neff_tot_lo  gq_neff_tot_hi  sg_neff_tot_lo_cm2  sg_neff_tot_hi_cm2
%          gq_neff_chi_lo  gq_neff_chi_hi  gq_neff_phi_lo  gq_neff_phi_hi
%          gq_neff_all_lo  gq_neff_all_hi  (the combined BBN-or-CMB band, plotted:
%          up to three excluded g_q intervals by role: *_all_lo/hi = dNeff<0
%          below any positive one, *_all_lo2/hi2 = dNeff>0, *_all_lo3/hi3 = the
%          rest above, usually open above)
%          gq_med_lo  gq_med_hi  sg_med_lo_cm2  sg_med_hi_cm2  (stable mediator
%          overcloses, Omega_A' h^2 > 0.12; freeze_coupled panels only)
% NaN = no solution / no exclusion at that mass point.
% -----------------------------------------------------------------------
if ~exist('dat_tag', 'var'), dat_tag = ''; end   % set by the batch block, if any
dat_fname = sprintf('output/compare_limits_g_%s_mphi%.4g_gx%.4g%s.dat', interaction, mphi, gx, dat_tag);
if ~isfolder('output'), mkdir('output'); end
fid = fopen(dat_fname, 'w');
fprintf(fid, '# compare_limits_g.m: interaction=%s  mphi=%.6g GeV  gx=%.6g  ghat=[%.6g %.6g %.6g]  T_RH=%.6g GeV  neff_bbn=%.3g+/-%.3g  neff_cmb=%.3g+/-%.3g  nsig=%.3g\n', ...
        interaction, mphi, gx, ghat(1), ghat(2), ghat(3), T_RH, neff_mu_bbn, neff_sig_bbn, neff_mu_cmb, neff_sig_cmb, neff_nsig);
fprintf(fid, ['# mx_GeV  gq_fi  gq_fo  gq_loop_fi  gq_loop_fo  sg_fi_cm2  sg_fo_cm2  sg_loop_fi_cm2  sg_loop_fo_cm2  sg_pert_cm2' ...
              '  gq_neff_tot_lo  gq_neff_tot_hi  sg_neff_tot_lo_cm2  sg_neff_tot_hi_cm2  gq_neff_chi_lo  gq_neff_chi_hi  gq_neff_phi_lo  gq_neff_phi_hi' ...
              '  gq_cmb_pos_lo  gq_cmb_pos_hi  gq_cmb_neg_lo  gq_cmb_neg_hi  sg_cmb_pos_lo_cm2  sg_cmb_pos_hi_cm2  sg_cmb_neg_lo_cm2  sg_cmb_neg_hi_cm2' ...
              '  gq_neff_all_lo  gq_neff_all_hi  sg_neff_all_lo_cm2  sg_neff_all_hi_cm2' ...
              '  gq_neff_all_lo2  gq_neff_all_hi2  sg_neff_all_lo2_cm2  sg_neff_all_hi2_cm2' ...
              '  gq_neff_all_lo3  gq_neff_all_hi3  sg_neff_all_lo3_cm2  sg_neff_all_hi3_cm2' ...
              '  gq_med_lo  gq_med_hi  sg_med_lo_cm2  sg_med_hi_cm2\n']);
for jj = 1:length(mx_vec)
    fprintf(fid, ['%.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e' ...
                  '  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e' ...
                  '  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e' ...
                  '  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e' ...
                  '  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n'], ...
        mx_vec(jj), gqvec_fi(jj), gqvec_fo(jj), gqvec_loop_fi(jj), gqvec_loop_fo(jj), ...
        sglim_fi(jj), sglim_fo(jj), sglim_loop_fi(jj), sglim_loop_fo(jj), sglim_pert(jj), ...
        gq_tot_lo(jj), gq_tot_hi(jj), sg_tot_lo(jj), sg_tot_hi(jj), ...
        gq_chi_lo(jj), gq_chi_hi(jj), gq_phi_lo(jj), gq_phi_hi(jj), ...
        gq_pos_lo_c(jj), gq_pos_hi_c(jj), gq_neg_lo_c(jj), gq_neg_hi_c(jj), ...
        sg_pos_lo_c(jj), sg_pos_hi_c(jj), sg_neg_lo_c(jj), sg_neg_hi_c(jj), ...
        gq_all_lo(jj), gq_all_hi(jj), sg_all_lo(jj), sg_all_hi(jj), ...
        gq_all_lo2(jj), gq_all_hi2(jj), sg_all_lo2(jj), sg_all_hi2(jj), ...
        gq_all_lo3(jj), gq_all_hi3(jj), sg_all_lo3(jj), sg_all_hi3(jj), ...
        gq_med_lo(jj), gq_med_hi(jj), sg_med_lo(jj), sg_med_hi(jj));
end
fclose(fid);
fprintf('[compare_limits_g] Wrote %s\n', dat_fname);

% Cache the computed arrays (see the "Computation cache" header above).
CACHE_SKIP = {'SHOW','BATCH','batch_params','REPLOT','REFRESH_EXT','REFRESH_NEFF','RELIC_CACHED','cache_fname','dat_tag', ...
              'dat_fname','fid','ans','plot_sigma_chiN','plot_gx','plot_svgg', ...
              'plot_gq','plot_kaon','plot_cmb','plot_integral','plot_kphi', ...
              'ext_specs','ext_show','nE','ext_excl','s_mix','SIGMA_YCAP','sigma_ycap_run_','FLOAT_ALIGN','CACHE_SKIP','cache_vars','wv_','k_','v_'};
wv_ = who; cache_vars = {};
for k_ = 1:numel(wv_)
    if any(strcmp(wv_{k_}, CACHE_SKIP)), continue; end
    v_ = eval(wv_{k_});
    if isnumeric(v_) || islogical(v_) || ischar(v_) || isstring(v_)
        cache_vars{end+1} = wv_{k_};  %#ok<SAGROW>
    end
end
save(cache_fname, cache_vars{:});
fprintf('[compare_limits_g] Cached %d arrays -> %s\n', numel(cache_vars), cache_fname);
else
% SIGMA_YCAP is a presentation setting of THIS run, but older caches carry the
% value they were built with, and a bare load would silently put it back.  It
% is parked in BATCH, a struct the cache never holds: any plain variable can be
% clobbered by the load (this file is a script, so a leftover from the previous
% panel ends up in the next panel's cache).
BATCH.sigma_ycap_run = SIGMA_YCAP;
load(cache_fname);
SIGMA_YCAP = BATCH.sigma_ycap_run;
fprintf('[compare_limits_g] Replot from %s (no recomputation)\n', cache_fname);
% A cache from before the two-interval dNeff band has no second layer.
if ~exist('gq_all_lo2', 'var')
    [gq_all_lo2, gq_all_hi2, sg_all_lo2, sg_all_hi2, vg_all_lo2, vg_all_hi2] = deal(nan(size(mx_vec)));
end
% ...and one from before the three-role band has no third (it put the rest in
% layer 2); REPLOT = 'neff' rebuilds all three.
if ~exist('gq_all_lo3', 'var')
    [gq_all_lo3, gq_all_hi3, sg_all_lo3, sg_all_hi3, vg_all_lo3, vg_all_hi3] = deal(nan(size(mx_vec)));
end
% ...and one from before 2026-11-17 has no mediator-relic band.
if ~exist('gq_med_lo', 'var')
    [gq_med_lo, gq_med_hi, sg_med_lo, sg_med_hi, vg_med_lo, vg_med_hi] = deal(nan(size(mx_vec)));
end
% Re-derive what the load may have clobbered with stale values: ext_show and nE
% follow the CURRENT SHOW / ext_specs, and a cache written when there were
% fewer external constraints carries shorter versions of both.
nE = size(ext_specs, 1);
ext_show = [SHOW.sn1987a, SHOW.neutronstar, SHOW.stellar, SHOW.fifthforce];
% External bands: a cache written before a constraint was added carries fewer
% rows than ext_specs has, and every ext_*(ie) read below would be out of
% bounds.  REPLOT = 'ext' rebuilds them (and re-saves); a plain replot says so.
if ~exist('le_ext_lo', 'var'), le_ext_lo = -20; end
if ~exist('le_ext_hi', 'var'), le_ext_hi = 8;   end
if ~exist('n_ext', 'var'), n_ext = 1 + round((le_ext_hi - le_ext_lo) / 0.04); end
if REFRESH_EXT
    t_ext = tic;
    [ext_gq_lo, ext_gq_hi, ext_sg_lo, ext_sg_hi, ext_vg_lo, ext_vg_hi] = ...
        ext_bands(ext_specs, mx_vec, mphi, gx, ghat, interaction, interaction_sig, ...
                  hbarc, le_ext_lo, le_ext_hi, n_ext, plot_svgg);
    save(cache_fname, 'ext_gq_lo', 'ext_gq_hi', 'ext_sg_lo', 'ext_sg_hi', ...
         'ext_vg_lo', 'ext_vg_hi', 'le_ext_lo', 'le_ext_hi', 'n_ext', '-append');
    fprintf('[compare_limits_g] External bands recomputed (%.1f s) -> %s\n', ...
            toc(t_ext), cache_fname);
    t_ext = tic;
    [gq_rat, sg_rat, vg_rat, gq_rat_hi, sg_rat_hi, vg_rat_hi] = ratio_curves(mx_vec, ...
        mphi, gx, ghat, interaction, interaction_sig, hbarc, rat_show_of(SHOW), plot_svgg);
    save(cache_fname, 'gq_rat', 'sg_rat', 'vg_rat', 'gq_rat_hi', 'sg_rat_hi', ...
         'vg_rat_hi', '-append');
    fprintf('[compare_limits_g] Ratio constraints recomputed (%.1f s) -> %s\n', ...
            toc(t_ext), cache_fname);
elseif size(gq_rat, 1) ~= numel(rat_show_of(SHOW))
    error('compare_limits_g:staleRat', ...
          ['%s predates the current ratio-constraint set (%d rows, now %d) -- ' ...
           'rerun with replot = ''ext'' to refresh them.'], ...
          cache_fname, size(gq_rat, 1), numel(rat_show_of(SHOW)));
elseif ~exist('ext_gq_lo', 'var') || numel(ext_gq_lo) ~= nE
    error('compare_limits_g:staleExt', ...
          ['%s predates the current external-constraint set (ext_specs now has ' ...
           '%d) -- rerun with replot = ''ext'' to refresh the bands.'], ...
          cache_fname, nE);
end
% A cache from before the kphi ceiling: no ceiling, i.e. excluded to the top,
% as it was drawn then.  replot = 'ext' computes it.
if ~exist('gq_rat_hi', 'var')
    gq_rat_hi = nan(size(gq_rat));  sg_rat_hi = gq_rat_hi;  vg_rat_hi = gq_rat_hi;
end
% The toggles are deliberately not cached, so re-apply what the compute block
% would have derived from the current SHOW.
if ~SHOW.neff, neff_color_mode = 'off'; end
end

% Bullet Cluster: sigma_T/m_chi depends on (mx, mphi, gx) only, so at fixed gx it
% excludes whole m_chi intervals at every g_q.  Closed form, so it is evaluated
% on every run (replots included) and never cached; written to its own sidecar.
[bullet_som, bullet_iv, bullet_bound] = bullet_mx_columns(mx_vec, mphi, gx, interaction);
bullet_fname = sprintf('output/compare_limits_g_%s_mphi%.4g_gx%.4g%s_bullet.dat', ...
                       interaction, mphi, gx, dat_tag);
fid = fopen(bullet_fname, 'w');
fprintf(fid, '# Bullet Cluster: interaction=%s  mphi=%.6g GeV  gx=%.6g  bound=%.3g cm^2/g\n', ...
        interaction, mphi, gx, bullet_bound);
for r_ = 1:size(bullet_iv, 1)
    fprintf(fid, '# excluded mx in [%.6e, %.6e] GeV\n', bullet_iv(r_, 1), bullet_iv(r_, 2));
end
fprintf(fid, '# mx_GeV  sigma_over_m_cm2_per_g  excluded\n');
fprintf(fid, '%.6e  %.6e  %d\n', [mx_vec(:)'; bullet_som(:)'; (bullet_som(:)' > bullet_bound)]);
fclose(fid);
if ~isempty(bullet_iv)
    fprintf('[compare_limits_g] Bullet Cluster excludes %d m_chi interval(s) -> %s\n', ...
            size(bullet_iv, 1), bullet_fname);
end


% -----------------------------------------------------------------------
%% Plot -- shared styling (used by all figures)
% -----------------------------------------------------------------------
use_latex_fonts();          % LaTeX (Computer Modern) font for all text
c_fi   = [0    0.45 0.74];
c_fo   = [0.85 0.33 0.10];
c_loop = [0.20 0.65 0.30];
c_med  = [0.55 0.35 0.10];   % stable-mediator overclosure band
c_bbn = [0.95 0.70 0.10];   % total |Delta-Neff| (magnitude mode)
c_pos = [0.95 0.70 0.10];   % BBN dNeff_total > 0 (dark radiation wins)
c_neg = [0.20 0.45 0.85];   % BBN dNeff_total < 0 (photon heating wins)
c_cmb_pos = [0.80 0.15 0.15];  % CMB dNeff_total > 0
c_cmb_neg = [0.10 0.60 0.70];  % CMB dNeff_total < 0 (photon heating)
c_phi = [0.60 0.40 0.10];   % phi-only edge lines
c_chi = [0.55 0.20 0.65];   % chi-only edge lines
% Meson/indirect-detection overlay colours + labels (kaon/kphi/cmb/integral/
% mesoninv; each excluded ABOVE its curve).  Nucleon externals use ext_specs
% colours above.
ext_rat_cols = {[0.60 0.00 0.30], [0.65 0.46 0.11], [0.09 0.75 0.81], [0.46 0.44 0.70], ...
                [0.05 0.15 0.50]};
ext_rat_lbls = {'$K^+\!\to\pi^+\chi\chi$', sprintf('$K^+\\!\\to\\pi^+%s$', phi_tex(interaction)), 'CMB', 'INTEGRAL', ...
                '$\pi^0,\eta,\eta''\to$ inv.'};
c_bullet = [0.45 0.50 0.05];  % Bullet Cluster m_chi columns

% -----------------------------------------------------------------------
%% Figure 1: sigma_chiN vs m_chi
% -----------------------------------------------------------------------
% Perturbativity gets ONE legend entry, not two: the shaded region owns it
% whenever it is drawn, and the boundary line falls back to owning it only if
% the shading is off.  (Both labelled put the same statement in the legend
% twice; neither labelled left an unexplained dotted line.)
pert_line_args = {'HandleVisibility', 'off'};
if SHOW.pert_line && ~SHOW.nonpert
    pert_line_args = {'DisplayName', 'Perturbativity ($g_q = 1$)'};
end

% Low-mass island walls (see the loop block).  A cache from before 2026-11-18
% has none: find them now with the same test, which costs one solve per
% candidate mass (plus a bisection where a wall exists).  Without the test, a
% run of lower-branch-only masses ahead of the first two-branch mass -- an
% island open-topped there, its upper branch above g_q = 1 -- was taken for a
% bump tip and erased (SxS, mphi = 1: the whole Oh2 = 0.12 line below 0.13 GeV).
if ~exist('relic_tip_mx', 'var') || numel(relic_tip_mx) ~= numel(mx_vec)
    tip_dec_ = 2;  if exist('COALESCE_MAX_DEC', 'var'), tip_dec_ = COALESCE_MAX_DEC; end
    relic_tip_mx = relic_island_tips(mx_vec, gqvec_loop_fi, gqvec_loop_fo, tip_dec_, ...
        relic_oh2_at(solver_name, mphi, gx, ghat, interaction, T_RH, production), Oh2_target);
    clear tip_dec_
end
% The loop branches on a grid with a column at each wall; an inserted column
% carries the couplings of the mass above it.  A lone lower root just below a
% wall (the bump tip, see the loop block; a cache from before 2026-11-18 still
% holds it) is not drawn: its line overshoots the wall and its fill runs to
% g_q = 1.
lo_d_ = gqvec_loop_fi;  sglo_d_ = sglim_loop_fi;
for j = find(isfinite(relic_tip_mx))
    k = j - 1;
    while k >= 1 && ~isnan(lo_d_(k)) && isnan(gqvec_loop_fo(k))
        lo_d_(k) = NaN;  sglo_d_(k) = NaN;  k = k - 1;
    end
end
[mx_rel, ia_rel, ins_rel] = relic_tip_columns(mx_vec, relic_tip_mx);
gq_rel_lo = lo_d_(ia_rel);           gq_rel_hi = gqvec_loop_fo(ia_rel);
sg_rel_lo  = sglo_d_(ia_rel);         sg_rel_hi  = sglim_loop_fo(ia_rel);
clear lo_d_ sglo_d_
sg_rel_pert = sglim_pert(ia_rel);
for k = find(ins_rel)
    sg_rel_lo(k)   = DMnucxsec(mx_rel(k), mphi, gx, gq_rel_lo(k) * ghat, interaction_sig) * hbarc^2;
    sg_rel_hi(k)   = DMnucxsec(mx_rel(k), mphi, gx, gq_rel_hi(k) * ghat, interaction_sig) * hbarc^2;
    sg_rel_pert(k) = DMnucxsec(mx_rel(k), mphi, gx, ghat, interaction_sig) * hbarc^2;
end

fig_sigma = gobjects(0); fig_gq = gobjects(0); fig_svgg = gobjects(0);
if plot_sigma_chiN
fig_sigma = figure; hold on;

mask = ~isnan(sg_rel_lo);
if any(mask) && SHOW.overabundant
    y_lo = sg_rel_lo;
    y_hi = sg_rel_hi;
    use_pert = isnan(y_hi) | (y_hi > sg_rel_pert);
    y_hi(use_pert) = sg_rel_pert(use_pert);
    mx_shade = mx_rel(mask);
    y_lo_m   = y_lo(mask);
    y_hi_m   = y_hi(mask);
    patch([mx_shade, fliplr(mx_shade)], [y_lo_m, fliplr(y_hi_m)], c_loop, ...
          'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off', ...
          'UserData', struct('label_float', true));   % label floats in the fill
    legend_swatch(c_loop, 0.15, 'Overabundant');
end

% Per-branch visibility: [FI, FO, loop-FIMP, loop-WIMP].  The two loop
% branches share one toggle (and one legend entry).
relic_show = [SHOW.relic_fi, SHOW.relic_fo, SHOW.relic_loop, SHOW.relic_loop];
for ii = 1:4
    if ~relic_show(ii), continue; end
    switch ii
        case 1, xv = mx_vec; yvec = sglim_fi;  gqc = gqvec_fi;  col = c_fi;   lbl = 'Freeze-in';           showleg = true;
        case 2, xv = mx_vec; yvec = sglim_fo;  gqc = gqvec_fo;  col = c_fo;   lbl = 'Freeze-out';          showleg = true;
        case 3, xv = mx_rel; yvec = sg_rel_lo; gqc = gq_rel_lo; col = c_loop; lbl = '$\Omega h^2 = 0.12$'; showleg = true;
        case 4, xv = mx_rel; yvec = sg_rel_hi; gqc = gq_rel_hi; col = c_loop; lbl = '';                    showleg = false;
    end
    pert = gqc <= 1;
    ys = yvec; ys(~pert) = NaN;
    yd = yvec; yd( pert) = NaN;
    if showleg
        loglog(xv, ys, '-',  'LineWidth', 2, 'Color', col, 'DisplayName', lbl);
    else
        loglog(xv, ys, '-',  'LineWidth', 2, 'Color', col, 'HandleVisibility', 'off');
    end
    loglog(xv, yd, '--', 'LineWidth', 2, 'Color', col, 'HandleVisibility', 'off');
end
if SHOW.relic_loop
    draw_relic_walls(mx_rel, ins_rel, gq_rel_lo, gq_rel_hi, c_loop, ...
        @(m, e) DMnucxsec(m, mphi, gx, e * ghat, interaction_sig) * hbarc^2);
end

set(gca, 'yscale', 'log', 'xscale', 'log');

yl_for_bbn = ylim; y_bottom = yl_for_bbn(1); y_top_bbn = yl_for_bbn(2);
[y_bottom, y_top_bbn] = row_frame('sigma', y_bottom, y_top_bbn, BATCH);
sg_cap = @(m) DMnucxsec(m, mphi, gx, ceiling_lookup(m, mx_vec, gq_ceiling, 10^le_scan_hi)*ghat, interaction_sig) * hbarc^2;

% Delta-Neff exclusion.  TOGGLE neff_color_mode (set above):
%   'combined'  -> ONE band: excluded by BBN or CMB, either sign (default)
%   'signed'    -> separate colours for dNeff>+max (dark) and <-max (gamma-heating)
%   'magnitude' -> one band per epoch for |dNeff-mu|>max
if strcmp(neff_color_mode, 'off')
    % Delta-Neff drawing disabled (SHOW.neff = false).
elseif strcmp(neff_color_mode, 'combined')
    draw_neff_layers(mx_vec, {sg_all_lo, sg_all_lo2, sg_all_lo3}, ...
        {sg_all_hi, sg_all_hi2, sg_all_hi3}, allon_all, c_bbn, 0.18, ...
        '$\Delta N_{\rm eff}$', y_bottom, y_top_bbn, sg_cap, 'shade');
elseif strcmp(neff_color_mode, 'signed')
    draw_neff_band(mx_vec, sg_pos_lo, sg_pos_hi, allon_pos, c_pos, 0.18, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm BBN}>%.2g$', thr_pos_bbn), y_bottom, y_top_bbn, sg_cap, 'shade');
    draw_neff_band(mx_vec, sg_neg_lo, sg_neg_hi, allon_neg, c_neg, 0.18, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm BBN}<%.2g$', thr_neg_bbn), y_bottom, y_top_bbn, sg_cap, 'shade');
    draw_neff_band(mx_vec, sg_pos_lo_c, sg_pos_hi_c, allon_pos_c, c_cmb_pos, 0.16, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm CMB}>%.2g$', thr_pos_cmb), y_bottom, y_top_bbn, sg_cap, 'shade');
    draw_neff_band(mx_vec, sg_neg_lo_c, sg_neg_hi_c, allon_neg_c, c_cmb_neg, 0.16, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm CMB}<%.2g$', thr_neg_cmb), y_bottom, y_top_bbn, sg_cap, 'shade');
else
    draw_neff_band(mx_vec, sg_tot_lo, sg_tot_hi, allon_tot, c_bbn, 0.18, ...
        sprintf('$|\\Delta N_{\\rm eff}^{\\rm BBN}-\\mu|>%.2g$', neff_nsig*neff_sig_bbn), y_bottom, y_top_bbn, sg_cap, 'shade');
    draw_neff_band(mx_vec, sg_tot_lo_c, sg_tot_hi_c, allon_tot_c, c_cmb_neg, 0.16, ...
        sprintf('$|\\Delta N_{\\rm eff}^{\\rm CMB}-\\mu|>%.2g$', neff_nsig*neff_sig_cmb), y_bottom, y_top_bbn, sg_cap, 'shade');
end

% phi-only / chi-only contributions as edge lines (where each species ALONE
% exceeds the bound).  chi (Dirac) appears only for light m_chi, where it is
% still relativistic at BBN; phi only where it survives to BBN.  These are a
% decomposition of the band above, not extra constraints, so 'combined' drops
% them.
if ~any(strcmp(neff_color_mode, {'combined', 'off'}))
    draw_neff_band(mx_vec, sg_phi_lo, sg_phi_hi, allon_phi, c_phi, 0.0, ...
        sprintf('$\\Delta N_{\\rm eff}\\,(%s)$', phi_tex(interaction)), y_bottom, y_top_bbn, sg_cap, 'edges');
    draw_neff_band(mx_vec, sg_chi_lo, sg_chi_hi, allon_chi, c_chi, 0.0, ...
        '$\Delta N_{\rm eff}\,(\chi)$', y_bottom, y_top_bbn, sg_cap, 'edges');
end

if SHOW.pert_line
    loglog(mx_vec, sglim_pert, ':', 'LineWidth', 2, 'Color', [0.3 0.3 0.3], pert_line_args{:});
end
if SHOW.nonpert
    yl = ylim;
    y_top = row_top('sigma', max(yl(2), max(sglim_pert)*1e3), BATCH);
    patch([mx_vec, fliplr(mx_vec)], [sglim_pert, y_top*ones(size(mx_vec))], ...
          [0.3 0.3 0.3], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    legend_swatch([0.3 0.3 0.3], 0.15, 'Non-perturbative');
    ylim(yl);
end

% Unitarity / EFT-validity ceiling (axialVector): dashed line + shade above.
if any(isfinite(sg_ceiling)) && SHOW.unitarity
    c_eft = [0.49 0.18 0.56];
    loglog(mx_vec, sg_ceiling, '-.', 'LineWidth', 2, 'Color', c_eft, ...
           'HandleVisibility', 'off');
    me = isfinite(sg_ceiling);
    yl_e = ylim; ytop_e = row_top('sigma', max(yl_e(2), max(sg_ceiling(me))*1e3), BATCH);
    patch([mx_vec(me), fliplr(mx_vec(me))], ...
          [sg_ceiling(me), ytop_e*ones(1, sum(me))], c_eft, ...
          'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    legend_swatch(c_eft, 0.12, "Unitarity ($\Gamma_{A'} \ge m_A$ )");
    ylim(yl_e);
end

% Small-mixing-breakdown (axialVector): dashed line + shade above (overlay only).
if any(isfinite(sg_mix)) && SHOW.large_mixing
    c_mix = [0.90 0 0.50];
    loglog(mx_vec, sg_mix, '-.', 'LineWidth', 2, 'Color', c_mix, ...
           'HandleVisibility', 'off');
    mm = isfinite(sg_mix);
    yl_m = ylim; ytop_m = row_top('sigma', max(yl_m(2), max(sg_mix(mm))*1e3), BATCH);
    patch([mx_vec(mm), fliplr(mx_vec(mm))], ...
          [sg_mix(mm), ytop_m*ones(1, sum(mm))], c_mix, ...
          'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    legend_swatch(c_mix, 0.10, sprintf('Large mixing ($\\sin\\theta_{A''}>%.2g$)', s_mix));
    ylim(yl_m);
end

% Stable mediator overcloses (freeze_coupled's Omega_A' h^2 > 0.12).
if SHOW.mediator_relic && any(isfinite(sg_med_lo))
    draw_neff_band(mx_vec, sg_med_lo, sg_med_hi, false(size(mx_vec)), c_med, 0.15, ...
        sprintf('$\\Omega_{%s}h^2 > 0.12$', phi_tex(interaction)), y_bottom, y_top_bbn, sg_cap, 'shade');
end

% External g_q-bands remapped to sigma_chiN (SN1987A / neutron star / axial).
for ie = 1:nE
    if isfinite(ext_gq_lo(ie)) && ext_show(ie)
        draw_neff_band(mx_vec, ext_sg_lo(ie,:), ext_sg_hi(ie,:), false(size(mx_vec)), ...
            ext_specs{ie,3}, 0.12, ext_specs{ie,2}, y_bottom, y_top_bbn, sg_cap, 'shade');
    end
end
% Meson/indirect exclusions (kaon/kphi/cmb/integral) remapped to sigma_chiN.
[eylo, eyhi] = frame_yrange(sglim_fi, sglim_fo, sglim_loop_fi, sglim_loop_fo, ...
    sglim_pert, sg_ceiling, sg_rat);
% The nucleon bands (SN1987A / neutron star / 5th force) ARE drawn just above,
% so the frame has to reach them.  Omitting ext_sg_* here silently CLIPPED
% whatever ran below the floor: on the mphi=0.1 panels the SN1987A band sits
% ~25 decades under this panel's own floor and survived only because share_y
% borrows the mphi=1 partner's -- regroup the row and it vanished with no
% warning.  They lower the FLOOR only; the top stays on frame_yrange's 98th
% percentile, which is what rejects the relic-solver spikes (pouring ~200 deep
% band values into that pool would drag the percentile down with them).  The
% g_q panel has always framed its ext_gq_*; this gives sigma the same deal.
eylo = min(eylo, ext_frame_floor(ext_sg_lo, ext_sg_hi, ...
                                 isfinite(ext_gq_lo(:)') & ext_show));
overlay_rat_externals(mx_vec, sg_rat, sg_rat_hi, ext_rat_cols, ext_rat_lbls, row_top('sigma', eyhi*1e3, BATCH));
apply_ylim('sigma', [eylo/5, eyhi*5], BATCH, SIGMA_YCAP);  % shared-row override wins
if SHOW.bullet, draw_mx_columns(bullet_iv, mx_vec, c_bullet, 0.15, 'Bullet Cluster'); end
xlim([mx_vec(1), mx_vec(end)]);   % truncate to the scanned mx range (no empty pad on the right)

xlabel('$m_\chi$ (GeV)', 'FontSize', FS.label, 'Interpreter', 'latex');
% The y label and numbers are drawn even on a panel that hides them, and only
% emptied afterwards by hide_y_labels -- see there for why.
ylabel('$\sigma_{\chi N}$ (cm$^2$)', 'FontSize', FS.label, 'Interpreter', 'latex');
if SHOW.legend
    legend('show', 'Interpreter', 'latex', 'FontSize', FS.legend, 'Location', 'northwest');
else
    lg = legend(gca); if ~isempty(lg), set(lg, 'Visible', 'off'); end
end
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', FS.tick);
% Full box on all four sides.  'Layer','top' is the load-bearing half: Box is
% on by default, but the default Layer='bottom' draws the axis lines UNDER the
% data, so the filled exclusion regions painted over the top and right spines
% and the box read as open on two sides.
set(gca, 'Box', 'on', 'Layer', 'top');
if ~STYLE_classic
    decade_ticks(gca);   % major tick on every decade
end
corner_caption({sprintf('$%s$', model_tex(interaction)); ...
                sprintf('$%s = %s$ GeV', mphi_tex(interaction), human_tex(mphi)); ...
                sprintf('$%s = %s$', gx_tex(interaction), human_tex(gx))}, FS.legend);
% Inline labels replace the legend.  Placed AFTER corner_caption so the caption
% box is already in the axes and the placement engine can avoid it.
if SHOW.inline_labels
    % The figure MUST already be at its final size: label extent is in points,
    % so the fraction of the axes a label covers -- and therefore every fit and
    % collision test -- depends on it.  save_batch_figure resizes at SAVE time,
    % which is far too late; laying labels out against the on-screen axes and
    % then shrinking the page silently invalidates every placement decision.
    set(gcf, 'Units', 'centimeters');
    fpos_ = get(gcf, 'Position');
    set(gcf, 'Position', [fpos_(1), fpos_(2), BATCH.figsize(1), BATCH.figsize(2)]);
    drawnow;
    if ~STYLE_classic
        tighten_axes(gca);   % reclaim MATLAB's default loose inset
        drawnow;
    end
    % No legend fallback: inline_labels always places a drawn object, taking a
    % crowded spot in preference to leaving a legend box on the panel.
    lbl_h = inline_labels(gca, 'FontSize', FS.inline, 'FloatAlign', FLOAT_ALIGN);
    lbl_h = place_labels_manually(lbl_h, LABEL_POS);
    set(lbl_h, 'Tag', 'inline_label');   % tagged for later lookup
end
if ~SHOW.ylabel, hide_y_labels(gca); end
set(gca,'yscale','log')
set(gca,'xscale','log')
if STYLE_classic, grid on; else, grid off; end
end  % plot_sigma_chiN

% =======================================================================
% Figure 2: g_q vs m_chi (reusing the
% g_q vectors already computed above -- no additional physics runs).
% =======================================================================
if plot_gq
fig_gq = figure; hold on;

% Overabundant band between FIMP and WIMP loop branches, capped at g_q = 1.
mask2 = ~isnan(gq_rel_lo);
if any(mask2) && SHOW.overabundant
    y_lo = gq_rel_lo;
    y_hi = gq_rel_hi;
    use_pert = isnan(y_hi) | (y_hi > 1);
    y_hi(use_pert) = 1;
    mx_shade = mx_rel(mask2);
    patch([mx_shade, fliplr(mx_shade)], [y_lo(mask2), fliplr(y_hi(mask2))], c_loop, ...
          'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off', ...
          'UserData', struct('label_float', true));   % label floats in the fill
    legend_swatch(c_loop, 0.15, 'Overabundant');
end

% Per-branch visibility: [FI, FO, loop-FIMP, loop-WIMP].  The two loop
% branches share one toggle (and one legend entry).
relic_show = [SHOW.relic_fi, SHOW.relic_fo, SHOW.relic_loop, SHOW.relic_loop];
for ii = 1:4
    if ~relic_show(ii), continue; end
    switch ii
        case 1, xv = mx_vec; yvec = gqvec_fi;  col = c_fi;   lbl = 'Freeze-in';           showleg = true;
        case 2, xv = mx_vec; yvec = gqvec_fo;  col = c_fo;   lbl = 'Freeze-out';          showleg = true;
        case 3, xv = mx_rel; yvec = gq_rel_lo; col = c_loop; lbl = '$\Omega h^2 = 0.12$'; showleg = true;
        case 4, xv = mx_rel; yvec = gq_rel_hi; col = c_loop; lbl = '';                    showleg = false;
    end
    pert = yvec <= 1;
    ys = yvec; ys(~pert) = NaN;
    yd = yvec; yd( pert) = NaN;
    if showleg
        loglog(xv, ys, '-',  'LineWidth', 2, 'Color', col, 'DisplayName', lbl);
    else
        loglog(xv, ys, '-',  'LineWidth', 2, 'Color', col, 'HandleVisibility', 'off');
    end
    loglog(xv, yd, '--', 'LineWidth', 2, 'Color', col, 'HandleVisibility', 'off');
end
if SHOW.relic_loop
    draw_relic_walls(mx_rel, ins_rel, gq_rel_lo, gq_rel_hi, c_loop, @(m, e) e);
end

set(gca, 'yscale', 'log', 'xscale', 'log');

% Delta-Neff bands in g_q directly (see the sigma panel for the mode toggle).
yl_for_bbn = ylim; y_bottom = yl_for_bbn(1); y_top_bbn = yl_for_bbn(2);
[y_bottom, y_top_bbn] = row_frame('gq', y_bottom, y_top_bbn, BATCH);
gq_cap = @(m) ceiling_lookup(m, mx_vec, gq_ceiling, 10^le_scan_hi);
if strcmp(neff_color_mode, 'off')
    % Delta-Neff drawing disabled (SHOW.neff = false).
elseif strcmp(neff_color_mode, 'combined')
    draw_neff_layers(mx_vec, {gq_all_lo, gq_all_lo2, gq_all_lo3}, ...
        {gq_all_hi, gq_all_hi2, gq_all_hi3}, allon_all, c_bbn, 0.18, ...
        '$\Delta N_{\rm eff}$', y_bottom, y_top_bbn, gq_cap, 'shade');
elseif strcmp(neff_color_mode, 'signed')
    draw_neff_band(mx_vec, gq_pos_lo, gq_pos_hi, allon_pos, c_pos, 0.18, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm BBN}>%.2g$', thr_pos_bbn), y_bottom, y_top_bbn, gq_cap, 'shade');
    draw_neff_band(mx_vec, gq_neg_lo, gq_neg_hi, allon_neg, c_neg, 0.18, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm BBN}<%.2g$', thr_neg_bbn), y_bottom, y_top_bbn, gq_cap, 'shade');
    draw_neff_band(mx_vec, gq_pos_lo_c, gq_pos_hi_c, allon_pos_c, c_cmb_pos, 0.16, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm CMB}>%.2g$', thr_pos_cmb), y_bottom, y_top_bbn, gq_cap, 'shade');
    draw_neff_band(mx_vec, gq_neg_lo_c, gq_neg_hi_c, allon_neg_c, c_cmb_neg, 0.16, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm CMB}<%.2g$', thr_neg_cmb), y_bottom, y_top_bbn, gq_cap, 'shade');
else
    draw_neff_band(mx_vec, gq_tot_lo, gq_tot_hi, allon_tot, c_bbn, 0.18, ...
        sprintf('$|\\Delta N_{\\rm eff}^{\\rm BBN}-\\mu|>%.2g$', neff_nsig*neff_sig_bbn), y_bottom, y_top_bbn, gq_cap, 'shade');
    draw_neff_band(mx_vec, gq_tot_lo_c, gq_tot_hi_c, allon_tot_c, c_cmb_neg, 0.16, ...
        sprintf('$|\\Delta N_{\\rm eff}^{\\rm CMB}-\\mu|>%.2g$', neff_nsig*neff_sig_cmb), y_bottom, y_top_bbn, gq_cap, 'shade');
end
if ~any(strcmp(neff_color_mode, {'combined', 'off'}))
    draw_neff_band(mx_vec, gq_phi_lo, gq_phi_hi, allon_phi, c_phi, 0.0, ...
        sprintf('$\\Delta N_{\\rm eff}\\,(%s)$', phi_tex(interaction)), y_bottom, y_top_bbn, gq_cap, 'edges');
    draw_neff_band(mx_vec, gq_chi_lo, gq_chi_hi, allon_chi, c_chi, 0.0, ...
        '$\Delta N_{\rm eff}\,(\chi)$', y_bottom, y_top_bbn, gq_cap, 'edges');
end

if SHOW.pert_line
    yline(1, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 2, pert_line_args{:});
end
if SHOW.nonpert
    yl = ylim;
    y_top = row_top('gq', max(yl(2), 1e3), BATCH);
    patch([mx_vec(1), mx_vec(end), mx_vec(end), mx_vec(1)], ...
          [1, 1, y_top, y_top], [0.3 0.3 0.3], ...
          'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    legend_swatch([0.3 0.3 0.3], 0.15, 'Non-perturbative');
    ylim(yl);
end

% Unitarity / EFT-validity ceiling (axialVector): dashed line + shade above.
if any(isfinite(gq_ceiling)) && SHOW.unitarity
    c_eft = [0.49 0.18 0.56];
    loglog(mx_vec, gq_ceiling, '-.', 'LineWidth', 2, 'Color', c_eft, ...
           'HandleVisibility', 'off');
    yl_e = ylim; ytop_e = row_top('gq', max(yl_e(2), 1e3), BATCH);
    me = isfinite(gq_ceiling);
    patch([mx_vec(me), fliplr(mx_vec(me))], ...
          [gq_ceiling(me), ytop_e*ones(1, sum(me))], c_eft, ...
          'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    legend_swatch(c_eft, 0.12, "Unitarity ($\Gamma_{A'} \ge m_A$ )");
    ylim(yl_e);
end

% Small-mixing-breakdown (axialVector): horizontal line at gq_mix + shade above.
if isfinite(gq_mix) && SHOW.large_mixing
    c_mix = [0.90 0 0.50];
    yl_m = ylim; ytop_m = row_top('gq', max(yl_m(2), 1e3), BATCH);
    loglog([mx_vec(1) mx_vec(end)], [gq_mix gq_mix], '-.', ...
           'LineWidth', 2, 'Color', c_mix, 'HandleVisibility', 'off');
    patch([mx_vec(1), mx_vec(end), mx_vec(end), mx_vec(1)], ...
          [gq_mix, gq_mix, ytop_m, ytop_m], c_mix, ...
          'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    legend_swatch(c_mix, 0.10, sprintf('Large mixing ($\\sin\\theta_{A''}>%.2g$)', s_mix));
    ylim(yl_m);
end

% Stable mediator overcloses (freeze_coupled's Omega_A' h^2 > 0.12).
if SHOW.mediator_relic && any(isfinite(gq_med_lo))
    draw_neff_band(mx_vec, gq_med_lo, gq_med_hi, false(size(mx_vec)), c_med, 0.15, ...
        sprintf('$\\Omega_{%s}h^2 > 0.12$', phi_tex(interaction)), y_bottom, y_top_bbn, gq_cap, 'shade');
end

% External constraints as horizontal g_q-bands (mphi & ghat fixed here).
for ie = 1:nE
    if isfinite(ext_gq_lo(ie)) && ext_show(ie)
        yl_sn = ylim;
        patch([mx_vec(1), mx_vec(end), mx_vec(end), mx_vec(1)], ...
              [ext_gq_lo(ie), ext_gq_lo(ie), ext_gq_hi(ie), ext_gq_hi(ie)], ...
              ext_specs{ie,3}, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
              'HandleVisibility', 'off');
        legend_swatch(ext_specs{ie,3}, 0.12, ext_specs{ie,2});
        ylim(yl_sn);
    end
end
% Meson/indirect exclusions (kaon/kphi/cmb/integral) as g_q curves.
[eylo, eyhi] = frame_yrange(gqvec_fi, gqvec_fo, gqvec_loop_fi, gqvec_loop_fo, ...
    gq_ceiling, gq_rat, ext_gq_lo, ext_gq_hi);
overlay_rat_externals(mx_vec, gq_rat, gq_rat_hi, ext_rat_cols, ext_rat_lbls, row_top('gq', eyhi*1e3, BATCH));
apply_ylim('gq', [eylo/5, eyhi*5], BATCH);   % shared-row override wins if present
if SHOW.bullet, draw_mx_columns(bullet_iv, mx_vec, c_bullet, 0.15, 'Bullet Cluster'); end
xlim([mx_vec(1), mx_vec(end)]);   % truncate to the scanned mx range (no empty pad on the right)

xlabel('$m_\chi$ (GeV)', 'FontSize', FS.label, 'Interpreter', 'latex');
% The y label and numbers are drawn even on a panel that hides them, and only
% emptied afterwards by hide_y_labels -- see there for why.
ylabel('$g_q$', 'FontSize', FS.label, 'Interpreter', 'latex');
if SHOW.legend
    legend('show', 'Interpreter', 'latex', 'FontSize', FS.legend, 'Location', 'northwest');
else
    lg = legend(gca); if ~isempty(lg), set(lg, 'Visible', 'off'); end
end
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', FS.tick);
% Full box on all four sides.  'Layer','top' is the load-bearing half: Box is
% on by default, but the default Layer='bottom' draws the axis lines UNDER the
% data, so the filled exclusion regions painted over the top and right spines
% and the box read as open on two sides.
set(gca, 'Box', 'on', 'Layer', 'top');
if ~STYLE_classic
    decade_ticks(gca);   % major tick on every decade
end
corner_caption({sprintf('$%s$', model_tex(interaction)); ...
                sprintf('$%s = %s$ GeV', mphi_tex(interaction), human_tex(mphi)); ...
                sprintf('$%s = %s$', gx_tex(interaction), human_tex(gx))}, FS.legend);
% Inline labels replace the legend.  Placed AFTER corner_caption so the caption
% box is already in the axes and the placement engine can avoid it.
if SHOW.inline_labels
    % The figure MUST already be at its final size: label extent is in points,
    % so the fraction of the axes a label covers -- and therefore every fit and
    % collision test -- depends on it.  save_batch_figure resizes at SAVE time,
    % which is far too late; laying labels out against the on-screen axes and
    % then shrinking the page silently invalidates every placement decision.
    set(gcf, 'Units', 'centimeters');
    fpos_ = get(gcf, 'Position');
    set(gcf, 'Position', [fpos_(1), fpos_(2), BATCH.figsize(1), BATCH.figsize(2)]);
    drawnow;
    if ~STYLE_classic
        tighten_axes(gca);   % reclaim MATLAB's default loose inset
        drawnow;
    end
    % No legend fallback: inline_labels always places a drawn object, taking a
    % crowded spot in preference to leaving a legend box on the panel.
    inline_labels(gca, 'FontSize', FS.inline, 'FloatAlign', FLOAT_ALIGN);
end
if ~SHOW.ylabel, hide_y_labels(gca); end
set(gca,'yscale','log'); set(gca,'xscale','log');
if STYLE_classic, grid on; else, grid off; end
end  % plot_gq

% =======================================================================
% Figure 3: <sigma v>_{gamma gamma -> chi chi} (at T = m_chi/20) vs m_chi.
% Same structure as Figure 1, with sigma_chiN -> <sigma v>_gg.  Curves solid
% where perturbative (g_q <= 1), dashed otherwise; the Delta-Neff band edges
% are the same exclusions remapped from g_q to <sigma v>_gg.
% =======================================================================
if plot_svgg
fig_svgg = figure; hold on;

% Overabundant band: between the FIMP (lo) and WIMP (hi) loop branches,
% capped at the perturbativity line where the WIMP branch is absent/above it.
mask = ~isnan(svgg_loop_fi);
if any(mask) && SHOW.overabundant
    y_lo = svgg_loop_fi;
    y_hi = svgg_loop_fo;
    use_pert = isnan(y_hi) | (y_hi > svgg_pert);
    y_hi(use_pert) = svgg_pert(use_pert);
    mx_shade = mx_vec(mask);
    patch([mx_shade, fliplr(mx_shade)], [y_lo(mask), fliplr(y_hi(mask))], c_loop, ...
          'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off', ...
          'UserData', struct('label_float', true));   % label floats in the fill
    legend_swatch(c_loop, 0.15, 'Overabundant');
end

% Per-branch visibility: [FI, FO, loop-FIMP, loop-WIMP].  The two loop
% branches share one toggle (and one legend entry).
relic_show = [SHOW.relic_fi, SHOW.relic_fo, SHOW.relic_loop, SHOW.relic_loop];
for ii = 1:4
    if ~relic_show(ii), continue; end
    switch ii
        case 1, yvec = svgg_fi;      gqc = gqvec_fi;      col = c_fi;   lbl = 'Freeze-in';        showleg = true;
        case 2, yvec = svgg_fo;      gqc = gqvec_fo;      col = c_fo;   lbl = 'Freeze-out';       showleg = true;
        case 3, yvec = svgg_loop_fi; gqc = gqvec_loop_fi; col = c_loop; lbl = '$\Omega h^2 = 0.12$'; showleg = true;
        case 4, yvec = svgg_loop_fo; gqc = gqvec_loop_fo; col = c_loop; lbl = '';                   showleg = false;
    end
    pert = gqc <= 1;
    ys = yvec; ys(~pert) = NaN;   % solid where perturbative
    yd = yvec; yd( pert) = NaN;   % dashed where non-perturbative
    if showleg
        loglog(mx_vec, ys, '-',  'LineWidth', 2, 'Color', col, 'DisplayName', lbl);
    else
        loglog(mx_vec, ys, '-',  'LineWidth', 2, 'Color', col, 'HandleVisibility', 'off');
    end
    loglog(mx_vec, yd, '--', 'LineWidth', 2, 'Color', col, 'HandleVisibility', 'off');
end

set(gca, 'yscale', 'log', 'xscale', 'log');
yl_for_bbn = ylim; y_bottom = yl_for_bbn(1); y_top_bbn = yl_for_bbn(2);
[y_bottom, y_top_bbn] = row_frame('svgg', y_bottom, y_top_bbn, BATCH);
svgg_cap = @(m) svgg_gq(m, mphi, gx, ceiling_lookup(m, mx_vec, gq_ceiling, 10^le_scan_hi), ghat, interaction, hbarc);
if strcmp(neff_color_mode, 'off')
    % Delta-Neff drawing disabled (SHOW.neff = false).
elseif strcmp(neff_color_mode, 'combined')
    draw_neff_layers(mx_vec, {vg_all_lo, vg_all_lo2, vg_all_lo3}, ...
        {vg_all_hi, vg_all_hi2, vg_all_hi3}, allon_all, c_bbn, 0.18, ...
        '$\Delta N_{\rm eff}$', y_bottom, y_top_bbn, svgg_cap, 'shade');
elseif strcmp(neff_color_mode, 'signed')
    draw_neff_band(mx_vec, vg_pos_lo, vg_pos_hi, allon_pos, c_pos, 0.18, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm BBN}>%.2g$', thr_pos_bbn), y_bottom, y_top_bbn, svgg_cap, 'shade');
    draw_neff_band(mx_vec, vg_neg_lo, vg_neg_hi, allon_neg, c_neg, 0.18, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm BBN}<%.2g$', thr_neg_bbn), y_bottom, y_top_bbn, svgg_cap, 'shade');
    draw_neff_band(mx_vec, vg_pos_lo_c, vg_pos_hi_c, allon_pos_c, c_cmb_pos, 0.16, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm CMB}>%.2g$', thr_pos_cmb), y_bottom, y_top_bbn, svgg_cap, 'shade');
    draw_neff_band(mx_vec, vg_neg_lo_c, vg_neg_hi_c, allon_neg_c, c_cmb_neg, 0.16, ...
        sprintf('$\\Delta N_{\\rm eff}^{\\rm CMB}<%.2g$', thr_neg_cmb), y_bottom, y_top_bbn, svgg_cap, 'shade');
else
    draw_neff_band(mx_vec, vg_tot_lo, vg_tot_hi, allon_tot, c_bbn, 0.18, ...
        sprintf('$|\\Delta N_{\\rm eff}^{\\rm BBN}-\\mu|>%.2g$', neff_nsig*neff_sig_bbn), y_bottom, y_top_bbn, svgg_cap, 'shade');
    draw_neff_band(mx_vec, vg_tot_lo_c, vg_tot_hi_c, allon_tot_c, c_cmb_neg, 0.16, ...
        sprintf('$|\\Delta N_{\\rm eff}^{\\rm CMB}-\\mu|>%.2g$', neff_nsig*neff_sig_cmb), y_bottom, y_top_bbn, svgg_cap, 'shade');
end
if ~any(strcmp(neff_color_mode, {'combined', 'off'}))
    draw_neff_band(mx_vec, vg_phi_lo, vg_phi_hi, allon_phi, c_phi, 0.0, ...
        sprintf('$\\Delta N_{\\rm eff}\\,(%s)$', phi_tex(interaction)), y_bottom, y_top_bbn, svgg_cap, 'edges');
    draw_neff_band(mx_vec, vg_chi_lo, vg_chi_hi, allon_chi, c_chi, 0.0, ...
        '$\Delta N_{\rm eff}\,(\chi)$', y_bottom, y_top_bbn, svgg_cap, 'edges');
end

if SHOW.pert_line
    loglog(mx_vec, svgg_pert, ':', 'LineWidth', 2, 'Color', [0.3 0.3 0.3], pert_line_args{:});
end

% Shade above the perturbativity line: g_q > 1 required -> ruled out.
ppmask = ~isnan(svgg_pert);
if any(ppmask) && SHOW.nonpert
    yl = ylim;
    y_top = row_top('svgg', max(yl(2), max(svgg_pert(ppmask))*1e3), BATCH);
    mxp = mx_vec(ppmask); svp = svgg_pert(ppmask);
    patch([mxp, fliplr(mxp)], [svp, y_top*ones(size(svp))], ...
          [0.3 0.3 0.3], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    legend_swatch([0.3 0.3 0.3], 0.15, 'Non-perturbative');
    ylim(yl);
end

% Stable mediator overcloses (freeze_coupled's Omega_A' h^2 > 0.12).
if plot_svgg && SHOW.mediator_relic && any(isfinite(vg_med_lo))
    draw_neff_band(mx_vec, vg_med_lo, vg_med_hi, false(size(mx_vec)), c_med, 0.15, ...
        sprintf('$\\Omega_{%s}h^2 > 0.12$', phi_tex(interaction)), y_bottom, y_top_bbn, svgg_cap, 'shade');
end

% External g_q-bands remapped to <sigma v>_gg (SN1987A / neutron star / axial).
if plot_svgg
    for ie = 1:nE
        if isfinite(ext_gq_lo(ie)) && ext_show(ie)
            draw_neff_band(mx_vec, ext_vg_lo(ie,:), ext_vg_hi(ie,:), false(size(mx_vec)), ...
                ext_specs{ie,3}, 0.12, ext_specs{ie,2}, y_bottom, y_top_bbn, svgg_cap, 'shade');
        end
    end
    % Meson/indirect exclusions (kaon/kphi/cmb/integral) remapped to <sigma v>_gg.
    [eylo, eyhi] = frame_yrange(svgg_fi, svgg_fo, svgg_loop_fi, svgg_loop_fo, ...
        svgg_pert, vg_rat);
    % Same clipping fix as the sigma panel above -- floor only.
    eylo = min(eylo, ext_frame_floor(ext_vg_lo, ext_vg_hi, ...
                                     isfinite(ext_gq_lo(:)') & ext_show));
    overlay_rat_externals(mx_vec, vg_rat, vg_rat_hi, ext_rat_cols, ext_rat_lbls, row_top('svgg', eyhi*1e3, BATCH));
    apply_ylim('svgg', [eylo/5, eyhi*5], BATCH);   % shared-row override wins if present
    if SHOW.bullet, draw_mx_columns(bullet_iv, mx_vec, c_bullet, 0.15, 'Bullet Cluster'); end
xlim([mx_vec(1), mx_vec(end)]);   % truncate to the scanned mx range (no empty pad on the right)
end

xlabel('$m_\chi$ (GeV)', 'FontSize', FS.label, 'Interpreter', 'latex');
% The y label and numbers are drawn even on a panel that hides them, and only
% emptied afterwards by hide_y_labels -- see there for why.
ylabel('$\langle \sigma v \rangle_{\gamma\gamma \to \chi\chi}$ (cm$^3$/s) at $T = m_\chi/20$', ...
       'FontSize', FS.label, 'Interpreter', 'latex');
if SHOW.legend
    legend('show', 'Interpreter', 'latex', 'FontSize', FS.legend, 'Location', 'northwest');
else
    lg = legend(gca); if ~isempty(lg), set(lg, 'Visible', 'off'); end
end
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', FS.tick);
% Full box on all four sides.  'Layer','top' is the load-bearing half: Box is
% on by default, but the default Layer='bottom' draws the axis lines UNDER the
% data, so the filled exclusion regions painted over the top and right spines
% and the box read as open on two sides.
set(gca, 'Box', 'on', 'Layer', 'top');
if ~STYLE_classic
    decade_ticks(gca);   % major tick on every decade
end
corner_caption({sprintf('$%s$', model_tex(interaction)); ...
                sprintf('$%s = %s$ GeV', mphi_tex(interaction), human_tex(mphi)); ...
                sprintf('$%s = %s$', gx_tex(interaction), human_tex(gx))}, FS.legend);
% Inline labels replace the legend.  Placed AFTER corner_caption so the caption
% box is already in the axes and the placement engine can avoid it.
if SHOW.inline_labels
    % The figure MUST already be at its final size: label extent is in points,
    % so the fraction of the axes a label covers -- and therefore every fit and
    % collision test -- depends on it.  save_batch_figure resizes at SAVE time,
    % which is far too late; laying labels out against the on-screen axes and
    % then shrinking the page silently invalidates every placement decision.
    set(gcf, 'Units', 'centimeters');
    fpos_ = get(gcf, 'Position');
    set(gcf, 'Position', [fpos_(1), fpos_(2), BATCH.figsize(1), BATCH.figsize(2)]);
    drawnow;
    if ~STYLE_classic
        tighten_axes(gca);   % reclaim MATLAB's default loose inset
        drawnow;
    end
    % No legend fallback: inline_labels always places a drawn object, taking a
    % crowded spot in preference to leaving a legend box on the panel.
    inline_labels(gca, 'FontSize', FS.inline, 'FloatAlign', FLOAT_ALIGN);
end
if ~SHOW.ylabel, hide_y_labels(gca); end
set(gca,'yscale','log'); set(gca,'xscale','log');
if STYLE_classic, grid on; else, grid off; end
end  % plot_svgg

% -----------------------------------------------------------------------
% Batch figure saving: only when a driver asked for it (batch_params.figdir).
% Each figure produced this run is written as <figdir>/<figtag>_<panel>.<ext>,
% and closed unless close_figs was set false, so a many-model sweep does not
% accumulate open windows.  Never runs on an interactive invocation.
% -----------------------------------------------------------------------
if BATCH.save
    save_batch_figure(fig_sigma, BATCH, 'sigma');
    save_batch_figure(fig_gq,   BATCH, 'gq');
    save_batch_figure(fig_svgg,  BATCH, 'svgg');
end
toc

% =======================================================================
% Local helpers.
% =======================================================================

function save_batch_figure(fh, B, panel)
% Write figure fh to <B.figdir>/<B.figtag>_<panel>.<ext> for every ext in
% B.figfmt, then close it unless B.close is false.  No-op for an empty/deleted
% handle (that panel was toggled off).
%
% Sizing: the figure is resized to B.figsize (cm) first.  Without it the page is
% whatever the on-screen window happened to be, which on a Retina display is
% ~3.4 x 2.7 in -- half the intended size, and far too small for the 16/14/13 pt
% fonts.  PDFs then go through printpdf, which sets PaperSize from the figure so
% the page is exactly B.figsize; exportgraphics is the fallback when printpdf is
% not on the path and tight-crops instead.  '.fig' goes through savefig, other
% raster formats through exportgraphics at 300 dpi.
if isempty(fh) || ~isgraphics(fh), return; end
if ~isfolder(B.figdir), mkdir(B.figdir); end
if isfield(B, 'figsize') && numel(B.figsize) == 2
    set(fh, 'Units', 'centimeters');
    pos = get(fh, 'Position');
    set(fh, 'Position', [pos(1), pos(2), B.figsize(1), B.figsize(2)]);
    drawnow;                     % let the axes relayout before the paper size is read
end
% PDF pages are cut to an explicit box, not to whatever ink happens to land on
% the page.  Two reasons:
%
%   * exportgraphics tight-crops to the ink, so a panel with less ink above its
%     axes came out SHORTER, and panels meant to sit side by side in one row did
%     not line up.  panel_layout pins the vertical geometry; this keeps the page
%     that pins it, by exporting the FULL figure height every time.
%   * which exporter ran used to depend on whether printpdf happened to be on
%     the MATLAB path -- printpdf writes the full page, exportgraphics crops --
%     so the same code produced two different geometries on two machines, and
%     whichever ran last silently changed how the figures had to be placed.
%     No path sniffing now.
%
% Horizontally the box still drops the y-label gutter from a panel whose labels
% hide_y_labels emptied: that gutter is blank, and cropping it is what lets the
% shared pair butt together in the document.
enforce_panel_layout(fh);
[bx, by, bw, bh] = panel_export_box(fh);
for k = 1:numel(B.figfmt)
    ext   = B.figfmt{k};
    fname = fullfile(B.figdir, sprintf('%s_%s.%s', B.figtag, panel, ext));
    try
        if strcmpi(ext, 'fig')
            savefig(fh, fname);
        elseif strcmpi(ext, 'pdf')
            set(fh, 'PaperUnits', 'centimeters', 'PaperSize', [bw bh], ...
                    'PaperPositionMode', 'manual', ...
                    'PaperPosition', [-bx, -by, B.figsize(1), B.figsize(2)]);
            print(fh, '-dpdf', '-vector', fname);
        else
            exportgraphics(fh, fname, 'Resolution', 300);
        end
        fprintf('[compare_limits_g] Wrote %s\n', fname);
    catch err
        fprintf(2, '[compare_limits_g] save %s failed: %s\n', fname, err.message);
    end
end
if B.close
    close(fh);
else
    % Kept open (close_figs = false); the appdata lets a caller re-save it
    % through this same function with the identical page box.
    setappdata(fh, 'label_edit', struct('B', B, 'panel', panel, 'save', @save_batch_figure));
end
end

function enforce_panel_layout(fh)
% Last word on the exported axes geometry, applied to the figure about to be
% printed.  tighten_axes already does this, but it is reached only through the
% SHOW.inline_labels branch, and at least two panels (lightScalarPP and
% lightScalarSP at mphi=1e-9) come out of the plotting block without it -- their
% axes was never the gca that branch tightened.  They kept
% PositionConstraint='outerposition', under which MATLAB RE-DERIVES the geometry
% at print time against the paper size set below, so what was measured on screen
% is not what lands in the PDF (Position read 0.7059 where the exported frame
% was 0.7640).
L  = panel_layout();
ax = findobj(fh, 'Type', 'axes');
for k = 1:numel(ax)
    % The whole box comes from panel_layout, never from the current TightInset:
    % after hide_y_labels that inset no longer includes the y gutter, and
    % re-measuring it would hand the unlabelled panel a WIDER plot box than the
    % labelled one it shares an axis with -- the regression 127099a fixed.
    p = [L.ti_left, L.ti_bottom, 1 - L.ti_left - L.ti_right, 1 - L.ti_bottom - L.ti_top];
    if p(3) > 0 && p(4) > 0
        set(ax(k), 'PositionConstraint', 'innerposition');
        set(ax(k), 'Position', p);
    end
end
end

function [x0, y0, w, h] = panel_export_box(fh)
% The rectangle of the figure that becomes the PDF page, in centimetres from the
% figure's bottom-left.  Full height always (panel_layout pins what sits where
% vertically, so every panel exports the same height); full width for a panel
% that still has its y label, and everything right of the now-blank gutter for
% one whose labels hide_y_labels emptied.
PAD = 0.02;                       % cm of air outside the leftmost ink
set(fh, 'Units', 'centimeters');
fp = get(fh, 'Position');
W  = fp(3); H = fp(4);
x0 = 0; y0 = 0; w = W; h = H;
ax = findobj(fh, 'Type', 'axes');
if isempty(ax), return; end
ax = ax(1);
ylab = get(get(ax, 'YLabel'), 'String');
ytl  = get(ax, 'YTickLabel');
if ~iscell(ytl), ytl = cellstr(ytl); end
hidden = isempty(ylab) && (isempty(ytl) || all(cellfun(@isempty, ytl)));
if hidden
    % Cut just left of the leftmost INK, not at the axes edge.  Once the y
    % labels are emptied the axes' own TightInset measures exactly what still
    % pokes out to the left -- the half of the first x tick label that overhangs
    % the axes.  Cutting at the axes edge slices that label in half ("10^-4"
    % came out as "-4").
    p  = get(ax, 'Position');
    ti = get(ax, 'TightInset');   % re-measured AFTER hide_y_labels
    x0 = max(0, (p(1) - ti(1))*W - PAD);
    w  = W - x0;
end
end

function overlay_rat_externals(x, Y, Yhi, cols, lbls, ytop_fill)
% Draw the meson/indirect ratio externals (rows of Y) on the current axes: each
% as a line at Y(ie,:) + translucent shade ABOVE it (the excluded side), filled
% up to the ceiling YHI(ie,:) where it is finite (kphi: the phi decays visibly
% above it; also drawn dash-dot) and to YTOP_FILL elsewhere (chosen above the
% figure's final frame top so the excluded region reaches the top edge).  NaN
% rows (constraint closed / not applicable) are skipped.  The caller sets the
% final y-limits (see frame_yrange).
%
% Each contiguous finite run of Y is its own region: a NaN gap is a real hole,
% not something to bridge.  Where a run stops inside the scan (kaon past NA62's
% m_miss^2 acceptance, INTEGRAL at the ends of its energy window), a vertical
% dash-dot edge joins floor to ceiling (or to the frame top), so the region
% reads as closed.  No edge at the ends of the scan, which are the frame.  As
% compare_limits_g_mphi.m.
for ie = 1:size(Y, 1)
    y = Y(ie, :);
    fin = isfinite(y);
    if ~any(fin), continue; end
    d = diff([false, fin, false]);
    run_lo = find(d == 1);  run_hi = find(d == -1) - 1;
    for r = 1:numel(run_lo)
        k = run_lo(r):run_hi(r);
        yt = ytop_fill*ones(1, numel(k));
        yh = Yhi(ie, k);
        yt(isfinite(yh)) = yh(isfinite(yh));
        % The legend row is the SHADE, not the boundary line: a filled patch
        % means "the model cannot lie here", and every exclusion has to read the
        % same way (a line is reserved for the relic loci).  The dash-dot
        % boundary is still drawn, it just does not claim the legend entry.
        loglog(x(k), y(k), '-.', 'LineWidth', 2, 'Color', cols{ie}, ...
               'HandleVisibility', 'off');
        if any(isfinite(yh))
            loglog(x(k), yh, '-.', 'LineWidth', 2, 'Color', cols{ie}, ...
                   'HandleVisibility', 'off');
        end
        for e = unique([k(1), k(end)])
            if e == 1 || e == numel(x), continue; end   % scan edge = frame edge
            loglog([x(e) x(e)], [y(e) yt(k == e)], '-.', 'LineWidth', 2, ...
                   'Color', cols{ie}, 'HandleVisibility', 'off');
        end
        patch([x(k), fliplr(x(k))], [y(k), fliplr(yt)], cols{ie}, ...
              'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    legend_swatch(cols{ie}, 0.15, lbls{ie}, '-.');
end
end

function lo = ext_frame_floor(v_lo, v_hi, drawn)
% Lowest positive finite value of the external nucleon bands that are actually
% DRAWN -- same condition as the draw loop (a finite g_q band AND its SHOW
% toggle), so hiding a band also stops it stretching the axis.  Returns Inf
% when nothing is drawn, which makes the min() at the call site a no-op.
lo = Inf;
if ~any(drawn), return; end
v = [v_lo(drawn, :), v_hi(drawn, :)];
v = v(isfinite(v) & v > 0);
if ~isempty(v), lo = min(v(:)); end
end

function [lo, hi] = frame_yrange(varargin)
% Robust y-range over all finite POSITIVE values of the supplied arrays.  The
% low edge is the true minimum (frame reaches the bottom of the plotted lines --
% in the g_q scan the freeze-in curves dive to very small sigma/<sv> as g_q->0).
% The high edge is the 98th percentile, which rejects the spurious relic-solver
% spikes (occasional huge g_q roots -> g_q^2-scaled sigma/<sv> outliers) that
% would otherwise blow the top of the log-axis up by ~20 decades.
v = [];
for k = 1:numel(varargin)
    a = varargin{k};
    v = [v, a(:).'];   %#ok<AGROW>
end
v = sort(v(isfinite(v) & v > 0));
if isempty(v), lo = 1e-30;  hi = 1;  return;  end
lo = v(1);
hi = v(min(numel(v), round(0.98*numel(v))));
if lo >= hi, hi = v(end);  end
end

function r = nd_cached(cache, mphi, mx, gx, g, interaction, T_RH)
% Memoized Delta-Neff (neff_simple.m), returned as r.bbn / r.cmb, each the row
% [dNeff_chi, dNeff_phi, dNeff_total] consumed by neff_band_edges.  neff_simple
% is CMB-only and does not decompose, so r.bbn and the chi/phi entries are NaN.
% Keyed on the (rounded) log10 of the q-coupling norm, which tracks g_q, so the
% band sweeps share one evaluation per point.  A failure gives NaN (no band)
% rather than killing the sweep.
key = round(log10(norm(g) + realmin) * 1e9);
if isKey(cache, key)
    r = cache(key);
    return;
end
r.bbn = [NaN NaN NaN];
try
    dN    = neff_simple(mphi, mx, gx, g, interaction, struct('T_RH', T_RH));
    r.cmb = [NaN, NaN, dN];
catch
    r.cmb = [NaN NaN NaN];
end
cache(key) = r;   %#ok<NASGU> handle object, mutated in place
end

function draw_neff_layers(mx, los, his, allon, col, alpha, label, yb, yt, cap, style)
% The combined dNeff band as up to three interval layers (see the band block).
% allon (the whole scan excluded) rides on layer 1.  The legend row / inline
% label goes on the first layer that has any geometry, so a panel whose only
% band is the open-above layer 3 still gets one.
done = false;
for k = 1:numel(los)
    al = false(size(mx));
    if k == 1, al = allon; end
    lab = '';
    if ~done && (any(isfinite(los{k})) || any(al))
        lab = label; done = true;
    end
    draw_neff_band(mx, los{k}, his{k}, al, col, alpha, lab, yb, yt, cap, style);
end
end

function z = neff_pull(r, mu_b, sig_b, mu_c, sig_c, use_b, use_c, signed)
% Worst-epoch pull of the SIGNED total Delta-Neff: max(|dNeff-mu|/sigma) over
% {BBN, CMB}.  z > nsig is the union of the per-epoch, per-sign exclusions,
% which is what the 'combined' band draws.  A NaN from either epoch propagates
% (two-input max would silently drop it and report the other epoch as the worst).
% signed = true returns the worst pull WITH its sign (|z| unchanged), which
% neff_band_edges uses to find allowed strips between a + and a - interval.
if nargin < 6, use_b = true; end
if nargin < 7, use_c = true; end
if nargin < 8, signed = false; end
z = [];
if use_b, z(end+1) = (r.bbn(3) - mu_b) / sig_b; end
if use_c, z(end+1) = (r.cmb(3) - mu_c) / sig_c; end
if isempty(z) || any(isnan(z)), z = NaN; return; end
[~, iw] = max(abs(z));
z = z(iw);
if ~signed, z = abs(z); end
end

function E = neff_edges_off()
% Shape neff_band_edges returns, for an epoch that has been switched off:
% every edge NaN, every allon false, so the drawing code leaves it out.
nan2 = [NaN NaN];
E = struct('chi',nan2, 'phi',nan2, 'tot',nan2, 'tot_pos',nan2, 'tot_neg',nan2, ...
           'allon', struct('chi',false,'phi',false,'tot',false, ...
                           'tot_pos',false,'tot_neg',false));
end


function [ext_gq_lo, ext_gq_hi, ext_sg_lo, ext_sg_hi, ext_vg_lo, ext_vg_hi] = ...
    ext_bands(ext_specs, mx_vec, mphi, gx, ghat, interaction, interaction_sig, ...
              hbarc, le_ext_lo, le_ext_hi, n_ext, plot_svgg)
% External nucleon constraints as g_q-BANDS -- see the call site for the
% reasoning.  A function so that a replot can refresh JUST these bands
% (REPLOT = 'ext') without re-solving anything else.
nE = size(ext_specs, 1);
ext_opts = struct('enable_bullet', false, 'enable_kaon', false, ...
                  'enable_kphi', false, 'enable_cmb', false, ...
                  'enable_integral', false, 'enable_mesoninv', false);
le_sn = linspace(le_ext_lo, le_ext_hi, n_ext);
excl_grid = false(numel(le_sn), nE);
for i = 1:numel(le_sn)
    Li = external_loglikes(mx_vec(1), mphi, gx, 10^le_sn(i)*ghat, interaction, ext_opts);
    for ie = 1:nE
        fi = Li.info.(ext_specs{ie, 1});
        excl_grid(i, ie) = isfield(fi, 'excluded') && fi.excluded;
    end
end
ext_gq_lo = nan(nE, 1);  ext_gq_hi = nan(nE, 1);
for ie = 1:nE
    if any(excl_grid(:, ie))
        ext_gq_lo(ie) = 10^min(le_sn(excl_grid(:, ie)));
        ext_gq_hi(ie) = 10^max(le_sn(excl_grid(:, ie)));
        warning('compare_limits_g:ext', '%s excludes gq in [%.3g, %.3g].', ...
            ext_specs{ie, 2}, ext_gq_lo(ie), ext_gq_hi(ie));
    end
end
% Remap each g_q-band onto sigma_chiN (always) and <sigma v>_gg (if plotted),
% mx-dependent because sigma/<sigma v> depend on mx at fixed g_q.
ext_sg_lo = nan(nE, numel(mx_vec));  ext_sg_hi = nan(nE, numel(mx_vec));
ext_vg_lo = nan(nE, numel(mx_vec));  ext_vg_hi = nan(nE, numel(mx_vec));
for ie = 1:nE
    if ~isfinite(ext_gq_lo(ie)), continue; end
    for j = 1:length(mx_vec)
        ext_sg_lo(ie, j) = DMnucxsec(mx_vec(j), mphi, gx, ext_gq_lo(ie)*ghat, interaction_sig) * hbarc^2;
        ext_sg_hi(ie, j) = DMnucxsec(mx_vec(j), mphi, gx, ext_gq_hi(ie)*ghat, interaction_sig) * hbarc^2;
    end
    if plot_svgg
        for j = 1:length(mx_vec)
            ext_vg_lo(ie, j) = svgg_gq(mx_vec(j), mphi, gx, ext_gq_lo(ie), ghat, interaction, hbarc);
            ext_vg_hi(ie, j) = svgg_gq(mx_vec(j), mphi, gx, ext_gq_hi(ie), ghat, interaction, hbarc);
        end
    end
end
end

function sv = svgg_gq(m, mphi, gx, gqv, ghat, interaction, hbarc)
% g_q scan wrapper for svgg_at_fo: gx is fixed, g_q scales ghat.  NaN g_q
% (no relic/edge solution) -> NaN out, so the parfor isn't fed a NaN g-vector.
if isnan(gqv), sv = NaN; return; end
sv = svgg_at_fo(m, mphi, gx, gqv*ghat, interaction, hbarc);
end

function sv = svgg_at_fo(m, mphi, gxv, g, interaction, hbarc)
% <sigma v>_{gamma gamma -> chi chi} at the freeze-out scale T = m/20,
% obtained from svavg('gg') (chi chi -> gamma gamma) by detailed balance.
% Returns cm^3/s; NaN in -> NaN out.  Uses the SAME chi equilibrium-density
% convention as freeze_loop's Y_eq (g_chi = 4):
%   n_chi_eq / n_gamma_eq = x^2 K2(x) / zeta(3),   x = m/T = 20.
if isnan(gxv), sv = NaN; return; end
T     = m / 20;                                  % freeze-out scale
x     = m / T;                                   % = 20
zeta3 = 1.2020569031595942;
try
    sv_xx = svavg(m, mphi, gxv, g, T, interaction, 'gg');   % GeV^-2, chi chi -> gg
catch err
    fprintf(2, '[svgg] mx=%.4g: %s\n', m, err.message);
    sv = NaN; return;
end
ratio = x^2 * besselk(2, x) / zeta3;             % n_chi_eq / n_gamma_eq
sv_gg = sv_xx * ratio^2;                         % GeV^-2, gg -> chi chi
% GeV^-2 -> cm^3/s:  sigma[cm^2] = sigma[GeV^-2]*hbarc^2; times c = 2.998e10 cm/s.
sv    = sv_gg * hbarc^2 * 2.998e10;
end

function gxc = ceiling_lookup(m, mx_vec, ceil_vec, fallback)
% EFT-ceiling value at mass m by log-log interpolation of ceil_vec(mx_vec).
% Falls back (no cap) where the ceiling is NaN/Inf -- non-axialVector (all NaN)
% or m_chi above all meson->chichi thresholds -- returning the supplied value.
fin = isfinite(ceil_vec) & ceil_vec > 0;
if ~any(fin), gxc = fallback; return; end
lg = interp1(log10(mx_vec(fin)), log10(ceil_vec(fin)), log10(m), 'linear', NaN);
if isnan(lg), gxc = fallback; else, gxc = 10^lg; end
end

function draw_neff_band(mx, y_lo, y_hi, allon, col, alpha, lbl, y_bottom, y_top, cap_fn, mode)
% Draw a Delta-Neff exclusion band given its lower/upper edges (y_lo, y_hi)
% in the current axes' y-units (sigma_chiN or g_q).  mode 'shade' fills the
% band (NaN upper edge -> cap_fn(m)); mode 'edges' draws the two edge lines
% only.  allon marks mass points excluded across the whole coupling scan
% (filled full-height in 'shade' mode).
shade   = strcmp(mode, 'shade');
labeled = false;
if shade
    loglog(mx, y_lo, '-.', 'LineWidth', 1.5, 'Color', col, 'HandleVisibility', 'off');
    loglog(mx, y_hi, '-.', 'LineWidth', 1.5, 'Color', col, 'HandleVisibility', 'off');
    % The band exists wherever EITHER edge is defined.  A missing lower edge
    % means "open below" (fill from y_bottom); a missing upper edge means
    % "open above" (fill up to cap_fn).  Keying the fill on the lower edge alone
    % left open-below bands (e.g. the CMB band, excluded at small coupling) with
    % an orphaned dash-dot upper edge and no legend entry.
    mask = ~isnan(y_lo) | ~isnan(y_hi);
    if any(mask)
        % ONE PATCH PER CONTIGUOUS RUN, not one polygon over every masked point.
        % A single polygon built from mx(mask) BRIDGES any gap in the band,
        % joining two disjoint excluded regions across a mass range that is not
        % excluded at all.  Splitting also gives each run a well-defined pair of
        % ends, which is what the closing edges below need.
        d    = diff([false, mask(:)', false]);
        runs = [find(d == 1); find(d == -1) - 1];        % [start; stop] per run
        for r = 1:size(runs, 2)
            idx = runs(1,r):runs(2,r);
            mb  = mx(idx);  ylo = y_lo(idx);  yhi = y_hi(idx);
            ylo(isnan(ylo)) = y_bottom;
            for kk = find(isnan(yhi(:))'), yhi(kk) = cap_fn(mb(kk)); end
            yhi = max(yhi, ylo);
            patch([mb, fliplr(mb)], [ylo, fliplr(yhi)], col, ...
                  'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            % CLOSE THE ENDS.  The two horizontal edges are drawn as dash-dot
            % lines above, but a band that stops partway across the mass range
            % had no drawn boundary there -- the fill simply ended in mid-air.
            % Draw the vertical only where the run genuinely terminates INSIDE
            % the plotted range; a run reaching the first or last mass point is
            % open, not closed, and a vertical there would be a fiction.
            if idx(1) > 1
                loglog(mb([1 1]), [ylo(1) yhi(1)], '-.', 'LineWidth', 1.5, ...
                       'Color', col, 'HandleVisibility', 'off');
            end
            if idx(end) < numel(mx)
                loglog(mb([end end]), [ylo(end) yhi(end)], '-.', 'LineWidth', 1.5, ...
                       'Color', col, 'HandleVisibility', 'off');
            end
        end
        if ~isempty(lbl), legend_swatch(col, alpha, lbl); end
        labeled = true;
    end
else
    if any(~isnan(y_lo))
        loglog(mx, y_lo, '-.', 'LineWidth', 1.8, 'Color', col, 'DisplayName', lbl);
        labeled = true;
    end
    if any(~isnan(y_hi))
        if labeled
            loglog(mx, y_hi, '-.', 'LineWidth', 1.8, 'Color', col, 'HandleVisibility', 'off');
        else
            loglog(mx, y_hi, '-.', 'LineWidth', 1.8, 'Color', col, 'DisplayName', lbl);
            labeled = true;
        end
    end
end
if shade && any(allon)
    e0 = diff([0, allon(:)', 0]);
    starts = find(e0 == 1); stops = find(e0 == -1) - 1;
    for k = 1:length(starts)
        idx = starts(k):stops(k);
        xr  = mx(idx);
        % Top of the full-height column follows cap_fn (the EFT ceiling for
        % axialVector; ~y_top elsewhere, where cap_fn is huge) so the column
        % does not claim exclusion above the ceiling.
        topr = arrayfun(@(mm) min(cap_fn(mm), y_top), xr);
        topr = max(topr, y_bottom);
        if numel(idx) == 1   % single-mass run: give the patch a finite width
            xr = xr * [1/1.01, 1.01]; topr = [topr, topr];
        end
        px = [xr, fliplr(xr)];
        py = [y_bottom * ones(size(xr)), fliplr(topr)];
        if ~labeled && k == 1
            patch(px, py, col, 'FaceAlpha', alpha, 'EdgeColor', 'none', ...
                  'HandleVisibility', 'off');
            legend_swatch(col, alpha, lbl);
            labeled = true;
        else
            patch(px, py, col, 'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        % BOUNDARY.  Every other exclusion in the panel is bounded by a dash-dot
        % line (the legend convention), but the full-height column was drawn
        % edgeless, so a region made ENTIRELY of allon columns -- which is what
        % the Delta-Neff band is once the answer stops depending on the coupling
        % -- arrived with no boundary anywhere.  Same rule the band fill uses:
        % draw an edge only where the region genuinely ends INSIDE the plotted
        % range.  A run reaching the first or last mass is open, and the top is
        % a boundary only where cap_fn cuts it below the axis.
        if idx(1) > 1
            loglog(xr([1 1]), [y_bottom topr(1)], '-.', 'LineWidth', 1.5, ...
                   'Color', col, 'HandleVisibility', 'off');
        end
        if idx(end) < numel(mx)
            loglog(xr([end end]), [y_bottom topr(end)], '-.', 'LineWidth', 1.5, ...
                   'Color', col, 'HandleVisibility', 'off');
        end
        if any(topr < y_top * (1 - 1e-9))
            loglog(xr, topr, '-.', 'LineWidth', 1.5, 'Color', col, ...
                   'HandleVisibility', 'off');
        end
    end
end
end



function draw_mx_columns(iv, mx, col, alpha, lbl)
% Full-height shading of the excluded m_chi intervals iv (n x 2), drawn after
% the panel's y range is final so it cannot move it.  Dash-dot verticals mark
% the edges that fall inside the plotted mass range.
if isempty(iv), return; end
yl = ylim;
nv = size(iv, 1);
V = zeros(4*nv, 2);  F = reshape(1:4*nv, 4, nv)';
for r = 1:nv
    V(4*r-3:4*r, :) = [iv(r,1) yl(1); iv(r,2) yl(1); iv(r,2) yl(2); iv(r,1) yl(2)];
end
patch('Faces', F, 'Vertices', V, 'FaceColor', col, 'FaceAlpha', alpha, ...
      'EdgeColor', 'none', 'HandleVisibility', 'off');
for r = 1:nv
    for xe = iv(r, :)
        if xe > mx(1) && xe < mx(end)
            loglog([xe xe], yl, '-.', 'LineWidth', 1.5, 'Color', col, 'HandleVisibility', 'off');
        end
    end
end
legend_swatch(col, alpha, lbl, '-.');
ylim(yl);
end

function [lo, hi] = row_frame(panel, lo, hi, B)
% Widen a panel-local fill range [LO, HI] to the shared-row y range, if one is
% being imposed.  Fills ("excluded above/below this line") are drawn before
% apply_ylim, against the panel's own frame; when share_y then hands the panel a
% taller row union, every fill stopped at the old top, short of the axis edge.
if isstruct(B) && isfield(B, 'ylims') && isstruct(B.ylims) && ...
        isfield(B.ylims, panel) && numel(B.ylims.(panel)) == 2
    lo = min(lo, B.ylims.(panel)(1));
    hi = max(hi, B.ylims.(panel)(2));
end
end

function hi = row_top(panel, hi, B)
[~, hi] = row_frame(panel, hi, hi, B);
end

function apply_ylim(panel, computed, B, cap)
% Set the y range for PANEL, and publish what was used.
%
% A row of panels that shares one y axis must also share one y RANGE -- inner
% panels have no tick numbers, so their scale is read off the left-hand panel
% and a differing range would be silently misread.  run_figure_sets_g measures
% each row first, unions the ranges, and passes them back in as B.ylims.
lims = computed;
if isstruct(B) && isfield(B, 'ylims') && isstruct(B.ylims) && ...
        isfield(B.ylims, panel) && numel(B.ylims.(panel)) == 2
    lims = B.ylims.(panel);
end
% Ceiling, applied AFTER the shared-row override so a capped panel cannot be
% re-expanded by the union -- and BEFORE the publish below, so the measure pass
% unions the CAPPED range.  That ordering is the whole point: it is the light
% mediator whose top is tens of decades above everyone else's, so clamping only
% at draw time would leave the union (and therefore the row) uncapped.
%
% A panel can sit entirely above the cap -- at mphi = 1e-9 the massless-mediator
% sigma_chiN runs 1e-20 .. 1e-6 cm^2, all of it off the top.  Rather than invert
% the axis, drop the floor one decade under the cap: the panel goes (correctly)
% blank, and its one-decade range cannot drag the row's union downwards either.
if nargin >= 4 && isfinite(cap)
    lims(2) = min(lims(2), cap);
    if lims(1) >= lims(2), lims(1) = cap / 10; end
end
ylim(lims);
% Publish into the caller (this file is a SCRIPT, so 'caller' is the driver's
% workspace) for the measure pass to collect.
pv = struct();
if evalin('caller', 'exist(''PANEL_YLIM'', ''var'')')
    pv = evalin('caller', 'PANEL_YLIM');
end
pv.(panel) = lims;
assignin('caller', 'PANEL_YLIM', pv);
end

function decade_ticks(ax)
% Major tick on EVERY decade of a log axis, with x labelled at every decade and
% y labelled sparsely (a 20-decade sigma axis cannot carry 20 labels).  MATLAB's
% automatic log ticks drop whole decades when the range is wide, which loses the
% reader's sense of scale between two labelled powers.
for a = 'XY'
    if ~strcmp(get(ax, [a 'Scale']), 'log'), continue; end
    lim = get(ax, [a 'Lim']);
    e0 = ceil(log10(lim(1))); e1 = floor(log10(lim(2)));
    if ~isfinite(e0) || ~isfinite(e1) || e1 < e0, continue; end
    ex = e0:e1;
    set(ax, [a 'Tick'], 10.^ex);
    if a == 'X'
        step = 1;                        % label every decade
    else
        step = 1;                        % ... then thin to ~6 labels
        while numel(1:step:numel(ex)) > 6, step = step + 1; end
    end
    lbl = repmat({''}, 1, numel(ex));
    for k = 1:step:numel(ex)
        lbl{k} = sprintf('$10^{%d}$', ex(k));
    end
    set(ax, [a 'TickLabel'], lbl);
end
% No minor ticks on y: with a major tick on every decade already, the minor
% decade subdivisions just fur the axis.  X keeps its own default.
set(ax, 'YMinorTick', 'off', 'TickLabelInterpreter', 'latex');
end

function hide_y_labels(ax)
% Empty the y label and the y tick numbers on an inner panel of a shared-y row.
%
% Called AFTER tighten_axes, never before, and that is the whole point.  The
% axes is sized from its TightInset, which shrinks when the y label and numbers
% are gone -- so a panel that hid them early would end up with a WIDER plot box
% than the labelled panel it shares an axis with, and the two could not be
% aligned in the document at any single \includegraphics width.  Drawing them,
% measuring, and only then emptying the strings leaves every panel of the row
% with the same plot box and the same blank gutter, which the LaTeX side trims
% off (trim + clip) to butt the panels flush.
set(get(ax, 'YLabel'), 'String', '');
set(ax, 'YTickLabel', repmat({''}, 1, numel(get(ax, 'YTick'))));
end

function L = panel_layout()
% FIXED vertical geometry shared by every exported panel, in normalized figure
% units.  Panels sit side by side in the document, so their plot boxes have to
% be the same size and at the same height -- which they are not if each is sized
% from its own TightInset.
%
% Measured over all 15 sigma panels (2026-08-19): the left and bottom insets are
% identical everywhere (0.2250 / 0.2202 -- same fonts, same tick labels), but the
% TOP inset ranges 0.0000 to 0.0468 depending on whether the topmost y tick label
% lands on the axis top and protrudes above it.  Sizing the axes from that spread
% the plot-box height over 6% (0.7059 to 0.7660), and the export then baked the
% difference into the page.
%
% So the vertical layout is pinned to the WORST case instead of measured per
% panel.  Panels needing less carry a little blank above the axes -- the same
% amount on every panel, so they still line up.
%
% The horizontal layout is pinned the same way (2026-10-21).  Sized per panel,
% the left gutter followed the widest y tick label: 0.1637 on S/PS/SP/PP but
% 0.1820 on axialVector, whose plot box therefore came out 2% narrower and whose
% row, scaled by \panelrow to fill the line, sat at a different size.  And the
% right inset measured 0, putting the frame on the page edge, so half its line
% width was clipped where the panels abut.  ti_left is the axialVector worst
% case; ti_right is ~2 pt of air for the frame.
L.ti_bottom = 0.2202;   % room for the x tick labels and the x axis label
L.ti_top    = 0.0468;   % room for a top y tick label that protrudes
L.ti_left   = 0.1825;   % room for the y tick numbers and the y axis label
L.ti_right  = 0.0040;   % keeps the right frame line inside the page
end

function tighten_axes(ax)
% Grow the axes to fill the figure, leaving only the room the tick labels and
% axis labels actually need.  MATLAB's default LooseInset reserves a fixed
% margin sized for a title and a colorbar that this figure does not have, so
% the axes occupy ~62%% of the page where ~73%% is available.  Must run BEFORE
% inline_labels: it changes the axes aspect ratio, and label rotations are
% computed from that.
%
% All four insets are fixed by panel_layout so every panel's plot box is the
% same size and in the same place -- see there.  This panel's own TightInset is
% only checked against them.
L  = panel_layout();
ti = get(ax, 'TightInset');
if ti(1) > L.ti_left + 1e-3 || ti(2) > L.ti_bottom + 1e-3 || ti(4) > L.ti_top + 1e-3
    warning('compare_limits_g:panelInset', ...
        ['panel needs more room than panel_layout allows (left %.4f vs %.4f, ' ...
         'bottom %.4f vs %.4f, top %.4f vs %.4f): labels WILL be clipped.  Raise ' ...
         'the constants in panel_layout and re-render the whole set, or the row ' ...
         'stops aligning.'], ti(1), L.ti_left, ti(2), L.ti_bottom, ti(4), L.ti_top);
end
p = [L.ti_left, L.ti_bottom, 1 - L.ti_left - L.ti_right, 1 - L.ti_bottom - L.ti_top];
if p(3) > 0 && p(4) > 0
    % PositionConstraint has to be pinned too, or this Position does not stick.
    % It defaults to 'outerposition', under which MATLAB keeps OuterPosition
    % fixed and RE-DERIVES Position from the current TightInset on the next
    % relayout -- and inline_labels triggers one merely by flipping ax.Units to
    % points to read the axes size.  It doubles as the marker telling
    % enforce_panel_layout which axes have already been through here.
    set(ax, 'PositionConstraint', 'innerposition');
    set(ax, 'Position', p);
end
end


function h = place_labels_manually(h, spec)
% Move the inline labels named in spec = {label, [x y] [, rot]; ...} (data
% units; rot in degrees, default 0) to those positions, centred, at the
% inline size (every placed label shares it; max() guards a mixed set).  A label the
% placer dropped (no room, e.g. a narrow column) is created from its legend
% swatch.  A name that matches neither is reported, so a renamed constraint
% cannot silently drop its override.  A 4th column hide=true makes the label
% invisible (it is kept, so it can be un-hidden).  Returns h plus any
% labels created here.
if isempty(spec) || isempty(h), return; end
fs_full = max([h.FontSize]);
for k = 1:size(spec, 1)
    hit = h(strcmp(get(h, 'String'), spec{k, 1}));
    if isempty(hit)
        sw = findobj(gca, 'Type', 'patch', 'DisplayName', spec{k, 1});
        if isempty(sw)
            warning('compare_limits_g:labelPos', 'label_pos: no inline label "%s".', spec{k, 1});
            continue
        end
        c = sw(1).FaceColor;
        L = 0.2126*c(1) + 0.7152*c(2) + 0.0722*c(3);   % as inline_labels' clamp_col
        if L > 0.40, c = c * (0.40 / L); end
        hit = text(gca, 0, 0, spec{k, 1}, 'Interpreter', 'latex', 'Color', c, ...
                   'Clipping', 'on', 'HandleVisibility', 'off');
        h = [h(:); hit];
    end
    rot = 0;
    if size(spec, 2) >= 3 && ~isempty(spec{k, 3}), rot = spec{k, 3}; end
    set(hit, 'Units', 'data', 'Position', [spec{k, 2}, 0], 'Rotation', rot, ...
             'FontSize', fs_full, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hide = size(spec, 2) >= 4 && ~isempty(spec{k, 4}) && logical(spec{k, 4});
    set(hit, 'Visible', matlab.lang.OnOffSwitchState(~hide));
end
end

function show = rat_show_of(SHOW)
% Row order of the ratio constraints (gq_rat / sg_rat / vg_rat and the
% ext_rat_* colours and labels): kaon, kphi, cmb, integral, mesoninv.
show = [SHOW.kaon, SHOW.kphi, SHOW.cmb, SHOW.integral, SHOW.mesoninv];
end

function [gq_rat, sg_rat, vg_rat, gq_rat_hi, sg_rat_hi, vg_rat_hi] = ratio_curves( ...
    mx_vec, mphi, gx, ghat, interaction, interaction_sig, hbarc, rat_show, plot_svgg)
% gq_crit(mx) for each ratio constraint (row order: rat_show_of), and its image
% on the sigma_chiN and <sigma v>_gg axes.  Hidden rows skip the root-find.
% *_hi: the ceiling of the excluded band (NaN = excluded all the way up); only
% kphi (row 2) has one, so only it is searched.
nmx = numel(mx_vec);
nr  = numel(rat_show);
gq_rat = nan(nr, nmx);  gq_rat_hi = nan(nr, nmx);
parfor j = 1:nmx
    m = mx_vec(j);
    col = nan(nr, 1);  col_hi = nan(nr, 1);
    try
        rfun = { @(e) kaon_invisible_ratio(m, mphi, gx, e*ghat, interaction), ...
                 @(e) mesonandCMB_ratio(m, mphi, gx, e*ghat, interaction, 'kphi'), ...
                 @(e) mesonandCMB_ratio(m, mphi, gx, e*ghat, interaction, 'cmb'), ...
                 @(e) mesonandCMB_ratio(m, mphi, gx, e*ghat, interaction, 'integral'), ...
                 @(e) mesonandCMB_ratio(m, mphi, gx, e*ghat, interaction, 'mesoninv') };
        for ie = 1:nr
            if ~rat_show(ie), continue; end   %#ok<PFBNS> hidden -> skip the root-find
            % g_q^2 rescale as the seed, then root-find: the mediator width
            % breaks g_q^2 scaling wherever a pole is resonant (see
            % docs/kaon_constraint.md).  Off resonance the root-find returns
            % the seed to solver tolerance.
            r1 = rfun{ie}(1);
            if isfinite(r1) && r1 > 0
                col(ie) = ratio_crit_coupling(rfun{ie}, 1/sqrt(r1));
            end
            if ~isfinite(col(ie))
                % kphi falls to ~0 at g_q = 1 once the phi decays well inside the
                % detector (exp(-L/ell); for SxS from mphi ~ 0.013 GeV), while
                % still excluding smaller g_q.  Underflowed to 0, or tiny but
                % nonzero, so that 1/sqrt(ratio) lands tens of decades above the
                % band and the root-find runs out of range (was a spurious gap at
                % mphi = 14-17 MeV, SxS): seed from a tiny g_q instead, where the
                % phi escapes and the ratio is ~g_q^2.
                e0 = 1e-12;  r0 = rfun{ie}(e0);
                if isfinite(r0) && r0 > 0
                    col(ie) = ratio_crit_coupling(rfun{ie}, e0/sqrt(r0));
                end
            end
            if ie == 2
                col_hi(ie) = ratio_ceiling_coupling(rfun{ie}, col(ie), 1e8);
            end
        end
    catch err
        fprintf(2, '[ext-ratio] mx=%.4g: %s\n', m, err.message);
    end
    gq_rat(:, j) = col;  gq_rat_hi(:, j) = col_hi;
end
[sg_rat, vg_rat] = ratio_images(gq_rat, mx_vec, mphi, gx, ghat, interaction, ...
                                interaction_sig, hbarc, plot_svgg);
[sg_rat_hi, vg_rat_hi] = ratio_images(gq_rat_hi, mx_vec, mphi, gx, ghat, interaction, ...
                                      interaction_sig, hbarc, plot_svgg);
end

function [sg, vg] = ratio_images(E, mx_vec, mphi, gx, ghat, interaction, ...
                                 interaction_sig, hbarc, plot_svgg)
% sigma_chiN and <sigma v>_gg at each g_q in E (rows = ratio constraints).
sg = nan(size(E));   vg = nan(size(E));
for ie = 1:size(E, 1)
    for j = 1:numel(mx_vec)
        e = E(ie, j);
        if ~isfinite(e), continue; end
        sg(ie, j) = DMnucxsec(mx_vec(j), mphi, gx, e*ghat, interaction_sig) * hbarc^2;
        if plot_svgg, vg(ie, j) = svgg_gq(mx_vec(j), mphi, gx, e, ghat, interaction, hbarc); end
    end
end
end

function [lo, hi] = overclosed_interval(le, OA, target)
% g_q interval on the coarse loop grid where a stable mediator's own relic
% exceeds target.  Edges are log-interpolated in Omega against the neighbouring
% grid point (grid step ~0.75 dex); an edge at the end of the grid stays at
% that grid point.  Disjoint pieces are merged into one span.  NaN if none.
lo = NaN;  hi = NaN;
ov = isfinite(OA) & OA > target;
if ~any(ov), return, end
i1 = find(ov, 1, 'first');  i2 = find(ov, 1, 'last');
lo = 10^oc_edge(le, OA, target, i1, i1 - 1);
hi = 10^oc_edge(le, OA, target, i2, i2 + 1);
end

function l = oc_edge(le, OA, target, i, n)
l = le(i);
if n < 1 || n > numel(le) || ~(isfinite(OA(n)) && OA(n) > 0), return, end
f = log(OA(i)/target) / log(OA(i)/OA(n));
l = le(i) + f*(le(n) - le(i));
end

function [mxa, ia, ins] = relic_tip_columns(mx, tip)
% The mx grid with a column inserted just below each mass j whose low-mass
% island wall relic_tip_mx(j) is finite.  ia maps every column to the grid mass
% whose loop-branch couplings it carries; ins flags the inserted columns.
n = numel(mx) + sum(isfinite(tip));
mxa = zeros(1, n);  ia = zeros(1, n);  ins = false(1, n);
k = 0;
for j = 1:numel(mx)
    if isfinite(tip(j))
        k = k + 1;  mxa(k) = tip(j);  ia(k) = j;  ins(k) = true;
    end
    k = k + 1;  mxa(k) = mx(j);  ia(k) = j;
end
end

function draw_relic_walls(mxa, ins, gq_lo, gq_hi, col, yof)
% The vertical Oh2 = 0.12 segment closing an island at each inserted column,
% from the lower to the upper loop branch; solid up to g_q = 1, dashed above,
% as the branches themselves.  yof(m, g_q) is the plotted quantity.
for k = find(ins)
    lo = gq_lo(k);  hi = gq_hi(k);
    if ~(isfinite(lo) && isfinite(hi)), continue; end
    m = mxa(k);  e1 = min(max(1, lo), hi);
    if lo < e1
        loglog([m m], [yof(m, lo) yof(m, e1)], '-',  'LineWidth', 2, 'Color', col, 'HandleVisibility', 'off');
    end
    if e1 < hi
        loglog([m m], [yof(m, e1) yof(m, hi)], '--', 'LineWidth', 2, 'Color', col, 'HandleVisibility', 'off');
    end
end
end


function tips = relic_island_tips(mx_vec, gq_lo, gq_hi, max_dec, oh_at, target)
% Low-mass walls of the relic islands; see the loop block.  At the first mass j
% of a two-branch run whose branches are more than MAX_DEC decades apart, bisect
% (6 steps in log mx) between mx(j-1) and mx(j) for where Oh2 at the branches'
% mid coupling crosses TARGET.  No wall if Oh2 there is not below TARGET at
% mx(j-1): the lower branch carries on below, open-topped.  OH_AT(m, g_q) is
% Oh2; TIPS(j) is the wall mass, NaN where there is none.
tips = nan(size(mx_vec));
both = ~isnan(gq_lo) & ~isnan(gq_hi);
for j = 2:numel(mx_vec)
    if ~both(j) || both(j-1), continue; end
    if abs(log10(gq_hi(j) / gq_lo(j))) <= max_dec, continue; end
    gq_mid = sqrt(gq_lo(j) * gq_hi(j));
    la = log(mx_vec(j-1));  lb = log(mx_vec(j));
    try
        fa = oh_at(exp(la), gq_mid) - target;
        if ~(isfinite(fa) && fa < 0), continue; end    % no wall below mx(j)
        for it = 1:6
            lc = 0.5 * (la + lb);
            fc = oh_at(exp(lc), gq_mid) - target;
            if ~isfinite(fc), break; end
            if fc < 0, la = lc; else, lb = lc; end
        end
    catch err
        fprintf(2, '[Loop tip] mx=%.4g: %s\n', mx_vec(j), err.message);
        continue
    end
    tips(j) = exp(0.5 * (la + lb));
end
end

function f = relic_oh2_at(solver_name, mphi, gx, ghat, interaction, T_RH, production)
% Oh2(m, g_q) with the relic solver this panel uses.
if strcmp(solver_name, 'freeze_coupled')
    f = @(m, gq) freeze_coupled(m, mphi, gx, gq * ghat, interaction, T_RH, false, production);
else
    f = @(m, gq) freeze_loop(m, mphi, gx, gq * ghat, interaction, T_RH, false, production);
end
end
