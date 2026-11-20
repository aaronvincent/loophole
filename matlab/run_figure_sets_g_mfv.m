% run_figure_sets_g_mfv.m
% MFV mirror of run_figure_sets_g.m: the same eta-scan sweep, but along the
% MINIMAL FLAVOUR VIOLATION coupling direction rather than the democratic one,
% and restricted to the two families where MFV is meaningful.
%
% Under MFV the quark couplings are proportional to the Yukawas,
%   g = eta * ghat,   ghat = [m_u, m_d, m_s] / (m_u + m_d + m_s),
% so the scan's single free normalization eta means exactly what it means in the
% democratic sweep (ghat sums to 1), and only the DIRECTION differs.  This is the
% same ansatz the GAMBIT MFV scans use (gu:gd:gs tied by `same_as` to these
% weights; see README.md).
%
% ONLY SS AND PS.  MFV is a statement about the SCALAR quark current, so it
% applies to the two families whose phi-SM vertex is scalar (heavy/lightScalar
% and heavy/lightScalarPS).  The SP and PP families couple to the PSEUDOSCALAR
% quark current, where "couplings proportional to the Yukawas" is not the MFV
% statement, and axialVector has its own structure -- all three are excluded
% here for the same reason gen_mfv_jobs.sh excludes axialVector.
%
%   models:  the SS and PS families only (RUNSET.families), with the mediator
%            regime chosen by mass -- light* at the eV-scale mphi, heavy* at the
%            two GeV-scale ones (see the job-list block for what that does/doesn't
%            change), so a family is 3 runs, not 6
%   mphi:    1e-9, 0.1, 1 GeV   (three figure sets per family)
%   gx:      fixed (RUNSET.gx),  g_ref: ghat (MFV direction, RUNSET.g_ref)
%
% Each run produces up to three compare_limits_g panels (toggled individually by
% RUNSET.plot_sigma / plot_eta / plot_svgg), named
% <figtag> = compare_limits_g_<interaction>_mphi<mphi>_gx<gx>_mfv --
%   <figtag>_sigma.pdf   sigma_chiN vs m_chi
%   <figtag>_eta.pdf     eta       vs m_chi
%   <figtag>_svgg.pdf    <sigma v>_{gg -> chi chi} vs m_chi
% gx is in the name so runs at different RUNSET.gx accumulate side by side.  The
% _mfv suffix is what keeps this sweep from colliding with the democratic one:
% the .dat/.mat stem compare_limits_g.m builds is
% <interaction>_mphi<mphi>_gx<gx> and carries NO record of the coupling
% direction, so without RUNSET.dat_tag = '_mfv' below an MFV run at the same
% (model, mphi, gx) would silently overwrite the democratic run's data.
% -- plus the .dat file compare_limits_g.m always writes, under output/.  The
% .dat holds every column regardless of the toggles; only <sigma v>_gg, whose
% block costs an extra svavg call per point, is skipped when its panel is off.
%
% compare_limits_g.m is invoked as a SCRIPT (calling it as a function trips the
% MATLAB parfor "source file not found" error), so it shares this workspace and
% overwrites names like `interaction`, `mphi`, `ii`, `j`, `jj`.  All driver state
% therefore lives in the single struct RUNSET, which compare_limits_g.m never
% touches, and the loop is driven by RUNSET.idx rather than a loop variable.
% Parameters reach compare_limits_g.m through `batch_params`.
%
% Runtime warning: a full compare_limits_g scan is minutes, so 2 families x 3
% masses is under an hour.  Trim RUNSET.families / RUNSET.mphi_vec for a subset, or
% set RUNSET.skip_existing = true to resume an interrupted sweep.

% -----------------------------------------------------------------------
% Settings -- edit these
% -----------------------------------------------------------------------
RUNSET = struct();
RUNSET.gx       = .01;                 % FIXED dark coupling

% MFV coupling direction: g = eta * ghat with ghat = [m_u, m_d, m_s]/sum, so
% ghat sums to 1 and eta keeps the meaning it has in the democratic sweep.
% The three masses are the ones nucleon_couplings.m:28 uses (2.2, 4.69, 95 MeV);
% keep them in sync with that file and with an MFV YAML, which hardcodes the
% resulting [0.02159191, 0.04603003, 0.93237805].
RUNSET.mfv_quark_masses = [0.0022, 0.00469, 0.095];      % [m_u m_d m_s], GeV
RUNSET.g_ref    = RUNSET.mfv_quark_masses / sum(RUNSET.mfv_quark_masses);
RUNSET.mphi_vec = [1e-9, 0.1, 1];    % mediator masses [GeV] -- 3 sets per model

% One row per physics family: {light variant, heavy variant}.  The variant used
% at a given mphi is chosen by mass, NOT run as two separate models: the eV-scale
% mediator gets the light* string, the two GeV-scale ones get heavy*.  Families
% with no heavy/light pair (axialVector) repeat the same string.
% NB the bare 'scalarPP' is absent because it is not a separate family but the
% name heavy*/light*PP already RESOLVE to: the heavy/light prefix is inert for
% the relic curves, full_sigma_name maps both to 'scalarPP' for sigma_chiN, and
% svavg whitelists and dispatches all three identically.  A 'scalarPP' row
% reproduced the ScalarPP row column-for-column (verified 2026-08-19: all 30
% .dat columns bit-identical, panels byte-identical).  Do not re-add it.
% NB the bare 'scalarPS' is deliberately absent: svavg's 'gg' channel whitelist
% omits it, which would leave the <sigma v> panel empty.  Use heavy*/light*PS.
% SS and PS only -- see the ONLY SS AND PS note in the header for why the SP,
% PP and axialVector rows of run_figure_sets_g.m are absent here.
RUNSET.families = { ...
    {'lightScalar',   'heavyScalar'  }, ...
    {'lightScalarPS', 'heavyScalarPS'} ...
};
RUNSET.mphi_light_max = 1e-6;        % mphi <= this -> light*; above -> heavy*
RUNSET.mx_vec  = logspace(-5, log10(0.200), 100);   % DM mass grid [GeV]
RUNSET.T_RH    = 0.01;               % reheat temperature [GeV]
RUNSET.production = 'sm';         % 'sm' (gamma gamma + meson pairs) or 'gg' (photons only)
RUNSET.figdir  = 'figures_g_mfv';    % output directory for the figures
RUNSET.figfmt  = {'pdf'};            % e.g. {'pdf','png'} or add 'fig'
RUNSET.figsize = [18 13];            % saved page size [w h] in cm (~7 x 5 in).  The
                                     % PDF page IS this size (printpdf sets PaperSize
% Font sizes in FINAL document points -- what the reader measures once
% \includegraphics has scaled the page.  panel_cm is the width ONE panel is
% placed at (0.32\textwidth on a 2 cm-margin letter page = 5.63 cm; use
% 0.48\textwidth = 8.45 cm for two-across).  Getting panel_cm wrong is the
% only way these come out the wrong size.
RUNSET.style = 'compact';            % 'compact' (packed: decade ticks, no grid,
                                     %   tight axes, doc-scaled fonts) or 'classic'
                                     %   (grid on, auto ticks, loose axes, absolute
                                     %   16/14/13/11 pt -- what the inline labels
                                     %   were tuned against).
RUNSET.sigma_ycap = 1e-25;           % ceiling on the sigma panel's y axis [cm^2],
                                     %   'compact' only -- the light-mediator
                                     %   panel otherwise sets a shared row range
                                     %   tens of decades taller than the curves.
                                     %   Inf disables; ignored under 'classic'.
RUNSET.font = struct('label', 7, 'tick', 6.5, 'legend', 6, 'inline', 4.7, ...
                     'panel_cm', 5.63);
                                     % from it); leave the font sizes in the plotter
                                     % alone unless you change this a lot.
RUNSET.dat_tag = '_mfv';             % suffix on the .dat/.mat names.  REQUIRED
                                     % here, not cosmetic: the stem
                                     % compare_limits_g.m builds records the
                                     % model, mphi and gx but NOT the coupling
                                     % direction, so an MFV run would otherwise
                                     % overwrite the democratic sweep's data at
                                     % the same (model, mphi, gx).  Append to it
                                     % (e.g. '_mfv_smoke') for test runs.
RUNSET.float_y_max_mphi = 1e-6;      % panels with mphi <= this are EXCLUDED from
                                     %   the shared range: they keep their own y
                                     %   axis, their own label and numbers, and no
                                     %   sigma_ycap.  The massless-mediator
                                     %   sigma_chiN sits tens of decades above the
                                     %   rest, so sharing with it either empties
                                     %   the other panels or (once capped) empties
                                     %   itself.  Inf -> nothing floats.
RUNSET.share_y = true;               % one y RANGE per model row, and y label +
                                     %   tick numbers only on the leftmost panel.
                                     %   Costs a measure pass per row (cheap: it
                                     %   builds the cache the render pass reuses).
RUNSET.replot = false;               % true -> reuse each job's cached .mat and only
                                     %   re-render (seconds, not hours).  Use after a
                                     %   legend / colour / SHOW change.  Errors if a
                                     %   job has no cache yet.
                                     % 'ext' -> same, but RE-RUN the external nucleon
                                     %   bands (SN1987A / neutron star / stellar /
                                     %   5th force) and write them back to the .mat.
                                     %   Use after adding or changing one of those
                                     %   constraints: the relic + Delta-Neff blocks
                                     %   are reused, the band scan costs seconds.
RUNSET.skip_existing = false;        % true -> skip a job whose figures already exist
RUNSET.close_figs    = true;         % close each figure after saving it.  Leave true
                                     % for a full sweep (6 families x 3 masses x 3
                                     % panels = 54 windows); set false when running a
                                     % short job list and you want to inspect them.

% Constraint overlays.  Any subset of compare_limits_g.m's SHOW fields; anything
% omitted keeps that script's default (all true).  E.g. to drop the indirect
% bounds from every figure in the sweep:
%   RUNSET.show = struct('cmb', false, 'integral', false);
RUNSET.show = struct();

% Manual inline-label positions, keyed by figtag; see run_figure_sets_g.m.
RUNSET.label_pos = label_pos_load();   % label_pos_g.json

% Which panels to draw.  All three come from the same physics scan, so turning
% one off saves plotting time, not solver time -- except plot_svgg, which gates
% compare_limits_g.m's extra <sigma v>_gg evaluation block.
RUNSET.plot_sigma = true;            % sigma_chiN vs m_chi
RUNSET.plot_eta   = false;            % eta       vs m_chi
RUNSET.plot_svgg  = false;            % <sigma v>_{gg -> chi chi} vs m_chi

% Overrides from a caller: any RUNSET field, applied after the
% settings above.  Consumed and cleared, like batch_params.
if exist('RUNSET_OVERRIDE', 'var') && isstruct(RUNSET_OVERRIDE)
    RUNSET.ovf = fieldnames(RUNSET_OVERRIDE);
    RUNSET.k = 1;
    while RUNSET.k <= numel(RUNSET.ovf)
        RUNSET.(RUNSET.ovf{RUNSET.k}) = RUNSET_OVERRIDE.(RUNSET.ovf{RUNSET.k});
        RUNSET.k = RUNSET.k + 1;
    end
    clear RUNSET_OVERRIDE
end

% -----------------------------------------------------------------------
% Build the job list: one job per (family, mphi), with the heavy/light variant
% picked by mphi.  NB with full_sigma_name in the loop (compare_limits_g.m
% rewrites heavy*/light* -> the full 'scalar*' sibling for sigma_chiN, and the
% relic solvers treat the prefix as an alias), this choice currently affects only
% the model NAME in the filenames and titles, not the numbers.  It is the right
% pairing to record, and becomes load-bearing if the full-sigma rewrite is ever
% dropped in favour of the regime approximations.
% (Plain loop variables are safe here: this runs before compare_limits_g.m is
% ever invoked, and nothing below reuses them.)
% -----------------------------------------------------------------------
RUNSET.jobs = {};
for ia = 1:numel(RUNSET.families)
    for im = 1:numel(RUNSET.mphi_vec)
        if RUNSET.mphi_vec(im) <= RUNSET.mphi_light_max
            nm = RUNSET.families{ia}{1};   % eV-scale mediator -> light variant
        else
            nm = RUNSET.families{ia}{2};   % GeV-scale mediator -> heavy variant
        end
        RUNSET.jobs{end+1} = {nm, RUNSET.mphi_vec(im)};
    end
end
% Panel suffixes actually being produced -- drives the skip_existing check.
RUNSET.panels = {};
if RUNSET.plot_sigma, RUNSET.panels{end+1} = 'sigma'; end
if RUNSET.plot_eta,   RUNSET.panels{end+1} = 'eta';   end
if RUNSET.plot_svgg,  RUNSET.panels{end+1} = 'svgg';  end
if isempty(RUNSET.panels)
    warning('run_figure_sets_g:noPanels', ...
        'plot_sigma/plot_eta/plot_svgg are all false: the sweep will write .dat files only.');
end

RUNSET.n_jobs = numel(RUNSET.jobs);
RUNSET.failed = {};
RUNSET.skipped = {};
RUNSET.t_all = tic;

if ~isfolder(RUNSET.figdir), mkdir(RUNSET.figdir); end

% -----------------------------------------------------------------------
% Run.  while + RUNSET.idx: a for-loop counter would be clobbered by
% compare_limits_g.m, which uses ii/j/jj/kk in its own workspace.
% -----------------------------------------------------------------------
% Jobs are family-major (3 mphi per family, in RUNSET.mphi_vec order), so one
% family is one contiguous group -- and one ROW of the figure layout, left to
% right.  With share_y that group is rendered twice: a measure pass to learn
% each panel's own y range, then a render pass with the union applied and the y
% label + tick numbers dropped on every panel that shares an axis to its left.
% Panels at mphi <= RUNSET.float_y_max_mphi are excluded from the union: they
% keep their own range, their own cap and their own labels, so a row is one
% floating panel plus a shared pair, not three panels on one unusable range.
%
% The measure pass is not wasted work.  It runs with RUNSET.replot as
% configured, so it is what builds the .mat cache; the render pass then forces
% replot = true and costs seconds.  The expensive solve still happens once.
RUNSET.group_size = max(1, numel(RUNSET.mphi_vec));
RUNSET.g0 = 1;
while RUNSET.g0 <= RUNSET.n_jobs
    RUNSET.gend  = min(RUNSET.g0 + RUNSET.group_size - 1, RUNSET.n_jobs);
    RUNSET.gidx  = RUNSET.g0:RUNSET.gend;
    RUNSET.gok   = true(1, numel(RUNSET.gidx));
    RUNSET.shared = struct();

    % ---------------- measure pass (share_y only) ----------------------
    RUNSET.gi = 1;
    while RUNSET.share_y && RUNSET.gi <= numel(RUNSET.gidx)
        RUNSET.idx = RUNSET.gidx(RUNSET.gi);
        RUNSET.this_interaction = RUNSET.jobs{RUNSET.idx}{1};
        RUNSET.this_mphi        = RUNSET.jobs{RUNSET.idx}{2};
        RUNSET.floats = RUNSET.this_mphi <= RUNSET.float_y_max_mphi;
        if RUNSET.floats
            RUNSET.what = 'building cache; y axis floats';
        else
            RUNSET.what = 'measuring y range';
        end
        fprintf('\n=== [%d/%d] %s  mphi=%.4g GeV  (%s) ===\n', ...
                RUNSET.idx, RUNSET.n_jobs, RUNSET.this_interaction, ...
                RUNSET.this_mphi, RUNSET.what);
        batch_params = struct();
        batch_params.interaction = RUNSET.this_interaction;
        batch_params.mphi        = RUNSET.this_mphi;
        batch_params.gx          = RUNSET.gx;
        batch_params.g_ref       = RUNSET.g_ref;
        batch_params.T_RH        = RUNSET.T_RH;
        batch_params.production  = RUNSET.production;
        batch_params.mx_vec      = RUNSET.mx_vec;
        batch_params.figsize     = RUNSET.figsize;
        batch_params.font        = RUNSET.font;
        batch_params.style       = RUNSET.style;
        if RUNSET.floats, batch_params.sigma_ycap = Inf;
        else,             batch_params.sigma_ycap = RUNSET.sigma_ycap; end
        batch_params.dat_tag     = RUNSET.dat_tag;
        batch_params.replot      = RUNSET.replot;
        batch_params.close_figs  = true;
        batch_params.show        = RUNSET.show;
        batch_params.plot_sigma_chiN = RUNSET.plot_sigma;
        batch_params.plot_eta        = RUNSET.plot_eta;
        batch_params.plot_svgg       = RUNSET.plot_svgg;
        % no figdir -> nothing is written; this pass only reads PANEL_YLIM back
        PANEL_YLIM = struct();
        try
            compare_limits_g
        catch RUNSET_err
            fprintf(2, '[run_figure_sets_g_mfv] MEASURE FAILED %s mphi=%.4g: %s\n', ...
                    RUNSET.this_interaction, RUNSET.this_mphi, RUNSET_err.message);
            RUNSET.gok(RUNSET.gi) = false;
        end
        close_unkept_figs
        % A floating panel still runs here -- this pass is what builds its .mat
        % cache -- but its range is deliberately left out of the union.
        if ~RUNSET.floats && RUNSET.gok(RUNSET.gi) && ...
                exist('PANEL_YLIM', 'var') && isstruct(PANEL_YLIM)
            RUNSET.pf = fieldnames(PANEL_YLIM);
            RUNSET.kf = 1;
            while RUNSET.kf <= numel(RUNSET.pf)
                RUNSET.nm = RUNSET.pf{RUNSET.kf};
                RUNSET.v  = PANEL_YLIM.(RUNSET.nm);
                if isfield(RUNSET.shared, RUNSET.nm)
                    RUNSET.shared.(RUNSET.nm) = ...
                        [min(RUNSET.shared.(RUNSET.nm)(1), RUNSET.v(1)), ...
                         max(RUNSET.shared.(RUNSET.nm)(2), RUNSET.v(2))];
                else
                    RUNSET.shared.(RUNSET.nm) = RUNSET.v;
                end
                RUNSET.kf = RUNSET.kf + 1;
            end
        end
        RUNSET.gi = RUNSET.gi + 1;
    end

    % ---------------- render pass --------------------------------------
    % The y label + numbers go on the leftmost panel of each y axis: every
    % floating panel (it owns its axis), plus the first panel of the shared set.
    RUNSET.gi_first_shared = 0;
    RUNSET.gi = 1;
    while RUNSET.gi <= numel(RUNSET.gidx)
        if RUNSET.jobs{RUNSET.gidx(RUNSET.gi)}{2} > RUNSET.float_y_max_mphi
            RUNSET.gi_first_shared = RUNSET.gi;
            break;
        end
        RUNSET.gi = RUNSET.gi + 1;
    end

    RUNSET.gi = 1;
    while RUNSET.gi <= numel(RUNSET.gidx)
        RUNSET.idx = RUNSET.gidx(RUNSET.gi);
        RUNSET.this_interaction = RUNSET.jobs{RUNSET.idx}{1};
        RUNSET.this_mphi        = RUNSET.jobs{RUNSET.idx}{2};
        RUNSET.floats = RUNSET.this_mphi <= RUNSET.float_y_max_mphi;
        % gx is in the name so runs at different RUNSET.gx do not overwrite each
        % other; RUNSET.dat_tag ('_mfv') separates this sweep from the democratic
        % one, whose figures share the model/mphi/gx stem.  Same convention as the
        % .dat/.mat the run writes beside them.
        RUNSET.figtag = sprintf('compare_limits_g_%s_mphi%.4g_gx%.4g%s', ...
                                RUNSET.this_interaction, RUNSET.this_mphi, ...
                                RUNSET.gx, RUNSET.dat_tag);
        RUNSET.have_all = ~isempty(RUNSET.panels);
        for kp = 1:numel(RUNSET.panels)
            RUNSET.have_all = RUNSET.have_all && isfile(fullfile(RUNSET.figdir, ...
                sprintf('%s_%s.%s', RUNSET.figtag, RUNSET.panels{kp}, RUNSET.figfmt{1})));
        end
        if RUNSET.skip_existing && RUNSET.have_all
            fprintf('=== [%d/%d] %s  mphi=%.4g GeV -- SKIPPED (figures exist) ===\n', ...
                    RUNSET.idx, RUNSET.n_jobs, RUNSET.this_interaction, RUNSET.this_mphi);
            RUNSET.skipped{end+1} = RUNSET.figtag;
            RUNSET.gi = RUNSET.gi + 1;
            continue;
        end
        if ~RUNSET.gok(RUNSET.gi)
            RUNSET.failed{end+1} = RUNSET.figtag;   % measure pass already failed
            RUNSET.gi = RUNSET.gi + 1;
            continue;
        end

        fprintf('\n=== [%d/%d] %s  mphi=%.4g GeV ===\n', RUNSET.idx, RUNSET.n_jobs, ...
                RUNSET.this_interaction, RUNSET.this_mphi);

        batch_params = struct();
        batch_params.interaction = RUNSET.this_interaction;
        batch_params.mphi        = RUNSET.this_mphi;
        batch_params.gx          = RUNSET.gx;
        batch_params.g_ref       = RUNSET.g_ref;
        batch_params.T_RH        = RUNSET.T_RH;
        batch_params.production  = RUNSET.production;
        batch_params.mx_vec      = RUNSET.mx_vec;
        batch_params.figdir      = RUNSET.figdir;
        batch_params.figtag      = RUNSET.figtag;
        batch_params.figfmt      = RUNSET.figfmt;
        batch_params.figsize     = RUNSET.figsize;
        batch_params.font        = RUNSET.font;
        batch_params.style       = RUNSET.style;
        batch_params.dat_tag     = RUNSET.dat_tag;
        batch_params.close_figs  = RUNSET.close_figs;
        batch_params.plot_sigma_chiN = RUNSET.plot_sigma;
        batch_params.plot_eta        = RUNSET.plot_eta;
        batch_params.plot_svgg       = RUNSET.plot_svgg;
        batch_params.show            = RUNSET.show;
        if isKey(RUNSET.label_pos, RUNSET.figtag)
            batch_params.label_pos   = RUNSET.label_pos(RUNSET.figtag);
        end
        if RUNSET.share_y
            batch_params.replot      = true;              % cache built by the measure pass
            if RUNSET.floats
                % Own axis: no shared range, no cap, and it keeps its numbers --
                % they cannot be read off the shared panels' scale.
                batch_params.sigma_ycap  = Inf;
                batch_params.show.ylabel = true;
            else
                batch_params.sigma_ycap  = RUNSET.sigma_ycap;
                batch_params.ylims       = RUNSET.shared;
                batch_params.show.ylabel = (RUNSET.gi == RUNSET.gi_first_shared);
            end
        else
            batch_params.replot      = RUNSET.replot;
            batch_params.sigma_ycap  = RUNSET.sigma_ycap;
        end

        try
            compare_limits_g        % invoked as a script, not a function
        catch RUNSET_err
            fprintf(2, '[run_figure_sets_g_mfv] FAILED %s mphi=%.4g: %s\n', ...
                    RUNSET.this_interaction, RUNSET.this_mphi, RUNSET_err.message);
            RUNSET.failed{end+1} = RUNSET.figtag;
            close_unkept_figs
        end
        RUNSET.gi = RUNSET.gi + 1;
    end

    RUNSET.g0 = RUNSET.gend + 1;
end

% Belt and braces: compare_limits_g.m already consumes (and clears) batch_params
% at its top, so it cannot leak into a later interactive run.
clear batch_params

fprintf('\n=== run_figure_sets_g_mfv done: %d/%d run, %d skipped, %d failed, %.1f s ===\n', ...
        RUNSET.n_jobs - numel(RUNSET.skipped) - numel(RUNSET.failed), RUNSET.n_jobs, ...
        numel(RUNSET.skipped), numel(RUNSET.failed), toc(RUNSET.t_all));
if ~isempty(RUNSET.failed)
    fprintf(2, 'Failed jobs:\n');
    for kf = 1:numel(RUNSET.failed)
        fprintf(2, '  %s\n', RUNSET.failed{kf});
    end
end
fprintf('Figures in %s/\n', RUNSET.figdir);
