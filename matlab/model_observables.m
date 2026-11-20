function obs = model_observables(mx, mphi, gx, g, interaction, T_RH, s_mix)
% model_observables  Map a single model point to its observables.
%
%   obs = model_observables(mx, mphi, gx, g, interaction [, T_RH, s_mix])
%
% One-stop, single-point evaluator intended for a likelihood call: takes the
% model parameters (masses, couplings, vertex structure) and returns the
% observables (relic abundance, Delta-Neff at the CMB, DM-nucleon cross
% section).  Each observable is computed in its own try/catch so one failure
% (e.g. an ode15s non-convergence at a stiff point, or the axialVector
% ultralight-mediator hard error) degrades to NaN for that field rather than
% aborting the whole point -- a sampler should map a NaN observable to a
% rejected / zero-likelihood draw.
%
% Inputs:
%   mx          DM mass m_chi [GeV]
%   mphi        mediator mass m_phi (or m_A' for axialVector) [GeV]
%   gx          DM-mediator coupling (rescaled internally by nothing here --
%               this is the *physical* coupling at which observables are read)
%   g           [g_u g_d g_s] quark couplings (sets the chiral loop D)
%   interaction vertex string ('heavyScalar', 'heavyScalarPS', 'axialVector', ...)
%   T_RH        reheat temperature [GeV].  Default 0.01 (10 MeV).
%   s_mix       axialVector small-mixing-breakdown threshold sin(theta): above
%               this the leading-order A'-meson mixing expansion is suspect.
%               Default 0.3 (matches compare_limits_g.m).  Ignored for scalars.
%
% Output struct obs:
%   .Oh2          relic abundance Omega_chi h^2 (authoritative loop solver)
%   .Oh2_mediator Omega h^2 of a mediator that is stable on the age of the
%                 Universe (freeze_coupled's info.Omega_phi_h2; in practice the
%                 axialVector A' below min(2 mx, 3 m_pi)); NaN from freeze_loop
%   .mediator_overclosed  stable mediator with Oh2_mediator > 0.12 (not folded
%                 into .Oh2 or .eft.valid -- the caller decides)
%   .dNeff_bbn    NaN (neff_simple has no BBN epoch)
%   .dNeff_cmb    dark-sector Delta-Neff at the CMB epoch  (signed; can be < 0)
%   .sigma_chiN   DM-nucleon cross section [cm^2]
%   .solver       'freeze_loop' or 'freeze_coupled' (which relic solver ran)
%   .neff         neff_simple info struct (contact weights, heating and
%                 phi-radiation pieces) -- [] if the Neff computation errored
%   .eft          EFT-validity flags for this point.  Always present; trivially
%                 valid for scalar mediators (no such boundary), meaningful only
%                 for axialVector:
%                   .valid          both boundaries satisfied (gate for a fit)
%                   .unitarity_ok   gx <= gx_ceiling
%                   .gx_ceiling     largest EFT-valid gx: min(axial_gx_ceiling,
%                                   sqrt(pi/2) mphi/mx); Inf if not axialVector
%                   .small_mixing_ok  sin(theta) < s_mix
%                   .sin_mix        sin(theta) = max_j|m_{jA'}|/mphi  (gx-INdependent)
%                   .s_mix          the threshold used
%   .errors       struct of caught error messages, one field per failed
%                 observable ('Oh2'/'neff'/'sigma'/'eft'); empty if all succeeded
%   .params       echo of the input point (mx, mphi, gx, g, interaction, T_RH, s_mix)
%
% Convention notes:
%   * gx is the physical coupling.  The compare_*.m scripts instead root-find
%     gx such that Oh2 = 0.12; here you supply gx and read what comes out --
%     the natural primitive for a likelihood that varies gx as a free
%     parameter.
%   * The relic value is the self-consistent loop solver (freeze_loop /
%     freeze_coupled), authoritative across the FIMP<->WIMP crossing.  It is
%     NOT the linearized Oh2_FI or the equilibrium-IC freeze_out (those are
%     each valid only on one side; see docs/boltzmann_solvers.md).
%   * Solver dispatch: freeze_loop treats phi as a virtual s-channel state and
%     is used for mphi > 50*T_RH; below that the on-shell phi is a bath species
%     and freeze_coupled tracks (Y_chi, Y_phi).
%
% Example:
%   o = model_observables(0.1, 10, 1, [1 1 1]/3, 'heavyScalarPS');
%   o.Oh2, o.dNeff_bbn, o.dNeff_cmb, o.sigma_chiN

if nargin < 6 || isempty(T_RH), T_RH = 0.01; end   % 10 MeV default
if nargin < 7 || isempty(s_mix), s_mix = 0.3;  end   % small-mixing threshold

alltheconstants;        % loads hbarc, Mpl, s0, rhoc, ... into this workspace

obs = struct();
obs.params = struct('mx', mx, 'mphi', mphi, 'gx', gx, 'g', g, ...
                    'interaction', interaction, 'T_RH', T_RH, 's_mix', s_mix);
obs.errors = struct();

% ------------------------------------------------------------------
% Relic abundance -- authoritative self-consistent loop solver.
% Same mphi-vs-T_RH dispatch as compare_limits_g.m.
% ------------------------------------------------------------------
if mphi > 50 * T_RH
    obs.solver = 'freeze_loop';
    relic_fn   = @freeze_loop;
else
    obs.solver = 'freeze_coupled';
    relic_fn   = @freeze_coupled;
end
obs.Oh2_mediator        = NaN;     % freeze_loop: phi is virtual, no relic
obs.mediator_overclosed = false;
try
    if strcmp(obs.solver, 'freeze_coupled')
        [obs.Oh2, finfo] = relic_fn(mx, mphi, gx, g, interaction, T_RH);
        obs.Oh2_mediator        = finfo.Omega_phi_h2;
        obs.mediator_overclosed = finfo.overclosed;
    else
        obs.Oh2 = relic_fn(mx, mphi, gx, g, interaction, T_RH);
    end
catch ME
    obs.Oh2          = NaN;
    obs.errors.Oh2   = ME.message;
end

% ------------------------------------------------------------------
% Delta-Neff (neff_simple): CMB epoch only.  dNeff_bbn is NaN, which the
% scan's neff_loglike counts as no contribution.
% ------------------------------------------------------------------
try
    [dNeff_cmb, ninfo] = neff_simple(mphi, mx, gx, g, interaction, struct('T_RH', T_RH));
    obs.dNeff_bbn = NaN;
    obs.dNeff_cmb = dNeff_cmb;
    obs.neff      = ninfo;
catch ME
    obs.dNeff_bbn   = NaN;
    obs.dNeff_cmb   = NaN;
    obs.neff        = [];
    obs.errors.neff = ME.message;
end

% ------------------------------------------------------------------
% DM-nucleon cross section [cm^2] (GeV^-2 from DMnucxsec, * hbarc^2).  Use the
% FULL angular-integrated sigma (full_sigma_name maps a heavy* relic model to its
% 'scalar*' full sibling) so the recorded observable is the accurate full-range
% value and matches the Julia/scan side, regardless of the relic approximation.
% ------------------------------------------------------------------
try
    obs.sigma_chiN = DMnucxsec(mx, mphi, gx, g, full_sigma_name(interaction)) * hbarc^2;
catch ME
    obs.sigma_chiN   = NaN;
    obs.errors.sigma = ME.message;
end

% ------------------------------------------------------------------
% EFT-validity flags.  For axialVector, two independent boundaries
% (mirrors compare_limits_g.m):
%   (1) unitarity ceiling, gx-dependent: gx <= min(axial_gx_ceiling,
%       sqrt(pi/2) mphi/mx).  axial_gx_ceiling keeps each meson BW width below
%       its mass and exists only while j -> chi chi is open (Inf for
%       mx > m_eta'/2).  The second term is perturbative unitarity of
%       chi chi -> chi chi via the longitudinal A', lambda_a mx/m_A' < sqrt(pi/2)
%       (Kahlhoefer, Schmidt-Hoberg, Schwetz & Vogl, JHEP 02 (2016) 016); it
%       bounds the 1/m_A'^4 of the J = 0 chi chi -> 3pi channel.
%   (2) small-mixing breakdown -- sin(theta)=max_j|m_{jA'}|/mphi >= s_mix;
%       gx-INdependent (set by g and mphi only).  axial_mixing_floor returns
%       mxmix = max_j|m_{jA'}|, from which sin(theta) = mxmix/mphi.
% The SP/PP families rest on the SAME kind of expansion -- the pseudoscalar
% mediator's mixing with pi0/eta/eta' -- so they get the same treatment, with
% sin_mix = max_j |theta_j| from sppp_mixing.m (see that file for the draft's
% statement of the condition).  There is no unitarity ceiling on that side.
% S(x)S and P(x)S go through the chiral loop factor D, not a mixing angle, so
% they really are trivially valid.
% ------------------------------------------------------------------
obs.eft = struct('valid', true, 'unitarity_ok', true, 'gx_ceiling', Inf, ...
                 'small_mixing_ok', true, 'sin_mix', NaN, 's_mix', s_mix);
if ismember(interaction, {'heavyScalarSP','lightScalarSP','scalarSP', ...
                          'scalarPP','heavyScalarPP','lightScalarPP'})
    try
        [~, ~, theta_max]       = sppp_mixing(mphi, g);
        obs.eft.sin_mix         = theta_max;
        obs.eft.small_mixing_ok = theta_max < s_mix;
        obs.eft.valid           = obs.eft.small_mixing_ok;
    catch ME
        obs.eft.valid           = false;
        obs.eft.small_mixing_ok = false;
        obs.errors.eft          = ME.message;
    end
elseif strcmp(interaction, 'axialVector')
    try
        gx_ceil = axial_gx_ceiling(mphi, mx, g);
        if ~isnan(gx_ceil)               % keep NaN (invalid): MATLAB's min drops it
            gx_ceil = min(gx_ceil, sqrt(pi/2)*mphi/mx);
        end
        [~, mxmix] = axial_mixing_floor(mphi, mx, g, s_mix);
        sin_mix = mxmix / mphi;
        obs.eft.gx_ceiling      = gx_ceil;
        obs.eft.unitarity_ok    = gx <= gx_ceil;     % NaN ceiling -> false (invalid)
        obs.eft.sin_mix         = sin_mix;
        obs.eft.small_mixing_ok = sin_mix < s_mix;
        obs.eft.valid           = obs.eft.unitarity_ok && obs.eft.small_mixing_ok;
    catch ME
        % axial_couplings hard-errors when the mediator is too light (the
        % unitarity boundary itself): the point is EFT-invalid, not unknown.
        obs.eft.valid        = false;
        obs.eft.unitarity_ok = false;
        obs.errors.eft       = ME.message;
    end
end

end
