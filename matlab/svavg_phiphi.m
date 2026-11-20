function sv = svavg_phiphi(mx, mphi, gx, g, T, interaction, moment)
% svavg_phiphi  Thermally averaged <sigma v> for chi chibar -> phi phi.
%
%   sv = svavg_phiphi(mx, mphi, gx, g, T, interaction)
%
% Used by neff_simple.m for the chi chi <-> phi phi rate
% (chi's link to the plasma through phi, and phi made by chi).  This is the t-channel phi-phi piece only,
% pulled from make_xxsigma(...,'phiphi') so the cross section stays the same
% one xxggxsec.m feeds the relic solvers.
%
% The phi-phi cross section is a smooth t-channel (no s-channel pole), so the
% Breit-Wigner waypoint / narrow-width machinery in svavg.m is not needed --
% a single log-s Gauss-Kronrod pass suffices.  Returns 0 only when the channel
% is unimplemented for this interaction.
%
% Nonzero for mx < mphi: the threshold is s >= 4 mphi^2, and s runs over
% [max(4 mx^2, 4 mphi^2), ...), so the channel costs exp(-2(mphi-mx)/T) there
% rather than closing.  s_lo below sets that range.
%
% Normalization matches svavg.m exactly (Gondolo-Gelmini, scaled Bessels):
%   <sigma v> = [1/(8 mx^4 T K2(mx/T)^2)]
%               * int_{s_lo}^{s_max} sigma(s) (s-4mx^2) sqrt(s) K1(sqrt(s)/T) ds.

% moment (optional, default 'number'): 'energy' returns <sigma v E_tot> in
% GeV^-1 instead of <sigma v> in GeV^-2, weighting each reaction by the lab
% energy of the annihilating pair.  Same kernel swap as svavg.m and ggxxrate.m.
if nargin < 7 || isempty(moment), moment = 'number'; end
if ~ismember(moment, {'number','energy'})
    error('svavg_phiphi:moment', 'moment must be ''number'' or ''energy''');
end

sig = make_xxsigma(mx, mphi, gx, g, interaction, 'phiphi');

K2scaled = besselk(2, mx/T, 1);
if ~isfinite(K2scaled) || K2scaled <= 0
    sv = 0; return;          % mx/T so large the average underflows -> negligible
end

integrand    = @(s)  sig(s).*(s - 4*mx^2).*sqrt(s) ...
                     .* ppkernel(sqrt(s), T, moment) .* exp((2*mx - sqrt(s))/T);
logintegrand = @(ls) integrand(exp(ls)).*exp(ls);

s_lo  = max(4*mx^2, 4*mphi^2);
s_max = max([100*T^2 + s_lo, 100*mx^2, 1e2]);
if s_max <= s_lo
    sv = 0; return;
end

I  = integral(logintegrand, log(s_lo), log(s_max), ...
              'RelTol', 1e-4, 'AbsTol', 1e-30);
sv = I / (8*mx^4*T*K2scaled.^2);
end

% =========================================================================
function K = ppkernel(sq, T, moment)
% See svavg.m's svkernel -- same construction, kept local so svavg_phiphi.m
% stays independently readable.
switch moment
    case 'number'
        K = besselk(1, sq./T, 1);
    case 'energy'
        K = sq .* besselk(0, sq./T, 1) + 2*T .* besselk(1, sq./T, 1);
    otherwise
        error('svavg_phiphi:moment', 'moment must be ''number'' or ''energy''');
end
end
