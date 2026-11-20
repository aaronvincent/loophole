function eta_max = axial_eta_ceiling(mphi, mx, gx, g)
% axial_eta_ceiling  EFT-validity ceiling on the eta-scan normalization for the
% axialVector interaction (companion of axial_gx_ceiling, for compare_limits_g).
%
%   eta_max = axial_eta_ceiling(mphi, mx, gx, g)
%
% In the eta-scan gx is FIXED and eta scales the quark couplings (g_eff = eta*g),
% so BOTH axial_couplings.m validity gates move with eta and the gx-solving
% axial_gx_ceiling cannot be reused (at eta=1 the small-mixing gate is already
% tripped for an ultralight m_A', so its sub-floor probe never clears and it
% returns NaN).  The two gates are:
%   * mediatorTooLight -- meson->chichi width reaches its mass; the width
%     ~ lambda_{jchi}^2 ~ (gx*eta)^2, so it grows with eta.
%   * mixingTooLarge   -- the A'-meson mixing mass m_{jA'} ~ eta reaches m_A',
%     breaking the sin(theta) ~ m_{jA'}/m_A' expansion; also grows with eta.
% Both are monotonic in eta, so the valid set is eta < eta_max and we bisect on
% log10(eta) for the boundary.  Returns Inf if eta=1e12 is still valid (no
% ceiling, e.g. heavy m_A'); NaN if even eta=1e-30 is invalid.

valid = @(le) local_valid(mphi, mx, gx, 10^le, g);

lo = -30; hi = 12;
if valid(hi),  eta_max = Inf;  return; end   % no ceiling in range
if ~valid(lo), eta_max = NaN;  return; end    % invalid even at vanishing eta

for it = 1:60
    mid = 0.5 * (lo + hi);
    if valid(mid), lo = mid; else, hi = mid; end
end
eta_max = 10^lo;
end

function tf = local_valid(mphi, mx, gx, eta, g)
tf = true;
try
    axial_couplings(mphi, mx, gx, eta * g);
catch ME
    if any(strcmp(ME.identifier, {'axial_couplings:mediatorTooLight', ...
                                  'axial_couplings:mixingTooLarge'}))
        tf = false;
    else
        rethrow(ME);
    end
end
end
