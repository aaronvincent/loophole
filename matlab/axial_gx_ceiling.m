function gx_max = axial_gx_ceiling(mphi, mx, g)
% axial_gx_ceiling  Unitarity / EFT-validity ceiling on the axialVector
% coupling gx (= lambda_a).
%
%   gx_max = axial_gx_ceiling(mphi, mx, g)
%
% Returns the largest gx for which the axial-vector meson-mixing EFT is
% self-consistent at (mphi = m_A', mx = m_chi, g = [a_u a_d a_s]).  It is the
% companion of the hard error thrown by axial_couplings.m: that routine refuses
% once any meson Breit-Wigner width reaches its mass, max_j Gamma_j/m_j >= 1;
% here we instead solve max_j Gamma_j/m_j = 1 FOR gx and return the boundary.
%
% Physics: the meson->chi chi width is
%   Gamma(j->chi chi) = 2 lambda_{j chi}^2 mx^2 sqrt(mj^2/4 - mx^2)/pi,
% with the longitudinally enhanced mixing coupling lambda_{j chi} = gx m_{jA'}/mphi^2
% (m_{jA'} the gx-independent A'-meson mixing mass).  Requiring
% Gamma_j = Gamma_SM,j + Gamma(j->chi chi) < m_j gives, per meson,
%   gx_max,j = (mphi^2/m_{jA'}) sqrt( (m_j - Gamma_SM,j) pi / (2 mx^2 sqrt(mj^2/4 - mx^2)) ),
% and the binding (most restrictive) meson sets gx_max = min_j gx_max,j.  Mesons
% with the chi chi channel closed (mj <= 2mx) impose no constraint (gx_max,j = Inf);
% if all three are closed the ceiling is Inf.  Scalings: gx_max ~ mphi^2 and, for
% mx << mj/2, gx_max ~ 1/mx.
%
% axial_couplings.m is the single source of truth for the mixing masses and SM
% widths, so rather than duplicate the chiral expressions we probe it at a small
% gx (below its validity floor) and recover the gx-independent m_{jA'} and
% Gamma_SM,j from the returned arrays.

% Probe axial_couplings at a gx small enough to clear its validity floor.
% lambda_{j chi} is linear in gx and the mixing masses are gx-independent, so
% any sub-floor probe recovers them exactly; Gamma(j->chi chi) ~ gx^2 -> 0.
gx_probe = 1;
mj = []; lam_chi = []; Gam_j = [];
ok = false;
for it = 1:200
    try
        [mj, lam_chi, ~, Gam_j] = axial_couplings(mphi, mx, gx_probe, g);
        ok = true;
        break;
    catch ME
        if strcmp(ME.identifier, 'axial_couplings:mediatorTooLight')
            gx_probe = gx_probe * 1e-6;   % step well below the floor
            continue;
        else
            rethrow(ME);
        end
    end
end
if ~ok
    gx_max = NaN;   % floor unreachable even at vanishing gx (absurdly light m_A')
    return;
end

% Recover the gx-independent inputs from the probe.  Use |m_{jA'}|: the width
% depends only on lambda_{j chi}^2, so the sign of the mixing mass is irrelevant
% (for g = [1 1 1]/3 the pi0/eta mixing masses are tiny and sign-flipped by the
% L5 quark-mass term; without the abs they produce a spurious negative ceiling).
m_jAp  = abs(lam_chi) * mphi^2 / gx_probe;                     % |lambda_{j chi}| = gx |m_{jA'}|/mphi^2
kin    = mx^2 .* sqrt(max(0, mj.^2/4 - mx^2)) / pi;        % Gamma(j->xx) = (gx m_{jA'}/mphi^2)^2 * kin
Gam_xx_probe = (lam_chi).^2 .* kin;
Gam_SM = Gam_j - Gam_xx_probe;

% Per-meson ceiling; chi chi-closed mesons (kin == 0) impose no bound.
open = kin > 0;
gx_max_j = inf(size(mj));
headroom = max(0, mj - Gam_SM);                                % m_j - Gamma_SM,j (~ m_j)
gx_max_j(open) = (mphi^2 ./ m_jAp(open)) .* ...
                 sqrt(headroom(open) ./ kin(open));
gx_max = min(gx_max_j);
end
