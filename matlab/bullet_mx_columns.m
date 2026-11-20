function [som, iv, bound] = bullet_mx_columns(mx_vec, mphi, gx, interaction)
% bullet_mx_columns  Bullet-cluster self-interaction bound on an m_chi grid.
%
%   [som, iv, bound] = bullet_mx_columns(mx_vec, mphi, gx, interaction)
%
% sigma_T/m_chi depends on (mx, mphi, gx) and the DM-side vertex only, so at
% fixed (mphi, gx) it excludes whole m_chi intervals, independent of the quark
% couplings.  Evaluated through external_loglikes (bullet only).
%
%   som    sigma_T/m_chi [cm^2/g] at each mx_vec point (NaN: vertex not covered)
%   iv     n x 2 excluded intervals [mx_lo mx_hi]; interior edges are solved
%          for sigma_T/m_chi = bound, edges at the grid ends are left there
%   bound  the sigma/m bound [cm^2/g]
opts = struct('enable_sn1987a', false, 'enable_neutronstar', false, ...
              'enable_stellar', false, 'enable_fifthforce', false, ...
              'enable_kaon', false, 'enable_kphi', false, 'enable_cmb', false, ...
              'enable_integral', false, 'enable_mesoninv', false, 'enable_beamdump', false);
som_at = @(m) bullet_info(m, mphi, gx, interaction, opts);
som = arrayfun(som_at, mx_vec);
L = external_loglikes(mx_vec(1), mphi, gx, [1 1 1], interaction, opts);
bound = L.info.bullet.bound;

ex = isfinite(som) & som > bound;
d  = diff([false, ex(:)', false]);
runs = [find(d == 1); find(d == -1) - 1];
iv = nan(size(runs, 2), 2);
f = @(lm) log(som_at(exp(lm))) - log(bound);
for r = 1:size(runs, 2)
    a = runs(1, r);  b = runs(2, r);
    iv(r, :) = mx_vec([a b]);
    if a > 1
        iv(r, 1) = exp(fzero(f, log(mx_vec([a-1 a]))));
    end
    if b < numel(mx_vec)
        iv(r, 2) = exp(fzero(f, log(mx_vec([b b+1]))));
    end
end
end

function s = bullet_info(m, mphi, gx, interaction, opts)
L = external_loglikes(m, mphi, gx, [1 1 1], interaction, opts);
s = L.info.bullet.sigma_over_m;
end
