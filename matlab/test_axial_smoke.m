% Smoke test for the axialVector interaction wiring.
mphi = 0.5; mx = 0.1; gx = 1; g = [1 1 1]/3; intr = 'axialVector';

fprintf('--- axial_couplings ---\n');
[mj, lam_chi, lam_gg, Gam_j] = axial_couplings(mphi, mx, gx, g);
fprintf('mj      = [%g %g %g]\n', mj);
fprintf('lam_chi = [%g %g %g]\n', lam_chi);
fprintf('lam_gg  = [%g %g %g]\n', lam_gg);
fprintf('Gam_j   = [%g %g %g]\n', Gam_j);
assert(all(isfinite([mj lam_chi lam_gg Gam_j])), 'couplings not finite');

fprintf('--- phi_partial_widths ---\n');
G = phi_partial_widths(mphi, mx, gx, g, intr);
fprintf('Gamma_A.total = %g (xx=%g)\n', G.total, G.xx);
assert(isfinite(G.total) && G.total>0, 'width bad');

% vector of s spanning all thresholds: below/at meson poles, above 4mphi^2
s = logspace(log10(4*mx^2*1.001), 3, 400).';

fprintf('--- ggxxxsec (gg->xx) ---\n');
sgg = ggxxxsec(mx, mphi, gx, g, s, intr);
fprintf('finite=%d  nonneg=%d  max=%g\n', all(isfinite(sgg)), all(real(sgg)>=-1e-30), max(real(sgg)));
assert(all(isfinite(sgg)), 'ggxxxsec not finite');

fprintf('--- xxggxsec (xx->all) ---\n');
% mx>mphi case so A''A'' + meson-pair channels open
mx2 = 0.6; mphi2 = 0.5;
s2 = logspace(log10(4*mx2^2*1.001), 3, 400).';
sall = xxggxsec(mx2, mphi2, gx, g, s2, intr);
fprintf('finite=%d  real-part max=%g  max|imag|=%g\n', all(isfinite(sall)), max(real(sall)), max(abs(imag(sall))));
assert(all(isfinite(sall)), 'xxggxsec not finite');

fprintf('--- DMnucxsec ---\n');
sN_h = DMnucxsec(mx, 1.0, gx, g, intr);    % heavy
sN_l = DMnucxsec(mx, 1e-7, gx, g, intr);   % light
fprintf('heavy=%g  light=%g  ratio=%g (expect ~3 vs 1 scaled by mphi^-4)\n', sN_h, sN_l, sN_h/sN_l);
assert(isfinite(sN_h) && isfinite(sN_l), 'DMnucxsec bad');

fprintf('--- svavg gg / all ---\n');
sv_gg  = svavg(mx, mphi, gx, g, 0.05, intr, 'gg');
sv_all = svavg(mx, mphi, gx, g, 0.05, intr, 'all');
fprintf('sv_gg=%g  sv_all=%g\n', sv_gg, sv_all);
assert(isfinite(sv_gg) && isfinite(sv_all) && sv_gg>=0 && sv_all>=0, 'svavg bad');

fprintf('--- ggxxrate ---\n');
Tg = [0.005 0.01 0.02];
R = ggxxrate(Tg, mx, mphi, gx, g, intr);
fprintf('rate = [%g %g %g]\n', R);
assert(all(isfinite(R)) && all(R>=0), 'ggxxrate bad');

fprintf('\nALL AXIAL SMOKE TESTS PASSED\n');
