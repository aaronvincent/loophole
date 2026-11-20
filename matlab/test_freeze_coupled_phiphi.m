function test_freeze_coupled_phiphi()
% test_freeze_coupled_phiphi  Checks for freeze_coupled's chi chibar <-> phi phi
% and fast-phi treatment, built for the regimes they change.
%
%   1. xxggxsec: 'all' = 'sm' + 'dark' to round-off, every family.
%   2. sv_darkpair_grid: forward <sv> = svavg_phiphi (independent adaptive
%      quadrature) to < 1%, S⊗S / P⊗S / axialVector, including mx < mphi; the
%      reverse rate is O(1) deep in the tail where mphi > mx (phi phi -> chi chibar
%      is exothermic), where the forward one is exponentially small.
%   3. Thermal limit: phi held in equilibrium by its SM link (light phi, large eta)
%      => freeze_coupled = freeze_loop with production 'all' (phi phi -> chi chibar
%      at equilibrium phi), on both relic branches.
%   4. Fast-phi elimination = full 2x2 integration (LOOP_FC_FULL2D=1) at a point
%      where both converge; the dark-exchange shortcut = full 2x2 (LOOP_FC_NOEXCH=1)
%      in a closed dark sector, and finite where the 2x2 system cannot start.
%   5. No step in Omega h^2(eta) across the Y_phi initial-condition and fast-phi
%      switches: finite everywhere, bounded log-slope.
% Run from matlab/:  matlab -batch "test_freeze_coupled_phiphi"   (~2 min)

warning('off', 'all');
np = 0;  nf = 0;
    function check(name, ok, val)
        if ok, np = np + 1; tag = 'PASS'; else, nf = nf + 1; tag = 'FAIL'; end
        fprintf('  %s  %-62s %s\n', tag, name, val);
    end

% --- 1 -------------------------------------------------------------------
s = logspace(log10(4*0.05^2*1.0001), 0, 60);
for ic = {'heavyScalar','heavyScalarPS','heavyScalarSP','heavyScalarPP','axialVector'}
    a = xxggxsec(0.05, 0.03, 0.5, 1e-2*[2 -1 -1], s, ic{1}, 'all');
    b = xxggxsec(0.05, 0.03, 0.5, 1e-2*[2 -1 -1], s, ic{1}, 'sm');
    d = xxggxsec(0.05, 0.03, 0.5, 1e-2*[2 -1 -1], s, ic{1}, 'dark');
    e = max(abs(a - b - d) ./ a);
    check(sprintf('all = sm + dark, %s', ic{1}), e < 1e-12, sprintf('%.1e', e));
end

% --- 2 -------------------------------------------------------------------
M = [1 3; 5 8; 20 5]*1e-3;
for ic = {'heavyScalar','heavyScalarPS','axialVector'}
    gphi = 1 + 2*strcmp(ic{1}, 'axialVector');
    for k = 1:3
        T = M(k,1) ./ [1 5 20 60];
        [lf, lr] = sv_darkpair_grid(M(k,1), M(k,2), 0.5, 1e-3*[2 -1 -1], T, ic{1}, gphi);
        ref = arrayfun(@(t) svavg_phiphi(M(k,1), M(k,2), 0.5, 1e-3*[2 -1 -1], t, ic{1}), T);
        e = max(abs(exp(lf)./ref - 1));
        check(sprintf('fwd = svavg_phiphi, %s mx=%g mphi=%g', ic{1}, M(k,1), M(k,2)), e < 1e-2, ...
              sprintf('%.1e', e));
        if M(k,2) > M(k,1)
            check(sprintf('rev O(1) at mphi/T = %g, %s', M(k,2)/T(end), ic{1}), ...
                  all(isfinite(lr)) && abs(lr(end) - lr(1)) < 10, mat2str(lr, 3));
        end
    end
end

% --- 3 -------------------------------------------------------------------
% lightScalarPS, mphi = 1e-9: phi is in the bath (Primakoff, pi gamma) at these eta.
for pt = {[0.0299 3.64], [0.03 2e-2]}
    mx = pt{1}(1);  eta = pt{1}(2);  g = eta*[2 -1 -1];
    fc = freeze_coupled(mx, 1e-9, 0.01, g, 'lightScalarPS', 0.01);
    fl = freeze_loop(mx, 1e-9, 0.01, g, 'lightScalarPS', 0.01, false, 'all');
    check(sprintf('thermal phi: coupled = loop(all), mx=%g eta=%g', mx, eta), ...
          abs(fc/fl - 1) < 2e-2, sprintf('%.4g / %.4g', fc, fl));
end

% --- 4 -------------------------------------------------------------------
g = 1e-4*[2 -1 -1];
setenv('LOOP_FC_FULL2D', '');   a = freeze_coupled(0.01, 0.1, 0.01, g, 'heavyScalarPS', 0.01);
setenv('LOOP_FC_FULL2D', '1');  b = freeze_coupled(0.01, 0.1, 0.01, g, 'heavyScalarPS', 0.01);
setenv('LOOP_FC_FULL2D', '');
check('fast-phi elimination = full 2x2 (heavyScalarPS mx=10 MeV mphi=0.1)', ...
      abs(a/b - 1) < 1e-3, sprintf('%.6g / %.6g', a, b));

% --- 4b: dark-exchange shortcut (closed sector, axialVector Gamma(A') = 0) --
for pt = {[0.006 0.01 0.03 1e-5], [0.02 0.01 0.03 1e-3]}
    q = pt{1};  g = q(4)*[2 -1 -1];
    setenv('LOOP_FC_NOEXCH', '');   a = freeze_coupled(q(1), q(2), q(3), g, 'axialVector', 0.01);
    setenv('LOOP_FC_NOEXCH', '1');  b = freeze_coupled(q(1), q(2), q(3), g, 'axialVector', 0.01);
    setenv('LOOP_FC_NOEXCH', '');
    check(sprintf('exchange shortcut = full 2x2, axialVector mx=%g mphi=%g gx=%g', q(1), q(2), q(3)), ...
          abs(a/b - 1) < 2e-3, sprintf('%.6g / %.6g', a, b));
end
v = freeze_coupled(0.006, 0.01, 1, 1e-3*[2 -1 -1], 'axialVector', 0.01);
check('closed sector at gx = 1 is finite (2x2 alone cannot start)', isfinite(v) && v > 0, sprintf('%.4g', v));

% --- 5 -------------------------------------------------------------------
le = linspace(-9, -1, 49);
for pt = {{'heavyScalarPS', 0.01, 0.1, 0.01}, {'heavyScalar', 0.001, 0.01, 1}, {'lightScalarSP', 0.06, 1e-9, 0.01}}
    [ic, mx, mp, gx] = pt{1}{:};
    v = arrayfun(@(l) freeze_coupled(mx, mp, gx, 10^l*[2 -1 -1], ic, 0.01), le);
    d = diff(log(v)) ./ diff(le*log(10));
    check(sprintf('Oh2(eta) continuous, %s mx=%g mphi=%g gx=%g', ic, mx, mp, gx), ...
          all(isfinite(v)) && max(abs(diff(d))) < 1.5, ...
          sprintf('max slope change %.2f', max(abs(diff(d)))));
end

fprintf('\ntest_freeze_coupled_phiphi: %d passed, %d failed\n', np, nf);
if nf > 0, error('test_freeze_coupled_phiphi: %d failed', nf); end
end
