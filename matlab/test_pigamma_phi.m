% test_pigamma_phi  Checks for pigamma_phi_xsec.m / pigamma_phi_rate.m.
%
%   1. sigma against an independent Python angular integration of the
%      explicit-4-vector |M|^2 (mu = 1, i.e. g_u + g_d = 1/B0).
%   2. Closed form vs threshold quadrature agree across the switch point.
%   3. sigma and the rate scale as (g_u + g_d)^2; g_s does not enter.
%   4. Massless-phi rate is continuous with the massive branch (z = 1e-4 switch).
B0 = 0.13957^2/(0.0022 + 0.00469);
gmu1 = [1/B0 0 0];                       % mu = 1
nfail = 0;
chk = @(name, ok) fprintf('%-58s %s\n', name, ternary(ok));

% 1. reference values: scratch Python check (dsigma/dcos from explicit
%    J = p/(p.k) - p'/(p'.k), quad over cos theta), alpha = 1/137.036.
ref = [0.03 1e-5 1.0213973319804588;
       0.05 0.05 0.29741317317170285;
       0.2  0.1  0.07901347176364891];
for r = 1:size(ref, 1)
    v = pigamma_phi_xsec(ref(r,1), ref(r,2), gmu1);
    ok = abs(v/ref(r,3) - 1) < 1e-9;  nfail = nfail + ~ok;
    chk(sprintf('sigma(s=%g, mphi=%g) = %.10g vs %.10g', ref(r,1), ref(r,2), v, ref(r,3)), ok);
end

% 2. continuity across THR_SWITCH = 1e-2 (closed form above, quadrature below)
for mphi = [1e-9 1e-3 0.1]
    sth = (0.13957 + mphi)^2;
    v = pigamma_phi_xsec(sth*(1 + 1e-2*[1-1e-9, 1+1e-9]), mphi, gmu1);
    ok = abs(v(1)/v(2) - 1) < 1e-6;  nfail = nfail + ~ok;
    chk(sprintf('switch continuity, mphi = %g: rel diff %.1e', mphi, abs(v(1)/v(2)-1)), ok);
end

% 3. coupling scaling
g1 = [2 -1 -1];  g2 = [2 -1 5];  g3 = 3*[2 -1 -1];
s = 0.05;  mphi = 1e-3;  T = 0.01;
ok = abs(pigamma_phi_xsec(s, mphi, g1)/pigamma_phi_xsec(s, mphi, g2) - 1) < 1e-12;
nfail = nfail + ~ok;  chk('sigma independent of g_s', ok);
ok = abs(pigamma_phi_rate(T, mphi, g3)/pigamma_phi_rate(T, mphi, g1) - 9) < 1e-6;
nfail = nfail + ~ok;  chk('rate scales as (g_u + g_d)^2', ok);

% 4. massless branch continuity
Ga = pigamma_phi_rate(T, 0.99e-4*T, g1);  Gb = pigamma_phi_rate(T, 1.01e-4*T, g1);
ok = abs(Ga/Gb - 1) < 1e-3;  nfail = nfail + ~ok;
chk(sprintf('rate continuous at mphi/T = 1e-4: rel diff %.1e', abs(Ga/Gb-1)), ok);

if nfail == 0
    fprintf('test_pigamma_phi: all passed\n');
else
    error('test_pigamma_phi:fail', '%d check(s) failed', nfail);
end

function s = ternary(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
