% test_threepi.m  Checks on the 3pi and pi pi phi channels.
% Run from matlab/.  Checks the default is off, then turns both channels on.

clear threepi_enabled pipiphi_enabled
setenv('LOOP_3PI', '');  setenv('LOOP_PIPIPHI', '');
nfail = 0;
off = ~threepi_enabled() && ~pipiphi_enabled() ...
      && threepi_weight(1, [2 -1 -1]) == 0 && pipiphi_weight(1, 1e-6, [2 -1 -1]) == 0 ...
      && axial3pi_weight(1, [1 -1 0]) == 0;
nfail = nfail + ~off;
fprintf('%-58s %s\n', 'both channels off by default', passfail(off));
setenv('LOOP_3PI', '1');  setenv('LOOP_PIPIPHI', '1');
clear threepi_enabled pipiphi_enabled

check = @(name, ok) fprintf('%-58s %s\n', name, passfail(ok));

g  = [1 0 -1];                       % g_u-g_d = 1, g_u+g_d-2g_s = 3: eta present
g0 = [2 -1 -1]/3;                    % g_u-g_d = 1, g_u+g_d-2g_s = 1

% 1. eta decoupled: W against the closed-form pi0/phi amplitudes, adaptive quad.
mpi0 = 0.13498; mpc = 0.13957; B0 = mpc^2/(0.0022 + 0.00469); f = 0.0924;
Yf = @(E, rs, m1, m2) sqrt(max((E.^2 - m1^2).*(1 - 4*m2^2./(rs^2 + m1^2 - 2*E*rs)), 0));
gnoeta = [0.5 -0.5 0];               % D = 1, g_u+g_d-2g_s = 0: no eta, no eta'
for rs = [0.45 0.8 1.6]
    s = rs^2;
    pm = integral(@(E) Yf(E, rs, mpi0, mpc).*((B0/f)*(s + mpi0^2 - 2*E*rs - mpi0^2)/(s - mpi0^2)).^2, ...
                  mpi0, (s + mpi0^2 - 4*mpc^2)/(2*rs), 'RelTol', 1e-10);
    z  = integral(@(E) Yf(E, rs, mpi0, mpi0), mpi0, (s - 3*mpi0^2)/(2*rs), 'RelTol', 1e-10) ...
         * ((B0/f)*s/(s - mpi0^2))^2 / 6;
    W  = threepi_weight(s, gnoeta);
    ok = abs(W/(pm + z) - 1) < 1e-6;
    nfail = nfail + ~ok;
    check(sprintf('W(s) vs adaptive quad at rs = %.2f (ratio %.8f)', rs, W/(pm + z)), ok);
end

% 2. Below threshold nothing; between the two thresholds only 3pi0.
ok = threepi_weight(0.15, g0) == 0;  nfail = nfail + ~ok;
check('zero below 3 mpi0', ok);
sb = (0.5*(3*mpi0 + 2*mpc + mpi0))^2;
ok = threepi_weight(sb, gnoeta) > 0 && sb < (2*mpc + mpi0)^2;  nfail = nfail + ~ok;
check('3pi0 open between the two thresholds', ok);

% 3. Vectorized = scalar.
sv = linspace(0.2, 3, 7);
wv = threepi_weight(sv, g);  ws = arrayfun(@(x) threepi_weight(x, g), sv);
ok = max(abs(wv./ws - 1)) < 1e-13 && isequal(size(wv), size(sv));  nfail = nfail + ~ok;
check('vectorized matches scalar, shape kept', ok);

% 4. Width: mphi = 0 is finite, below threshold nothing, above it dominant.
Gz = phi_partial_widths(0, 0.01, 1e-3, g0, 'lightScalarPP');
Gl = phi_partial_widths(0.3, 0.2, 1e-3, g0, 'heavyScalarPP');
Gh = phi_partial_widths(1, 0.6, 1e-3, g0, 'heavyScalarPP');
ok = isfinite(Gz.total) && Gl.pipi == 0 && Gh.pipi > 100*Gh.gg;  nfail = nfail + ~ok;
check(sprintf('widths: mphi=0 finite, mphi=1: 3pi/gg = %.3g', Gh.pipi/Gh.gg), ok);

% 5. Cross section is in the SM channels, with the S(x)P suppression.
mx = 0.25;  s = 0.8^2;
PP = xxggxsec(mx, 1, 1e-3, 1e-4*g0, s, 'heavyScalarPP', 'had');
SP = xxggxsec(mx, 1, 1e-3, 1e-4*g0, s, 'heavyScalarSP', 'had');
ok = PP > 0 && abs(SP/PP/((s - 4*mx^2)/s) - 1) < 1e-3;  nfail = nfail + ~ok;
check(sprintf('xxggxsec: SP/PP = (s-4mx^2)/s (ratio %.5f)', SP/PP/((s - 4*mx^2)/s)), ok);

% 6. pi pi phi: zero below threshold, finite at mphi = 0, grows as g^4.
ok = pipiphi_weight(0.06, 1e-6, g0) == 0;  nfail = nfail + ~ok;
check('pi pi phi zero below 2 mpi0 + mphi', ok);
ok = isfinite(pipiphi_weight(0.2, 0, g0));  nfail = nfail + ~ok;
check('pi pi phi finite at mphi = 0', ok);
r = pipiphi_weight(0.25, 1e-6, 2*g0)/pipiphi_weight(0.25, 1e-6, g0);
ok = abs(r/16 - 1) < 1e-12;  nfail = nfail + ~ok;
check(sprintf('pi pi phi scales as g^4 (ratio %.6f)', r), ok);

% 7. Axial vector: WT, WL against the note's s_{+-} integrals (adaptive quad),
%    a different route from the E_0 moments in axial3pi_weight.
lam = @(a, b, c) a.^2 + b.^2 + c.^2 - 2*(a.*b + b.*c + c.*a);
ga  = [0.5 -0.5 0];                  % a_u - a_d = 1
mjA = axial_mixing(ga);  Da2 = (mjA(1)/f)^2;   % current normalization m_{pi0 A'}/F
for rs = [0.45 0.8 1.6]
    s   = rs^2;  hi = (rs - mpi0)^2;
    bpm = @(x) sqrt(1 - 4*mpc^2./x);  b00 = @(x) sqrt(1 - 4*mpi0^2./x);
    T   = Da2*integral(@(x) bpm(x).*lam(s, x, mpi0^2).^1.5, 4*mpc^2, hi, 'RelTol', 1e-10) / (4*f^2*s^2);
    L   = Da2*mpi0^4/(4*f^2*s*(s - mpi0^2)^2) * ( ...
          integral(@(x) (x - mpi0^2).^2.*bpm(x).*sqrt(lam(s, x, mpi0^2)), 4*mpc^2, hi, 'RelTol', 1e-10) ...
          + s^2/6*integral(@(x) b00(x).*sqrt(lam(s, x, mpi0^2)), 4*mpi0^2, hi, 'RelTol', 1e-10));
    [WT, WL] = axial3pi_weight(s, ga);
    ok = abs(WT/T - 1) < 1e-6 && abs(WL/L - 1) < 1e-6;  nfail = nfail + ~ok;
    check(sprintf('axial WT, WL vs s+- quad at rs = %.2f (%.8f %.8f)', rs, WT/T, WL/L), ok);
end

% 8. PCAC: the longitudinal part is the pseudoscalar pi0 path times (mpi^2/B0)^2.
sv = linspace(0.2, 3, 7);
[~, WL] = axial3pi_weight(sv, ga);
r  = WL ./ (Da2*mpi0^4/B0^2 * threepi_weight(sv, gnoeta));
ok = max(abs(r(sv > 9*mpi0^2) - 1)) < 1e-12;  nfail = nfail + ~ok;
check('axial WL = (mpi^2/B0)^2 x pseudoscalar W (PCAC)', ok);

% 9. Only m_{pi0 A'} enters (isoscalar a: the L5 isospin-breaking remnant);
%    no transverse 3pi0 between the thresholds.
[WT, WL]   = axial3pi_weight(1, [1 1 1]);
[WT1, WL1] = axial3pi_weight(1, ga);
m111 = axial_mixing([1 1 1]);  r2 = (m111(1)/mjA(1))^2;
[WTb, WLb] = axial3pi_weight(sb, ga);
ok = abs(WT/WT1/r2 - 1) < 1e-12 && abs(WL/WL1/r2 - 1) < 1e-12 && r2 < 1e-4 ...
     && WTb == 0 && WLb > 0;  nfail = nfail + ~ok;
check(sprintf('axial: scales as m_{pi0 A''}^2 (isoscalar/isovector %.2g); 3pi0 J = 0 only', r2), ok);

% 10. On the A' pole the J = 1 piece is the spin-1 Breit-Wigner peak
%     16 pi/(s beta^2) (3/4) B_xx B_3pi, with Gamma(A' -> chi chi) and
%     Gamma(A' -> 3pi) from phi_partial_widths: fixes the relative normalization.
mA = 0.8;  mx = 0.2;  la = 1e-2;  ga2 = 0.3*ga;
GA = phi_partial_widths(mA, mx, la, ga2, 'axialVector');
s  = mA^2;  beta = sqrt(1 - 4*mx^2/s);
sJ1 = axial3pi_sigma(s, mx, mA, la, GA.total, ga2) - axial3pi_sigma(s, mx, mA, la, Inf, ga2);
sBW = 16*pi/(s*beta^2) * 3/4 * GA.xx*GA.pipi/GA.total^2;
ok = GA.pipi > 0 && abs(sJ1/sBW - 1) < 1e-10;  nfail = nfail + ~ok;
check(sprintf('axial J = 1 peak = spin-1 BW (ratio %.10f)', sJ1/sBW), ok);

% 11. Width finite at mA = 0, zero below 3 mpi0; xxggxsec carries the channel.
G0 = phi_partial_widths(0, 0.01, 1e-3, ga, 'axialVector');
Gb = phi_partial_widths(0.3, 0.1, 1e-3, ga, 'axialVector');
ok = isfinite(G0.total) && G0.pipi == 0 && Gb.pipi == 0;  nfail = nfail + ~ok;
check('axial width: mA = 0 finite, zero below 3 mpi0', ok);
mx = 0.3;  mA = 1;  s = [0.5 0.9 1.5];  ga3 = 1e-2*[2 -1 -1];
on  = xxggxsec(mx, mA, 1e-2, ga3, s, 'axialVector', 'had');
GA  = phi_partial_widths(mA, mx, 1e-2, ga3, 'axialVector').total;
s3  = axial3pi_sigma(s, mx, mA, 1e-2, GA, ga3);
setenv('LOOP_3PI', '');  clear threepi_enabled
off = xxggxsec(mx, mA, 1e-2, ga3, s, 'axialVector', 'had');
setenv('LOOP_3PI', '1');  clear threepi_enabled
ok = max(abs((on - off)./s3 - 1)) < 1e-12 && all(s3 > 0);  nfail = nfail + ~ok;
check('xxggxsec axialVector: on - off = axial3pi_sigma', ok);

% 12. pi+pi-A': Dalitz weights against a brute-force 2D integral of the
%     polarization-summed current (P = -g + k k/mA^2, CM frame), not the X moments.
mA = 0.05;  gq = [1.5 -0.5 0.2];  mjA2 = axial_mixing(gq);
c4 = 4*(mjA2(1)/0.0924)^2*(gq(1) - gq(2))^2;
ok = true;
for s = [0.2 0.6 2]
    rs = sqrt(s);
    EAhi = (s + mA^2 - 4*mpc^2)/(2*rs);
    fT = @(EA, Ep) c4*(3 + (EA.^2 - mA^2)/mA^2);
    fL = @(EA, Ep) c4*mpi0^4/(s - mpi0^2)^2*(-s + s*EA.^2/mA^2);
    lo = @(EA) min_lim(EA, s, mA, mpc);  hi = @(EA) max_lim(EA, s, mA, mpc);
    bT = integral2(fT, mA, EAhi, lo, hi, 'RelTol', 1e-10, 'AbsTol', 0);
    bL = integral2(fL, mA, EAhi, lo, hi, 'RelTol', 1e-10, 'AbsTol', 0);
    [WT, WL] = pipiA_weight(s, mA, gq);
    ok = ok && abs(WT/bT - 1) < 1e-6 && abs(WL/bL - 1) < 1e-6;
end
nfail = nfail + ~ok;
check(sprintf('pi pi A'': WT, WL = brute-force Dalitz (last %.2e, %.2e)', WT/bT - 1, WL/bL - 1), ok);

% 13. pi+pi-A': zero below threshold and for a_u = a_d; scales as
%     (a_u - a_d)^2 m_{pi0 A'}^2; xxggxsec carries it (on - off) with LOOP_3PI off.
[a1, b1] = pipiA_weight((2*mpc + mA)^2*0.999, mA, gq);
[a2, b2] = pipiA_weight(0.6, mA, [1 1 0.3]);
[a3, b3] = pipiA_weight(0.6, mA, 2*gq);  [a4, b4] = pipiA_weight(0.6, mA, gq);
ok = a1 == 0 && b1 == 0 && a2 == 0 && b2 == 0 && abs(a3/a4/16 - 1) < 1e-12 && abs(b3/b4/16 - 1) < 1e-12;
setenv('LOOP_3PI', '');  clear threepi_enabled
mx = 0.3;  mA = 0.1;  s = [0.5 0.9 1.5];  ga3 = 1e-2*[2 -1 -1];
on  = xxggxsec(mx, mA, 1e-2, ga3, s, 'axialVector', 'sm');
sa  = axial3pi_sigma(s, mx, mA, 1e-2, phi_partial_widths(mA, mx, 1e-2, ga3, 'axialVector').total, ga3);
setenv('LOOP_PIPIPHI', '');  clear pipiphi_enabled
off = xxggxsec(mx, mA, 1e-2, ga3, s, 'axialVector', 'sm');
[z1, z2] = pipiA_weight(0.6, mA, gq);
setenv('LOOP_3PI', '1');  setenv('LOOP_PIPIPHI', '1');  clear threepi_enabled pipiphi_enabled
ok = ok && max(abs((on - off)./sa - 1)) < 1e-12 && all(sa > 0) && z1 == 0 && z2 == 0;
nfail = nfail + ~ok;
check('pi pi A'': threshold, (a_u - a_d)^2 m_{pi0A''}^2 scaling, switch, xxggxsec', ok);

setenv('LOOP_3PI', '');  setenv('LOOP_PIPIPHI', '');
clear threepi_enabled pipiphi_enabled
fprintf('\ntest_threepi: %d failure(s)\n', nfail);

% E_+ range at fixed E_A (pi+pi- recoiling against the A')
function v = min_lim(EA, s, mA, m)
[v, ~] = ep_lims(EA, s, mA, m);
end
function v = max_lim(EA, s, mA, m)
[~, v] = ep_lims(EA, s, mA, m);
end
function [lo, hi] = ep_lims(EA, s, mA, m)
rs = sqrt(s);  w = s + mA^2 - 2*EA*rs;  kA = sqrt(max(EA.^2 - mA^2, 0));
b  = sqrt(max(1 - 4*m^2./w, 0));
lo = ((rs - EA) - kA.*b)/2;  hi = ((rs - EA) + kA.*b)/2;
end

function s = passfail(ok)
if ok, s = 'ok'; else, s = 'FAIL'; end
end
