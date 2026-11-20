% test_primakoff_phi  Checks for primakoff_phi_rate.m (docs/phi_thermalisation_check.md).
%   1. sigma_C closed form vs direct quadrature of the screened |M|^2 over -t
%   2. continuity across the X/mD^2 = 1e-2 series switch
%   3. large-s asymptote (alpha g2/8)(ln(s/mD^2) - 7/4)
%   4. rate exactly ∝ g2; massless rate ~ T (G/T^... ) and mphi-independent at mphi << T
%   5. annihilation piece vs a direct evaluation
npass = 0;  nfail = 0;
alpha = 1/137.036;
sc.sig = @(s, mphi, D, ~) primakoff_phi_xsec(s, mphi, D);
% 1
D = 1e-7;  mphi = 3e-3;
for s = [1.2e-5, 4e-5, 1e-3, 1e-1]
    X = s - mphi^2;
    f = @(x) x.*(s^2 + (X - x).^2)./(x + D).^2;
    ref = alpha*integral(f, 0, X, 'RelTol', 1e-12, 'AbsTol', 0)/(16*s^2);
    got = sc.sig(s, mphi, D, alpha);
    c = abs(got/ref - 1) < 1e-9;  [npass, nfail] = tally(c, npass, nfail);
    fprintf('%s  sigma_C closed form, s = %g: rel err %.2e\n', ternary(c), s, abs(got/ref-1));
end
% 2
D = 1e-4;  mphi = 0;  s1 = 1e-2*D*(1 - 1e-12);  s2 = 1e-2*D*(1 + 1e-12);
X = s1;  f = @(x) x.*(s1^2 + (X - x).^2)./(x + D).^2;
ref = alpha*integral(f, 0, X, 'RelTol', 1e-13, 'AbsTol', 0)/(16*s1^2);
a = sc.sig(s1, mphi, D, alpha);  b = sc.sig(s2, mphi, D, alpha);
c = abs(a/ref - 1) < 1e-5 && abs(b/a - 1) < 1e-4;  [npass, nfail] = tally(c, npass, nfail);
fprintf('%s  series branch vs quadrature %.1e, jump at switch %.1e\n', ternary(c), abs(a/ref-1), abs(b/a-1));
% 3
D = 1e-8;  s = 1e2;
asym = alpha/8*(log(s/D) - 7/4);  got = sc.sig(s, 0, D, alpha);
c = abs(got/asym - 1) < 1e-6;  [npass, nfail] = tally(c, npass, nfail);
fprintf('%s  large-s asymptote: rel err %.2e\n', ternary(c), abs(got/asym-1));
% 4
T = logspace(-4, -2, 5);
G1 = primakoff_phi_rate(T, 1e-9, 1);  G2 = primakoff_phi_rate(T, 1e-9, 7.5);
c = max(abs(G2./(7.5*G1) - 1)) < 1e-14;  [npass, nfail] = tally(c, npass, nfail);
fprintf('%s  rate exactly ∝ g2\n', ternary(c));
G0 = primakoff_phi_rate(T, 0, 1);  Gs = primakoff_phi_rate(T, 1e-7, 1);
c = max(abs(Gs./G0 - 1)) < 1e-3;  [npass, nfail] = tally(c, npass, nfail);
fprintf('%s  mphi << T independent of mphi (max dev %.1e)\n', ternary(c), max(abs(Gs./G0-1)));
% high-T, fe = 1: G/T^3 = (8 T^4/(32 pi^4) ∫ 2 y^4 sigma K1 dy)/(T^3/pi^2); compare with
% the asymptotic-sigma estimate (alpha/8)(ln(s/mD^2)-7/4) at <s> -- just check ballpark
Th = 1;  [~, o] = primakoff_phi_rate(Th, 0, 1);
yy = @(y) 2*y.^4.*besselk(1, y).*(alpha/8).*(log(y.^2/((4*pi*alpha/3))) - 7/4);
est = 8*Th/(32*pi^4)*integral(yy, 1e-6, 80)*Th^4/(Th^3/pi^2);
c = abs(o.compton/est - 1) < 0.05;  [npass, nfail] = tally(c, npass, nfail);
fprintf('%s  Compton rate vs asymptotic-sigma estimate at T = 1 GeV: ratio %.3f\n', ternary(c), o.compton/est);
% 5: annihilation, mphi = 0: sigma = alpha/24 const, R = 4 T/(32pi^4) (alpha/24) ∫ s^{3/2} K1 ds
%    = 4 T/(32 pi^4)(alpha/24) * 2 T^6 ∫ y^4 K1 dy = ... ∫_0^∞ y^4 K1(y) dy = 16
Ga_ref = 4*Th/(32*pi^4)*(alpha/24)*2*16*Th^5/(Th^3/pi^2);
c = abs(o.annihilation/Ga_ref - 1) < 1e-6;  [npass, nfail] = tally(c, npass, nfail);
fprintf('%s  annihilation, massless, T = 1 GeV: rel err %.1e\n', ternary(c), abs(o.annihilation/Ga_ref-1));
fprintf('\ntest_primakoff_phi: %d passed, %d failed\n', npass, nfail);
if nfail > 0, error('test_primakoff_phi:fail', '%d failures', nfail); end

function [p, f] = tally(c, p, f)
if c, p = p + 1; else, f = f + 1; end
end
function s = ternary(c)
if c, s = 'PASS'; else, s = 'FAIL'; end
end
