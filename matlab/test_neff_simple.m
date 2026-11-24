function test_neff_simple()
% test_neff_simple  F(y) normalisation and the
% massless limits of the Boehm-Dolan-McCabe N_eff formula in neff_simple.m.
%
%   test_neff_simple        % errors on failure

g = [2 -1 -1];  gx = 1;  I = 'lightScalar';
nfail = 0;

% F(y) is read back through info with T_nu_dec = 1, so y = m.  The coupling
% flags are irrelevant here; F is computed regardless.
F = @(y, which) Fvia(y, which, gx, g, I);

% --- normalisation and limits ---
nfail = nfail + check('F_fermion(0) = 1',     F(1e-8,'chi'), 1,   1e-6);
nfail = nfail + check('F_boson(0) = 8/7',     F(1e-8,'phi'), 8/7, 1e-6);
nfail = nfail + check('F_fermion(50) ~ 0',    F(50,'chi'),   0,   1e-15, 'abs');
ys = [0.01 0.1 0.5 1 2 5 10];
Ff = arrayfun(@(y) F(y,'chi'), ys);  Fb = arrayfun(@(y) F(y,'phi'), ys);
nfail = nfail + assert_true('F decreasing (fermion)', all(diff(Ff) < 0));
nfail = nfail + assert_true('F decreasing (boson)',   all(diff(Fb) < 0));

% --- massless limits, at points where the rates couple both species at
% T_nu_dec = 2 MeV (mphi = 1e-5: y_phi = 5e-3; mx = 1e-7: y_chi = 5e-5) ---
[d, In] = neff_simple(1e-5, 1e-7, gx, g, I);
nfail = nfail + assert_true('limit point: phi and chi coupled', In.phi_coupled && In.chi_coupled);
nfail = nfail + check('phi+chi massless: X = 9/11', In.X, 9/11, 1e-4);
nfail = nfail + check('phi+chi massless: dNeff', d, 3.046*((20/11)^(-4/3) - 1), 1e-4);
[d, In] = neff_simple(1e-5, 1e-7, gx, g, I, struct('species','phi'));
nfail = nfail + check('phi only massless: X = 2/11', In.X, 2/11, 1e-4);
nfail = nfail + check('phi only massless: dNeff', d, 3.046*((13/11)^(-4/3) - 1), 1e-4);

% --- phi thermalised through chi, left as free radiation when chi annihilates:
% bath {gamma, e, chi, phi} at T_nu_dec (S = 10), {gamma, phi} after, so
% T_phi = T_gamma = (10/3)^(1/3) T_nu.  P_free divides out phi decays.
% phi must have NO photon link of its own: the gamma gamma <-> phi phi meson
% loop goes as g_q^4 and couples phi directly from g_q ~ 1e-2, so take
% g_q = 2e-3 (w_phi ~ 1e-3) and a large gx for full chi <-> phi contact. ---
[~, In] = neff_simple(1e-11, 1e-4, 30, 2e-3*[1 1 1]/3, 'heavyScalar');
nfail = nfail + check('phi via chi: T_phi/T_nu = (10/3)^(1/3)', ...
    In.Tphi_over_Tnu / In.P_free^(1/4), (10/3)^(1/3), 2e-2);

% --- uncoupled -> 0 ---
[d, In] = neff_simple(1e-9, 1e-4, gx, 1e-9*g, I);
nfail = nfail + assert_true('tiny gq: nothing coupled, dNeff = 0', ...
    ~In.phi_coupled && ~In.chi_coupled && abs(d) < 1e-5);

if nfail > 0
    error('test_neff_simple:fail', '%d check(s) failed', nfail);
end
fprintf('test_neff_simple: all checks PASS\n');
end

function v = Fvia(y, which, gx, g, I)
[~, In] = neff_simple(y, y, gx, g, I, struct('T_nu_dec', 1));
v = In.(['F_' which]);
end

function f = check(name, got, want, tol, mode)
if nargin < 5, mode = 'rel'; end
if strcmp(mode, 'abs'), err = abs(got - want); else, err = abs(got - want)/abs(want); end
f = ~(err <= tol);
fprintf('%-40s got %.8g  want %.8g  err %.2e  %s\n', name, got, want, err, ifs(f));
end

function f = assert_true(name, c)
f = ~c;
fprintf('%-40s %s\n', name, ifs(f));
end

function s = ifs(f)
if f, s = 'FAIL'; else, s = 'ok'; end
end
