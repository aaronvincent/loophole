% Regenerate the M7 reference dump for the Julia validation.
% Run from anywhere:  matlab -batch "run('.../test/matlab/dump_m7.m')"
% Writes test/ref/m7_ref.csv  (key,value -- one scalar per row).
%
% M7 covers the physics merged from `main` into the port: the pseudoscalar
% nucleon coupling h~_P for SP/PP (DMnucxsec), the axialVector EFT-validity
% ceilings/floors and the isolated chi chibar -> A'A' piece.
%
% NOTE on sigma_nucleon: routed through full_sigma_name, as dump_m5.m does.
% The Julia port's sigma_nucleon is ALWAYS the full angular-integrated form
% (observables.jl:40) -- HeavyScalarPP/LightScalarPP/ScalarPP share one method,
% as do HeavyScalarSP/LightScalarSP -- because the prefix is inert everywhere
% except a direct DMnucxsec call (docs/dmnucxsec.md).  The five names below therefore
% collapse to TWO distinct values on both sides, and that aliasing contract is
% what this block tests.  Until 2026-08-05 M7 called the raw heavy*/light*
% approximations instead and asserted they differ; the port never implemented
% those, so lightScalarPP/lightScalarSP failed permanently -- a test/port
% contract mismatch, not a physics defect.
here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'ref');
comb   = fullfile(here, '..', '..', '..', '..', 'matlab');
addpath(comb);

fid = fopen(fullfile(refdir, 'm7_ref.csv'), 'w');
fprintf(fid, 'key,value\n');
w = @(k, v) fprintf(fid, '%s,%.17g\n', k, v);

% --- h~_P nucleon sigma (SP/PP), g=[2 -1 -1]/3, mx=0.05, mphi=0.1, gx=1 -----
g  = [2 -1 -1]/3;
mx = 0.05;  mphi = 0.1;  gx = 1;
for ic = {'heavyScalarPP','lightScalarPP','heavyScalarSP','lightScalarSP','scalarPP'}
    w(['signuc_' ic{1}], DMnucxsec(mx, mphi, gx, g, full_sigma_name(ic{1})));
end
% degenerate g -> eta'-singlet projection zeroes h~_P -> sigma_N = 0
w('signuc_deg_heavyScalarPP', ...
  DMnucxsec(mx, mphi, gx, [1 1 1]/3, full_sigma_name('heavyScalarPP')));

% --- axialVector EFT ceilings / floors (mphi = mx = 0.05, g = [1 1 1]) ------
w('axial_gx_ceiling', axial_gx_ceiling(0.05, 0.05, [1 1 1]));
[gq_mix, mxmix, jbind] = axial_mixing_floor(0.05, 0.05, [1 1 1], 0.3);
w('axial_mixing_floor_gq_mix', gq_mix);
w('axial_mixing_floor_mxmix',   mxmix);
w('axial_mixing_floor_jbind',   jbind);
w('axial_gq_ceiling', axial_gq_ceiling(0.05, 0.05, 1e-3, [1 1 1]));

% --- chi chibar -> A'A' thermal average -------------------------------------
% Julia: AxialVector(mx=0.05, mphi=0.02, gx=1e-3, g=(1,1,1)) at T = mx/20.
w('svavg_phiphi', svavg_phiphi(0.05, 0.02, 1e-3, [1 1 1], 0.05/20, 'axialVector'));

fclose(fid);
disp('dump_m7 done');
