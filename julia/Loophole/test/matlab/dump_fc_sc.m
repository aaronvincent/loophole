% Regenerate test/ref/fc_sc_ref.csv for validate_fc_sc.jl: freeze_coupled where
% phi's 2 -> 2 conversions (Primakoff, pi gamma -> pi phi) fill a phi that
% gamma gamma -> phi cannot (mphi = 10 mx, small gx, small g_q), plus points
% where they are inert.  Each row is evaluated with the conversions on and off
% (LOOP_PHI_SC=0).  A second block covers chi chibar <-> phi phi in every family.
% Columns: interaction, mx, mphi, gx, g1, g2, g3, T_RH,
% Oh2_on, Oh2_off
here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'ref');
addpath(fullfile(here, '..', '..', '..', '..', 'matlab'));
warning('off', 'all');
fid = fopen(fullfile(refdir, 'fc_sc_ref.csv'), 'w');
fprintf(fid, 'interaction,mx,mphi,gx,g1,g2,g3,T_RH,Oh2_on,Oh2_off\n');
T_RH = 0.01;
for ic = {'heavyScalar', 'heavyScalarPP'}
  for mphi = [1e-4, 1e-3, 1e-2]
    for gx = [1e-9, 1e-5]
      for gq = [1e-6, 1e-5, 1e-3]
        g = gq*[2 -1 -1];  mx = mphi/10;
        setenv('LOOP_PHI_SC', '');   on  = freeze_coupled(mx, mphi, gx, g, ic{1}, T_RH);
        setenv('LOOP_PHI_SC', '0');  off = freeze_coupled(mx, mphi, gx, g, ic{1}, T_RH);
        fprintf(fid, '%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
                ic{1}, mx, mphi, gx, g, T_RH, on, off);
      end
    end
  end
end
% chi chibar <-> phi phi block (every family, incl. axialVector's g_phi = 3):
% mx < mphi < 2 mx (phi cannot decay to chi chibar), mx > mphi, mphi > 2 mx.
P = [0.006 0.01 1    1e-3;  0.006 0.01 0.01 1e-5;  0.02 0.01 0.1 1e-3;
     0.001 0.01 1    1e-3;  0.03  0.1  0.3  1e-4];
for ic = {'heavyScalar', 'heavyScalarPS', 'heavyScalarSP', 'heavyScalarPP', 'axialVector'}
  for k = 1:size(P, 1)
    mx = P(k,1);  mphi = P(k,2);  gx = P(k,3);  g = P(k,4)*[2 -1 -1];
    setenv('LOOP_PHI_SC', '');   on  = freeze_coupled(mx, mphi, gx, g, ic{1}, T_RH);
    setenv('LOOP_PHI_SC', '0');  off = freeze_coupled(mx, mphi, gx, g, ic{1}, T_RH);
    fprintf(fid, '%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
            ic{1}, mx, mphi, gx, g, T_RH, on, off);
  end
end
setenv('LOOP_PHI_SC', '');
fclose(fid);
disp('dump_fc_sc done');
