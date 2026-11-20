% Regenerate the M3 relic-solver reference dumps.
% Columns: solver, interaction, mx, mphi, gx, T_RH, Oh2
here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'ref');
comb   = fullfile(here, '..', '..', '..', '..', 'matlab');
addpath(comb);
alltheconstants;                 % Mpl, s0, rhoc, T0
warning('off', 'all');

fid = fopen(fullfile(refdir, 'm3_ref.csv'), 'w');
fprintf(fid, 'solver,interaction,mx,mphi,gx,T_RH,Oh2\n');

g = [1 1 1]/3;
T_RH = 0.01;

% --- freeze_loop / freeze_out / Oh2_FI : heavy mediator (mphi=10) ---
mphi = 10;
for interaction = {'heavyScalar','heavyScalarPS'}
  ic = interaction{1};
  for mx = [0.01, 0.1]
    for gx = [1e-4, 1e-2, 1, 5]
      v = freeze_loop(mx, mphi, gx, g, ic, T_RH);
      fprintf(fid, 'freeze_loop,%s,%.16e,%.16e,%.16e,%.16e,%.16e\n', ic, mx, mphi, gx, T_RH, v);
      v = freeze_out(mx, mphi, gx, g, ic);
      fprintf(fid, 'freeze_out,%s,%.16e,%.16e,%.16e,%.16e,%.16e\n', ic, mx, mphi, gx, T_RH, v);
      v = Oh2_FI(mx, mphi, gx, g, ic, T_RH, Mpl, s0, rhoc, T0);
      fprintf(fid, 'Oh2_FI,%s,%.16e,%.16e,%.16e,%.16e,%.16e\n', ic, mx, mphi, gx, T_RH, v);
    end
  end
end

% --- freeze_coupled : light mediator (mphi <= 10 T_RH) ---
for mphi = [1e-4, 0.05]
  for mx = [0.005, 0.05]
    for gx = [1e-3, 0.0035, 0.1]
      v = freeze_coupled(mx, mphi, gx, g, 'heavyScalar', T_RH);
      fprintf(fid, 'freeze_coupled,heavyScalar,%.16e,%.16e,%.16e,%.16e,%.16e\n', mx, mphi, gx, T_RH, v);
    end
  end
end

fclose(fid);
disp('dump_m3 done');
