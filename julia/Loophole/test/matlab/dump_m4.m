% Regenerate the M4 reference dump: model_observables (relic + dNeff bbn/cmb +
% sigma_chiN) over a (interaction, mx, mphi, gx) grid spanning the freeze_loop
% (heavy mediator) and freeze_coupled (sub-MeV mediator) regimes, light/heavy
% chi (the dNeff_chi gate), and a sub-MeV phi (the CMB photon-heating branch).
% Columns: interaction, mx, mphi, gx, g1, g2, g3, T_RH, Oh2, dNeff_bbn,
%          dNeff_cmb, sigma_chiN
here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'ref');
comb   = fullfile(here, '..', '..', '..', '..', 'matlab');
addpath(comb);
warning('off', 'all');

fid = fopen(fullfile(refdir, 'm4_ref.csv'), 'w');
fprintf(fid, 'interaction,mx,mphi,gx,g1,g2,g3,T_RH,Oh2,dNeff_bbn,dNeff_cmb,sigma_chiN\n');

T_RH = 0.01;
g    = [1 1 1]/3;

interactions = {'heavyScalar','lightScalar','heavyScalarPS','lightScalarPS'};
mx_vec   = [0.005, 0.05];
mphi_vec = [1e-4, 0.05, 10];
gx_vec   = [1e-3, 1];

rows = {};
for ii = 1:numel(interactions)
  ic = interactions{ii};
  for mx = mx_vec
    for mphi = mphi_vec
      for gx = gx_vec
        rows(end+1,:) = {ic, mx, mphi, gx, g}; %#ok<AGROW>
      end
    end
  end
end

% explicit anchor points (the two verified MATLAB points)
rows(end+1,:) = {'heavyScalarPS', 0.1,  10,   1,      [1 1 1]/3};
rows(end+1,:) = {'heavyScalar',   0.01, 1e-4, 0.0035, [1 1 1]};

for r = 1:size(rows,1)
  ic = rows{r,1}; mx = rows{r,2}; mphi = rows{r,3}; gx = rows{r,4}; gg = rows{r,5};
  o = model_observables(mx, mphi, gx, gg, ic, T_RH);
  fprintf(fid, '%s,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n', ...
    ic, mx, mphi, gx, gg(1), gg(2), gg(3), T_RH, ...
    o.Oh2, o.dNeff_bbn, o.dNeff_cmb, o.sigma_chiN);
end

fclose(fid);
disp('dump_m4 done');
