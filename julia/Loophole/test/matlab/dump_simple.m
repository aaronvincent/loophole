% Regenerate the neff_simple reference dump (test/ref/simple_ref.csv) for
% validate_simple.jl.  Grid spans: phi via its photon link (light phi, large
% g_q), chi via gamma gamma <-> chi chi, phi made by chi (ultralight phi, gx = 1),
% phi decays before recombination, phi via gamma gamma <-> phi phi (light phi), heavy phi (chi only), chi decoupled,
% the pi gamma -> pi phi turn-on (decoupled-before-T_nu_dec phi as free radiation),
% and Primakoff e gamma -> e phi keeping light phi in the bath to T_late < T_nu_dec.
% Columns: interaction, mx, mphi, gx, g1, g2, g3, T_RH, dNeff, dNeff_heat,
%          dNeff_phi_rad, w_phi, w_chi
here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'ref');
comb   = fullfile(here, '..', '..', '..', '..', 'matlab');
addpath(comb);
warning('off', 'all');

fid = fopen(fullfile(refdir, 'simple_ref.csv'), 'w');
fprintf(fid, 'interaction,mx,mphi,gx,g1,g2,g3,T_RH,dNeff,dNeff_heat,dNeff_phi_rad,w_phi,w_chi\n');
T_RH = 0.01;
rows = {};
for ic = {'heavyScalar', 'heavyScalarPS'}
  for mphi = [1e-10, 1e-4, 0.1, 1]
    for mx = [1e-4, 3e-3]
      for gq = [1e-3, 1e-1, 1]
        for gx = [1e-5, 1]
          rows(end+1,:) = {ic{1}, mx, mphi, gx, gq*[1 1 1]/3}; %#ok<AGROW>
        end
      end
    end
  end
end
rows(end+1,:) = {'axialVector', 1e-3, 0.1, 1e-3, 1e-3*[1 1 1]/3};
rows(end+1,:) = {'lightScalar', 3.5e-2, 1e-10, 1, 1e-3*[1 1 1]/3};
for mphi = [1e-9, 1e-6]                    % gamma gamma <-> phi phi via the SP/PP vertex
  rows(end+1,:) = {'heavyScalarPP', 1e-3, mphi, 1e-4, [1 1 1]/3}; %#ok<AGROW>
end
for mphi = [1e-9, 1e-3]                    % pi gamma -> pi phi turn-on (gq_th ~ 2e-6)
  for gq = [1e-6, 3e-6, 1e-5]
    rows(end+1,:) = {'heavyScalar', 1e-2, mphi, 1e-5, gq*[2 -1 -1]}; %#ok<AGROW>
  end
end
for mphi = [1e-9, 1e-4]                    % Primakoff + late release (gq_th ~ 1e-5)
  for gq = [3e-5, 1e-4, 1e-3]
    rows(end+1,:) = {'heavyScalar', 1e-2, mphi, 1e-5, gq*[2 -1 -1]}; %#ok<AGROW>
  end
end
rows(end+1,:) = {'heavyScalarSP', 1e-3, 1e-9, 1e-4, 1e-2*[1 1 1]/3};
for r = 1:size(rows,1)
    [ic, mx, mphi, gx, g] = rows{r,:};
    try
        [d, I] = neff_simple(mphi, mx, gx, g, ic, struct('T_RH', T_RH));
        v = [d, I.dNeff_heat, I.dNeff_phi_rad, I.w_phi, I.w_chi];
    catch
        v = nan(1, 5);
    end
    fprintf(fid, '%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
            ic, mx, mphi, gx, g(1), g(2), g(3), T_RH, v);
end
fclose(fid);
fprintf('wrote %d rows\n', size(rows,1));
