% Regenerate the M5 reference dump for the multi-pole models (SP/PP family +
% axialVector): the cross sections sigma_ggxx (ggxxxsec), sigma_xxgg
% (make_xxsigma 'gg'), sigma_xxall (xxggxsec) over an s-grid, plus the mediator
% widths (phi_partial_widths), the meson poles (sppp/axial_couplings), and
% sigma_nucleon (DMnucxsec) -- for both degenerate g=[1 1 1]/3 (the {eta',phi'}
% singlet switch) and generic g (the {pi0,phi',eta} set).
% Long format columns: interaction, gtag, mx, mphi, gx, quantity, arg, value
here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'ref');
comb   = fullfile(here, '..', '..', '..', '..', 'matlab');
addpath(comb);
warning('off', 'all');

fid = fopen(fullfile(refdir, 'm5_ref.csv'), 'w');
fprintf(fid, 'interaction,gtag,mx,mphi,gx,quantity,arg,value\n');

w = @(ic, gt, mx, mphi, gx, q, a, v) fprintf(fid, ...
    '%s,%s,%.16e,%.16e,%.16e,%s,%.16e,%.16e\n', ic, gt, mx, mphi, gx, q, a, v);

interactions = {'heavyScalarSP','lightScalarSP','scalarPP','heavyScalarPP', ...
                'lightScalarPP','axialVector'};
gpats = {'deg', [1 1 1]/3; 'gen', [0.5 0.3 0.2]};
pts   = [0.05 0.8 0.1; 0.3 0.5 0.01; 0.7 0.9 0.1];   % mx mphi gx

for ii = 1:numel(interactions)
  ic = interactions{ii};
  for gi = 1:size(gpats,1)
    gt = gpats{gi,1}; g = gpats{gi,2};
    for pi = 1:size(pts,1)
      mx = pts(pi,1); mphi = pts(pi,2); gx = pts(pi,3);

      % cross sections over an s-grid (just above 4mx^2 up to 60)
      s_lo  = 4*mx^2*1.05;
      sgrid = logspace(log10(s_lo), log10(60), 10);
      sig_gg = make_xxsigma(mx, mphi, gx, g, ic, 'gg');
      for s = sgrid
        w(ic, gt, mx, mphi, gx, 'ggxx',  s, ggxxxsec(mx,mphi,gx,g,s,ic));
        w(ic, gt, mx, mphi, gx, 'xxgg',  s, sig_gg(s));
        w(ic, gt, mx, mphi, gx, 'xxall', s, xxggxsec(mx,mphi,gx,g,s,ic));
      end

      % mediator widths
      G = phi_partial_widths(mphi, mx, gx, g, ic);
      w(ic, gt, mx, mphi, gx, 'wgg',  0, G.gg);
      w(ic, gt, mx, mphi, gx, 'wxx',  0, G.xx);
      w(ic, gt, mx, mphi, gx, 'wtot', 0, G.total);

      % meson poles (s_p = mj^2, w_p = mj*Gam)
      if strcmp(ic, 'axialVector')
        [mj, ~, ~, Gam] = axial_couplings(mphi, mx, gx, g);
      else
        [mj, ~, ~, Gam] = sppp_couplings(mphi, mx, gx, g, ic);
      end
      for j = 1:numel(mj)
        w(ic, gt, mx, mphi, gx, 'pole_s', j, mj(j)^2);
        w(ic, gt, mx, mphi, gx, 'pole_w', j, mj(j)*Gam(j));
      end

      % DM-nucleon cross section [GeV^-2] -- the FULL angular-integrated form, to
      % match the Julia sigma_nucleon (full-only). full_sigma_name maps a heavy*
      % name to its 'scalar*' full sibling (scalarPP/axialVector pass through).
      w(ic, gt, mx, mphi, gx, 'signuc', 0, DMnucxsec(mx,mphi,gx,g,full_sigma_name(ic)));
    end
  end
end

% --- end-to-end model_observables for a few multi-pole points: exercises the
%     generic multi-pole thermal_pole_integral path (svavg/ggxxrate) + the relic
%     solvers, which the per-s sigma checks above don't cover. T_RH=0.01.
T_RH = 0.01;
obs_pts = { 'heavyScalarSP', 'gen', 0.05, 0.8,  0.1;
            'scalarPP',      'deg', 0.05, 0.8,  0.1;
            'heavyScalarPP', 'gen', 0.05, 0.5,  0.1;
            'axialVector',   'gen', 0.05, 0.8,  0.1;
            'axialVector',   'deg', 0.02, 0.05, 0.01 };   % freeze_coupled (mphi<=0.1)
for k = 1:size(obs_pts,1)
  ic = obs_pts{k,1}; gt = obs_pts{k,2};
  if strcmp(gt,'deg'), g = [1 1 1]/3; else, g = [0.5 0.3 0.2]; end
  mx = obs_pts{k,3}; mphi = obs_pts{k,4}; gx = obs_pts{k,5};
  o = model_observables(mx, mphi, gx, g, ic, T_RH);
  w(ic, gt, mx, mphi, gx, 'Oh2',   0, o.Oh2);
  w(ic, gt, mx, mphi, gx, 'dNbbn', 0, o.dNeff_bbn);
  w(ic, gt, mx, mphi, gx, 'dNcmb', 0, o.dNeff_cmb);
end

fclose(fid);
disp('dump_m5 done');
