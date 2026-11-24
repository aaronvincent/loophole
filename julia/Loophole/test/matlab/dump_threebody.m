% Reference dump for the SP/PP three-body final states (threepi_weight.m,
% pipiphi_weight.m, threepi_sigma.m, Gamma(phi -> 3pi), svavg) and the
% axialVector 3pi and pi+pi-A' channels (axial3pi_weight.m, pipiA_weight.m,
% axial3pi_sigma.m, Gamma(A' -> 3pi)).  Must run with
% LOOP_3PI=1 LOOP_PIPIPHI=1 in the environment:
%   LOOP_3PI=1 LOOP_PIPIPHI=1 matlab -batch "dump_threebody"
% Long format columns: interaction, gtag, mx, mphi, gx, quantity, arg, value
assert(threepi_enabled_env() && pipiphi_enabled_env(), ...
    'dump_threebody: set LOOP_3PI=1 and LOOP_PIPIPHI=1');
here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'ref');
comb   = fullfile(here, '..', '..', '..', '..', 'matlab');
addpath(comb);
warning('off', 'all');

fid = fopen(fullfile(refdir, 'threebody_ref.csv'), 'w');
fprintf(fid, 'interaction,gtag,mx,mphi,gx,quantity,arg,value\n');
w = @(ic, gt, mx, mphi, gx, q, a, v) fprintf(fid, ...
    '%s,%s,%.16e,%.16e,%.16e,%s,%.16e,%.16e\n', ic, gt, mx, mphi, gx, q, a, v);

interactions = {'heavyScalarSP', 'heavyScalarPP'};
gpats = {'deg', [1 1 1]/3; 'gen', [0.5 0.3 0.2]};
pts   = [0.2 1.0 1; 0.15 0.6 0.5; 0.25 0.02 0.01; 0.3 1e-6 0.01];   % mx mphi gx

for ii = 1:numel(interactions)
  ic = interactions{ii};  isSP = contains(ic, 'SP');
  for gi = 1:size(gpats,1)
    gt = gpats{gi,1}; g = gpats{gi,2};
    for pi = 1:size(pts,1)
      mx = pts(pi,1); mphi = pts(pi,2); gx = pts(pi,3);
      G = phi_partial_widths(mphi, mx, gx, g, ic);
      w(ic, gt, mx, mphi, gx, 'w3pi', 0, G.pipi);
      w(ic, gt, mx, mphi, gx, 'wtot', 0, G.total);
      [~, ~, ~, Gam] = sppp_couplings(mphi, mx, gx, g, ic);
      sgrid = logspace(log10(4*mx^2*1.02), log10(4), 12);
      for s = sgrid
        w(ic, gt, mx, mphi, gx, 'W3pi',  s, threepi_weight(s, g));
        w(ic, gt, mx, mphi, gx, 'Wppf',  s, pipiphi_weight(s, mphi, g));
        w(ic, gt, mx, mphi, gx, 'sig3',  s, threepi_sigma(s, mx, mphi, gx, Gam(2), g, isSP));
        w(ic, gt, mx, mphi, gx, 'xxall', s, xxggxsec(mx, mphi, gx, g, s, ic));
        w(ic, gt, mx, mphi, gx, 'xxhad', s, xxggxsec(mx, mphi, gx, g, s, ic, 'had'));
      end
    end
  end
end
% one <sigma v> point per family, near T_RH where the 3-body channels matter
for ii = 1:numel(interactions)
  ic = interactions{ii};
  mx = 0.2; mphi = 1.0; gx = 1; g = [0.5 0.3 0.2];
  for T = [0.01 0.03]
    w(ic, 'gen', mx, mphi, gx, 'sv_all', T, svavg(mx, mphi, gx, g, T, ic));
    w(ic, 'gen', mx, mphi, gx, 'sv_sm',  T, svavg(mx, mphi, gx, g, T, ic, 'sm'));
  end
end
% axialVector: chi chi -> 3pi (J = 0 + J = 1), pi+pi-A', and Gamma(A' -> 3pi)
ic = 'axialVector';  g = 0.01*[2 -1 -1];
pts = [0.3 1.0 0.01; 0.25 0.6 0.01; 0.3 0.1 0.01; 0.2 0.3 1e-3];   % mx mA lambda_a
for pi = 1:size(pts,1)
  mx = pts(pi,1); mphi = pts(pi,2); gx = pts(pi,3);
  G = phi_partial_widths(mphi, mx, gx, g, ic);
  w(ic, 'ax', mx, mphi, gx, 'w3pi', 0, G.pipi);
  w(ic, 'ax', mx, mphi, gx, 'wtot', 0, G.total);
  sgrid = [logspace(log10(4*mx^2*1.02), log10(4), 12), mphi^2*[0.999 1 1.001]];
  for s = sgrid(sgrid > 4*mx^2)
    [WT, WL] = axial3pi_weight(s, g);
    w(ic, 'ax', mx, mphi, gx, 'WT',    s, WT);
    w(ic, 'ax', mx, mphi, gx, 'WL',    s, WL);
    [WT, WL] = pipiA_weight(s, mphi, g);
    w(ic, 'ax', mx, mphi, gx, 'WTa',   s, WT);
    w(ic, 'ax', mx, mphi, gx, 'WLa',   s, WL);
    w(ic, 'ax', mx, mphi, gx, 'sig3',  s, axial3pi_sigma(s, mx, mphi, gx, G.total, g));
    w(ic, 'ax', mx, mphi, gx, 'xxall', s, xxggxsec(mx, mphi, gx, g, s, ic));
    w(ic, 'ax', mx, mphi, gx, 'xxhad', s, xxggxsec(mx, mphi, gx, g, s, ic, 'had'));
  end
  T = mx/20;
  w(ic, 'ax', mx, mphi, gx, 'sv_all', T, svavg(mx, mphi, gx, g, T, ic));
  w(ic, 'ax', mx, mphi, gx, 'sv_sm',  T, svavg(mx, mphi, gx, g, T, ic, 'sm'));
end
fclose(fid);
fprintf('wrote %s\n', fullfile(refdir, 'threebody_ref.csv'));

function tf = threepi_enabled_env(), tf = strcmp(getenv('LOOP_3PI'), '1'); end
function tf = pipiphi_enabled_env(), tf = strcmp(getenv('LOOP_PIPIPHI'), '1'); end
