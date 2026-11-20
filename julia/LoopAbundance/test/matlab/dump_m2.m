% Regenerate the M2 thermal-average reference dumps.
% Columns: tag, kind, T, val   (kind in svgrid_all/svgrid_gg/sv_all/sv_gg/ggrate)
here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'ref');
comb   = fullfile(here, '..', '..', '..', '..', 'matlab');
addpath(comb);

fid = fopen(fullfile(refdir, 'm2_ref.csv'), 'w');
fprintf(fid, 'tag,kind,T,val\n');

pts = { ...
  struct('tag','hs_res',  'int','heavyScalar',  'mx',0.1,'mphi',10,'gx',1,'g',[1 1 1]/3), ... % phi pole in range
  struct('tag','hsps_res','int','heavyScalarPS','mx',0.1,'mphi',10,'gx',1,'g',[1 1 1]/3), ...
  struct('tag','hs_nopole','int','heavyScalar', 'mx',8,  'mphi',10,'gx',1,'g',[1 1 1]/3), ... % pole below threshold
  struct('tag','hs_pp',   'int','heavyScalar',  'mx',20, 'mphi',10,'gx',1,'g',[1 1 1]/3) ...  % phiphi open
};

Tgrid = logspace(-3, 0, 12);

for ii = 1:numel(pts)
    p = pts{ii};
    svg_all = svavg_grid(p.mx, p.mphi, p.gx, p.g, Tgrid, p.int);
    svg_gg  = svavg_grid(p.mx, p.mphi, p.gx, p.g, Tgrid, p.int, 'gg');
    grate   = ggxxrate(Tgrid, p.mx, p.mphi, p.gx, p.g, p.int);
    for kk = 1:numel(Tgrid)
        T = Tgrid(kk);
        fprintf(fid, '%s,svgrid_all,%.16e,%.16e\n', p.tag, T, svg_all(kk));
        fprintf(fid, '%s,svgrid_gg,%.16e,%.16e\n',  p.tag, T, svg_gg(kk));
        fprintf(fid, '%s,sv_all,%.16e,%.16e\n', p.tag, T, svavg(p.mx,p.mphi,p.gx,p.g,T,p.int));
        fprintf(fid, '%s,sv_gg,%.16e,%.16e\n',  p.tag, T, svavg(p.mx,p.mphi,p.gx,p.g,T,p.int,'gg'));
        fprintf(fid, '%s,ggrate,%.16e,%.16e\n', p.tag, T, grate(kk));
    end
end
fclose(fid);
disp('dump_m2 done');
