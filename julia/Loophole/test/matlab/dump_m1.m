% Regenerate the M1 reference dumps for the Julia validation.
% Run from anywhere:  matlab -batch "run('.../test/matlab/dump_m1.m')"
% Writes test/ref/m1_ref.csv and test/ref/m1_widths.csv.
here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'ref');
comb   = fullfile(here, '..', '..', '..', '..', 'matlab');
addpath(comb);

fid = fopen(fullfile(refdir, 'm1_ref.csv'), 'w');
fprintf(fid, 'tag,s,ggxx,xxgg,xxall\n');
wfid = fopen(fullfile(refdir, 'm1_widths.csv'), 'w');
fprintf(wfid, 'tag,gg,xx,pipi,KK,total\n');

pts = { ...
  struct('tag','hs_a',   'int','heavyScalar',  'mx',0.1,'mphi',10,'gx',1,'g',[1 1 1]/3), ...
  struct('tag','hsps_a', 'int','heavyScalarPS','mx',0.1,'mphi',10,'gx',1,'g',[1 1 1]/3), ...
  struct('tag','hs_pp',  'int','heavyScalar',  'mx',20, 'mphi',10,'gx',1,'g',[1 1 1]/3), ...
  struct('tag','hsps_pp','int','heavyScalarPS','mx',20, 'mphi',10,'gx',1,'g',[1 1 1]/3) ...
};

for ii = 1:numel(pts)
    p = pts{ii};
    G = phi_partial_widths(p.mphi, p.mx, p.gx, p.g, p.int);
    fprintf(wfid, '%s,%.16e,%.16e,%.16e,%.16e,%.16e\n', p.tag, G.gg, G.xx, G.pipi, G.KK, G.total);

    s_lo = 4*p.mx^2;
    if p.mx >= p.mphi
        svec = sort([logspace(log10(s_lo*1.0001), log10(1e6), 50), 4*p.mphi^2*1.01]);
    else
        svec = sort([logspace(log10(s_lo*1.0001), log10(1e4), 50), 99, 100, 101]);
    end
    [siggg, ~] = make_xxsigma(p.mx, p.mphi, p.gx, p.g, p.int, 'gg');
    for s = svec
        a = ggxxxsec(p.mx, p.mphi, p.gx, p.g, s, p.int);
        b = siggg(s);
        c = xxggxsec(p.mx, p.mphi, p.gx, p.g, s, p.int);
        fprintf(fid, '%s,%.16e,%.16e,%.16e,%.16e\n', p.tag, s, a, b, c);
    end
end
fclose(fid); fclose(wfid);
disp('dump_m1 done');
