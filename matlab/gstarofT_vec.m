function gstar = gstarofT(T,sflag)
%Aaron Vincent 
%input is GeV
%sflag: use gstars; no sflag, use gstar not s
%returns gstar as computed by the interpolation in http://arxiv.org/pdf/0910.1066.pdf
if isrow(T)
    rowflag = 1;
    T = T';
else
    rowflag = 0;
end
L = length(T);
if ~sflag
    a0 = 1.21;
  a1 =repmat([0.572   0.330   0.579   0.138   0.108],L,1);
  a2=repmat([ -8.77   -2.95   -1.80   -0.162   3.76],L,1);
  a3=repmat([ 0.682   1.01   0.165   0.934   0.869],L,1);
  T = log([T T T T T]);
  II = ones(L,5);
  
  gstar = a0.*ones(L,1) + sum(a1.*(II+tanh((T-a2)./a3)),2);

else
  aS0 =1.36;
  aS1 = repmat([ 0.498   0.327   0.579   0.140   0.109],L,1);
  aS2 = repmat([ -8.74   -2.89   -1.79   -0.102   3.82],L,1);
  aS3 = repmat([ 0.693   1.01   0.155   0.963   0.907],L,1);
  
  T = log([T T T T T]);
    II = ones(L,5);
  
    gstar = aS0.*ones(L,1) + sum(aS1.*(II+tanh((T-aS2)./aS3)),2);
end
  gstar = exp(gstar);
  if rowflag
      gstar = gstar';
  end
