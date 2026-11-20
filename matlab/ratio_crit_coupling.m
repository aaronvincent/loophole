function gc = ratio_crit_coupling(rfun, g0)
%RATIO_CRIT_COUPLING  Coupling at which an exclusion ratio first reaches 1.
%
%   GC = RATIO_CRIT_COUPLING(RFUN, G0) returns the smallest coupling GC with
%   RFUN(GC) = 1 -- the lower edge of the excluded region.  G0 is a starting
%   guess, normally the gx^2-rescaled estimate 1/sqrt(RFUN(1)).
%
%   Why this is not just G0: the rescale assumes ratio ~ gx^2.  That fails for
%   K+ -> pi+ chi chi whenever the phi pole sits inside the Dalitz region
%   (2 mx < mphi < mK+ - mpi+).  There the E_pi integral is saturated by the
%   resonance, so F ~ 1/Gamma_phi, and Gamma(phi -> chi chibar) ~ gx^2 cancels
%   the explicit gx^2: the rate plateaus and G0 is orders of magnitude too weak
%   (~200x in gx at mphi = 0.1 GeV).  Off resonance the ratio really is ~ gx^2
%   and the root-find returns G0 to solver tolerance.
%
%   The search runs in log10(coupling).  Returns NaN if no root can be
%   bracketed (RFUN non-finite, or never crossing 1 in ~60 decades).

if ~(isfinite(g0) && g0 > 0)
    gc = NaN;  return;
end

f = @(L) local_logratio(rfun, L);      % zero when ratio = 1

% Upper bracket: need ratio >= 1.  A resonance can only ENHANCE the rate over
% the gx^2 extrapolation, so g0 is already excluded and L0 usually suffices.
% Quarter-decade steps: kphi is excluded only in a band when the phi decays
% visibly inside the detector at larger coupling, and it can be under a decade
% wide.
L0 = log10(g0);
Lhi = L0;  fhi = f(Lhi);  k = 0;
while ~(isfinite(fhi) && fhi >= 0) && k < 160
    Lhi = Lhi + 0.25;  fhi = f(Lhi);  k = k + 1;
end
% Not found above g0: the ratio may fall again at large coupling (kphi: the phi
% decays visibly inside the detector), making g0 = 1/sqrt(RFUN(1)) far too
% high.  Look below g0 instead.
if ~(isfinite(fhi) && fhi >= 0)
    Lhi = L0 - 0.25;  fhi = f(Lhi);  k = 0;
    while ~(isfinite(fhi) && fhi >= 0) && k < 240
        Lhi = Lhi - 0.25;  fhi = f(Lhi);  k = k + 1;
    end
end
if ~(isfinite(fhi) && fhi >= 0), gc = NaN;  return; end
if fhi == 0, gc = 10^Lhi;  return; end

% Lower bracket: need ratio < 1.  Well below the resonance plateau the width is
% gamma-gamma dominated, F goes flat, and ratio ~ gx^2 again -> a root exists.
Llo = Lhi - 1;  flo = f(Llo);  k = 0;
while ~(isfinite(flo) && flo < 0) && k < 60
    Llo = Llo - 1;  flo = f(Llo);  k = k + 1;
end
if ~(isfinite(flo) && flo < 0), gc = NaN;  return; end

gc = 10^fzero(f, [Llo Lhi]);
end

% -------------------------------------------------------------------------
function y = local_logratio(rfun, L)
% log(ratio), with non-finite / non-positive ratios mapped to NaN so the
% bracket expanders above treat them as "unusable" rather than as a sign.
r = rfun(10^L);
if isfinite(r) && r > 0
    y = log(r);
else
    y = NaN;
end
end
