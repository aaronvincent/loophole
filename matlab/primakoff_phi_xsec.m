function [sig, sig_ann] = primakoff_phi_xsec(s, mphi, D)
% primakoff_phi_xsec  sigma(e gamma -> e phi) and sigma(e+ e- -> gamma phi) at
% g_phigamgam = 1 [GeV^-2], massless electrons, vectorised in s [GeV^2].
%
%   [sig, sig_ann] = primakoff_phi_xsec(s, mphi, D)      D = mD^2 [GeV^2]
%
% sigma_C = alpha I/(16 s^2), I as in primakoff_phi_rate.m (closed form for
% X/D >= 1e-2, X = s - mphi^2; a three-term series in X/D below, where the
% closed form cancels).  sigma_A = alpha (1 - mphi^2/s)^3 / 24.
alpha = 1/137.036;
X = max(s - mphi^2, 0);
I = zeros(size(s));
r = X ./ D;
big = r >= 1e-2;
if any(big(:))
    Xb = X(big);  sb = s(big);
    B = Xb + D;  C = sb.^2 + B.^2;
    % F(D+X) - F(D), with the D C / y term combined as -C X/(D+X)
    I(big) = (Xb.^2 + 2*D*Xb)/2 - (2*B + D).*Xb + (C + 2*B*D).*log1p(Xb/D) ...
             - C.*Xb./(D + Xb);
end
sm = ~big & X > 0;
if any(sm(:))
    % X << D: expand 1/(x+D)^2 = (1 - 2x/D + 3x^2/D^2)/D^2 under the integral
    Xs = X(sm);  ss = s(sm);
    m0 = ss.^2.*Xs.^2/2 + Xs.^4/12;              % ∫ x (s^2 + (X-x)^2)
    m1 = ss.^2.*Xs.^3/3 + Xs.^5/30;              % ∫ x^2 (...)
    m2 = ss.^2.*Xs.^4/4 + Xs.^6/60;              % ∫ x^3 (...)
    I(sm) = (m0 - 2*m1/D + 3*m2/D^2)/D^2;
end
sig = alpha * I ./ (16 * s.^2);
sig(s <= 0) = 0;
sig_ann = alpha * max(0, 1 - mphi^2./s).^3 / 24;
sig_ann(s <= mphi^2) = 0;
end
