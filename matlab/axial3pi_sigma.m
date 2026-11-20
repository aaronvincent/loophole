function sig = axial3pi_sigma(s, mx, mA, gx, GA, g)
% axial3pi_sigma  chi chibar -> 3pi + pi+pi-A' for the axial-vector mediator A'.
% s vector [GeV^2], all s > 4 mx^2.  mA = m_A', gx = lambda_a, GA = Gamma(A')
% total.  J = 1 (transverse, resonant at s = m_A'^2) plus J = 0 (longitudinal,
% no A' pole); they do not interfere.  Weights and switches: axial3pi_weight.m
% (LOOP_3PI), pipiA_weight.m (LOOP_PIPIPHI).
[WT, WL]   = axial3pi_weight(s, g);
[WTa, WLa] = pipiA_weight(s, mA, g);
WT = WT + WTa;  WL = WL + WLa;
beta = sqrt(1 - 4*mx^2./s);
sig  = gx^2/(32*pi^3) * (beta.*WT ./ (6*((s - mA^2).^2 + mA^2*GA^2)) ...
                         + mx^2*WL ./ (beta.*s.^2*mA^4));
end
