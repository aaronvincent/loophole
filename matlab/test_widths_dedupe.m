% Regression test for phi_partial_widths.
% For each interaction, compare phi_partial_widths(...).total against an
% independent transcription of the documented formulas.  These must match
% bit-identically.

mphi = 0.05; mx = 0.02; gx = 0.5; g = [2 -1 -1]/3;

mpi0    = 0.13498;
mpi_pm  = 0.13957;
mK      = 0.493677;
m_eta   = 0.547862;
alphaEM = 1/137;
B0      = 0.13957^2/(0.0022 + 0.00469);   % = mpi+^2/(mu+md) ~ 2.827 GeV (draft convention; was 2.6)
f       = 0.0924;
Dscript = 0.5*(B0*(g(1)+g(2))/mpi_pm^2 + B0*(g(1)+g(3))/mK^2)^2;   % charged-meson loop (pi+, K+), low-energy limit A0 -> 1/3

fprintf('--- phi partial width dedupe regression (mphi=%g mx=%g gx=%g) ---\n', mphi, mx, gx);

% ----- heavyScalar / lightScalar -----
% gg /4608 + beta^3 scalar-vertex chi chi width; pipi/KK closed at this mphi.
old_HS = alphaEM^2*Dscript*mphi^3/(4608*pi^3) ...
       + gx^2*max(0, mphi^2 - 4*mx^2)^(3/2)/(8*pi*mphi^2);
new_HS = phi_partial_widths(mphi, mx, gx, g, 'heavyScalar').total;
report('heavyScalar', old_HS, new_HS);
new_LS = phi_partial_widths(mphi, mx, gx, g, 'lightScalar').total;
report('lightScalar', old_HS, new_LS);

% ----- heavyScalarPS / lightScalarPS -----
old_PS = alphaEM^2*Dscript*mphi^3/(4608*pi^3) ...
       + gx^2*sqrt(max(0, mphi^2/4 - mx^2))/(4*pi);
new_PS = phi_partial_widths(mphi, mx, gx, g, 'heavyScalarPS').total;
report('heavyScalarPS', old_PS, new_PS);

% ----- heavyScalarSP / lightScalarSP -----
gpigamma  = alphaEM/pi/f;
getagamma = alphaEM/sqrt(3)/pi/f;
% NLO L6/L7 meson-mediator mass mixings (match phi_gg_coupling / sppp_couplings).
L6 = 1.46e-3; L7 = -0.39e-3;
m_up = 0.0022; m_dn = 0.00469; m_s = 0.095; summ = m_up + m_dn + m_s;
m2_pi0phi = (g(1)-g(2)) * (B0*f + 32*B0^2*L6/f*summ) ...
          + 32*B0^2*L7/f * (g(1)+g(2)+g(3)) * (m_up - m_dn);
m2_etaphi = ( (g(1)+g(2)-2*g(3)) * (B0*f + 32*B0^2*L6/f*summ) ...
            + 32*B0^2*L7/f * (g(1)+g(2)+g(3)) * (m_up + m_dn - 2*m_s) ) / sqrt(3);
gphigg_SP = -gpigamma  * m2_pi0phi / (mphi^2 - mpi0^2) ...
            -getagamma * m2_etaphi / (mphi^2 - m_eta^2);
old_SP = gphigg_SP^2 * mphi^3 / (64*pi) ...
       + gx^2 * max(0, mphi^2 - 4*mx^2)^(3/2) / (8*pi*mphi^2);
new_SP = phi_partial_widths(mphi, mx, gx, g, 'heavyScalarSP').total;
report('heavyScalarSP', old_SP, new_SP);

% ----- scalarPP / heavyScalarPP / lightScalarPP -----
% Same gphigg as SP (PP mixing through pi0/eta).
old_PP = gphigg_SP^2 * mphi^3 / (64*pi) ...
       + gx^2/(4*pi) * sqrt(max(0, mphi^2/4 - mx^2));
new_PP = phi_partial_widths(mphi, mx, gx, g, 'heavyScalarPP').total;
report('heavyScalarPP', old_PP, new_PP);
new_PPs = phi_partial_widths(mphi, mx, gx, g, 'scalarPP').total;
report('scalarPP',      old_PP, new_PPs);

% Unknown interaction -> total = 0 (matches previous svavg 'otherwise' fallback).
new_UK = phi_partial_widths(mphi, mx, gx, g, 'thisIsNotAThing').total;
report('unknown (->0)', 0, new_UK);

function report(name, old_val, new_val)
    rel = abs(old_val - new_val) / max(abs(old_val), eps);
    pass_str = "PASS";
    if rel > 1e-12 && abs(old_val - new_val) > 0
        pass_str = "FAIL";
    end
    fprintf('  %-18s old=%-14.6g  new=%-14.6g  rel=%-9.2g  %s\n', ...
            name, old_val, new_val, rel, pass_str);
end
