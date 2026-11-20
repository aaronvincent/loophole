function sig = threepi_sigma(s, mx, mphi, gx, Gphi, g, isSP)
% threepi_sigma  chi chibar -> 3pi + pi pi phi for the SP/PP sector.
% s vector [GeV^2], all s > 4 mx^2.  isSP selects the scalar chi vertex.  The
% amplitudes and the on/off switches are in threepi_weight.m, pipiphi_weight.m.
sig  = (threepi_weight(s, g) ./ ((s - mphi^2).^2 + mphi^2*Gphi^2) ...
        + pipiphi_weight(s, mphi, g)) * gx^2 / (128*pi^3) ./ sqrt(1 - 4*mx^2./s);
if isSP
    sig = sig .* (s - 4*mx^2)./s;
end
end
