function [I0, I1, I2] = dalitz_moments(s, m1, m2)
% dalitz_moments  I_n(s) = int dE Y(E) X(E)^n over the energy of the m1 particle.
%
% s is a column vector.  The other two particles have mass m2.
% Y = sqrt((E^2-m1^2)(1 - 4 m2^2/w)) is the E2 range at fixed E,
% w = s + m1^2 - 2 E sqrt(s), X = w - mpi0^2.  Y has square-root zeros at both
% ends, so t = (1 - cos th)/2 and the midpoint rule in th converge spectrally.
n  = 32;
th = pi*((1:n) - 0.5)/n;
t  = (1 - cos(th))/2;
dt = sin(th)/2*pi/n;
rs = sqrt(s);
Elo = m1*ones(size(s));
Ehi = (s + m1^2 - 4*m2^2)./(2*rs);
E  = Elo + (Ehi - Elo).*t;
w  = s + m1^2 - 2*E.*rs;
Y  = sqrt(max((E.^2 - m1^2).*(1 - 4*m2^2./w), 0));
X  = w - 0.13498^2;
wt = (Ehi - Elo).*dt.*Y;
I0 = sum(wt, 2);
I1 = sum(wt.*X, 2);
I2 = sum(wt.*X.^2, 2);
end
