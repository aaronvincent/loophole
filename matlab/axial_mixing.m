function mjA = axial_mixing(g)
% axial_mixing  A'-meson mixing masses [m_{pi0 A'}, m_{eta A'}, m_{eta' A'}] [GeV].
%
%   mjA = axial_mixing(g)     g = [a_u a_d a_s]
%
% Single source of truth for axial_couplings.m (lambda_{j chi} = lambda_a
% m_{jA'}/m_A'^2) and apipi_phi_rate.m (the A' pi -> pi pi current).
%
% NOTE (unresolved sign): the NLO L4 term carries a "1 +/- 16 B0 L4/F^2 sum_i m_i"
% in the source expressions; the sign is not fixed by the note.  We use "+"
% here.  See docs/axialvector.md -- flip `corr4` below to test.

B0    = 0.13957^2/(0.0022 + 0.00469);   % GeV, = mpi+^2/(mu+md) ~ 2.827 (was 2.6)
F     = 0.0924;    % GeV, pion decay constant
L4    = 1.87e-3;   % NLO chiral low-energy constants
L5    = 1.22e-3;
m_up  = 0.0022;    % GeV
m_dn  = 0.00469;
m_s   = 0.095;

summ  = m_up + m_dn + m_s;
corr4 = 1 + 16*B0*L4/F^2 * summ;        % "+" branch (sign ambiguity, see above)

au = g(1); ad = g(2); as = g(3);
mjA = [ F        * ( (au - ad)        *corr4 + 16*B0*L5*(au*m_up - ad*m_dn)            /F^2 ), ...
        sqrt(3)*F/3 * ( (au + ad - 2*as)*corr4 + 16*B0*L5*(au*m_up + ad*m_dn - 2*as*m_s)/F^2 ), ...
        sqrt(6)*F/3 * ( (au + ad + as)  *corr4 + 16*B0*L5*(au*m_up + ad*m_dn + as*m_s)  /F^2 ) ];
end
