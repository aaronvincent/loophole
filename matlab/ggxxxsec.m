
function sig = ggxxxsec(mx,mphi,gx,g,s,interaction)
%gamma gamma -> x x cross section
%scalar cross section
%heavy mediator limit
mpi0 = 0.13498; mpi_pm = 0.13957; Gpi = 7.7e-9; alphaEM = 1/137;
m_eta = 0.547862;
mK= 0.493677;
mN = 0.938;


m_up = 0.0022; %GeV
m_dn = 0.00469;
m_s = 0.095;
B0 = 0.13957^2/(0.0022 + 0.00469); % GeV, = mpi+^2/(mu+md) ~ 2.827 (draft convention; was 2.6)
L4 = 1.87e-3;
L5 = 1.22e-3;
f = 0.0924;
fpu = 0.02;
fpd = 0.026;
fps = 0.118;
m_up = 0.0022; %GeV
m_dn = 0.00469;
m_s = 0.095;

Dscript = 0.5*(B0*(g(1)+g(2))/mpi_pm^2 + B0*(g(1)+g(3))/mK^2)^2;   % charged-meson loop (pi+, K+), low-energy limit A0 -> 1/3
if s > 4*mx^2
    if strcmp(interaction,'heavyScalar') || strcmp(interaction,'lightScalar')
        % Full s-channel scalar mediator with BW width.  Matches xxggxsec.m
        % so that detailed balance is preserved between gg->xx and xx->gg.
        % Single formula now covers both 'heavy' (s << mphi^2) and 'light'
        % (s ~ mphi^2) regimes; the BW propagator handles both.
        Gphi = phi_partial_widths(mphi, mx, gx, g, interaction).total;
        prop2 = (s - mphi^2).^2 + mphi^2*Gphi^2;
        % Source formula is dsigma/dOmega; multiply by 4*pi for total sigma.
        % Prefactor /4608 (not /9216) restores detailed balance with
        % xxggxsec.m: S_gamma_gamma = 1/2 in the chi chi -> gamma gamma
        % final state requires sigma(gg->xx)/sigma(xx->gg) = 2(s-4mx^2)/s,
        % see docs/cross_sections.md.
        sig = gx^2*alphaEM^2*Dscript/(4608*pi^3) ...
              .* s.*(s - 4*mx^2)./prop2 .* sqrt((s - 4*mx^2)./s);
    elseif strcmp(interaction,'heavyScalarPS') || strcmp(interaction,'lightScalarPS')
        % Pseudoscalar phi-chi-chi vertex i*g_chi*phi*chibar*gamma5*chi.
        % chi-chi vertex has no beta^2 = (s-4mx^2)/s suppression, unlike
        % the scalar case.  Phi width matches xxggxsec.m heavyScalarPS.
        Gphi = phi_partial_widths(mphi, mx, gx, g, interaction).total;
        prop2 = (s - mphi^2).^2 + mphi^2*Gphi^2;
        % Source formula is dsigma/dOmega; multiply by 4*pi for total sigma.
        % Prefactor /4608 (not /9216) restores detailed balance with
        % xxggxsec.m, exactly as in the SS branch above: S_gamma_gamma = 1/2
        % in the chi chi -> gamma gamma final state requires
        % sigma(gg->xx)/sigma(xx->gg) = 2(s-4mx^2)/s.
        sig = gx^2*alphaEM^2*Dscript/(4608*pi^3) ...
              .* s.^2 ./ prop2 .* sqrt((s - 4*mx^2)./s);

    elseif (strcmp(interaction,'lightScalarPP') | strcmp(interaction,'heavyScalarPP') | strcmp(interaction,'scalarPP'))
        % Pseudoscalar DM + pseudoscalar mediator.  Intermediate resonances
        % {pi0, phi', eta} generically, or the {eta', phi'} singlet set for
        % g ~ [1 1 1] (see sppp_couplings.m).  Couplings/widths/propagators
        % mirror xxggxsec.m scalarPP exactly so detailed balance holds.
        [mj_arr, gcx_arr, ggx_arr, Gam_arr] = sppp_couplings(mphi, mx, gx, g, interaction);
        sum_prop = zeros(size(s));
        for jj = 1:numel(mj_arr)
            sum_prop = sum_prop + gcx_arr(jj)*ggx_arr(jj) ...
                       ./ (mj_arr(jj)^2 - s - 1i*mj_arr(jj)*Gam_arr(jj));
        end
        % Source formula is dsigma/dOmega; multiply by 4*pi for total sigma.
        % /64pi (not /128pi) carries the identical-photon factor 2 relative to
        % the xx->gg piece of xxggxsec.m, as in the SS/axialVector branches.
        sig = s.^2/(64*pi) .* abs(sum_prop).^2 .* sqrt((s - 4*mx^2)./s);

    elseif strcmp(interaction,'heavyScalarSP') || strcmp(interaction,'lightScalarSP')
        % Scalar DM, pseudoscalar mediator (PP-mixing).  Intermediate
        % resonances {pi0, phi', eta} generically, or {eta', phi'} for
        % g ~ [1 1 1] (see sppp_couplings.m), each BW-regulated so the
        % propagator is integrable across the pole.  Matches xxggxsec.m
        % heavyScalarSP -- detailed balance preserved.
        [mj_arr, gcx_arr, ggx_arr, Gam_arr] = sppp_couplings(mphi, mx, gx, g, interaction);
        sum_prop = zeros(size(s));
        for jj = 1:numel(mj_arr)
            sum_prop = sum_prop + gcx_arr(jj)*ggx_arr(jj) ...
                       ./ (mj_arr(jj)^2 - s - 1i*mj_arr(jj)*Gam_arr(jj));
        end
        % Source formula is dsigma/dOmega; multiply by 4*pi for total sigma.
        % abs(...)^2 (not (...).^2) now that the propagators are complex.
        % /64pi (not /128pi) carries the identical-photon factor 2 relative to
        % the xx->gg piece of xxggxsec.m, as in the SS/axialVector branches.
        sig = sqrt(s).*(s - 4*mx^2).^(3/2)/(64*pi) .* abs(sum_prop).^2;

    elseif strcmp(interaction,'axialVector')
        % gamma gamma -> chi chibar (axial-vector mediator A'), via the three
        % SM-meson s-channel mixing poles pi0, eta, eta' (NOT an A' pole).
        % Detailed-balance crossing of the chi chibar -> gamma gamma piece in
        % xxggxsec.m: sigma(gg->xx) = 2(s-4mx^2)/s * sigma(xx->gg).  The factor
        % 2 is identical-photon statistics S_gg = 1/2 (same convention as the
        % heavyScalar /4608-vs-/9216 split, see docs/cross_sections.md); it turns the /32pi
        % of xx->gg into /16pi here with the beta velocity factor flipped.
        [mj, lam_chi, lam_gg, Gam_j] = axial_couplings(mphi, mx, gx, g);
        sum_prop = zeros(size(s));
        for jj = 1:3
            sum_prop = sum_prop + lam_chi(jj)*lam_gg(jj) ...
                       ./ (s - mj(jj)^2 + 1i*mj(jj)*Gam_j(jj));
        end
        sig = mx^2 * s.^2/(16*pi) .* abs(sum_prop).^2 .* sqrt((s - 4*mx^2)./s);

    end
else
    sig = 0;
end