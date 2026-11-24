function sig = xxggxsec(mx,mphi,gx,g,s,interaction,parts)
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
% parts (optional, default 'all') selects final states:
%   'all'  every open channel
%   'sm'   SM final states only: gamma gamma + meson pairs
%   'had'  meson pairs only (pi pi, K K; SP/PP/axialVector jj)
%   'dark' the dark mediator pair only (= 'all' - 'sm')
% 'sm' leaves out the dark pairs: phi phi, A'A' and, for SP/PP, the phi'phi'
% term of the jj sum (index 2 in sppp_couplings' ordering).  It is the
% detailed-balance partner of chi production from the equilibrium SM bath.
% 'dark' feeds freeze_coupled's chi chibar <-> phi phi rates (sv_darkpair_grid).
if nargin < 7 || isempty(parts), parts = 'all'; end
use_dark = any(strcmp(parts, {'all', 'dark'}));
use_gg   = ~any(strcmp(parts, {'had', 'dark'}));
use_had  = ~strcmp(parts, 'dark');
keep_j   = @(jj) (jj ~= 2 && use_had) || (jj == 2 && use_dark);

if s > 4*mx^2
    if strcmp(interaction,'heavyScalar') || strcmp(interaction,'lightScalar')
        % Breit-Wigner total width of phi (regulates s-channel resonance).
        % heavy / light alias: the BW propagator (s - mphi^2)^2 + (mphi*Gphi)^2
        % handles both s << mphi^2 and s ~ mphi^2; no separate code path needed.
        Gphi = phi_partial_widths(mphi, mx, gx, g, interaction).total;

        % chi chi -> phi phi  (t-channel; open for s > 4*mphi^2).
        sig = zeros(size(s));
    % s > 4*mphi^2 is the only gate: it is the physical threshold, and s ranges
    % over [4*mx^2, ...) in a thermal bath, so the channel is reachable for
    % mx < mphi at a cost of exp(-2(mphi-mx)/T).  The integration range crosses
    % the branch point at s = 4*mphi^2 when mx < mphi; svavg.m passes 4*mphi^2 as
    % a continuum waypoint so the adaptive rule resolves the cusp.
        mask_pp = (s > 4*mphi^2);
        if use_dark && any(mask_pp)
            sp    = s(mask_pp);
            kin   = sqrt((sp - 4*mx^2).*(sp - 4*mphi^2));
            A     = -32*mx^4 + 16*mx^2*(sp - mphi^2) + 6*mphi^4 - 4*mphi^2*sp + sp.^2;
            atarg = sqrt((sp - 4*mx^2).*(4*mphi^2 - sp)) ./ (2*mphi^2 - sp);
            t1    = 2*(mx^2*(sp - 4*mphi^2) + mphi^4) .* A .* atan(atarg) ...
                    ./ ((sp - 4*mx^2).^(3/2) .* (2*mphi^2 - sp) .* sqrt(4*mphi^2 - sp));
            t2    = (16*mx^4 + 2*mx^2*(sp - 8*mphi^2) + 3*mphi^4) ./ (4*mx^2 - sp);
            % /32pi (was /16pi): identical-particle S = 1/2 for the phi phi
            % final state.  Verified against a direct tree-level trace
            % calculation, which reproduces this form exactly.
            sig(mask_pp) = real(gx^4 * kin .* (t1 + t2) ...
                           ./ (32*pi * sp .* (mx^2*(sp - 4*mphi^2) + mphi^4)));
        end

        % chi chi -> gamma gamma  (s-channel phi, Breit-Wigner regulated)
        % /9216 (was /4608): the identical-photon S_gg = 1/2 belongs here, not
        % only in the prod/ann ratio.  This is now exactly
        % pi s Gamma_xx Gamma_gg / (p_chi^2 [(s-mphi^2)^2 + mphi^2 Gphi^2])
        % with the same Gamma(phi->gg) = D alpha^2 mphi^3 / 4608 pi^3 that
        % phi_partial_widths.m uses.
        prop2 = (s - mphi^2).^2 + mphi^2*Gphi^2;
        sig = sig + use_gg*gx^2*alphaEM^2*Dscript/pi^3/9216 ...
              .* s.*(s - 4*mx^2) ./ prop2 .* sqrt(s./(s - 4*mx^2));

        % chi chi -> meson pairs (meson_pairs.m): S mu^2 times the pi+pi-
        % form, with sqrt(s - 4 m^2) -> sqrt(lambda(s, mA^2, mB^2)/s)
        if use_had
            P = meson_pairs(g);
            for j = 1:numel(P)
                mask = (s > (P(j).mA + P(j).mB)^2);
                if any(mask)
                    sp  = s(mask);
                    tp  = sqrt((sp - (P(j).mA+P(j).mB)^2).*(sp - (P(j).mA-P(j).mB)^2)./sp);
                    sig(mask) = sig(mask) + P(j).S*P(j).mu^2*gx^2/(32*pi) ...
                                .* tp.*sqrt(sp - 4*mx^2) ...
                                ./ (sp.*((sp - mphi^2).^2 + mphi^2*Gphi^2));
                end
            end
        end
    elseif strcmp(interaction,'heavyScalarPS') || strcmp(interaction,'lightScalarPS')
        % Dirac fermion DM, scalar mediator, pseudoscalar coupling igx*phi*chibar*gamma5*chi.
        % chi and chibar are distinct; no identical-particle factor.
        % Gamma(phi->chibar chi): pseudoscalar Dirac vertex gives sqrt(mphi^2/4 - mx^2)/(4*pi).
        Gphi = phi_partial_widths(mphi, mx, gx, g, interaction).total;
        prop2 = (s - mphi^2).^2 + mphi^2*Gphi^2;

        % chibar chi -> phi phi  (t-channel; open for s > 4*mphi^2).
        sig = zeros(size(s));
        mask_pp = (s > 4*mphi^2);
        if use_dark && any(mask_pp)
            % Positive and ~v_rel in the non-relativistic regime: a pseudoscalar
            % chi vertex has no s-wave -- at threshold the two diagrams sum to
            % vbar(pslash_1 + pslash_2)u = 0.  Matches a direct tree-level trace
            % calculation to all digits.
            % The 1/16pi already carries the identical-particle S = 1/2.
            sp    = s(mask_pp);
            kin   = sqrt((sp - 4*mx^2).*(sp - 4*mphi^2));
            t1    = -sqrt((sp - 4*mphi^2)./(sp - 4*mx^2)) ...
                    .* (2*mx^2*(sp - 4*mphi^2) + 3*mphi^4) ...
                    ./ (2*(mx^2*(sp - 4*mphi^2) + mphi^4));
            t2    = (6*mphi^4 - 4*mphi^2*sp + sp.^2) ...
                    ./ ((sp - 4*mx^2).*(sp - 2*mphi^2)) .* atanh(kin ./ (sp - 2*mphi^2));
            sig(mask_pp) = gx^4*(t1 + t2) ./ (16*pi*sp);
        end

        % chibar chi -> gamma gamma  (s-channel phi, Breit-Wigner regulated)
        % /9216 (was /4608): identical-photon S_gg = 1/2, as in the
        % heavyScalar branch above.
        sig = sig + use_gg*gx^2*alphaEM^2*Dscript/(9216*pi^3) .* s.^2 ./ prop2 ...
              .* sqrt(s./(s - 4*mx^2));

        % chibar chi -> meson pairs (meson_pairs.m)
        if use_had
            P = meson_pairs(g);
            for j = 1:numel(P)
                mask = (s > (P(j).mA + P(j).mB)^2);
                if any(mask)
                    sp  = s(mask);
                    tp  = sqrt((sp - (P(j).mA+P(j).mB)^2).*(sp - (P(j).mA-P(j).mB)^2)./sp);
                    sig(mask) = sig(mask) + P(j).S*P(j).mu^2*gx^2/(32*pi) ...
                                .* tp ./ (sqrt(sp - 4*mx^2) .* prop2(mask));
                end
            end
        end

    elseif strcmp(interaction,'heavyScalarSP') || strcmp(interaction,'lightScalarSP')
        % scalar DM, pseudoscalar mediator; annihilation is fully relativistic in both cases
        % Couplings/widths from sppp_couplings: {pi0, phi', eta} generically,
        % or {eta', phi'} for g ~ [1 1 1].  Matches ggxxxsec.m SP branch.
        [mj_arr, gcx_arr, ggx_arr, Gamma_arr] = sppp_couplings(mphi, mx, gx, g, interaction);

        % chi chibar -> gamma gamma -- BW-regulated complex propagators
        % match ggxxxsec.m SP branch (detailed balance preserved).
        sum_prop = zeros(size(s));
        for jj = 1:numel(mj_arr)
            sum_prop = sum_prop + gcx_arr(jj)*ggx_arr(jj) ...
                       ./ (mj_arr(jj)^2 - s - 1i*mj_arr(jj)*Gamma_arr(jj));
        end
        % /128pi (was /64pi): identical-photon S_gg = 1/2, as in the SS/PS
        % branches above.
        sig = use_gg * s.^(3/2) .* sqrt(s - 4*mx^2) / (128*pi) .* abs(sum_prop).^2;

        % chi chibar -> jj  (t-channel; sum over the same poles: pi0/phi'/eta,
        % or eta'/phi' for g ~ [1 1 1]).  s > 4*mj^2 is the only gate: it is
        % the physical threshold, and s runs over [4*mx^2, ...) in a thermal
        % bath, so the channel is reachable for mx < mj at a cost of
        % exp(-2(mj-mx)/T).  svavg.m passes 4*mj^2 as a continuum waypoint
        % (extra_wp covers pi0/eta/eta' and mphi = m_phi') so the adaptive
        % rule resolves the threshold cusp.
        for jj = 1:numel(mj_arr)
            mj  = mj_arr(jj);
            gcx = gcx_arr(jj);
            mask = (s > 4*mj^2);
            if keep_j(jj) && any(mask)
                sp    = s(mask);
                kin   = sqrt((sp - 4*mx^2).*(sp - 4*mj^2));
                A     = -32*mx^4 + 16*mx^2*(sp - mj^2) + 6*mj^4 - 4*mj^2*sp + sp.^2;
                t1    = 2*(mx^2*(sp - 4*mj^2) + mj^4) .* A .* atanh(kin ./ (2*mj^2 - sp)) ...
                        ./ ((sp - 4*mx^2).^(3/2) .* (2*mj^2 - sp) .* sqrt(sp - 4*mj^2));
                t2    = (16*mx^4 + 2*mx^2*(sp - 8*mj^2) + 3*mj^4) ./ (4*mx^2 - sp);
                % /32pi (was /16pi): identical-particle S = 1/2 for the jj
                % final state, as in the heavyScalar phi-phi branch.
                sig(mask) = sig(mask) + real(gcx^4 * kin .* (t1 + t2) ...
                            ./ (32*pi * sp .* (mx^2*(sp - 4*mj^2) + mj^4)));
            end
        end

        % chi chi -> 3pi (threepi_weight.m); Gamma_arr(2) is the mediator's width.
        if use_had
            sig = sig + threepi_sigma(s, mx, mphi, gx, Gamma_arr(2), g, true);
        end


    elseif strcmp(interaction,'scalarPP') || strcmp(interaction,'heavyScalarPP') || strcmp(interaction,'lightScalarPP')
        % chi chibar -> all, via pseudoscalar DM + pseudoscalar mediator.
        % Couplings/widths from sppp_couplings: {pi0, phi', eta} generically,
        % or {eta', phi'} for g ~ [1 1 1].
        [mj_arr, gcx_arr, ggx_arr, Gamma_arr] = sppp_couplings(mphi, mx, gx, g, interaction);

        % gamma gamma channel (complex propagators regulate on-shell resonances)
        % Loop over mediators so s can be a vector (as passed by integral())
        sum_prop = zeros(size(s));
        for jj = 1:numel(mj_arr)
            sum_prop = sum_prop + gcx_arr(jj)*ggx_arr(jj) ...
                       ./ (mj_arr(jj)^2 - s - 1i*mj_arr(jj)*Gamma_arr(jj));
        end
        % /128pi (was /64pi): identical-photon S_gg = 1/2, as in the SS/PS
        % branches above.
        sig = use_gg * s.^2/(128*pi) .* abs(sum_prop).^2 .* sqrt(s./(s - 4*mx^2));

        % meson pair channels (sum over kinematically open j: pi0/phi'/eta,
        % or eta'/phi' for g ~ [1 1 1]).  The atanh(arg)/(sp-2mj^2) pole sits
        % at sp = 2*mj^2, below the sp > 4*mj^2 threshold, so it is never in
        % range; the closed form is finite and vanishes as beta at threshold.
        for jj = 1:numel(mj_arr)
            mj  = mj_arr(jj);
            gcx = gcx_arr(jj);
            mask = (s > 4*mj^2);
            if keep_j(jj) && any(mask)
                % Same form as the heavyScalarPS phi-phi branch (positive,
                % ~v_rel); verified against a direct trace calculation.
                % The 1/16pi already carries the identical-particle S = 1/2.
                sm   = s(mask);
                kin  = sqrt((sm - 4*mx^2).*(sm - 4*mj^2));
                num1 = -sqrt((sm - 4*mj^2)./(sm - 4*mx^2)) ...
                       .* (2*mx^2*(sm - 4*mj^2) + 3*mj^4) ...
                       ./ (2*(mx^2*(sm - 4*mj^2) + mj^4));
                num2 = (6*mj^4 - 4*mj^2*sm + sm.^2) ...
                       ./ ((sm - 4*mx^2).*(sm - 2*mj^2)) .* atanh(kin ./ (sm - 2*mj^2));
                sig(mask) = sig(mask) + gcx^4*(num1 + num2) ./ (16*pi*sm);
            end
        end

        % chi chi -> 3pi (threepi_weight.m); Gamma_arr(2) is the mediator's width.
        if use_had
            sig = sig + threepi_sigma(s, mx, mphi, gx, Gamma_arr(2), g, false);
        end

    elseif strcmp(interaction,'axialVector')
        % chi chibar -> all for the axial-vector mediator A'.  Three pieces
        % (full chi chibar -> all expression):
        %   gamma gamma : 3 SM-meson s-channel mixing poles (pi0/eta/eta'),
        %                 BW-regulated; this piece is the detailed-balance
        %                 partner of ggxxxsec.m's gg->xx (factor-2 S_gg there).
        %   A'A'        : t/u-channel chi exchange (open for s > 4 m_A'^2).
        %   jj          : meson pairs (j = pi0/eta/eta', open for s > 4 mj^2).
        % A'A' and jj reuse the heavyScalar-phiphi / scalarPP-meson-pair
        % closed forms' structure.  Both are gated on the physical threshold
        % only (s > 4 m^2); sub-threshold mx costs exp(-2(m-mx)/T), which is
        % not zero for mx within a few T of m.
        [mj, lam_chi, lam_gg, Gam_j] = axial_couplings(mphi, mx, gx, g);

        % --- gamma gamma --------------------------------------------------
        sum_prop = zeros(size(s));
        for jj = 1:3
            sum_prop = sum_prop + lam_chi(jj)*lam_gg(jj) ...
                       ./ (s - mj(jj)^2 + 1i*mj(jj)*Gam_j(jj));
        end
        % /32pi (was /16pi): identical-photon S_gg = 1/2, as in the 
        % other four
        % families.
        sig = use_gg * mx^2 * s.^2/(32*pi) .* abs(sum_prop).^2 .* sqrt(s./(s - 4*mx^2));

        % --- A'A' (t/u-channel chi exchange) ------------------------------
        mask_pp = (s > 4*mphi^2);
        if use_dark && any(mask_pp)
            sp    = s(mask_pp);
            sqneg = sqrt((4*mphi^2 - sp).*(sp - 4*mx^2));   % imaginary for sp>4mphi^2
            atrm  = atan(sqneg ./ (2*mphi^2 - sp));
            base  = mphi^4 - 4*mphi^2*mx^2 + mx^2*sp;
            poly  = 16*mphi^8 - 16*mphi^6.*(4*mx^2 + sp) + 2*mphi^4.*sp.*(8*mx^2 + 5*sp) ...
                    - 4*mphi^2.*sp.*(-16*mx^4 - 2*mx^2*sp + sp.^2) - 16*mx^4*sp.^2 + sp.^4;
            num   = -8*mphi^8 + 8*mphi^6.*(2*mx^2 + sp) ...
                    + mphi^4.*(64*mx^4 - 28*mx^2*sp - 3*sp.^2) ...
                    + 4*mphi^2*mx^2.*sp.*(3*sp - 16*mx^2) ...
                    + 2*base.*poly.*atrm ./ ((2*mphi^2 - sp).*sqneg) ...
                    - 2*mx^2*sp.^2.*(sp - 4*mx^2);
            % /32pi (was /16pi): identical-particle S = 1/2 for A'A'.  NOTE:
            % unlike the other pair channels this one is NOT independently
            % verified -- my trace check of the axial channels reproduced
            % neither normalization.  Applied for
            % consistency with the rest of the draft's symmetry-factor pass.
            sig(mask_pp) = sig(mask_pp) + real( gx^4 * sqrt((sp - 4*mphi^2)./(sp - 4*mx^2)) ...
                           .* num ./ (32*pi*mphi^4.*sp.*base) );
        end

        % --- meson pairs jj (j = pi0, eta, eta') --------------------------
        for jj = 1:3
            mjj = mj(jj);
            mask = (s > 4*mjj^2);
            if use_had && any(mask)
                sp  = s(mask);
                kin = sqrt((sp - 4*mjj^2).*(sp - 4*mx^2));
                den = mjj^4 + mx^2*(sp - 4*mjj^2);
                t1  = kin .* (mx^2*(4*mjj^4 + 8*mjj^2*sp - 2*sp.^2) - 3*mjj^4*sp) ./ den;
                t2  = 2*(sp.*(6*mjj^4 - 4*mjj^2*sp + sp.^2) ...
                      + 4*mx^2*(2*mjj^4 - 4*mjj^2*sp + sp.^2)) ...
                      .* atanh(kin ./ (2*mjj^2 - sp)) ./ (2*mjj^2 - sp);
                % /32pi (was /16pi): identical-particle S = 1/2 for jj.  Same
                % caveat as the A'A' channel above -- not independently
                % verified.
                sig(mask) = sig(mask) + real( lam_chi(jj)^4 * (t1 + t2) ./ (32*pi*sp) );
            end
        end

        % --- chi chi -> 3pi, pi+pi-A' (axial3pi_sigma.m); J = 1 part of 3pi
        % resonant at m_A'^2
        if use_had && (threepi_enabled() || pipiphi_enabled())
            GA  = phi_partial_widths(mphi, mx, gx, g, interaction).total;
            sig = sig + axial3pi_sigma(s, mx, mphi, gx, GA, g);
        end

    end
else
    sig = 0;
end