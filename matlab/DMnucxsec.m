function sig = DMnucxsec(mx,mphi,gx,g,interaction)
%Scalar, heavy mediator
%ignoring heavy quarks
%gx is the dark matter coupling
%g is a 3x1 of the u, d, s couplings
%quark composition https://arxiv.org/pdf/1202.1292 table II
%proton 
mN = 0.938;   % nucleon mass [GeV]
% Effective nucleon couplings: scalar hp (SS/PS), pseudoscalar htildeP (SP/PP),
% axial a_N (axialVector).  Single source of truth -> nucleon_couplings.m.
[hp, htildeP, a_N] = nucleon_couplings(g);

mu = mN*mx/(mN+mx);
vrel = 220/3e5;

% --- sigma_chiN regime is selected by the interaction NAME -----------------
% The suffix sets the vertex structure (SS/PS/SP/PP); the prefix selects the regime:
%   'heavyScalar*' -> heavy-mediator analytic limit (1/mphi^4 contact);
%   'lightScalar*' -> light-mediator closed forms;
%   'scalar*'      -> the FULL angular-integrated form, valid across the whole
%                     mediator-mass range.
% 'axialVector' is always full (its smooth heavy<->light crossover is built in).
% full_sigma_name.m maps a heavy*/light* name to its 'scalar*' full sibling so the
% canonical-observable pipeline (model_observables, compare_*) reports the full form;
% the heavy*/light* approximations are for explicit DMnucxsec(...,'heavy/lightScalar*').
scalar_family = {'scalar','heavyScalar','lightScalar', ...
                 'scalarPS','heavyScalarPS','lightScalarPS', ...
                 'scalarSP','heavyScalarSP','lightScalarSP', ...
                 'scalarPP','heavyScalarPP','lightScalarPP'};
if any(strcmp(interaction, scalar_family))
    % Vertex structure: a pseudoscalar chi vertex (PS/PP) zeroes its mass term;
    % a pseudoscalar nucleon vertex (SP/PP) zeroes its term AND swaps hp -> htildeP.
    isPSchi = any(strcmp(interaction,{'heavyScalarPS','lightScalarPS','scalarPS', ...
                                      'heavyScalarPP','lightScalarPP','scalarPP'}));
    isPSnuc = any(strcmp(interaction,{'heavyScalarSP','lightScalarSP','scalarSP', ...
                                      'heavyScalarPP','lightScalarPP','scalarPP'}));
    if isPSnuc, hN = htildeP; else, hN = hp; end       % pseudoscalar nucleon -> htildeP
    sCM = (mx+mN)^2;
    if strncmp(interaction,'heavy',5)
        % Heavy-mediator limit (mphi^2 >> T, mx^2, mN^2): the leading term of the
        % full angular integral (verify_dmnucxsec_fullrange.m checks it against the
        % full form at large mphi). SS s-wave; PS/SP ~vrel^2 (p-wave); PP ~vrel^4.
        %   SS -> gx^2 hp^2      mu^2/(pi mphi^4)
        %   PS -> gx^2 hp^2      mN^2 mu^2 vrel^2/(2 pi s mphi^4)
        %   SP -> gx^2 htildeP^2 mx^2 mu^2 vrel^2/(2 pi s mphi^4)
        %   PP -> gx^2 htildeP^2      mu^4 vrel^4/(3 pi s mphi^4)
        if     ~isPSchi && ~isPSnuc        % SS
            sig = gx^2*hN^2*mu^2/(pi*mphi^4);
        elseif  isPSchi && ~isPSnuc        % PS
            sig = gx^2*hN^2*mN^2*mu^2*vrel^2/(2*pi*sCM*mphi^4);
        elseif ~isPSchi &&  isPSnuc        % SP
            sig = gx^2*hN^2*mx^2*mu^2*vrel^2/(2*pi*sCM*mphi^4);
        else                               % PP
            sig = gx^2*hN^2*mu^4*vrel^4/(3*pi*sCM*mphi^4);
        end
    elseif strncmp(interaction,'light',5)
        % Light-mediator closed forms. SS keeps the propagator with an IR regulator
        % q0=mx*vrel; PS/SP/PP are mphi-independent.
        %   SS -> gx^2 hp^2      mu^2/(pi (q0^2+mphi^2)^2)
        %   PS -> gx^2 hp^2      /(8 pi mx^2 vrel^2)
        %   SP -> gx^2 htildeP^2 /(8 pi mN^2 vrel^2)
        %   PP -> gx^2 htildeP^2 /(16 pi (mx+mN)^2)
        q0sq = mx^2*vrel^2;
        if     ~isPSchi && ~isPSnuc        % SS
            sig = gx^2*hN^2*mu^2/(pi*(q0sq+mphi^2)^2);
        elseif  isPSchi && ~isPSnuc        % PS
            sig = gx^2*hN^2/(8*pi*mx^2*vrel^2);
        elseif ~isPSchi &&  isPSnuc        % SP
            sig = gx^2*hN^2/(8*pi*mN^2*vrel^2);
        else                               % PP
            sig = gx^2*hN^2/(16*pi*sCM);
        end
    else
        % Full angular-integrated form, valid across the whole mediator-mass range.
        % Integrate dsigma/dOmega = gx^2 hN^2/(64 pi^2) (a-t)(b-t)/(s (mphi^2-t)^2)
        % over solid angle at reference vrel: with t=-q^2, p_cm=mu*vrel,
        % dsigma/dt=(pi/p_cm^2)dsigma/dOmega, t in [-T,0], T=4 mu^2 vrel^2,
        %   sigma = gx^2 hN^2/(64 pi mu^2 vrel^2 s) * mphi^2 * J(x,al,be),
        %   x=T/mphi^2, al=a/mphi^2, be=b/mphi^2.  Each pseudoscalar vertex zeroes its
        %   own mass term: SS=(4mx^2,4mN^2), PS=(0,4mN^2), SP=(4mx^2,0), PP=(0,0).
        T  = 4*mu^2*vrel^2;
        if isPSchi, al = 0; else, al = 4*mx^2/mphi^2; end
        if isPSnuc, be = 0; else, be = 4*mN^2/mphi^2; end
        x  = T/mphi^2;
        sig = gx^2*hN^2/(64*pi*mu^2*vrel^2*sCM) * (mphi^2 * scalar_angint(x,al,be));
    end
elseif strcmp(interaction,'axialVector')
    % Axial-vector mediator A': spin-dependent DM-nucleon scattering, full
    % mediator-mass range.  a_N = sum_q a_q Delta_q^p (proton spin content; the
    % quark coupling g plays the role of a_q). Integrate the supplied
    %   dsigma/dOmega = lambda_a^2 a_N^2/(32 pi^2 m_A'^4 s (m_A'^2-t)^2) * G(s,t)
    % over solid angle (t=-q^2, p_cm=mu*vrel, t in [-T,0], T=4 mu^2 vrel^2):
    %   sigma = lambda_a^2 a_N^2/(32 pi mu^2 vrel^2 m_A'^4 s) * I.
    % Limits (both old forms recovered exactly, smooth crossover replaces the old
    % hard mphi>q0 switch): heavy (m_A'^2>>T) -> 3 lambda_a^2 a_N^2 mu^2/(pi m_A'^4)
    % (G0 bracket -> 24 mx^2 mN^2); light (m_A'->0) -> lambda_a^2 a_N^2 mu^2/(pi m_A'^4)
    % (only G2*m_A'^2*J0 survives).  a_N from nucleon_couplings (top of file).
    sCM = (mx+mN)^2;
    T   = 4*mu^2*vrel^2;
    x   = T/mphi^2;
    [J0,J1,J2] = jmoments(x);
    % G(s,t)=G0+G1 t+G2 t^2, collected from the supplied G(s,t) (mn=mN).
    G0 = mphi^4*( 2*mN^4 + 4*mN^2*(5*mx^2 - sCM) + 2*mx^4 - 4*mx^2*sCM + 2*sCM^2 );
    G1 = -16*mN^2*mx^2*mphi^2 + mphi^4*( 2*sCM - 4*mN^2 - 4*mx^2 );
    G2 = 8*mN^2*mx^2 + mphi^4;
    % I = int_{-T}^0 G/(mphi^2-t)^2 dt, with the J0,J1,J2 t-moments (x=T/mphi^2).
    Iax = G0*J2/mphi^2 - G1*J1 + G2*mphi^2*J0;
    sig = gx^2*a_N^2/(32*pi*mu^2*vrel^2*mphi^4*sCM) * Iax;
end

end

function J = scalar_angint(x,al,be)
% Dimensionless angular integral for the SS/SP/PS/PP DM-nucleon sigma (a pseudoscalar
% chi vertex sets al=0, a pseudoscalar nucleon vertex sets be=0):
%   J(x,al,be) = x + (al+be-2) ln(1+x) + (1-al-be+al*be) x/(1+x)
%              = J0 + (al+be)*J1 + al*be*J2,
% with J0,J1,J2 the cancellation-safe t-moments from jmoments(x).
    [J0,J1,J2] = jmoments(x);
    J = J0 + (al+be).*J1 + al.*be.*J2;
end

function [J0,J1,J2] = jmoments(x)
% Cancellation-safe moments of 1/(mphi^2-t)^2 over t in [-T,0] (x=T/mphi^2). With
% u=-t, w=u/mphi^2:
%   J0 = int_0^x w^2/(1+w)^2 dw = x-2ln(1+x)+x/(1+x)   (~x^3/3),
%   J1 = int_0^x w  /(1+w)^2 dw = ln(1+x)-x/(1+x)       (~x^2/2),
%   J2 = int_0^x      1/(1+w)^2 dw = x/(1+x).
% J0,J1 cancel to O(x^3)/O(x^2); evaluating them naively loses all precision for
% heavy mphi (x~1e-9), so use their convergent series below x=0.05. J2 is
% cancellation-free. Matches direct quadrature to ~5e-16. Shared by the scalar
% family (scalar_angint) and the axialVector DM-nucleon cross section.
    J2 = x./(1+x);
    if x < 0.05
        n  = 3:12;  J0 = sum((-1).^(n+1).*(n-2)./n .* x.^n);   % x^3/3 - x^4/2 + ...
        m  = 2:12;  J1 = sum((-1).^m    .*(m-1)./m .* x.^m);   % x^2/2 - 2x^3/3 + ...
    else
        L  = log1p(x);
        J0 = x - 2*L + J2;
        J1 = L - J2;
    end
end

