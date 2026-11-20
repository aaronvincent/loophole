function verify_dmnucxsec_fullrange()
% Smoke test: DMnucxsec.m closed forms vs brute-force quadrature of the supplied
% dsigma/dOmega for every vertex structure (SS/SP/PS/PP/axialVector), across the
% full mediator-mass range. sigma = int_{-T}^0 (pi/(mu^2 vrel^2)) dsigma/dOmega dt,
% T = 4 mu^2 vrel^2.  Reports max relative error (should be ~1e-10 or better).
    mN = 0.938; vrel = 220/3e5;
    Du = 0.842; Dd = -0.427; Ds = -0.085;
    g  = [0.7, -0.3, 0.15];                 % non-degenerate so htildeP ~= 0
    gx = 0.8;
    % nucleon couplings (mirror DMnucxsec.m)
    fpu=0.02; fpd=0.026; fps=0.118; mu_=0.0022; md_=0.00469; ms_=0.095;
    hp = mN*(g(1)*fpu/mu_ + g(2)*fpd/md_ + g(3)*fps/ms_);
    invm = 1/mu_+1/md_+1/ms_; gpr = (g(1)/mu_+g(2)/md_+g(3)/ms_)/invm;
    Dnu=0.84; Dnd=-0.44; Dns=-0.03;
    htp = mN*(Dnu/mu_*(g(1)-gpr)+Dnd/md_*(g(2)-gpr)+Dns/ms_*(g(3)-gpr));
    aN  = g(1)*Du+g(2)*Dd+g(3)*Ds;

    % FULL angular-integrated form (the unprefixed 'scalar*' names + axialVector).
    its = {'scalar','scalarSP','scalarPS','scalarPP','axialVector'};
    mxs = [1e-3 1e-2 0.1 1];
    mps = [1e-4 1e-2 0.1 1 10 100];

    fprintf('%-16s %10s %10s %14s %14s %10s\n','interaction','mx','mphi','closed','quad','rel.err');
    worst = 0;
    for ii = 1:numel(its)
        it = its{ii};
        for mx = mxs
            mu = mN*mx/(mN+mx); sCM = (mx+mN)^2; T = 4*mu^2*vrel^2;
            for mphi = mps
                sig_c = DMnucxsec(mx,mphi,gx,g,it);
                % dsigma/dOmega(t) for this vertex
                switch it
                  case 'scalar'        % SS
                    dsdO = @(t) gx^2*hp^2/(64*pi^2).*(4*mx^2-t).*(4*mN^2-t)./(sCM*(mphi^2-t).^2);
                  case 'scalarSP'      % SP: -htp^2 t(4mx^2-t)/(s(mphi^2-t)^2), * gx^2
                    dsdO = @(t) -gx^2*htp^2/(64*pi^2).*t.*(4*mx^2-t)./(sCM*(mphi^2-t).^2);
                  case 'scalarPS'      % PS: -hp^2 t(4mN^2-t)/(s(mphi^2-t)^2), * gx^2
                    dsdO = @(t) -gx^2*hp^2/(64*pi^2).*t.*(4*mN^2-t)./(sCM*(mphi^2-t).^2);
                  case 'scalarPP'      % PP: htp^2 t^2/(s(mphi^2-t)^2), * gx^2
                    dsdO = @(t) gx^2*htp^2/(64*pi^2).*t.^2./(sCM*(mphi^2-t).^2);
                  case 'axialVector'
                    G = @(t) 2*mN^4*mphi^4 ...
                        + 4*mN^2.*(mx^2.*(5*mphi^4-4*mphi^2.*t+2*t.^2) - mphi^4.*(sCM+t)) ...
                        + mphi^4.*(2*mx^4 - 4*mx^2.*(sCM+t) + 2*sCM^2 + 2*sCM.*t + t.^2);
                    dsdO = @(t) gx^2*aN^2./(32*pi^2*mphi^4*sCM*(mphi^2-t).^2).*G(t);
                end
                dsdt = @(t) (pi/(mu^2*vrel^2)).*dsdO(t);
                sig_q = integral(dsdt, -T, 0, 'RelTol',1e-13,'AbsTol',0);
                re = abs(sig_c - sig_q)/max(abs(sig_q),realmin);
                worst = max(worst, re);
                fprintf('%-16s %10.1e %10.1e %14.6e %14.6e %10.1e\n', ...
                        it, mx, mphi, sig_c, sig_q, re);
            end
        end
    end
    fprintf('\nWORST relative error (full form vs quadrature): %.2e\n', worst);

    % --- heavy-mediator approximation: DMnucxsec(...,'heavyScalar*') is the leading
    %     term of the full form and must converge to it as mphi grows.  Check the
    %     heavy/full pairs at a large mphi where the limit is excellent.
    pairs = {'heavyScalar','scalar'; 'heavyScalarPS','scalarPS'; ...
             'heavyScalarSP','scalarSP'; 'heavyScalarPP','scalarPP'};
    mphi_heavy = 1e3;
    fprintf('\nHeavy-limit check at mphi=%g:\n', mphi_heavy);
    fprintf('%-16s %10s %14s %14s %10s\n','heavy name','mx','heavy','full','rel.err');
    worst_heavy = 0;
    for pp = 1:size(pairs,1)
        for mx = mxs
            sig_h = DMnucxsec(mx,mphi_heavy,gx,g,pairs{pp,1});
            sig_f = DMnucxsec(mx,mphi_heavy,gx,g,pairs{pp,2});
            re = abs(sig_h - sig_f)/max(abs(sig_f),realmin);
            worst_heavy = max(worst_heavy, re);
            fprintf('%-16s %10.1e %14.6e %14.6e %10.1e\n', pairs{pp,1}, mx, sig_h, sig_f, re);
        end
    end
    fprintf('\nWORST heavy-limit rel.err (heavy approx vs full): %.2e\n', worst_heavy);

    if worst < 1e-8 && worst_heavy < 1e-4
        fprintf('PASS: full form matches quadrature and heavy approx converges to full.\n');
    else
        fprintf('FAIL: full-vs-quad %.2e (<1e-8?) or heavy-vs-full %.2e (<1e-4?).\n', ...
                worst, worst_heavy);
    end
end
