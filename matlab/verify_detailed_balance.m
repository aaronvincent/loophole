% verify_detailed_balance.m
%
% Cross-check the gg->xx (FI pipeline) and xx->gg (loop pipeline) cross
% sections satisfy detailed balance for `heavyScalar` after the BW
% propagator fix in ggxxxsec.m.
%
% Two independent tests:
%   (A) Pointwise s ratio test:
%         R_sig(s) = sigma_gg->xx(s) / sigma_xx->gg(s)
%       Expected from crossing for gg<->xx-bar with Dirac chi:
%         R_expected = 2 * (s - 4 m_chi^2) / s
%       (factor of 2 = identical-photon-pair statistics in final state for
%        the chi chi-bar -> gamma gamma direction).
%
%   (B) Equilibrium collision-rate test:
%         G_FI(T) from ggxxrate.m should equal
%         n_chi_eq(T)^2 * <sigma v>_xx->gg(T) * (overall normalization)
%       up to a fixed prefactor that does not depend on T or mphi.

clear; close all;

interaction = 'heavyScalar';
mphi = 10;
mx   = 0.1;
gx   = 1;
g    = [1 1 1]/3;

% --- (A) pointwise ratio ----------------------------------------------------
s_test = logspace(log10(4*mx^2*1.001), log10(1e4*mphi^2), 60);
sig_gx = arrayfun(@(s) ggxxxsec(mx,mphi,gx,g,s,interaction), s_test);
sig_xg = arrayfun(@(s) xxggxsec(mx,mphi,gx,g,s,interaction), s_test);

% xxggxsec returns chi-chi -> all; isolate the gamma gamma piece by
% recomputing it standalone here so we don't mix in pi-pi / K-K / phi-phi.
mpi0 = 0.13498; mpi_pm = 0.13957; mK = 0.493677; alphaEM = 1/137; B0 = 0.13957^2/(0.0022 + 0.00469);
Dscript = 0.5*(B0*(g(1)+g(2))/mpi_pm^2 + B0*(g(1)+g(3))/mK^2)^2;   % charged-meson loop (pi+, K+), low-energy limit A0 -> 1/3
Gphi = phi_partial_widths(mphi, mx, gx, g, interaction).total;
prop2 = (s_test - mphi^2).^2 + mphi^2*Gphi^2;
sig_xg_gg_only = gx^2*alphaEM^2*Dscript/pi^3/9216 ...
                .* s_test.*(s_test - 4*mx^2)./prop2 ...
                .* sqrt(s_test./(s_test - 4*mx^2));

R_obs = sig_gx ./ sig_xg_gg_only;
R_exp = 2*(s_test - 4*mx^2)./s_test;

fprintf('\n=== (A) Pointwise sigma ratio (gg->xx / xx->gg, gg-only) ===\n');
fprintf('%12s %14s %14s %10s\n','sqrt(s) [GeV]','R_observed','R_expected','ratio');
idx = round(linspace(1,length(s_test),8));
for k = idx
    fprintf('%12.4g %14.4e %14.4e %10.4f\n', sqrt(s_test(k)), R_obs(k), R_exp(k), R_obs(k)/R_exp(k));
end

% --- (B) equilibrium-rate test ---------------------------------------------
T_vec = logspace(log10(0.3*mx), log10(3*mphi), 8);
G_FI    = ggxxrate(T_vec, mx, mphi, gx, g, interaction);
sv_xg   = arrayfun(@(T) svavg(mx,mphi,gx,g,T,interaction), T_vec);

% n_chi_eq for Dirac (g_spin=2), counting particle only
n_chi_eq = arrayfun(@(T) (2/(2*pi^2)) * T * mx^2 * besselk(2, mx/T), T_vec);

% G_FI from ggxxrate is gamma-bath collision integral. By detailed balance
% with identical-photon-pair statistics, equilibrium balance reads
%   (1/2) n_gamma_eq^2 <sv>_gg->xx  =  n_chi_eq * n_chibar_eq <sv>_xx->gg
% so the FI rate / 2 should match n_chi_eq^2 <sv>_xx->gg up to the
% normalization that ggxxrate factors into G(T).
R_rate = (G_FI/2) ./ (n_chi_eq.^2 .* sv_xg);

fprintf('\n=== (B) (G_FI/2) / (n_chi_eq^2 <sv>_xg) vs T ===\n');
fprintf('%10s %14s %14s %14s %14s\n','T [GeV]','G_FI','n_chi_eq^2*sv','ratio','log10(ratio)');
for k = 1:length(T_vec)
    fprintf('%10.3g %14.4e %14.4e %14.4e %14.4f\n', ...
        T_vec(k), G_FI(k), n_chi_eq(k)^2*sv_xg(k), R_rate(k), log10(R_rate(k)));
end

fprintf('\nIf BW fix worked, ratio (B) should be T-independent (constant).\n');
fprintf('If detailed balance factor of 2 is missing in code, ratio (A) shows it.\n');

% --- (C) all five families ---------------------------------------------------
% (A) covers heavyScalar only, which is how the missing identical-photon
% factor 2 in the PS/SP/PP production channels went unnoticed.  Here the
% parameters are chosen so xxggxsec returns the gamma gamma piece alone:
%   mx = 0.05 < mphi and below every meson pole -> phi phi / jj channels gated
%   s in (4mx^2, 4mpi0^2)                       -> every meson pair closed
% so no standalone re-derivation of the gg piece is needed.
mxC = 0.05; mphiC = 10; gxC = 1; gC = [1 1 1]/3;
fams = {'heavyScalar','heavyScalarPS','heavyScalarSP','heavyScalarPP','axialVector'};
sC = logspace(log10(4*mxC^2*1.05), log10(4*0.13498^2*0.95), 12);
R_expC = 2*(sC - 4*mxC^2)./sC;

fprintf('\n=== (C) sigma(gg->xx)/sigma(xx->gg) = 2(s-4mx^2)/s, all families ===\n');
fprintf('%16s %16s %8s\n','family','max |obs/exp - 1|','verdict');
allpass = true;
for k = 1:numel(fams)
    a = arrayfun(@(s) ggxxxsec(mxC,mphiC,gxC,gC,s,fams{k}), sC);
    b = arrayfun(@(s) xxggxsec(mxC,mphiC,gxC,gC,s,fams{k}), sC);
    dev = max(abs((a./b)./R_expC - 1));
    pass = dev < 1e-10;
    allpass = allpass && pass;
    tags = {'FAIL','PASS'};
    fprintf('%16s %16.3e %8s\n', fams{k}, dev, tags{pass+1});
end
verdicts = {'SOME FAILED','ALL PASS'};
fprintf('\n(C) %s\n', verdicts{allpass+1});
