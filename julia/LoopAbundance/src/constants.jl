# Physical constants — single source of truth (port of alltheconstants.m plus
# the constants that the MATLAB physics files (ggxxxsec.m, xxggxsec.m,
# phi_partial_widths.m, ...) each redefine locally).  All energies in GeV.

# --- cosmology / EM (alltheconstants.m) ---------------------------------
const HBARC   = 1.97e-14      # GeV·cm
const GEV_TO_G = 1.783e-24    # grams per GeV (1 GeV/c² in g), for σ/m in cm²/g
const KB      = 8.62e-14      # GeV/K
const HUBBLE_h = 0.7
const ME      = 511e-6        # electron mass, GeV
const RHOC    = 1.05e-5 * HUBBLE_h^2 * HBARC^3   # critical density, GeV^4
const MPL     = 1.22e19       # standard (non-reduced) Planck mass, GeV
const ALPHA_EM = 1/137
const T0      = 2.7 * KB      # CMB temperature today, GeV

# entropy density today; depends on g_s(T0) so computed from gstarofT
const S0_ENTROPY = gstarofT(T0, 1) * 2 * pi^2 / 45 * T0^3   # GeV^3

# relic-abundance conversion constants (freeze_loop.m:135-136, hard-coded there)
const S0_CM3      = 2891.2      # CMB entropy density today, cm^-3
const RHO_CRIT_H2 = 1.054e-5    # GeV·cm^-3

# --- meson / hadron masses ---------------------------------------------
const MPI0    = 0.13498
const MPI_PM  = 0.13957
const MK      = 0.493677
const MK0     = 0.497611
const M_ETA   = 0.547862
const M_ETAP  = 0.95778
const MN      = 0.938          # nucleon mass

# --- meson widths (used as BW regulators) ------------------------------
const GPI     = 7.7e-9
const G_ETA   = 1.31e-6
const G_ETAP  = 1.88e-4        # eta' total width
const G_ETAP_GG = 4.34e-6      # eta' -> gamma gamma partial width
const G_ETA_GG  = 5.16e-7      # eta  -> gamma gamma partial width (1.31 keV * 39.36%)

# --- chiral / quark parameters -----------------------------------------
const B0      = 0.13957^2 / (0.0022 + 0.00469)  # chiral condensate, GeV: mpi+^2/(mu+md) ~ 2.827 (draft convention; was 2.6)
const FPI     = 0.0924         # pion decay constant, GeV
const L4      = 1.87e-3        # NLO chiral LECs
const L5      = 1.22e-3
const L6      = 1.46e-3
const L7      = -0.39e-3
const M_UP    = 0.0022
const M_DN    = 0.00469
const M_S     = 0.095
const FPU     = 0.02           # nucleon scalar couplings
const FPD     = 0.026
const FPS     = 0.118

"""
    dscript(g) -> D

Chiral φγγ loop factor `D = ½ (B0(g_u+g_d)/mπ+² + B0(g_u+g_s)/m_K+²)²`.
"""
dscript(g) = 0.5*(B0 * (g[1] + g[2]) / MPI_PM^2 + B0 * (g[1] + g[3]) / MK^2)^2  # charged-meson loop (pi+, K+), low-energy limit A0 -> 1/3
