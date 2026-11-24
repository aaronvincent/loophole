# External limits

All curves: excluded above. Mass column is log10(m_φ / eV) unless stated.

## Mediator limits

| file | bound | coupling constrained | source |
|---|---|---|---|
| `SN1987AlogmphiloggN_pseudo.dat` | SN1987A | pseudoscalar nucleon coupling h̃_N | Ramani & Woolley, arXiv:1905.04319, Table 1 |
| `SN1987AScalarlogmphiloggN.dat` | SN1987A | scalar nucleon coupling h_N (no-decay limit) | Hardy, Sokolov & Stubbs, JHEP 04 (2025) 013, arXiv:2410.17347, Fig. 6 |
| `NeutronStarLoggNLogmphi.dat` | neutron-star cooling | scalar neutron coupling `hn` | Fiorillo et al., PRL 135 (2025) 211003, arXiv:2506.19906, Fig. 1 |
| `alphaphiNNvsmphieV.dat` | stellar cooling (TRGB + HB), key `stellar` | log10 α_φNN = g²/4π, g = (hp + hn)/2; S⊗S, P⊗S | Hardy & Lasenby, JHEP 02 (2017) 033, arXiv:1611.05852 |
| `axialConstr.dat` | fifth force, axial vector | log10 g_A² vs log10 λ (m), λ = 1/m_A′ | Cong et al., Rev. Mod. Phys. 97 (2025) 025005, arXiv:2408.15691, Fig. 10(b) |
| `PProtonProtonFifthForce.dat` | fifth force, spin–spin | log10 g_P², g_P = `htildeP` | Cong et al. (as above), Fig. 13 |
| `PProtonNeutronFifthForce.dat` | fifth force, spin–spin | log10 \|g_P g_N\|, g_N = `htildeN` | Cong et al. (as above), Fig. 13 |

## Direct-detection curves (`projections/`)

Columns: log10(m_χ / MeV), log10(σ / cm²). Used by `analysis/dd/` for the m_φ = 10 GeV overlay.

| file | content | source |
|---|---|---|
| `Combined Existing Spin Independent Limits.csv` | existing SI limits (envelope) | CRESST (arXiv:1707.06749, 1904.00498), LUX (1811.11241), CDEX (1905.00354, 2111.11243), XENON1T Migdal (1907.12771), EDELWEISS Migdal (2203.03993) |
| `Herald Spin Independent curve 1 100 m.csv` | HeRALD projection, 1 kg·day, 40 eV threshold | Hertel et al., PRD 100 (2019) 092007, arXiv:1810.06283, Fig. 9 |
| `Herald Spin Independent 1 MeV 100 m.csv` | HeRALD projection, 100 kg·yr, 1 meV threshold | as above |
| `QUESTDMC Spin Independent Projection.csv` | QUEST-DMC SI projection | Darvishi et al., JCAP 10 (2025) 044, arXiv:2505.17995 |
| `QUEST DMC Neutron Only Spin Dependent.csv` | QUEST-DMC SD-neutron projection | as above |
| `SuperCDMS Spin Independent Projection.csv` | SuperCDMS SNOLAB SI projection | Agnese et al., PRD 95 (2017) 082002, arXiv:1610.00006 |
| `Spin Dependent Limits.csv` | existing SD-proton limits | |
