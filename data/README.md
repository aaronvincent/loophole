# External limits

All curves: excluded above. Mass column is log10(m_φ / eV) unless stated.

| file | bound | coupling constrained | source |
|---|---|---|---|
| `SN1987A` | SN1987A | pseudoscalar | Hari |
| `SN1987AScalar` | SN1987A | scalar | Hardy, Sokolov, Stubbs, arXiv:2410.17347 |
| `NeutronStar` | neutron star | scalar; code uses neutron coupling `hn` | Damiano, arXiv:2506.19906 |
| `alphaphiNNvsmphieV` | stellar cooling (TRGB + HB), key `stellar` | log10 α_φNN = g²/4π, g = (hp + hn)/2; S⊗S, P⊗S only | Hardy & Lasenby (`HardyLasenbyalphaphiNN.png`) |
| `Axial` | axial-vector bounds | log10 g_A² vs log10 λ (m) | arXiv:2408.15691 |
| `PProtonProtonFifthForce` | spin–spin fifth force | log10 g_P², g_P = `htildeP` | `ppfifthforce.png` |
| `PProtonNeutronFifthForce` | spin–spin fifth force | log10 \|g_P g_N\|, g_N = `htildeN` | `ppfifthforce.png` |

