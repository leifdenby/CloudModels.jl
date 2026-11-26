# Thesis Ch 3 sensitivity figures (Figs 3.2–3.8)

Scripts that redo the cloud-model sensitivity figures of thesis Ch 3
(L. Denby 2017, PhD thesis, Cambridge,
<https://www.repository.cam.ac.uk/handle/1810/269850>, Sect. 3.4, printed
pp. 46–59) with CloudModels.jl. They are not part of the package.

Run all of them from the repository root:

```
julia --project=. paper/thesis_ch3/make_all.jl [outdir] [dRHdz_1]
```

or a single one, e.g. `julia --project=. paper/thesis_ch3/fig3_4.jl`. Each
script writes `<name>.png` and `<name>_summary.csv` (cloud base, cloud top,
peak w and why the integration stopped, per run) to `outdir`
(default `paper/thesis_ch3/output/`, which git ignores).

| script | thesis figure | varied |
|---|---|---|
| `fig3_2.jl` | 3.2, drag without entrainment | C_D ∈ {1, 0.5, 0.25, 0} |
| `fig3_3.jl` | 3.3, entrainment vs drag | (C_D, β) with β + 3/8 C_D fixed; r0 ∈ {200, 300} m |
| `fig3_4.jl` | 3.4, entrainment rate | β ∈ {0.2, 0.1, 0}, plus cloud top against β (`fig3_4_beta_sweep.png`) |
| `fig3_5.jl` | 3.5, no entrainment of MSE or hydrometeors | β ∈ {0.1, 0.2, 0.3, 0.4, 10} |
| `fig3_6.jl` | 3.6, ambient RH gradient | dRH/dz ∈ {0, -10, -20} %/km |
| `fig3_7.jl` | 3.7, rain-out length-scale | l_pr ∈ {10, 100, 1000, ∞} m |
| `fig3_8.jl` | 3.8, boundary-layer moisture perturbation | Δq_v ∈ {0, 0.2, 0.4} g/kg; r0 ∈ {100, 200} m |

## Set-up

Shared by all scripts (`common.jl`):

- **Environment:** `TwoLayerMoistIsentropicPBL` with RH0 = 70 %,
  T0 = 299.3 K, p0 = 101325 Pa (thesis Sect. 5.1.2), and z_BL and z_INV
  as in each figure's title.
- **RH gradient above the boundary layer:** −0.5e-3 1/m (−50 %/km), except
  in Fig 3.6. The thesis doesn't state this value. It is inferred from the
  environment lines in the thesis figures: q_v,e ≈ 13.5 g/kg at 700 m and 0
  near 2.5 km, and s_e ≈ 310 kJ/kg at z_INV. The package default,
  −0.1e-3 1/m, gives a much moister environment that doesn't match them.
  Override it with the second argument.
- **Cloud base:** the lifting condensation level of the unperturbed surface
  parcel (`compute_LCL`), at 742 m. The thesis used CCFM's cloud trigger,
  which isn't ported, and its bases are 700–880 m.
- **Parameters:** γ = 0.5 (thesis value, passed explicitly; the package
  default is 1.0). Each script sets C_D, β and l_pr as given in the thesis
  legends.
- **Read off the thesis plots, not stated in the text:** r0 and w0 (noted
  in each script), and C_D = 0 for Fig 3.5.
- **Integration:** Tsit5 with rtol 1e-8, atol 1e-10 and steps of at most
  5 m (thesis Sect. 3.3). It stops when w crosses 0 or r exceeds 10 km
  (pyclouds' default stopping conditions), or at the top of the domain. The
  "cloud top" is the height where the integration stopped.
- **Tendencies:** `plume_equations!` is called through a wrapper that zeroes
  dF/dz first, so the scripts don't depend on the `dFdz .= 0.0` fix in PR
  #19.

Diagnostics follow pyclouds' plotting code (`pyclouds/plot/parcel.py`). The
moist static energy is c̄_l T + q_v L_v + g z with constant L_v. The density
difference is ρ_c − ρ_e, so negative means buoyant. Supersaturation is
(q_v/q_v,sat − 1)·100 %.
