# Thesis Fig 3.8: boundary-layer moisture perturbation of the surface parcel,
# Δq_v ∈ {0, 0.2, 0.4} g/kg, at cloud-base radius r0 ∈ {100, 200} m (thesis
# Sect. 3.4.5); C_D = 0.506, β = 0.2, l_pr = 100 m; z_BL = 600 m,
# z_INV = 2400 m. Cloud base is the LCL of the perturbed surface parcel.
# Read off the thesis figure: w0 ≈ 2 m/s.
@isdefined(THESIS_GAMMA) || include(joinpath(@__DIR__, "common.jl"))

function fig3_8()
    env = thesis_profile(z_BL=600.0, z_INV=2400.0)
    runs = []
    for dqv in (0.0, 0.2e-3, 0.4e-3), r0 in (100.0, 200.0)
        z0, F0 = cloudbase_state(env; r0=r0, w0=2.0, dqv=dqv)
        push!(runs, run_plume(env, F0, z0, (C_D=0.506, β=0.2, l_pr=100.0); z_max=4000.0,
                              label="Δq_v=$(dqv * 1e3)g/kg, r0=$(Int(r0))m"))
    end
    summarise("fig3_8", runs)
    plt = plot_profiles(
        runs, [:r, :w, :T, :qv, :ql, :qr, :drho, :mse, :ss];
        title="Fig 3.8: moisture perturbation (C_D=0.506, β=0.2, l_pr=100m)",
        xlims=Dict(:r => (0, 1000)),
    )
    save_figure(plt, "fig3_8")
    return runs
end

fig3_8()
