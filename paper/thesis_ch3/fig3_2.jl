# Thesis Fig 3.2: effect of drag without entrainment.
# β = 0, l_pr = 100 m, C_D ∈ {1, 0.5, 0.25, 0}; profile RH0 = 70 %,
# T0 = 299.3 K, z_BL = 600 m, z_INV = 2400 m; plume started at the LCL.
# Read off the thesis figure: cloud base ≈ 730 m, r0 ≈ 200 m, w0 ≈ 2 m/s.
@isdefined(THESIS_GAMMA) || include(joinpath(@__DIR__, "common.jl"))

function fig3_2()
    env = thesis_profile(z_BL=600.0, z_INV=2400.0)
    z0, F0 = cloudbase_state(env; r0=200.0, w0=2.0)
    runs = [
        run_plume(env, F0, z0, (C_D=C_D, β=0.0, l_pr=100.0); z_max=4000.0,
                  label="C_D=$(C_D), β=0")
        for C_D in (1.0, 0.5, 0.25, 0.0)
    ]
    summarise("fig3_2", runs)
    plt = plot_profiles(
        runs, [:r, :w, :T, :qv, :ql, :qr, :mse, :drho, :ss];
        title="Fig 3.2: drag, no entrainment (RH0=70%, T0=299.3K, z_BL=600m, z_INV=2400m, l_pr=100m)",
        ylims=(500.0, 4000.0), xlims=Dict(:r => (0, 600), :T => (276, 294)),
    )
    save_figure(plt, "fig3_2")
    return runs
end

fig3_2()
