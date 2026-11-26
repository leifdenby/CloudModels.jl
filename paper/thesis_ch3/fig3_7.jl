# Thesis Fig 3.7: rain-out length-scale. l_pr ∈ {10, 100, 1000, ∞} m;
# C_D = 0.5, β = 0.2; z_BL = 600 m, z_INV = 2400 m.
# Read off the thesis figure: cloud base ≈ 700 m, r0 ≈ 200 m, w0 ≈ 2 m/s.
@isdefined(THESIS_GAMMA) || include(joinpath(@__DIR__, "common.jl"))

function fig3_7()
    env = thesis_profile(z_BL=600.0, z_INV=2400.0)
    z0, F0 = cloudbase_state(env; r0=200.0, w0=2.0)
    runs = [
        run_plume(env, F0, z0, (C_D=0.5, β=0.2, l_pr=l_pr); z_max=4000.0,
                  label="l_pr=$(l_pr)m")
        for l_pr in (10.0, 100.0, 1000.0, Inf)
    ]
    summarise("fig3_7", runs)
    plt = plot_profiles(
        runs, [:r, :w, :T, :qv, :ql, :qr, :drho, :mse, :ss];
        title="Fig 3.7: rain-out length-scale (C_D=0.5, β=0.2, z_INV=2400m)",
        xlims=Dict(:r => (0, 1400)),
    )
    save_figure(plt, "fig3_7")
    return runs
end

fig3_7()
