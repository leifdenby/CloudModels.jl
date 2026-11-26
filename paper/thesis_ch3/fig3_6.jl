# Thesis Fig 3.6: ambient relative humidity above the boundary layer.
# dRH/dz ∈ {0, -10, -20} %/km; z_BL = 700 m, z_INV = 6000 m; C_D = 0.5,
# β = 0.2, l_pr = ∞.
# Read off the thesis figure: cloud base ≈ 700 m, r0 ≈ 100 m, w0 ≈ 1 m/s.
@isdefined(THESIS_GAMMA) || include(joinpath(@__DIR__, "common.jl"))

function fig3_6()
    runs = []
    envs = []
    for dRHdz in (0.0, -0.1e-3, -0.2e-3)
        env = thesis_profile(z_BL=700.0, z_INV=6000.0, dRHdz_1=dRHdz)
        z0, F0 = cloudbase_state(env; r0=100.0, w0=1.0)
        push!(runs, run_plume(env, F0, z0, (C_D=0.5, β=0.2, l_pr=Inf); z_max=6000.0,
                              label="dRH/dz=$(dRHdz * 1e5)%/km"))
        push!(envs, env)
    end
    summarise("fig3_6", runs)
    plt = plot_profiles(
        runs, [:r, :w, :T, :drho, :mse, :ql, :dqv, :ss, :rh];
        title="Fig 3.6: ambient RH gradient (z_BL=700m, z_INV=6000m, C_D=0.5, β=0.2, l_pr=∞)",
        ylims=(500.0, 6200.0), xlims=Dict(:r => (0, 1000), :rh => (65, 105)),
        envs=envs,
    )
    save_figure(plt, "fig3_6")
    return runs
end

fig3_6()
