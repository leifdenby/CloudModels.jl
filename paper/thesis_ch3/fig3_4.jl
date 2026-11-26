# Thesis Fig 3.4: effect of the entrainment rate. C_D = 0,
# β ∈ {0.2, 0.1, 0}, l_pr = 100 m, z_INV = 2400 m.
# Read off the thesis figure: cloud base ≈ 700 m, r0 ≈ 200 m, w0 ≈ 2 m/s.
#
# The caption says the cloud-top height is "significantly below the isotherm"
# only with entrainment and that "doubling the entrainment rate strongly
# affects the cloud-top height". To check that claim this script also writes
# `fig3_4_beta_sweep.png`: cloud top against β, with and without rain-out.
@isdefined(THESIS_GAMMA) || include(joinpath(@__DIR__, "common.jl"))

function fig3_4()
    env = thesis_profile(z_BL=600.0, z_INV=2400.0)
    z0, F0 = cloudbase_state(env; r0=200.0, w0=2.0)
    runs = [
        run_plume(env, F0, z0, (C_D=0.0, β=β, l_pr=100.0); z_max=4000.0,
                  label="C_D=0, β=$(β)")
        for β in (0.2, 0.1, 0.0)
    ]
    summarise("fig3_4", runs)
    plt = plot_profiles(
        runs, [:r, :w, :T, :qv, :ql, :qr, :mse, :drho, :ss];
        title="Fig 3.4: entrainment rate (C_D=0, l_pr=100m, z_INV=2400m)",
        ylims=(500.0, 4500.0), xlims=Dict(:r => (0, 600), :T => (276, 294)),
    )
    save_figure(plt, "fig3_4")

    # cloud top against β
    βs = collect(0.0:0.0125:0.4)
    sweep = Dict{Any,Any}()
    for l_pr in (100.0, Inf)
        sweep[l_pr] = [
            run_plume(env, F0, z0, (C_D=0.0, β=β, l_pr=l_pr); z_max=4000.0,
                      label="β=$(β), l_pr=$(l_pr)")
            for β in βs
        ]
    end
    summarise("fig3_4_beta_sweep", vcat(sweep[100.0], sweep[Inf]))
    plt2 = plot(; xlabel="entrainment coefficient β", ylabel="cloud top [m]",
                legend=:topright, size=(650, 450), ylims=(1000, 4200), xlims=(-0.01, 0.41),
                title="Cloud top vs β (Fig 3.4 set-up: C_D=0, r0=200m, w0=2m/s)",
                titlefontsize=10)
    hline!(plt2, [env.z_INV]; color=:gray50, ls=:dash, label="z_INV (isothermal above)")
    hline!(plt2, [4000.0]; color=:gray80, ls=:dot, label="integration top")
    plot!(plt2, βs, [r.z_top for r in sweep[100.0]]; marker=:circle, ms=3,
          color=PALETTE[1], label="CloudModels.jl, l_pr = 100 m (as Fig 3.4)")
    plot!(plt2, βs, [r.z_top for r in sweep[Inf]]; marker=:diamond, ms=3,
          color=PALETTE[2], label="CloudModels.jl, l_pr = ∞ (no rain-out)")
    # cloud tops read off thesis Fig 3.4 (β = 0 reaches the plotted top)
    scatter!(plt2, [0.0, 0.1, 0.2], [4000.0, 2250.0, 2430.0]; marker=:star5, ms=9,
             color=PALETTE[3], label="thesis Fig 3.4 (read off plot)")
    save_figure(plt2, "fig3_4_beta_sweep")
    return runs, sweep
end

fig3_4()
