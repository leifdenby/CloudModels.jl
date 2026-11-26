# Thesis Fig 3.3: entrainment vs drag. (C_D, β) varied together so that
# β + 3/8 C_D (the drag + entrainment term in the momentum equation) is
# constant: (0.7, 0.125), (0.5, 0.2), (0.3, 0.275); r0 ∈ {200, 300} m;
# l_pr = ∞ (no rain-out); z_INV raised to 4000 m.
# Read off the thesis figure: cloud base ≈ 700 m, w0 ≈ 1.6 m/s.
@isdefined(THESIS_GAMMA) || include(joinpath(@__DIR__, "common.jl"))

function fig3_3()
    env = thesis_profile(z_BL=600.0, z_INV=4000.0)
    runs = []
    for (C_D, β) in ((0.7, 0.125), (0.5, 0.2), (0.3, 0.275)), r0 in (200.0, 300.0)
        z0, F0 = cloudbase_state(env; r0=r0, w0=1.6)
        push!(runs, run_plume(env, F0, z0, (C_D=C_D, β=β, l_pr=Inf); z_max=4000.0,
                              label="C_D=$(C_D), β=$(β), r0=$(Int(r0))m"))
    end
    summarise("fig3_3", runs)
    plt = plot_profiles(
        runs, [:r, :w, :drho, :ql, :qr, :ss];
        title="Fig 3.3: entrainment vs drag (z_INV=4000m, l_pr=∞)",
        xlims=Dict(:r => (0, 1000)),
    )
    save_figure(plt, "fig3_3")
    return runs
end

fig3_3()
