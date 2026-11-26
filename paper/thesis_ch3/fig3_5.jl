# Thesis Fig 3.5: entrainment of mass and momentum only. Moist static energy
# and hydrometeors are not entrained (`entrain_moist_static_energy=false`,
# `entrain_hydrometeors=false`: dT/dz = -g/c_p,d plus latent heating from the
# microphysics, as in pyclouds), β ∈ {0.1, 0.2, 0.3, 0.4}; the text also
# mentions β = 10, which is added.
# Not given in the thesis: C_D and l_pr. In the thesis figure the peak w
# falls as β rises (7.2 m/s at β = 0.1 to 6.3 m/s at β = 0.4). With
# C_D = 0.5 it rises instead (faster-growing radius, less drag), so C_D = 0 is
# used; l_pr = 100 m (the thesis default). Read off the figure: r0 ≈ 100 m,
# w0 ≈ 1 m/s (cloud base there is ≈ 830 m, higher than the LCL used here).
@isdefined(THESIS_GAMMA) || include(joinpath(@__DIR__, "common.jl"))

function fig3_5()
    env = thesis_profile(z_BL=600.0, z_INV=2400.0)
    z0, F0 = cloudbase_state(env; r0=100.0, w0=1.0)
    switches = (entrain_moist_static_energy=false, entrain_hydrometeors=false)
    runs = [
        run_plume(env, F0, z0, (C_D=0.0, β=β, l_pr=100.0, switches...); z_max=4000.0,
                  label="β=$(β)")
        for β in (0.1, 0.2, 0.3, 0.4, 10.0)
    ]
    summarise("fig3_5", runs)
    plt = plot_profiles(
        runs, [:r, :w, :T, :drho, :mse, :ss];
        title="Fig 3.5: no entrainment of MSE or hydrometeors (C_D=0, l_pr=100m)",
        ylims=(500.0, 4000.0), xlims=Dict(:r => (0, 2000), :T => (278, 294)),
    )
    save_figure(plt, "fig3_5")
    return runs
end

fig3_5()
