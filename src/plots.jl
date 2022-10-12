function plot_profile(sol)
    desc = Dict(
        :q_v => "water vapour concentation [g/kg]",
        :q_l => "cloud liquid concentation [g/kg]",
        :T => "temperature [g/kg]",
        :r => "cloud radius [m]",
        :w => "vertical velocity [m/s]"
    )
    g(sol, v) = getindex.(sol.u, v)
    pg(sol, v) = plot(g(sol, v), g(sol, :z), xlabel=desc[v], label="")
    sg(sol, v, s) = g(sol, v) .* s

    params = sol.prob.p

    z_prof = g(sol, :z)
    p_env = uconvert.(u"Pa", params.environment.(z_prof * u"m", :p))
    T_env = uconvert.(u"K", params.environment.(z_prof * u"m", :T))
    qv_env = params.environment.(z_prof * u"m", :qv)
    T_cld = g(sol, :T) * u"K"
    qv_sat_env = CloudModels.calc_qv_sat.(T_env, p_env)
    qv_sat_cld = CloudModels.calc_qv_sat.(T_cld, p_env)
    rh_prof = g(sol, :q_v) ./ qv_sat_cld
    ρ_env = params.environment.(z_prof * u"m", :rho)

    qv_cld, ql_cld, qr_cld, qi_cld = g(sol, :q_v), g(sol, :q_l), g(sol, :q_r), g(sol, :q_i)
    qd_cld = 1.0 .- (qv_cld- ql_cld - qr_cld - qi_cld)
    ρ_cld = CloudModels.calc_mixture_density.(p_env, T_cld, qd_cld, qv_cld, ql_cld, qr_cld, qi_cld)

    p_temp = pg(sol, :T)
    plot!(p_temp, ustrip.(T_env), ustrip.(z_prof), label="env")

    plt_qv = plot(sg(sol, :q_v, 1.0e3), z_prof, label="qv", xlabel="water vapour conc [g/kg]")
    plot!(plt_qv, qv_sat_cld .* 1.0e3, z_prof, label="qv_sat", color="red")
    # plot!(plt_qv, qv_env .* 1.0e3, z_prof, label="qv_sat", color=:green, linestyle=:dash)

    plot_rh = plot(rh_prof, g(sol, :z), label="", xlabel="relative humidity [1]")
    vline!(plot_rh, [1.0], linestyle=:dash, color=:black)

    plot(
        pg(sol, :r),
        pg(sol, :w),
        p_temp,
        plt_qv,
        plot(sg(sol, :q_l, 1.0e3), z_prof, xlabel=desc[:q_l], label=""),
        plot_rh,
        plot(ustrip.(ρ_cld - ρ_env), z_prof, xlabel="Δρ [kg/m^3]"),
        layout=(4, 2),
        size=(600, 1000),
        margin=20Plots.px
    )
end