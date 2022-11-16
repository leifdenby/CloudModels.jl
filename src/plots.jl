
function plot_profile(sols...; vars=[:r, :w, :T, :qv, :ql, :rh, :Δrho], wrap=3, size=400)
    plot_parts = []
    for var_name in vars
        push!(plot_parts, plot_profile_var(sols..., var_name=var_name))
    end
    
    ncols = minimum([wrap, length(plot_parts)])
    nrows = length(plot_parts) ÷ wrap
    if ncols * nrows < length(plot_parts)
        nrows += 1
    end

    plot(
        plot_parts...,
        layout=(nrows, ncols),
        size=(ncols*size, nrows*size),
        margin=20Plots.px,
    )
end


function plot_profile_var(sols...; var_name=:r)
    env_linestyle = :dash

    desc = Dict(
        :q_v => "water vapour concentation [g/kg]",
        :q_l => "cloud liquid concentation [g/kg]",
        :T => "temperature [g/kg]",
        :r => "cloud radius [m]",
        :w => "vertical velocity [m/s]"
    )
    
    p_var = plot()
    
    g(sol, v) = getindex.(sol.u, v)
    sg(sol, v, s) = g(sol, v) .* s

    colors = palette(:tab10)
    for (sol, color) in zip(sols, colors)
        pg(sol, v) = plot!(p_var, g(sol, v), g(sol, :z), xlabel=desc[v], label="", ylabel="alt [m]", color=color)
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

        if var_name == :r
            pg(sol, :r)
        elseif var_name == :w
            pg(sol, :w)
        elseif var_name == :T
            pg(sol, :T)
            plot!(p_var, ustrip.(T_env), ustrip.(z_prof), label=nothing, linestyle=env_linestyle, color=color)
        elseif var_name == :qv
            p = plot!(p_var, sg(sol, :q_v, 1.0e3), z_prof, label=nothing, xlabel="water vapour conc [g/kg]", color=color)
            plot!(p_var, qv_env .* 1.0e3, z_prof, label=nothing, color=color, linestyle=:dashdot)
            plot!(p_var, qv_sat_cld .* 1.0e3, z_prof, label=nothing, color=color, linestyle=:dashdot)
        elseif var_name == :rh
            plot!(p_var, rh_prof, g(sol, :z), label="", xlabel="relative humidity [1]", color=color)
            if length(p_var.series_list) == 1
                vline!(p_var, [1.0], linestyle=:dash, color=:black, label=nothing)
            end
        elseif var_name == :Δrho
            Δrho = ustrip.(ρ_cld - ρ_env) .* 1.0e3
            plot!(p_var, Δrho, z_prof, xlabel="Δρ [g/m^3]", color=color)
        elseif var_name == :ql
            plot!(p_var, sg(sol, :q_l, 1.0e3), z_prof, xlabel=desc[:q_l], label="", color=color)
        else
            throw("plotting of $(var_name) not implemented")
        end
    end
    
    return p_var
end