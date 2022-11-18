
function plot_profile(sols...; vars=[:r, :w, :T, :qv, :ql, :rh, :Δrho], wrap=3, size=400)
    N = length(vars)
    ncols = minimum([wrap, N])
    nrows = N ÷ wrap
    if ncols * nrows < N
        nrows += 1
    end

    plot_parts = []
    for (i, var_name) in enumerate(vars)
        if i == N
            add_legend = true
        else
            add_legend = false
        end
        push!(plot_parts, plot_profile_var(sols..., var_name=var_name, legend=add_legend))
    end

    plot(
        plot_parts...,
        layout=(nrows, ncols),
        size=(ncols*size, nrows*size),
        margin=20Plots.px,
    )
end

function make_name(sol)
    params = Dict(pairs(sol.prob.p))
    pop!(params, :environment)
    params_s = join(["$(k)=$(v)" for (k,v) in params], " ")
    return params_s
end

function plot_profile_var(sols...; var_name=:r, legend=:outertopright)
    env_linestyle = :dash

    desc = Dict(
        :q_v => "water vapour concentation [g/kg]",
        :q_l => "cloud liquid concentation [g/kg]",
        :T => "temperature [g/kg]",
        :r => "cloud radius [m]",
        :w => "vertical velocity [m/s]"
    )
    
    p_var = plot(;legend=legend)
    
    g(sol, v) = getindex.(sol.u, v)
    sg(sol, v, s) = g(sol, v) .* s

    colors = palette(:tab10)
    for (sol, color) in zip(sols, colors)
        label = make_name(sol)
        pg(sol, v) = plot!(p_var, g(sol, v), g(sol, :z), xlabel=desc[v], label=label, ylabel="alt [m]", color=color)
        params = sol.prob.p

        z_prof = g(sol, :z)
        p_env = uconvert.(u"Pa", params.environment.(z_prof * u"m", :p))
        T_env = uconvert.(u"K", params.environment.(z_prof * u"m", :T))
        qv_env = params.environment.(z_prof * u"m", :qv)
        ρ_env = params.environment.(z_prof * u"m", :rho)
        rh_env = params.environment.(z_prof * u"m", :rh)

        T_cld = g(sol, :T) * u"K"
        qv_sat_cld = CloudModels.calc_qv_sat.(T_cld, p_env)
        rh_c = g(sol, :q_v) ./ qv_sat_cld

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
            p = plot!(p_var, sg(sol, :q_v, 1.0e3), z_prof, label=label, xlabel="water vapour conc [g/kg]", color=color)
            plot!(p_var, qv_env .* 1.0e3, z_prof, label=label, color=color, linestyle=env_linestyle)
            plot!(p_var, qv_sat_cld .* 1.0e3, z_prof, label=label, color=color, linestyle=:dashdot)
        elseif var_name == :rh
            plot!(p_var, rh_c, g(sol, :z), label=label, xlabel="relative humidity [1]", color=color)
            if length(p_var.series_list) == 1
                vline!(p_var, [1.0], linestyle=:dash, color=:black, label=label)
            end
            plot!(p_var, rh_env, g(sol, :z), label="", color=color, linestyle=env_linestyle)
        elseif var_name == :Δrho
            Δrho = ustrip.(ρ_cld - ρ_env) .* 1.0e3
            plot!(p_var, Δrho, z_prof, xlabel="Δρ [g/m^3]", color=color, label=label)
        elseif var_name == :ql
            plot!(p_var, sg(sol, :q_l, 1.0e3), z_prof, xlabel=desc[:q_l], label=label, color=color)
        else
            throw("plotting of $(var_name) not implemented")
        end
    end
    
    return p_var
end