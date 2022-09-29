using Test
using CloudModels
using Unitful
using OrdinaryDiffEq
using ComponentArrays
using Plots
using Revise

@testset "eqns" begin
    F = ComponentArray(r=100u"m", w=1.0u"m/s", T=300u"K", q_v=13.0e-3, 
    q_l=0.0, q_r=0.0, q_i=0.0, p=100e3u"Pa", q_pr=0.0)
    dFdz = ComponentArray(Dict([ v => 0.0 * unit(F[v] * u"1/m") for v in keys(F)]))
    env_profile = CloudModels.StandardIsentropicAtmosphere()
    p = (environment=env_profile, β=0.2)
    CloudModels.parcel_equations!(dFdz, F, 0.0u"m", p)
end

#@testset "integrate equations" begin

#F = ComponentArray(r=100u"m", w=1.0u"m/s", T=300u"K", q_v=13.0e-3, q_l=0.0, q_r=0.0, q_i=0.0, p=100e3u"Pa", q_pr=0.0)


function setup_callbacks()
    condition(u, t, integrator) = u.z - 2000.0
    cb_ztop = ContinuousCallback(condition,terminate!)

    mphys_zero_cbs = [
        DiscreteCallback((u, t, i) -> u[v] < 0.0, terminate!) for v in [:q_v, :q_r, :q_l]
    ]
    CallbackSet(cb_ztop, mphys_zero_cbs...)
end

env_profile = CloudModels.ProfileRICO.RICO_profile()

function make_initial_condition(env_profile, z0)
    dT = 0.0
    dqv = 3.0e-3
    T0 = ustrip(env_profile(z0 * u"m", :T)) + dT
    qv0 = ustrip(env_profile(z0 * u"m", :qv)) + dqv
    F = ComponentArray(r=300, w=1.0, T=T0, q_v=qv0, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0, z=z0)
    return F
end

function nounit_parcel_equations!(dFdt, F, params, t)
    F_units = ComponentArray(
        r=F.r * u"m",
        w=F.w * u"m/s",
        T=F.T * u"K",
        q_v=F.q_v,
        q_l=F.q_l,
        q_r=F.q_r,
        q_i=F.q_i,
        q_pr=F.q_pr,
        p=F.p * u"Pa",
        z=F.z * u"m",
    )
    dFdz_units = ComponentArray(Dict([ v => 0.0 * unit(F_units[v] * u"1/m") for v in keys(F)]))

    CloudModels.parcel_equations!(dFdz_units, F_units, F.z * u"m", params)
    dFdz = ustrip.(dFdz_units)
    dzdt = F.w
    for v in keys(F)
        dFdt[v] = dFdz[v] / dzdt
    end
    dFdt[:z] = dzdt
end

#env_profile = CloudModels.StandardIsentropicAtmosphere()
#env_profile = CloudModels.StableAtmosphereTest()
params = (environment=env_profile, β=0.2)
F = make_initial_condition(env_profile, 300.0)

prob = ODEProblem(nounit_parcel_equations!, F, [0.0, 1000.0], params)
sol = solve(prob, Tsit5(), saveat=0.1, callback=setup_callbacks())


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

    z_prof = g(sol, :z)
    p_env = params.environment.(z_prof * u"m", :p) .|> u"Pa"
    T_env = params.environment.(z_prof * u"m", :T) .|> u"K"
    qv_env = params.environment.(z_prof * u"m", :qv)
    T_cld = g(sol, :T) * u"K"
    qv_sat_env = CloudModels.calc_qv_sat.(T_env, p_env)
    qv_sat_cld = CloudModels.calc_qv_sat.(T_cld, p_env)
    rh_prof = g(sol, :q_v) ./ qv_sat_cld

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
        layout=(3, 2),
        size=(600, 1000),
        margin=20Plots.px
    )
end

plot_profile(sol)