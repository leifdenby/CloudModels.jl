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

F = ComponentArray(r=300, w=1.0, T=300, q_v=10.0e-3, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0, z=0.0)

qd0 = 1.0 - F.q_v - F.q_l - F.q_r - F.q_i
p0 = params.environment(F.z * u"m", :p)
rho0 = CloudModels.calc_mixture_density(p0, F.T * u"K", qd0, F.q_v, F.q_l, F.q_r, F.q_i)
rho0_env = params.environment(F.z * u"m", :rho)


#env_profile = CloudModels.StandardIsentropicAtmosphere()
env_profile = CloudModels.StandardIsothermalAtmosphere()
params = (environment=env_profile, β=0.0)
prob = ODEProblem(nounit_parcel_equations!, F, [0.0, 1000], params)
sol = solve(prob, Tsit5())

g(sol, v) = getindex.(sol.u, v)
pg(sol, v) = plot(g(sol, v), g(sol, :z), label=string(v))

g(sol, :q_v)

p_env = params.environment.(g(sol, :z) * u"m", :p) .|> u"Pa"
qv_sat_env = CloudModels.calc_qv_sat.(g(sol, :T) * u"K", p_env)
rh_prof = g(sol, :q_v) ./ qv_sat_env
z_prof = g(sol, :z)

plot(ustrip.(p_env), z_prof)

function plot_qv(sol)
    p = pg(sol, :q_v)
    plot!(p, qv_sat_env, g(sol, :z), label="qv_sat")
end

plot_qv(sol)

plot(
    pg(sol, :r),
    pg(sol, :w),
    pg(sol, :T),
    pg(sol, :q_v),
    pg(sol, :q_l),
    plot(rh_prof, g(sol, :z), label="rh")
)

plot(
    plot(sol, vars=[:r]),
    plot(sol, vars=[:w])
)
