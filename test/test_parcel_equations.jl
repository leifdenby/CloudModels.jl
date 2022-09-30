using Test
using CloudModels
using Unitful
using OrdinaryDiffEq
using ComponentArrays


@testset "eqns" begin
    F = ComponentArray(r=100u"m", w=1.0u"m/s", T=300u"K", q_v=13.0e-3, 
    q_l=0.0, q_r=0.0, q_i=0.0, p=100e3u"Pa", q_pr=0.0)
    dFdz = ComponentArray(Dict([ v => 0.0 * unit(F[v] * u"1/m") for v in keys(F)]))
    env_profile = CloudModels.StandardIsentropicAtmosphere()
    p = (environment=env_profile, β=0.2)
    CloudModels.parcel_equations!(dFdz, F, 0.0u"m", p)
end


function setup_callbacks()
    condition(u, t, integrator) = u.z - 2000.0
    cb_ztop = ContinuousCallback(condition,terminate!)

    mphys_zero_cbs = [
        DiscreteCallback((u, t, i) -> u[v] < 0.0, x-> terminate!(x, :NegativeMphys)) for v in [:q_v, :q_r, :q_l]
    ]
    CallbackSet(cb_ztop, mphys_zero_cbs...)
end


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

@testset "parcel-integration" begin
    env_profile = CloudModels.ProfileRICO.RICO_profile()
    params = (environment=env_profile, β=0.2)
    F = make_initial_condition(env_profile, 300.0)

    prob = ODEProblem(nounit_parcel_equations!, F, [0.0, 1000.0], params)
    sol = solve(prob, Euler(), saveat=0.1, callback=setup_callbacks(), dt=1.0)
    CloudModels.plot_profile(sol)

    g(sol, v) = getindex.(sol.u, v)
    #check that some condensation has occoured :)
    @test sum(g(sol, :q_l)) > 0.0
end