using Test
using CloudModels
# using Unitful
using OrdinaryDiffEq
using ComponentArrays


@testset "eqns" begin
    ICs = Dict(
        "simple" => ComponentArray(r=100u"m", w=1.0u"m/s", T=300u"K", q_v=13.0e-3, q_l=0.0, q_r=0.0, q_i=0.0, p=100e3u"Pa", q_pr=0.0, z=0.0u"m"),
        "rain cond/evap" => ComponentArray(r = 26.1u"m", w = 3.0u"m/s", T = 326.8u"K", q_v = 0.64, q_l = 0.019, q_r = 1.16e-7, q_i = 0.0, q_pr = 0.0, p = 98077.78139210862u"Pa", z = 303.8320580602206u"m")
    )
    @testset "$(name) state" for (name, F) in ICs
        dFdt = ComponentArray(Dict([ v => 0.0 * unit(F[v] * u"1/s") for v in keys(F)]))
        env_profile = CloudModels.StandardIsentropicAtmosphere()
        p = (environment=env_profile, β=0.2)
        dFdt_initial = copy(dFdt)
        CloudModels.parcel_equations!(dFdt, F, p, 0.0u"s")
        @test dFdt != dFdt_initial
    end
end


function setup_callbacks(;z_max=2000.0)
    condition(u, t, integrator) = u.z - z_max
    cb_ztop = ContinuousCallback(condition,terminate!)
    cb_zbottom = ContinuousCallback((u, t, integrator) -> u.z - 200.0,terminate!)
    radius_negative = ContinuousCallback((u, t, integrator) -> u.r, terminate!)
    # TODO: improve this cloud-top prediction to take into account the radius
    neutral_buoyancy_cloudtop = DiscreteCallback((u, t, integrator) -> u.w < 0 || u.r > 3000, terminate!)

    mphys_zero_cbs = [
        DiscreteCallback((u, t, i) -> u[v] < 0.0, x-> terminate!(x, :NegativeMphys)) for v in [:q_v, :q_r, :q_l]
    ]
    CallbackSet(cb_ztop, cb_zbottom, radius_negative, mphys_zero_cbs...)
end


function make_initial_condition(env_profile; z0=300, dT=0.0, dqv=1.0e-3, r0=300, w0=1.0)
    T0 = ustrip(env_profile(z0 * u"m", :T)) + dT
    qv0 = ustrip(env_profile(z0 * u"m", :qv)) + dqv
    F = ComponentArray(r=r0, w=w0, T=T0, q_v=qv0, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0, z=z0)
    return F
end


getvar(sol, v) = getindex.(sol.u, v)

@testset "parcel-integration" begin
    env_profile = CloudModels.ProfileRICO.RICO_profile()
    params = (environment=env_profile, β=0.2)
    
    @testset "changing initial velocity" begin
        sols = []
        w0 = 1.0
        for w0 in [0.5, 1.0, 2.0]
            U0 = make_initial_condition(env_profile, r0=400, w0=w0, dqv=1.0e-3, z0=500.0)
            prob = ODEProblem(CloudModels.parcel_equations!, U0, [0.0, 700.0], params)
            sol = solve(prob, Euler(), dt=1.0, saveat=10.0, callback=setup_callbacks())
            @test sum(getvar(sol, :q_l)) > 0.0
            push!(sols, sol)
        end
        
        # higher initial velocity should lead to higher cloud-top height
        zmax = [maximum(getvar(sol, :z)) for sol in sols]
        @test sort(zmax) == zmax
    end

    @testset "changing entrainment rate" begin
        sols = []
        for β in [0.2, 0.1, 0.0]
            params = (environment=env_profile, β=β)
            U0 = make_initial_condition(env_profile, r0=400, w0=1.0, dqv=1.0e-3, z0=500.0)
            prob = ODEProblem(CloudModels.parcel_equations!, U0, [0.0, 700.0], params)
            sol = solve(prob, Euler(), dt=1.0, saveat=10.0, callback=setup_callbacks(z_max=00.0))
            push!(sols, sol)
        end
        
        # lower entrainment rate should lead to higher cloud-top height
        zmax = [maximum(getvar(sol, :z)) for sol in sols]
        @test sort(zmax) == zmax
    end
end