using Test
using CloudModels
# using Unitful
using OrdinaryDiffEq
using ComponentArrays


function setup_callbacks()
    cb_zbottom = ContinuousCallback((u, z, integrator) -> z - 200.0,terminate!)
    radius_negative = ContinuousCallback((u, z, integrator) -> u.r, terminate!)
    # TODO: improve this cloud-top prediction to take into account the radius
    neutral_buoyancy_cloudtop = DiscreteCallback((u, z, integrator) -> u.w < 0 || u.r > 3000, terminate!)

    mphys_zero_cbs = [
        DiscreteCallback((u, z, i) -> u[v] < 0.0, x-> terminate!(x, :NegativeMphys)) for v in [:q_v, :q_r, :q_l]
    ]
    CallbackSet(cb_zbottom, radius_negative, mphys_zero_cbs...)
end


function make_initial_condition(env_profile; z0=300, dT=0.0, dqv=1.0e-3, r0=300, w0=1.0)
    T0 = ustrip(env_profile(z0 * u"m", :T)) + dT
    qv0 = ustrip(env_profile(z0 * u"m", :qv)) + dqv
    F = ComponentArray(r=r0, w=w0, T=T0, q_v=qv0, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0)
    return F
end


@testset "eqns" begin
    ICs = Dict(
        "simple" => ComponentArray(r=100u"m", w=1.0u"m/s", T=300u"K", q_v=13.0e-3, q_l=0.0, q_r=0.0, q_i=0.0, p=100e3u"Pa", q_pr=0.0),
        "rain cond/evap" => ComponentArray(r = 26.1u"m", w = 3.0u"m/s", T = 326.8u"K", q_v = 0.64, q_l = 0.019, q_r = 1.16e-7, q_i = 0.0, q_pr = 0.0, p = 98077.78139210862u"Pa", z = 303.8320580602206u"m")
    )
    @testset "$(name) state" for (name, F) in ICs
        dFdz = ComponentArray(Dict([ v => 0.0 * unit(F[v] * u"1/m") for v in keys(F)]))
        env_profile = CloudModels.StandardIsentropicAtmosphere()
        p = (environment=env_profile, β=0.2)
        dFdz_initial = copy(dFdz)
        CloudModels.plume_equations!(dFdz, F, 0.0u"m", p)
        @test dFdz != dFdz_initial
    end
end


getvar(sol, v) = getindex.(sol.u, v)

@testset "plume-integration" begin
    env_profile = CloudModels.ProfileRICO.RICO_profile()
    params = (environment=env_profile, β=0.2)
    
    @testset "changing initial velocity" begin
        sols = []
        for w0 in [0.5, 1.0, 2.0]
            z0 = 500
            U0 = make_initial_condition(env_profile, r0=400, w0=w0, dqv=1.0e-3, z0=z0)
            prob = ODEProblem(CloudModels.plume_equations!, U0, [z0, 1000.0], params)
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
            z0 = 500
            U0 = make_initial_condition(env_profile, r0=400, w0=1.0, dqv=1.0e-3, z0=z0)
            prob = ODEProblem(CloudModels.plume_equations!, U0, [z0, 1000.0], params)
            sol = solve(prob, Euler(), dt=1.0, saveat=10.0, callback=setup_callbacks(z_max=00.0))
            push!(sols, sol)
        end
        
        # lower entrainment rate should lead to higher cloud-top height
        zmax = [maximum(getvar(sol, :z)) for sol in sols]
        @test sort(zmax) == zmax
    end
end