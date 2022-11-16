using Test
using CloudModels
using Unitful
using OrdinaryDiffEq
using ComponentArrays


include("integration_common.jl")


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
        CloudModels.parcel_equations!(dFdz, F, 0.0u"m", p)
        @test dFdz != dFdz_initial
    end
end


getvar(sol, v) = getindex.(sol.u, v)

@testset "parcel-integration" begin
    env_profile = CloudModels.ProfileRICO.RICO_profile()
    params = (environment=env_profile, β=0.2)
    
    @testset "changing initial velocity" begin
        sols = []
        for w0 in [0.5, 1.0, 2.0]
            U0 = make_initial_condition(env_profile, r0=400, w0=w0, dqv=1.0e-3, z0=500.0)
            prob = ODEProblem(parcel_equations!, U0, [0.0, 700.0], params)
            sol = solve(prob, Euler(), dt=1.0, saveat=10.0, callback=setup_callbacks(z_max=00.0))
            @test sum(getvar(sol, :q_l)) > 0.0
            push!(sols, sol)
        end
        
        # higher initial velocity should lead to higher cloud-top height
        zmax = [maximum(getvar(sol, :z)) for sol in sols]
        @test sort(zmax) == zmax
    end
end