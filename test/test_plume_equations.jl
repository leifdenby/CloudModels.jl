using Test
using CloudModels
# using Unitful
using OrdinaryDiffEq
using ComponentArrays

# NB: in the plume formulation height `z` is the independent variable, so
# (unlike the parcel equations) `z` is not part of the state vector

function setup_plume_callbacks()
    # don't allow radius to become negative
    radius_negative = ContinuousCallback((u, z, integrator) -> u.r, terminate!)
    # cloud-top: stop when the vertical velocity reaches zero (dw/dz ~ 1/w
    # diverges there)
    w_zero = ContinuousCallback((u, z, integrator) -> u.w, terminate!)

    mphys_zero_cbs = [
        DiscreteCallback((u, z, i) -> u[v] < 0.0, x-> terminate!(x, :NegativeMphys)) for v in [:q_v, :q_r, :q_l]
    ]
    CallbackSet(radius_negative, w_zero, mphys_zero_cbs...)
end


function make_plume_initial_condition(env_profile; z0=300, dT=0.0, dqv=1.0e-3, r0=300, w0=1.0)
    T0 = ustrip(env_profile(z0 * u"m", :T)) + dT
    qv0 = ustrip(env_profile(z0 * u"m", :qv)) + dqv
    F = ComponentArray(r=r0, w=w0, T=T0, q_v=qv0, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0)
    return F
end


@testset "plume eqns" begin
    ICs = Dict(
        "simple" => ComponentArray(r=100u"m", w=1.0u"m/s", T=300u"K", q_v=13.0e-3, q_l=0.0, q_r=0.0, q_i=0.0, p=100e3u"Pa", q_pr=0.0),
        "rain cond/evap" => ComponentArray(r = 26.1u"m", w = 3.0u"m/s", T = 326.8u"K", q_v = 0.64, q_l = 0.019, q_r = 1.16e-7, q_i = 0.0, q_pr = 0.0, p = 98077.78139210862u"Pa")
    )
    @testset "$(name) state" for (name, F) in ICs
        dFdz = zero(F)
        env_profile = CloudModels.StandardIsentropicAtmosphere()
        p = (environment=env_profile, β=0.2)
        dFdz_initial = copy(dFdz)
        CloudModels.plume_equations!(dFdz, F, p, 0.0u"m")
        @test dFdz != dFdz_initial
        @test all(isfinite, dFdz)
    end

    @testset "tendencies don't depend on contents of output array" begin
        # ODE solvers reuse the derivative array between calls, so the
        # tendencies must overwrite (not accumulate into) it
        env_profile = CloudModels.ProfileRICO.RICO_profile()
        p = (environment=env_profile, β=0.2)
        F = make_plume_initial_condition(env_profile, r0=400, w0=1.0, z0=500)
        dFdz_a = zero(F)
        CloudModels.plume_equations!(dFdz_a, copy(F), p, 500.0)
        dFdz_b = zero(F) .+ 1.0
        CloudModels.plume_equations!(dFdz_b, copy(F), p, 500.0)
        @test dFdz_a == dFdz_b
    end
end


getvar(sol, v) = getindex.(sol.u, v)

@testset "plume-integration" begin
    env_profile = CloudModels.ProfileRICO.RICO_profile()
    z0 = 500.0
    z_max = 3500.0

    @testset "changing initial velocity" begin
        params = (environment=env_profile, β=0.2)
        sols = []
        for w0 in [0.5, 1.0, 2.0]
            U0 = make_plume_initial_condition(env_profile, r0=400, w0=w0, dqv=1.0e-3, z0=z0)
            prob = ODEProblem(CloudModels.plume_equations!, U0, [z0, z_max], params)
            sol = solve(prob, Euler(), dt=1.0, saveat=10.0, callback=setup_plume_callbacks())
            @test sum(getvar(sol, :q_l)) > 0.0
            push!(sols, sol)
        end

        # higher initial velocity should lead to higher cloud-top height
        ztop = [sol.t[end] for sol in sols]
        @test all(ztop .< z_max)
        @test sort(ztop) == ztop
    end

    @testset "changing entrainment rate" begin
        sols = []
        for β in [0.2, 0.1, 0.0]
            params = (environment=env_profile, β=β)
            U0 = make_plume_initial_condition(env_profile, r0=400, w0=1.0, dqv=1.0e-3, z0=z0)
            prob = ODEProblem(CloudModels.plume_equations!, U0, [z0, z_max], params)
            sol = solve(prob, Euler(), dt=1.0, saveat=10.0, callback=setup_plume_callbacks())
            push!(sols, sol)
        end

        # lower entrainment rate should lead to higher cloud-top height
        ztop = [sol.t[end] for sol in sols]
        @test sort(ztop) == ztop
    end
end
