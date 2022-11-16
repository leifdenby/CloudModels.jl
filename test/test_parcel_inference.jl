using Test
using CloudModels
#using Unitful
using OrdinaryDiffEq
using ComponentArrays
using Plots
using Revise
using Turing
using StatsPlots


include("integration_common.jl")


#prob = ODEProblem(nounit_parcel_equations!, F, [0.0, 1000.0], params)

env_profile = CloudModels.ProfileRICO.RICO_profile();

sols = []
for w0 in [0.5, 1.0, 2.0]
    params = (environment=env_profile, β=0.2)
    U0 = make_initial_condition(env_profile, r0=400, w0=w0, dqv=1.0e-3, z0=500.0)
    prob = ODEProblem(parcel_equations!, U0, [0.0, 700.0], params)
    sol = solve(prob, Euler(), dt=1.0, saveat=10.0, callback=setup_callbacks(z_max=00.0))
    push!(sols, sol)
end

CloudModels.plot_profile(sols..., wrap=4, size=300)

CloudModels.plot_profile(sols..., vars=[:r, :w, :qv, :rh, :Δrho, :T], wrap=3, size=300)

CloudModels.plot_profile_var(sols..., var_name=:ql)
CloudModels.plot_profile_var(sols..., var_name=:r)


plot!(p, rand(Float32, (100,)), sin)

z_ = 500:5:1e3
plot(env_profile.(z_, :T), z_)


p = plot(rand(Float32, (10)), sin)
typeof(p)
fieldnames(typeof(p))

@model function discrete_cloud_evolution(data)
    dT = 0.0
    dqv = 3.0e-3
    z0 = 300.0
    T0 = ustrip(env_profile(z0 * u"m", :T)) + dT
    qv0_μ = ustrip(env_profile(z0 * u"m", :qv)) + dqv
    #qv0 ~ truncated(Normal(qv0_μ, 1.0e-3); lower=0.0)
    #r0 ~ truncated(Normal(280.0, 20.0), lower=0.0)
    r0 ~ truncated(Normal(280.0, 100.0), lower=100.0)
    F = ComponentArray(r=r0, w=1.0, T=T0, q_v=qv0_μ, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0, z=z0)
    
    σ ~ InverseGamma(2, 3)

    prob = ODEProblem(parcel_equations!, F, [0.0, 1000.0], params)
    pred = solve(prob, Euler(), saveat=10.0, callback=setup_callbacks(), dt=1.0)
    #pred = solve(prob, RK4(), u0=F; saveat=10.0, callback=setup_callbacks())
    
    @show length(pred) r0
    
    # fit only on vertical velocity
    for i in 1:length(data)
        if data[i].z > sol.u[end].z
            break
        end
        # interpolate into solution at heights where we have "observations"
        data[i].w ~ Normal(pred(data[i].z).w, σ)
    end
end

g(sol, v) = getindex.(sol.u, v)
g(sol, :z)

model = discrete_cloud_evolution(sol.u);
chain = sample(model, MH(), 50)
chain = sample(model, HMC(0.01, 5), 100)
chain = sample(model, NUTS(), MCMCSerial(), 100, 1)

plot(chain)

sol_w = getindex.(sol.u, :w)
plot(sol.t, sol_w, marker=:circle)

posterior_samples = sample(chain[[:r0, ]], 10; replace=false)
p2 = plot()
for p in eachrow(Array(posterior_samples))
    dT = 0.0
    dqv = 3.0e-3
    z0 = 300.0
    T0 = ustrip(env_profile(z0 * u"m", :T)) + dT
    qv0_μ = ustrip(env_profile(z0 * u"m", :qv)) + dqv
    r0 = p[1]

    F = ComponentArray(r=r0, w=1.0, T=T0, q_v=qv0_μ, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0, z=z0)
    prob = ODEProblem(parcel_equations!, F, [0.0, 1000.0], params)
    sol_p = solve(prob, Euler(), u0=F, saveat=10.0, callback=setup_callbacks(), dt=1.0)
    @show sol_p.u[1]
    plot!(p2, g(sol_p, :z), g(sol_p, :w); alpha=0.1, color="#BBBBBB")
end
plot!(p2, g(sol, :z), g(sol, :w), color=:red)
plot!(p2)