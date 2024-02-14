using Test
using CloudModels
#using Unitful
using OrdinaryDiffEq
using ComponentArrays
using Plots
using Turing
using StatsPlots
using Distributions
using LinearAlgebra

function setup_callbacks(; stop_on_descent=true)
    cbs = []
    # don't go below z=0m
    push!(cbs, ContinuousCallback((u, z, integrator) -> z, x-> terminate!(x, :NegativeAltitude)))
    # don't allow radius to become negative
    push!(cbs, ContinuousCallback((u, z, integrator) -> u.r, x-> terminate!(x, :NegativeRadius)))

    if stop_on_descent
        # stop of vertical velocity becomes negative
        push!(cbs, ContinuousCallback((u, z, integrator) -> u.w, x -> terminate!(x, :NegativeVerticalVelocity)))
    end
    
    push!(cbs,[
        DiscreteCallback((u, z, i) -> u[v] < 0.0, x-> terminate!(x, :NegativeMphys))
        for v in [:q_v, :q_r, :q_l]
    ]...)

    CallbackSet(cbs...)
end


function make_initial_condition(env_profile; z0=300, dT=0.0, dqv=1.0e-3, r0=300, w0=1.0)
    T0 = ustrip(env_profile(z0 * u"m", :T)) + dT
    qv0 = ustrip(env_profile(z0 * u"m", :qv)) + dqv
    p0 = ustrip(env_profile(z0 * u"Pa", :p))
    F = ComponentArray(r=r0, w=w0, T=T0, q_v=qv0, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=p0, z=z0)
    return F
end

function setup_problem(r0, env_profile; t_max=2000.0)
    prob_params = (; environment=env_profile)
    U0 = make_initial_condition(env_profile, r0=r0)
    prob = ODEProblem(CloudModels.parcel_equations!, U0, [0.0, t_max], prob_params)
    return prob
end

function create_fake_observations(r0, env_profile)
    prob = setup_problem(r0, env_profile)
    sol = solve(prob, Euler(), saveat=10.0, callback=setup_callbacks(), dt=1.0)

    getvar(sol, v) = getindex.(sol.u, v)

    t_obs = sol.t
    z_obs = getvar(sol, :z)
    w_true = getvar(sol, :w)
    
    w_noisy = w_true + rand(Normal(0.0, 1.0e-2), size(w_true))
    return t_obs, z_obs, w_noisy, w_true
end


@model function discrete_cloud_evolution(t_obs, w_obs, env_profile)
    # σ ~ InverseGamma(2, 3)
    # σ_alt ~ InverseGamma(2, 3)

    # r0 ~ Normal(300.0, 100.0)
    # r1 ~ Normal(3.0, 1.0)
    # r1 ~ TruncatedNormal(3.0, 1.0, 1.0, 10.0)
    # r0 = r1 * 100.0
    r0 ~ TruncatedNormal(200.0, 10.0, 100.0, 1000.0)
    
    prob = setup_problem(r0, env_profile)
    pred = nothing
    try
        pred = solve(prob, Euler(), saveat=10.0, callback=setup_callbacks(), dt=10.0)
    catch
        Turing.@addlogprob! -Inf
        return
    end
    
    # if length(pred) != length(t_obs)
        # Turing.@addlogprob! -Inf
        # return
    # end

    # w_pred = getindex.(pred.(t_obs), :w)
    # w_obs ~ MvNormal(w_pred, σ * I)
    t_max = pred.t[end]
    n = 0
    
    for i in 1:length(t_obs)
        n += 1
        # @show i t_obs[i] pred.t[end]
        if t_obs[i] > pred.t[end]
            break
        end
        w_obs[i] ~ Normal(pred(t_obs[i]).w, 1.0e-2)
    end
    
    #@show n r0 t_max

    return nothing
end

# using DynamicHMC

env_profile = CloudModels.ProfileRICO.RICO_profile();
t_obs, z_obs, w_obs, w_true = create_fake_observations(300, env_profile)
scatter(w_obs, t_obs)

model = discrete_cloud_evolution(t_obs, w_obs, env_profile);
chain = sample(model, MH(), MCMCThreads(), 3000, 8)
plot(chain)

# TODO work how to get other samplers to work...
chain = sample(model, HMC(0.0005, 5), 500)
plot(chain)
# chain = sample(model, NUTS(), 200)
# chain = sample(model, DynamicNUTS(), 1000)
chain = sample(model, NUTS(), MCMCThreads(), 1000, 8)
# chain = sample(model, HMCDA(200, 0.65, 0.3), 100)
plot(chain)
# chain = sample(model, HMC(0.0001, 5), MCMCThreads(), 200, 8)
# chain = sample(model, NUTS(), MCMCThreads(), 100, 10)
# chain = sample(model, NUTS(0.45), MCMCThreads(), 100, 3; progress=false)



posterior_samples = sample(chain[[:r0, ]], 50; replace=false)
plot(posterior_samples)
posterior_samples = [250., 300.0]
sols = []
for p in eachrow(Array(posterior_samples))
    r0 = p[1]
    @info r0

    prob = setup_problem(r0, env_profile)
    sol_p = solve(prob, Euler(), saveat=10.0, callback=setup_callbacks(), dt=1.0)
    push!(sols, sol_p)
end

length(sols)
CloudModels.plot_profile_var(sols..., var_name=:w, color=:grey)
plot!(w_true, z_obs, color=:red, linewidth=5)

getvar(sol, v) = getindex.(sol.u, v)
scatter(w_obs, t_obs)
for i in 1:length(sols)
    plot!(getvar(sols[i], :w), sols[i].t, color=:grey, label=nothing, alpha=0.1)
end
plot!(w_true, t_obs, color=:red, linewidth=5)




prob = setup_problem(450, env_profile)
sol_p = solve(prob, Euler(), saveat=10.0, callback=setup_callbacks(stop_on_descent=true), dt=1.0)

CloudModels.plot_profile(sol_p)