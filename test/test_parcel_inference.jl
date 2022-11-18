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


include("integration_common.jl")


function create_fake_observations(env_profile)
    prob_params = (; environment=env_profile)
    U0 = make_initial_condition(env_profile)
    prob = ODEProblem(parcel_equations!, U0, [0.0, 700.0], prob_params)
    sol = solve(prob, Euler(), saveat=10.0, callback=setup_callbacks(), dt=1.0)

    getvar(sol, v) = getindex.(sol.u, v)

    z_obs = getvar(sol, :z)
    w_true = getvar(sol, :w)
    
    w_noisy = w_true + rand(Normal(0.0, 1.0e-2), size(w_true))
    return z_obs, w_noisy, w_true
end


@model function discrete_cloud_evolution(z_obs, w_obs, env_profile)
    σ ~ InverseGamma(2, 3)
    σ_alt ~ InverseGamma(2, 3)

    z0 = z_obs[1]
    r0 ~ truncated(Normal(200.0, 100.0), lower=0.0)

    U0 = make_initial_condition(env_profile, r0=r0)
    prob_params = (; environment=env_profile)
    prob = ODEProblem(parcel_equations!, U0, [0.0, 700.0], prob_params)
    # pred = solve(prob, Euler(), saveat=10.0, callback=setup_callbacks(), dt=1.0)
    pred = solve(prob, RK4(), saveat=10.0, callback=setup_callbacks())

    # only compare up to the height of the max observation or prediction
    zmax_comp = minimum([pred.u[end].z, z_obs[end]])
    k = argmax(filter(z -> z <= zmax_comp, z_obs))
    z_comp = z_obs[1:k]
    w_pred = getindex.(pred.(z_comp), :w)
    w_obs[1:k] ~ MvNormal(w_pred, σ * I)

    # try to make prediction top center on obs top
    z_obs[end] ~ Normal(pred.u[end].z, σ_alt)

    return nothing
end

env_profile = CloudModels.ProfileRICO.RICO_profile();
z_obs, w_obs, w_true = create_fake_observations(env_profile)

plot(w_obs, z_obs)
plot!(w_true, z_obs)

model = discrete_cloud_evolution(z_obs, w_obs, env_profile);
chain = sample(model, MH(), 100)

# TODO work how to get other samplers to work...
# chain = sample(model, HMC(0.01, 5), 100)
# chain = sample(model, NUTS(), MCMCSerial(), 100, 1)
# chain = sample(model, NUTS(), MCMCThreads(), 100, 10)
# chain = sample(model, NUTS(0.45), MCMCThreads(), 100, 3; progress=false)

plot(chain)


posterior_samples = sample(chain[[:r0, ]], 50; replace=false)
plot(posterior_samples)
sols = []
for p in eachrow(Array(posterior_samples))
    r0 = p[1]
    @show r0

    U0 = make_initial_condition(env_profile, r0=r0)
    prob_params = (; environment=env_profile)
    prob = ODEProblem(parcel_equations!, U0, [0.0, 700.0], prob_params)
    sol_p = solve(prob, Euler(), saveat=10.0, callback=setup_callbacks(), dt=1.0)
    @show sol_p.u[end].z
    push!(sols, sol_p)
end

length(sols)
CloudModels.plot_profile_var(sols..., var_name=:w)
plot!(w_true, z_obs, color=:red, linewidth=5)