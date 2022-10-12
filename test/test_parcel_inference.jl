using Test
using CloudModels
#using Unitful
using OrdinaryDiffEq
using ComponentArrays
using Plots
using Revise
using Turing
using StatsPlots

function setup_callbacks()
    condition(u, t, integrator) = u.z - 2000.0
    cb_ztop = ContinuousCallback(condition,terminate!)
    cb_zbottom = ContinuousCallback((u, t, integrator) -> u.z - 200.0,terminate!)

    mphys_zero_cbs = [
        DiscreteCallback((u, t, i) -> u[v] < 0.0, x-> terminate!(x, :NegativeMphys)) for v in [:q_v, :q_r, :q_l]
    ]
    CallbackSet(cb_ztop, cb_zbottom, mphys_zero_cbs...)
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

function parcel_equations!(dFdt, F, params, t)
    dFdt .= 0.0
    dzdt = F.w
    CloudModels.parcel_equations!(dFdt, F, F.z * u"m", params)
    dFdt ./= dzdt
    dFdt[:z] = dzdt
end


env_profile = CloudModels.ProfileRICO.RICO_profile();
params = (environment=env_profile, β=0.2)
F = make_initial_condition(env_profile, 300.0)

#prob = ODEProblem(nounit_parcel_equations!, F, [0.0, 1000.0], params)
prob = ODEProblem(parcel_equations!, F, [0.0, 1000.0], params)
#sol = solve(prob, Euler(), saveat=10.0, callback=setup_callbacks(), dt=1.0)
sol = solve(prob, RK4(), saveat=10.0, callback=setup_callbacks())
CloudModels.plot_profile(sol)
length(sol)


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