# common functionality for integrating the equations, maybe these should go
# in the actual module somewhere...


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
    CallbackSet(cb_ztop, cb_zbottom, radius_negative, neutral_buoyancy_cloudtop, mphys_zero_cbs...)
end


function make_initial_condition(env_profile; z0=300, dT=0.0, dqv=1.0e-3, r0=300, w0=1.0)
    T0 = ustrip(env_profile(z0 * u"m", :T)) + dT
    qv0 = ustrip(env_profile(z0 * u"m", :qv)) + dqv
    F = ComponentArray(r=r0, w=w0, T=T0, q_v=qv0, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0, z=z0)
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
    # dF/dz is stored in `dFdt` vector, but need dF/dt: dF/dt = dF/dz * dz/dt
    dFdt .*= dzdt
    dFdt[:z] = dzdt
end