gamma=0.5
C_D=0.506
l_pr=100.0u"m"


# droplet-size distribution constant
N0 = 1.0e7u"m^-4"  # [m^-4]

"""
Momentum equation

State variables:
    w: verticaly velocity
    r: cloud radius
    rho_c: cloud density
"""
function calc_dw_dz(p, w_c, r_c, T_c, qd_c, qv_c, ql_c, qi_c, qr_c, rho_e, mu)
    rho_c = calc_mixture_density(p, T_c, qd_c, qv_c, ql_c, qi_c, qr_c)
    B = (rho_e - rho_c) / rho_e

    return (
        1.0
        / w_c
        * (
            g / (1.0 + gamma) * B
            - mu * w_c ^ 2.0
            - 3.0 / 8.0 * C_D * w_c ^ 2.0 / r_c
        )
    )
end

function calc_dT_dz(
    r_c,
    T_c,
    qd_c,
    qv_c,
    ql_c,
    qi_c,
    dql_c__dz,
    dqi_c__dz,
    dqv_c__dz,
    T_e,
    qv_e,
    mu,
)
    """
    State variables:
        qd_c, qv_c, ql_c, qi_c: in-cloud dry air, water vapour, liquid water, ice
        T_c: in-cloud absolute temp

        dql_c__dz, qdi_c__dz: vertical in

        qd_c, qv_c, ql_c, qi_c: environment (constant) dry air, water vapour, liquid water, ice
        T_e: environment absolute temp
    """

    # latent-heat corrected for change in temperature
    L_v_ = L_v + (cp_v - cp_l) * (T_c - T00)
    L_s_ = L_s + (cp_v - cp_i) * (T_c - T00)
    L_f_ = L_v + (cp_l - cp_i) * (T_c - T00)

    # equations below don't work with ice yet
    if qi_c > 0.0
        throw("ice not implemented yet")
    end

    c_cm_p = cp_d * qd_c + cp_l * (qv_c + ql_c + qi_c)

    qd_e = 1.0 - qv_e
    ql_e, qi_e, qr_e = 0.0, 0.0, 0.0

    c_em_p = cp_d * qd_e + cp_l * (qv_e + ql_e + qi_e)

    # difference in environment and cloud moist static energy
    Ds = (c_em_p * T_e + qv_e * L_v_ - qi_e * L_f_) - (
        c_cm_p * T_c + qv_c * L_v_ - qi_c * L_f_
    )

    dqd_c__dz = -(dqv_c__dz + dql_c__dz + dqi_c__dz)

    # heat changes due to phase changes (and entrainment of species)
    dedz_q = (
        T_c * cp_d * dqd_c__dz
        + (T_c * cp_l + L_v_) * dqv_c__dz
        + T_c * cp_l * dql_c__dz
        + (T_c * cp_i - L_f_) * dqi_c__dz
    )

    # *actual* mixture heat capacity (if temperature dependency of
    # latent heats are considered in derivation then the *actual*
    # mixture heat capacity comes out)
    c_cm = cp_d * qd_c + qv_c * cp_v + ql_c * cp_l + qi_c * cp_i

    return -g / c_cm + mu * Ds / c_cm - 1.0 / c_cm * dedz_q
end

function calc_dr_dz(
    p,
    w_c,
    r_c,
    T_c,
    qd_c,
    qv_c,
    ql_c,
    qr_c,
    qi_c,
    dqv_c__dz,
    dql_c__dz,
    dqi_c__dz,
    dTc_dz,
    dw_dz,
    mu,
)
    """
    Mass conservation equation

    State variables:
        qc_g: specific concentration of gas species in-cloud
        T_c: absolute temperature in cloud
        rho_c: density in-cloud
        rho_cg: density of gas mixture (wrt gas volume)

        (assumed constant)
        rho_i: in-cloud density of ice
        rho_l: in-cloud density of liquid water
    """

    # Total specific concentration stays unchanged
    dqd_c__dz = -(dqv_c__dz + dql_c__dz + dqi_c__dz)

    # in-cloud mixture density
    rho_c = calc_mixture_density(
        p, T_c, qd_c, qv_c, ql_c, qi_c, qr_c
    )

    # in-cloud gas density
    rho_cg = calc_mixture_density(p, T_c, qd_c, qv_c, 0.0, 0.0, 0.0)

    # effective gas constant
    Rs_c = (R_v * qv_c + R_d * qd_c) / (qv_c + qd_c)

    # in-cloud specific constant of gas constituents
    qg_c = qd_c + qv_c

    return (
        r_c
        / 2.0
        * (
            qg_c * rho_c / rho_cg * rho_c / rho_cg * g / (Rs_c * T_c)
            + qg_c * rho_c / rho_cg * 1.0 / T_c * dTc_dz
            + rho_c / (rho_cg * Rs_c) * (dqv_c__dz * R_v + dqd_c__dz * R_d)
            + rho_c / rho_i * dqi_c__dz
            + rho_c / rho_l * dql_c__dz
            + mu  #  1./M*dM_dz\
            - 1.0 / w_c * dw_dz
        )
    )
end

"""
Estimate rate at which rain-droplets leave the cloudy air parcel. Based
on the relative velocity of the cloud parcel and the fall-speed of the
rain-droplets
"""
function calc_dqr_dz__rainout(rho_c, q_r, w)
    if q_r <= 0.0
        return 0.0u"1/m"
    end

    # size-distribtion length-scale
    l = (8.0 * rho_l * pi * N0 / (q_r * rho_c)) ^ 0.25

    # fall-speed coefficient taken from the r > 0.5mm expression for
    # fall-speed from Herzog '98
    a_r = 201.0u"m^.5 / s"
    # reference density
    rho0 = 1.12u"kg/m^3"

    # charateristic velocity, XXX: this is definitely wrong, but couldn't
    # work out how to do the integral at the time
    w_r = a_r * sqrt(1.0 / l * rho0 / rho_c)

    # rainout fraction
    f = w_r / w / l_pr

    return f * q_r
end


function parcel_equations!(dFdz, F, z, params)
#    @info "what"
#    @show dFdz
#    @show F
    try
    environment = params[:environment]
    r = F[:r]
    w = F[:w]
    T = F[:T]
    q_v = F[:q_v]
    q_l = F[:q_l]
    q_i = F[:q_i]
    q_r = F[:q_r]
    q_d = 1.0 - q_v - q_l - q_i - q_r

    # cloud is assumed to be at same pressure as in environment
    p = environment(z, :p)
    # NB: need to make sure that pressure is set since the microphysics
    # needs it and we don't integrate a pressure gradient, instead we
    # assume pressure balance with environment
    F[:p] = p

    # rho_e = self.environment.rho(z)
    T_e = environment(z, :T)
    qv_e = environment(z, :qv)

    qd_e = 1.0 - qv_e
    rho_e = calc_mixture_density(p, T_e, qd_e, qv_e, 0.0, 0.0, 0.0)

    # calculate entrainment rate
    rho_c = calc_mixture_density(p, T, q_d, q_v, q_l, q_i, q_r)
    B = (rho_e - rho_c) / rho_e
    mu = params.β / r

    # 1. Estimate change in vertical velocity with initial state
    dwdz = calc_dw_dz(p, w, r, T, q_d, q_v, q_l, q_r, q_i, rho_e, mu)

    dFdz[:w] = dwdz

    # assume no condesates present in environment
    ql_e, qr_e, qi_e = 0.0, 0.0, 0.0

    # as well as tracer (water) changes from microphysics we have to consider entrainment
    # dFdz_entrain__q = zero(dFdz)
    dFdz_entrain__q = copy(dFdz) .* 0.0

    # dqd_dz__ent = mu/rho_c*(qd_e*rho_e - q_d*rho_c)
    dqd_dz__ent = mu * (qd_e - q_d)
    # dFdz_entrain__q[:q_v] = -q_v*dqd_dz__ent
    # dFdz_entrain__q[:q_l] = -q_l*dqd_dz__ent
    # dFdz_entrain__q[:q_r] = -q_r*dqd_dz__ent
    # dFdz_entrain__q[:q_i] = -q_i*dqd_dz__ent

    # this formulation gives a strong impact of entrainment on
    # hydrometeors, seems best I think
    dFdz_entrain__q[:q_v] = mu * (qv_e - q_v)
    dFdz_entrain__q[:q_l] = mu * (ql_e - q_l)
    dFdz_entrain__q[:q_r] = mu * (qr_e - q_r)
    dFdz_entrain__q[:q_i] = mu * (qi_e - q_i)

    # Old formulation: not exactly sure why this doesn't work, but the
    # above is from the perspective of entraining dry air instead of
    # detraining hydrometeors, maybe there's something important here?

    # dFdz_entrain__q[:q_v] = mu/rho_c*(qv_e*rho_e - q_v*rho_c)
    # dFdz_entrain__q[:q_l] = mu/rho_c*(ql_e*rho_e - q_l*rho_c)
    # dFdz_entrain__q[:q_r] = mu/rho_c*(qr_e*rho_e - q_r*rho_c)
    # dFdz_entrain__q[:q_i] = mu/rho_c*(qi_e*rho_e - q_i*rho_c)

    # 2. estimate new state from phase changes predicted by microphysics

    # dFdt_micro = ComponentArray(Dict([ v => 0.0 * unit(F[v] * u"1/s") for v in keys(F)]))
    dFdt_micro = copy(dFdz) .* 0.0
    dFdt_microphysics!(dFdt_micro, F, 0.0u"s")

    dFdz_micro = dFdt_micro / w  # w = dz/dt
    # temperature effect will be determined from temperature equation,
    # microphysics only provides changes due to microphysics itself over short
    # time-increment of rise, so we only include the changes to the microphysics
    # species
    for v in [:q_v, :q_l, :q_r, :q_i]
        dFdz[v] += dFdz_micro[v]
        dFdz[v] += dFdz_entrain__q[v]
    end

    # 3. Estimate temperature change forgetting about phase-changes for now
    # (i.e. considering only adiabatic adjustment and entrainment)

    # in terms of the thermodynamics cloud-water and rain-water are treated identically
    dql_c__dz = dFdz[:q_l] + dFdz[:q_r]
    dqi_c__dz = dFdz[:q_i]
    dqv_c__dz = dFdz[:q_v]
    ql_c = q_l + q_r

    dTdz_s = calc_dT_dz(r, T, q_d, q_v, ql_c, q_i, dql_c__dz, dqi_c__dz, dqv_c__dz, T_e, qv_e, mu)

    dFdz[:T] += dTdz_s
    dTdz = dFdz[:T]

    # 4. Use post microphysics state (and phase changes from microphysics) to estimate radius change
    drdz = calc_dr_dz(p, w, r, T, q_d, q_v, q_l, q_r, q_i, dqv_c__dz, dql_c__dz, dqi_c__dz, dTdz, dwdz, mu)

    dFdz[:r] = drdz

    # 5. Estimate fraction of rain that leaves parcel
    dqr_dz__rainout = calc_dqr_dz__rainout(rho_c, q_r, w)
    dFdz[:q_r] -= dqr_dz__rainout
    dFdz[:q_pr] = dqr_dz__rainout
    catch
        @warn F
        rethrow()
    end
end