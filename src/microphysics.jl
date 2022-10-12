N0i = 200 * 1.0e6u"1/m^3"  # initial aerosol number concentration [m-3]
r0 = 0.1e-6u"m"  # cloud droplet initial radius

model_constraint = "isobaric"

"""
Collection of microphysics routines for use with cloud-model integration.
"""
function calc_cp_m(F)
    qv = F[:q_v]
    ql = F[:q_l]
    qr = F[:q_r]
    qi = F[:q_i]
    qd = 1.0 - qv - ql - qr - qi

    if qi > 0.0
        throw("not implemented")
    end

    return cp_d * qd + (ql + qr) * cp_l + qv * cp_v
end

function calc_cv_m(F)
    qv = F[:q_v]
    ql = F[:q_l]
    qr = F[:q_r]
    qi = F[:q_i]
    qd = 1.0 - qv - ql - qr - qi

    if qi > 0.0
        throw("not implemented")
    end

    return cv_d * qd + (ql + qr) * cv_l + qv * cv_v
end

"""
Create rain droplets through collision and coalescence of cloud
droplets

from Cloud and Precipitation Microphysics by J. Straka 2009, eqn 9.2, p. 256
"""
function _dqr_dt__autoconversion(ql, qg, rho_g)
    k_c = 1.0e-3u"1/s"
    a_c = 5.0e-4

    ## XXX: I think these expressions are dimensionally incorrect compared to Straka's book
    # dqr_dt = k_c * (ql - qg / rho_g * a_c)
    # TODO: Is it reasonable to limit this? It seems the autoconversion
    # equation doesn't really make sense to describe breakup
    # dqr_dt = max(0.0, dqr_dt)
    # return dqr_dt
    H(x) = x > 0 ? x : 0

    if ql == 0.0
        return 0.0u"1/s"
    else
        return k_c * H(ql - a_c)
    end
end

function _dqr_dt__accretion(ql, qg, qr, rho_g)
    # if there is no rain there is nothing to accrete onto
    if qr <= 0.0
        return 0.0u"1/s"
    end

    # TODO: rederive these equations to verify that they are correct
    G3p5 = 3.32399614155  # = Gamma(3.5)
    N0r = 1.0e7u"1/m^4"
    a_r = 201.0u"m^(1/2)/s"
    rho0 = 1.12u"kg/m^3"

    λr = (pi * (qg * rho_l) / (qr * rho_g) * N0r) ^ (1.0 / 4.0)

    dqr_dt = (
        pi
        / 4.0
        * N0r
        * a_r
        * sqrt(rho0 / rho_g)
        * G3p5
        * λr ^ (-3.5)
        * ql
    )

    return max(dqr_dt, 0.0u"1/s")
end

function _dql_dt__cond_evap(qv, ql, rho, p, T)
    # condensation evaporation of cloud droplets (given number of droplets
    # and droplet radius calculated above)
    qv_sat = calc_qv_sat(T, p)
    Sw = qv / qv_sat
    
    if ql < 1.0e-12
        return 0.0
    end

    # number of aerosols stays constant (to add aerosol activation only a
    # fraction if the original present aerosols would be "activated" at
    # this point)
    Nc = N0i

    if ql == 0.0
        if Sw > 1.0
            r_c = r0
        else
            r_c = 0.0u"m"
        end
    else
        @show ql rho
        r_c = (ql * rho / (4.0 / 3.0 * pi * N0i * rho_l)) ^ (1.0 / 3.0)
        # droplet's should at least be as big as their initial (aerosol) size
        r_c = max(r0, r_c)
    end

    Ka = calc_thermal_conductivity_coefficient(T)
    Fk = (L_v / (R_v * T) - 1.0) * L_v / (Ka * T) * rho_l

    pv_sat = calc_pv_sat(T)
    Dv = calc_water_vapour_diffusion_coefficient(T, p)
    Fd = R_v * T / (pv_sat * Dv) * rho_l

    # compute rate of change of condensate from diffusion
    dql_dt = 4 * pi * rho_l / rho * Nc * r_c * (Sw - 1.0) / (Fk + Fd)

    return dql_dt
end

"""
Condensation and evaporation of rain. Similar to cloud-water droplet
condensation/evaporation but includes corrections for "ventilation" and
and droplet-size is assumed to follow a Marshall-Palmer distribution:

    N(r)dr = N0 exp(-l*r) dr

Arguments:
    rho: cloud-mixture density
    qv: water vapour specific mass
    qr: rain water specific mass
    T: cloud-mixture temperature
    p: pressure
"""
function _dqr_dt__cond_evap(qv, qr, rho, p, T)
    G2p75 = 1.608359421985546  # = Gamma(2.75)

    # droplet-size distribution constant
    N0 = 1.0e7u"1/m^4"  # [m^-4]

    # fall-speed coefficient taken from the r > 0.5mm expression for
    # fall-speed from Herzog '98
    a_r = 201.0u"m^(1/2)/s"
    # reference density
    rho0 = 1.12u"kg/m^3"

    # can't do cond/evap without any rain-droplets present
    if qr <= 0.0
        return 0.0u"1/s"
    end

    # computer super/sub-saturation
    qv_sat = calc_qv_sat(T, p)
    Sw = qv / qv_sat

    # size-distribtion length-scale [1/m]
    λ = (8.0 * rho_l * pi * N0 / (qr * rho)) ^ 0.25

    # w_r = a_r * sqrt(1.0 / λ * rho0 / rho)
    # Nr = N0 / λ

    # air condutivity and diffusion effects
    Ka = calc_thermal_conductivity_coefficient(T)
    Fk = (L_v / (R_v * T) - 1) * L_v / (Ka * T) * rho_l

    pv_sat = calc_pv_sat(T)
    Dv = calc_water_vapour_diffusion_coefficient(T, p)
    Fd = R_v * T / (pv_sat * Dv) * rho_l

    # compute the ventilation coefficient `f`
    # fall-velocity
    # dynamic viscosity
    mu = calc_dynamic_viscosity(T)

    f = 1.0 + 0.22 * (2.0 * a_r * rho / mu) ^ 0.5 * (rho0 / rho) ^ 0.25 * G2p75 / (λ ^ (3/4))

    # compute rate of change of condensate from diffusion
    dqr_dt = 4 * pi * rho_l / rho * N0 / λ ^ 2.0 * (Sw - 1.0) / (Fk + Fd) * f

    return dqr_dt
end


function dFdt_microphysics!(dFdt, F, t)
    qv = F[:q_v]
    ql = F[:q_l]
    qr = F[:q_r]
    qi = F[:q_i]

    # the integrator may have put us below zero, this is non-conservative,
    # but I don't have any other options here since I can't modify the
    # integrator's logic
    if ql < 0.0
        ql = 0.0
        F[:q_l] = ql
    end

    qd = 1.0 - qv - ql - qr - qi
    T = F[:T]
    p = F[:p]

    # mixture density
    rho = calc_mixture_density(p, T, qd, qv, ql, qi, qr)
    # gas density
    rho_g = calc_mixture_density(p, T, qd, qv, 0.0, 0.0, 0.0)

    dql_dt = _dql_dt__cond_evap(qv, ql, rho, p, T)

    qg = qv + qd
    dqr_dt_autoc = _dqr_dt__autoconversion(ql, qg, rho_g)
    dqr_dt_accre = _dqr_dt__accretion(ql, qg, qr, rho_g)
    dqr_dt_condevap = _dqr_dt__cond_evap(qv, qr, rho, p, T)

    dqr_dt = dqr_dt_autoc + dqr_dt_accre
    
    dFdt[:q_l] = dql_dt - dqr_dt
    dFdt[:q_v] = -dql_dt - dqr_dt_condevap
    dFdt[:q_r] = dqr_dt + dqr_dt_condevap


    for v in [:q_l, :q_v, :q_r]
        if F[v] < 0.0
            @warn dFdt v
        end
    end

    if model_constraint == "isometric"
        c_m = calc_cv_m(F)
    elseif model_constraint == "isobaric"
        c_m = calc_cp_m(F)
    else
        throw("Model constraint mode '$(model_constraint)%s' not implemented")
    end

    dFdt[:T] = L_v / c_m * dql_dt
end