const p0vs = 611.2u"Pa"
const a0_lq = 17.67
const a1_lq = -32.19u"K"
const a0_ice = 22.587
const a1_ice = 0.7u"K"
const T00 = 273.15u"K"


function calc_pv_sat_liquid(T)
    return p0vs * exp((a0_lq * (T - T00) / (T + a1_lq)))
end

function calc_pv_sat_ice(T)
    return p0vs * exp((a0_ice * (T - T00) / (T + a1_ice)))
end

function calc_pv_sat(T)
    if T > T00
        return calc_pv_sat_liquid(T)
    else
        return calc_pv_sat_ice(T)
    end
end

function calc_qv_sat(T, p)
    pv_sat = calc_pv_sat(T)
    epsilon = R_d / R_v
    qv_sat = (epsilon * pv_sat) / (p - (1.0 - epsilon) * pv_sat)

    return qv_sat
end


"""
from G. Thompson '07 microphysics scheme
"""
function calc_dynamic_viscosity_thompson(T)
    Tc = T - T00

    if T > T00
        return (1.718 + 0.0049 * Tc) * 1.0e-5
    else
        return (1.718 + 0.0049 * Tc - 1.2e-5 * Tc ^ 2.0) * 1.0e-5
    end
end

"""
Sources:
    Rogers & Yau 1989
"""
function calc_dynamic_viscosity(T)
    return 1.72e-5 * (393.0u"K" / (T + 120.0u"K")) * (T / T00) ^ (3.0 / 2.0) * u"kg/m/s"
end



"""
Linear model of thermal conductivity in [J/m/s/K]

takes in temperature temperature in [K]
"""
function calc_thermal_conductivity_coefficient(T)
    a_K = 8.0e-5u"J/m/s/K^2"
    b_K = 2.4e-2u"J/m/s/K"

    return a_K * (T - T00) + b_K
end


"""
Coefficient of water vapour diffusion [m^2/s]
"""
function calc_water_vapour_diffusion_coefficient(T, p)
    # ATHAM_constants = {
    #     "a": 2.11e-5,
    #     "b": 1.94,
    # }

    # these constants were obtained by fitting against data in Rogers & Yau
    a = 2.20e-5u"m^2/s"
    b = 1.92

    # tabulated values in Rogers & Yau are given for p=100kPa=100000Pa reference pressure,
    # have to scale by pressure get the correct diffusivity
    p0 = 100e3u"Pa"

    return a * (T / T00) ^ b * p0 / p
end