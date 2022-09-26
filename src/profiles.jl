struct ReferenceAtmosphere
    rho0::typeof(1.0u"kg/m^3")
    p0::typeof(1.0u"Pa")
    dTdz::typeof(1.0u"K/km")
end

rho0 = 1.205 * u"kg/m^3"
p0 = 101325.0 * u"Pa"
T0 = p0 / (rho0 * R_d) |> u"K"
κ = R_d / cp_d

StandardIsothermalAtmosphere = ReferenceAtmosphere(rho0, p0, 0.0u"K/km")
StandardIsentropicAtmosphere = ReferenceAtmosphere(rho0, p0, -g/cp_d)
ConstantDensityAtmospere = ReferenceAtmosphere(rho0, p0, -g/R_d)

function calc_temperature(z, prof::ReferenceAtmosphere)
    return T0 + prof.dTdz * z
end

function calc_density(z, prof::ReferenceAtmosphere)
    if prof.dTdz == 0.0
        return rho0 * exp(-z * g * R_d * T0)
    else
        α = g / (prof.dTdz * R_d)
        T = calc_temperature(z, prof)
        return rho0 * T0 ^ (α + 1.0) * T ^ (- α - 1.0)
    end
end

function calc_pressure(z, prof::ReferenceAtmosphere)
    ρ = calc_density(z, prof)
    T = calc_temperature(z, prof)
    return ρ * R_d * T
end

function calc_potential_temperature(z, prof::ReferenceAtmosphere)
    T = calc_temperature(z, prof)
    p = calc_pressure(z, prof)
    return T * (p / p0) ^ (-κ)
end