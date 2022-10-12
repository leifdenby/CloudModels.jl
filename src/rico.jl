module ProfileRICO

using ComponentArrays
using LinearInterpolations
using OrderedCollections

macro u_str(unit)
    return 1.0
end

uconvert(_, v) = v
unit(v) = 1.0

"""
Based on KNMI's synthesis of the RICO field compaign for a LES intercomparison study

http://projects.knmi.nl/rico/setup3d.html

OBS: The sea surface temperature is different from the temperature at z=0m,
which is necessary to have surface fluxes.
"""

# surface conditions
ps = 101540.0u"Pa"  # [Pa], surface pressure
p0 = 100000.0u"Pa"  # [Pa], reference pressure
Ts = 299.8u"K"  # [K], sea surface temperature

# constants
Lv = 2.5e6u"J/kg"  # [J/kg]
cp_d = 1005.0u"J/kg/K" # [J/kg/K]
g = 9.81u"m/s^2"  # [m/s^2]
R_d = 287.0u"J/kg/K"  # [J/kg/K]

# XXX: R_v and cp_v are not given in the RICO test definition on the
# the KNMI site I will use what I believe are standard values here
R_v= 461.51u"J/kg/K"
cp_v= 1859.0u"J/kg/K"

reference_points = OrderedDict(
    0.0u"m" => ComponentArray(qt=16.0e-3, theta_l=297.9u"K"),
    740.0u"m" => ComponentArray(qt=13.8e-3, theta_l=297.9u"K"),
    3260.0u"m" => ComponentArray(qt=2.4e-3, theta_l=317.0u"K"),
    4000.0u"m" => ComponentArray(qt=2.4e-3, theta_l=317.0u"K"),
)

"""Create a vertical profile that we can interpolate into later.
Integrating with the hydrostatic assumption.
"""
function integrate_profile(ref_profile, dz)
    z_max = 4e3u"m"

    z = 0.0u"m"
    p = ps

    # Cathy suggested using the liquid water potential temperature as the
    # temperature in the first model level
    T = ref_profile(0.0u"m")[:theta_l]

    profile = []

    points = OrderedDict()

    while z < z_max
        qt, theta_l = ref_profile(z)

        # assume no liquid water
        qv = qt
        qd = 1.0 - qt

        R_l = R_d * qd + R_v * qv
        c_l = cp_d * qd + cp_v * qv

        T = uconvert(u"K", theta_l / ((p0 / p) ^ (R_l / c_l)))
        rho = uconvert(u"kg/m^3", 1.0 / ((qd * R_d + qv * R_v) * T / p))  # + 1.0/(ql/rho_l), ql = 0.0

        points[z] = ComponentArray(rho=rho, p=p, T=T)

        # integrate pressure
        z += dz
        p += uconvert(u"Pa", -rho * g * dz)
    end

    return points
end

struct RICO_profile
    ref_profile
    itp2
end


function RICO_profile()
    ref_profile = Interpolate(collect(keys(reference_points)), collect(values(reference_points)))
    points_integrated = integrate_profile(ref_profile, 100u"m");
    itp2 = Interpolate(collect(keys(points_integrated)), collect(values(points_integrated)))
    return RICO_profile(ref_profile, itp2)
end


function (prof::RICO_profile)(z, var_name::Symbol)
    if var_name in [:theta_l, :qt]
        return prof.ref_profile(z)[var_name]
    end
    if var_name in [:rho, :p, :T]
        return prof.itp2(z)[var_name]
    end
    if var_name == :qv
        return prof.ref_profile(z)[:qt]
    end
    throw("Can't compute $(var_name)")
end

end