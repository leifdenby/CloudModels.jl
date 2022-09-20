using ComponentArrays
using LinearInterpolations
using OrderedCollections

"""
Based on KNMI's synthesis of the RICO field compaign for a LES intercomparison study

http://projects.knmi.nl/rico/setup3d.html

OBS: The sea surface temperature is different from the temperature at z=0m,
which is necessary to have surface fluxes.
"""

# surface conditions
ps = 101540.0  # [Pa], surface pressure
p0 = 100000.0  # [Pa], reference pressure
Ts = 299.8  # [K], sea surface temperature

# constants
Lv = 2.5e6  # [J/kg]
cp_d = 1005.0  # [J/kg/K]
g = 9.81  # [m/s^2]
R_d = 287.0  # [J/kg/K]

# XXX: R_v and cp_v are not given in the RICO test definition on the
# the KNMI site I will use what I believe are standard values here
R_v= 461.51
cp_v= 1859.0

reference_points = OrderedDict(
    0.0 => ComponentArray(qt=16.0, theta_l=297.9),
    740.0 => ComponentArray(qt=13.8, theta_l=297.9),
    3260.0 => ComponentArray(qt=2.4, theta_l=317.0),
    4000.0 => ComponentArray(qt=2.4, theta_l=317.0),
)

"""Create a vertical profile that we can interpolate into later.
Integrating with the hydrostatic assumption.
"""
function integrate_profile(ref_profile, dz)
    z_max = 4e3

    z = 0.0
    p = ps

    # Cathy suggested using the liquid water potential temperature as the
    # temperature in the first model level
    T = ref_profile(0.0)[:theta_l]

    profile = []

    points = OrderedDict()

    n = 0
    while z < z_max
        qt, theta_l = ref_profile(z)

        # assume no liquid water
        ql = 0.0
        qv = qt

        qd = 1.0 - qt

        R_l = R_d * qd + R_v * qv
        c_l = cp_d * qd + cp_v * qv

        T = theta_l / ((p0 / p) ^ (R_l / c_l))
        # T = self.iteratively_find_temp(theta_l=theta_l, p=p, q_t=qt, q_l=ql, T_initial=T)

        rho = 1.0 / ((qd * R_d + qv * R_v) * T / p)  # + 1.0/(ql/rho_l), ql = 0.0

        points[z] = ComponentArray(rho=rho, p=p, T=T)
        #profile.append((z, rho, p, T))

        # integrate pressure
        z += dz
        p += -rho * g * dz

        n += 1
    end

    return points
end


function create_profile()
    ref_profile = Interpolate(collect(keys(reference_points)), collect(values(reference_points)))
    points_integrated = integrate_profile(ref_profile, 100);
    itp2 = Interpolate(collect(keys(points_integrated)), collect(values(points_integrated)))
    return ref_profile, itp2
end

ref_profile, itp2 = create_profile()


function get_rico_var(z, var_name)
    if var_name in [:theta_l, :qt]
        return ref_profile(z)[var_name]
    end
    if var_name in [:rho, :p, :T]
        return itp2(z)[var_name]
    end
    if var_name == :qv
        return ref_profile(z)[:qt]
    end
    throw("Can't compute $(var_name)")
end

export get_rico_var