"""
Idealised two-layer moist atmosphere: a sub-saturated well-mixed boundary
layer (isentropic, constant water vapour specific concentration) below `z_BL`,
a layer with a fixed lapse rate `dTdz_mid` and linearly decreasing relative
humidity up to the inversion at `z_INV`, and an isothermal layer above that.

Ported from `TwoLayerMoistIsentropicPBL` in pyclouds
(`pyclouds/reference/atmos/stratification_profiles.py`,
github.com/leifdenby/python-cloudmodel). The implementation follows pyclouds
step by step, including:

- the boundary-layer lapse rate is `-g/c_p` with `c_p` of the moist surface
  air,
- pressure and density are found by integrating the hydrostatic balance in
  pressure steps of `dp = -100 Pa` up to `p_min`, and linearly interpolated
  (clamped at the ends) in between,
- above `z_INV` relative humidity decreases at `dRHdz_2` (hard-coded to
  -0.4e-3 1/m in pyclouds) without a lower bound, so it becomes negative
  far enough above the inversion.

Use as `prof(z, var)` with `var` one of `:T`, `:p`, `:rho`, `:qv`, `:rh`.
"""
struct TwoLayerMoistIsentropicPBL
    z_BL::Float64
    z_INV::Float64
    T0::Float64
    p0::Float64
    RH0::Float64
    q_v0::Float64
    rho0::Float64
    dTdz_BL::Float64
    dTdz_mid::Float64
    dRHdz_1::Float64
    dRHdz_2::Float64
    # relative humidity at top of boundary layer, used for linear RH profile
    # above
    RH_BL_top::Float64
    # hydrostatically integrated profile
    z_pts::Vector{Float64}
    p_pts::Vector{Float64}
    rho_pts::Vector{Float64}
end

"""
Linear interpolation of `ys(xs)` at `x`, holding the end values constant
outside `xs` (like `numpy.interp`). `xs` must be increasing.
"""
function _interp_clamped(x, xs, ys)
    if x <= xs[1]
        return ys[1]
    elseif x >= xs[end]
        return ys[end]
    end
    i = searchsortedlast(xs, x)
    t = (x - xs[i]) / (xs[i + 1] - xs[i])
    return ys[i] + t * (ys[i + 1] - ys[i])
end

function _two_layer_temperature(z, T0, z_BL, z_INV, dTdz_BL, dTdz_mid)
    if z == 0.0
        return T0
    elseif z <= z_BL
        return T0 + z * dTdz_BL
    else
        T_BL_top = T0 + z_BL * dTdz_BL
        if z <= z_INV
            return T_BL_top + (z - z_BL) * dTdz_mid
        else
            return T_BL_top + (z_INV - z_BL) * dTdz_mid
        end
    end
end

function _two_layer_rel_humidity_above_BL(z, RH_BL_top, z_BL, z_INV, dRHdz_1, dRHdz_2)
    if z <= z_INV
        return RH_BL_top + (z - z_BL) * dRHdz_1
    else
        RH_top = RH_BL_top + (z_INV - z_BL) * dRHdz_1
        return RH_top + (z - z_INV) * dRHdz_2
    end
end

function TwoLayerMoistIsentropicPBL(;
    z_BL,
    RH0,
    T0,
    z_INV,
    dRHdz_1=-0.1e-3,
    p0=101325.0,
    dTdz_mid=-6.0e-3,
    dRHdz_2=-0.4e-3,
    p_min=500e2,
)
    qv_sat_0 = calc_qv_sat(T0, p0)
    q_v0 = RH0 * qv_sat_0
    q_d0 = 1.0 - q_v0

    cp = q_v0 * cp_v + q_d0 * cp_d
    R = q_v0 * R_v + q_d0 * R_d
    rho0 = p0 / (R * T0)
    dTdz_BL = -g / cp

    T_f(z) = _two_layer_temperature(z, T0, z_BL, z_INV, dTdz_BL, dTdz_mid)
    rho_f(T, p, qv) = 1.0 / (((1.0 - qv) * R_d + qv * R_v) * T / p)

    # numerical hydrostatic integration in pressure steps (to take into
    # account that the gas mixture changes with height)
    dp = -100.0
    p_pts = [Float64(p0)]
    rho_pts = [Float64(rho0)]
    z_pts = [0.0]

    while p_pts[end] > p_min
        p_new = p_pts[end] + dp
        # dp/dz = -rho * g
        dz = -dp / (rho_pts[end] * g)
        z_new = z_pts[end] + dz

        T = T_f(z_new)
        if z_new <= z_BL
            qv = q_v0
        else
            # NB: pyclouds evaluates the boundary-layer top relative humidity
            # from the pressure profile integrated so far (excluding the
            # current level, and clamped at its top), so we do the same here
            p_BL_top = _interp_clamped(z_BL, z_pts, p_pts)
            RH_BL_top = q_v0 / calc_qv_sat(T_f(z_BL), p_BL_top)
            RH = _two_layer_rel_humidity_above_BL(
                z_new, RH_BL_top, z_BL, z_INV, dRHdz_1, dRHdz_2
            )
            qv = RH * calc_qv_sat(T, p_new)
        end
        push!(p_pts, p_new)
        push!(z_pts, z_new)
        push!(rho_pts, rho_f(T, p_new, qv))
    end

    p_BL_top = z_BL == 0.0 ? p0 : _interp_clamped(z_BL, z_pts, p_pts)
    RH_BL_top = q_v0 / calc_qv_sat(T_f(z_BL), p_BL_top)

    return TwoLayerMoistIsentropicPBL(
        z_BL, z_INV, T0, p0, RH0, q_v0, rho0, dTdz_BL, dTdz_mid, dRHdz_1,
        dRHdz_2, RH_BL_top, z_pts, p_pts, rho_pts,
    )
end

function calc_temperature(z, prof::TwoLayerMoistIsentropicPBL)
    return _two_layer_temperature(
        z, prof.T0, prof.z_BL, prof.z_INV, prof.dTdz_BL, prof.dTdz_mid
    )
end

function calc_pressure(z, prof::TwoLayerMoistIsentropicPBL)
    return z == 0.0 ? prof.p0 : _interp_clamped(z, prof.z_pts, prof.p_pts)
end

function calc_density(z, prof::TwoLayerMoistIsentropicPBL)
    return z == 0.0 ? prof.rho0 : _interp_clamped(z, prof.z_pts, prof.rho_pts)
end

function calc_rel_humidity(z, prof::TwoLayerMoistIsentropicPBL)
    if z <= prof.z_BL
        T = calc_temperature(z, prof)
        p = calc_pressure(z, prof)
        return prof.q_v0 / calc_qv_sat(T, p)
    else
        return _two_layer_rel_humidity_above_BL(
            z, prof.RH_BL_top, prof.z_BL, prof.z_INV, prof.dRHdz_1, prof.dRHdz_2
        )
    end
end

function calc_qv(z, prof::TwoLayerMoistIsentropicPBL)
    T = calc_temperature(z, prof)
    p = calc_pressure(z, prof)
    return calc_rel_humidity(z, prof) * calc_qv_sat(T, p)
end

function (prof::TwoLayerMoistIsentropicPBL)(z, var_name::Symbol)
    if var_name == :T
        return calc_temperature(z, prof)
    elseif var_name == :p
        return calc_pressure(z, prof)
    elseif var_name == :rho
        return calc_density(z, prof)
    elseif var_name == :qv
        return calc_qv(z, prof)
    elseif var_name == :rh
        return calc_rel_humidity(z, prof)
    else
        throw("Not implemented: `$(var_name)`")
    end
end
