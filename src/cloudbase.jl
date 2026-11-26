"""
Raised when no condensation level is found below the maximum search height.
"""
struct CloudbaseNotFoundException <: Exception
    z_max::Float64
end

Base.showerror(io::IO, e::CloudbaseNotFoundException) = print(
    io, "CloudbaseNotFoundException: parcel didn't saturate below z=$(e.z_max)m"
)

"""
    compute_LCL(environment, T0, qv0; dz=1.0, z_max=10e3)

Compute the lifting condensation level: the height at which a sub-saturated
parcel, lifted adiabatically from `z=0` with temperature `T0` and water vapour
specific concentration `qv0` (and no condensate), first becomes saturated.
The parcel is at the environment's pressure and its temperature decreases at
the dry-adiabatic lapse rate of the moist gas mixture, `-g/c_p` with
`c_p = q_d c_p,d + q_v c_p,v`. Returns `(z, T)` at the condensation level.

Ported from `compute_LCL` in pyclouds
(`pyclouds/integration/parcel_initiation.py`), including its fixed
integration step `dz` and its search ceiling of `z_max`.
"""
function compute_LCL(environment, T0, qv0; dz=1.0, z_max=10e3)
    qd = 1.0 - qv0
    cp_g = qd * cp_d + qv0 * cp_v

    # temperature decreases at dry adiabatic lapse rate
    dTdz = -g / cp_g

    z = 0.0
    T = T0
    while true
        p = environment(z, :p)
        qv_sat = calc_qv_sat(T, p)

        Sw = qv0 / qv_sat

        if Sw > 1.0
            return z, T
        end
        if z > z_max
            throw(CloudbaseNotFoundException(z_max))
        end

        z += dz
        T += dz * dTdz
    end
end
