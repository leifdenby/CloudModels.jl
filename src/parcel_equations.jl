function parcel_equations!(dFdt, F, params, t)
    dFdt .= 0.0
    dzdt = F.w
    plume_equations!(dFdt, F, params, F.z * u"m")
    # dF/dz is stored in `dFdt` vector, but need dF/dt: dF/dt = dF/dz * dz/dt
    dFdt .*= dzdt
    dFdt[:z] = dzdt
end