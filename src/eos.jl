

function calc_mixture_density(p, T, qd, qv, ql, qr, qi)
    """
    Compute the mixture density from the import full equation of state

    Constants:
        R_d: specific gas constant of dry air
        R_v: specific gas constant of water vapour
        rho_l: density of liquid water
        rho_i: density of ice
    """
    rho_inv = (
        (qd * R_d + qv * R_v) * T / p + (ql + qr) / rho_l + qi / rho_i
    )

    return 1.0 / rho_inv
end
