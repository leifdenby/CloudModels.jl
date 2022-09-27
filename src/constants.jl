"""
Constants:
    cp_d: heat capacity of dry air at constant pressure
    cp_v: heat capacity of liquid water at constant pressure
    L_v: latent heat of vapourisation (vapour -> liquid)
    L_s: latent heat sublimation (vapour -> solid)
"""

L_v = 2.5008e6 * u"J/kg"
L_s = 2.8345e6 * u"J/kg"
R_d = 287.05 * u"J/K/kg"
R_v = 461.51 * u"J/K/kg"
cp_d = 1005.46 * u"J/K/kg"
cv_d = 717.60 * u"J/K/kg"
cp_v = 1859.0 * u"J/K/kg"
cv_v = 1402.5 * u"J/K/kg"
cp_l = 4183.0 * u"J/K/kg"
cv_l = cp_l # same as cp as liquid is assumed incompressible
cp_i = 2103.0 * u"J/K/kg"  # taken from ATHAM
rho_l = 1000.0 * u"kg/m^3"
rho_i = 500.0 * u"kg/m^3"
g = 9.80665 * u"m/s^2"

export R_d, R_v, L_v, L_s, cp_d, cv_d, cp_v, cv_v, rho_l, rho_i, g