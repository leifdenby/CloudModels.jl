"""
Constants:
    cp_d: heat capacity of dry air at constant pressure
    cp_v: heat capacity of liquid water at constant pressure
    L_v: latent heat of vapourisation (vapour -> liquid)
    L_s: latent heat sublimation (vapour -> solid)
"""

R_d= 287.05
R_v= 461.51
L_v= 2.5008e6
L_s= 2.8345e6
cp_d= 1005.46
cv_d= 717.60  # J/kg/K
cp_v= 1859.0
cv_v= 1402.5  # J/kg/K
cp_l= 4183.0
cv_l= 4183.0  # same as cp as liquid is assumed incompressible
cp_i= 2103.0  # taken from ATHAM
rho_l= 1000.0
rho_i= 500.0
g= 9.80665

export R_d, R_v, L_v, L_s, cp_d, cv_d, cp_v, cv_v, rho_l, rho_i, g