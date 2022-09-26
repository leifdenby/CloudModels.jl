using Test
using CloudModels

@testset "eqsn" begin
    F0 = ComponentArray(r=100, w=1.0, T=300, q_v=13.0e-3, q_l=0.0, q_r=0.0, q_i=0.0, p=100e3, q_pr=0.0)
    dFdz = zero(F0)
    p = (environment=get_rico_var, β=0.2)
    parcel_equations(dFdz, F0, 0.0, p)
end