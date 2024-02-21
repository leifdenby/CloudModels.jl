using Test
using CloudModels
using ComponentArrays
using BenchmarkTools


@testset "heat-capacities" begin
    qv0 = 0.0  # [kg/kg]
    qr0 = 0.0
    ql0 = 0.0
    qi0 = 0.0
    qd0 = 1.0 - qv0 - ql0 - qr0 - qi0
    
    F0 = ComponentArray(q_v=qv0, q_l=ql0, q_r=qr0, q_i=qi0)
    cp_m0 = CloudModels.calc_cp_m(F0)
    cv_m0 = CloudModels.calc_cv_m(F0)

    F1 = ComponentArray(q_v=qv0+1.0e-3, q_l=ql0, q_r=qr0, q_i=qi0)
    cp_m1 = CloudModels.calc_cp_m(F1)
    cv_m1 = CloudModels.calc_cv_m(F1)

    F2 = ComponentArray(q_v=qv0, q_l=ql0+1.0e-3, q_r=qr0, q_i=qi0)
    cp_m2 = CloudModels.calc_cp_m(F2)
    cv_m2 = CloudModels.calc_cv_m(F2)

    F3 = ComponentArray(q_v=qv0, q_l=ql0, q_r=qr0+1.0e-3, q_i=qi0)
    cp_m3 = CloudModels.calc_cp_m(F3)
    cv_m3 = CloudModels.calc_cv_m(F3)
    
    # water vapour has higher heat capacity than dry air
    @test cp_m0 < cp_m1
    @test cv_m0 < cv_m1
    # liquid water has higher heat capacity than dry air and water vapour
    @test cp_m1 < cp_m2
    @test cv_m1 < cv_m2
    # liquid water and rain have same heat capacity
    @test cp_m2 == cp_m3
    @test cv_m2 == cv_m3
    
    trial = @benchmark CloudModels.calc_cp_m($F0)
    @test trial.allocs == 0
    trial = @benchmark CloudModels.calc_cv_m($F0)
    @test trial.allocs == 0
end