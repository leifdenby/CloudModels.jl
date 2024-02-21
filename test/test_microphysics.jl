using Test
using CloudModels
using ComponentArrays
using BenchmarkTools
# using Unitful


@testset "microphysics" begin
    p0 = 101.0e3u"Pa"
    T0 = 300.0u"K"
    qv0_sat = CloudModels.calc_qv_sat(T0, p0)

    @testset "processes" begin
        @testset "cloud cond-evap Sw=$(Sw) ql=$(ql)" for Sw in [0.9, 1.0, 1.1], ql in [0.0, 0.1e-3]
            qv = qv0_sat*Sw
            qr = 0.0
            qi = 0.0
            qd = 1.0 - qv - ql - qr - qi
            rho = CloudModels.calc_mixture_density(p0, T0, qd, qv, ql, qr, qi)
            dql_dt = CloudModels._dql_dt__cond_evap(qv, ql, rho, p0, T0)
            if Sw > 1.0
                @test dql_dt > 0.0u"1/s"
            elseif Sw < 1.0 && ql > 0.0
                @test dql_dt < 0.0u"1/s"
            else
                @test dql_dt == 0.0u"1/s"
            end
            trial = @benchmark CloudModels._dql_dt__cond_evap($qv, $ql, $rho, $p0, $T0)
            @test trial.allocs == 0
        end

        @testset "rain cond-evap qr=$(qr) Sw=$(Sw)" for Sw in [0.9, 1.0, 1.1], qr in [0.0, 0.1e-3]
            qv = qv0_sat*Sw
            ql = 0.0
            qi = 0.0
            qd = 1.0 - qv - ql - qr - qi
            rho = CloudModels.calc_mixture_density(p0, T0, qd, qv, ql, qr, qi)
            dqr_dt = CloudModels._dqr_dt__cond_evap(qv, qr, rho, p0, T0)
            if Sw > 1.0 && qr > 0.0
                @test dqr_dt > 0.0u"1/s"
            elseif Sw < 1.0 && qr > 0.0
                @test dqr_dt < 0.0u"1/s"
            else
                @test dqr_dt == 0.0u"1/s"
            end
            trial = @benchmark CloudModels._dqr_dt__cond_evap($qv, $qr, $rho, $p0, $T0)
            @test trial.allocs == 0
        end

        @testset "accretion qr=$(qr) ql=$(ql)" for qr in [0.0, 0.1e-3], ql in [0.0, 0.1e-3]
            qv = qv0_sat
            qi = 0.0
            qd = 1.0 - qv - ql - qr - qi
            rho_g = CloudModels.calc_mixture_density(p0, T0, qd, 0.0, 0.0, qr, 0.0)
            qg = qv + qd
            dqr_dt = CloudModels._dqr_dt__accretion(ql, qg, qr, rho_g)
            if qr == 0.0 || ql == 0.0
                @test dqr_dt == 0.0u"1/s"
            else
                @test dqr_dt > 0.0u"1/s"
            end
            trial = @benchmark CloudModels._dqr_dt__accretion($ql, $qg, $qr, $rho_g)
            @test trial.allocs == 0
        end
    end

    @testset "super-saturated" begin
        # start of super-saturated
        F = ComponentArray(T=T0, q_v=1.1*qv0_sat, q_l=0.0, q_r=0.0, q_i=0.0, p=p0)
        dFdt = ComponentArray(Dict([ v => 0.0 * unit(F[v] * u"1/s") for v in keys(F)]))
        CloudModels.dFdt_microphysics!(dFdt, F, nothing, 0)
        # because we are super-saturated we expect condensation to occour which
        # should release heat and increase the amount of cloud condensate
        @test dFdt.T > 0u"K/s"
        @test dFdt.q_v < 0u"1/s"
        @test dFdt.q_l > 0u"1/s"
        
        trial = @benchmark CloudModels.dFdt_microphysics!($dFdt, $F, nothing, 0)
        @test trial.allocs == 0
    end

    @testset "super-saturated /w condensate" begin
        # start of super-saturated
        F = ComponentArray(T=T0, q_v=1.1*qv0_sat, q_l=0.1, q_r=0.0, q_i=0.0, p=p0)
        dFdt = ComponentArray(Dict([ v => 0.0 * unit(F[v] * u"1/s") for v in keys(F)]))
        CloudModels.dFdt_microphysics!(dFdt, F, nothing, 0)
        # because we are super-saturated we expect condensation to occour which
        # should release heat and increase the amount of cloud condensate
        @test dFdt.T > 0u"K/s"
        @test dFdt.q_v < 0u"1/s"
        @test dFdt.q_l > 0u"1/s"
        @test dFdt.q_r > 0u"1/s"

        trial = @benchmark CloudModels.dFdt_microphysics!($dFdt, $F, nothing, 0)
        @test trial.allocs == 0
    end

    @testset "sub-saturated" begin
        # start of sub-saturated with cloud-condensate present
        F = ComponentArray(T=T0, q_v=qv0_sat*0.9, q_l=1.0e-3, q_r=0.0, q_i=0.0, p=p0)
        dFdt = ComponentArray(Dict([ v => 0.0 * unit(F[v] * u"1/s") for v in keys(F)]))
        CloudModels.dFdt_microphysics!(dFdt, F, nothing, 0)
        # because we are sub-saturated we expect condensation top occour which
        # should take heat and decrease the amount of cloud condensate
        @test dFdt.T < 0u"K/s"
        @test dFdt.q_v > 0u"1/s"
        @test dFdt.q_l < 0u"1/s"

        trial = @benchmark CloudModels.dFdt_microphysics!($dFdt, $F, nothing, 0)
        @test trial.allocs == 0
    end
end


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