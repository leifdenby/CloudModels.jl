using Test
using CloudModels
using Unitful


@testset "microphysics" begin
    p0 = 101.0e3u"Pa"
    T0 = 300.0u"K"
    qv0_sat = CloudModels.calc_qv_sat(T0, p0)

    @testset "processes" begin
        @test_skip @testset "cloud cond-evap Sw=$(Sw) ql=$(ql)" for Sw in [0.9, 1.0, 1.1], ql in [0.0, 0.1e-3]
            qv = qv0_sat*Sw
            qr = 0.0
            qi = 0.0
            qd = 1.0 - qv - ql - qr - qi
            rho = CloudModels.calc_mixture_density(p0, T0, qd, qv, ql, qr, qi)
            dql_dt = CloudModels._dql_dt__cond_evap(qv, qr, rho, T0, p0)
            if Sw > 1.0
                @test dql_dt > 0.0u"1/s"
            elseif Sw < 1.0 && ql > 0.0
                @test dql_dt < 0.0u"1/s"
            else
                @test dql_dt == 0.0u"1/s"
            end
        end

        @testset "rain cond-evap qr=$(qr) Sw=$(Sw)" for Sw in [0.9, 1.0, 1.1], qr in [0.0, 0.1e-3]
            qv = qv0_sat*Sw
            ql = 0.0
            qi = 0.0
            qd = 1.0 - qv - ql - qr - qi
            rho = CloudModels.calc_mixture_density(p0, T0, qd, qv, ql, qr, qi)
            dqr_dt = CloudModels._dqr_dt__cond_evap(qv, qr, rho, T0, p0)
            if Sw > 1.0 && qr > 0.0
                @test dqr_dt > 0.0u"1/s"
            elseif Sw < 1.0 && qr > 0.0
                @test dqr_dt < 0.0u"1/s"
            else
                @test dqr_dt == 0.0u"1/s"
            end
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
        end
    end

    @testset "super-saturated" begin
        # start of super-saturated
        F0 = ComponentArray(T=T0, q_v=1.1*qv0_sat, q_l=0.0, q_r=0.0, q_i=0.0, p=p0)
        dFdt0 = CloudModels.dFdt_microphysics(F0, 0)
        # because we are super-saturated we expect condensation to occour which
        # should release heat and increase the amount of cloud condensate
        @test dFdt0.T > 0u"K/s"
        @test dFdt0.q_v < 0u"1/s"
        @test dFdt0.q_l > 0u"1/s"
    end

    @testset "super-saturated /w condensate" begin
        # start of super-saturated
        F0 = ComponentArray(T=T0, q_v=1.1*qv0_sat, q_l=0.1, q_r=0.0, q_i=0.0, p=p0)
        dFdt0 = CloudModels.dFdt_microphysics(F0, 0)
        # because we are super-saturated we expect condensation to occour which
        # should release heat and increase the amount of cloud condensate
        @test dFdt0.T > 0u"K/s"
        @test dFdt0.q_v < 0u"1/s"
        @test dFdt0.q_l > 0u"1/s"
        @test dFdt0.q_r > 0u"1/s"
    end

    @test_skip @testset "sub-saturated" begin
        # start of sub-saturated with cloud-condensate present
        F1 = ComponentArray(T=T0, q_v=qv0_sat*0.9, q_l=1.0e-3, q_r=0.0, q_i=0.0, p=p0)
        dFdt1 = CloudModels.dFdt_microphysics(F1, 0)
        # because we are sub-saturated we expect condensation top occour which
        # should take heat and decrease the amount of cloud condensate
        @test dFdt1.T < 0u"K/s"
        @test dFdt1.q_v > 0u"1/s"
        @test dFdt1.q_l < 0u"1/s"
    end
end