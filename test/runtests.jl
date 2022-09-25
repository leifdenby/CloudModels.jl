using Test
using CloudModels
using ComponentArrays

@testset "EoS" begin
    p0 = 101.325e3  # [hPa]
    T0 = 273.15 + 15.0  # [K]
    qv0 = 0.0  # [kg/kg]
    qr0 = 0.0
    ql0 = 0.0
    qi0 = 0.0
    rho0 = CloudModels.calc_mixture_density(p0, T0, qv0, ql0, qr0, qi0)
    @test isapprox(rho0, 1.225; atol=0.001)

    # more water vapour makes the mixture lighter
    rho1 = CloudModels.calc_mixture_density(p0, T0, qv0 + 1.0e-3, ql0, qr0, qi0)
    @test rho1 < rho0

    # more liquid makes the mixture heavier
    rho2 = CloudModels.calc_mixture_density(p0, T0, qv0, ql0 + 1.0e-3, qr0, qi0)
    rho3 = CloudModels.calc_mixture_density(p0, T0, qv0, ql0, qr0 + 1.0e-3, qi0)
    @test rho2 > rho0
    @test rho3 > rho0
end


@testset "microphysics" begin
    p0 = 101.0e3
    T0 = 300.0
    qv0_sat = CloudModels.calc_qv_sat(T0, p0)

    @testset "super-saturated" begin
        # start of super-saturated
        F0 = ComponentArray(T=T0, q_v=1.1*qv0_sat, q_l=0.0, q_r=0.0, q_i=0.0, p=p0)
        dFdt0 = CloudModels.dFdt_microphysics(F0, 0)
        # because we are super-saturated we expect condensation to occour which
        # should release heat and increase the amount of cloud condensate
        @test dFdt0.T > 0
        @test dFdt0.q_v < 0
        @test dFdt0.q_l > 0
    end

    @testset "super-saturated /w condensate" begin
        # start of super-saturated
        F0 = ComponentArray(T=T0, q_v=1.1*qv0_sat, q_l=0.1, q_r=0.0, q_i=0.0, p=p0)
        dFdt0 = CloudModels.dFdt_microphysics(F0, 0)
        # because we are super-saturated we expect condensation to occour which
        # should release heat and increase the amount of cloud condensate
        @test dFdt0.T > 0
        @test dFdt0.q_v < 0
        @test dFdt0.q_l > 0
        @test dFdt0.q_r > 0
    end

    @testset "sub-saturated" begin
        # start of sub-saturated with cloud-condensate present
        F1 = ComponentArray(T=T0, q_v=qv0_sat*0.9, q_l=1.0e-3, q_r=0.0, q_i=0.0, p=p0)
        dFdt1 = CloudModels.dFdt_microphysics(F1, 0)
        # because we are sub-saturated we expect condensation top occour which
        # should take heat and decrease the amount of cloud condensate
        @test dFdt1.T < 0
        @test dFdt1.q_v > 0
        @test dFdt1.q_l < 0
    end
end