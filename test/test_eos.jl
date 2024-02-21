using Test
using CloudModels
using BenchmarkTools
# using Unitful

@testset "EoS" begin
    p0 = 101.325e3u"Pa"
    T0 = (273.15 + 15.0)u"K"
    qv0 = 0.0  # [kg/kg]
    qr0 = 0.0
    ql0 = 0.0
    qi0 = 0.0
    qd0 = 1.0 - qv0 - ql0 - qr0 - qi0
    rho0 = CloudModels.calc_mixture_density(p0, T0, qd0, qv0, ql0, qr0, qi0)
    @test isapprox(rho0, 1.225u"kg/m^3"; atol=0.001u"kg/m^3")

    dq = 1.0e-3

    # more water vapour makes the mixture lighter
    rho1 = CloudModels.calc_mixture_density(p0, T0, qd0-dq, qv0+dq, ql0, qr0, qi0)
    @test rho1 < rho0

    # more liquid makes the mixture heavier
    rho2 = CloudModels.calc_mixture_density(p0, T0, qd0-dq, qv0, ql0+dq, qr0, qi0)
    rho3 = CloudModels.calc_mixture_density(p0, T0, qd0-dq, qv0, ql0, qr0+dq, qi0)
    @test rho2 > rho0
    @test rho3 > rho0

    trial = @benchmark CloudModels.calc_mixture_density($p0, $T0, $qd0, $qv0, $ql0, $qr0, $qi0)
    @test trial.allocs == 0
end