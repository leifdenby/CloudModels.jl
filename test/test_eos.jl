using Test
using CloudModels

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