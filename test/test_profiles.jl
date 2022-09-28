using Test
using CloudModels
using Unitful

profiles = [
    CloudModels.StandardIsentropicAtmosphere(),
    CloudModels.StandardIsothermalAtmosphere(),
    CloudModels.ConstantDensityAtmospere()
]


@testset "profiles" for prof in profiles
    ρ0 = calc_density(0.0u"m", prof)
    p0 = calc_pressure(0.0u"m", prof)

    @test ρ0 == CloudModels.rho0
    if prof == CloudModels.ConstantDensityAtmospere()
        @test calc_density(1u"km", prof) == ρ0
    else
        # density should decrease with height
        @test calc_density(1u"km", prof) < ρ0
    end

    @test p0 == CloudModels.p0
    # pressure should decrease with height
    @test calc_pressure(1u"km", prof) < p0
end