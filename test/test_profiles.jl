using Test
using CloudModels
using Unitful

profile_types = [
    CloudModels.StandardIsentropicAtmosphere,
    CloudModels.StandardIsothermalAtmosphere,
    CloudModels.ConstantDensityAtmospere,
    CloudModels.ProfileRICO.RICO_profile,
]


@testset "profiles $(p_type)" for p_type in profile_types
    prof = p_type()
    ρ0 = prof(0.0u"m", :rho)
    p0 = prof(0.0u"m", :p)

    if prof == CloudModels.ConstantDensityAtmospere()
        @test prof(1u"km", :rho) == ρ0
    elseif prof == CloudModels.ProfileRICO.RICO_profile()
        # pass
    else
        # density should decrease with height
        @test prof(1u"km", :rho) < ρ0
    end

    # TODO expose reference values in RICO profile
    if p_type != CloudModels.ProfileRICO.RICO_profile
        @test p0 == prof.p0
        @test ρ0 == CloudModels.rho0
    end

    # pressure should decrease with height
    @test prof(1u"km", :p) < p0
end

using Plots

z_ = collect(0:10:3500)u"m"
prof = CloudModels.ProfileRICO.RICO_profile()
plot(ustrip.(prof.(z_, :rho)), ustrip.(z_))

prof(0.0u"m", :rho)
prof.itp2