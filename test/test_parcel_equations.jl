using Test
using CloudModels
using Unitful

@testset "eqns" begin
    F = ComponentArray(r=100u"m", w=1.0u"m/s", T=300u"K", q_v=13.0e-3, 
    q_l=0.0, q_r=0.0, q_i=0.0, p=100e3u"Pa", q_pr=0.0)
    dFdz = ComponentArray(Dict([ v => 0.0 * unit(F[v] * u"1/m") for v in keys(F)]))
    env_profile = CloudModels.StandardIsentropicAtmosphere()
    p = (environment=env_profile, β=0.2)
    CloudModels.parcel_equations!(dFdz, F, 0.0u"m", p)
end