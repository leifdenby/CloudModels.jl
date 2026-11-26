using Test
using CloudModels
using ComponentArrays

# `calc_mixture_density(p, T, qd, qv, ql, qr, qi)` treats rain like cloud
# water (density rho_l) and ice with density rho_i (thesis eqn 3.5). These
# tests check that the callers pass rain and ice in that order, i.e. that rain
# is never treated as ice. Ice is not implemented in the model, so q_i = 0 in
# all of them; with rain swapped into the ice slot the density would be
# computed with rho_i instead of rho_l.

@testset "rain is liquid in the mixture density" begin
    p = 95.0e3u"Pa"
    T = 290.0u"K"
    q_v = 14.0e-3
    q_l = 1.0e-3
    q_r = 1.0e-3
    q_i = 0.0
    q_d = 1.0 - q_v - q_l - q_r - q_i

    rho_c = CloudModels.calc_mixture_density(p, T, q_d, q_v, q_l, q_r, q_i)
    # sanity check: with rain as ice the density differs (by about 1e-6
    # relative), which is what the tests below resolve
    rho_c_rain_as_ice = CloudModels.calc_mixture_density(p, T, q_d, q_v, q_l, q_i, q_r)
    @test !isapprox(rho_c, rho_c_rain_as_ice; rtol=1.0e-10)

    @testset "calc_dw_dz and calc_dr_dz only see total liquid" begin
        # moving condensate between cloud water and rain must not change the
        # density-dependent dynamics
        w, r = 2.0u"m/s", 300.0u"m"
        rho_e = CloudModels.calc_mixture_density(p, T - 1.0u"K", 1.0 - q_v, q_v, 0.0, 0.0, 0.0)
        mu, gamma, C_D = 0.2 / r, 0.5, 0.506
        dw_a = CloudModels.calc_dw_dz(p, w, r, T, q_d, q_v, q_l, q_r, q_i, rho_e, mu, gamma, C_D)
        dw_b = CloudModels.calc_dw_dz(p, w, r, T, q_d, q_v, q_l + q_r, 0.0, q_i, rho_e, mu, gamma, C_D)
        @test isapprox(dw_a, dw_b; rtol=1.0e-12)

        dqv, dql, dqi, dT, dw = -1.0e-6u"1/m", 1.0e-6u"1/m", 0.0u"1/m", -6.0e-3u"K/m", 1.0e-3u"1/s"
        dr_a = CloudModels.calc_dr_dz(p, w, r, T, q_d, q_v, q_l, q_r, q_i, dqv, dql, dqi, dT, dw, mu)
        dr_b = CloudModels.calc_dr_dz(p, w, r, T, q_d, q_v, q_l + q_r, 0.0, q_i, dqv, dql, dqi, dT, dw, mu)
        @test isapprox(dr_a, dr_b; rtol=1.0e-12)
    end

    @testset "microphysics uses the mixture density with rain as liquid" begin
        F = ComponentArray(T=T, q_v=q_v, q_l=q_l, q_r=q_r, q_i=q_i, p=p)
        dFdt = zero(F)
        CloudModels.dFdt_microphysics!(dFdt, F, nothing, 0.0u"s")

        dql_dt = CloudModels._dql_dt__cond_evap(q_v, q_l, rho_c, p, T)
        dqr_dt_condevap = CloudModels._dqr_dt__cond_evap(q_v, q_r, rho_c, p, T)
        @test isapprox(dFdt.q_v, -dql_dt - dqr_dt_condevap; rtol=1.0e-12)
    end

    @testset "plume rain-out uses the mixture density with rain as liquid" begin
        env_profile = CloudModels.ProfileRICO.RICO_profile()
        z = 1000.0u"m"
        w, r = 2.0u"m/s", 300.0u"m"
        F = ComponentArray(r=r, w=w, T=T, q_v=q_v, q_l=q_l, q_r=q_r, q_i=q_i, q_pr=0.0, p=0.0)
        params = (environment=env_profile, β=0.2)
        dFdz = zero(F)
        CloudModels.plume_equations!(dFdz, F, params, z)

        # plume_equations! uses the environment pressure in the cloud
        p_env = env_profile(z, :p)
        rho_c_env = CloudModels.calc_mixture_density(p_env, T, q_d, q_v, q_l, q_r, q_i)
        # default N0 and l_pr from plume_equations!
        dqr_dz__rainout = CloudModels.calc_dqr_dz__rainout(rho_c_env, q_r, w, 1.0e7u"m^-4", 100.0u"m")
        @test isapprox(dFdz.q_pr, dqr_dz__rainout; rtol=1.0e-12)
    end
end
