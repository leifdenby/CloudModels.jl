using Test
using CloudModels
using ComponentArrays

@testset "entrainment switches" begin
    env = CloudModels.TwoLayerMoistIsentropicPBL(
        z_BL=600.0, RH0=0.7, T0=299.3, p0=101686.0, z_INV=2.6e3
    )
    z = 1000.0
    T_e = env(z, :T)
    p_e = env(z, :p)
    qv_e = env(z, :qv)

    # slightly supersaturated cloudy parcel, warmer than the environment
    T_c = T_e + 0.5
    F = ComponentArray(
        r=300.0, w=1.0, T=T_c, q_v=1.001 * CloudModels.calc_qv_sat(T_c, p_e),
        q_l=1.0e-4, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0
    )
    β = 0.2
    mu = β / F.r

    function eval_dFdz(; kwargs...)
        params = (environment=env, β=β, kwargs...)
        dFdz = zero(F)
        CloudModels.plume_equations!(dFdz, copy(F), params, z)
        return dFdz
    end

    dFdz_default = eval_dFdz()

    @testset "defaults unchanged" begin
        dFdz_explicit = eval_dFdz(
            entrain_hydrometeors=true,
            entrain_moist_static_energy=true,
            entrain_liquid_static_energy=false,
            temperature_dependent_latent_heats=true,
        )
        @test dFdz_explicit == dFdz_default
    end

    @testset "entrain_hydrometeors" begin
        dFdz = eval_dFdz(entrain_hydrometeors=false)
        # only the entrainment tendencies of the water species are removed
        @test dFdz[:q_v] ≈ dFdz_default[:q_v] - mu * (qv_e - F.q_v)
        @test dFdz[:q_l] ≈ dFdz_default[:q_l] - mu * (0.0 - F.q_l)
        @test dFdz[:w] == dFdz_default[:w]
    end

    @testset "no static energy entrainment" begin
        dFdz = eval_dFdz(entrain_moist_static_energy=false)
        # dry-adiabatic cooling plus latent heating from the microphysics
        F_m = copy(F)
        F_m.p = p_e
        dFdt_micro = zero(F)
        CloudModels.dFdt_microphysics!(dFdt_micro, F_m, nothing, 0.0)
        @test dFdz[:T] ≈ -CloudModels.g / CloudModels.cp_d + dFdt_micro[:T] / F.w
        @test dFdz[:T] != dFdz_default[:T]
    end

    @testset "liquid static energy" begin
        dFdz = eval_dFdz(
            entrain_moist_static_energy=false, entrain_liquid_static_energy=true
        )
        @test dFdz[:T] != dFdz_default[:T]

        # entraining both moist and liquid static energy isn't allowed
        @test_throws "Can only entrain" eval_dFdz(entrain_liquid_static_energy=true)

        # without entrainment or phase changes: dry adiabat of the mixture
        # with all water as vapour
        q_d = 1.0 - F.q_v - F.q_l
        args = (F.r, T_c, q_d, F.q_v, F.q_l, 0.0)
        dTdz = CloudModels.calc_dT_dz(
            args..., 0.0, 0.0, 0.0, T_e, qv_e, 0.0;
            entrain_moist_static_energy=false, entrain_liquid_static_energy=true,
        )
        c_cm_p = CloudModels.cp_d * q_d + CloudModels.cp_v * (F.q_v + F.q_l)
        @test dTdz ≈ -CloudModels.g / c_cm_p

        # condensation heats with the temperature-corrected L_v
        dql = 1.0e-6
        dTdz_cond = CloudModels.calc_dT_dz(
            args..., dql, 0.0, -dql, T_e, qv_e, 0.0;
            entrain_moist_static_energy=false, entrain_liquid_static_energy=true,
        )
        L_v_T = CloudModels.L_v + (CloudModels.cp_v - CloudModels.cp_l) * (T_c - CloudModels.T00)
        @test dTdz_cond - dTdz ≈ L_v_T / c_cm_p * dql
    end

    @testset "temperature_dependent_latent_heats" begin
        dFdz = eval_dFdz(temperature_dependent_latent_heats=false)
        @test dFdz[:T] != dFdz_default[:T]

        q_d = 1.0 - F.q_v - F.q_l
        dqv = -1.0e-6
        function dTdz_(T; T_dep)
            CloudModels.calc_dT_dz(
                F.r, T, q_d, F.q_v, F.q_l, 0.0, -dqv, 0.0, dqv, T_e, qv_e, 0.0;
                temperature_dependent_latent_heats=T_dep,
            )
        end
        # identical at the reference temperature of the latent heats
        @test dTdz_(CloudModels.T00; T_dep=true) ≈ dTdz_(CloudModels.T00; T_dep=false)
        # otherwise differ by the change in L_v
        c_cm = (
            CloudModels.cp_d * q_d + CloudModels.cp_v * F.q_v
            + CloudModels.cp_l * F.q_l
        )
        dL_v = (CloudModels.cp_v - CloudModels.cp_l) * (T_c - CloudModels.T00)
        @test dTdz_(T_c; T_dep=true) - dTdz_(T_c; T_dep=false) ≈ -dL_v * dqv / c_cm
    end
end
