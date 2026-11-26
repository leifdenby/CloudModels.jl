using Test
using CloudModels

@testset "lifting condensation level" begin
    # same two-layer setup as in test_two_layer_profile.jl
    prof = CloudModels.TwoLayerMoistIsentropicPBL(
        z_BL=600.0, RH0=0.7, T0=299.3, p0=101686.0, z_INV=2.6e3
    )
    T_surf = prof(0.0, :T)
    p_surf = prof(0.0, :p)

    function lcl(dT, RH)
        T0 = T_surf + dT
        qv0 = RH * CloudModels.calc_qv_sat(T0, p_surf)
        return CloudModels.compute_LCL(prof, T0, qv0)
    end

    @testset "matches pyclouds" begin
        # reference values from pyclouds' `compute_LCL`
        # (github.com/leifdenby/python-cloudmodel, commit 4223631) for a
        # parcel with surface temperature perturbation dT and relative
        # humidity RH
        cases = [
            # (dT, RH, z_LCL, T_LCL)
            (0.1, 0.6, 1048.0, 289.2891042926539),
            (1.0, 0.6, 1056.0, 290.11806219673855),
            (0.1, 0.8, 471.0, 294.87222560505063),
            (0.2, 0.7, 743.0, 292.345141438113),
        ]
        for (dT, RH, z_ref, T_ref) in cases
            z, T = lcl(dT, RH)
            @test z == z_ref
            @test T ≈ T_ref rtol = 1e-10
        end
    end

    @testset "behaviour" begin
        # same checks as pyclouds' tests/test_LCL.py
        z1, _ = lcl(0.1, 0.6)
        @test 0.0 < z1 < 1300.0
        # higher initial temperature raises the condensation level
        z2, _ = lcl(1.0, 0.6)
        @test z1 < z2
        # higher initial relative humidity lowers the condensation level
        z3, _ = lcl(0.1, 0.8)
        @test z3 < z1

        # parcel is saturated at the LCL but not one step below it
        T0 = T_surf + 0.1
        qv0 = 0.6 * CloudModels.calc_qv_sat(T0, p_surf)
        z, T = CloudModels.compute_LCL(prof, T0, qv0)
        @test qv0 > CloudModels.calc_qv_sat(T, prof(z, :p))
        cp_g = (1.0 - qv0) * CloudModels.cp_d + qv0 * CloudModels.cp_v
        T_below = T + CloudModels.g / cp_g
        @test qv0 < CloudModels.calc_qv_sat(T_below, prof(z - 1.0, :p))
        # and has followed the dry adiabat of the moist gas mixture
        @test T ≈ T0 - CloudModels.g / cp_g * z rtol = 1e-8

        # a saturated parcel condenses straight away
        z_sat, T_sat = lcl(0.0, 1.01)
        @test z_sat == 0.0
        @test T_sat == T_surf

        # no condensation level below the search ceiling (z=1048m for this
        # parcel, see above)
        @test_throws CloudModels.CloudbaseNotFoundException CloudModels.compute_LCL(
            prof, T0, qv0; z_max=500.0
        )
    end
end
