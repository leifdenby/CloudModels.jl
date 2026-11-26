using Test
using CloudModels

# setup used in the pyclouds notebook "Entrainment contribution terms - moist plume"
two_layer_setup = (z_BL=600.0, RH0=0.7, T0=299.3, p0=101686.0, z_INV=2.6e3)

@testset "two-layer moist isentropic PBL profile" begin
    prof = CloudModels.TwoLayerMoistIsentropicPBL(; two_layer_setup...)

    @testset "matches pyclouds" begin
        # reference values from pyclouds' `TwoLayerMoistIsentropicPBL`
        # (github.com/leifdenby/python-cloudmodel, commit 4223631) with
        # `constants=default_constants` and the setup above
        z = [0.0, 100.0, 300.0, 600.0, 601.0, 1000.0, 2000.0, 2600.0, 3000.0, 4000.0]
        T_ref = [
            299.3, 298.33688444562443, 296.4106533368732, 293.52130667374644,
            293.51530667374647, 291.12130667374646, 285.12130667374646,
            281.52130667374644, 281.52130667374644, 281.52130667374644,
        ]
        p_ref = [
            101686.0, 100540.01166531892, 98275.76975077731, 94948.13711680176,
            94937.18229649315, 90645.30224039452, 80566.63261211374,
            74965.90757295533, 71424.50124004515, 63274.2583068766,
        ]
        rho_ref = [
            1.172919786497075, 1.163444913075372, 1.144633632879301,
            1.1167388345156328, 1.1166398975745677, 1.0762875827783978,
            0.9792194143580478, 0.9238405044547741, 0.8808781147007567,
            0.7820945049414023,
        ]
        rh_ref = [
            0.7, 0.7336391723082781, 0.8067192312683142, 0.9328095021604256,
            0.9327095021604256, 0.8928095021604255, 0.7928095021604256,
            0.7328095021604255, 0.5728095021604255, 0.17280950216042545,
        ]
        qv_ref = [
            0.014951305049351979, 0.014951305049351979, 0.014951305049351979,
            0.014951305049351979, 0.014945797689987779, 0.012876015995879002,
            0.008690985281441433, 0.006763486789394193, 0.005550436791150945,
            0.0018916217108976415,
        ]

        @test prof.(z, :T) ≈ T_ref rtol = 1e-10
        @test prof.(z, :p) ≈ p_ref rtol = 1e-10
        @test prof.(z, :rho) ≈ rho_ref rtol = 1e-10
        @test prof.(z, :rh) ≈ rh_ref rtol = 1e-10
        @test prof.(z, :qv) ≈ qv_ref rtol = 1e-10
        # number of levels in the hydrostatic integration (pyclouds: 518)
        @test length(prof.z_pts) == 518
    end

    @testset "structure" begin
        # well-mixed boundary layer: constant water vapour, dry-adiabatic
        # lapse rate for the moist air
        @test prof(0.0, :rh) ≈ two_layer_setup.RH0
        @test prof(300.0, :qv) ≈ prof(0.0, :qv)
        @test prof(600.0, :qv) ≈ prof(0.0, :qv)
        q_v0 = prof(0.0, :qv)
        cp = (1.0 - q_v0) * CloudModels.cp_d + q_v0 * CloudModels.cp_v
        @test (prof(500.0, :T) - prof(100.0, :T)) / 400.0 ≈ -CloudModels.g / cp

        # above the boundary layer: prescribed lapse rate and RH gradient up
        # to the inversion, isothermal above
        @test (prof(2000.0, :T) - prof(1000.0, :T)) / 1000.0 ≈ -6.0e-3
        @test (prof(2000.0, :rh) - prof(1000.0, :rh)) / 1000.0 ≈ -0.1e-3
        @test prof(3500.0, :T) == prof(2600.0, :T)
        @test (prof(3500.0, :rh) - prof(3000.0, :rh)) / 500.0 ≈ -0.4e-3

        # hydrostatic balance, dp/dz = -rho g
        dz = 1.0
        for z in [200.0, 1500.0, 3000.0]
            dpdz = (prof(z + dz, :p) - prof(z - dz, :p)) / (2dz)
            @test dpdz ≈ -prof(z, :rho) * CloudModels.g rtol = 1e-3
        end

        # pressure and density decrease with height
        z_ = collect(0.0:100.0:5000.0)
        @test issorted(prof.(z_, :p), rev=true)
        @test issorted(prof.(z_, :rho), rev=true)
    end
end
