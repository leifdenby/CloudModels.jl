using Test
using CloudModels
using ComponentArrays
using OrdinaryDiffEq

@testset "prescribed entrainment rate μ(z)" begin
    env = CloudModels.TwoLayerMoistIsentropicPBL(
        z_BL=600.0, RH0=0.7, T0=299.3, p0=101686.0, z_INV=2.6e3
    )

    @testset "calc_entrainment_rate" begin
        # default: β/r with β=0.2
        @test CloudModels.calc_entrainment_rate((;), 400.0, 1000.0) == 0.2 / 400.0
        @test CloudModels.calc_entrainment_rate((β=0.1,), 400.0, 1000.0) == 0.1 / 400.0
        # prescribed constant or profile, independent of radius
        @test CloudModels.calc_entrainment_rate((μ=1.0e-3,), 400.0, 1000.0) == 1.0e-3
        μ_z(z) = 1.0e-3 * z / 1000.0
        @test CloudModels.calc_entrainment_rate((μ=μ_z,), 400.0, 1500.0) == 1.5e-3
        @test CloudModels.calc_entrainment_rate((μ=μ_z,), 100.0, 1500.0) == 1.5e-3
        @test_throws ArgumentError CloudModels.calc_entrainment_rate(
            (μ=μ_z, β=0.2), 400.0, 1500.0
        )
    end

    z = 1000.0
    T_c = env(z, :T) + 0.5
    F = ComponentArray(
        r=300.0, w=1.0, T=T_c,
        q_v=1.001 * CloudModels.calc_qv_sat(T_c, env(z, :p)),
        q_l=1.0e-4, q_r=0.0, q_i=0.0, q_pr=0.0, p=0.0
    )
    function eval_dFdz(params)
        dFdz = zero(F)
        CloudModels.plume_equations!(dFdz, copy(F), (environment=env, params...), z)
        return dFdz
    end

    @testset "plume equations" begin
        # a prescribed μ(z) equal to β/r gives the same tendencies as β
        β = 0.2
        @test eval_dFdz((μ=_ -> β / F.r,)) ≈ eval_dFdz((β=β,))
        @test eval_dFdz((μ=0.0,)) ≈ eval_dFdz((β=0.0,))
        # μ(z) is evaluated at the plume height
        @test eval_dFdz((μ=z_ -> z_ == z ? 1.0e-3 : NaN,)) ≈ eval_dFdz((μ=1.0e-3,))
        @test eval_dFdz((μ=1.0e-3,)) != eval_dFdz((β=β,))
    end

    @testset "parcel integration" begin
        # start at the lifting condensation level of a surface parcel
        T0 = env(0.0, :T) + 0.2
        qv0 = 0.7 * CloudModels.calc_qv_sat(T0, env(0.0, :p))
        z_clb, T_clb = CloudModels.compute_LCL(env, T0, qv0)
        U0 = ComponentArray(
            r=300.0, w=1.0, T=T_clb, q_v=qv0, q_l=0.0, q_r=0.0, q_i=0.0,
            q_pr=0.0, p=0.0, z=z_clb
        )
        function z_top(params)
            stop = DiscreteCallback((u, t, i) -> u.w < 0.05 || u.r <= 0.0, terminate!)
            prob = ODEProblem(
                CloudModels.parcel_equations!, U0, [0.0, 3000.0],
                (environment=env, params...)
            )
            sol = solve(prob, Euler(), dt=1.0, callback=stop)
            return maximum(getindex.(sol.u, :z))
        end
        # μ(z) == β/r(z) isn't expressible as a profile, but zero is
        @test z_top((μ=z -> 0.0,)) == z_top((β=0.0,))
        # stronger prescribed entrainment above cloud base lowers cloud top
        z_tops = [z_top((μ=z -> μ0,)) for μ0 in [0.0, 1.0e-3, 3.0e-3]]
        @test z_tops[1] > z_tops[2] > z_tops[3]
        @test z_tops[1] > z_clb
    end
end
