using Test
using CloudModels
using ComponentArrays
using OrdinaryDiffEq

# Regression test against pyclouds, the thesis-era Python implementation
# (github.com/leifdenby/python-cloudmodel). Reference values are in
# `pyclouds_reference.jl`: tendencies dF/dz at fixed states and a plume
# integrated from the LCL on the two-layer profile, with γ, C_D, β and l_pr set
# to pyclouds' defaults.
#
# Known differences between the two codes (found by patching each Julia
# formulation into pyclouds, after which all tendencies agree exactly). They
# are marked `@test_broken` below until it is decided which version is right:
#   A. autoconversion: pyclouds `k_c*max(q_l - q_g/ρ_g*a_c, 0)`, Julia
#      `k_c*max(q_l - a_c, 0)` (rain starts at a slightly higher q_l in Julia)
#   B. rain condensation/evaporation ventilation factor: pyclouds
#      `G2p75/λ^(1/4)`, Julia `G2p75/λ^(3/4)`; also T/273.0 vs T/T00 in the
#      dynamic viscosity
#   C. in-cloud gas density ρ_cg in dr/dz: pyclouds gas mass / gas volume
#      `p/(R_s T)`, Julia mixture mass / gas volume (smaller by q_d + q_v)
#   D. rain given the density of ice in the mixture density (q_r and q_i
#      swapped in the `calc_mixture_density` calls; cancels in dw/dz)
include("pyclouds_reference.jl")

const REF = PYCLOUDS_REF

function pyclouds_regression_setup()
    env = CloudModels.TwoLayerMoistIsentropicPBL(; REF.profile...)
    params = (environment=env, REF.params...)
    return env, params
end

const REGRESSION_SWITCHES = Dict(
    "default" => NamedTuple(),
    "no_hydrometeor_entr" => (entrain_hydrometeors=false,),
    "LSE_entr" => (entrain_moist_static_energy=false, entrain_liquid_static_energy=true),
    "no_static_energy_entr" => (entrain_moist_static_energy=false,),
    "const_latent_heats" => (temperature_dependent_latent_heats=false,),
)

function regression_dFdz!(dF, F, params, z)
    dF .= 0.0
    # pass a copy, since the microphysics may clip q_l in the state
    CloudModels.plume_equations!(dF, copy(F), params, z)
    return nothing
end

@testset "pyclouds regression" begin
    env, params = pyclouds_regression_setup()

    @testset "initial state at the LCL" begin
        T_s = env(0.0, :T) + 0.2
        qv0 = 0.7 * CloudModels.calc_qv_sat(T_s, env(0.0, :p))
        z_clb, T_clb = CloudModels.compute_LCL(env, T_s, qv0)
        @test qv0 ≈ REF.initial.q_v rtol = 1e-12
        @test z_clb ≈ REF.initial.z rtol = 1e-12
        @test T_clb ≈ REF.initial.T rtol = 1e-10
    end

    @testset "tendencies: $(sw_name)" for (sw_name, sw) in REGRESSION_SWITCHES
        p = merge(params, sw)
        for st in REF.states
            F = ComponentArray(; st.F...)
            dF = zero(F)
            regression_dFdz!(dF, F, p, st.z)
            ref = REF.tendencies[sw_name][st.name]
            has_rain = st.F.q_r > 0.0

            # momentum is identical in all states
            @test dF.w ≈ ref.w rtol = 1e-12

            if !has_rain
                # without rain the thermodynamics and microphysics are identical
                for v in (:T, :q_v, :q_l, :q_r, :q_pr)
                    @test isapprox(dF[v], ref[v], rtol=1e-12, atol=1e-20)
                end
                # dr/dz differs by difference C only
                @test dF.r ≈ ref.r rtol = 1e-4
                if st.F.q_l > 0.0
                    @test_broken dF.r ≈ ref.r rtol = 1e-10
                else
                    @test dF.r ≈ ref.r rtol = 1e-12
                end
            else
                # with rain, differences A, B and D enter all water and
                # temperature tendencies
                @test_broken all(
                    isapprox(dF[v], ref[v], rtol=1e-6) for v in (:T, :q_v, :q_l, :q_r, :r)
                )
            end
        end
    end

    @testset "plume from the LCL" begin
        F0 = ComponentArray(
            r=REF.initial.r, w=REF.initial.w, T=REF.initial.T, q_v=REF.initial.q_v,
            q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0, p=env(REF.initial.z, :p),
        )
        # stop where pyclouds' default stopping functions would: w crosses 0
        # or r exceeds 10 km
        cbs = CallbackSet(
            ContinuousCallback((u, z, i) -> u.w, nothing, terminate!),
            ContinuousCallback((u, z, i) -> u.r - 10e3, terminate!, nothing),
        )
        prob = ODEProblem(regression_dFdz!, F0, (REF.initial.z, 6000.0), params)
        sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-12, callback=cbs)

        z_top = sol.t[end]
        # pyclouds: 2843.18 m, Julia currently 0.54 m lower (A-D)
        @test abs(z_top - REF.plume.z_top) < 2.0
        @test_broken abs(z_top - REF.plume.z_top) < 0.01

        below_rain = REF.plume.z .< REF.plume.z_rain_onset
        for (k, z) in enumerate(REF.plume.z)
            u = sol(z)
            if below_rain[k]
                # only difference C acts here (relative differences ~3e-5)
                for v in (:w, :r, :T, :q_v, :q_l)
                    @test u[v] ≈ REF.plume[v][k] rtol = 1e-4
                end
            else
                # rain physics differs (A, B, D): currently within 0.8% in w
                # and r, 1.7% in q_l, 1e-5 in T, 2e-4 in q_v, and up to
                # 2.7e-6 kg/kg in q_r (rain starts later in Julia)
                @test u.w ≈ REF.plume.w[k] rtol = 1e-2
                @test u.r ≈ REF.plume.r[k] rtol = 1e-2
                @test u.T ≈ REF.plume.T[k] rtol = 1e-5
                @test u.q_v ≈ REF.plume.q_v[k] rtol = 2e-4
                @test u.q_l ≈ REF.plume.q_l[k] rtol = 2.5e-2
                @test u.q_r ≈ REF.plume.q_r[k] atol = 4e-6
            end
        end
        @test_broken all(
            isapprox(sol(z)[v], REF.plume[v][k], rtol=1e-6)
            for (k, z) in enumerate(REF.plume.z) for v in (:w, :r, :T, :q_v, :q_l, :q_r)
        )
    end
end
