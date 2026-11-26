using Test
using CloudModels

# Thesis Ch 3 (eqn 3.28 and below eqn 3.33): L_v = h_v - h_l, L_f = h_l - h_i,
# and L_f(T) = L_f(T00) + (c_l - c_i) (T - T00). With L_s = h_v - h_i this
# means L_s(T) = L_v(T) + L_f(T) at every temperature.

@testset "latent heats" begin
    T00 = 273.15u"K"
    L_v0, L_s0, L_f0 = CloudModels.calc_latent_heats(T00)

    # at the reference temperature the constants are returned
    @test L_v0 == CloudModels.L_v
    @test L_s0 == CloudModels.L_s
    # latent heat of fusion at 0 C is about 334 kJ/kg
    @test isapprox(L_f0, 3.34e5u"J/kg"; rtol=0.01)

    @testset "L_s = L_v + L_f at T=$(T)" for T in [233.15u"K", 263.15u"K", T00, 300.0u"K"]
        L_v_, L_s_, L_f_ = CloudModels.calc_latent_heats(T)
        @test isapprox(L_s_, L_v_ + L_f_; rtol=1.0e-12)
        # Kirchhoff: dL_f/dT = c_l - c_i
        @test isapprox(L_f_ - L_f0, (CloudModels.cp_l - CloudModels.cp_i) * (T - T00); rtol=1.0e-9, atol=1.0e-6u"J/kg")
    end
end
