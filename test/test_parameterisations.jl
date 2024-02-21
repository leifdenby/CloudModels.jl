using Test
using CloudModels
using BenchmarkTools
#using Unitful

@testset "parameterisations" begin
    @testset "viscosity" begin
        T0 = 300.0u"K"
        μ0 = CloudModels.calc_dynamic_viscosity(T0)
        @test unit(μ0) === u"kg/m/s"
        μ1 = CloudModels.calc_dynamic_viscosity(T0 + 5.0u"K")
        # dynamic viscosity should go up with temperature
        @test μ1 > μ0

        trial = @benchmark CloudModels.calc_dynamic_viscosity($T0)
        @test trial.allocs == 0
    end
end