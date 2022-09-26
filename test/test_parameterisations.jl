using Test
using CloudModels

@testset "parameterisations" begin
    @testset "viscosity" begin
        T0 = 300.0
        μ0 = CloudModels.calc_dynamic_viscosity(T0)
        μ1 = CloudModels.calc_dynamic_viscosity(T0 + 5.0)
        # dynamic viscosity should go up with temperature
        @test μ1 > μ0
    end

end