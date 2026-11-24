using Test
using CloudModels
using ComponentArrays

include("test_parcel_equations.jl")
include("test_parameterisations.jl")
include("test_heat_capacities.jl")
include("test_eos.jl")
include("test_microphysics.jl")
include("test_profiles.jl")
include("test_plume_equations.jl")

# Not part of the suite (exploratory MCMC scripts, slow, need Turing/StatsPlots):
# - test_parcel_inference.jl
# - test_plume_inference.jl (also currently broken, see note in the file)