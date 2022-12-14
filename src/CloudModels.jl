module CloudModels

#using Unitful
using ComponentArrays
using LinearInterpolations
using Plots


macro u_str(unit)
    return 1.0
end

uconvert(_, v) = v
unit(v) = 1.0
ustrip(v) = v

export @u_str, uconvert, unit, ustrip


include("constants.jl")
include("eos.jl")
include("parameterisations.jl")
include("microphysics.jl")
include("profiles.jl")
include("parcel_equations.jl")
include("plots.jl")

include("Reference/Reference.jl")
include("rico.jl")

export calc_density, calc_temperature, calc_pressure


end # module CloudModels