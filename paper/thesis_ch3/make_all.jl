# Run all thesis Ch 3 figure scripts:
#
#   julia --project=. paper/thesis_ch3/make_all.jl [outdir] [dRHdz_1]
include(joinpath(@__DIR__, "common.jl"))

for name in ("fig3_2", "fig3_3", "fig3_4", "fig3_5", "fig3_6", "fig3_7", "fig3_8")
    include(joinpath(@__DIR__, "$(name).jl"))
end
