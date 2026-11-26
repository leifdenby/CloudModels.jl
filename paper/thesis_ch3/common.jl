# Shared set-up, integration and plotting for the thesis Ch 3 figure scripts
# (Denby 2017, PhD thesis, Figs 3.2-3.8). See `README.md` in this folder.
#
# Every figure script `include`s this file, builds a list of runs with
# `run_plume` and plots them with `plot_profiles`.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")  # headless GR

using CloudModels
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Printf

const CM = CloudModels

# Thesis Ch 3 values (Sect. 3.1.2: "γ = 0.5, C_D = 0.5"); the package default
# for γ is 1.0, so it is always passed explicitly here.
const THESIS_GAMMA = 0.5

# Command-line arguments (all optional):
#   1. output directory (default `output/` next to this file)
#   2. RH gradient above the boundary layer for the figures that use the
#      default profile, in 1/m (default -0.5e-3, see `thesis_profile`)
outdir() = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "output")
const DRHDZ_DEFAULT = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : -0.5e-3

"""
    thesis_profile(; z_BL=600, z_INV=2400, dRHdz_1=DRHDZ_DEFAULT)

The idealised two-layer profile used throughout thesis Sect. 3.4 (figure
titles: RH0 = 70 %, T0 = 299.3 K). The surface pressure is p0 = 101325 Pa
(thesis Sect. 5.1.2).

The relative-humidity gradient above the boundary layer is not stated in the
thesis. The environment q_v and moist static energy lines in Figs 3.2, 3.4,
3.5, 3.7 and 3.8 (q_v,e ≈ 13.5 g/kg at 700 m, ≈ 3.5 g/kg at 1850 m, 0 near
2500 m; s_e ≈ 310 kJ/kg at z_INV) match `dRHdz_1 = -0.5e-3` 1/m (-50 %/km),
not the package's `TwoLayerMoistIsentropicPBL` default of -0.1e-3, so that
is the default here. Fig 3.6 sets its own values. The cloud profiles in Fig
3.2 lie between runs with -0.5e-3 and -0.1e-3, so the value can be changed
with the second command-line argument.
"""
function thesis_profile(; z_BL=600.0, z_INV=2400.0, dRHdz_1=DRHDZ_DEFAULT)
    return CM.TwoLayerMoistIsentropicPBL(;
        z_BL=z_BL, RH0=0.7, T0=299.3, z_INV=z_INV, dRHdz_1=dRHdz_1, p0=101325.0,
    )
end

"""
    cloudbase_state(env; r0, w0, dT=0, dqv=0)

Initial state at cloud base: a surface parcel with the environment's surface
temperature and water vapour plus the perturbation (`dT` [K], `dqv` [kg/kg])
is lifted dry-adiabatically to its lifting condensation level
(`CloudModels.compute_LCL`). Returns `(z_base, F0)`.
"""
function cloudbase_state(env; r0, w0, dT=0.0, dqv=0.0)
    T_s = env(0.0, :T) + dT
    qv_s = env(0.0, :qv) + dqv
    z_b, T_b = CM.compute_LCL(env, T_s, qv_s)
    F0 = ComponentArray(
        r=r0, w=w0, T=T_b, q_v=qv_s, q_l=0.0, q_r=0.0, q_i=0.0, q_pr=0.0,
        p=env(z_b, :p),
    )
    return z_b, F0
end

# `plume_equations!` on this branch accumulates into `dFdz` (the
# `dFdz .= 0.0` fix is in PR #19), and the microphysics clips q_l in the state
# it is given, so zero the tendencies and pass a copy.
function plume_rhs!(dF, F, params, z)
    dF .= 0.0
    CM.plume_equations!(dF, copy(F), params, z)
    return nothing
end

"""
    run_plume(env, F0, z0, params; z_max, label)

Integrate the plume equations in height from `z0` until the updraft stops
(w crosses 0), the radius exceeds 10 km (pyclouds' default stopping
conditions) or `z_max` is reached. Returns a NamedTuple with the solution,
the cloud top and why the integration stopped.
"""
function run_plume(env, F0, z0, params; z_max=4000.0, label="", dtmax=5.0)
    stop = Ref("z_max reached")
    cbs = CallbackSet(
        ContinuousCallback(
            (u, z, i) -> u.w, nothing,
            i -> (stop[] = "w = 0"; terminate!(i)),
        ),
        ContinuousCallback(
            (u, z, i) -> u.r - 10e3,
            i -> (stop[] = "r > 10 km"; terminate!(i)), nothing,
        ),
    )
    p = (environment=env, γ=THESIS_GAMMA, params...)
    prob = ODEProblem(plume_rhs!, copy(F0), (z0, z_max), p)
    sol = solve(
        prob, Tsit5(); reltol=1e-8, abstol=1e-10, dtmax=dtmax, callback=cbs,
        maxiters=10^6,
    )
    if !(sol.retcode in (ReturnCode.Success, ReturnCode.Terminated))
        stop[] = "solver: $(sol.retcode)"
    end
    return (sol=sol, z_base=z0, z_top=sol.t[end], stop=stop[], label=label, params=p, env=env)
end

# ---------------------------------------------------------------------------
# diagnostics

"""
Moist static energy (thesis eq. 3.31, no ice) [J/kg]:
s = c̄_l T + q_v L_v + g z with c̄_l = q_d c_p,d + (q_v + q_l + q_r) c_l and
constant L_v, as in pyclouds' plotting code (`pyclouds/plot/parcel.py`).
"""
function moist_static_energy(T, q_v, q_l, q_r, z)
    q_d = 1.0 - q_v - q_l - q_r
    c_l = q_d * CM.cp_d + (q_v + q_l + q_r) * CM.cp_l
    return c_l * T + q_v * CM.L_v + CM.g * z
end

"""
Diagnostic profiles of a run, sampled at the solver's steps. Densities use
`calc_mixture_density(p, T, qd, qv, ql, qr, qi)` in its own argument order
(rain as liquid); the environment density is the same expression with the
environment's T and q_v (as in `plume_equations!`).
"""
function diagnostics(run)
    sol, env = run.sol, run.env
    z = sol.t
    u = sol.u
    r = [x.r for x in u]
    w = [x.w for x in u]
    T = [x.T for x in u]
    qv = [x.q_v for x in u]
    ql = [x.q_l for x in u]
    qr = [x.q_r for x in u]
    p = [env(zz, :p) for zz in z]
    qd = 1.0 .- qv .- ql .- qr
    rho_c = CM.calc_mixture_density.(p, T, qd, qv, ql, qr, 0.0)
    qv_e = [env(zz, :qv) for zz in z]
    T_e = [env(zz, :T) for zz in z]
    rho_e = CM.calc_mixture_density.(p, T_e, 1.0 .- qv_e, qv_e, 0.0, 0.0, 0.0)
    qv_sat = CM.calc_qv_sat.(T, p)
    return (
        z=z, r=r, w=w, T=T, qv=qv, ql=ql, qr=qr,
        mse=moist_static_energy.(T, qv, ql, qr, z),
        drho=rho_c .- rho_e,
        ss=(qv ./ qv_sat .- 1.0) .* 100.0,
        dqv=qv .- qv_e,
        rh_c=qv ./ qv_sat .* 100.0,
    )
end

function env_profiles(env, z)
    T = [env(zz, :T) for zz in z]
    qv = [env(zz, :qv) for zz in z]
    return (
        z=z, T=T, qv=qv,
        mse=moist_static_energy.(T, qv, 0.0, 0.0, z),
        rh=[env(zz, :rh) for zz in z] .* 100.0,
    )
end

# panel definitions: (diagnostic field, scale, x-label, environment field)
const PANELS = Dict(
    :r => (:r, 1.0, "radius [m]", nothing),
    :w => (:w, 1.0, "vertical velocity [m/s]", nothing),
    :T => (:T, 1.0, "temperature [K]", :T),
    :qv => (:qv, 1e3, "water vapour [g/kg]", :qv),
    :ql => (:ql, 1e3, "liquid water [g/kg]", nothing),
    :qr => (:qr, 1e3, "rain water [g/kg]", nothing),
    :mse => (:mse, 1e-3, "moist static energy [kJ/kg]", :mse),
    :drho => (:drho, 1e3, "ρ_c - ρ_e [g/m³]", nothing),
    :ss => (:ss, 1.0, "supersaturation [%]", nothing),
    :dqv => (:dqv, 1e3, "q_v,c - q_v,e [g/kg]", nothing),
    :rh => (:rh_c, 1.0, "relative humidity [%]", :rh),
)

const PALETTE = [
    "#4c72b0", "#55a868", "#c44e52", "#8172b2", "#ccb974", "#64b5cd",
]

"""
    plot_profiles(runs, panels; title, ylims, xlims=Dict(), envs=nothing)

One subplot per panel symbol (see `PANELS`), one line per run, with the
environment profile where one applies (`envs`: one environment per run, or a
single one shared by all).
"""
function plot_profiles(runs, panels; title="", ylims=nothing, xlims=Dict(),
                       ncols=3, colors=PALETTE, envs=nothing)
    diags = diagnostics.(runs)
    z_lo = minimum(r.z_base for r in runs) - 100.0
    z_hi = isnothing(ylims) ? maximum(r.z_top for r in runs) + 200.0 : ylims[2]
    yl = isnothing(ylims) ? (z_lo, z_hi) : ylims
    envs = isnothing(envs) ? [runs[1].env] : envs
    env_line_styles = [:solid, :dash, :dot, :dashdot]

    plts = []
    for (k, pn) in enumerate(panels)
        field, scale, xlabel, env_field = PANELS[pn]
        plt = plot(; xlabel=xlabel, ylabel=(k - 1) % ncols == 0 ? "height [m]" : "",
                   ylims=yl, legend=false, guidefontsize=8, tickfontsize=7)
        if !isnothing(env_field)
            z_env = collect(range(max(yl[1], 0.0), yl[2]; length=300))
            for (ie, env) in enumerate(unique(envs))
                ep = env_profiles(env, z_env)
                plot!(plt, getfield(ep, env_field) .* scale, z_env;
                      color=:gray50, lw=1.2, ls=env_line_styles[mod1(ie, 4)])
            end
        end
        if pn in (:ql, :qr, :drho, :ss, :dqv, :qv)
            vline!(plt, [0.0]; color=:gray70, ls=:dot, lw=1)
        end
        for (i, d) in enumerate(diags)
            plot!(plt, getfield(d, field) .* scale, d.z;
                  color=colors[mod1(i, length(colors))], lw=2)
        end
        if haskey(xlims, pn)
            xlims!(plt, xlims[pn]...)
        end
        push!(plts, plt)
    end
    # legend panel
    leg = plot(; framestyle=:none, legend=:inside, legendfontsize=7,
               foreground_color_legend=nothing, xlims=(0, 1), ylims=(0, 1))
    for (i, r) in enumerate(runs)
        plot!(leg, [NaN], [NaN]; color=colors[mod1(i, length(colors))], lw=2,
              label=@sprintf("%s (top %.0f m, %s)", r.label, r.z_top, r.stop))
    end
    plot!(leg, [NaN], [NaN]; color=:gray50, lw=1.2, label="environment")
    push!(plts, leg)

    nrows = cld(length(plts), ncols)
    while length(plts) < nrows * ncols
        push!(plts, plot(; framestyle=:none))
    end
    return plot(plts...; layout=(nrows, ncols), size=(330 * ncols, 300 * nrows),
                plot_title=title, plot_titlefontsize=10, left_margin=4Plots.mm,
                bottom_margin=3Plots.mm)
end

"""
Print a one-line summary per run and write it as CSV to `<name>_summary.csv`.
"""
function summarise(name, runs)
    lines = String[]
    csv = ["figure,run,z_base,z_top,max_w,max_r,stop,dRHdz_default"]
    push!(lines, @sprintf("%-40s %8s %8s %8s %8s  %s", "run", "z_base", "z_top", "max w", "max r", "stopped by"))
    for r in runs
        d = diagnostics(r)
        push!(lines, @sprintf("%-40s %8.0f %8.0f %8.2f %8.0f  %s",
                              r.label, r.z_base, r.z_top, maximum(d.w), maximum(d.r), r.stop))
        push!(csv, @sprintf("%s,\"%s\",%.1f,%.1f,%.3f,%.1f,%s,%g", name, r.label,
                            r.z_base, r.z_top, maximum(d.w), maximum(d.r), r.stop, DRHDZ_DEFAULT))
    end
    txt = join(lines, "\n")
    println("== $(name) (dRHdz_1 default $(DRHDZ_DEFAULT) 1/m)\n", txt)
    mkpath(outdir())
    write(joinpath(outdir(), "$(name)_summary.csv"), join(csv, "\n") * "\n")
    return txt
end

function save_figure(plt, name)
    mkpath(outdir())
    path = joinpath(outdir(), "$(name).png")
    savefig(plt, path)
    println("wrote ", path)
    return path
end
