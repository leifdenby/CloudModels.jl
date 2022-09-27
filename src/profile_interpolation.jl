
y0 = ComponentArray(r=100, w=2)
y1 = ComponentArray(r=200, w=-1)

ys = [y0, y1]
xs = [0, 10]

itp = Interpolate(xs, ys)
itp(5.0)
itp

reference_points = Dict(
    0.0 => ComponentArray(qt=16.0, theta_l=797.9),
    740.0 => ComponentArray(qt=13.8, theta_l=797.9),
    3260.0 => ComponentArray(qt=2.4, theta_l=317.0),
    4000.0 => ComponentArray(qt=2.4, theta_l=317.0),
)

values(reference_points)

itp = Interpolate(collect(keys(reference_points)), collect(values(reference_points)))

itp(300.0)