# Convenience functions to plot gradients of the ODEs.
function meshgrid(xrange, yrange)
    xx = [x for y in yrange, x in xrange]
    yy = [y for y in yrange, x in xrange]
    return xx, yy
end

# Get the gradients of the vector field of an ODE system at a grid of points in the state space
function get_gradient(prob, xsym, ysym, xx, yy; t=nothing)
    ## The order of state variables (unknowns) in the ODE system is not guaranteed. So we may need to swap the order of x and y when calling ∂F.
    swap_or_not(x, y; xidx=1) = xidx == 1 ? [x, y] : [y, x]
    ∂F = prob.f
    ps = prob.p
    sys = prob.f.sys
    xidx = SymbolicIndexingInterface.variable_index(sys, xsym)
    yidx = SymbolicIndexingInterface.variable_index(sys, ysym)
    dx = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[xidx], xx, yy)
    dy = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[yidx], xx, yy)
    return (dx, dy)
end

function normalize_gradient(dx, dy, xscale, yscale; transform=cbrt)
    maxdx = maximum(abs, dx)
    maxdy = maximum(abs, dy)
    dx_norm = @. transform(dx / maxdx) * xscale
    dy_norm = @. transform(dy / maxdy) * yscale
    return (dx_norm, dy_norm)
end
