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

# Get the gradients of a hand-written ODE right-hand side at a grid of points in the state
# space. `f!` is the mutating function, `u0` a prototype state whose components are addressed
# by name (typically a `ComponentArray`) and `p` its parameter array. Only `xname` and `yname`
# are overwritten with the grid values, so this also draws a 2-D slice of a larger system with
# the remaining states fixed at their `u0` values.
function get_gradient(f!, u0, p, xname, yname, xx, yy; t=0.0)
    u = copy(u0)     # never mutate the prototype passed by the caller
    du = similar(u)
    v = map(xx, yy) do x, y
        u[xname] = x
        u[yname] = y
        f!(du, u, p, t)
        (du[xname], du[yname])
    end
    dx = map(first, v)
    dy = map(last, v)
    return (dx, dy)
end

function normalize_gradient(dx, dy, xscale, yscale; transform=cbrt)
    maxdx = maximum(abs, dx)
    maxdy = maximum(abs, dy)
    dx_norm = @. transform(dx / maxdx) * xscale
    dy_norm = @. transform(dy / maxdy) * yscale
    return (dx_norm, dy_norm)
end

# Helpers for writing the ODE right-hand side of a reaction network by hand.

# The flux of a mass-action reaction: its rate constant times the mass-action monomial of its reactants. A reaction without reactants (e.g. `k, 0 --> X`) contributes its rate constant as-is.
mass_action(k) = k                 # zero-order flux, e.g. constant production
mass_action(k, x) = k * x          # first-order flux
mass_action(k, x, y) = k * x * y   # second-order flux

# Hill-type rate laws, named `hil`/`hilr` so as not to clash with Catalyst's `mm`, `hill` and `hillr`.
hil(x, k=one(x)) = x / (x + k)     # activation, saturating from 0 to 1
hil(x, k, n) = hil(x^n, k^n)       # Hill coefficient n
hilr(x, k=one(x)) = k / (x + k)    # repression, saturating from 1 to 0
hilr(x, k, n) = hilr(x^n, k^n)     # Hill coefficient n

# Replace one parameter and return a copy of the parameter array.
# Only parameters that `p` already has can be replaced this way; to introduce a new one, build the array from scratch, e.g. `ComponentArray(merge(NamedTuple(p), (new_k = 1.0,)))`.
function with_param(p, name::Symbol, val)
    q = copy(p)
    q[name] = val
    return q
end
