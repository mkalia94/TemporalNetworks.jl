using Statistics: middle
import Plots: plot
"""
    bisection(f, a, b; fa = f(a), fb = f(b), ftol, wtol)

Bisection algorithm for finding the root ``f(x) ≈ 0`` within the initial bracket
`[a,b]`.

Returns a named tuple
`(x = x, fx = f(x), isroot = ::Bool, iter = ::Int, ismaxiter = ::Bool)`.

Terminates when either

1. `abs(f(x)) < ftol` (`isroot = true`),
2. the width of the bracket is `≤wtol` (`isroot = true`, to account for discontinuities),
3. `maxiter` number of iterations is reached. (`isroot = false, maxiter = true`).

which are tested for in the above order. Therefore, care should be taken not to make `wtol` too large.

"""
function bisection(f, a::Real, b::Real; fa::Real = f(a), fb::Real = f(b),
                   ftol = √eps(), wtol = 1e-6, maxiter = 1000)
    @assert fa * fb ≤ 0 "initial values don't bracket zero"
    @assert isfinite(a) && isfinite(b)
    _bisection(f, float.(promote(a, b, fa, fb, ftol, wtol))..., maxiter)
end

function _bisection(f, a, b, fa, fb, ftol, wtol, maxiter)
    iter = 0
    abs(fa) < ftol && return (x = a, fx = fa, isroot = true, iter = iter, ismaxiter = false)
    abs(fb) < ftol && return (x = b, fx = fb, isroot = true, iter = iter, ismaxiter = false)
    while true
        iter += 1
        m = middle(a, b)
        fm = f(m)
        abs(fm) < ftol && return (x = m, fx = fm, isroot = true, iter = iter, ismaxiter = false)
        if abs(b-a) ≤ wtol
            println("Bisection done: f($(a)) = $(f(a)), f($(b)) = $(f(b))")
            return (x = m, fx = fm, isroot = true, iter = iter, ismaxiter = false) 
        end
        if fa * fm > 0
            a, fa = m, fm
        else
            b, fb = m, fm
        end
        iter == maxiter && return (x = m, fx = fm, isroot = false, iter = iter, ismaxiter = true)
    end
end


