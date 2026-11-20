# Fixed Gauss-Legendre nodes and graded log-s panel edges for the vectorized
# fast-path thermal integrators (svavg_grid, ggxxrate).  Ports of
# gl_nodes_weights.m (Golub-Welsch) and loggrid_edges.m.

using LinearAlgebra: SymTridiagonal, eigen

const _GL_CACHE = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()

"""
    gl_nodes_weights(n) -> (x, w)

n-point Gauss-Legendre nodes/weights on [-1,1] via Golub-Welsch
(eigendecomposition of the Legendre Jacobi matrix).  Cached per n.
"""
function gl_nodes_weights(n::Integer)
    haskey(_GL_CACHE, n) && return _GL_CACHE[n]
    k = collect(1:n-1)
    beta = k ./ sqrt.(4 .* k .^ 2 .- 1)
    E = eigen(SymTridiagonal(zeros(n), beta))
    x = E.values                       # ascending
    w = 2 .* (E.vectors[1, :]) .^ 2
    res = (collect(x), collect(w))
    _GL_CACHE[n] = res
    return res
end

"""
    loggrid_edges(L_lo, L_hi, L_extra) -> edges

Sorted panel edges in log(s).  Geometrically graded above L_lo (the thermal
kernel concentrates in a thin layer of width ~ T/mx just above threshold), with
extra breakpoints (meson thresholds, BW sub-panels, notch edges) merged in.
"""
function loggrid_edges(L_lo::Real, L_hi::Real, L_extra=Float64[])
    span = L_hi - L_lo
    d = exp10.(range(-4, log10(span); length=14))      # distances above L_lo
    thr = L_lo .+ d[(d .> 0) .& (d .< span)]
    return sort(unique(vcat(L_lo, L_hi, collect(Float64.(L_extra)), thr)))
end
