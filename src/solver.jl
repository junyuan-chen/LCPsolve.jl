# The algorithm for solving the complementarity problem
# follows a Matlab routine written by Yuval Tassa.
# The original Matlab script can be accessed via MATLAB Central File Exchange at
# https://www.mathworks.com/matlabcentral/fileexchange/20952-lcp-mcp-solver-newton-based

# It is intended to generate computation results that are
# identical to what would be produced by the original Matlab script.

const DenseOrSparseMatrix{T} = Union{Matrix{T},SparseMatrixCSC{T,Int}}

struct FBCache{TF<:AbstractFloat}
    a::Vector{TF}
    b::Vector{TF}
    s::Vector{TF}
    da::SparseMatrixCSC{TF,Int}
    db::SparseMatrixCSC{TF,Int}
    at::Vector{TF}
    bt::Vector{TF}
    st::Vector{TF}
end

function FBCache{TF}(n::Int, nt::Int) where {TF<:AbstractFloat}
    a = Vector{TF}(undef, n)
    b = Vector{TF}(undef, n)
    s = Vector{TF}(undef, n)
    da = spdiagm(0=>Vector{TF}(undef, n))
    db = spdiagm(0=>Vector{TF}(undef, n))
    at = Vector{TF}(undef, nt)
    bt = Vector{TF}(undef, nt)
    st = Vector{TF}(undef, nt)
    return FBCache(a, b, s, da, db, at, bt, st)
end

# Hold the evaluation of Fischer-Burmeister function and its Jacobian
mutable struct FB{TF<:AbstractFloat}
    φ::Vector{TF}
    J::SparseMatrixCSC{TF,Int}
end

function FB{TF}(n::Int) where {TF<:AbstractFloat}
    φ = Vector{TF}(undef, n)
    J = spzeros(TF, n, n)
    return FB{TF}(φ, J)
end

# Norm of the residuals
ψ(fb::FB) = 0.5*dot(fb.φ, fb.φ)

"""
    LCP{TF<:AbstractFloat, TM<:DenseOrSparseMatrix{TF}}

A type for storing the specification of the linear complementarity problem
and arrays used in the process of solving the problem.
"""
struct LCP{TF<:AbstractFloat, TM<:DenseOrSparseMatrix{TF}}
    n::Int
    nt::Int
    M::TM
    x::Vector{TF}
    nx::Vector{TF}
    q::Vector{TF}
    l::Vector{TF}
    u::Vector{TF}
    bl::BitVector
    bu::BitVector
    blu::BitVector
    bf::BitVector
    il::BitVector
    iu::BitVector
    keep::BitVector
    ca::FBCache{TF}
    fb::FB{TF}
    nfb::FB{TF}

    function LCP{TF,TM}(M::TM, q::Vector{TF},
        l::Vector{TF}=fill(0.0, size(q)), u::Vector{TF}=fill(Inf, size(q))
        ) where {TF<:AbstractFloat,TM<:DenseOrSparseMatrix{TF}}
        n = length(q)
        size(M) != (n,n) &&
            throw(DimensionMismatch("matrix M must be square"))
        length(l) != n &&
            throw(DimensionMismatch("length of vector l does not match matrix M"))
        length(u) != n &&
            throw(DimensionMismatch("length of vector u does not match matrix M"))
        x = Vector{TF}(undef, n)
        nx = Vector{TF}(undef, n)
        bl = (l.>-Inf) .& (u.== Inf)
        bu = (l.==-Inf) .& (u.<Inf)
        blu = (l.>-Inf) .& (u.<Inf)
        bf = (l.==-Inf) .& (u.==Inf)
        il = BitVector(undef, n)
        iu = BitVector(undef, n)
        nt = sum(blu)
        keep = trues(n)
        ca = FBCache{TF}(n, nt)
        fb = FB{TF}(n)
        nfb = FB{TF}(n)
        return new(n, nt, M, x, nx, q, l, u, bl, bu, blu, bf, il, iu, keep, ca, fb, nfb)
    end
end

"""
    LCP(M::DenseOrSparseMatrix, q::Vector, l=fill(0.0, size(q)), u=fill(Inf, size(q)))

Specify a complementarity problem that finds a real vector `x` such that for each element indexed by `i`
```math
    l[i] < x[i] < u[i]   =>   (Mx)[i] + q[i] = 0;
    (Mx)[i] + q[i] < 0   =>   x[i] = u[i];
    (Mx)[i] + q[i] > 0   =>   x[i] = l[i].
```
The matrix `M` can be sparse but must be square.

If `l` and `u` are not provided,
the problem is specified as finding an `x` such that
```math
            x >= 0;
       Mx + q >= 0;
    x'(Mx + q) = 0.
```
"""
function LCP(M::DenseOrSparseMatrix, q::Vector,
    l::Vector=fill(0.0, size(q)), u::Vector=fill(Inf, size(q)))
    T = promote_type(eltype(M), eltype(q), eltype(l), eltype(u))
    M = typeof(M) <: Matrix ? Matrix{T}(M) : SparseMatrixCSC{T}(M)
    q = Vector{T}(q)
    l = Vector{T}(l)
    u = Vector{T}(u)
    return LCP{T,DenseOrSparseMatrix{T}}(M, q, l, u)
end

# Evaluate Fischer-Burmeister function and its Jacobian
function update!(fb::FB, lcp::LCP, x::Vector)
    M, q, l, u = lcp.M, lcp.q, lcp.l, lcp.u
    n, nt, bl, bu, blu, bf, ca = lcp.n, lcp.nt, lcp.bl, lcp.bu, lcp.blu, lcp.bf, lcp.ca
    a, b, s, da, db, at, bt, st = ca.a, ca.b, ca.s, ca.da, ca.db, ca.at, ca.bt, ca.st
    a .= x
    b .= M*x+q

    @inbounds a[bl] .= a[bl] .- l[bl]
    @inbounds a[bu] .= u[bu] .- a[bu]
    @inbounds b[bu] .= -b[bu]

    if nt > 0
        @inbounds at .= u[blu] .- x[blu]
        @inbounds bt .= -b[blu]
        st .= sqrt.(at.^2 .+ bt.^2)
        @inbounds a[blu] .= x[blu] .- l[blu]
        @inbounds b[blu] .= st .- at .- bt
        @inbounds M[blu,:] = -sparse(1:nt, blu, at./st.-ones(nt), nt, n) - sparse(1:nt,1:nt,bt./st.-ones(nt))*M[blu,:]
    end

    s .= sqrt.(a.^2 .+ b.^2)
    fb.φ .= s .- a .- b
    @inbounds fb.φ[bu] .= -fb.φ[bu]
    @inbounds fb.φ[bf] .= -b[bf]

    da.nzval .= a./s .- ones(n)
    db.nzval .= b./s .- ones(n)
    @inbounds da.nzval[bf] .= 0.0
    @inbounds db.nzval[bf] .= -1.0
    fb.J = da + db*M
end    

"""
    solve!(lcp, x0; tol, μ, μ_step, μ_min, max_iter, b_tol, store_trace, show_trace)

# Arguments
- `lcp::LCP{TF,DenseOrSparseMatrix{TF}}`: the complementarity problem to be solved.
- `x0::Vector{TF} = min.(max.(ones(TF,lcp.n), lcp.l), lcp.u)`: the initial value of `x`.
- `tol::TF = 1e-12`: tolerance level for convergence based on the norm of the residual.
- `μ::TF = 1e-3`: initial value of Levenberg-Marquardt μ coefficient.
- `μ_step::TF = 5.0`: constant by which μ is multiplied or divided.
- `μ_min::TF = 1e-5`: value below which μ is set to zero (pure Gauss-Newton).
- `max_iter::Int = 20`: maximum number of iterations.
- `b_tol::TF = 1e-6`: tolerance of degenerate complementarity.
- `store_trace::Bool = false`: store a trace of the solver's state in each iteration.
- `show_trace::Bool = false`: print a trace of the solver's state in each iteration.

Solve the complementarity problem specified by `lcp`.
Return an object of type `SolverResults` that contains the results.
In particular, the field `sol` of `SolverResults` contains the solution;
the field `converged` indicates whether convergence is reached.
"""
function solve!(lcp::LCP{TF,DenseOrSparseMatrix{TF}},
    x0::Vector{TF} = min.(max.(ones(TF,lcp.n), lcp.l), lcp.u);
    tol::TF = 1e-12,
    μ::TF = 1e-3,
    μ_step::TF = 5.0,
    μ_min::TF = 1e-5,
    max_iter::Int = 20,
    b_tol::TF = 1e-6,
    store_trace::Bool = false,
    show_trace::Bool = false) where TF<:AbstractFloat

    M, l, x, nx, u, n, nt = lcp.M, lcp.l, lcp.x, lcp.nx, lcp.u, lcp.n, lcp.nt
    fb, nfb, il, iu, keep = lcp.fb, lcp.nfb, lcp.il, lcp.iu, lcp.keep
    x .= x0
    update!(fb, lcp, x)
    ψ_old = ψ(fb)
    new_x = true
    tr = SolverTrace()

    for i = 1:max_iter
        if new_x
            # Separate out dimensions that are almost constrained
            il .= (abs.(x.-l) .< b_tol) .& (abs.(x.-l) .< abs.(x.-u)) .& (abs.(fb.φ) .< b_tol)
            iu .= (abs.(x.-u) .< b_tol) .& (.~ il) .& (abs.(fb.φ) .< b_tol)
            keep .= .~(il .| iu)
            nx[il] .= l[il]
            nx[iu] .= u[iu] 
        end

        Jk = fb.J[keep,keep]
        H = Jk'*Jk + μ*I
        Jphi = Jk'*fb.φ[keep]
        d = -(cholesky(H)\Jphi)
        nx[keep] .= x[keep] .+ d

        update!(nfb, lcp, nx)
        ψ_new = ψ(nfb)
        r = (ψ_old-ψ_new)/-(Jphi'*d + 0.5*(d'*H*d))

        # Increase/decrease μ if r too small/large
        r < 0.3 && (μ = max(μ*μ_step, μ_min))
        r > 0.8 && (μ = μ/μ_step * (μ > μ_min))
        if r > 0
            x .= nx
            fb.φ .= nfb.φ
            fb.J = nfb.J
            ψ_old = ψ_new
            new_x = true
        end
        (store_trace||show_trace) && update!(tr, i, ψ_new, r, μ, store_trace, show_trace)
        ψ_new < tol && return SolverResults(x0, x, i, ψ_new, r, μ, tol, true, tr)
        i == max_iter && return SolverResults(x0, x, max_iter, ψ_new, r, μ, tol, false, tr)
    end
end
