# An illustrative usage of LCPsolver for solving a version of the Hopenhayn Model

# For model details, please see
# https://benjaminmoll.com/wp-content/uploads/2020/06/hopenhayn.pdf
# The implementation follows the Matlab script
# https://benjaminmoll.com/wp-content/uploads/2020/06/hopenhayn.m

using LinearAlgebra
using SparseArrays
using LCPsolver

# Discretized firm productivity
struct Grid{TF<:AbstractFloat}
    nZgrid::Int
    vZ::Vector{TF}
    dz::TF
    dz2::TF
end

function Grid(;
    nZgrid = 1000,
    Zmin = 0.0,
    Zmax = 1.0)

    vZ = collect(LinRange(Zmin, Zmax, nZgrid))
    dz = vZ[2]-vZ[1]
    dz2 = dz^2
    return Grid(nZgrid, vZ, dz, dz2)
end

# Parameters for the problem
struct Param{TF<:AbstractFloat}
    ρ::TF
    α::TF
    ϵ::TF
    φ::TF
    vstar::TF
    c_f::TF
    c_e::TF
    m_bar::TF
    η::TF
    vψ::Vector{TF}
    μ::TF
    σ_bar::TF
    vσ2::Vector{TF}
end

function Param(g::Grid;
    ρ = 0.05,
    α = 0.5,
    ϵ = 0.5,
    φ = 0.5,
    vstar = 0.0,
    c_f = 0.05,
    c_e = 0.6,
    m_bar = 0.1,
    η = 1000.0,
    ψ_lb = 0.7,
    μ = -0.01,
    σ_bar = 0.01)

    iψ_lb = ceil(Int, ψ_lb*g.nZgrid)
    vψ = zeros(g.nZgrid)
    vψ[iψ_lb:end] .= 1/((g.nZgrid-iψ_lb+1)*g.dz)
    vσ2 = (g.vZ.*σ_bar).^2
    return Param(ρ, α, ϵ, φ, vstar, c_f, c_e, m_bar, η, vψ, μ, σ_bar, vσ2)
end

# Quantities for each firm indexed by productivity
struct Firm{TF<:AbstractFloat}
    vN::Vector{TF}  # Labor demand
    vF::Vector{TF}  # Output
    vΠ::Vector{TF}  # Profit
    vV::Vector{TF}  # Value
    vG::Vector{TF}  # Distribution
end

function Firm(g::Grid)
    vN = zeros(g.nZgrid)
    vF = zeros(g.nZgrid)
    vΠ = zeros(g.nZgrid)
    vV = zeros(g.nZgrid)
    vG = zeros(g.nZgrid)
    return Firm(vN, vF, vΠ, vV, vG)
end

function update!(f::Firm, p::Float64, w::Float64, g::Grid, pa::Param)
    f.vN .= (pa.α.*g.vZ.*p./w).^(1.0 /(1.0-pa.α))
    f.vF .= g.vZ.*f.vN.^pa.α
    f.vΠ .= p.*f.vF .- w.*f.vN .- pa.c_f
end

struct UpwindArrays{TF<:AbstractFloat}
    A::SparseMatrixCSC{TF, Int}
    B::SparseMatrixCSC{TF, Int}
    q::Vector{TF}
end

function UpwindArrays(g::Grid, p::Param)
    X = -min.(p.μ,0)./g.dz .+ p.vσ2./(2.0*g.dz2)
    Y = -max.(p.μ,0)./g.dz .+ min.(p.μ,0)./g.dz .- p.vσ2./g.dz2
    Z = max.(p.μ,0)./g.dz .+ p.vσ2./(2.0*g.dz2)

    Y[1] = Y[1] + X[1]
    Y[end] = Y[end] + Z[end]
    A = spdiagm(-1 => X[2:end], 0 => Y, 1 => Z[1:end-1])
    B = sparse(p.ρ*I, g.nZgrid, g.nZgrid) - A
    q = zeros(g.nZgrid)
    return UpwindArrays(A, B, q)
end

function vfi!(g::Grid, pa::Param, f::Firm; p = 0.6, w = 1.0, relax_w = 0.2, relax_p = 0.001, max_iter = 100, tol = 1e-5)
    ua = UpwindArrays(g, pa)
    AA = copy(ua.A) # KF array
    vS = zeros(g.nZgrid)
    lcp = LCP(ua.B, ua.q)
    Bvstar = reshape(sum(ua.B, dims = 2), g.nZgrid).*pa.vstar
    x0 = zeros(g.nZgrid)
    
    for i_w = 1:max_iter
        for i_p = 1:max_iter
            update!(f, p, w, g, pa)
            lcp.q .= -f.vΠ .+ Bvstar
            s = solve!(lcp, x0)            
            vS .= s.sol
            f.vV .= vS .+ pa.vstar

            # Find the threshold for firm exit
            ix = findfirst(x->x!=0.0, vS) - 1

            m = pa.m_bar*exp(pa.η*(g.dz*dot(f.vV, pa.vψ)-pa.c_e))
            # Reconstruct AA with the new threshold
            AA.nzval .= ua.A.nzval
            AA.nzval[1:3*ix-1] .= 0.0
            AA.nzval[1:3:3*ix-2] .= 1.0
            f.vG .= -AA'\(m.*pa.vψ)

            Q = max(dot(f.vF, f.vG.*g.dz), 1e-4)
            p_implied = Q^(-pa.ϵ)
            diff_p = abs(p-p_implied)
            # println(" p iteration ", i_p, " diff_p = ", diff_p, " p = ", p)
            p = relax_p*p_implied + (1.0-relax_p)*p
            diff_p < tol && break
        end
        N = max(dot(f.vN, f.vG.*g.dz), 1e-4)
        w_implied = N^pa.φ

        diff_w = abs(w-w_implied)
        # println(" w iteration ", i_w, " diff_w = ", diff_w, " w = ", w)
        w = relax_w*w_implied + (1.0-relax_w)*w
        if diff_w < tol
            println("Equilibrium found, w = ", w, " p = ", p)
            break
        end
    end
end

function main()
    g = Grid()
    pa = Param(g)
    f = Firm(g)
    @time vfi!(g, pa, f)
    return f, g
end

f, g = main()
