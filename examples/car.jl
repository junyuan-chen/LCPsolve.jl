# An illustrative usage of LCPsolve for solving a consumption-saving model
# with an indivisible durable

# For model details, please see
# https://benjaminmoll.com/wp-content/uploads/2020/06/car.pdf
# The implementation follows the Matlab script
# https://benjaminmoll.com/wp-content/uploads/2020/06/car.m

using LinearAlgebra
using SparseArrays
using LCPsolve

mutable struct Household{TF<:AbstractFloat, TI<:Integer, TFN1<:Function, TFN2<:Function, TFN3<:Function}
    ρ::TF
    r::TF
    κ::TF
    p0::TF
    p1::TF
    vA::Vector{TF}
    nAgrid::TI
    da::TF
    nZgrid::TI
    vY::Vector{TF}
    u::TFN1
    u_prime::TFN2
    mU::Matrix{TF}
    c_policy::TFN3
    mCf::Matrix{TF}
    mCb::Matrix{TF}
    mC::Matrix{TF}
    Δ::TF
    mV::Matrix{TF}
    mVstar::Matrix{TF}
    mDaVf::Matrix{TF}
    mDaVb::Matrix{TF}
    mSf::Matrix{TF}
    mSb::Matrix{TF}
    smA::SparseMatrixCSC{TF, TI}
    smB::SparseMatrixCSC{TF, TI}
    smD::SparseMatrixCSC{TF, TI}
end

function Household(; ρ = 0.05,
    r = 0.045,
    κ = 0.25,
    p0 = 0.2,
    p1 = 0.1,
    aMin = -0.02,
    aMax = 3.0,
    nAgrid = 500,
    nZgrid = 2,
    vY = [0.1, 0.1],
    s = 2.0,
    Δ = 1000.0)

    vA = collect(LinRange(aMin, aMax, nAgrid))
    da = (aMax-aMin)/(nAgrid-1)
    mU = zeros(nAgrid,nZgrid)
    mCf = zeros(nAgrid,nZgrid)
    mCb = zeros(nAgrid,nZgrid)
    mC = zeros(nAgrid,nZgrid)

    u = c::Float64 -> c^(1.0-s)/(1.0-s)
    u_prime = c::Float64 -> c^(-s)
    c_policy = mDaV::Matrix{Float64} -> max.(mDaV, 1e-10).^(-1/s)
    mV = u.(vY'.+r.*vA)./ρ
    mVstar = Matrix{Float64}(undef, size(mV))
    mDaVf = Matrix{Float64}(undef, size(mV))
    mDaVb = Matrix{Float64}(undef, size(mV))
    mSf = Matrix{Float64}(undef, size(mV))
    mSb = Matrix{Float64}(undef, size(mV))
    smA = spzeros(length(mV), length(mV))
    smB = spzeros(length(mV), length(mV))
    smD = sparse((1.0/Δ+ρ)*I, length(mV), length(mV))
    return Household(ρ, r, κ, p0, p1, vA, nAgrid, da, nZgrid, vY, u, u_prime, mU, c_policy, mCf, mCb, mC, Δ, mV, mVstar, mDaVf, mDaVb, mSf, mSb, smA, smB, smD)
end

function diff_forward!(H::Household)
    @views @inbounds H.mDaVf[1:end-1,:] .= (H.mV[2:end,:] .- H.mV[1:end-1,:])./H.da
    # Just to fill the value, only used when aMax is too small
    H.mDaVf[end,:] = H.u_prime.(H.vY' .+ H.r*H.vA[end])
end

function diff_backward!(H::Household)
    @views @inbounds H.mDaVb[2:end,:] .= (H.mV[2:end,:] .- H.mV[1:end-1,:])./H.da
    # State constraint boundary condition
    H.mDaVb[1,:] = H.u_prime.(H.vY' .+ H.r*H.vA[1])
end

savings(H::Household, mC::Matrix) = H.vY' .+ H.r.*H.vA .- mC

function upwind_implicit!(H::Household)
    H.mC .= H.vY' .+ H.r.*H.vA
    H.mCf = H.c_policy(H.mDaVf)
    H.mCb = H.c_policy(H.mDaVb)
    H.mSf = savings(H, H.mCf)
    H.mSb = savings(H, H.mCb)
    for i in eachindex(H.mC)
        @inbounds H.mC[i] = H.mSf[i]>0.0 ? H.mCf[i] : H.mSb[i]<0.0 ? H.mCb[i] : H.mC[i]
    end
    H.mU .= H.u.(H.mC)
    H.mU[:,2] .= H.mU[:,2].+ H.κ
    mX = -min.(H.mSb, 0.0)./H.da
    mZ = max.(H.mSf, 0.0)./H.da
    mY = -mZ.-mX
    H.smA = spdiagm(0 => mY[:], -1 => mX[2:end], 1 => mZ[1:end-1])
    H.smB .= -H.smA .+ H.smD
end

function update_mVstar!(H::Household)
    i_buy = ceil(Int, H.p0/H.da)
    H.mVstar[i_buy+1:end, 1] .= H.mV[1:end-i_buy, 2]
    slope = (H.mVstar[i_buy+2,1]-H.mVstar[i_buy+1,1])/H.da
    H.mVstar[1:i_buy,1] .= H.mVstar[i_buy+1,1].+slope.*(H.vA[1:i_buy].-H.vA[i_buy+1])
    i_sell = ceil(Int, H.p1/H.da)
    H.mVstar[1:end-i_sell, 2] .= H.mV[i_sell+1:end, 1]
    H.mVstar[end-i_sell+1:end, 2] .= H.mV[end,1]
end

function vfi!(H::Household; tol = 1e-6, maxit = 100)
    diff = 10.0
    vz0 = zeros(size(H.smB,1))
    vq = zeros(size(H.smB,1))
    mVnew = Matrix{Float64}(undef, size(H.mV))
    lcp = LCP(H.smB, vq)
    for i = 1:maxit
        diff_forward!(H)
        diff_backward!(H)
        upwind_implicit!(H)
        update_mVstar!(H)
        vz0 .= H.mV[:] .- H.mVstar[:]
        lcp.M .= H.smB
        lcp.q .= -H.mU[:] .- H.mV[:]./H.Δ .+ H.smB*H.mVstar[:]
        res = solve!(lcp, vz0)
        z = res.sol
        mVnew .= reshape(z, H.nAgrid, 2) .+ H.mVstar
        diff = maximum(abs.(mVnew.-H.mV))
        H.mV .= mVnew
        if diff < tol
            println(" Iteration = ", i, " Max Diff = ", diff)
            println(" Converge!")
            break
        end
    end 
end

function main()
    H = Household()
    @time vfi!(H)
return H
end

H = main()
