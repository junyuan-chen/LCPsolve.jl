struct SolverState{T<:Real}
    iter::Int
    ψ::T
    r::T
    μ::T
end

struct SolverTrace
    states::Vector{SolverState}
end

SolverTrace() = SolverTrace(Array{SolverState}(undef, 0))

function show(io::IO, st::SolverState)
    @printf io "iter = %2d, ψ = %3.0e, r = %3.1f, μ = %3.0e\n" st.iter st.ψ st.r st.μ
end

push!(tr::SolverTrace, st::SolverState) = push!(tr.states, st)

getindex(tr::SolverTrace, i::Integer) = getindex(tr.states, i)

function show(io::IO, tr::SolverTrace)
    for state in tr.states
        show(io, state)
    end
end

function update!(tr::SolverTrace,
                 iter::Int,
                 ψ::Real,
                 r::Real,
                 μ::Real,
                 store_trace::Bool,
                 show_trace::Bool)
    st = SolverState(iter, ψ, r, μ)
    store_trace && push!(tr, st)
    show_trace && show(st)
end

struct SolverResults{T<:Real}
    initial_x::Vector{T}
    sol::Vector{T}
    iter::Int
    ψ::T
    r::T
    μ::T
    tol::T
    converged::Bool
    trace::SolverTrace
end

function show(io::IO, re::SolverResults)
    @printf io "Results of Solving LCP\n"
    @printf io " * Convergence: %s\n" re.converged
    @printf io " * iter = %2d, ψ = %3.0e, r = %3.1f, μ = %3.0e" re.iter re.ψ re.r re.μ
end
