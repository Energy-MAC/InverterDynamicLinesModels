abstract type Perturbation end
mutable struct Simulation
    problem::DiffEqBase.ODEProblem
    perturbations::Vector{<:Perturbation}
end

struct CircuitTrip <: Perturbation
    time::Float64
    droppped_circuits::Float64
    callback::DiffEqBase.DiscreteCallback
        function CircuitTrip(;time::Float64, droppped_circuits::Float64)
            condition(u, t, integrator) = t == time
            function affect!(integrator)
                integrator.p[4] = integrator.p[4]*droppped_circuits+1
                integrator.p[5] = integrator.p[5]*droppped_circuits+1
                integrator.p[4] = integrator.p[6]/(droppped_circuits+1)
                integrator.p[5] = integrator.p[7]/(droppped_circuits+1)
            end
            return CircuitTrip(
                   time,
                   droppped_circuits,
                   DiffEqBase.DiscreteCallback(condition,affect!)
            )
        end
end
CirtuitTrip(t::Float64) = CirtuitTrip(t, 1.0)

struct PowerOutputIncrease <: Perturbation
    time::Float64
    new_value::Float64
    callback::DiffEqBase.DiscreteCallback
        function PowerOutputIncrease(;time::Float64, droppped_circuits::Float64)
            condition(u, t, integrator) = t == time
            function affect!(integrator)
                integrator.p[10] = new_value
            end
            return PowerOutputIncrease(
                   time,
                   droppped_circuits,
                   DiffEqBase.DiscreteCallback(condition,affect!)
            )
        end
end

function instantiate_ode(M::ModelOperatingPoint{T, N}, perturbations::Perturbation...; tspan::Tuple) where {T <: InverterModel, N <: NetworkModel}
    model = get_ode_system(T, N)
    initial_conditions = Array{Pair}(undef, length(M.u0))
    parameters_ = Array{Pair}(undef, length(M.parameters))
    for (ix, val) in enumerate(M.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    for (ix, val) in enumerate(M.parameters)
        parameters_[ix] = MTK.parameters(model)[ix] => val
    end
    problem = DiffEqBase.ODEProblem(
        model,
        initial_conditions,
        tspan,
        parameters_,
        # Don't use jac = true https://github.com/SciML/ModelingToolkit.jl/issues/390
        jac = false,
    )
    return Simulation(problem, [perturbations...])
end

function instantiate_ode(system::PSY.System, ::Type{T}, ::Type{N}, perturbations::Perturbation...; tspan::Tuple) where {T <: InverterModel, N <: NetworkModel}
    model = get_ode_system(T, N)
    steady_state = instantiate_model(T, N, system; solve_powerflow = true)
    initial_conditions = Array{Pair}(undef, length(steady_state.u0))
    parameters_ = Array{Pair}(undef, length(steady_state.parameters))
    for (ix, val) in enumerate(steady_state.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    for (ix, val) in enumerate(steady_state.parameters)
        parameters_[ix] = MTK.parameters(model)[ix] => val
    end
    problem = DiffEqBase.ODEProblem(
        model,
        initial_conditions,
        tspan,
        parameters_,
        # Don't use jac = true https://github.com/SciML/ModelingToolkit.jl/issues/390
        jac = false,
    )
    Simulation(problem, [perturbations...])
end
