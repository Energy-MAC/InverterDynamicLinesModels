function instantiate_ode(M::ModelOperatingPoint{T, N}; tspan::Tuple) where {T <: InverterModel, N <: NetworkModel}
    model = get_ode_system(T, N)
    initial_conditions = Array{Pair}(undef, length(M.u0))
    parameters_ = Array{Pair}(undef, length(M.parameters))
    for (ix, val) in enumerate(M.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    for (ix, val) in enumerate(M.parameters)
        parameters_[ix] = MTK.parameters(model)[ix] => val
    end
    return DiffEqBase.ODEProblem(
        model,
        initial_conditions,
        tspan,
        parameters_,
        # Don't use jac = true https://github.com/SciML/ModelingToolkit.jl/issues/390
        jac = false,
    )
end

function instantiate_ode(system::PSY.System, ::Type{T}, ::Type{N}; tspan::Tuple) where {T <: InverterModel, N <: NetworkModel}
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
    return DiffEqBase.ODEProblem(
        model,
        initial_conditions,
        tspan,
        parameters_,
        # Don't use jac = true https://github.com/SciML/ModelingToolkit.jl/issues/390
        jac = false,
    )
end
