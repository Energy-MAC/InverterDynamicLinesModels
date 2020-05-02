function instantiate_ode(M::ModelOperatingPoint; tspan::Tuple)
    model = get_ode_system()
    initial_conditions = Array{Pair}(undef, length(M.u0))
    for (ix, val) in enumerate(M.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    return DiffEqBase.ODEProblem(
        model,
        initial_conditions,
        tspan,
        M.parameters,
        jac = true,
    )
end

function instantiate_ode(system::PSY.System; tspan::Tuple, kwargs...)
    model = get_ode_system()
    steady_state = instantiate_model(system; solve_powerflow = true)
    initial_conditions = Array{Pair}(undef, length(steady_state.u0))
    for (ix, val) in enumerate(steady_state.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    return DiffEqBase.ODEProblem(
        model,
        initial_conditions,
        tspan,
        steady_state.parameters,
        jac = true,
    )
end
