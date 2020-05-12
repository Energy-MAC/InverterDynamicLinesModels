function instantiate_ode_vsm(M::ModelOperatingPoint; tspan::Tuple)
    model = get_ode_system_vsm()
    initial_conditions = Array{Pair}(undef, length(M.u0))
    for (ix, val) in enumerate(M.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    return DiffEqBase.ODEProblem(model, initial_conditions, tspan, M.parameters, jac = false)
end

function instantiate_ode_vsm(system::PSY.System; tspan::Tuple, kwargs...)
    model = get_ode_system_vsm()
    steady_state = instantiate_model_vsm(system; solve_powerflow = true)
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


function instantiate_ode_droop(M::ModelOperatingPoint; tspan::Tuple)
    model = get_ode_system_droop()
    initial_conditions = Array{Pair}(undef, length(M.u0))
    for (ix, val) in enumerate(M.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    return DiffEqBase.ODEProblem(model, initial_conditions, tspan, M.parameters, jac = true)
end

function instantiate_ode_droop(system::PSY.System; tspan::Tuple, kwargs...)
    model = get_ode_system_droop()
    steady_state = instantiate_model_droop(system; solve_powerflow = true)
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
