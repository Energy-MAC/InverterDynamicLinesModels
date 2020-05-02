function instantiate_ode(system::PSY.System; tspan::Tuple, kwargs...)
    model = get_ode_system()
    parameter_values = instantiate_parameters(model, system)
    solve_pf = get(kwargs, :solve_powerflow, false)
    initial_conditions = instantiate_initial_conditions(
        model,
        parameter_values,
        system;
        solve_powerflow = solve_pf,
    )
    return DiffEqBase.ODEProblem(
        model,
        initial_conditions,
        tspan,
        parameter_values,
        jac = true,
    )
end
