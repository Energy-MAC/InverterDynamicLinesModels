
function solve_steady_state(initial_guess, parameter_values)
    _, model_rhs, _, variables, params = get_internal_model(nothing)
    @assert length(initial_guess) == length(model_rhs) == length(variables)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    _nl_system = MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
    nlsys_func = MTK.generate_function(_nl_system, expression = Val{false})[2]
    _parameter_values = [x.second for x in parameter_values]
    # Requires https://github.com/SciML/ModelingToolkit.jl/issues/323 to work
    #nlsys_jac = MTK.generate_jacobian(_nl_system, expression = Val{false})[2] # second is in-place
    sol = NLsolve.nlsolve(
        (out, x) -> nlsys_func(out, x, _parameter_values),
        #(out, x) -> nlsys_jac(out, x, _parameter_values),
        initial_guess,
    )
    println(sol)
    return sol.zero
end
function instantiate_initial_conditions(model, parameter_values)# system::PSY.System)
    #TODO: SolvePowerFlow here for eg_d, eg_q and others if needed.
    _initial_guess = [
        1.0,    #eg_d
        0.0,    #eg_q
        0.5,    #is_d
        0,      #is_q
        0.5,    #ig_d
        0,      #ig_q
        0.5,
        0.0,
        0.0,    #ξ_d
        0.0,    #ξ_q
        0.0,    #γ_d
        0.0,    #γ_q
        1.0,    #vdc
        0.5,    #ibat
        0.0,    #η
        0.0,    #κ
        0.0,    #M
        0.0,    #L
    ]
    _initial_conditions = solve_steady_state(_initial_guess, parameter_values)
    initial_conditions = Array{Pair}(undef, length(_initial_conditions))
    for (ix, val) in enumerate(_initial_guess)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    return initial_conditions
end
