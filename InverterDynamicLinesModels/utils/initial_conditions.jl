
function solve_steady_state(initial_guess, parameter_values)
    _, model_rhs, _, variables, params = get_internal_model(nothing)
    @assert length(initial_guess) == length(model_rhs) == length(variables)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    _nl_system = MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
    nlsys_func = MTK.generate_function(_nl_system, expression = Val{false})[2]
    _parameter_values = [x.second for x in parameter_values]
    nlsys_jac = MTK.generate_jacobian(_nl_system, expression = Val{false})[2] # second is in-place
    sol = NLsolve.nlsolve(
        (out, x) -> nlsys_func(out, x, _parameter_values),
        (out, x) -> nlsys_jac(out, x, _parameter_values),
        initial_guess,
    )
    println(sol)
    return sol.zero
end
function instantiate_initial_conditions(model, parameter_values) #, system::PSY.System)
    #TODO: SolvePowerFlow here for eg_d, eg_q and others if needed.
    _initial_guess = [
        0.0,    #il_r
        0.0,    #il_i
        1.0,    #vg_from_r
        0.0,    #vg_from i
        0.95,   #ef_d
        -0.1,   #ef_q
        0.5,    #ic_d
        0.0,    #ic_q
        0.49,   #if_d
        -0.1,   #if_q
        0.0015, #ξ_d
        -0.07,  #ξ_q
        0.05,    #γ_d
        -0.001,    #γ_q
        0.2,    #θ
        1.0,    #ω
        0.025,   #qf
    ]
    _initial_conditions = solve_steady_state(_initial_guess, parameter_values)
    initial_conditions = Array{Pair}(undef, length(_initial_conditions))
    for (ix, val) in enumerate(_initial_guess)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    return initial_conditions
end
