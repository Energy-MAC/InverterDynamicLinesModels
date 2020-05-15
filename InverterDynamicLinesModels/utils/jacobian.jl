struct ModelJacobian
    J_func::Function
    J_Matrix::Matrix{Float64}
end

function get_jacobian_function_vsm()
    _, model_rhs, _, variables, params = get_internal_model_vsm(nothing)
    @assert length(model_rhs) == length(variables)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    _nl_system = MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
    nlsys_jac = MTK.generate_jacobian(_nl_system, expression = Val{false})[2] # second is in-place
    return nlsys_jac
end

function instantiate_jacobian_vsm(M::ModelOperatingPoint)
    jac = get_jacobian_function_vsm()
    jac_eval = (out, u0, params) -> jac(out, u0, params)
    param_eval = (out, params) -> jac(out, M.u0, params)
    n = length(M.u0)
    J = zeros(n, n)
    param_eval(J, M.parameters)
    return ModelJacobian(jac_eval, J)
end

function get_jacobian_function_droop()
    _, model_rhs, _, variables, params = get_internal_model_droop(nothing)
    @assert length(model_rhs) == length(variables)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    _nl_system = MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
    nlsys_jac = MTK.generate_jacobian(_nl_system, expression = Val{false})[2] # second is in-place
    return nlsys_jac
end

function instantiate_jacobian_droop(M::ModelOperatingPoint)
    jac = get_jacobian_function_droop()
    jac_eval = (out, u0, params) -> jac(out, u0, params)
    param_eval = (out, params) -> jac(out, M.u0, params)
    n = length(M.u0)
    J = zeros(n, n)
    param_eval(J, M.parameters)
    return ModelJacobian(jac_eval, J)
end

function (J::ModelJacobian)(M::ModelOperatingPoint)
    J.J_func(J.J_Matrix, M.u0, M.parameters)
    return J.J_Matrix
end

## Do the functor to see if you want to check 
"""
Only Lines Algebraic:
ix = trues(length(M.u0))
ix[1:4] .= false
states_ix = ix
vars_ix = .!ix
gy = jac[vars_ix, vars_ix]
gx = jac[vars_ix, states_ix]
fy = jac[states_ix, vars_ix]
fx = jac[states_ix, states_ix]
Jred = fx - fy*inv(gy)*gx 

Filter+Lines Algebraic:
ix = trues(length(M.u0))
ix[1:10] .= false
states_ix = ix
vars_ix = .!ix
gy = jac[vars_ix, vars_ix]
gx = jac[vars_ix, states_ix]
fy = jac[states_ix, vars_ix]
fx = jac[states_ix, states_ix]
Jred = fx - fy*inv(gy)*gx 
"""