function get_ode_system(::Type{T}, ::Type{N}) where {T <: InverterModel, N <: NetworkModel}
    model_lhs, model_rhs, states, _, params = get_internal_model(T, N)
    t = params[1]
    _eqs = model_lhs .~ model_rhs
    return MTK.ODESystem(_eqs, t, [states...], [params...][2:end])
end

function get_nonlinear_system(
    ::Type{T},
    ::Type{N},
) where {T <: InverterModel, N <: NetworkModel}
    _, model_rhs, _, variables, params = get_internal_model(T, N)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    return MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
end
