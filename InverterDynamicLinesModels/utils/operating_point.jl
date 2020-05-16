struct ModelOperatingPoint{T <: InverterModel, N <: NetworkModel}
    sys_func::Function
    u0::Vector{Float64}
    parameters::Vector{Float64}
end

function instantiate_model(
    ::Type{T},
    ::Type{N},
    system::PSY.System;
    solve_powerflow = false,
) where {T <: InverterModel, N <: NetworkModel}
    nl_sys = get_nonlinear_system(T, N)
    variable_count = length(states(nl_sys))
    nlsys_func = MTK.generate_function(nl_sys, expression = Val{false})[2]
    sys_f = (out, x, param) -> nlsys_func(out, x, param)

    if solve_powerflow
        PSY.solve_powerflow!(system, NLsolve.nlsolve)
    end

    inverter = collect(PSY.get_components(PSY.DynamicInverter, system))
    isempty(inverter) && @error("There are no inverters in the system")
    bus_voltage = PSY.get_voltage(PSY.get_bus(inverter[1]))
    bus_angle = PSY.get_angle(PSY.get_bus(inverter[1]))

    initial_guess = [
        0.0,    #il_r
        0.0,    #il_i
        bus_voltage * cos(bus_angle),    #vg_from_r
        bus_voltage * sin(bus_angle),    #vg_from_r
        0.95,   #ef_d
        -0.1,   #ef_q
        0.5,    #ic_d
        0.0,    #ic_q
        0.49,   #if_d
        -0.1,   #if_q
        0.0015, #ξ_d
        -0.07,  #ξ_q
        0.05,   #γ_d
        -0.001, #γ_q
        0.2,    #θ
        1.0,    #ω
        0.025,   #qf
    ]
    parameter_mapping = instantiate_parameters(T, system)
    parameter_values = [x.second for x in parameter_mapping]
    M = ModelOperatingPoint{T, N}(sys_f, initial_guess, parameter_values)
    M()
    return M
end

function (M::ModelOperatingPoint)(parameters::Vector{Float64})
    res = NLsolve.nlsolve((out, x) -> M.sys_func(out, x, parameters), M.u0)
    !NLsolve.converged(res) && @error("NLsolve failed to converge")
    M.parameters .= parameters
    M.u0 .= res.zero
    return M.u0
end

function (M::ModelOperatingPoint)(
    parameters::Array{Pair{Variable{ModelingToolkit.Parameter{Number}}, Float64}, 1},
)
    parameter_values = [x.second for x in parameters]
    return M(parameter_values)
end

(M::ModelOperatingPoint)() = M(M.parameters)
