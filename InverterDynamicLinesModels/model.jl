function get_internal_model(::Nothing)
    # Model Parameters
    params = MTK.@parameters begin
        t
        # AC side quantities
        ωb      # Base Frequency
        # Line impedance
        lg      # Line reactance
        rg      # Line resistance
        cg      # Line capacitance
        #Reference set-point input
        pʳ      # Active Power Setpoint
        qʳ      # Reactive Power Setpoint
        vʳ      # Voltage Setpoint
        ωʳ      # Frequency Setpoint
        # Filter parameters
        lf      # Filter reactance
        cf      # Filter capacitance
        rf      # Filter Resistance
        rt      # Transformer resistance
        lt      # Transformer reactance
        # OuterControl Loops
        Ta
        Kq
        kω
        kq
        ωf
        # SRF Current Control
        kpc     # Current control propotional gain
        kic     # Current control integral gain
        kffi    # FeedForward Control for current
        ωad     # Active damping low pass filter cut-off frequency
        kad
        # SRF Voltage Control
        kpv     # Voltage control propotional gain
        kiv     # Voltage control integral gain
        kffv    # FeedForward Control for voltage
        # Virtual Impedance
        rv
        lv
    end

    MTK.@derivatives d'~t

    # Definition of the states
    states = MTK.@variables begin
        eg_d(t) #d-axis capacitor filter voltage
        eg_q(t) #q-axis capacitor filter voltage
        is_d(t) #d-axis current flowing into filter
        is_q(t) #q-axis current flowing into filter
        ig_d(t) #d-axis current flowing into grid
        ig_q(t) #q-axis current flowing into grid
        pf(t)   #Filtered active power measurement
        qf(t)   #Filtered reactive power measurement
        ξ_d(t)  #d-axis integrator term for outer AC/DC PI controller
        ξ_q(t)  #q-axis integrator term for outer AC/DC PI controller
        γ_d(t)  #d-axis integrator term for inner AC/DC PI controller
        γ_q(t)  #d-axis integrator term for inner AC/DC PI controller
    end

    # Definition of the variables for non-linear system. Requires https://github.com/SciML/ModelingToolkit.jl/issues/322 to eliminate
    variables = MTK.@variables begin
        eg_d #d-axis capacitor filter voltage
        eg_q #q-axis capacitor filter voltage
        is_d #d-axis current flowing into filter
        is_q #q-axis current flowing into filter
        ig_d #d-axis current flowing into grid
        ig_q #q-axis current flowing into grid
        pf   #Filtered active power measurement
        qf   #Filtered reactive power measurement
        ξ_d  #d-axis integrator term for outer AC/DC PI controller
        ξ_q  #q-axis integrator term for outer AC/DC PI controller
        γ_d  #d-axis integrator term for inner AC/DC PI controller
        γ_q  #d-axis integrator term for inner AC/DC PI controller
    end

    # Expressions
    pm = eg_d * ig_d + eg_q * ig_q  # AC Active Power Calculation
    qm = -eg_d * ig_q + eg_q * ig_d # AC Reactive Power Calculation
    ω_a = ωʳ + Dp * (pʳ - pf)  # Active Power Drop
    v_hat = vʳ + Dq * (qʳ - qf) # Reactive Power Drop
    v_iref_d = v_hat - rv * ig_d + ω_a * lv * ig_q # Inner voltage controller d PI
    v_iref_q = -rv * ig_q- ω_a * lv * ig_d # Inner voltage controller q PI
    i_hat_d = kvp * (v_iref_d - eg_d) + kvi * ξ_d - ω_a * cf * eg_q # Inner current controller d PI
    i_hat_q = kvp * (v_iref_d - eg_q) + kvi * ξ_q + ω_a * cf * eg_d # Inner current controller q PI
    v_md = kip * (i_hat_d - is_d) + kii * γ_d - ω_a * lf * is_q
    v_mq = kip * (i_hat_q - is_q) + kii * γ_q + ω_a * lf * is_d
    p_inv = v_md * is_d + v_mq * is_q
    q_inv = -v_md * is_q + v_mq * is_d

    model_rhs = [
        ### Grid forming equations
        #𝜕eg_d/𝜕t
        (ωb / cf) * (is_d - ig_d) + ω_a * ωb * eg_q
        #𝜕eg_q/𝜕t
        (ωb / cf) * (is_q - ig_q) - ω_a * ωb * eg_d
        #𝜕is_d/𝜕t
        (ωb / lf) * (v_md - eg_d) - (rf * ωb / lf) * is_d + ωb * ω_a * is_q
        #𝜕is_q/𝜕t
        (ωb / lf) * (v_mq - eg_q) - (rf * ωb / lf) * is_q - ωb * ω_a * is_d
        #𝜕ig_d/𝜕t
        (ωb / lt) * (eg_d - v_gd) - (rt * ωb / lt) * ig_d + ωb * ω_a * ig_q
        #𝜕ig_q/𝜕t
        (ωb / lt) * (eg_q - v_gq) - (rt * ωb / lt) * ig_q - ωb * ω_a * ig_d
        #𝜕pf/𝜕t
        ωz * (pm - pf)
        #𝜕qf/𝜕t
        ωz * (qm - qf)
        #𝜕ξ_d/𝜕t
        v_iref_d - eg_d
        #𝜕ξ_q/𝜕t
        v_iref_q - eg_q
        #𝜕γ_d/𝜕t
        i_hat_d - is_d
        #𝜕γ_q/𝜕t
        i_hat_q - is_q
        ## Line equations

    ]

    # Temporary until SteadyState problems are resolved.
    model_lhs = [
        d(eg_d)
        d(eg_q)
        d(is_d)
        d(is_q)
        d(ig_d)
        d(ig_q)
        d(pf)
        d(qf)
        d(ξ_d)
        d(ξ_q)
        d(γ_d)
        d(γ_q)

    ]

    return model_lhs, model_rhs, states, variables, params
end

function get_model()
    model_lhs, model_rhs, states, _, params = get_internal_model(nothing)
    t = params[1]
    return MTK.ODESystem(model_lhs .~ model_rhs, t, [states...], [params...][2:end])
end

function instantiate_model(
    model,
    tspan::Tuple,
    system::PSY.System,
)
    parameter_values = instantiate_parameters(model, system)
    initial_conditions = instantiate_initial_conditions(model, parameter_values, system)
    return DiffEqBase.ODEProblem(
        model,
        initial_conditions,
        tspan,
        parameter_values,
        jac = true,
    )
end
