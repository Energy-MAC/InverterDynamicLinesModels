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
        # Grid States
        ig_r(t)         # real current flowing into grid
        ig_i(t)         # imaginary current flowing into grid
        vg_from_r(t)    # real voltage from-side of line
        vg_from_i(t)    # imaginary voltage from-side of line
        vg_to_r(t)      # real voltage to-side of line
        vg_to_i(t)      # imaginary voltage to-side of line
        # Filter States
        eg_r(t) # real capacitor filter voltage
        eg_i(t) # imaginary capacitor filter voltage
        ic_r(t) # real current flowing into filter
        ic_i(t) # imaginary current flowing into filter
        if_r(t) # real current flowing out of the filter
        if_i(t) # imaginary current flowing out of the filter
        # Inner Loop Control States
        ξ_d(t)  #d-axis integrator term for outer AC/DC PI controller
        ξ_q(t)  #q-axis integrator term for outer AC/DC PI controller
        γ_d(t)  #d-axis integrator term for inner AC/DC PI controller
        γ_q(t)  #d-axis integrator term for inner AC/DC PI controller
        # VSM Control States
        θ(t)    # Outer Control Angle
        ω(t)    # Outer Control Frequency
    end

    # Definition of the variables for non-linear system. Requires https://github.com/SciML/ModelingToolkit.jl/issues/322 to eliminate
    variables = MTK.@variables begin
        # Grid States
        ig_r        # real current flowing into grid
        ig_i        # imaginary current flowing into grid
        vg_from_r   # real voltage from-side of line
        vg_from_i   # imaginary voltage from-side of line
        vg_to_r     # real voltage to-side of line
        vg_to_i     # imaginary voltage to-side of line
        # Filter States
        ef_r        # real capacitor filter voltage
        ef_i        # imaginary capacitor filter voltage
        ic_r        # real current flowing into filter
        ic_i        # imaginary current flowing into filter
        if_r        # real current flowing out of the filter
        if_i        # imaginary current flowing out of the filter
        # Inner Loop Control States
        ξ_d         # d-axis integrator term for outer AC/DC PI controller
        ξ_q         # q-axis integrator term for outer AC/DC PI controller
        γ_d         # d-axis integrator term for inner AC/DC PI controller
        γ_q         # d-axis integrator term for inner AC/DC PI controller
        # VSM Control States
        θ           # Outer Control Angle
        ω           # Outer Control Frequency
    end

    # Expressions
    pm = ef_r * ic_r + ef_q * ig_q  # AC Active Power into the filter
    qm = -eg_d * ig_q + eg_q * ig_d # AC Reactive Power into the filter
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
        # Line Equations
        (ω_b / L) * ((V_r_from[1] - V_r_to[1]) - (R * Il_r - L * Il_i))
        (ω_b / L) * ((V_i_from[1] - V_i_to[1]) - (R * Il_i + L * Il_r))
        ωb / cf * Ir_cnv - ωb / cf * Ir_filter + ωb * ω_sys * Vi_filter
        ωb / cf * Ii_cnv - ωb / cf * Ii_filter - ωb * ω_sys * Vr_filter
        ωb / cf * Ir_cnv - ωb / cf * Ir_filter + ωb * ω_sys * Vi_filter
        ωb / cf * Ii_cnv - ωb / cf * Ii_filter - ωb * ω_sys * Vr_filter
        #Filter Equations
        ωb / lf * Vr_cnv - ωb / lf * Vr_filter - ωb * rf / lf * Ir_cnv + ωb * ω_sys * Ii_cnv
        ωb / lf * Vi_cnv - ωb / lf * Vi_filter - ωb * rf / lf * Ii_cnv - ωb * ω_sys * Ir_cnv
        ωb / cf * Ir_cnv - ωb / cf * Ir_filter + ωb * ω_sys * Vi_filter
        ωb / cf * Ii_cnv - ωb / cf * Ii_filter - ωb * ω_sys * Vr_filter
        ωb / lg * Vr_filter - ωb / lg * V_tR - ωb * rg / lg * Ir_filter + ωb * ω_sys * Ii_filter
        ωb / lg * Vi_filter - ωb / lg * V_tI - ωb * rg / lg * Ii_filter - ωb * ω_sys * Ir_filter

        ### Inner Control Equations
        #𝜕ξ_d/𝜕t
        v_iref_d - eg_d
        #𝜕ξ_q/𝜕t
        v_iref_q - eg_q
        #𝜕γ_d/𝜕t
        i_hat_d - is_d
        #𝜕γ_q/𝜕t
        i_hat_q - is_q

        ### Outer Control Equations

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
