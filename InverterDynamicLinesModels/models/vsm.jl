function get_internal_model_vsm(::Nothing)
    # Model Parameters
    params = MTK.@parameters begin
        t
        # AC side quantities
        Ωb      # Base Frequency
        ω_sys   # Grid frequency in pu
        Sb      # Base Power of system
        # Line impedance
        lg      # Line reactance
        rg      # Line resistance
        cg      # Line capacitance
        # Infinite bus voltage
        vg_to_r     #real voltage to-side of line
        vg_to_i     #imaginary voltage to-side of line
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
        M    # Virtual Inertia Constant
        kd   # Active Power PLL Frequency Damping
        kω   # Active Power Frequency Setpoint Damping
        kq   # Reactive Power Droop
        ωf   # Cut-Off frequency Low-Pass Filter (both Active and Reactive)
        # SRF Voltage Control
        kpv     # Voltage control propotional gain
        kiv     # Voltage control integral gain
        kffi    # FeedForward Control for current in voltage control
        # SRF Current Control
        kpc     # Current control propotional gain
        kic     # Current control integral gain
        kffv    # FeedForward Control for voltage in current control
        # Virtual Impedance
        rv
        lv
        # Base Power
        Sinv    # Base Power for Inverter
    end

    MTK.@derivatives d'~t

    # Definition of the states
    states = MTK.@variables begin
        # Grid States
        il_r(t)         # real current flowing in line
        il_i(t)         # imaginary current flowing in line
        vg_from_r(t)    # real voltage from-side of line
        vg_from_i(t)    # imaginary voltage from-side of line
        #vg_to_r(t)      # real voltage to-side of line
        #vg_to_i(t)      # imaginary voltage to-side of line
        # Filter States
        ef_d(t) # d-axis capacitor filter voltage
        ef_q(t) # q-axis capacitor filter voltage
        ic_d(t) # d-axis current flowing into filter
        ic_q(t) # q-axis current flowing into filter
        if_d(t) # d-axis flowing out of the filter
        if_q(t) # q-axis current flowing out of the filter
        # Inner Loop Control States
        ξ_d(t)  #d-axis integrator term for outer AC/DC PI controller
        ξ_q(t)  #q-axis integrator term for outer AC/DC PI controller
        γ_d(t)  #d-axis integrator term for inner AC/DC PI controller
        γ_q(t)  #d-axis integrator term for inner AC/DC PI controller
        # VSM Control States
        θ(t)    # Outer Control Angle
        ω(t)    # Outer Control Frequency
        qf(t)          # Filtered Reactive Power
    end

    # Definition of the variables for non-linear system. Requires https://github.com/SciML/ModelingToolkit.jl/issues/322 to eliminate
    variables = MTK.@variables begin
        # Grid States
        il_r        # real current flowing into grid
        il_i        # imaginary current flowing into grid
        vg_from_r   # real voltage from-side of line
        vg_from_i   # imaginary voltage from-side of line
        #vg_to_r     # real voltage to-side of line
        #vg_to_i     # imaginary voltage to-side of line
        # Filter States
        ef_d # d-axis capacitor filter voltage
        ef_q # q-axis capacitor filter voltage
        ic_d # d-axis current flowing into filter
        ic_q # q-axis current flowing into filter
        if_d # d-axis flowing out of the filter
        if_q # q-axis current flowing out of the filter
        # Inner Loop Control States
        ξ_d         # d-axis integrator term for outer AC/DC PI controller
        ξ_q         # q-axis integrator term for outer AC/DC PI controller
        γ_d         # d-axis integrator term for inner AC/DC PI controller
        γ_q         # d-axis integrator term for inner AC/DC PI controller
        # VSM Control States
        θ           # Outer Control Angle
        ω           # Outer Control Frequency
        qf          # Filtered Reactive Power
    end

    # Expressions
    pm = ef_d * if_d + ef_q * if_q  # AC Active Power into the filter
    qm = -ef_d * if_q + ef_q * if_d # AC Reactive Power into the filter
    #ω_a = ωʳ + Dp * (pʳ - pf)  # Active Power Droop: Not applicable in VSM
    v_hat = vʳ + kq * (qʳ - qf) # Reactive Power Droop
    v_iref_d = v_hat - rv * if_d + ω * lv * if_q # Inner voltage controller d PI
    v_iref_q = -rv * if_q - ω * lv * if_d # Inner voltage controller q PI
    i_hat_d = kpv * (v_iref_d - ef_d) + kiv * ξ_d - ω * cf * ef_q + kffi * if_d # Inner current controller d PI
    i_hat_q = kpv * (v_iref_q - ef_q) + kiv * ξ_q + ω * cf * ef_d + kffi * if_q # Inner current controller q PI
    v_md = kpc * (i_hat_d - ic_d) + kic * γ_d - ω * lf * ic_q + kffv * ef_d
    v_mq = kpc * (i_hat_q - ic_q) + kic * γ_q + ω * lf * ic_d + kffv * ef_q
    p_inv = v_md * ic_d + v_mq * ic_q
    q_inv = -v_md * ic_q + v_mq * ic_d
    if_r = (Sinv / Sb) * (cos(θ) * if_d - sin(θ) * if_q)
    if_i = (Sinv / Sb) * (sin(θ) * if_d + cos(θ) * if_q)
    vg_from_d = cos(θ) * vg_from_r + sin(θ) * vg_from_i
    vg_from_q = -sin(θ) * vg_from_r + cos(θ) * vg_from_i

    model_rhs = [
        # Line Equations
        #𝜕il_r/𝜕t
        (Ωb / lg) * ((vg_from_r - vg_to_r) - (rg * il_r - lg * ω_sys * il_i))
        #𝜕il_i/𝜕t
        (Ωb / lg) * ((vg_from_i - vg_to_i) - (rg * il_i + lg * ω_sys * il_r))
        #𝜕vg_from_r/𝜕t
        (Ωb / cg) * (if_r - il_r) + Ωb * ω_sys * vg_from_i
        ##𝜕vg_from_i/𝜕t
        (Ωb / cg) * (if_i - il_i) - Ωb * ω_sys * vg_from_r
        #Filter Equations
        #𝜕ef_d/𝜕t
        (Ωb / cf) * (ic_d - if_d) + Ωb * ω_sys * ef_q
        #𝜕ef_q/𝜕t
        (Ωb / cf) * (ic_q - if_q) - Ωb * ω_sys * ef_d
        #𝜕ic_d/𝜕t
        (Ωb / lf) * (v_md - ef_d) - Ωb * (rf / lf) * ic_d + Ωb * ω_sys * ic_q
        #𝜕ic_q/𝜕t
        (Ωb / lf) * (v_mq - ef_q) - Ωb * (rf / lf) * ic_q - Ωb * ω_sys * ic_d
        #𝜕if_d/𝜕t
        (Ωb / lt) * (ef_d - vg_from_d) - Ωb * (rt / lt) * if_d + Ωb * ω_sys * if_q
        #𝜕if_q/𝜕t
        (Ωb / lt) * (ef_q - vg_from_q) - Ωb * (rt / lt) * if_q - Ωb * ω_sys * if_d
        ### Inner Control Equations
        #𝜕ξ_d/𝜕t
        v_iref_d - ef_d
        #𝜕ξ_q/𝜕t
        v_iref_q - ef_q
        #𝜕γ_d/𝜕t
        i_hat_d - ic_d
        #𝜕γ_q/𝜕t
        i_hat_q - ic_q
        ### Outer Control Equations
        #𝜕θ/𝜕t
        Ωb * (ω - ω_sys)
        #𝜕ω/𝜕t
        (1 / M) * ((pʳ - pm) + kω * (ωʳ - ω))
        #𝜕qf/𝜕t
        ωf * (qm - qf)
    ]

    # Temporary until SteadyState problems are resolved.
    model_lhs = [
        d(il_r)
        d(il_i)
        d(vg_from_r)
        d(vg_from_i)
        d(ef_d)
        d(ef_q)
        d(ic_d)
        d(ic_q)
        d(if_d)
        d(if_q)
        d(ξ_d)
        d(ξ_q)
        d(γ_d)
        d(γ_q)
        d(θ)
        d(ω)
        d(qf)
    ]

    return model_lhs, model_rhs, states, variables, params
end

function get_ode_system_vsm()
    model_lhs, model_rhs, states, _, params = get_internal_model_vsm(nothing)
    t = params[1]
    _eqs = model_lhs .~ model_rhs
    return MTK.ODESystem(_eqs, t, [states...], [params...][2:end])
end

function get_nonlinear_system_vsm()
    _, model_rhs, _, variables, params = get_internal_model_vsm(nothing)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    return MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
end
