function instantiate_parameters(::Type{VInertia}, system::PSY.System)
    # TODO: Automate better with PSY getter functions
    # AC side quantities
    # AC side quantities
    Ωb,      # Base Frequency
    ω_sys,   # Grid frequency in pu
    Sb,      # Base Power of system
    # Line impedance
    lg,      # Line reactance
    rg,      # Line resistance
    cg,      # Line capacitance
    # Infinite bus voltage
    vg_to_r,     #real voltage to-side of line
    vg_to_i,     #imaginary voltage to-side of line
    #Reference set-point input
    pʳ,      # Active Power Setpoint
    qʳ,      # Reactive Power Setpoint
    vʳ,      # Voltage Setpoint
    ωʳ,      # Frequency Setpoint
    # Filter parameters
    lf,      # Filter reactance
    cf,      # Filter capacitance
    rf,      # Filter Resistance
    rt,      # Transformer resistance
    lt,      # Transformer reactance
    # OuterControl Loops
    M,    # Virtual Inertia Constant
    kd,   # Active Power PLL Frequency Damping
    kp,   # Active Power Frequency Setpoint Damping (1 / kω)
    kq,   # Reactive Power Droop
    ωf,   # Cut-Off frequency Low-Pass Filter (both Active and Reactive)
    # SRF Voltage Control
    kpv,     # Voltage control propotional gain
    kiv,     # Voltage control integral gain
    kffi,    # FeedForward Control for current in voltage control
    # SRF Current Control
    kpc,     # Current control propotional gain
    kic,     # Current control integral gain
    kffv,   # FeedForward Control for voltage in current control
    # Virtual Impedance
    rv,
    lv,
    # Base Power
    Sinv = MTK.parameters(get_nonlinear_system(VInertia, DynamicLines)) # Sampling time = parameters(model)

    fb = PowerSystems.get_frequency(system) # Get using PSY. Not updating json file.
    _Ωb = 2 * pi * fb
    _Sb = 100.0 # Get using PSY
    _Sinv = 2.75 # Get using PSY

    p = [
        Ωb => _Ωb # Get using PSY
        ω_sys => 1.0
        Sb => _Sb # Get using PSY
        # Line impedance: Kiwi Line 42km
        lg => 0.0025289*42
        rg => 0.0002429*42
        cg => 0.0006648*42
        # Infinite bus voltage
        vg_to_r => 1.0
        vg_to_i => 0.0
        #Reference set-point inputs
        pʳ => 0.5 # Get using PSY
        qʳ => 0.0 # Get using PSY
        vʳ => 1.02 # Get using PSY
        ωʳ => 1.0
        # Filter parameters
        lf => 0.08  # Get using PSY
        cf => 0.074 # Get using PSY
        rf => 0.003 # Get using PSY
        # Transformer Parameters
        rt => 0.2 # Get using PSY
        lt => 0.075  # Get using PSY
        # Outer Control Loops
        # Active Power Droop Control
        M => 0.05 #2.0
        kd => 400.0 # Get using PSY
        kp => 0.02 #20.0
        kq => 0.02
        ωf => 1000.0
        # Inner Control Loops
        # SRF Voltage Control
        kpv => 0.59 # Get using PSY
        kiv => 736.0 # Get using PSY
        kffi => 0.0 # Get using PSY
        # SRF Current Control
        kpc => 1.27 # Get using PSY
        kic => 14.3  # Get using PSY
        kffv => 0.0 # Get using PSY
        # Virtual Impedance
        rv => 0     # Get using PSY
        lv => 0.2   # Get using PSY
        # Base Power
        Sinv => _Sinv * 10
    ]
    return p
end

function instantiate_parameters(::Type{DroopModel}, system::PSY.System)
    # TODO: Automate better with PSY getter functions
    # AC side quantities
    # AC side quantities
    Ωb,      # Base Frequency
    ω_sys,   # Grid frequency in pu
    Sb,      # Base Power of system
    # Line impedance
    lg,      # Line reactance
    rg,      # Line resistance
    cg,      # Line capacitance
    # Infinite bus voltage
    vg_to_r,     #real voltage to-side of line
    vg_to_i,     #imaginary voltage to-side of line
    #Reference set-point input
    pʳ,      # Active Power Setpoint
    qʳ,      # Reactive Power Setpoint
    vʳ,      # Voltage Setpoint
    ωʳ,      # Frequency Setpoint
    # Filter parameters
    lf,      # Filter reactance
    cf,      # Filter capacitance
    rf,      # Filter Resistance
    rt,      # Transformer resistance
    lt,      # Transformer reactance
    # OuterControl Loops
    kp,   # Active Power Droop
    kq,   # Reactive Power Droop
    kα,   # Frequency Power Droop for Reactive Power
    kβ,   # Voltage Power Droop for Active Power
    ωf,   # Cut-Off frequency Low-Pass Filter (both Active and Reactive)
    # SRF Voltage Control
    kpv,     # Voltage control propotional gain
    kiv,     # Voltage control integral gain
    kffi,    # FeedForward Control for current in voltage control
    # SRF Current Control
    kpc,     # Current control propotional gain
    kic,     # Current control integral gain
    kffv,   # FeedForward Control for voltage in current control
    # Virtual Impedance
    rv,
    lv,
    # Base Power
    Sinv = MTK.parameters(get_nonlinear_system(DroopModel,  DynamicLines)) # Sampling time = parameters(model)

    fb = PowerSystems.get_frequency(system) # Get using PSY. Not updating json file.
    _Ωb = 2 * pi * fb
    _Sb = 100.0 # Get using PSY
    _Sinv = 2.75 # Get using PSY

    p = [
        Ωb => _Ωb # Get using PSY
        ω_sys => 1.0
        Sb => _Sb # Get using PSY
        # Line impedance: Kiwi Line 42km
        lg => 0.0025289*42
        rg => 0.0002429*42
        cg => 0.0006648*42
        # Infinite bus voltage
        vg_to_r => 1.0
        vg_to_i => 0.0
        #Reference set-point inputs
        pʳ => 0.5 # Get using PSY
        qʳ => 0.0 # Get using PSY
        vʳ => 1.02 # Get using PSY
        ωʳ => 1.0
        # Filter parameters
        lf => 0.08  # Get using PSY
        cf => 0.074 # Get using PSY
        rf => 0.003 # Get using PSY
        # Transformer Parameters
        rt => 0.2 # Get using PSY
        lt => 0.075  # Get using PSY
        # Outer Control Loops
        # Active Power Droop control
        kp => 0.02
        kq => 0.02
        kα => 0.005
        kβ => 0.005
        ωf => 1000.0
        # Inner Control Loops
        # SRF Voltage Control
        kpv => 0.59 # Get using PSY
        kiv => 736.0 # Get using PSY
        kffi => 0.0 # Get using PSY
        # SRF Current Control
        kpc => 1.27 # Get using PSY
        kic => 14.3  # Get using PSY
        kffv => 0.0 # Get using PSY
        # Virtual Impedance
        rv => 0     # Get using PSY
        lv => 0.2   # Get using PSY
        # Base Power
        Sinv => _Sinv * 10
    ]
    return p
end
