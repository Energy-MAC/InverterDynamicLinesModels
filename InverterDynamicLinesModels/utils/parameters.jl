
function instantiate_parameters(model) #, system::PSY.System)
    # TODO: Automate better with PSY getter functions
    # AC side quantities
    ωb,      # Base Frequency
    # Grid impadance
    lg,      # Grid reactance
    rg,      # Grid resistance
    #Reference set-point input
    pʳ,      # Active Power Setpoint
    qʳ,      # Reactive Power Setpoint
    vʳ,      # Voltage Setpoint
    ωʳ,      # Frequency Setpoint
    # Load at rated voltage
    vl,      # Load rated voltage
    pl,      # Load rated power
    # Filter parameters
    lf,      # Filter reactance
    cf,      # Filter capacitance
    rf,      # Filter Resistance
    ωz,      # Filtering frequency
    # Transformer Parameters
    rt,      # Transformer resistance
    lt,      # Transformer reactance
    # OuterControl Loops
    Dp,      # Active Power Droop
    Dq,      # Reactive Power Droop
    # SRF Current Control
    kip,     # Current control propotional gain
    kii,     # Current control integral gain
    kffi,    # Current control differential gain
    # SRF Voltage Control
    kvp,     # Voltage control propotional gain
    kvi,     # Voltage control integral gain
    kffv,    # Voltage control differential gain
    # Virtual Impedance
    rv,
    lv,
    # DC Source Parameters
    leq,
    req,
    vb,
    cdc,
    # DC/DC converter controller parameters
    vdcʳ,    # DC Voltage reference
    kivb,
    kpvb,
    kpib,
    kiib,
    Ts = MTK.parameters(model) # Sampling time = parameters(model)

    Ub = 690 # Get using PSY
    fb = 60 # Get using PSY.
    _ωb = 2 * pi * fb
    Sb = 0.8e6 # Get using PSY
    Vb = 690 # Get using PSY
    # System base for AC side
    Ib = Sb / Ub
    Zb = Vb / Ib
    Lb = Zb / _ωb
    Cb = 1 / (Zb * _ωb)
    # System base for DC side
    Vb_dc = 2 * Vb
    Ib_dc = Sb / Vb_dc
    Zb_dc = Vb_dc / Ib_dc
    Lb_dc = Zb_dc / _ωb
    Cb_dc = 1 / (Zb_dc * _ωb)

    p = [
        ωb => _ωb
        # Grid impadance (Currently not used)
        lg => 0.2 # Get using PSY
        rg => 0.01 # Get using PSY
        #Reference set-point inputs
        pʳ => 0.5 # Get using PSY
        qʳ => 0 # Get using PSY
        vʳ => 1 # Get using PSY
        ωʳ => 1
        # TODO: what is 2.02 ?
        # Load at rated voltage
        vl => 1 # Get using PSY
        pl => 0.5 # Get using PSY
        # Filter parameters
        lf => 0.08  # Get using PSY
        cf => 0.074 # Get using PSY
        rf => 0.003 # Get using PSY
        # Filtering frequency
        ωz => 0.1 * _ωb
        # Transformer Parameters
        rt => 0.01 # Get using PSY
        lt => 0.2  # Get using PSY
        # Outer Control Loops
        # Active Power Droop Control
        Dp => 0.02 # Get using PSY
        # Reactive Power Droop Control
        Dq => 0.001 # Get using PSY
        # Inner Control Loops
        # SRF Current Control
        kip => 1.27 # Get using PSY
        kii => 14.3 # Get using PSY
        kffi => 0.0 # Get using PSY
        # SRF Voltage Control
        kvp => 0.59 # Get using PSY
        kvi => 736  # Get using PSY
        kffv => 1.0 # Get using PSY
        # Virtual Impedance
        rv => 0     # Get using PSY
        lv => 0.2   # Get using PSY
        # DC source Parameters
        leq => 5e-3 / Lb_dc # Get using PSY
        req => 0.0 / Zb_dc  # Get using PSY
        vb => 1000 / Vb_dc  # Get using PSY
        vdcʳ => 2.02 * Ub / Vb_dc # Get using PSY
        cdc => 100e-6 / Cb_dc # Get using PSY
        # DC/DC converter controller parameters
        kpvb => 0.6 # Get using PSY
        kivb => 4   # Get using PSY
        kpib => 0.3863 # Get using PSY
        kiib => 10.34 # Get using PSY
        Ts => 1 / (3.2e3) # Get using PSY
    ]
    return p
end
