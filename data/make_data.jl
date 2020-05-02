using PowerSystems
using NLsolve
const PSY = PowerSystems

omib_file_dir= joinpath(pwd(), "data/OMIB.raw")
omib_sys = System(PowerModelsData(omib_file_dir), time_series_in_memory = true, runchecks=false)
slack_bus = get_components_by_name(Component, omib_sys, "BUS 1")[1]
inf_source = Source(
    name = "InfBus",
    available = true,
    activepower = -0.,
    reactivepower = 0.0,
    bus = slack_bus,
    X_th = 0.000005
)
PSY.add_component!(omib_sys, inf_source)

inv_bus = get_components_by_name(Component, omib_sys, "BUS 2")[1]
inv_bus.bustype = BusTypes.PQ
battery = GenericBattery(
    name = "Battery",
    primemover = PrimeMovers.BA,
    available = true,
    bus = inv_bus,
    energy = 5.0,
    capacity = (min = 5.0, max = 100.0),
    rating = 0.0275, #Value in per_unit of the system
    activepower = 0.01375,
    inputactivepowerlimits = (min = 0.0, max = 50.0),
    outputactivepowerlimits = (min = 0.0, max = 50.0),
    reactivepower = 0.0,
    reactivepowerlimits = (min = -50.0, max = 50.0),
    efficiency = (in = 0.80, out = 0.90),
)
add_component!(omib_sys, battery)
res = solve_powerflow!(omib_sys, nlsolve)

###### Converter Data ######
converter() = AverageConverter(
    v_rated = 690.0,
    s_rated = 2.75,
)
###### DC Source Data ######
dc_source() = FixedDCSource(voltage = 600.0) #Not in the original data, guessed.

###### Filter Data ######
filter() = LCLFilter(
    lf = 0.08,
    rf = 0.003,
    cf = 0.074,
    lg = 0.2,
    rg = 0.01,
)

###### PLL Data ######
pll() = KauraPLL(
    ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
    kp_pll = 0.084,  #PLL proportional gain
    ki_pll = 4.69,   #PLL integral gain
)

###### Outer Control ######
function outer_control()
    function virtual_inertia()
        return VirtualInertia(
        Ta = 2.0,
        kd = 400.0,
        kω = 20.0,
        ωb = 2 * pi * 50.0,
    )
    end
    function reactive_droop()
        return ReactivePowerDroop(
        kq = 0.2,
        ωf = 1000.0,
    )
    end
    return OuterControl(virtual_inertia(), reactive_droop())
end

######## Inner Control ######
inner_control() = CurrentControl(
    kpv = 0.59,     #Voltage controller proportional gain
    kiv = 736.0,    #Voltage controller integral gain
    kffv = 0.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
    rv = 0.0,       #Virtual resistance in pu
    lv = 0.2,       #Virtual inductance in pu
    kpc = 1.27,     #Current controller proportional gain
    kic = 14.3,     #Current controller integral gain
    kffi = 0.0,     #Binary variable enabling the current feed-forward in output of current controllers
    ωad = 50.0,     #Active damping low pass filter cut-off frequency
    kad = 0.2,
)

inverter = PSY.DynamicInverter(
        number = 1,
        name = "VSM",
        bus = get_bus(battery),
        ω_ref = 1.0,
        V_ref = get_voltage(get_bus(battery)),
        P_ref = get_activepower(battery),
        Q_ref = get_reactivepower(battery),
        MVABase = get_rating(battery),
        converter = converter(),
        outer_control = outer_control(),
        inner_control = inner_control(),
        dc_source = dc_source(),
        freq_estimator = pll(),
        filter = filter(),
    )

add_component!(omib_sys, inverter)
to_json(omib_sys, joinpath(pwd(), "data/OMIB_inverter.json"))
