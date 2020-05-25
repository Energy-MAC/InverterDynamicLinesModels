using OrdinaryDiffEq #Gets the solvers
using PowerSystems
using Plots
plotlyjs()

include(joinpath(pwd(), "InverterDynamicLinesModels", "InverterDynamicLinesModels.jl"))
# Only need to run this line to re-generate the system data
#include(joinpath(pwd(), "data","make_data.jl"))

# Load Data with PF solution from file
omib_sys = System(joinpath(pwd(), "data", "OMIB_inverter.json"))

#### VSM model with Full Lines ####
#Instantiate VSM Model with DynamicLines
M_vsm = instantiate_model(VInertia, DynamicLines, omib_sys)
#Instantiate Jacobian
jac_exp_vsm = get_jacobian_function(VInertia, DynamicLines);
fjac_vsm = eval(jac_exp_vsm);
jac_vsm = instantiate_jacobian(M_vsm, fjac_vsm)
#Instantiate Small Signal Object
ss_vsm = instantiate_small_signal(M_vsm, jac_vsm)
#Update Small Signal Object
ss_vsm(M_vsm, jac_vsm)
#Report Eigenvalues
ss_vsm.eigen_vals

#Run ODE problem
ode_vsm1 = instantiate_ode(M_vsm, CircuitTrip(time = 1.0); tspan = (0.0, 5))
sol_vsm1 = solve(ode_vsm1, GRK4T())
plot(sol_vsm1)

ode_vsm = instantiate_ode(M_vsm, PowerOutputIncrease(time = 1.0, new_value = 0.65); tspan = (0.0, 5))
sol_vsm = solve(ode_vsm, GRK4T())
plot(sol_vsm)

#### VSM model with Algebraic Lines ####
M_vsm_slines = instantiate_model(VInertia, StaticLines, omib_sys)
jac_vsm_slines = instantiate_jacobian(M_vsm_slines)
ss_vsm_slines = instantiate_small_signal(M_vsm_slines, jac_vsm_slines)
ss_vsm_slines(M_vsm_slines, jac_vsm_slines)

# Report Eigenvalues
ss_vsm_slines.eigen_vals

#Run ODE problem
ode_vsm_slines = instantiate_ode(M_vsm_slines, CircuitTrip(time = 1.0); tspan = (0.0, 5))
sol_vsm_slines = solve(ode_vsm_slines, GRK4T())
plot(sol_vsm_slines)

#### VSM Model with Filter+Lines Algebraic ####
M_vsm_static = instantiate_model(VInertia, ACStatic, omib_sys)
jac_vsm_static = instantiate_jacobian(M_vsm_static)
ss_vsm_static = instantiate_small_signal(M_vsm_static, jac_vsm_static)
ss_vsm_static(M_vsm_static, jac_vsm_static)

#Report Eigenvalues
ss_vsm_static.eigen_vals

#Run ODE problem
ode_vsm_static = instantiate_ode(M_vsm_static; tspan = (0.0, 5))
sol_vsm_static = solve(ode_vsm_static, GRK4T())
plot(sol_vsm_static)

#= Droop Model not functional
#### Droop Model with Full Dynamic Lines ####
M_droop = instantiate_model(DroopModel, DynamicLines, omib_sys)
jac_droop = instantiate_jacobian(M_droop)
ss_droop = instantiate_small_signal(M_droop, jac_droop)
ss_droop(M_droop, jac_droop)

#Report Eigenvalues
ss_droop.eigen_vals

#Run ODE problem
ode_droop = instantiate_ode(M_droop; tspan = (0.0, 5))
sol_droop = solve(ode_droop, GRK4T())
plot(sol_droop)
=#
