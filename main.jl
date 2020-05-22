#using OrdinaryDiffEq #Gets the solvers
using PowerSystems
#using Plots

include(joinpath(pwd(), "InverterDynamicLinesModels", "InverterDynamicLinesModels.jl"))
# Only need to run this line to re-generate the system data
#include(joinpath(pwd(), "data","make_data.jl"))

# Load Data with PF solution from file
omib_sys = System(joinpath(pwd(), "data", "OMIB_inverter.json"))

#### VSM model with Full Lines ####

#Instantiate Parameters for VSM model
parameter_mapping = instantiate_parameters(VInertia, omib_sys)
#Instantiate VSM Model with DynamicLines
M_vsm = instantiate_model(VInertia, DynamicLines, omib_sys)
#Instantiate x0 for parameters
u0 = M_vsm(parameter_mapping) # works as a test, not really necessary to call
#Instantiate Jacobian
jac_vsm = instantiate_jacobian(M_vsm)
#Instantiate Small Signal Object
ss_vsm = instantiate_small_signal(M_vsm, jac_vsm)
#Update Small Signal Object
ss_vsm(M_vsm, jac_vsm)
#Report Eigenvalues
ss_vsm.eigen_vals

#### VSM model with Algebraic Lines ####

parameter_mapping = instantiate_parameters(VInertia, omib_sys)
M_vsm_slines = instantiate_model(VInertia, StaticLines, omib_sys)
u0 = M_vsm_slines(parameter_mapping) # works as a test, not really necessary to call
jac_vsm_slines = instantiate_jacobian(M_vsm_slines)
ss_vsm_slines = instantiate_small_signal(M_vsm_slines, jac_vsm_slines)
ss_vsm_slines(M_vsm_slines, jac_vsm_slines)

# Report Eigenvalues
ss_vsm_slines.eigen_vals

#### VSM Model with Filter+Lines Algebraic ####

parameter_mapping = instantiate_parameters(VInertia, omib_sys)
M_vsm_static = instantiate_model(VInertia, ACStatic, omib_sys)
u0 = M_vsm_static(parameter_mapping) # works as a test, not really necessary to call
jac_vsm_static = instantiate_jacobian(M_vsm_static)
ss_vsm_static = instantiate_small_signal(M_vsm_static, jac_vsm_static)
ss_vsm_static(M_vsm_static, jac_vsm_static)

#Report Eigenvalues
ss_vsm_static.eigen_vals

#### Droop Model with Full Dynamic Lines ####
parameter_mapping = instantiate_parameters(DroopModel, omib_sys)
M_droop = instantiate_model(DroopModel, DynamicLines, omib_sys)
u0 = M_droop(parameter_mapping) # works as a test, not really necessary to call
jac_droop = instantiate_jacobian(M_droop)
ss_droop = instantiate_small_signal(M_droop, jac_droop)
ss_droop(M_droop, jac_droop)

#Report Eigenvalues
ss_droop.eigen_vals

# Returns Generic ODE system and solves
#ode_prob = instantiate_ode(omib_sys; tspan = (0.0, 5))
#ode_prob = instantiate_ode(M; tspan = (0.0, 5))
#sol1 = solve(ode_prob, Rosenbrock23())
#plot(sol1, vars = (0, 13))

#=

parameters.pl = 0.6;

tspan = (0.0,1)
prob = ODEProblem(ode_system!,sol1.u[end],tspan, parameters)
sol2 = solve(prob)

plot(sol2,vars=(0,13),title = "DC Voltage After Load Step")

function condition(u,t,integrator)
    t == 0.2
end

function affect!(integrator)
  parameters.pl = 0.6;
end
cb = DiscreteCallback(condition, affect!)

parameters = get_params(Ub,fb,ωb,Sb,Vb)
const tstop = [0.2]

parameters.pl

prob = ODEProblem(ode_system!,ic.zero,tspan, parameters)
sol3 = solve(prob,Tsit5(),callback = cb, tstops=tstop)

plot(sol1,vars=(0,13),title = "DC Voltage Before Load Step")

parameters.pl
=#
