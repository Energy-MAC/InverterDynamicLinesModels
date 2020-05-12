using OrdinaryDiffEq #Gets the solvers
using PowerSystems
#using Plots

include(joinpath(pwd(), "InverterDynamicLinesModels", "InverterDynamicLinesModels.jl"))
# Only need to run this line to re-generate the system data
#include(joinpath(pwd(), "data","make_data.jl"))

# Load Data with PF solution from file
omib_sys = System(joinpath(pwd(), "data", "OMIB_inverter.json"))

# Instantiate analysis objects
parameter_mapping = instantiate_parameters_vsm(omib_sys)
M = instantiate_model_vsm(omib_sys)
u0 = M(parameter_mapping) # works as a test, not really necessary to call
jac = instantiate_jacobian_vsm(M)
ss = instantiate_small_signal(M, jac)
ss(M, jac)


# Test of parameter sweep for the gain of the integral gain of voltage
println("$(parameter_mapping[24])")
param_space = range(0.1,10000, length=5000)
res = Vector{Number}(undef, 5000)
for (i, val) in enumerate(param_space)
    parameter_values[24] = val
    M(parameter_values)
    ss(M, jac)
    res[i] = ss.eigen_vals[end]
end
plot(param_space, res)

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
