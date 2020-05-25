using PowerSystems
using Plots

include(joinpath(pwd(), "InverterDynamicLinesModels", "InverterDynamicLinesModels.jl"))
# Only need to run this line to re-generate the system data
#include(joinpath(pwd(), "data","make_data.jl"))

# Load Data with PF solution from file
omib_sys = System(joinpath(pwd(), "data", "OMIB_inverter.json"))

#Instantiate Parameters for VSM model
parameter_mapping = instantiate_parameters(VInertia, omib_sys)
#Instantiate Model for VSM + Dynamic Lines
M_vsm = instantiate_model(VInertia, DynamicLines, omib_sys)
u0 = M_vsm(parameter_mapping) # works as a test, not really necessary to call
#Instantiate Jacobian
jac_exp_vsm = get_jacobian_function(VInertia, DynamicLines);
fjac_vsm = eval(jac_exp_vsm);
jac_vsm = instantiate_jacobian(M_vsm, fjac_vsm)
#Instantiate Small Signal Analysis
ss_vsm = instantiate_small_signal(M_vsm, jac_vsm)
ss_vsm(M_vsm, jac_vsm)

# Create Parameter Vector to modify it later
parameter_values = [x.second for x in parameter_mapping];
# Update p_ref to 1.0
parameter_values[10] = 1.0;

#Initiate Sweep over kp and kq
N = 75
res = Matrix{Complex}(undef, N, 5*N) #Max Eigenvalue
damp = Matrix{Float64}(undef, N, 5*N) #Damping Ratio of Max Eigenvalue
kq_crit = Vector{Float64}(undef, N) #Critical kq for each kq
p_fact = Array{Float64}(undef, (N, 5*N, 19)) #Participation Factors for each kp, kq
damp_crit = Vector{Float64}(undef, N) #Damping Ratio for each kp at kq_crit
p_fact_crit = Matrix{Float64}(undef, N, 19) #Participation Factors for each kp at kq_crit
eig_crit = Vector{Complex}(undef, N) #Eigenvalue at kq_crit for each kp
kp_range = range(0.001, 0.42, length = N) #kp_range
kq_range = range(0.001, 10, length = 5*N) #kq_range

#Initiate Sweep
for (i, kpval) in enumerate(kp_range)
    parameter_values[21] = kpval
    for (j, kqval) in enumerate(kq_range)
        parameter_values[22] = kqval
        M_vsm(parameter_values)
        jac_vsm(M_vsm)
        ss_vsm(M_vsm, jac_vsm)
        eig = ss_vsm.eigen_vals[end]
        res[i, j] = eig
        damp[i, j] = round(ss_vsm.damping_vector[end], digits = 4)
        p_fact[i, j, :] = round.(ss_vsm.participation_factors[:, end], digits = 4)
    end
    idx = findall(x -> real(x) >= 0, res[i, :]);
    if isempty(idx)
        kq_crit[i] = -1
        damp_crit[i] = damp[i, 1]
        p_fact_crit[i, :] = p_fact[i, 1, :]
    else
        idx_crit = idx[1]
        eig_crit[i] = res[i, idx_crit]
        kq_crit[i] = round(kq_range[idx_crit], digits = 4)
        damp_crit[i] = damp[i, idx_crit]
        p_fact_crit[i, :] = p_fact[i, idx_crit, :]
    end
end 


#Instantiate Parameters for VSM model + Static Lines
parameter_mapping = instantiate_parameters(VInertia, omib_sys)
M_vsm = instantiate_model(VInertia, StaticLines, omib_sys)
u0 = M_vsm(parameter_mapping) # works as a test, not really necessary to call
jac_exp_vsm = get_jacobian_function(VInertia, StaticLines);
fjac_vsm = eval(jac_exp_vsm);
jac_vsm = instantiate_jacobian(M_vsm, fjac_vsm)
ss_vsm = instantiate_small_signal(M_vsm, jac_vsm)
ss_vsm(M_vsm, jac_vsm)

# Create Parameter Vector to modify it later
parameter_values = [x.second for x in parameter_mapping];
parameter_values[10] = 1.0;

#Initiate Sweep over kp and kq
N = 75
res = Matrix{Complex}(undef, N, 5*N)
damp_alg = Matrix{Float64}(undef, N, 5*N)
p_fact_alg = Array{Float64}(undef, (N, 5*N, 13))
damp_crit_alg = Vector{Float64}(undef, N)
p_fact_crit_alg = Matrix{Float64}(undef, N, 13)
eig_crit_alg = Vector{Complex}(undef, N)
kq_crit_alg = Vector{Float64}(undef, N)
kp_range = range(0.001, 0.42, length = N)
kq_range = range(0.001, 20, length = 5*N)
for (i, kpval) in enumerate(kp_range)
    parameter_values[21] = kpval
    for (j, kqval) in enumerate(kq_range)
        parameter_values[22] = kqval
        M_vsm(parameter_values)
        jac_vsm(M_vsm)
        ss_vsm(M_vsm, jac_vsm)
        eig = ss_vsm.eigen_vals[end]
        res[i, j] = eig     
        damp_alg[i, j] = round(ss_vsm.damping_vector[end], digits = 4)
        p_fact_alg[i, j, :] = round.(ss_vsm.participation_factors[:, end], digits = 4)
    end
    idx = findall(x -> real(x) >= 0, res[i, :]);
    if isempty(idx)
        kq_crit_alg[i] = -1
        damp_crit_alg[i] = damp_alg[i, 2] #Store for kq = 0.02
        p_fact_crit_alg[i, :] = p_fact_alg[i, 2, :]
    else
        idx_crit = idx[1]
        eig_crit_alg[i] = res[i, idx_crit]
        kq_crit_alg[i] = round(kq_range[idx_crit], digits = 4)
        damp_crit_alg[i] = damp_alg[i, idx_crit]
        p_fact_crit_alg[i, :] = p_fact_alg[i, idx_crit, :]
    end
end 

#Instantiate Parameters for VSM model + Full AC side
parameter_mapping = instantiate_parameters(VInertia, omib_sys)
M_vsm = instantiate_model(VInertia, ACStatic, omib_sys)
u0 = M_vsm(parameter_mapping) # works as a test, not really necessary to call
jac_exp_vsm = get_jacobian_function(VInertia, ACStatic);
fjac_vsm = eval(jac_exp_vsm);
jac_vsm = instantiate_jacobian(M_vsm, fjac_vsm)
ss_vsm = instantiate_small_signal(M_vsm, jac_vsm)
ss_vsm(M_vsm, jac_vsm)

# Create Parameter Vector to modify it later
parameter_values = [x.second for x in parameter_mapping];
parameter_values[10] = 1.0;

#Initiate Sweep over kp and kq
N = 75
res = Matrix{Complex}(undef, N, 7*N)
damp_st = Matrix{Float64}(undef, N, 7*N)
p_fact_st = Array{Float64}(undef, (N, 7*N, 7))
damp_crit_st = Vector{Float64}(undef, N)
p_fact_crit_st = Matrix{Float64}(undef, N, 7)
eig_crit_st = Vector{Complex}(undef, N)
kq_crit_st = Vector{Float64}(undef, N)
kp_range = range(0.001, 0.42, length = N)
kq_range = range(0.001, 20, length = 7*N)
for (i, kpval) in enumerate(kp_range)
    parameter_values[21] = kpval
    for (j, kqval) in enumerate(kq_range)
        parameter_values[22] = kqval
        M_vsm(parameter_values)
        jac_vsm(M_vsm)
        ss_vsm(M_vsm, jac_vsm)
        eig = ss_vsm.eigen_vals[end]
        res[i, j] = eig
        damp_st[i, j] = round(ss_vsm.damping_vector[end], digits = 4)
        p_fact_st[i, j, :] = round.(ss_vsm.participation_factors[:, end], digits = 4)
    end
    idx = findall(x -> real(x) >= 0, res[i, :]);
    if isempty(idx)
        kq_crit_st[i] = 25 #Store at 25 if no kq_crit is find
        damp_crit_st[i] = damp_st[i, 2] #Store for kq = 0.02
        p_fact_crit_st[i, :] = p_fact_st[i, 2, :]
    else
        idx_crit = idx[1]
        eig_crit_st[i] = res[i, idx_crit]
        kq_crit_st[i] = round(kq_range[idx_crit], digits = 4)
        damp_crit_st[i] = damp_st[i, idx_crit]
        p_fact_crit_st[i, :] = p_fact_st[i, idx_crit, :]
    end
end 

plot(kp_range, kq_crit)
plot!(kp_range, kq_crit_alg)
plot!(kp_range, kq_crit_st)