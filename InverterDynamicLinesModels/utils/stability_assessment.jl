struct SmallSignal
    eigen_vals::Vector{Complex{Float64}}
    R_eigen_vect::Matrix{Complex{Float64}}
    L_eigen_vect::Matrix{Complex{Float64}}
    participation_factors::Matrix{Float64}
    damping_vector::Vector{Float64}
    max_eigenvalue::Dict{Symbol, Float64}
end

function instantiate_small_signal(M::ModelOperatingPoint, J::ModelJacobian)
    eigen_vals, R_eigen_vect = LinearAlgebra.eigen(J(M))
    L_eigen_vect = inv(R_eigen_vect)
    damping_vector = Vector{Float64}(undef, length(eigen_vals))
    max_eigenvalue = Dict(:val => -Inf, :loc => -1)
    participation_factors = zeros(size(L_eigen_vect))
    for (ix, eigen_val) in enumerate(eigen_vals)
        den = sum(abs.(L_eigen_vect[ix, :]).*abs.(R_eigen_vect[:,ix]))
        participation_factors[ix,:] = abs.(L_eigen_vect[ix, :]).*abs.(R_eigen_vect[:,ix])./den
        current_eigenvalue = real(eigen_val)
        damping_vector[ix] = -1*real(eigen_val) / sqrt(real(eigen_val)^2 + imag(eigen_val)^2)
        if max_eigenvalue[:val] < current_eigenvalue
           max_eigenvalue[:val] = current_eigenvalue
           max_eigenvalue[:loc] = ix
        end
    end
    return SmallSignal(eigen_vals,
                     R_eigen_vect,
                     L_eigen_vect,
                     participation_factors,
                     damping_vector,
                     max_eigenvalue)
end

function (S::SmallSignal)(M::ModelOperatingPoint, J::ModelJacobian)
    _J = J(M)
    S.eigen_vals .= LinearAlgebra.eigvals(_J)
    S.L_eigen_vect .= LinearAlgebra.eigvecs(_J)
    S.L_eigen_vect .= LinearAlgebra.inv(S.R_eigen_vect)
    for (ix, eigen_val) in enumerate(S.eigen_vals)
        den = sum(abs.(S.L_eigen_vect[ix, :]).*abs.(S.R_eigen_vect[:,ix]))
        S.participation_factors[ix,:] .= abs.(S.L_eigen_vect[ix, :]).*abs.(S.R_eigen_vect[:,ix])./den
        current_eigenvalue = real(eigen_val)
        S.damping_vector[ix] = -1*real(eigen_val) / sqrt(real(eigen_val)^2 + imag(eigen_val)^2)
        if S.max_eigenvalue[:val] < current_eigenvalue
           S.max_eigenvalue[:val] = current_eigenvalue
           S.max_eigenvalue[:loc] = ix
        end
    end
