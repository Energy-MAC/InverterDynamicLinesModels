function max_eigenvalue(M::ModelOperatingPoint, J::ModelJacobian)
    eigen_vals, eigen_vect = LinearAlgebra.eigen(J(M))
    return maximum(real.(eigen_vals))
end
