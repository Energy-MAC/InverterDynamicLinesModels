function Base.show(io::IO, m::MIME"text/plain", ::ModelJacobian)
    println(io, "Jacobian")
end

function Base.show(io::IO, m::MIME"text/plain", ::ModelOperatingPoint)
    println(io, "ModelOperatingPoint")
end
