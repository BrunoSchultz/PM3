struct Cosmology{T <: AbstractFloat}
    Ω_m0::T
    Ω_Λ0::T
    Ω_k0::T
    A_init::T
    H0::T
end

# Defaults to Float64, but you can call ΛCDM(Float32) for GPU runs
function ΛCDM(::Type{T}=Float64) where {T <: AbstractFloat}
    return Cosmology{T}(
        T(0.315),   
        T(0.685),   
        T(0.0),     
        T(0.01),    
        T(67.4)     
    )
end

function Hubble(a::Float64, C::Cosmology)
    H = sqrt(C.H0^2 * (C.Ω_m0 / a^3 + C.Ω_k0 / a^2 + C.Ω_Λ0))
    return H
end

function f(a::Float64, C::Cosmology)
    f = (a^(-1) * (C.Ω_m0 + C.Ω_k0 * a +C.Ω_Λ0 * a^3))^(-0.5)
    return f
end

function Dt(C::Cosmology)
    Dt = 5 / 2 * C.Ω_m0 * (C.Ω_m0^(4 / 7) - C.Ω_Λ0 + (1 + C.Ω_m0 / 2)*(1 + C.Ω_Λ0 / 70))^(-1)
end
