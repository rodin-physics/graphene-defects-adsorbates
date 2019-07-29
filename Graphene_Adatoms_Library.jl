using QuadGK
using LinearAlgebra

## Parameters
t = 2.8;        # Hopping integral
ν = 1e-3;       # Small number for relative tolerance
η = 1e-7;       # Small number for imaginary part
NumEvals = 1e4; # Maximum number of evaluations in integrals

# Graphene lattice vectors in Angstroms
d1 = 2.46 .* [1 √(3)] ./ 2;
d2 = 2.46 .* [-1 √(3)] ./ 2;

# Sublattice spinors
A = [1, 0]
B = [0, 1]

# Graphene Coordinate type
struct GrapheneCoord
    u :: Int
    v :: Int
    sublattice :: Array{Int64, 1}
end

## Helper functions
# Auxiliary W function
function W(z :: ComplexF64, x :: Float64)
    return ((z / t)^2 - 1) / (4 * cos(x)) - cos(x)
end

# Auxiliary Y function used in the calculation of Ω
function Y(n :: Int, z :: ComplexF64, x :: Float64)
    W_ = W(z, x)
    return (W_ - √(W_ - 1) * √(W_ + 1))^n / (√(W_ - 1) * √(W_ + 1))
end

## Main Functions
function Ω(z :: ComplexF64, u :: Int, v :: Int) :: ComplexF64
    f_int(x) = exp(1im * (u - v) * x) / cos(x) * Y(abs.(u + v), z, x)
    res = quadgk(f_int, 0, 2 * π, rtol = ν, maxevals = NumEvals)
    return (res[1] / (8 * π * t^2))
end

# The propagator function picks out the correct element of the Ξ matrix based
# on the sublattices of the graphene coordinates
function Propagator(Imp_l :: GrapheneCoord, Imp_m :: GrapheneCoord, z) :: ComplexF64
    u = Imp_l.u - Imp_m.u
    v = Imp_l.v - Imp_m.v
    if Imp_l.sublattice == Imp_m.sublattice
        return (z * Ω(z, u, v))
    elseif (Imp_l.sublattice == A && Imp_m.sublattice ==  B)
        return (- t * (Ω(z, u, v) + Ω(z, u + 1, v) + Ω(z, u, v + 1)))
    elseif (Imp_l.sublattice == B && Imp_m.sublattice ==  A)
        return (- t * (Ω(z, u, v) + Ω(z, u - 1, v) + Ω(z, u, v - 1)))
    else
        error("Illegal sublattice parameter")
    end
end

# Scattering matrix Λ
function Λ(V, ϵ, z, ImpsT :: Array{GrapheneCoord, 1}) :: Array{ComplexF64, 2}
    nImps = length(ImpsT)
    ImpsT_Mat = repeat(ImpsT, 1, nImps)
    Imps_Mat = permutedims(ImpsT_Mat)
    Γ_Inv = (z - ϵ) .* Matrix{Int}(I, nImps, nImps)
    Prop = V^2 .* (map((x, y) -> Propagator(x, y, z), ImpsT_Mat, Imps_Mat))
    return (Γ_Inv .- Prop)
end

# Integrand used to calculate the local density
function Δρ_Integrand(Loc :: GrapheneCoord, V, ϵ, z, Imps :: Array{GrapheneCoord, 1})
    Λ_Inv = inv(Λ(V, ϵ, z, Imps))
    PropVector = reshape(map(x -> Propagator(Loc, x , z), Imps), (1, length(Imps)))
    return ((PropVector * Λ_Inv * permutedims(PropVector))[1])
end

# Local density function
function Δρ(Loc :: GrapheneCoord, V, ϵ, μ, Imps :: Array{GrapheneCoord, 1}) :: ComplexF64
    f_int(x) = Δρ_Integrand(Loc, V, ϵ, μ + 1im * x, Imps )
    res =  quadgk(f_int, -Inf, 0, Inf, maxevals  = NumEvals, rtol = ν)
    return (V^2 * res[1] / (2 * π))
end

# Impurity self-energy function
function Σ(z) :: ComplexF64
    return (z * Ω(z, 0, 0))
end

# Integrand used to compute the interaction energy
function F_I_Integrand(V, ϵ, z, Imps :: Array{GrapheneCoord, 1})
    Λ_ = Λ(V, ϵ, z, Imps)
    denom = z - ϵ - V^2 * Σ(z)
    return (-log(det(Λ_ ./ denom)) / (2 * π))
end

# Impurity interaction energy
function F_I(V, ϵ, μ, Imps :: Array{GrapheneCoord, 1})
    f_int(x) = F_I_Integrand(V, ϵ, μ + 1im * x, Imps)
    res = quadgk(f_int, -Inf, 0, Inf, maxevals = NumEvals, rtol = ν)
    return (res[1])
end

## Data Processing Functions

# In order to plot the calculated results as a lattice, one needs to create
# the correct coordinate arrays

function Data_Process(A_Lattice, B_Lattice)
    sz = size(A_Lattice)
    nPts = floor(Int,(sz[1] - 1)/2);
    # Lattice shift between the two sublattices
    latticeShift = -1/√(3) * 2.46;
    # Arrays
    d1s = -nPts : 1 : nPts; # d1 vectors
    d2s = -nPts : 1 : nPts; # d2 vectors

    D1S = repeat(d1s, 1, 2 * nPts + 1);
    D2S = repeat(d2s', 2 * nPts + 1, 1);

    # Coordinates of the carbon atoms
    XS = d1[1] .* D1S .+ d2[1] .* D2S;
    YS = d1[2] .* D1S .+ d2[2] .* D2S;

    # Flatten the coordinates and the data
    XS_A = reshape(XS, 1, sz[1]^2);
    YS_A = reshape(YS, 1, sz[1]^2);
    A_Lattice = reshape(A_Lattice, 1, sz[1]^2);

    XS_B = reshape(XS, 1, sz[1]^2);
    YS_B = reshape(YS, 1, sz[1]^2) .+ latticeShift;
    B_Lattice = reshape(B_Lattice, 1, sz[1]^2);

    # Combine all the coordinates and data for both sublattices
    XS = [XS_A XS_B];
    YS = [YS_A YS_B];
    dta = [A_Lattice B_Lattice];

    # Add values to make sure that the color scheme is centered around zero
    maxVal = maximum(dta);
    minVal = minimum(dta);
    maxAbs = max(abs.(maxVal), abs.(minVal))

    # Put the additional points away from the main data
    xMax = 1000;
    xMin = -1000;
    yMax = 1000;
    yMin = -1000;
    minV = - maxAbs;
    maxV = maxAbs;

    XS = [XS xMax xMin];
    YS = [YS yMax yMin];
    dta = [dta maxV minV];
    return([XS; YS; dta])
end
