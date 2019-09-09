using QuadGK
using LinearAlgebra

## Parameters
const ν = 1e-2;       # Small number for relative tolerance
const η = 1e-8;       # Small number for absolute tolerance
const NumEvals = 1e6; # Maximum number of evaluations in integrals

# Graphene hopping integral and lattice vectors in Angstroms
const t = 2.8;
const d1 = 2.46 .* [1 √(3)] ./ 2;
const d2 = 2.46 .* [-1 √(3)] ./ 2;

# Sublattice spinors
const A = [1, 0]
const B = [0, 1]

# Graphene Coordinate type: R = u * d1 + v * d2
struct GrapheneCoord
    u :: Int
    v :: Int
    sublattice :: Array{Int64, 1}
end

# Localized Impurity type
struct Impurity
    pos :: GrapheneCoord
    V :: Float64    # Coupling energy
    ϵ :: Float64    # On-site energy
end

## Helper functions
# Auxiliary W function
function W(z :: ComplexF64, x :: Float64)
    return (((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x))
end

# Auxiliary Y function used in the calculation of Ω
function Y(n :: Int, z, x)
    W_ = W(z, x)
    return (W_ - √(W_ - 1) * √(W_ + 1))^n / (√(W_ - 1) * √(W_ + 1))
end

# When computing Ω, occasionally the integrand becomes small enough to give NaN
# This helper functions is used to catch these instances
function Ω_Integrand(z, u, v, x :: Float64)
    res = exp(1.0im * (u - v) * x) / cos(x) * Y(abs.(u + v), z, x)
    if isnan(res) == true
        return 0.0 + 0.0im
    else
        return res
    end
end

## Main Functions
function Ω(z, u, v)
    f_int(x) = Ω_Integrand(z, u, v, x)
    res = quadgk(f_int, 0, 2 * π, maxevals = NumEvals)
    return (res[1] / (8.0 * π * t^2) )
end

# The propagator function picks out the correct element of the Ξ matrix based
# on the sublattices of the graphene coordinates
function Propagator(Imp_l :: GrapheneCoord, Imp_m :: GrapheneCoord, z)
    u = Imp_l.u - Imp_m.u
    v = Imp_l.v - Imp_m.v
    if Imp_l.sublattice == Imp_m.sublattice
        return (z * Ω(z, u, v))
    elseif ([Imp_l.sublattice, Imp_m.sublattice] == [A, B])
        return (- t * (Ω(z, u, v) + Ω(z, u + 1, v) + Ω(z, u, v + 1)))
    elseif ([Imp_l.sublattice, Imp_m.sublattice] == [B, A])
        return (- t * (Ω(z, u, v) + Ω(z, u - 1, v) + Ω(z, u, v - 1)))
    else
        error("Illegal sublattice parameter")
    end
end

# Inverse full impurity Green's function
function Λ_Inv(z, Imps :: Array{Impurity, 1})

    nImps = length(Imps)
    ϵ = map(x -> x.ϵ, Imps)         # Array of impurity on-site energies
    V = map(x -> x.V, Imps)         # Array of impurity coupling energies
    loc = map(x -> x.pos, Imps)     # Array of impurity positions

    ϵ_Mat = Diagonal(ϵ)             #  Diagonal matrix of on-site energies
    V_Mat = Diagonal(V)             # Diagonal matrix of coupling energies


    # Impurity Green's function for isolated impurities (as a diagonal matrix)
    Γ_Inv = z .* Matrix{Int}(I, nImps, nImps) .- ϵ_Mat

    locT_Mat = repeat(loc, 1, nImps)    # Impurity position matrix
    loc_Mat = permutedims(locT_Mat)     # Transpose of impurity position matrix


    Prop = V_Mat * (map((x, y) -> Propagator(x, y, z), locT_Mat, loc_Mat)) * V_Mat
    return (Γ_Inv .- Prop)
end

# Integrand used to calculate the local density
function Δρ_Integrand(Loc, z, Imps)
    Λ = inv(Λ_Inv(z, Imps))
    imp_loc = map(x -> x.pos, Imps)

    V = map(x -> x.V, Imps)
    V_Mat = Diagonal(V)

    PropVectorL = reshape(map(x -> Propagator(Loc, x , z), imp_loc), (1, length(Imps)))
    PropVectorR = reshape(map(x -> Propagator(x, Loc , z), imp_loc), (length(Imps), 1))

    return ((PropVectorL * V_Mat * Λ * V_Mat * PropVectorR)[1])
end

# Local density function
function Δρ(Loc, μ, Imps)
    f_int(x) = Δρ_Integrand(Loc, μ + 1im * x, Imps )
    res = quadgk(f_int, -Inf, 0, Inf, maxevals  = NumEvals, rtol = ν, atol = η)
    return (res[1] / (2 * π))
end

# Impurity self-energy function
function Σ(z)
    return (z * Ω(z, 0, 0))
end

# Integrand used to compute the interaction energy. We catch NaN and return 0
function F_I_Integrand(z, Imps)
    nImps = length(Imps)
    Λ_Inv_ = Λ_Inv(z, Imps)

    ϵ = map(x -> x.ϵ, Imps)
    V = map(x -> x.V, Imps)
    loc = map(x -> x.pos, Imps)

    Σ_ = Σ(z)

    ϵ_Mat = Diagonal(ϵ)
    V_Mat = Diagonal(V)

    Λ0_Inv = z .* Matrix{Int}(I, nImps, nImps) .- ϵ_Mat .- V_Mat .* Σ_ .* V_Mat

    res = -log(det(Λ_Inv_ * inv(Λ0_Inv))) / (2 * π)

    if isnan(res) == true
        return 0.0 + 0.0im
    else
        return res
    end
end

# Impurity interaction energy
function F_I(μ, Imps)
    f_int(x) = F_I_Integrand(μ + 1im * x, Imps)
    res = quadgk(f_int, -Inf, 0, Inf, maxevals = NumEvals, rtol = ν, atol = η)
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

    return([XS; YS; dta])
end
