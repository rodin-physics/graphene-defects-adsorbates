using Distributed
using DelimitedFiles

@everywhere include("Graphene_Adatoms_Library.jl")
@everywhere include("Computation_Settings.jl")

# Calculation
@everywhere function f_A_Sublattice(x1, x2)
    println("x1 = $(x1), x2 = $(x2)")
    Loc = GrapheneCoord(x1, x2, A)
    res = F_I(V, ϵ, μ, [Imp_Origin, Loc])
    println(res)
    return res[1]
end

@everywhere function f_B_Sublattice(x1, x2)
    println("x1 = $(x1), x2 = $(x2)")
    Loc = GrapheneCoord(x1, x2, B)
    res = F_I(V, ϵ, μ, [Imp_Origin, Loc])
    println(res)
    return res[1]
end

resA = pmap(f_A_Sublattice, D1S, D2S)
resB = pmap(f_B_Sublattice, D1S, D2S)

writedlm("Data/Interaction/Interaction_A.txt", real(resA))
writedlm("Data/Interaction/Interaction_B.txt", real(resB))
