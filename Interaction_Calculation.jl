using Distributed
using DelimitedFiles

@everywhere include("Computation_Settings.jl")

# Calculation
@everywhere function f_A_Sublattice(x1, x2)
    println("x1 = $(x1), x2 = $(x2)")
    Loc = GrapheneCoord(x1, x2, A)
    Imp = Impurity(Loc, V, ϵ)
    res = F_I(μ, [Imp_Origin, Imp])
    println(res)
    return res[1]
end

@everywhere function f_B_Sublattice(x1, x2)
    println("x1 = $(x1), x2 = $(x2)")
    Loc = GrapheneCoord(x1, x2, B)
    Imp = Impurity(Loc, V, ϵ)
    res = F_I(μ, [Imp_Origin, Imp])
    println(res)
    return res[1]
end

resA = pmap(f_A_Sublattice, D1S, D2S)
resB = pmap(f_B_Sublattice, D1S, D2S)

writedlm("Data/Interaction/F_mu_04_A_Shallow.txt", real(resA))
writedlm("Data/Interaction/F_mu_04_B_Shallow.txt", real(resB))
