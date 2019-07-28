include("Graphene_Adatoms_Library.jl")

# Parameters
ϵ = -2.4;      # On-site energy
μ = -0.0;      # Chemical potential
V = 4.0;       # Coupling strength
nPts = 3;      # Number of grid points away from zero

ω = -0.4;      # Energy parameter

# Arrays
d1s = -nPts : 1 : nPts; # d1 vectors
d2s = -nPts : 1 : nPts; # d2 vectors

D1S = repeat(d1s, 1, 2 * nPts + 1);
D2S = repeat(d2s', 2 * nPts + 1, 1);

# Impurity Arrangement
Imp_Origin = GrapheneCoord(0, 0, A)
Imp1 = GrapheneCoord(0, 0, A)
Imp2 = GrapheneCoord(3, -3, B)
Imp3 = GrapheneCoord(3, 3, A)

Imps = [Imp3, Imp2]
