using Profile
include("Graphene_Adatoms_Library.jl")

# Parameters
# ϵ = -1.0;      # On-site energy
ϵ = -7.2;      # On-site energy
μ = -0.2;      # Chemical potential
V = 3.0;       # Coupling strength
nPts = 25;     # Number of grid points away from zero

# Arrays
d1s = -nPts : 1 : nPts; # d1 vectors
d2s = -nPts : 1 : nPts; # d2 vectors

# Lattice vector matrices
D1S = repeat(d1s, 1, 2 * nPts + 1);
D2S = repeat(d2s', 2 * nPts + 1, 1);

# Impurity Arrangement
Loc_Origin = GrapheneCoord(0, 0, A)
Imp_Origin = Impurity(Loc_Origin, V, ϵ)

Loc_1 = GrapheneCoord(10, -10, A)
Loc_2 = GrapheneCoord(-10, 10, A)
Loc_3 = GrapheneCoord(4, 2, B)
Imp_1 = Impurity(Loc_1, V, ϵ)
Imp_2 = Impurity(Loc_2, V, ϵ)
Imp_3 = Impurity(Loc_3, V, ϵ)

Imps = [Imp_1, Imp_2]

# F_I(μ, Imps)
