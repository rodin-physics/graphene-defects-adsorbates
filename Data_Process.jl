using DelimitedFiles

include("Graphene_Adatoms_Library.jl")

# A_Lattice = readdlm("Data/Interaction/F_mu_04_A_Shallow.txt")
# B_Lattice = readdlm("Data/Interaction/F_mu_04_B_Shallow.txt")
#
# data = Data_Process(A_Lattice, B_Lattice)
#
# writedlm("Data/Interaction/Processed/F_mu_04_Shallow.dat", data')


A_Lattice = readdlm("Data/Density/rho_mu_02_Imp2_A_H.txt")
B_Lattice = readdlm("Data/Density/rho_mu_02_Imp2_B_H.txt")

data = Data_Process(A_Lattice, B_Lattice)

writedlm("Data/Density/Processed/rho_mu_02_Imp2_H.dat", data')
