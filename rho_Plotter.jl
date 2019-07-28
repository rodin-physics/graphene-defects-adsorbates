using Plots
using LaTeXStrings
using DelimitedFiles

include("Graphene_Adatoms_Library.jl")

A_Lattice = readdlm("Data/Density/rho_A.txt")
B_Lattice = readdlm("Data/Density/rho_B.txt")

data = Data_Process(A_Lattice, B_Lattice)

XS = data[1,:]
YS = data[2,:]
ρ = data[3,:]

pyplot();
plot(   leg = false,
        aspect_ratio=1,
        # markeralpha = 0.25,
        # xaxis = (L"\frac{U}{m-\epsilon}", font(20)),
        # yaxis = (L"\frac{\mu}{m}", font(20)),
        xtickfont = font(12),
        ytickfont = font(12),
        # yticks = 0.25:0.25:1,
        ylims = (-60,60),
        xlims = (-60,60),
        size = (400,400)
        )

scatter!(XS',YS', zcolor = ρ',
        markerstrokecolor = :white,
        markersize = 3,
        color = :RdBu)

savefig("Test.pdf")
