using Plots
using LaTeXStrings
using DelimitedFiles

include("Graphene_Adatoms_Library.jl")

data = readdlm("Data/Interaction/Processed/F_mu_04_H.dat")

XS = data[:, 1]
YS = data[:, 2]
FI = data[:, 3]

bound = 1e-5;
bound = bound * 1e3;

FI = FI * 1e3

pyplot();
scatter(XS,YS,
        marker_z = FI,
        markerstrokecolor = :white,
        markerstrokewidth = 0.001,
        markersize = 3,
        leg = false,
        aspect_ratio=1,
        # xaxis = (L"\AA", font(20, "Serif")),
        # yaxis = (L"\AA", font(20, "Serif")),
        xtickfont = font(12, "Serif"),
        ytickfont = font(12, "Serif"),
        ylims = (-60,60),
        xlims = (-60,60),
        size = (500,400),
        color = :coolwarm,
        clim = (-bound, bound),
        colorbar = true
        )
println("Plot Done")
savefig("F_mu_04_H.pdf")
println("Plot Saved")
