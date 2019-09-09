using Plots
using LaTeXStrings
using DelimitedFiles

include("Graphene_Adatoms_Library.jl")

data = readdlm("Data/Density/Processed/rho_mu_04_Imp1_H.dat")

XS = data[:, 1]
YS = data[:, 2]
ρ = data[:, 3]

ρ = ρ * 1e3;

bound = 1e-3;
bound = bound * 1e3;

pyplot();
scatter(XS,YS,
        marker_z = ρ,
        markerstrokecolor = :white,
        markerstrokewidth = 0.001,
        markersize = 9,
        leg = false,
        aspect_ratio = 1,
        # xaxis = (L"\AA", font(20, "Serif")),
        # yaxis = (L"\AA", font(20, "Serif")),
        xtickfont = font(12, "Serif"),
        ytickfont = font(12, "Serif"),
        ylims = (-20,20),
        xlims = (-20,20),
        size = (500,400),
        color = :coolwarm,
        clim = (-bound, bound),
        colorbar = true,
        colorbar_title = L"\Delta\rho\times 10^3"
        )
println("Plot Done")
savefig("rho_mu_04_Imp1_H_Zoom.pdf")
println("Plot Saved")

# 2 * π / (2*0.4/(t * 2.46e-8) * 1e-7)
