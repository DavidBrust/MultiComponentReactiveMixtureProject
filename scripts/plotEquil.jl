# script to plot STC efficiency from exported data (.scv)
using CSV, DataFrames, Plots, Colors

df = CSV.read("./data/RWGS_Methanation_equil/RWGS_Methanation_equil.csv", DataFrame)

p = Plots.plot(
     xguide="Temperature / °C",
     yguide="Mole fraction / -",
     title="Equilibrium Composition",
     legend=:outertopright,
     size=(450,300),
     )

cols = distinguishable_colors(5, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
cols[1]

plot!(p, df.T,df.MoleFrac_CO_out, c=cols[1], label="CO", lw=2,ls=:auto)
plot!(p, df.T,df.MoleFrac_H2_out, c=cols[2], label="H2", lw=2,ls=:auto)
plot!(p, df.T,df.MoleFrac_CH4_out, c=cols[3], label="CH4", lw=2,ls=:auto)
plot!(p, df.T,df.MoleFrac_H2O_out, c=cols[4], label="H2O", lw=2,ls=:auto)
plot!(p, df.T,df.MoleFrac_CO2_out, c=cols[5], label="CO2", lw=2,ls=:auto)

lens!([550, 650], [0.1, 0.2], inset = (1, bbox(0.52, 0.05, 0.25, 0.25)))

savefig(p,"./img/out/RWGS_equil.svg")