# script to plot STC efficiency from exported data (.scv)
using CSV, DataFrames, Plots

df = CSV.read("./data/out/STC_bak.csv", DataFrame)

p = Plots.plot(
    xguide="Feed Flow / sccm",
    yguide="Solar Concentration / kW m-2",
    title="Solar-To-Chemical Efficiency",
    size=(400,266)
    )

heatmap!(p,unique(df.Qflow)./ufac"cm^3/minute",unique(df.C),reshape(unique(df.STCs), 12,:),right_margin=5Plots.mm)
#savefig(p,"./img/out/STC_bak.svg")