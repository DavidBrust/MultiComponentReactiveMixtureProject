module EX3_CooledDuct

using Revise
using ColorSchemes
using GridVisualize
using CairoMakie
GridVisualize.default_plotter!(CairoMakie)



include("../notebooks/CooledDuct.jl")

nref = 3
grid = grid_2D(nref=nref)

solt, sys = setup_run_sim(grid, data);

sol_ss = run_ss(solt,sys);

function create_plot(data)
    (;iT) = data
	T = sol_ss[iT,:]

	levels = 41

	colormap = ColorSchemes.jet1[10:90]

	res = (650,250)
	vis = GridVisualizer(layout=(1,1), resolution=1 .*res, aspect=1.2, levels=0, colormap=colormap, colorbar=:horizontal, colorbarticks=6)
	
	scalarplot!(vis[1,1], grid, T, levels=0, linewidth=0, colorlevels=levels+2, title = "Temperature / K",
	xlabel = "Length / m", ylabel = "Height Y / m")
	
	sc = reveal(vis)
	fn = "img/out/2024-06-07/Cooled_Duct_nref3.png"
	GridVisualize.save(fn, sc)
end

end

