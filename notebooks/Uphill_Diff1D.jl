### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c21e1942-628c-11ee-2434-fd4adbdd2b93
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise, Test
	using VoronoiFVM
	using GridVisualize, ExtendableGrids

	using LessUnitful	
	using PlutoUI, PlutoVista, Plots
	using MultiComponentReactiveMixtureProject
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 2c4b16ab-8f60-467b-b608-2fea9fbc741c
md"""
# Introduction
Demo Notebook to investigate the isothermal 1D multi-component transport problem with reacting gas mixture consisting of 5 reacting and 1 inert species with non-constant physical properties.

For high convective flowrates a "hump" in CO2 molar fraction occurs upstream of the reaction zone.
A grid study is performed to assess if that "hump" corresponds to "uphill diffusion" or is a numerical artifact.
"""

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
PlutoUI.TableOfContents(title="Uphill Diffusion 1D")

# ╔═╡ 83fa22fa-451d-4c30-a4b7-834974245996
function grid1D(nref=0)
	nx = (10.0)*2^(nref)
    h=1/nx
	X=(0:h:1)*ufac"cm"
	grid=simplexgrid(X)
	cellmask!(grid,[0.4]*ufac"cm",[0.6]*ufac"cm",2)	# catalyst region
	grid
end

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
gridplot(grid1D(), resolution=(600,200))

# ╔═╡ 832f3c15-b75a-4afe-8cc5-75ff3b4704d6
begin
	const Γ_left = 1
	const Γ_right = 2
end;

# ╔═╡ e17a383e-4e57-4c4d-9e08-0a24ae53af4f
@doc MultiComponentReactiveMixtureProject.DMS_Info_isothermal()

# ╔═╡ 9d0f4f0f-2a09-4f6c-a102-e45c8d39bd3e
md"""
# Uphill diffusion
For high convective flow rates, a "humb" of CO2 molar fraction upstream of the reaction zone occurs. A grid refinement study is performed to see if this is a numerical artifact.
"""

# ╔═╡ 04acab16-0848-40a7-8d70-9c2ae6c1ca9b
md"""
Calculate the height of the "hump" of CO2 molar fraction ($x_{\text{CO}_2}$) that is observable at high flow rates via
```math
	\Delta x_{\text{CO}_2} = \max (x_{\text{CO}_2}) - x_{\text{CO}_2, \text{in}}
```
For low flow rates, the phenomenon does not occur and $x_{\text{CO}_2}$ takes on its maximum at the inlet, from which follows that $\Delta x_{\text{CO}_2} =0$.
"""

# ╔═╡ 98cee333-4702-4974-bbd1-5ef8eed9a672
md"""
Select flowrate: $(@bind flowrate Select([:high, :low], default=:high))
"""

# ╔═╡ cc7e137e-8e7c-4d7e-ba70-517c3da0a42d
nflowin = flowrate == :high ? 0.005 : 0.0001

# ╔═╡ a56afd30-76ad-404a-9510-5f2d1b5dcf5d
data = ReactorData(
		inlet_boundaries=[Γ_left],
		outlet_boundaries=[Γ_right],
		irradiated_boundaries=[],
		side_boundaries=[],
		solve_T_equation=false,
		nflowin=nflowin,
		Treac = 273.15+650
)

# ╔═╡ 3f56ffca-8ab8-4c32-b2b1-f2ec28ccf8b7
function calc(data,nref=0)
	(;gni) = data
	
	grid = grid1D(nref)
	inival, sys = PTR_init_system(1, grid, data)
	
	control = SolverControl(nothing, sys;)
		control.handle_exceptions=true
		control.Δu_opt=100.0

	times = [0,100]
	tsol = solve(sys;inival=inival,times,control,verbose="a")
	sol = tsol(tsol.t[end])
	
	MaxCO2 = maximum(sol[gni[:CO2],:])
	DeltaCO2 = MaxCO2 - sol[gni[:CO2],1]
	nx = length(grid[Coordinates])
	
	return sol,grid,sys,DeltaCO2,nx
end

# ╔═╡ bb2ef3c0-96fc-4a90-a714-a0cfa08ac178
function PerformRefinement(data;refmax = 6)

	aDeltaCO2 = []
	an = []
	asol = []
	agrid = []

	
	for nref=0:refmax
		sol,grid,sys,DeltaCO2,n = calc(data,nref)
		push!(asol, sol)
		push!(agrid, grid)
		push!(aDeltaCO2, DeltaCO2)
		push!(an, n)			
	end

	return aDeltaCO2, an, asol, agrid
 end

# ╔═╡ 5da2971b-8d60-4ccf-b4e7-a3696b111e32
aDeltaCO2, an, asol, agrid = PerformRefinement(data)

# ╔═╡ abedcbf9-c99c-4969-b97a-3bc0295061bb
let	
	sol=asol[end]
	grid=agrid[end]
	(;ip,p,ng,gn) = data
	cols = distinguishable_colors(ng, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
	pcols = map(col -> (red(col), green(col), blue(col)), cols)
	vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300))
	for i=1:ng
		scalarplot!(vis, grid, sol[i,:], clear=false, color=pcols[i],label=gn[i])
	end
	reveal(vis)
end

# ╔═╡ 05e5e167-7c1c-43f4-8a8d-e173f2fb7029
md"""
## Grid Refinement
"""

# ╔═╡ c5190db5-ed3d-4084-ba16-496cc825fa9d
let
	p=Plots.plot(xlabel="Number of gridpoints (1D)",ylabel="Delta xCO2")
	Plots.plot!(p,an,aDeltaCO2, m=:circle,label="$(flowrate) flow rate",ylim=(0,1.25*maximum(aDeltaCO2)))
end

# ╔═╡ 49033775-2358-4463-a5dd-cb5629f35142
@test isapprox(aDeltaCO2[end], 0.01702152459555617)

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─2c4b16ab-8f60-467b-b608-2fea9fbc741c
# ╟─d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═83fa22fa-451d-4c30-a4b7-834974245996
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═832f3c15-b75a-4afe-8cc5-75ff3b4704d6
# ╠═e17a383e-4e57-4c4d-9e08-0a24ae53af4f
# ╟─9d0f4f0f-2a09-4f6c-a102-e45c8d39bd3e
# ╟─04acab16-0848-40a7-8d70-9c2ae6c1ca9b
# ╟─98cee333-4702-4974-bbd1-5ef8eed9a672
# ╠═cc7e137e-8e7c-4d7e-ba70-517c3da0a42d
# ╠═a56afd30-76ad-404a-9510-5f2d1b5dcf5d
# ╠═3f56ffca-8ab8-4c32-b2b1-f2ec28ccf8b7
# ╠═bb2ef3c0-96fc-4a90-a714-a0cfa08ac178
# ╠═5da2971b-8d60-4ccf-b4e7-a3696b111e32
# ╠═abedcbf9-c99c-4969-b97a-3bc0295061bb
# ╟─05e5e167-7c1c-43f4-8a8d-e173f2fb7029
# ╠═c5190db5-ed3d-4084-ba16-496cc825fa9d
# ╠═49033775-2358-4463-a5dd-cb5629f35142
