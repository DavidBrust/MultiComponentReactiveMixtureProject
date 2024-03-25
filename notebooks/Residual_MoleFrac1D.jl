### A Pluto.jl notebook ###
# v0.19.40

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
	using Revise
	using VoronoiFVM
	using ExtendableGrids, GridVisualize,ExtendableSparse
	using NLsolve, LinearSolve
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots, Printf
	using PlutoUI, Colors
	using Test
	using MultiComponentReactiveMixtureProject
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 2c4b16ab-8f60-467b-b608-2fea9fbc741c
md"""
# Introduction
Demo Notebook with $n=3$ ideal gas species in isothermal case with constant properties to assess the discretization approach used for the Multi-component transport model where the last species $X_3$ is inert.

The system consists of the three species $X_1, X_2, X_3$. The feed composition on molar basis is:
-  $X_1 = 0.8$
-  $X_2 = 0.2$
-  $X_3 = 0.0$

For the problem formulated here the molar fraction of the last species $X_3$ should vanish, since it is not participating in any reactions and it is not in the feed flow. In the model, the molar fraction of species $X_3$, $x_3$, is computed from solving $\Sigma x_i = 1$.

Solutions show small but non-vanishing values for $x_3$. A grid study is performed to see the influence of the grid spacing on the residuals of the molar fraction of species $X_3$.

Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=true))
"""

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
PlutoUI.TableOfContents(title="Demo")

# ╔═╡ 83fa22fa-451d-4c30-a4b7-834974245996
function grid1D(;nref=0)
	nx = (10.0)*2^(nref)
    h=1/nx
	X=(0:h:1)*ufac"cm"
	grid=simplexgrid(X)
	# catalyst region
	cellmask!(grid,[0.4]*ufac"cm",[0.6]*ufac"cm",2)	
	grid
end

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
gridplot(grid1D(nref=1), resolution=(600,200))

# ╔═╡ 832f3c15-b75a-4afe-8cc5-75ff3b4704d6
begin
	const Γ_left = 1
	const Γ_right = 2
end;

# ╔═╡ 94b0c012-b92c-41b7-9fa2-63e966dabd77
@doc MultiComponentReactiveMixtureProject.DMS_Info_isothermal()

# ╔═╡ 3440d4d8-3e03-4ff3-93f1-9afd7aaf9c41
md"""
## Chemical Reaction
The system consists of the three species $X_1, X_2, X_3$. The feed composition on molar basis is:
-  $X_1 = 0.8$
-  $X_2 = 0.2$
-  $X_3 = 0.0$

In the interior of the domain a reaction leading to increase in total species number (moles) occurs:

$X_2 \rightarrow 3 X_1$
"""

# ╔═╡ 3f56ffca-8ab8-4c32-b2b1-f2ec28ccf8b7
function calce(data,nref=0)
	
	grid = grid1D(nref=nref)
	inival,sys = PTR_init_system(1, grid, data)
	
	control = SolverControl(nothing, sys;)
		control.handle_exceptions=true
		control.Δu_opt=100

	times=[0,20]
	tsol=solve(sys;inival=inival,times,control)
	sol = tsol(tsol.t[end]);
	
	abse = maximum(abs.(sol[3,:]))
	nx = length(grid[Coordinates])
	
	sol,grid,sys,abse,nx
	
end

# ╔═╡ 184f70a9-f049-4017-ad28-027ae606d0ca
md"""
# Residual Molar Fraction
Investigate molar fraction of species 3 ($x_3$), which is the inert species and which should have a vanishing molar fraction, since it is not part of the feed.

Its value takes on a small but not vanishing value which probably results as a form of discretization error. It seems $x_3$ reaches its maximum at the outlet.

A grid refinement study is performed and the maximum value of $x_3$ is plotted over the number of grid points in the 1D grid.
"""

# ╔═╡ a5e9ffdb-083a-4192-bd8e-6f0ccb79598e
function make_data()
	MinKin = KinData{}(;
		ng = 3,
		gnames = [:A, :B, :C],
		Fluids = Vector{FluidProps}(undef,3),
		gn = Dict(1:3 .=> [:A, :B, :C]),
		nr = 1,
		rnames = [:R1],
		rn = Dict(1:1 .=> [:R1]),
		rni = Dict(value => key for (key, value) in Dict(1:1 .=> [:R1])),
		nuij = vcat(
			[3, -1, 0], #R1 : B -> 3A
		),
		ki_ref = Dict( [:R1] .=>  log.([5.0]) ),
		Ei = Dict( [:R1] .=> [0.0]*ufac"kJ/mol"),
		Ki_ref = Dict( [:R1] .=> [Inf]),
	
		Kj_ref = Dict( [:R1] .=> 0.0 ),
		TKj_ref = Dict( [:R1] .=>  0.0 ),
		ΔHj = Dict( [:R1] .=> 0.0 ),
	)
	
	data = ReactorData(;
		dim = 1,
		inlet_boundaries = [Γ_left],
		outlet_boundaries = [Γ_right],
		irradiated_boundaries = [],
		side_boundaries = [],
		solve_T_equation = false,
		constant_properties = true,
		is_reactive = true,
		kinpar = MinKin,
		lcat = 1.0,
		ng = 3,
		Tamb = 273.15,
		m = [2.0,6.0,21.0]*ufac"g/mol",
		X0 = [0.2, 0.8, 0.0],
		mfluxin = 0.01*ufac"kg/(m^2*s)"
	)
	return data
end

# ╔═╡ bb2ef3c0-96fc-4a90-a714-a0cfa08ac178
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	refmax = 6
	ae = []
	an = []
	asol = []
	agrid = []
	asys = []

	data_res = make_data()
	
	for nref=0:refmax
	 	sol,grid,sys,e,n = calce(data_res,nref)
	 	push!(ae, e)
	 	push!(an, n)
	 	push!(asol, sol)
	 	push!(agrid, grid)
	 	push!(asys, sys)
	end	 	
end
  ╠═╡ =#

# ╔═╡ abedcbf9-c99c-4969-b97a-3bc0295061bb
# ╠═╡ skip_as_script = true
#=╠═╡
let
	sol=asol[end]
	grid=agrid[end]
	(;ip,p) = data_res
	vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300))
	scalarplot!(vis, grid, sol[1,:], clear=false, label="x1")
	scalarplot!(vis, grid, sol[2,:], clear=false, color=:red, label="x2")
	scalarplot!(vis, grid, sol[3,:], clear=false, color=:blue, label="x3")

	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ 7dd363bf-d6e0-4dfb-b29f-85aa1fb62429
md"""
## Grid refinement
"""

# ╔═╡ c5190db5-ed3d-4084-ba16-496cc825fa9d
# ╠═╡ skip_as_script = true
#=╠═╡
let
	p=Plots.plot(xlabel="Number of gridpoints (1D)",ylabel="Abs. error")
	Plots.plot!(p,an,ae, xaxis=:log, yaxis=:log,m=:circle,label="max(abs(x$(data_res.ng)))",)
	Plots.plot!(p,an,1.0e-4 ./an, xaxis=:log, yaxis=:log, label="1/nx")
end
  ╠═╡ =#

# ╔═╡ 6f02c9cf-7670-4224-9239-b825d8331174
function runtests()
	sol,grid,sys,err,n = calce(make_data(),6)
    @test isapprox(err, 1.8146659931244213e-7)	
end;

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─2c4b16ab-8f60-467b-b608-2fea9fbc741c
# ╠═d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═83fa22fa-451d-4c30-a4b7-834974245996
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═832f3c15-b75a-4afe-8cc5-75ff3b4704d6
# ╟─94b0c012-b92c-41b7-9fa2-63e966dabd77
# ╟─3440d4d8-3e03-4ff3-93f1-9afd7aaf9c41
# ╠═3f56ffca-8ab8-4c32-b2b1-f2ec28ccf8b7
# ╟─184f70a9-f049-4017-ad28-027ae606d0ca
# ╠═bb2ef3c0-96fc-4a90-a714-a0cfa08ac178
# ╠═a5e9ffdb-083a-4192-bd8e-6f0ccb79598e
# ╠═abedcbf9-c99c-4969-b97a-3bc0295061bb
# ╟─7dd363bf-d6e0-4dfb-b29f-85aa1fb62429
# ╠═c5190db5-ed3d-4084-ba16-496cc825fa9d
# ╠═6f02c9cf-7670-4224-9239-b825d8331174
