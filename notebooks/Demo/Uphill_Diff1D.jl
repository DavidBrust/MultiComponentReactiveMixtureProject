### A Pluto.jl notebook ###
# v0.19.36

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
	Pkg.activate(joinpath(@__DIR__,"../.."))
	using Revise
	using VoronoiFVM
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve, LinearSolve
	using StaticArrays
	using Test
	using LessUnitful
	
	using PlutoVista, Plots
	using PlutoUI, Colors, Printf

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 2c4b16ab-8f60-467b-b608-2fea9fbc741c
md"""
# Introduction
Demo Notebook to investigate the isothermal 1D multi-component transport problem with reacting gas mixture consisting of 5 reacting and 1 inert species with non-constant physical properties.

For high convective flowrates a "hump" in CO2 molar fraction occurs upstream of the reaction zone.
A grid study is performed to assess if that "hump" corresponds to "uphill diffusion" or is a numerical artifact.

Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=true))
"""

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
PlutoUI.TableOfContents(title="Demo")

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
@doc FixedBed.DMS_Info_isothermal()

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
Select flowrate: $(@bind flowrate Select([:high, :low]))
"""

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
function runSimRef(data,nref=0)
	(;p,ip,X0,ng,nflowin)=data
	
	mygrid=grid1D(nref)

	
	times=flowrate == :high ? [0,20.0] : [0,100.0]

	sys=VoronoiFVM.System( 	mygrid;
							data=data,
							flux=FixedBed.DMS_flux,
							reaction=FixedBed.DMS_reaction,
							storage=FixedBed.DMS_storage,
							bcondition=FixedBed.PCR_bcond,
							bflux=FixedBed.PCR_bflux,
							bstorage=FixedBed.PCR_bstorage,
							boutflow=FixedBed.DMS_boutflow,
							outflowboundaries=data.outlet_boundaries,
							assembly=:edgewise
							)
	
	enable_species!(sys; species=collect(1:(ng+1))) # gas phase species pi & ptotal	
	
	inival=unknowns(sys)

	inival[ip,:].=p
	for i=1:ng
		#inival[i,:] .= 1.0/ng
		inival[i,:] .= X0[i]
	end
	
	control = SolverControl(nothing, sys;)
		control.handle_exceptions=true
		control.Δu_opt=1_000.0
		control.maxiters=250

	tsol=solve(sys;inival=inival,times,control,verbose="nae")
	sol = tsol(tsol.t[end]);

	sol,mygrid,sys
end

# ╔═╡ 05e5e167-7c1c-43f4-8a8d-e173f2fb7029
md"""
## Grid Refinement
"""

# ╔═╡ 3f56ffca-8ab8-4c32-b2b1-f2ec28ccf8b7
function calc(data,nref=0)
	(;gni) = data
	sol,grid,sys = runSimRef(data,nref)
	
	MaxCO2 = maximum(sol[gni[:CO2],:])
	DeltaCO2 = MaxCO2 - sol[gni[:CO2],1]
	nx = length(grid[Coordinates])
	sol,grid,sys,DeltaCO2,nx
end

# ╔═╡ bb2ef3c0-96fc-4a90-a714-a0cfa08ac178
begin
	#if RunSim
	 	refmax = 6
	 	aDeltaCO2 = []
	 	an = []
		asol = []
		agrid = []
		asys = []
		
		nflowin = flowrate == :high ? 0.005 : 0.0001
		data_uh = ReactorData(
			inlet_boundaries=[Γ_left],
			outlet_boundaries=[Γ_right],
			irradiated_boundaries=[],
			side_boundaries=[],
			solve_T_equation=false,
			nflowin=nflowin,
			Treac = 273.15+650
		)
		
	 	for nref=0:refmax
	 		sol,grid,sys,DeltaCO2,n = calc(data_uh,nref)
	 		push!(asol, sol)
			push!(agrid, grid)
			push!(asys, sys)
			push!(aDeltaCO2, DeltaCO2)
	 		push!(an, n)			
	 	end	 	
	#end
 end

# ╔═╡ abedcbf9-c99c-4969-b97a-3bc0295061bb
let
	if RunSim
		sol=asol[end]
		grid=agrid[end]
		(;ip,p,ng,gn) = data_uh
		cols = distinguishable_colors(ng, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
		pcols = map(col -> (red(col), green(col), blue(col)), cols)
		vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300))
		for i=1:ng
			scalarplot!(vis, grid, sol[i,:], clear=false, color=pcols[i],label=gn[i])
		end
		reveal(vis)
	end
end

# ╔═╡ aee47c66-97ba-43fa-a8ec-dad7c0d4f20c
let
	if RunSim
		(;gn,gni,m,mfluxin,mmix0,X0,ng,ip) = data_uh
		sol=asol[end]
		sys=asys[end]
		in_,out_=FixedBed.DMS_checkinout(sol,sys,data_uh)
	
		nout(i) = out_[i]/m[i]
		nin(i) = mfluxin/mmix0 *1.0*ufac"m^2"*X0[i]
		RI=sum(integrate(sys,FixedBed.DMS_reaction,sol),dims=2) # reaction integral

		println("Total mass inflows and outflows:")
		@printf "IN: %2.6e \t OUT: %2.6e \t REACT: %2.6e kg/hr \nSUM: %2.6e kg/hr\n\n" in_[ip]/ufac"kg/hr" out_[ip]/ufac"kg/hr" RI[ip]/ufac"kg/hr" (in_[ip]+out_[ip]+RI[ip])/ufac"kg/hr"
		
		println("Molar species inflows, outflows and reaction integrals:")
		for i = 1:ng
			@printf "%s\tIN: %2.6e \t OUT: %2.6e \t REACT: %2.6e mol/hr \n\tSUM: %2.6e mol/hr\n" gn[i] nin(i)/ufac"mol/hr" nout(i)/ufac"mol/hr" -RI[i]/m[i]/ufac"mol/hr" (nin(i)+nout(i)-RI[i]/m[i])/ufac"mol/hr"
		end
	end
end

# ╔═╡ c5190db5-ed3d-4084-ba16-496cc825fa9d
if RunSim
	p=Plots.plot(xlabel="Number of gridpoints (1D)",ylabel="Delta xCO2")
	Plots.plot!(p,an,aDeltaCO2, m=:circle,label="$(flowrate) flow rate",ylim=(0,1.25*maximum(aDeltaCO2)))
end

# ╔═╡ 7a798fbc-472a-44fc-a058-f1db25b09459
@test aDeltaCO2[end] ≈ 0.01702152459546069

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─2c4b16ab-8f60-467b-b608-2fea9fbc741c
# ╠═d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═83fa22fa-451d-4c30-a4b7-834974245996
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═832f3c15-b75a-4afe-8cc5-75ff3b4704d6
# ╠═e17a383e-4e57-4c4d-9e08-0a24ae53af4f
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╟─9d0f4f0f-2a09-4f6c-a102-e45c8d39bd3e
# ╟─04acab16-0848-40a7-8d70-9c2ae6c1ca9b
# ╟─98cee333-4702-4974-bbd1-5ef8eed9a672
# ╟─abedcbf9-c99c-4969-b97a-3bc0295061bb
# ╟─aee47c66-97ba-43fa-a8ec-dad7c0d4f20c
# ╟─05e5e167-7c1c-43f4-8a8d-e173f2fb7029
# ╠═c5190db5-ed3d-4084-ba16-496cc825fa9d
# ╠═7a798fbc-472a-44fc-a058-f1db25b09459
# ╠═bb2ef3c0-96fc-4a90-a714-a0cfa08ac178
# ╠═3f56ffca-8ab8-4c32-b2b1-f2ec28ccf8b7
