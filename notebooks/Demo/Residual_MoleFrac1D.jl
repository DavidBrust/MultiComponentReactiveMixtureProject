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

	using LessUnitful
	
	using PlutoVista, Plots, Printf
	using PlutoUI, Colors

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 2c4b16ab-8f60-467b-b608-2fea9fbc741c
md"""
# Introduction
Demo Notebook with $\nu=3$ ideal gas species and constant properties to assess the discretization approach used for the Multi-component transport model where the $\nu$ (last) species is inert.

In the case setup in this notebook, the molar fraction of the last species $x_\nu$ should vanish, since it is not participating in any reactions and it is not in the feed flow. In the model, gas species $\nu$ is computed from $\Sigma x_i = 1$.

Solutions show small but non-vanishing values for $x_\nu$. A grid study is performed to see the influence of the grid spacing on the residuals of the molar fraction of species $\nu$.

Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=false))
"""

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
PlutoUI.TableOfContents(title="Demo")

# ╔═╡ 83fa22fa-451d-4c30-a4b7-834974245996
function grid1D(nref=0)
	nx = (10.0)*2^(nref)
    h=1/nx
	X=(0:h:1)*ufac"cm"
	grid=simplexgrid(X)
	# catalyst region
	cellmask!(grid,[0.4]*ufac"cm",[0.6]*ufac"cm",2)	
	grid
end

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
gridplot(grid1D(1), resolution=(600,200))

# ╔═╡ 832f3c15-b75a-4afe-8cc5-75ff3b4704d6
begin
	const Γ_left = 1
	const Γ_right = 2
end;

# ╔═╡ 0fadb9d2-1ccf-4d44-b748-b76d911784ca
md"""
## Mass Continuity
```math
\begin{align}
	\frac{\partial \rho}{\partial t} + \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\end{align}
```
"""

# ╔═╡ b94513c2-c94e-4bcb-9342-47ea48fbfd14
md"""
## Species Mass Transport
```math
\begin{align}
	\frac{\partial \rho_i}{\partial t} + \nabla \cdot \left( \vec \Phi_i + \rho_i \vec v \right ) - R_i &= 0 ~, \qquad i = 1 ... \nu \\
		\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= -\sum_{j=1 \atop j \neq i}^{\nu} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
		\sum_{i=1}^\nu x_i &= 1
\end{align}
```
"""

# ╔═╡ c886dd12-a90c-40ab-b9d0-32934c17baee
md"""
where $\rho$ is the (total) mixture density, $\vec v$ is the mass-averaged (barycentric)  mixture velocity calculated with the Darcy equation, $x_i$, $w_i$ and $M_i$ are the molar fraction, mass fraction and molar mass of species $i$ respectively, $\vec \Phi_i$ is the mass flux of species $i$ ($\frac{\text{kg}}{\text{m}^2 \text{s}}$) and $R_i$ is the species mass volumetric source/sink ($\frac{\text{mol}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
"""

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

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
function runSimRef(data,nref=0)

	mygrid=grid1D(nref)
	strategy = nothing
	times=[0,20]


	
	(;p,ip,X0)=data
	ng=ngas(data)
	
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

	control = SolverControl(strategy, sys;)
		control.handle_exceptions=true
		control.maxiters=250
		control.Δu_opt=1.0e5

	tsol=solve(sys;inival=inival,times,control)
	sol = tsol(tsol.t[end]);

	sol,mygrid,sys
end

# ╔═╡ 184f70a9-f049-4017-ad28-027ae606d0ca
md"""
# Residual Molar Fraction
Investigate molar fraction of species 3 ($x_3$), which is the inert species and which should have a vanishing molar fraction, since it is not part of the feed.

Its value takes on a small but not vanishing value which probably results as a form of discretization error. It seems $x_3$ reaches its maximum at the outlet.

A grid refinement study is performed and the maximum value of $x_3$ is plotted over the number of grid points in the 1D grid.
"""

# ╔═╡ 7dd363bf-d6e0-4dfb-b29f-85aa1fb62429
md"""
## Grid Refinement
"""

# ╔═╡ 3f56ffca-8ab8-4c32-b2b1-f2ec28ccf8b7
function calce(data,nref=0)	
	sol,grid,sys = runSimRef(data,nref)
	abse = maximum(abs.(sol[3,:]))
	nx = length(grid[Coordinates])
	sol,grid,sys,abse,nx
end

# ╔═╡ bb2ef3c0-96fc-4a90-a714-a0cfa08ac178
begin
	if RunSim
	 	refmax = 6
	 	ae = []
	 	an = []
		asol = []
		agrid = []
		asys = []

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
		
		data_res = ReactorData(;
			inlet_boundaries=[Γ_left],
			outlet_boundaries=[Γ_right],
			irradiated_boundaries=[],
			side_boundaries=[],
			solve_T_equation=false,
			constant_properties = true,
			is_reactive = true,
			kinpar=MinKin,
			lcat = 1.0,
			ng = 3,
			Tamb = 273.15,
			m = [2.0,6.0,21.0]*ufac"g/mol",
			X0 = [0.2, 0.8, 0.0],
			mfluxin = 0.01*ufac"kg/(m^2*s)"
		)
			
	 	for nref=0:refmax
	 		sol,grid,sys,e,n = calce(data_res,nref)
	 		push!(ae, e)
	 		push!(an, n)
			push!(asol, sol)
			push!(agrid, grid)
			push!(asys, sys)
	 	end	 	
	end
 end

# ╔═╡ abedcbf9-c99c-4969-b97a-3bc0295061bb
let
	if RunSim
		sol=asol[end]
		grid=agrid[end]
		(;ip,p) = data_res
		vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300))
		scalarplot!(vis, grid, sol[1,:], clear=false, label="x1")
		scalarplot!(vis, grid, sol[2,:], clear=false, color=:red, label="x2")
		scalarplot!(vis, grid, sol[3,:], clear=false, color=:blue, label="x3")
	
		reveal(vis)
	end
end

# ╔═╡ c5190db5-ed3d-4084-ba16-496cc825fa9d
if RunSim
	p=Plots.plot(xlabel="Number of gridpoints (1D)",ylabel="Abs. error")
	Plots.plot!(p,an,ae, xaxis=:log, yaxis=:log,m=:circle,label="max(abs(x$(data_res.ng))",)
	Plots.plot!(p,an,1.0e-4 ./an, xaxis=:log, yaxis=:log, label="1/nx")
end

# ╔═╡ e4fb2336-781d-4c9c-b038-d9325ced57b6
function checkinout(sys,sol)	
	tfact=TestFunctionFactory(sys)
	tf_in=testfunction(tfact,[Γ_right],[Γ_left])
	tf_out=testfunction(tfact,[Γ_left],[Γ_right])	
	(;in=integrate(sys,tf_in,sol),out=integrate(sys,tf_out,sol) )
end

# ╔═╡ 3ffbe1ba-d5aa-4bf0-a5c4-92cec7e0c199
let
	if RunSim
		(;ip,gn,gni,m,mfluxin,mmix0,X0,ng) = data_res
		sol=asol[end]
		sys=asys[end]
		grid=agrid[end]
		
		in_,out_=checkinout(sys,sol)
	
		nout(i) = out_[i]/m[i]
		nin(i) = mfluxin/mmix0 *bareas(Γ_left,sys,grid)*X0[i]
		RI=sum(integrate(sys,FixedBed.DMS_reaction,sol),dims=2) # reaction integral

		println("Total mass inflows and outflows:")
		@printf "IN: %2.6e \t OUT: %2.6e \t REACT: %2.6e kg/hr \nSUM: %2.6e kg/hr\n\n" in_[ip]/ufac"kg/hr" out_[ip]/ufac"kg/hr" RI[ip]/ufac"kg/hr" (in_[ip]+out_[ip]+RI[ip])/ufac"kg/hr"
		
		println("Molar species inflows, outflows and reaction integrals:")
		for i = 1:ng
			@printf "%i\tIN: %2.6e \t OUT: %2.6e \t REACT: %2.6e mol/hr \n\tSUM: %2.6e mol/hr\n" i nin(i)/ufac"mol/hr" nout(i)/ufac"mol/hr" -RI[i]/m[i]/ufac"mol/hr" (nin(i)+nout(i)-RI[i]/m[i])/ufac"mol/hr"
		end
	end
end

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─2c4b16ab-8f60-467b-b608-2fea9fbc741c
# ╠═d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═83fa22fa-451d-4c30-a4b7-834974245996
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═832f3c15-b75a-4afe-8cc5-75ff3b4704d6
# ╟─0fadb9d2-1ccf-4d44-b748-b76d911784ca
# ╟─b94513c2-c94e-4bcb-9342-47ea48fbfd14
# ╟─c886dd12-a90c-40ab-b9d0-32934c17baee
# ╟─3440d4d8-3e03-4ff3-93f1-9afd7aaf9c41
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╟─184f70a9-f049-4017-ad28-027ae606d0ca
# ╟─abedcbf9-c99c-4969-b97a-3bc0295061bb
# ╠═3ffbe1ba-d5aa-4bf0-a5c4-92cec7e0c199
# ╟─7dd363bf-d6e0-4dfb-b29f-85aa1fb62429
# ╟─c5190db5-ed3d-4084-ba16-496cc825fa9d
# ╠═bb2ef3c0-96fc-4a90-a714-a0cfa08ac178
# ╠═3f56ffca-8ab8-4c32-b2b1-f2ec28ccf8b7
# ╠═e4fb2336-781d-4c9c-b038-d9325ced57b6
