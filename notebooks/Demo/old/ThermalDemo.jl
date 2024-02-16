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
	using Revise, Test
	using VoronoiFVM
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve, LinearSolve,Pardiso, ILUZero
	using StaticArrays

	using LessUnitful
	using DataFrames
	using PlutoVista, Plots
	using PlutoUI, Colors
	using MultiComponentReactiveMixtureProject
	using Printf
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═╡ skip_as_script = true
#=╠═╡
PlutoUI.TableOfContents(title="Photo Thermal (PT) Reactor")
  ╠═╡ =#

# ╔═╡ b2791860-ae0f-412d-9082-bb2e27f990bc
md"""
# Introduction
Demonstration notebook for the photo thermal catalytic reactor (PTR) model. Solve energy equation alongside multicomponent species transport. Include reactive gas mixture (CO2,H2,CO,CH4,H2O,N2) with variable physical properties and a Ni based catalyst described with kinetics from published literature.

Select problem dimension: $(@bind dim Select([2, 3], default=2))

Select grid refinement level: $(@bind nref Select([0,1,2,3], default=0))

Check the box to __start the simulation__: $(@bind RunSim PlutoUI.CheckBox(default=true))
"""

# ╔═╡ 4e05ab31-7729-4a4b-9c14-145118477715
# ╠═╡ skip_as_script = true
#=╠═╡
if dim == 3
	@bind xcut Slider(linspace(0,16,17)*ufac"cm",show_value=true,default=8*ufac"cm")
end
  ╠═╡ =#

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	grid, inb,irrb,outb,sb,catr =  grid_boundaries_regions(dim, nref=nref)
	if dim == 2
		gridplot(grid, resolution=(660,300), aspect=4.0, zoom=2.8)
	else
		gridplot(grid,  xplane=xcut, resolution=(660,460), zoom=1.8, )
	end
end
  ╠═╡ =#

# ╔═╡ a078e1e1-c9cd-4d34-86d9-df4a052b6b96
md"""
# Introduction
Reactor Simulation of PC reactor for cylindrical symmetrical geometry (2D).

This notebook covers overall mass transport through porous material based on Darcy equation, multicomponent species transport based on Maxwell-Stefan equations for diffusion and superimposed by convective (bulk) transport with the velocity field obtained by Darcy equation.

Also the thermal energy equation is solved, taking into account convective-diffusive heat fluxes and irradiation boundary conditions.
"""

# ╔═╡ a1ea393e-f123-4ad0-affa-885db325cfd5
@doc MultiComponentReactiveMixtureProject.DMS_Info_isothermal()

# ╔═╡ 415f6fa7-d5b5-40a2-806e-3d8a61541c2e
@doc MultiComponentReactiveMixtureProject.DMS_Info_thermal()

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
function ThermalDemo(dim; nref=nref, times=nothing, mfluxin = nothing, verbose="aen")
	if dim == 2
		times = isnothing(times) ? [0,50.0] : times
	else
		times = isnothing(times) ? [0,5.0] : times
	end

	grid, inb,irrb,outb,sb,catr =  grid_boundaries_regions(dim,nref=nref)

	data=ReactorData(
		dim=dim,
		nflowin = 7.4*ufac"mol/hr",
		nom_flux = 100.0*ufac"kW/m^2",		
		#nom_flux = 0.0*ufac"kW/m^2",
		dt_hf_irrad = (2.0, 10.0),
		dt_hf_enth = (2.0, 3.0),
		T_gas_in = 273.15 + 25,
		#T_gas_in = 273.15 + 600,
		Nu = 0.0,
		X0 = [0,0.5,0,0,0.5,0.0], # H2 / CO2 = 1/1
		inlet_boundaries=inb,
		irradiated_boundaries=irrb,
		outlet_boundaries=outb,
		side_boundaries=sb,
		catalyst_regions=catr,
		rhos=5.0*ufac"kg/m^3" # set solid density to low value to reduce thermal inertia of system
		)
	(;p,ip,Tamb,iT,iTw,iTp,ng,gni,X0)=data

	inival,sys = init_system(dim, grid, data)

	if dim == 2
		control = SolverControl(nothing, sys;)
	else
		control = SolverControl(;
        method_linear = KrylovJL_GMRES(),
        precon_linear = VoronoiFVM.factorizationstrategy(
			MKLPardisoLU(), NoBlock(), sys),
   		)
	end
		control.handle_exceptions=true
		control.Δt_min=1.0e-6
		control.Δu_opt=100
		
	solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="a",log=true)
	
	return solt,grid,sys,data
end

# ╔═╡ fac7a69d-5d65-43ca-9bf3-7d9d0c9f2583
# ╠═╡ skip_as_script = true
#=╠═╡
if RunSim
	solt,grid,sys,data=ThermalDemo(dim);
end;
  ╠═╡ =#

# ╔═╡ 927dccb1-832b-4e83-a011-0efa1b3e9ffb
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Initialisation and Solve
The simulation is setup as a transient simulation. An initialisation strategy is employed where different physics are enabled step by step once a stationary state is established. Initially, no heat is transported and no chemical reactions take place. 

1. Velocity field (mass flow is ramped up from 1-100 % in T=$(data.dt_mf) s)
2. Temperature field through enthalpy flux (ramped up from 0-100 % in T=$(data.dt_hf_enth) s)
3. Temperature field through irradiation b.c. (ramped up from 0-100 % in T=$(data.dt_hf_irrad) s)

The mass flow boundary condition into the reactor domain is "ramped up" starting from a low value and linearly increasing until the final value is reached. A time delay is given to let the flow stabilize. Once the flow field is established, heat transport is ramped up until a stable temperature field is established. Finally, the reactivity of the catalyst is "ramped up" until its final reactivity value is reached.
"""
  ╠═╡ =#

# ╔═╡ 1cc9d6c4-e2d6-4501-ae4d-d7568dee1e8f
#=╠═╡
plothistory(solt)
  ╠═╡ =#

# ╔═╡ 3207839f-48a9-49b6-9861-e5e74bc593a4
# ╠═╡ skip_as_script = true
#=╠═╡
MultiComponentReactiveMixtureProject.DMS_print_summary_ext(solt,grid,sys,data)
  ╠═╡ =#

# ╔═╡ 5d5ac33c-f738-4f9e-bcd2-efc43b638109
# ╠═╡ skip_as_script = true
#=╠═╡
let
	inflow_rate, outflow_rate, reaction_rate, stored_amount, I_in, I_out, I_reac = BoundaryFlows_Integrals(solt, sys, data)
	(;ng, gn, gni, iT, ip) = data

	k=gni[:CO]
	#k=iT

	if k in 1:ng
		name = gn[k]
	elseif k == ip
		name = "Total Mass"
	elseif k == iT
		name = "Enthalpy Flow"
	end
	
	@printf "%s In: %2.2e \t Out: %2.2e \t React: %2.2e \nIn - Out: %2.4e \nStorage tEnd -t0: %2.4e" name I_in[k] I_out[k] I_reac[k] (I_in+I_out-I_reac)[k] (stored_amount[end]-stored_amount[1])[k]

	vis=GridVisualizer(resolution=(600,300), xlabel="Time / s", ylabel="Molar flow / Total Moles")

	scalarplot!(vis, solt.t[2:end], stack(inflow_rate, dims=1)[:,k], label="Inflow rate")
	scalarplot!(vis, solt.t[2:end], stack(outflow_rate, dims=1)[:,k], label="Outflow rate", color=:red, clear=false)	
	scalarplot!(vis, solt.t[2:end], stack(-reaction_rate, dims=1)[:,k], label="Reaction rate",  color=:blue, clear=false)
	#scalarplot!(vis, solt.t[2:end], stack(stored_amount, dims=1)[:,k], label="Stored amount", color=:green, clear=false, )
	reveal(vis)

end
  ╠═╡ =#

# ╔═╡ d5b816c3-f5c6-4762-a610-5a2efb77d4ff
#=╠═╡
function runtests()
	inflow_rate, outflow_rate, reaction_rate, stored_amount, I_in, I_out, I_reac = BoundaryFlows_Integrals(solt, sys, data)
	(;ng, gn, gni, iT, ip, m) = data

	RR = stack(-reaction_rate, dims=1)[end, gni[:CO]] / m[gni[:CO]] / ufac"mol/hr"
	
    @test isapprox(RR, 1.0439224917035148)	
end;
  ╠═╡ =#

# ╔═╡ 98468f9e-6dee-4b0b-8421-d77ac33012cc
md"""
### Temperature
1) Porous frit + catalyst layer domain
2) Window inner surface
3) Bottom plate
"""

# ╔═╡ f798e27a-1d7f-40d0-9a36-e8f0f26899b6
#=╠═╡
@bind t Slider(solt.t,show_value=true,default=solt.t[end])
  ╠═╡ =#

# ╔═╡ 5588790a-73d4-435d-950f-515ae2de923c
# ╠═╡ skip_as_script = true
#=╠═╡
sol = solt(t);
  ╠═╡ =#

# ╔═╡ 994d4a87-3f27-4a51-b061-6111c3346d60
#=╠═╡
MultiComponentReactiveMixtureProject.DMS_print_summary(sol,grid,sys,data)
  ╠═╡ =#

# ╔═╡ 99b59260-7651-45d0-b364-4f86db9927f8
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;iT,iTw,iTp,irradiated_boundaries,outlet_boundaries)=data
	#vis=GridVisualizer(layout=(3,1), resolution=(680,900))
	vis=GridVisualizer(layout=(1,1), resolution=(680,300))
	scalarplot!(vis[1,1],grid, sol[iT,:] .- 273.15, zoom = 2.8, aspect=4.0, show=true)
end
  ╠═╡ =#

# ╔═╡ 58c0b05d-bb0e-4a3f-af05-71782040c8b9
if dim == 2
md"""
- (1,1): T-profile at r=0
- (2,1): T-profile at z=0
- (1,2): Window T-profile
- (2,2): Bottom Plate T-profile
"""
end

# ╔═╡ 8de4b22d-080c-486f-a6a9-41e8a5489966
# ╠═╡ show_logs = false
#=╠═╡
let
	if dim == 2
		(;iT,iTw,iTp,irradiated_boundaries,outlet_boundaries) = data
		vis=GridVisualizer(layout=(2,2), resolution=(680,600))
		function _2to1(a,b)
			a[1]=b[2]
		end
		_grid,_,_,_,_,_ = grid_boundaries_regions(dim,nref=nref)
		newreg = num_bfaceregions(grid) + 1
		bfacemask!(_grid, [0.0,0.0].*ufac"cm",[0.0,0.5].*ufac"cm",newreg)
		grid1D = subgrid(_grid, [newreg]; boundary = true, transform = _2to1)
		sol1D=view(sol[iT, :], grid1D)
		scalarplot!(vis[1,1],grid1D, sol1D .-273.15, label="Temperature along Y-axis", clear=false)
	
		function __2to1(a,b)
			a[1]=b[1]
		end
		grid1D = subgrid(grid, outlet_boundaries; boundary = true, transform = __2to1)
		sol1D=view(sol[iT, :], grid1D)
		scalarplot!(vis[2,1],grid1D, sol1D .-273.15, label="Temperature along X-axis", clear=false)
		
	    # window
		bgridw = subgrid(grid, irradiated_boundaries; boundary = true, transform = __2to1)
		bsolw=view(sol[iTw, :], bgridw)
		scalarplot!(vis[1,2],bgridw, bsolw.-273.15,)
		# bottom plate
		bgridp = subgrid(grid, outlet_boundaries; boundary = true, transform = __2to1)
		bsolp=view(sol[iTp, :], bgridp)
		scalarplot!(vis[2,2],bgridp, bsolp.-273.15,show=true)
	end
end
  ╠═╡ =#

# ╔═╡ c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
md"""
### Molar fractions
1) CO
2) CO2
"""

# ╔═╡ 111b1b1f-51a5-4069-a365-a713c92b79f4
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;ip,p,gn,gni) = data
	ng=ngas(data)
	if dim == 2
		vis=GridVisualizer(layout=(4,1), resolution=(680,900))
		scalarplot!(vis[1,1], grid, sol[gni[:CO],:], aspect = 4.0,zoom = 2.8) # CO
		scalarplot!(vis[2,1], grid, sol[gni[:CO2],:], aspect = 4.0,zoom = 2.8) # CO2
		scalarplot!(vis[3,1], grid, sol[gni[:N2],:], aspect = 4.0,zoom = 2.8) # N2

		cols = distinguishable_colors(ng)
		# plot species molar fractions along frit thickness (along y axis)
		function _2to1(a,b)
			a[1]=b[2]
		end
		_grid,_,_,_,_,_ = grid_boundaries_regions(dim)
		bfacemask!(_grid, [3.0,0.0].*ufac"cm",[3.0,0.5].*ufac"cm",5)
	    grid1D = subgrid(_grid, [5]; boundary = true, transform = _2to1)
		for i=1:ng
			sol1D=view(sol[i, :], grid1D)
			scalarplot!(vis[4,1],grid1D, sol1D, label=gn[i], color=cols[i],clear=false)
		end
		reveal(vis)
	
	else
		vis=GridVisualizer(layout=(3,1), resolution=(400,1200), outlinealpha=0.0)
		scalarplot!(vis[1,1], grid, sol[1,:])
		scalarplot!(vis[2,1], grid, sol[2,:])
		scalarplot!(vis[3,1], grid, sol[3,:])
	end	
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ eb9dd385-c4be-42a2-8565-cf3cc9b2a078
md"""
### Flow field
1. Pressure
2. Density
3. Velocity X
4. Velocity Y
"""

# ╔═╡ de69f808-2618-4add-b092-522a1d7e0bb7
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;p,m,ip,iT,Tamb,mfluxin) = data
	ng = ngas(data)
	mmix = []
	for j in 1:length(sol[1,:])
		_mmix=0
		for i=1:ng
			_mmix += sol[i,j]*m[i]
		end
		push!(mmix, _mmix)
	end
	
	ps = sol[ip,:]
	Ts = sol[iT,:]
	#rho = @. ps * mmix /(ph"R"*T)
	rho = @. ps * mmix /(ph"R"*Ts)
	
	if dim == 2
		vis=GridVisualizer(layout=(4,1), resolution=(600,800))
		scalarplot!(vis[1,1], grid, ps, aspect=4.0, zoom=3.5) # Total pressure
		scalarplot!(vis[2,1], grid, rho, aspect=4.0, zoom=3.5) # Total Density
		nf = nodeflux(sys, sol)
		massflux = nf[:,ip,:]
		scalarplot!(vis[3,1], grid, massflux[1,:]./rho, aspect=4.0, zoom=3.5) # Velocity - X
		scalarplot!(vis[4,1], grid, massflux[2,:]./rho, aspect=4.0, zoom=3.5) # Velocity - Y
	else
		vis=GridVisualizer(layout=(2,1), resolution=(400,800), outlinealpha=0.0)
		scalarplot!(vis[1,1], grid, ps, title="Total Pressure")
		scalarplot!(vis[2,1], grid, rho, title="Total Density")
		
	end
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─d3278ac7-db94-4119-8efd-4dd18107e248
# ╟─b2791860-ae0f-412d-9082-bb2e27f990bc
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═4e05ab31-7729-4a4b-9c14-145118477715
# ╟─a078e1e1-c9cd-4d34-86d9-df4a052b6b96
# ╠═a1ea393e-f123-4ad0-affa-885db325cfd5
# ╠═415f6fa7-d5b5-40a2-806e-3d8a61541c2e
# ╟─927dccb1-832b-4e83-a011-0efa1b3e9ffb
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═fac7a69d-5d65-43ca-9bf3-7d9d0c9f2583
# ╠═5588790a-73d4-435d-950f-515ae2de923c
# ╠═1cc9d6c4-e2d6-4501-ae4d-d7568dee1e8f
# ╠═994d4a87-3f27-4a51-b061-6111c3346d60
# ╠═3207839f-48a9-49b6-9861-e5e74bc593a4
# ╠═5d5ac33c-f738-4f9e-bcd2-efc43b638109
# ╠═d5b816c3-f5c6-4762-a610-5a2efb77d4ff
# ╟─98468f9e-6dee-4b0b-8421-d77ac33012cc
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╟─99b59260-7651-45d0-b364-4f86db9927f8
# ╟─58c0b05d-bb0e-4a3f-af05-71782040c8b9
# ╟─8de4b22d-080c-486f-a6a9-41e8a5489966
# ╟─c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
# ╟─111b1b1f-51a5-4069-a365-a713c92b79f4
# ╟─eb9dd385-c4be-42a2-8565-cf3cc9b2a078
# ╟─de69f808-2618-4add-b092-522a1d7e0bb7
