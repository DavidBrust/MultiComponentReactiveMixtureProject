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
	using Revise, Test
	using VoronoiFVM, ExtendableGrids, GridVisualize, SparseArrays
	using LinearSolve, Pardiso, ExtendableSparse
	
	using LessUnitful
	using PlutoUI, PlutoVista, Plots
	using CSV, Tables, Dates, Printf
	using StaticArrays
	using MultiComponentReactiveMixtureProject
	
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

Select problem dimension: $(@bind dim Select([2,3], default=2))

Select grid refinement level: $(@bind nref Select([0,1,2,3], default=0))

Check the box to __start the simulation__: $(@bind RunSim PlutoUI.CheckBox(default=false))
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
	grid, inb,irrb,outb,sb,catr,accr =  PTR_grid_boundaries_regions(dim, nref=nref)
	if dim == 2
		gridplot(grid, resolution=(660,300), aspect=4.0, zoom=2.8)
	else
		gridplot(grid,  xplane=xcut, resolution=(660,460), zoom=1.8, )
	end
end
  ╠═╡ =#

# ╔═╡ a1ea393e-f123-4ad0-affa-885db325cfd5
@doc MultiComponentReactiveMixtureProject.DMS_Info_isothermal()

# ╔═╡ 415f6fa7-d5b5-40a2-806e-3d8a61541c2e
@doc MultiComponentReactiveMixtureProject.DMS_Info_thermal()

# ╔═╡ 1638178e-840b-4abe-9f46-8b0bbe3d606a
Wolf_rWGS

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
function ThermalDemo(dim; nref=nref)

	grid, inb, irrb, outb, sb, catr, permr =  PTR_grid_boundaries_regions(dim, nref=nref)

	data = ReactorData(
		dim=dim,
		kinpar=Wolf_rWGS,
		#kinpar=XuFroment,
		p = 3.5*ufac"bar",
		#nflowin = 7.4*ufac"mol/hr",
		nflowin = 3.5*7.4*ufac"mol/hr",
		nom_flux = 70.0*ufac"kW/m^2",
		mcat = 3000*ufac"mg",
		dt_hf_irrad = (2.0, 10.0),
		dt_hf_enth = (2.0, 3.0),
		T_gas_in = 273.15 + 25,
		
		X0 = [0,0.5,0,0,0.5,0.0], # H2 / CO2 = 1/1
		inlet_boundaries=inb,
		irradiated_boundaries=irrb,
		outlet_boundaries=outb,
		side_boundaries=sb,
		catalyst_regions=catr,
		permeable_regions=permr,

		include_dpdt=true
	)
	
	inival,sys = PTR_init_system(dim, grid, data)

	if dim == 2
		times = [0,1000.0]
		control = SolverControl(nothing, sys;)
	elseif dim == 3
		times = [0,200.0]
		control = SolverControl(
			GMRESIteration(MKLPardisoLU(), EquationBlock()),
			sys
		)
	end
	control.handle_exceptions=true
	control.Δu_opt=100
	control.Δt_max=100
		
	solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="aen",log=true)
	
	return solt,grid,sys,data
end

# ╔═╡ fac7a69d-5d65-43ca-9bf3-7d9d0c9f2583
if RunSim
	solt,grid,sys,data=ThermalDemo(dim);
end;

# ╔═╡ 927dccb1-832b-4e83-a011-0efa1b3e9ffb
md"""
# Initialisation and Solve
The simulation is setup as a transient simulation. An initialisation strategy is employed where different physics are enabled step by step once a stationary state is established. Initially, no heat is transported and no chemical reactions take place. 

1. Velocity field (mass flow is ramped up from 1-100 % in T=$(data.dt_mf) s)
2. Temperature field through enthalpy flux (ramped up from 0-100 % in T=$(data.dt_hf_enth) s)
3. Temperature field through irradiation b.c. (ramped up from 0-100 % in T=$(data.dt_hf_irrad) s)

The mass flow boundary condition into the reactor domain is "ramped up" starting from a low value and linearly increasing until the final value is reached. A time delay is given to let the flow stabilize. Once the flow field is established, heat transport is ramped up until a stable temperature field is established. Finally, the reactivity of the catalyst is "ramped up" until its final reactivity value is reached.
"""

# ╔═╡ 1cc9d6c4-e2d6-4501-ae4d-d7568dee1e8f
plothistory(solt)

# ╔═╡ 3207839f-48a9-49b6-9861-e5e74bc593a4
# ╠═╡ skip_as_script = true
#=╠═╡
MultiComponentReactiveMixtureProject.Print_summary_ext(solt,grid,sys,data)
  ╠═╡ =#

# ╔═╡ 5d5ac33c-f738-4f9e-bcd2-efc43b638109
# ╠═╡ skip_as_script = true
#=╠═╡
let
	inflow_rate, outflow_rate, reaction_rate, stored_amount, I_in, I_out, I_reac = BoundaryFlows_Integrals(solt, sys, data)
	(;ng, gn, gni, iT, ip) = data

	#k=gni[:H2]
	#k=gni[:CO]
	#k=iT
	k=ip

	if k in 1:ng
		name = gn[k]
	elseif k == ip
		name = "Total Mass"
	elseif k == iT
		name = "Enthalpy Flow"
	end
	
	@printf "%s In: %2.2e \t Out: %2.2e \t React: %2.2e \nIn - Out: %2.4e \nStorage tEnd -t0: %2.4e" name I_in[k] I_out[k] I_reac[k] (I_in+I_out-I_reac)[k] (stored_amount[end]-stored_amount[1])[k]

	vis=GridVisualizer(resolution=(500,600), layout=(3,1), xlabel="Time / s", ylabel="Inflow / Outflow / Reaction Rate")

	function plot_flows!(k,vis)
		scalarplot!(vis, solt.t[2:end], stack(inflow_rate, dims=1)[:,k], label="Inflow rate")
		scalarplot!(vis, solt.t[2:end], stack(outflow_rate, dims=1)[:,k], label="Outflow rate", color=:red, clear=false)	
		scalarplot!(vis, solt.t[2:end], stack(-reaction_rate, dims=1)[:,k], label="Reaction rate",  color=:blue, clear=false)
		#scalarplot!(vis, solt.t[2:end], stack(stored_amount, dims=1)[:,k], label="Stored amount", color=:green, clear=false, )
	end

	plot_flows!(ip,vis[1,1])
	plot_flows!(gni[:H2],vis[2,1])
	plot_flows!(gni[:CO],vis[3,1])
	
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ 560ad300-42fc-4528-a3ec-95bcd66cdbce
md"""
# Stationary solution
"""

# ╔═╡ 70cdb28c-4b23-4ea4-8cd4-5eb97a3b930a
function run_ss(solt,sys)	
	VoronoiFVM.solve(
		sys;
		time = solt.t[end],
		inival=solt(solt.t[end]),
	)
end

# ╔═╡ e148f083-4d4e-4fe8-960d-bccd00689c9b
sol_ss = run_ss(solt,sys);

# ╔═╡ e6828f65-fc35-4e2e-aedd-324ccfe4a22c
function write_sol(sol; desc="")
	
	_t = now()
	tm = "$(hour(_t))_$(minute(_t))_$(second(_t))"
	path = "../../data/out/$(Date(_t))/"
	fn = tm*"_sol_$(desc).csv"
	try
        mkpath(path)
    catch e
        println("Directory " * path * " already exists.")
    end
	CSV.write(path*fn, Tables.table(sol), decimal=',', delim=";")
end

# ╔═╡ f99203e7-e53e-4109-b6ff-7fb87d290324
#write_sol(solt(3.0), desc="include_dpdt=$(data.include_dpdt)")

# ╔═╡ dbb6346c-e08a-4ad0-a985-3052272cf6c7
function Test_RR(sol_ss, sys, data)
	(;gni, m) = data
	
	inflow_rate, outflow_rate, reaction_rate, = BoundaryFlows_Integrals(sol_ss, sys, data)

	return -reaction_rate[gni[:CO]] / m[gni[:CO]] / ufac"mol/hr"
end

# ╔═╡ 380c74fb-66c4-43fb-a3f5-9c942b13fa0d
if dim == 2
	@test isapprox(Test_RR(sol_ss, sys, data), 1.0401674474564733)
elseif dim == 3
	@test isapprox(Test_RR(sol_ss, sys, data), 0.7774951984340692)
end

# ╔═╡ 98468f9e-6dee-4b0b-8421-d77ac33012cc
md"""
### Temperature
1) Porous frit + catalyst layer domain
2) Window inner surface
3) Bottom plate
"""

# ╔═╡ 58c0b05d-bb0e-4a3f-af05-71782040c8b9
if dim == 2
md"""
- (1,1): T-profile at r=0
- (2,1): T-profile at z=0
- (1,2): Window T-profile
- (2,2): Bottom Plate T-profile
"""
end

# ╔═╡ c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
md"""
### Molar fractions
1) CO
2) CO2
3) CH4
"""

# ╔═╡ eb9dd385-c4be-42a2-8565-cf3cc9b2a078
md"""
### Flow field
1. Pressure
2. Density
3. Velocity X
4. Velocity Y
"""

# ╔═╡ f798e27a-1d7f-40d0-9a36-e8f0f26899b6
@bind t Slider(solt.t,show_value=true,default=solt.t[end])

# ╔═╡ 5588790a-73d4-435d-950f-515ae2de923c
sol = solt(t);

# ╔═╡ 994d4a87-3f27-4a51-b061-6111c3346d60
MultiComponentReactiveMixtureProject.Print_summary(sol,grid,sys,data)

# ╔═╡ 99b59260-7651-45d0-b364-4f86db9927f8
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;iT,iTw,iTp,irradiated_boundaries,outlet_boundaries)=data
	#vis=GridVisualizer(layout=(3,1), resolution=(680,900))
	vis=GridVisualizer(layout=(1,1), resolution=(680,300))
	scalarplot!(vis[1,1],grid, sol_ss[iT,:] .- 273.15, zoom = 2.8, aspect=4.0, show=true)
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

# ╔═╡ c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
md"""
### Molar fractions
1) CO
2) CO2
"""

# ╔═╡ eb9dd385-c4be-42a2-8565-cf3cc9b2a078
md"""
### Flow field
1. Pressure
2. Density
3. Velocity X
4. Velocity Y
"""

# ╔═╡ 107b390f-f9e6-4879-89a7-ec1373bafb52
md"""
### Source term in Enthalpy Eq
Visualize distribution of magnitude of source term from $\partial p / \partial t$ [W/m³]:
"""

# ╔═╡ f798e27a-1d7f-40d0-9a36-e8f0f26899b6
@bind t Slider(solt.t,show_value=true,default=solt.t[end])

# ╔═╡ b42ce84e-9f97-488a-9311-24c809437623
sol = solt(t);

# ╔═╡ 994d4a87-3f27-4a51-b061-6111c3346d60
MultiComponentReactiveMixtureProject.Print_summary(sol,grid,sys,data)

# ╔═╡ 8de4b22d-080c-486f-a6a9-41e8a5489966
# ╠═╡ show_logs = false
let
	if dim == 2
		(;iT,iTw,iTp,irradiated_boundaries,outlet_boundaries) = data
		vis=GridVisualizer(layout=(2,2), resolution=(680,600))
		function _2to1(a,b)
			a[1]=b[2]
		end
		_grid,_,_,_,_,_ = PTR_grid_boundaries_regions(dim,nref=nref)
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

# ╔═╡ 111b1b1f-51a5-4069-a365-a713c92b79f4
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;ip,p,gn,gni) = data
	ng=ngas(data)
	if dim == 2
		#vis=GridVisualizer(layout=(4,1), resolution=(680,900))
		vis=GridVisualizer(layout=(3,1), resolution=(680,900))
		scalarplot!(vis[1,1], grid, sol[gni[:CO],:], aspect = 4.0,zoom = 2.8) # CO
		scalarplot!(vis[2,1], grid, sol[gni[:CO2],:], aspect = 4.0,zoom = 2.8) # CO2
		#scalarplot!(vis[3,1], grid, sol[gni[:N2],:], aspect = 4.0,zoom = 2.8) # N2

		cols = distinguishable_colors(ng)
		# plot species molar fractions along frit thickness (along y axis)
		function _2to1(a,b)
			a[1]=b[2]
		end
		_grid,_,_,_,_,_ = PTR_grid_boundaries_regions(dim,nref=nref)
		n_max_reg = grid[NumBFaceRegions] + 1
		bfacemask!(_grid, [3.0,0.0].*ufac"cm",[3.0,0.5].*ufac"cm", n_max_reg)
	    grid1D = subgrid(_grid, [n_max_reg]; boundary = true, transform = _2to1)
		for i=1:ng
			sol1D=view(sol[i, :], grid1D)
			#scalarplot!(vis[4,1],grid1D, sol1D, label=gn[i], color=cols[i],clear=false)
			scalarplot!(vis[3,1],grid1D, sol1D, label=gn[i], color=cols[i],clear=false)
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

# ╔═╡ 5547b97e-5adf-48ec-9fb9-55d54c1503a4
let
	(;idpdt, include_dpdt) = data
	if include_dpdt
		vis=GridVisualizer(resolution=(680,300))
		scalarplot!(vis,grid, sol[idpdt,:], zoom = 1.5, aspect=4.0, show=true)
	end
end

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

# ╔═╡ bcaae53b-d58b-4e36-9b79-471b02acaea6
md"""
# Auxiliary
"""

# ╔═╡ b519541f-9ccd-4032-bc51-9a1abaecadba
md"""
### Peclet Number
Describes the ratio of advective transport to diffusive transport. Can be formulated  for mass- and heat transfer, $\text{Pe}^{\text m}$ and $\text{Pe}^{\text h}$, respectively.

```math
	\text{Pe}^{\text m}= \frac{L \vec v}{D}
```

```math
	\text{Pe}^{\text h}= \frac{L \vec v \rho c_p}{\lambda}
```

"""


# ╔═╡ 57b8f297-7616-41f3-bef3-5237c6cfded6
md"""
### Knudsen Number
```math
\text{Kn} = \frac{\lambda}{l} = \frac{k_{\text B}T}{\sqrt 2 \pi \sigma^2pl}
```
where $\lambda$ is the mean free path length of the fluid, ideal gas in this case and $l$ is the pore diameter of the porous medium. Determine the mean free path length for ideal gas from kinetic gas theory and kinetic collision cross-sections (diameter).

-  $\text{Kn} < 0.01$: Continuum flow
-  $0.01 <\text{Kn} < 0.1$: Slip flow
-  $0.1 < \text{Kn} < 10$: Transitional flow
-  $\text{Kn} > 10$: Free molecular flow
"""

# ╔═╡ eaaf24b5-fabd-4363-8dbf-ebcc1e6416d1
function RePrPeKn(T, p, data)

	(;mfluxin, T_gas_in, mmix0, Fluids, X0, dp) = data
	ng = ngas(data)

	#T = T_gas_in
    dens = p*mmix0/(ph"R"*T)
	
	#ρf = density_idealgas(Fluids, T, p, x)
	dyn_visc, therm_cond = dynvisc_thermcond_mix(data, T, X0)

	u0 = mfluxin/dens
    heat_cap = heatcap_mix(Fluids, T, X0)
	
	Re = u0*dens*dp/dyn_visc # Reynolds number
	Pr = heat_cap*dyn_visc/(therm_cond*mmix0) # Prandtl number

	L = 5*ufac"mm" # length scale of reactor in main flow direction
	Pe_h = L*u0*dens*heat_cap/(mmix0*therm_cond) # Peclet heat transport

	D = zeros(Float64,ng,ng,)
	D_matrix!(D, T, p, data)	
	
	Pe_m = L*u0/maximum(D)
	
	σs = Dict(:H2 => 289*ufac"pm", :CH4 => 380*ufac"pm", :H2O => 265*ufac"pm", :N2 => 364*ufac"pm", :CO => 376*ufac"pm", :CO2 => 330*ufac"pm")
	Kn = ph"k_B"*T/(sqrt(2)*pi*σs[:H2O]^2*p*dp)

	
	Re, Pr, Pe_h, Pe_m, Kn	
end

# ╔═╡ 1196e9ed-024a-4469-95cf-a8622ecaf413
Re, Pr, Pe_h, Pe_m, Kn = RePrPeKn(600+273.15, 1*ufac"bar", data)

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─d3278ac7-db94-4119-8efd-4dd18107e248
# ╟─b2791860-ae0f-412d-9082-bb2e27f990bc
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═4e05ab31-7729-4a4b-9c14-145118477715
# ╠═a1ea393e-f123-4ad0-affa-885db325cfd5
# ╠═415f6fa7-d5b5-40a2-806e-3d8a61541c2e
# ╠═1638178e-840b-4abe-9f46-8b0bbe3d606a
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═fac7a69d-5d65-43ca-9bf3-7d9d0c9f2583
# ╠═b42ce84e-9f97-488a-9311-24c809437623
# ╟─927dccb1-832b-4e83-a011-0efa1b3e9ffb
# ╠═1cc9d6c4-e2d6-4501-ae4d-d7568dee1e8f
# ╠═994d4a87-3f27-4a51-b061-6111c3346d60
# ╠═3207839f-48a9-49b6-9861-e5e74bc593a4
# ╟─5d5ac33c-f738-4f9e-bcd2-efc43b638109
# ╟─560ad300-42fc-4528-a3ec-95bcd66cdbce
# ╠═70cdb28c-4b23-4ea4-8cd4-5eb97a3b930a
# ╠═e148f083-4d4e-4fe8-960d-bccd00689c9b
# ╠═e6828f65-fc35-4e2e-aedd-324ccfe4a22c
# ╠═f99203e7-e53e-4109-b6ff-7fb87d290324
# ╠═dbb6346c-e08a-4ad0-a985-3052272cf6c7
# ╠═380c74fb-66c4-43fb-a3f5-9c942b13fa0d
# ╟─98468f9e-6dee-4b0b-8421-d77ac33012cc
# ╠═99b59260-7651-45d0-b364-4f86db9927f8
# ╟─58c0b05d-bb0e-4a3f-af05-71782040c8b9
# ╟─8de4b22d-080c-486f-a6a9-41e8a5489966
# ╟─c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
# ╠═111b1b1f-51a5-4069-a365-a713c92b79f4
# ╟─eb9dd385-c4be-42a2-8565-cf3cc9b2a078
# ╟─107b390f-f9e6-4879-89a7-ec1373bafb52
# ╠═5547b97e-5adf-48ec-9fb9-55d54c1503a4
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╟─de69f808-2618-4add-b092-522a1d7e0bb7
# ╟─bcaae53b-d58b-4e36-9b79-471b02acaea6
# ╠═1196e9ed-024a-4469-95cf-a8622ecaf413
# ╟─b519541f-9ccd-4032-bc51-9a1abaecadba
# ╟─57b8f297-7616-41f3-bef3-5237c6cfded6
# ╠═eaaf24b5-fabd-4363-8dbf-ebcc1e6416d1
