### A Pluto.jl notebook ###
# v0.19.35

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
	using NLsolve, LinearSolve,Pardiso, ILUZero
	using StaticArrays

	using LessUnitful
	using DataFrames
	using PlutoVista, Plots
	using PlutoUI, Colors
	using FixedBed
	using Printf
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═╡ skip_as_script = true
#=╠═╡
PlutoUI.TableOfContents(title="Photo Catalytic (PC) Reactor")
  ╠═╡ =#

# ╔═╡ b2791860-ae0f-412d-9082-bb2e27f990bc
md"""
# Introduction
Demonstration notebook for the photo thermal catalytic reactor (PCR) model. Solve energy equation alongside multicomponent species transport. Include reactive gas mixture (CO2,H2,CO,CH4,H2O,N2) with variable physical properties and a Ni based catalyst described with kinetics from published literature.

Select problem dimension: $(@bind dim Select([2, 3], default=2))

Check the box to __start the simulation__: $(@bind RunSim PlutoUI.CheckBox(default=false))
"""

# ╔═╡ 4e05ab31-7729-4a4b-9c14-145118477715
# ╠═╡ skip_as_script = true
#=╠═╡
if dim == 3
	@bind xcut Slider(linspace(0,16,17)*ufac"cm",show_value=true,default=8*ufac"cm")
end
  ╠═╡ =#

# ╔═╡ bc811695-8394-4c35-8ad6-25856fa29183
function grid_boundaries_regions(dim)
	Ω_catalyst = 2
	W=16
	H=0.5
	if dim == 2
		Γ_bottom = 1
		#Γ_bottom_insulating = 7
		Γ_side = 2
		Γ_sym = 4		
		Γ_top_permeable = 5
		Γ_top_irradiated = 6

		W=W/2 # axisymmetry, half domain is sufficient
		R1=(0:1:(W-1))*ufac"cm"
		R2=((W-1):(9/10):W-1/10)*ufac"cm"
		R3=((W-1/10):(1/10):W).*ufac"cm"
		R=glue(R1,R2)
		R=glue(R,R3)
		Z=(0:H/10:H)*ufac"cm"
		grid=simplexgrid(R,Z)
		circular_symmetric!(grid)
	
		cellmask!(grid,[0,9/10*H].*ufac"cm",[W,H].*ufac"cm",Ω_catalyst) # catalyst layer	
		bfacemask!(grid, [0,H].*ufac"cm",[W-1,H].*ufac"cm",Γ_top_permeable)
		bfacemask!(grid, [0,H].*ufac"cm",[W-2,0.5].*ufac"cm",Γ_top_irradiated) 
		#bfacemask!(grid, [W-1/10,0].*ufac"cm",[W,0].*ufac"cm",Γ_bottom_insulating) 
		
		inb = [Γ_top_permeable,Γ_top_irradiated]
		irrb = [Γ_top_irradiated]
		outb = [Γ_bottom]
		sb = [Γ_side]
	else
		Γ_side_1 = 1 
		Γ_side_2 = 2
		Γ_side_3 = 3
		Γ_side_4 = 4		
		Γ_bottom = 5
		Γ_top_permeable = 7
		Γ_top_irradiated = 8

		X=(0:1:W)*ufac"cm"
		Y=(0:1:W)*ufac"cm"
		Z=(0:H/10:H)*ufac"cm"	
		grid=simplexgrid(X,Y,Z)
	
		# catalyst region
		cellmask!(grid,[0,0,9/10*H].*ufac"cm",[W,W,H].*ufac"cm",Ω_catalyst) # catalyst layer	
		bfacemask!(grid, [1,1,H].*ufac"cm",[W-1,W-1,H].*ufac"cm",Γ_top_permeable)
		bfacemask!(grid, [2,2,H].*ufac"cm",[W-2,W-2,H].*ufac"cm",Γ_top_irradiated)

		inb = [Γ_top_permeable,Γ_top_irradiated]
		irrb = [Γ_top_irradiated]
		outb = [Γ_bottom]
		sb = [Γ_side_1,Γ_side_2,Γ_side_3,Γ_side_4]
	end

	return grid, inb, irrb, outb, sb, [Ω_catalyst]
end;

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	grid, inb,irrb,outb,sb,catr =  grid_boundaries_regions(dim)
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

# ╔═╡ 0fadb9d2-1ccf-4d44-b748-b76d911784ca
md"""
## Overall Mass Continuity
Mixture mass flow (overall mass flow) through the pore space of the porous medium. Mixture mass averaged velocity is calculated from Darcy equation. The void fraction (porosity) is given by $\epsilon$.

```math
\begin{align}
	\frac{\partial \epsilon \rho}{\partial t} + \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\end{align}
```
"""

# ╔═╡ b94513c2-c94e-4bcb-9342-47ea48fbfd14
md"""
## Species Mass Continuity and Transport
```math
\begin{align}
	\frac{\partial \epsilon \rho_i}{\partial t} + \nabla \cdot \left( \vec \Phi_i + \rho_i \vec v \right ) - r_i(\varrho) &= 0 ~, \qquad i = 1 ... \nu \\
		\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= -\sum_{j=1 \atop j \neq i}^{\nu} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
		\sum_{i=1}^\nu x_i &= 1
\end{align}
```
"""

# ╔═╡ c886dd12-a90c-40ab-b9d0-32934c17baee
md"""
where $\rho$ is the (total) mixture density, $\vec v$ is the mass-averaged (barycentric)  mixture velocity calculated with the Darcy equation, $x_i$, $w_i$ and $M_i$ are the molar fraction, mass fraction and molar mass of species $i$ respectively, $\vec \Phi_i$ is the __diffusive__ mass flux of species $i$ ($\frac{\text{kg}}{\text{m}^2 \text{s}}$) and $r_i$ is the species mass volumetric source/sink ($\frac{\text{kg}}{\text{m}^3 \text{s}}$) of gas phase species $i$.

The __convective__ mass flux of species $i$ is the product of the superficial mean velocity with the partial mass density $\rho_i \vec v$. The combined __convective-diffusive__ species mass flux $\vec \Psi_i = \vec \Phi_i + \rho_i \vec v$ can be introduced as an auxiliary variable.
"""

# ╔═╡ 8f2549f4-b0a6-440f-af94-6880e0814dc2
md"""
## Thermal Energy Transport
"""

# ╔═╡ 05ef2813-194a-4ba8-81aa-4ec96f68b01c
md"""
The thermal energy equation considers convective-diffusive transport of thermal energy and is formulated with effective transport parameters that are a consequence of the quasi-homogeneous phase approach:

```math
\begin{align}
\frac{\partial ( \varepsilon c c_p + (1-\varepsilon) \rho_{\text s} c_{\text s}) T}{\partial t} - \nabla \cdot \left(\lambda_{\text{eff}} \nabla T - h_{\text{mix}} \rho \vec v \right) = 0
\end{align}
```
The term for the time derivative of the temperature describes the heat storage capacity of the media in the domain and is a superposition of the heat capacities for gas and solid phase weighted by their respective volume fraction.

Here $\lambda_{\text{eff}}$ is the effective thermal conductivity in the quasi-homogeneous approach, which is a (complicated) function of the thermal conductivities of the solid and gas phases.

The convective heat transport ("thermal drift") within the porous medium is a result of gas phase mass flow carrying enthalpy and is expressed via the sum  of species (mass) fluxes and species mass-specific enthalpies $\sum_i^n \vec \Psi_i h_i$.
A simplification is made under the assumption that $\sum_i^n\vec \Psi_i h_i \approx (\rho \vec v) h_{\text{mix}}$, leading to above thermal energy transport equation.
Here the mixture mass specific enthalpy is defined by $h_{\text{mix}}=\sum_i^n x_i h_i / M_{\text{mix}}$.

The heat release from chemical reactions is considered as part of the species enthalpies.
"""

# ╔═╡ 87d3a394-47b1-4c16-aa67-2ad1f16e48ff
md"""
As part of the initialisation strategy (see next section) the convective contribution of heat flux is ramped up at the same time with the heat transport boundary conditions after the flow field has been established.
"""

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
function ThermalDemo(dim; times=nothing, mfluxin = nothing, verbose="aen")
	if dim == 2
		times = isnothing(times) ? [0,25.0] : times
	else
		times = isnothing(times) ? [0,15.0] : times
	end

	grid, inb,irrb,outb,sb,catr =  grid_boundaries_regions(dim)
	
	data=ReactorData(
		is_reactive = false,
		#G_lamp = 70.0*ufac"kW/m^2",
		G_lamp = 0.0*ufac"kW/m^2",
		dt_hf_irrad = (30.0,31.0),
		T_gas_in = 273.15 +300,
		X0 = [0,0,0,0,0,1.0], # N2
		k_nat_conv = 0.0,
		inlet_boundaries=inb,
		irradiated_boundaries=irrb,
		outlet_boundaries=outb,
		side_boundaries=sb,
		catalyst_regions=catr,
		rhos=5.0*ufac"kg/m^3" # set solid density to low value to reduce thermal inertia of system
		)
	(;p,ip,Tamb,iT,iTw,iTp,ng,gni,X0)=data
	ng=ngas(data)

	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=FixedBed.DMS_flux,
							reaction=FixedBed.DMS_reaction,
							storage=FixedBed.DMS_storage,
							bcondition=FixedBed.PCR_bcond,
							bflux=FixedBed.PCR_bflux,
							bstorage=FixedBed.PCR_bstorage,
							boutflow=FixedBed.DMS_boutflow,
							outflowboundaries=outb,
							assembly=:edgewise
							)

	enable_species!(sys; species=collect(1:(ng+2))) # gas phase species xi, ptotal & T
	enable_boundary_species!(sys, iTw, irrb) # window temperature as boundary species in upper chamber
	enable_boundary_species!(sys, iTp, outb) # plate temperature as boundary species in lower chamber
	inival=unknowns(sys)

	inival[ip,:].=p
	inival[[iT,iTw,iTp],:] .= Tamb
		
	for i=1:ng
		inival[i,:] .= X0[i]
	end
	
	#nd_ids = unique(grid[CellNodes][:,grid[CellRegions] .== Ω_catalyst])
	catalyst_nodes = []
	for reg in catr
		catalyst_nodes = vcat(catalyst_nodes, unique(grid[CellNodes][:,grid[CellRegions] .== reg]) )
	end
		
	cat_vol = sum(nodevolumes(sys)[unique(catalyst_nodes)])

	data.lcat = data.mcat/cat_vol
	local Ain = 0.0
	for boundary in inb
		Ain += bareas(boundary,sys,grid)
	end
	data.mfluxin = data.mflowin / Ain

	if dim == 2
		control = SolverControl(nothing, sys;)
	else
		control = SolverControl(;
        method_linear = KrylovJL_GMRES(
            gmres_restart = 10,
            restart = true,
            itmax = 100,
        ),
        precon_linear = VoronoiFVM.factorizationstrategy(
			MKLPardisoLU(), NoBlock(), sys),
   		)
	end
		control.handle_exceptions=true
		control.Δu_opt=10_000.0
		
	solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="nae")
	
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

# ╔═╡ 0e86c197-32ae-4cff-8ab2-fe7847e5514a
#=╠═╡
let
	HeatFluxes_EB_I(solt,grid,sys,data)
end
  ╠═╡ =#

# ╔═╡ 5d5ac33c-f738-4f9e-bcd2-efc43b638109
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;m,ip,iT,gn,gni,poros,mflowin,nflowin,W0,T_gas_in,Tamb,X0,outlet_boundaries,inlet_boundaries,dt_mf,dt_hf_enth)=data
	ng=ngas(data)
	vis=GridVisualizer(resolution=(600,300), xlabel="Time / s", ylabel="Molar flow / Total Moles")
	
	tfact=TestFunctionFactory(sys)	
	tf_out=testfunction(tfact,inlet_boundaries,outlet_boundaries)
	tf_in=testfunction(tfact,outlet_boundaries,inlet_boundaries)
		
	inflow_rate=Float64[]
	inflow_rate_manual=Float64[]
	outflow_rate=Float64[]
	reaction_rate=Float64[]
	stored_amount=Float64[]

	#k=gni[:N2]
	k=iT
	for i=2:length(solt)
		m_ = 1
		W_ = 1
		#fac = k in 1:ng ? ufac"mol/hr" : ufac"kg/hr"
		if k in 1:ng
			m_ = m[k]
			W_ = W0[k]
			ifr=mflowin*W_*ramp(solt.t[i]; du=(0.1,1), dt=dt_mf)			
		elseif k == iT
			ifr=integrate(sys,tf_in,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])[iT]
			
			#ifr_manual=nflowin*(enthalpy_mix(data, T_gas_in, X0)-enthalpy_mix(data, Tamb, X0)) * ramp(solt.t[i]; du=(0.0,1), dt=dt_hf_enth)
			ifr_manual=nflowin*enthalpy_mix(data, T_gas_in, X0) * ramp(solt.t[i]; du=(0.0,1), dt=dt_hf_enth)
			push!(inflow_rate_manual,ifr_manual)		
		end
		ofr=integrate(sys,tf_out,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])
		push!(inflow_rate,ifr/m_)
		
		push!(outflow_rate,ofr[k]/m_)		
		rr = integrate(sys,sys.physics.reaction,solt[i])[k,2]
		amount = sum(integrate(sys,sys.physics.storage,solt[i]), dims=2)[k]
		push!(reaction_rate, rr/m_)
		push!(stored_amount, amount/m_)
   	end

	
	# integrals
	I_in=0.0
	I_out=0.0
	I_reac=0.0
	for i=1:length(solt)-1
		I_in+=inflow_rate[i]*(solt.t[i+1]-solt.t[i])
		I_out+=outflow_rate[i]*(solt.t[i+1]-solt.t[i])
		I_reac+=reaction_rate[i]*(solt.t[i+1]-solt.t[i])
	end
	if k in 1:ng
		name = gn[k]
	elseif k == ip
		name = "Total Mass"
	elseif k == iT
		name = "Enthalpy Flow"
	end
	
	@printf "%s In: %2.2e \t Out: %2.2e \t React: %2.2e \nIn - Out: %2.4e \nStorage tEnd -t0: %2.4e" name I_in I_out I_reac I_in+I_out-I_reac stored_amount[end]-stored_amount[1]

	scalarplot!(vis, solt.t[2:end], inflow_rate, label="Inflow rate")
	scalarplot!(vis, solt.t[2:end], inflow_rate_manual, label="Inflow MANUAL", color=:pink, clear=false)
	scalarplot!(vis, solt.t[2:end], outflow_rate, label="Outflow rate", color=:red, clear=false)	
	scalarplot!(vis, solt.t[2:end], -reaction_rate, label="Reaction rate",  color=:blue, clear=false)
	#scalarplot!(vis, solt.t[2:end], stored_amount, label="Stored amount", color=:green, clear=false, )
	reveal(vis)

end
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

# ╔═╡ d6a073e4-f4f6-4589-918f-20b61a780dad
#=╠═╡
let
	(;inlet_boundaries,outlet_boundaries)=data
	tfact=TestFunctionFactory(sys)
	tf_out=testfunction(tfact,inlet_boundaries,outlet_boundaries)
	#tf_out=testfunction(tfact,[1,3,4,5,6,7],[2])
	out=integrate(sys,tf_out,sol)
	#scalarplot(grid,tf_out)
end
  ╠═╡ =#

# ╔═╡ 994d4a87-3f27-4a51-b061-6111c3346d60
#=╠═╡
FixedBed.DMS_print_summary(sol,grid,sys,data)
  ╠═╡ =#

# ╔═╡ 3207839f-48a9-49b6-9861-e5e74bc593a4
# ╠═╡ skip_as_script = true
#=╠═╡
FixedBed.DMS_print_summary_ext(sol,sys,data)
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
		_grid,_,_,_,_,_ = grid_boundaries_regions(dim)
		bfacemask!(_grid, [0.0,0.0].*ufac"cm",[0.0,0.5].*ufac"cm",8)
		grid1D = subgrid(_grid, [8]; boundary = true, transform = _2to1)
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
		vis=GridVisualizer(layout=(3,1), resolution=(680,900))
		scalarplot!(vis[1,1], grid, sol[gni[:CO],:], aspect = 4.0,zoom = 2.8) # CO
		scalarplot!(vis[2,1], grid, sol[gni[:CO2],:], aspect = 4.0,zoom = 2.8) # CO2

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

# ╔═╡ 68ca72ae-3b24-4c09-ace1-5e340c8be3d4
function len(grid)
	coord = grid[Coordinates]
	L=0.0
	if dim == 1
		L=coord[end]
	elseif dim == 2
		L=coord[1,end]
	else
		L=coord[2,end]
	end
	L*ufac"m"
end

# ╔═╡ db77fca9-4118-4825-b023-262d4073b2dd
md"""
### Peclet Number
```math
\text{Pe}_L= \frac{L \vec v}{D}

```
"""

# ╔═╡ ae8c7993-a89f-438a-a72a-d4a0c9a8ce57
# ╠═╡ skip_as_script = true
#=╠═╡
let
	L=maximum(grid[Coordinates][dim,:])
	(;mfluxin,mmix0,p,Tamb) = data
	ng = ngas(data)
	rho0 = p*mmix0/(ph"R"*Tamb)
	v0 = mfluxin / rho0
	D = DiffCoeffs(mydata)
	Pe = L*v0 / minimum(D[D.>0])
end
  ╠═╡ =#

# ╔═╡ e7497364-75ef-4bd9-87ca-9a8c2d97064c
md"""
### Damköhler Number
```math
\text{Da}= \frac{k_{\text{react}}}{k_{\text{conv}}} = \frac{\dot r}{c_0} \tau = \frac{\dot r L}{c_0 v}

```
"""

# ╔═╡ e000c100-ee46-454e-b049-c1c29daa9a56
# ╠═╡ skip_as_script = true
#=╠═╡
let
	L=maximum(grid[Coordinates][dim,:])
	(;mfluxin,mmix0,lcat,p,X0,gni) = data
	T = 650 + 273.15
	c0 = p*X0/(ph"R"*T)
	rho0 = p*mmix0/(ph"R"*T)
	v0 = mfluxin / rho0	
	RR = -lcat*ri(data,T,p*X0)
	tau = L/v0
	Da = maximum(RR)/c0[gni[:H2]] * tau
end
  ╠═╡ =#

# ╔═╡ f6e54602-b63e-4ce5-b8d5-f626fbe5ae7a
md"""
### Reynolds Number
```math
\text{Re} = \frac{\rho v D}{\eta}  
```
where $\rho$ is the density of the fluid, ideal gas in this case, $v$ is the magnitufe of the velocity, $D$ is a representative length scale, here it is the mean pore diameter, and $\eta$ is the dynamic viscosity of the fluid.
"""

# ╔═╡ f690c19f-22e3-4428-bc1c-3ed7d1646e71
# ╠═╡ skip_as_script = true
#=╠═╡
let	
	(;mfluxin,mmix0,X0,p,Tamb,m,dp,Fluids) = data
	ng = ngas(data)
	rho0 = p*mmix0/(ph"R"*Tamb)
	v0 = mfluxin / rho0

	ηf, λf = dynvisc_thermcond_mix(data, Tamb, X0)
    
    cf = heatcap_mix(Fluids, Tamb, X0)
	
	Re = v0*rho0*dp/ηf # Reynolds number
	#Pr = cf*ηf/λf # Prandtl number
	#Pe = u0*ρf*cf*d/λf # Peclet
	#Re,Pr,Pe
end
  ╠═╡ =#

# ╔═╡ 0d507c5e-eb0f-4094-a30b-a4a82fd5c302
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

# ╔═╡ 5bbe72b2-2f80-4dae-9706-7ddb0b8b6dbe
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;dp,Tamb,p) = data
	σs = Dict(:H2 => 289*ufac"pm", :CH4 => 380*ufac"pm", :H2O => 265*ufac"pm", :N2 => 364*ufac"pm", :CO => 376*ufac"pm", :CO2 => 330*ufac"pm")

	Kn = ph"k_B"*Tamb/(sqrt(2)*pi*σs[:H2O]^2*p*dp)	
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─d3278ac7-db94-4119-8efd-4dd18107e248
# ╟─b2791860-ae0f-412d-9082-bb2e27f990bc
# ╟─a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═d6a073e4-f4f6-4589-918f-20b61a780dad
# ╠═4e05ab31-7729-4a4b-9c14-145118477715
# ╠═bc811695-8394-4c35-8ad6-25856fa29183
# ╟─a078e1e1-c9cd-4d34-86d9-df4a052b6b96
# ╟─0fadb9d2-1ccf-4d44-b748-b76d911784ca
# ╟─b94513c2-c94e-4bcb-9342-47ea48fbfd14
# ╟─c886dd12-a90c-40ab-b9d0-32934c17baee
# ╟─8f2549f4-b0a6-440f-af94-6880e0814dc2
# ╟─05ef2813-194a-4ba8-81aa-4ec96f68b01c
# ╟─87d3a394-47b1-4c16-aa67-2ad1f16e48ff
# ╟─927dccb1-832b-4e83-a011-0efa1b3e9ffb
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═fac7a69d-5d65-43ca-9bf3-7d9d0c9f2583
# ╠═5588790a-73d4-435d-950f-515ae2de923c
# ╠═0e86c197-32ae-4cff-8ab2-fe7847e5514a
# ╠═994d4a87-3f27-4a51-b061-6111c3346d60
# ╠═3207839f-48a9-49b6-9861-e5e74bc593a4
# ╟─5d5ac33c-f738-4f9e-bcd2-efc43b638109
# ╟─98468f9e-6dee-4b0b-8421-d77ac33012cc
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╟─99b59260-7651-45d0-b364-4f86db9927f8
# ╟─58c0b05d-bb0e-4a3f-af05-71782040c8b9
# ╟─8de4b22d-080c-486f-a6a9-41e8a5489966
# ╟─c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
# ╟─111b1b1f-51a5-4069-a365-a713c92b79f4
# ╟─eb9dd385-c4be-42a2-8565-cf3cc9b2a078
# ╟─de69f808-2618-4add-b092-522a1d7e0bb7
# ╟─68ca72ae-3b24-4c09-ace1-5e340c8be3d4
# ╟─db77fca9-4118-4825-b023-262d4073b2dd
# ╠═ae8c7993-a89f-438a-a72a-d4a0c9a8ce57
# ╟─e7497364-75ef-4bd9-87ca-9a8c2d97064c
# ╠═e000c100-ee46-454e-b049-c1c29daa9a56
# ╟─f6e54602-b63e-4ce5-b8d5-f626fbe5ae7a
# ╠═f690c19f-22e3-4428-bc1c-3ed7d1646e71
# ╟─0d507c5e-eb0f-4094-a30b-a4a82fd5c302
# ╠═5bbe72b2-2f80-4dae-9706-7ddb0b8b6dbe