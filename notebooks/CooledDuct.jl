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

# ╔═╡ 349e7220-dc69-11ee-13d2-8f95e6ee5c96
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using PlutoUI
	using LessUnitful
	using StaticArrays
	using ExtendableSparse, Pardiso, LinearSolve
	using VoronoiFVM
	using ExtendableGrids, GridVisualize, CairoMakie, ColorSchemes, Plots
	using CSV, DataFrames
	using Revise
	using MultiComponentReactiveMixtureProject
	GridVisualize.default_plotter!(CairoMakie)
end;

# ╔═╡ 0f102f06-3ff3-4bcc-8892-8d9190a87849
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(aside=true)
  ╠═╡ =#

# ╔═╡ d14462c6-f63b-4a61-a1d9-4bcdb8e30e3d
md"""
Supplementary notebook for:$br
__Transport of heat and mass for reactive gas mixtures in porous media: modeling and application__.

This notebook is used to demonstrate the capabilites of the model for coupled heat- and mass transfer by reproducing simulation results for an example case presented in [1].
"""

# ╔═╡ 38c3ddb4-b44a-4981-9000-0a1d303bd9ac
md"""
# Introduction
The non-isothermal modeling framework is applied to a binary gas mixture (air) flowing through a rectangular duct that is filled with a porous material corresponding to the test case presented in [1].
"""

# ╔═╡ c8a3fd5c-91ae-407e-a564-caf4c3665fcc
# ╠═╡ skip_as_script = true
#=╠═╡
@doc MultiComponentReactiveMixtureProject.DMS_Info_isothermal()
  ╠═╡ =#

# ╔═╡ 192627ce-11f2-41af-850a-76adcba24168
# ╠═╡ skip_as_script = true
#=╠═╡
@doc MultiComponentReactiveMixtureProject.DMS_Info_thermal()
  ╠═╡ =#

# ╔═╡ b3a6fe03-be46-4159-96ab-477a42d0eec5
md"""
# Simulation setup
The simulation workflow comprises the following steps:
- Setup the modeling domain
- Define simulation parameters
- Define boundary conditions
- Compose the system and solve it
- Visualize the solution
"""

# ╔═╡ c3dbf8b3-fc5f-44ff-be2c-ca4200f5bd6c
md"""
## Modeling domain
"""

# ╔═╡ 9da82e5a-169e-4758-8e50-73fa8793aa80
begin
	const L = 80ufac"mm"
	const H = 10ufac"mm"
	const mflowin = 1.0e-3ufac"kg/m^3"
	#const mflowin = 5.0e-4ufac"kg/m^3"
end;

# ╔═╡ 138b4d12-af0f-4a1c-a4aa-4f55e6038664
begin
	const Γ_bottom = 1
	const Γ_right = 2
	const Γ_top = 3
	const Γ_left = 4
	const Γ_left_inflow = 5
	const Γ_right_outflow = 6

	const Γ_inlet_wall_bottom = 7
	const Γ_inlet_wall_top = 8
	
	const Ω_free = 2
	const Ω_porous = 3
	
end;

# ╔═╡ 6939978d-9590-407b-80dc-54721c3f672d
md"""
Select grid refinement level: $(@bind nref Select([0,1,2,3], default=0))

"""

# ╔═╡ 7be1c79e-08d8-493c-bce0-114c0c003dd7
function grid_2D(;nref=0, L=L, H=H, efix=0.01ufac"mm")

	nx = 20*2^nref
	ny = 10*2^nref

	lmin = L/nx/10
    lmax = L/nx

	hmin = H/ny/10
	hmax = H/ny

	#X_ = collect(-L/2:L/2/nx:0) # uniform spacing, inlet region
	
	#XLeft = collect(0:L/nx:(3/4*L)) # uniform spacing
	XLeft = collect(-L/2:L/nx:(3/4*L)) # uniform spacing
    XRight = geomspace(3/4*L, L, lmax, lmin)
	#XLeft = glue(X_, XLeft)
	X = glue(XLeft, XRight)
	#X = collect(0:L/nx:L)

	Y_ = [-efix, 0.0]
	#YBottom = geomspace(0, H/2, hmin, hmax)
	YBottom = collect(0:H/2/ny:H/2)
	YBottom = glue(Y_, YBottom)
	
	#YTop = geomspace(H/2, H, hmax, hmin)
	YTop = collect(H/2:H/2/ny:H)
	Y_ = [H, H + efix]
	YTop = glue(YTop, Y_)
	Y = glue(YBottom, YTop)
    #Y = collect(0:H/ny:H)
	
    grid = simplexgrid(X, Y)
	cellmask!(grid, [-L/2,0], [0,H], Ω_free) # free space inlet region
	cellmask!(grid, [0,0], [L,H], Ω_porous) # porous region
	bfacemask!(grid, [-L/2,0], [-L/2,H], Γ_left_inflow)
	bfacemask!(grid, [L,0], [L,H], Γ_right_outflow)
	bfacemask!(grid, [-L/2,-efix], [0,-efix], Γ_inlet_wall_bottom)
	bfacemask!(grid, [-L/2,H + efix], [0,H + efix], Γ_inlet_wall_top)

	grid
end

# ╔═╡ 869652e5-f15e-43d4-8fbc-724e866892b6
grid = grid_2D(nref=nref)

# ╔═╡ b55537bf-9982-4997-8a2a-1972127bdd86
# ╠═╡ skip_as_script = true
#=╠═╡
let
	vis = GridVisualizer(resolution=(660,250))
	gridplot!(vis, grid, linewidth=0.2, aspect = 2)
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ e9cb07eb-cfbb-4802-bc7f-6de7a6ad8ac6
md"""
## Simulation parameters
"""

# ╔═╡ 81099644-bd82-43b7-a433-b0415d8776f0
function fhp(data,x,y)

	(;mflowin, T_gas_in, p, mmix0) = data

	c0 = p/(ph"R"*T_gas_in)
	rho0 = mmix0 * c0
	v0 = mflowin / (H*1.0ufac"m") / rho0
	
	yh=y/H
	# 2D velocitx vector, x->hp , y -> 0
	return 6*v0*yh*(1.0-yh), 0
	
end

# ╔═╡ 3954d223-efb2-4aed-8560-263efe5a480e
function fhp_(f,u,edge,data)

	(;mflowin, T_gas_in, p, mmix0) = data

	c0 = p/(ph"R"*T_gas_in)
	rho0 = mmix0 * c0
	v0 = mflowin / (H*1.0ufac"m") / rho0
	
	x,y = edge.coord[:,edge.index]
	
	yh=y/H

	
	return v0,0
	
	# 2D velocitx vector, x component follow hp profile, y component is 0
	#return 6*v*yh*(1.0-yh),0
end

# ╔═╡ f008f30f-0137-4c54-8d4e-a6f589c4a952
md"""
Setup data structure corresponding to Air:
1) Oxygen (O2)
2) Nitrogen (N2)
"""

# ╔═╡ e7ca4902-0e14-48ca-bcc6-96b06c85a39d
Air = KinData{}(;
    ng = 2,
    gnames = [:O2, :N2],
    Fluids = [O2, N2],
    gn = Dict(1:2 .=> [:O2, :N2]),
    gni = Dict(value => key for (key, value) in Dict(1:2 .=> [:O2, :N2])),
	nr = 0,
	rnames = [],
)

# ╔═╡ 7f1d9cf8-7785-48c1-853c-74680188121f
data = ReactorData(
	dim = 2,
	inflow_boundaries = [Γ_left_inflow],
	outflow_boundaries = [Γ_right_outflow],
	permeable_regions = [Ω_free, Ω_porous],
	X0 = [0.21,0.79],
	#X0 = [0.0,1.0],
	mflowin = mflowin,
	mfluxin = mflowin/(H*1.0ufac"m"),
	kinpar = Air,
	Tamb = 50.0 + 273.15,
	T_gas_in = 100 + 273.15,
	
	p = 101.3*ufac"kPa",

	constant_properties = true,
	constant_species_viscosities = [
		dynvisc_gas(O2, 75.0+273.15),
		dynvisc_gas(N2, 75.0+273.15)
	],
	constant_species_thermal_conductivities = [
		thermcond_gas(O2, 75.0+273.15),
		thermcond_gas(N2, 75.0+273.15)
	],		
	solve_T_equation = true,
	is_reactive = false,
	include_Soret_Dufour = false,
	
	perm = [0.0,1.0,1.0]*1.23e-10*ufac"m^2",
	rhos = 5900.0*ufac"kg/m^3",
	cs = 500.0*ufac"J/(kg*K)",
	lambdas = 1.5*ufac"W/(m*K)",
	poros = [0.0,1.0,0.9],
	γ_τ = [1.0,1.0,1.1]

	# perm = 1.23e-10*ufac"m^2" * 1.0e6,
)

# ╔═╡ 0d3c46d6-d839-4969-8840-68d039e4ef9a
begin
	const evelo=edgevelocities(grid,(x,y) -> fhp(data,x,y))
	const bfvelo=bfacevelocities(grid,(x,y) -> fhp(data,x,y))
end;

# ╔═╡ bcfc3138-3cbb-4386-a33b-573b6c39caf9
md"""
## Boundary conditions
"""

# ╔═╡ 0c5e24c0-4aa5-44a2-b2fd-db78795485af
md"""
The modeling domain corresponds to a closed system, so no-flux boundary conditions for species mass apply on all domain boundaries. The homogeneous Neumann boundary conditions correspond to the default boundary conditon and thus nothing needs to be specified in particular.
"""

# ╔═╡ b5870a92-89e2-4eae-a79f-6033b1f3489e
function bcondition(f,u,bnode,data)
	(;p,ip,iT,inflow_boundaries,outflow_boundaries,dt_hf_enth,mfluxin,dt_mf,solve_T_equation,W0,X0,T_gas_in,Tamb,mmix0)=data

	ng = ngas(data)
	
	eps_ = 1/1e-4
	
	#boundary_robin!(f,u,bnode, species=iT,region=Γ_bottom,value=Tamb * eps_,factor=eps_)
	boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_bottom,value=Tamb)
	#boundary_robin!(f,u,bnode, species=iT,region=Γ_top,value=Tamb * eps_,factor=eps_)
	boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_top,value=Tamb)

	#boundary_dirichlet!(f,u,bnode, species=iT, region=[Γ_inlet_wall_bottom,Γ_inlet_wall_top],value=T_gas_in)
	#boundary_robin!(f,u,bnode, species=iT, region=[Γ_inlet_wall_bottom,Γ_inlet_wall_top],value=T_gas_in * eps_,factor=eps_)
	
	if bnode.region in inflow_boundaries
		r_mfluxin = mfluxin*ramp(bnode.time; du=(0.1,1), dt=dt_mf)
		
		f[ip] = -r_mfluxin # Neumann bc. for total mass flux
		for i=1:(ng-1)
			f[i] = -r_mfluxin*W0[i] # Neumann bc. for species mass flux
		end

	end

	if solve_T_equation
		if bnode.region in inflow_boundaries

			# heatflux from enthalpy inflow
			hin = zero(eltype(u))
			@inbounds for i=1:ng
				# hin += enthalpy_gas_thermal(Fluids[i],T_gas_in)
				hin += X0[i]*enthalpy_gas_thermal(data.Fluids[i],T_gas_in)
			end

			r_hf_enth = mfluxin/mmix0 * hin
			f[iT] = -r_hf_enth * ramp(bnode.time; du=(0.0,1), dt=dt_hf_enth)		
		end

	end

	
	# outflow boundary condition
	for boundary in outflow_boundaries
        boundary_dirichlet!(f,u,bnode, species=ip,region=boundary,value=p)
    end	
		
end

# ╔═╡ 13cfe122-eea8-4bbb-aaa0-5bcc74a247d1
md"""
# Solving and plotting
"""

# ╔═╡ 7b0a84b5-60d3-4dd9-89e9-29c88282cb25
md"""
## Transient solution
"""

# ╔═╡ 3a063481-d447-4bf5-9c49-ecde37a0fcea
md"""
Setup the system of equations as a transient system, solve the system and return the transient solution:
"""

# ╔═╡ bf036015-80a1-4b5a-9ddd-9bfc939979e0
md"""
## Stationary solution
"""

# ╔═╡ 9274b233-687d-471d-8bac-e14e9a0cb7c0
md"""
## Comparison with published literature
__Taken from [1]__: $br
Fig. 4. Temperature distribution (K) in the two-dimensional fluid domain partially
filled by a porous solid. Total mass flow rate of the gas mixture: 4E-3 kg/s.
"""

# ╔═╡ c9c931d9-0a2a-47de-9bee-493c472def48
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
$(LocalResource("../data/Goell2012/Goell2012_CooledDuct.png", :width => 600))
"""
  ╠═╡ =#

# ╔═╡ 665c5826-4516-4ffa-a9d8-def8e3985ddb
md"""
__This work__: $br
"""

# ╔═╡ 5abe4d49-f398-4d0a-8055-9c4b1014e74f
md"""
### Nusselt number

```math
	\mathrm{Nu}=\frac{2H}{T_{\mathrm w}-T_{\mathrm m}} \left. \frac{\partial T}{\partial y} \right|_{\mathrm w}
```
Where $T_{\mathrm m}$ is the mass flux weighted, mean temperature over the cross-section of the duct:

```math
	T_{\mathrm m} = \frac{\int_0^H (\rho v T) dy}{\int_0^H (\rho v) dy} 
```
"""

# ╔═╡ 65dbb492-4795-44ca-afcb-fb2a2c925d92
md"""
# References
1) __Göll, Stephan; Piesche, Manfred (2012)__: Multi-component gas transport in micro-porous domains: Multidimensional simulation at the macroscale. In: International Journal of Heat and Mass Transfer 55 (1-3), S. 480–487. DOI: 10.1016/j.ijheatmasstransfer.2011.09.049.
"""

# ╔═╡ e8425c71-666a-462e-9c4d-fc480810f922
md"""
# Function definitions
"""

# ╔═╡ 056119e3-ec74-4860-abb4-b75f6b16878d
function quad_trap(v,coord) 
    N = size(coord,1)
	int = 0.0
    for i=1:N-1
        xk = coord[i+1] - coord[i]
        int = int + (v[i]+v[i+1])/2*xk
    end
    int
end

# ╔═╡ 1e51701d-a893-4056-8336-a3772b85abe4
function setup_run_sim(grid, data)
	(;ng, ip, iT, Tamb, p, X0, inflow_boundaries, outflow_boundaries, permeable_regions) = data
	
	sys=VoronoiFVM.System(
		grid;
		data=data,
		flux=DMS_flux,
		reaction=DMS_reaction,
		storage=DMS_storage,
		bcondition=bcondition,
		boutflow=DMS_boutflow,
		outflowboundaries=outflow_boundaries,
		assembly=:edgewise,
		unknown_storage=:dense
	)

	#enable_species!(sys; species=collect(1:(ng+2))) # gas phase species xi, ptotal, T
	enable_species!(sys; species=collect(1:(ng+1)), regions=permeable_regions) # gas phase species xi, ptotal
	enable_species!(sys; species=iT) # T

	inival=unknowns(sys)
	inival .= zero(eltype(inival))
	
	for reg in permeable_regions
		perm_nodes = unique(grid[CellNodes][:,grid[CellRegions] .== reg]) 

		inival[ip,perm_nodes].=p
		for i=1:ng
			
			inival[i,perm_nodes] .= X0[i]
		end
	end
	
	inival[iT,:].=Tamb


	#control = SolverControl(
	#	GMRESIteration(MKLPardisoLU(), EquationBlock()),
	# 	sys
	#)

	
	control = SolverControl(nothing, sys;)
	control.handle_exceptions=true
	control.Δt_max=100.0
	control.Δt_grow=1.5
	control.Δu_opt=1000.0
	

	times=[0,1000.0]
	
	sol=VoronoiFVM.solve(sys;inival=inival,times,control,log=true,verbose="nae")
	return (sol,sys)
end

# ╔═╡ 035d4123-7092-4429-8cfd-1e5926e84493
# ╠═╡ skip_as_script = true
#=╠═╡
solt, sys = setup_run_sim(grid, data);
  ╠═╡ =#

# ╔═╡ 5a0900cc-df10-4176-b903-358b3e00415c
#=╠═╡
md"""
Move the slider to change $t$: $(if isa(solt, TransientSolution)
	@bind t PlutoUI.Slider(solt.t,show_value=false,default=solt.t[end])	
end)
"""
  ╠═╡ =#

# ╔═╡ 26bab6eb-7457-4fb7-b8e2-5148769891ff
#=╠═╡
sol = solt(t);
  ╠═╡ =#

# ╔═╡ ae6e4bb7-46e8-4f95-b337-0b4589c43cbf
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;iT,ip,gni) = data
	vis =GridVisualizer(layout=(3,1), resolution=(660,450), colorbarticks=4)
	x_O2_ = sol[gni[:O2],:]
	x_O2_ = x_O2_[x_O2_ .!= 0.0]
	x_N2_ = sol[gni[:N2],:]
	x_N2_ = x_N2_[x_N2_ .!= 0.0]
	
	scalarplot!(vis[1,1], grid, sol[gni[:O2],:], limits=(minimum(x_O2_), maximum(x_O2_)), title = "O2 molar fraction",)
	scalarplot!(vis[2,1], grid, sol[gni[:N2],:], limits=(minimum(x_N2_), maximum(x_N2_)), title = "N2 molar fraction",)
	scalarplot!(vis[3,1], grid, sol[iT,:], title = "Temperature")

	reveal(vis)
	
end
  ╠═╡ =#

# ╔═╡ 076b4a28-be0f-46f0-9857-e6f886c4b118
#=╠═╡
md"""
Plot the transient solution at time $t=$ $(round(t)):
"""
  ╠═╡ =#

# ╔═╡ 48366e85-cffb-4c5c-ac3a-807e47f858c7
function run_ss(solt,sys)	
	sol_steadystate = VoronoiFVM.solve(
		sys;
		time = solt.t[end],
		inival=solt(solt.t[end]),
		verbose="na"
	)
end

# ╔═╡ 9b13bc55-63eb-46bc-acea-f9aec90b340f
#=╠═╡
sol_ss = run_ss(solt,sys);
  ╠═╡ =#

# ╔═╡ 47c852e5-38eb-4c1f-8d60-d180e4826d05
#=╠═╡
nf = nodeflux(sys, sol_ss);
  ╠═╡ =#

# ╔═╡ 77f663a6-96e5-4a1c-843c-cd66fd1382b5
#=╠═╡
function Profiles_cross(sol, xpos, data)
	(;ip,iT) = data
	grid = grid_2D(nref=nref)

	N_MAX_REG = grid[NumBFaceRegions] + 1
	bfacemask!(grid, [xpos,0],[xpos,H],N_MAX_REG)

	grid_1D  = subgrid(grid, [N_MAX_REG], boundary=true, transform=(a,b)->a[1]=b[2]) # transform y coordinate of parent grid into x coordinate of subgrid
	
	T_profile_cross = view(sol[iT, :], grid_1D)

	MF = nf[:,ip,:]
	MF_x = MF[1,:]

	MF_profile_cross = view(MF[1,:], grid_1D)
			
	T_profile_cross, MF_profile_cross, grid_1D	
end
  ╠═╡ =#

# ╔═╡ 4acaadd4-6102-44f5-b602-00465bf3feca
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;iT,ip,gni,m) = data

	T = sol_ss[iT,:]
	p = sol_ss[ip,:]
	
	ng = ngas(data)
	mmix = []	
	for j in 1:length(sol_ss[1,:])
		_mmix=0
		for i=1:ng
			_mmix += sol_ss[i,j] * m[i]
		end
		push!(mmix, _mmix)
	end
	
	dens = @. p * mmix /(ph"R"*T)

	p_ = p[p .!= 0.0]
	dens_ = dens[dens .!= 0.0]
	
	#nf = nodeflux(sys, sol_ss)

	massflux = nf[:,ip,:]
	massflux_x_ = massflux[1,:]
	massflux_x_ = massflux_x_[dens .!= 0.0]
	
	vel_x = massflux[1,:]./dens
	vel_x_ = vel_x[dens .!= 0.0]
	vel_y = massflux[2,:]./dens
	vel_y_ = vel_y[dens .!= 0.0]
	
	vis =GridVisualizer(layout=(5,1), resolution=(600,750))
	
	scalarplot!(vis[1,1], grid, T, title = "Temperature / K")
	scalarplot!(vis[2,1], grid, p, limits=(minimum(p_),maximum(p_)), climits=(minimum(p_),maximum(p_)), title = "Pressure / Pa",)
	scalarplot!(vis[3,1], grid, dens, limits=(minimum(dens_),maximum(dens_)), title = "Density / kg m-3",)
		
	scalarplot!(vis[4,1], grid, massflux[1,:], limits=(minimum(massflux_x_),maximum(massflux_x_)), title = "Mass flux (X) / kg s-1 m-2")
	#scalarplot!(vis[4,1], grid, massflux[2,:], title = "Mass flux (Y) / kg s-1 m-2")
	
	scalarplot!(vis[5,1], grid, vel_x, limits=(minimum(vel_x_),maximum(vel_x_)), title = "Velocity (X) / m s-1")
	#scalarplot!(vis[6,1], grid, vel_y, title = "Velocity (Y) / m s-1")

	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ a21d9a07-de69-4884-8d7d-742413f9a95a
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;iT) = data
	T = sol_ss[iT,:]

	levels = 41

	colormap = ColorSchemes.jet1[20:90]

	res = (650,250)
	vis =GridVisualizer(layout=(1,1), resolution=res, aspect=1.2, levels=0, colormap=colormap, colorbar=:horizontal, colorbarticks=6)
	
	scalarplot!(vis[1,1], grid, T, levels=0, linewidth=0, colorlevels=levels+2, title = "Temperature / K")
	
	sc = reveal(vis)
	#fn = "../img/out/2024-06-07/Cooled_Duct_nref2.png"
	#GridVisualize.save(fn, sc)
end
  ╠═╡ =#

# ╔═╡ 5623ace0-4b62-4ec7-b54f-c110d38bc06b
#=╠═╡
function Nu(x, data)
	(;iT, Tamb) = data
		
	T_profile_cross,
	MF_profile_cross,
	grid_1D_profile_cross = Profiles_cross(sol_ss, x, data)

	y_coord = vec(grid_1D_profile_cross[Coordinates])
	
	MFlow = quad_trap(MF_profile_cross, y_coord)
	Tm = quad_trap(MF_profile_cross.*T_profile_cross, y_coord) / MFlow 

	# finite diff approx.
	dtdy_w = (T_profile_cross[end]-T_profile_cross[end-1])/(y_coord[end]-y_coord[end-1])

	Nu = 2*H/(Tamb - Tm) * dtdy_w
end
  ╠═╡ =#

# ╔═╡ 91980f61-4173-4840-a1ae-f1d743daa2c4
# ╠═╡ skip_as_script = true
#=╠═╡
Nu(60ufac"mm", data)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═349e7220-dc69-11ee-13d2-8f95e6ee5c96
# ╠═0f102f06-3ff3-4bcc-8892-8d9190a87849
# ╟─d14462c6-f63b-4a61-a1d9-4bcdb8e30e3d
# ╟─38c3ddb4-b44a-4981-9000-0a1d303bd9ac
# ╠═c8a3fd5c-91ae-407e-a564-caf4c3665fcc
# ╠═192627ce-11f2-41af-850a-76adcba24168
# ╟─b3a6fe03-be46-4159-96ab-477a42d0eec5
# ╟─c3dbf8b3-fc5f-44ff-be2c-ca4200f5bd6c
# ╠═9da82e5a-169e-4758-8e50-73fa8793aa80
# ╠═138b4d12-af0f-4a1c-a4aa-4f55e6038664
# ╠═869652e5-f15e-43d4-8fbc-724e866892b6
# ╠═b55537bf-9982-4997-8a2a-1972127bdd86
# ╟─6939978d-9590-407b-80dc-54721c3f672d
# ╠═7be1c79e-08d8-493c-bce0-114c0c003dd7
# ╟─e9cb07eb-cfbb-4802-bc7f-6de7a6ad8ac6
# ╠═81099644-bd82-43b7-a433-b0415d8776f0
# ╠═0d3c46d6-d839-4969-8840-68d039e4ef9a
# ╠═3954d223-efb2-4aed-8560-263efe5a480e
# ╟─f008f30f-0137-4c54-8d4e-a6f589c4a952
# ╠═e7ca4902-0e14-48ca-bcc6-96b06c85a39d
# ╠═7f1d9cf8-7785-48c1-853c-74680188121f
# ╟─bcfc3138-3cbb-4386-a33b-573b6c39caf9
# ╟─0c5e24c0-4aa5-44a2-b2fd-db78795485af
# ╠═b5870a92-89e2-4eae-a79f-6033b1f3489e
# ╟─13cfe122-eea8-4bbb-aaa0-5bcc74a247d1
# ╟─7b0a84b5-60d3-4dd9-89e9-29c88282cb25
# ╟─3a063481-d447-4bf5-9c49-ecde37a0fcea
# ╠═035d4123-7092-4429-8cfd-1e5926e84493
# ╠═26bab6eb-7457-4fb7-b8e2-5148769891ff
# ╠═47c852e5-38eb-4c1f-8d60-d180e4826d05
# ╟─5a0900cc-df10-4176-b903-358b3e00415c
# ╟─076b4a28-be0f-46f0-9857-e6f886c4b118
# ╟─ae6e4bb7-46e8-4f95-b337-0b4589c43cbf
# ╟─bf036015-80a1-4b5a-9ddd-9bfc939979e0
# ╟─4acaadd4-6102-44f5-b602-00465bf3feca
# ╠═9b13bc55-63eb-46bc-acea-f9aec90b340f
# ╟─9274b233-687d-471d-8bac-e14e9a0cb7c0
# ╠═c9c931d9-0a2a-47de-9bee-493c472def48
# ╟─665c5826-4516-4ffa-a9d8-def8e3985ddb
# ╠═a21d9a07-de69-4884-8d7d-742413f9a95a
# ╟─5abe4d49-f398-4d0a-8055-9c4b1014e74f
# ╠═91980f61-4173-4840-a1ae-f1d743daa2c4
# ╠═5623ace0-4b62-4ec7-b54f-c110d38bc06b
# ╟─65dbb492-4795-44ca-afcb-fb2a2c925d92
# ╟─e8425c71-666a-462e-9c4d-fc480810f922
# ╠═056119e3-ec74-4860-abb4-b75f6b16878d
# ╠═77f663a6-96e5-4a1c-843c-cd66fd1382b5
# ╠═1e51701d-a893-4056-8336-a3772b85abe4
# ╠═48366e85-cffb-4c5c-ac3a-807e47f858c7
