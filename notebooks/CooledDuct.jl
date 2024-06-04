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
	using VoronoiFVM
	using ExtendableGrids, GridVisualize, CairoMakie
	using CSV, DataFrames
	using Revise
	using MultiComponentReactiveMixtureProject
	GridVisualize.default_plotter!(CairoMakie)	
end;

# ╔═╡ 0f102f06-3ff3-4bcc-8892-8d9190a87849
TableOfContents(aside=true)

# ╔═╡ d14462c6-f63b-4a61-a1d9-4bcdb8e30e3d
md"""
Supplementary notebook for:$br
__Transport of heat and mass for reactive gas mixtures in porous media: modeling and application__.

This notebook is used as validation of the iso-thermal model implementation by reproducing experimental results in [1] and simulation results in [2].
"""

# ╔═╡ 38c3ddb4-b44a-4981-9000-0a1d303bd9ac
md"""
# Introduction
The isothermal modeling framework is applied to a ternary gas mixture diffusing in a Loschmidt tube. Experimental results are presented in [1] and simulation results using an alternative approach are presented in [2].
"""

# ╔═╡ 2628cb2d-c1ef-4ad0-8ee4-38e45f864838
md"""
## Model equations
"""

# ╔═╡ 75fddee6-e057-4e91-a239-2033370b00fc
md"""
Reiterating the presented modeling equations from Section 2.3, the species mass balances are solved for the isothermal conditions for the ternary gas mixture in the Loschmidt diffusion cell:
```math
\begin{align}
    \partial_t \rho_i + \nabla \cdot \vec J_i  = 0,\qquad i=1,\dots,n.
    \end{align}   
```
"""

# ╔═╡ 43148504-814c-46ec-985a-2d790e1265e4
md"""
To complete close the model, expressions for the diffusive species mass fluxes $J_i$ are introduced corresponding to equations in Section 2.3:

Diffusive speceis mass fluxes $J_i$:
```math
\begin{equation}
\begin{split}
	\frac{p}{RT}\frac{1}{M_{\rm mix}}\mathsf{\vec d}_i &=-\sum_{j:j\not=i}\frac{w_j\vec J_i-w_i\vec J_j}{M_iM_jD_{ij}},\qquad i=1,\dots,n,\\
	\mathsf{\vec d}_i&=\mathsf{\vec d}_i
	\\&=\nabla x_i+(x_i{-}w_i)\nabla \log p ,\qquad i=1,\dots,n. 
\end{split}
\end{equation}
```
"""

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

# ╔═╡ 138b4d12-af0f-4a1c-a4aa-4f55e6038664
begin
	const Γ_bottom = 1
	const Γ_right = 2
	const Γ_top = 3
	const Γ_left = 4
	const Γ_right_outflow = 5
	
	const Ω_permeable = 2
	# const Ω_free = 1
end;

# ╔═╡ 6939978d-9590-407b-80dc-54721c3f672d
md"""
Select grid refinement level: $(@bind nref Select([0,1,2,3], default=0))

"""

# ╔═╡ 7be1c79e-08d8-493c-bce0-114c0c003dd7
function grid_2D(;nref=0, L=80ufac"mm", H=10ufac"mm", efix=0.01ufac"mm")

	nx = 20*2^nref
	ny = 10*2^nref

	lmin = L/nx/10
    lmax = L/nx

	hmin = H/ny/10
	hmax = H/ny
	
	XLeft = collect(0:L/nx:(3/4*L)) # uniform spacing
    XRight = geomspace(3/4*L, L, lmax, lmin)	
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
	cellmask!(grid, [0,0], [L,H], Ω_permeable) # gas permeable region
	bfacemask!(grid, [L,0], [L,H], Γ_right_outflow)

	grid
end

# ╔═╡ 869652e5-f15e-43d4-8fbc-724e866892b6
grid = grid_2D(nref=nref)

# ╔═╡ b55537bf-9982-4997-8a2a-1972127bdd86
let
	vis = GridVisualizer(resolution=(600,150))
	gridplot!(vis, grid, linewidth=0.5, aspect=2.0)
	reveal(vis)
end

# ╔═╡ e9cb07eb-cfbb-4802-bc7f-6de7a6ad8ac6
md"""
## Simulation parameters
"""

# ╔═╡ f008f30f-0137-4c54-8d4e-a6f589c4a952
md"""
Setup data structure corresponding to the ernary gas mixture consisting of the noble gases:
1) Methane (CH4)
2) Argon (Ar)
3) Hydrogen (H2)
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
	inflow_boundaries = [Γ_left],
	outflow_boundaries = [Γ_right_outflow],
	permeable_regions = [Ω_permeable],
	X0 = [0.21,0.79],
	mflowin = 2.0e-3,
	kinpar = Air,
	Tamb = 50.0 + 273.15,
	#T_gas_in = 50.0 + 273.15,
	T_gas_in = 100 + 273.15,
	
	p = 101.3*ufac"kPa",

	solve_T_equation = true,
	is_reactive = false,
	include_Soret_Dufour = false,
	
	perm = [0.0,1.0]*1.23e-10*ufac"m^2",
	rhos = 5900.0*ufac"kg/m^3",
	cs = 500.0*ufac"J/(kg*K)",
	lambdas = 1.5*ufac"W/(m*K)",
	poros = [0.0,0.9],
	γ_τ = [1.0,1.1]

	# perm = 1.23e-10*ufac"m^2" * 1.0e6,
)

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
	
	boundary_robin!(f,u,bnode, species=iT,region=Γ_bottom,value=Tamb * eps_,factor=eps_)
	boundary_robin!(f,u,bnode, species=iT,region=Γ_top,value=Tamb * eps_,factor=eps_)

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
			# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative
			f[iT] = -r_hf_enth * ramp(bnode.time; du=(0.0,1), dt=dt_hf_enth)
						
		end

	end

	for boundary in outflow_boundaries
        boundary_dirichlet!(f,u,bnode, species=ip,region=boundary,value=p)
    end	
		
end

# ╔═╡ 13cfe122-eea8-4bbb-aaa0-5bcc74a247d1
md"""
# Solving and plotting
"""

# ╔═╡ a4fd9977-577e-45c9-a5c3-c0cd8a5dd012
md"""
Following the imposed boundray conditions for the thermal energy equation, a temperature gradient of 100 K establishes over the length of the domain that is the driving force for thermodiffusion, leading to a partial separation of the initially uniform, equimolar mixture.
"""

# ╔═╡ 7b0a84b5-60d3-4dd9-89e9-29c88282cb25
md"""
## Transient solution
"""

# ╔═╡ 3a063481-d447-4bf5-9c49-ecde37a0fcea
md"""
Setup the system of equations as a transient system, solve the system and return the transient solution:
"""

# ╔═╡ 9274b233-687d-471d-8bac-e14e9a0cb7c0
md"""
## Comparison with published literature
The time evolution of the species mole fraction averaged over the Bottom and Top parts of the Loschmidt diffusion cell as shown in [1], Figure 5.7 is reproduced below using the model and implementation presented in this work.
"""

# ╔═╡ af0a2719-3e0a-420e-9c8c-4cbbcb828cb1
begin
	path = "../data/Goell2012/Loschmidt/"
	xAr_Bottom = DataFrame(
		CSV.File(
			path*"xAr_Bottom.csv",
			delim=';',
			header=["time", "xAr_Bottom"]),		
	)
	xAr_Top = DataFrame(
		CSV.File(
			path*"xAr_Top.csv",
			delim=';',
			header=["time", "xAr_Top"]),		
	)
	xCH4_Bottom = DataFrame(
		CSV.File(
			path*"xCH4_Bottom.csv",
			delim=';',
			header=["time", "xCH4_Bottom"]),		
	)
	xCH4_Top = DataFrame(
		CSV.File(
			path*"xCH4_Top.csv",
			delim=';',
			header=["time", "xCH4_Top"]),		
	)
end;

# ╔═╡ 65dbb492-4795-44ca-afcb-fb2a2c925d92
md"""
# References
1) __Taylor, Ross; Krishna, Rajamani (1993)__: Multicomponent mass transfer. Wiley, New York.
1) __Göll, Stephan; Piesche, Manfred (2012)__: Multi-component gas transport in micro-porous domains: Multidimensional simulation at the macroscale. In: International Journal of Heat and Mass Transfer 55 (1-3), S. 480–487. DOI: 10.1016/j.ijheatmasstransfer.2011.09.049.
"""

# ╔═╡ e8425c71-666a-462e-9c4d-fc480810f922
md"""
# Function definitions
"""

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

	

	control = SolverControl(nothing, sys;)
	control.handle_exceptions=true
	control.Δt_max=100.0
	control.Δu_opt=200.0
	

	times=[0,3000.0]
	
	sol=VoronoiFVM.solve(sys;inival=inival,times,control,log=true)
	return (sol,sys)
end

# ╔═╡ 035d4123-7092-4429-8cfd-1e5926e84493
solt, sys = setup_run_sim(grid, data);

# ╔═╡ 5a0900cc-df10-4176-b903-358b3e00415c
md"""
Move the slider to change $t$: $(if isa(solt, TransientSolution)
	@bind t PlutoUI.Slider(solt.t,show_value=false,default=solt.t[end])	
end)
"""

# ╔═╡ 076b4a28-be0f-46f0-9857-e6f886c4b118
md"""
Plot the transient solution at time $t=$ $(round(t)):
"""

# ╔═╡ 26bab6eb-7457-4fb7-b8e2-5148769891ff
sol = solt(t);

# ╔═╡ ae6e4bb7-46e8-4f95-b337-0b4589c43cbf
let
	(;iT,ip,gni) = data
	vis =GridVisualizer(layout=(4,1), resolution=(600,600))
	scalarplot!(vis[1,1], grid, sol[gni[:O2],:], title = "O2 molar fraction",)
	scalarplot!(vis[2,1], grid, sol[gni[:N2],:], title = "N2 molar fraction",)
	scalarplot!(vis[3,1], grid, sol[iT,:], title = "Temperature",)
	scalarplot!(vis[4,1], grid, sol[ip,:], title = "Pressure",)

	reveal(vis)
end

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
sol_ss = run_ss(solt,sys);

# ╔═╡ 4acaadd4-6102-44f5-b602-00465bf3feca
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
	
	nf = nodeflux(sys, sol_ss)

	massflux = nf[:,ip,:]
	vel_x = massflux[1,:]./dens
	vel_x_ = vel_x[dens .== 0.0] .= 0.0
	vel_y = massflux[2,:]./dens
	
	vis =GridVisualizer(layout=(5,1), resolution=(600,750))
	
	scalarplot!(vis[1,1], grid, T, title = "Temperature / K")
	scalarplot!(vis[2,1], grid, p, limits=(minimum(p_),maximum(p_)), title = "Pressure / Pa",)
	scalarplot!(vis[3,1], grid, dens, limits=(minimum(dens_),maximum(dens_)), title = "Density / kg m-3",)
		
	scalarplot!(vis[4,1], grid, massflux[1,:], title = "Mass flux (X) / kg s-1 m-2")
	#scalarplot!(vis[4,1], grid, massflux[2,:], title = "Mass flux (Y) / kg s-1 m-2")
	
	scalarplot!(vis[5,1], grid, vel_x, limits=(minimum(vel_x_),maximum(vel_x_)), title = "Velocity (X) / m s-1")
	#scalarplot!(vis[6,1], grid, vel_y, title = "Velocity (Y) / m s-1")

	reveal(vis)
	#massflux[1,:]
	
end

# ╔═╡ 7e09011f-9bc7-4925-953c-803b6ac1869f
Print_summary(sol_ss,grid,sys,data)

# ╔═╡ 0aef8cc3-daea-4dd0-98bc-4188b1baffc9
function fcn_identity(f,u,node,data)
	(;ip) = data
	ng=ngas(data)

    for i=1:ng
	    f[i] = u[i]
    end
    f[ip] = u[ip]
end

# ╔═╡ Cell order:
# ╠═349e7220-dc69-11ee-13d2-8f95e6ee5c96
# ╠═0f102f06-3ff3-4bcc-8892-8d9190a87849
# ╟─d14462c6-f63b-4a61-a1d9-4bcdb8e30e3d
# ╟─38c3ddb4-b44a-4981-9000-0a1d303bd9ac
# ╟─2628cb2d-c1ef-4ad0-8ee4-38e45f864838
# ╟─75fddee6-e057-4e91-a239-2033370b00fc
# ╟─43148504-814c-46ec-985a-2d790e1265e4
# ╟─b3a6fe03-be46-4159-96ab-477a42d0eec5
# ╟─c3dbf8b3-fc5f-44ff-be2c-ca4200f5bd6c
# ╠═138b4d12-af0f-4a1c-a4aa-4f55e6038664
# ╠═869652e5-f15e-43d4-8fbc-724e866892b6
# ╠═b55537bf-9982-4997-8a2a-1972127bdd86
# ╟─6939978d-9590-407b-80dc-54721c3f672d
# ╠═7be1c79e-08d8-493c-bce0-114c0c003dd7
# ╠═7f1d9cf8-7785-48c1-853c-74680188121f
# ╟─e9cb07eb-cfbb-4802-bc7f-6de7a6ad8ac6
# ╟─f008f30f-0137-4c54-8d4e-a6f589c4a952
# ╠═e7ca4902-0e14-48ca-bcc6-96b06c85a39d
# ╟─bcfc3138-3cbb-4386-a33b-573b6c39caf9
# ╟─0c5e24c0-4aa5-44a2-b2fd-db78795485af
# ╠═b5870a92-89e2-4eae-a79f-6033b1f3489e
# ╟─13cfe122-eea8-4bbb-aaa0-5bcc74a247d1
# ╟─a4fd9977-577e-45c9-a5c3-c0cd8a5dd012
# ╟─7b0a84b5-60d3-4dd9-89e9-29c88282cb25
# ╟─3a063481-d447-4bf5-9c49-ecde37a0fcea
# ╠═035d4123-7092-4429-8cfd-1e5926e84493
# ╟─5a0900cc-df10-4176-b903-358b3e00415c
# ╟─076b4a28-be0f-46f0-9857-e6f886c4b118
# ╠═ae6e4bb7-46e8-4f95-b337-0b4589c43cbf
# ╠═26bab6eb-7457-4fb7-b8e2-5148769891ff
# ╟─9274b233-687d-471d-8bac-e14e9a0cb7c0
# ╠═9b13bc55-63eb-46bc-acea-f9aec90b340f
# ╠═4acaadd4-6102-44f5-b602-00465bf3feca
# ╟─af0a2719-3e0a-420e-9c8c-4cbbcb828cb1
# ╟─65dbb492-4795-44ca-afcb-fb2a2c925d92
# ╟─e8425c71-666a-462e-9c4d-fc480810f922
# ╠═7e09011f-9bc7-4925-953c-803b6ac1869f
# ╠═1e51701d-a893-4056-8336-a3772b85abe4
# ╠═48366e85-cffb-4c5c-ac3a-807e47f858c7
# ╠═0aef8cc3-daea-4dd0-98bc-4188b1baffc9
