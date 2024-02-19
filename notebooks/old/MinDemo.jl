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
	using NLsolve, LinearSolve, Pardiso, ILUZero
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots, Printf
	using PlutoUI, Colors

	using MultiComponentReactiveMixtureProject
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
PlutoUI.TableOfContents(title="Demo")

# ╔═╡ 2c4b16ab-8f60-467b-b608-2fea9fbc741c
md"""
# Introduction
Demo Notebook with 3 ideal gas species and constant properties to demonstrate isothermal Multicomponent species transport.

Select problem dimension: $(@bind dim Select([1, 2, 3], default=2))

Check the box to __start the simulation__: $(@bind RunSim PlutoUI.CheckBox(default=false))
"""

# ╔═╡ 4e05ab31-7729-4a4b-9c14-145118477715
if dim == 3
	@bind xcut Slider(linspace(0,1,21)*ufac"cm",show_value=true,default=0.5*ufac"cm")
end

# ╔═╡ 832f3c15-b75a-4afe-8cc5-75ff3b4704d6
begin
	const Ω_catalyst = 2
	if dim == 1
		const Γ_left = 1
		const Γ_right = 2
	elseif dim == 2
		const Γ_bottom = 1
		const Γ_right = 2
		const Γ_top = 3
		const Γ_left = 4
	else
		const Γ_left = 1
		const Γ_front = 2
		const Γ_right = 3
		const Γ_back = 4
		const Γ_bottom = 5
		const Γ_top = 6
	end
end;

# ╔═╡ 83fa22fa-451d-4c30-a4b7-834974245996
function grid1D()
	X=(0:0.02:1)*ufac"cm"
	grid=simplexgrid(X)
	# catalyst region
	cellmask!(grid,[0.4]*ufac"cm",[0.6]*ufac"cm",Ω_catalyst)	
	grid
end

# ╔═╡ 4dae4173-0363-40bc-a9ca-ce5b4d5224cd
function grid2D()
	X=(0:0.1:1)*ufac"cm"
	Y=(0:0.1:1)*ufac"cm"
	grid=simplexgrid(X,Y)

	# catalyst region
	cellmask!(grid,[0.4,0.3].*ufac"cm",[0.6,0.7].*ufac"cm",Ω_catalyst)
	
	grid
end

# ╔═╡ 561e96e2-2d48-4eb6-bb9d-ae167a622aeb
function grid3D()
	X=(0:0.1:1)*ufac"cm"
	Y=(0:0.1:1)*ufac"cm"
	Z=(0:0.1:1)*ufac"cm"
	grid=simplexgrid(X,Y,Z)

	# catalyst region
	cellmask!(grid,[0.3,0.4,0.3].*ufac"cm",[0.7,0.6,0.7].*ufac"cm",Ω_catalyst)
	#bfacemask!(grid, [0.2,1,0.2].*ufac"cm",[0.8,1,0.8].*ufac"cm",7)
	
	grid
end

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
let
	if dim == 1
		gridplot(grid1D(), resolution=(600,200))
	elseif dim == 2
		gridplot(grid2D())
	else
		gridplot(grid3D(); xplane=xcut, show=true, outlinealpha=0.0 )
	end
end

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
Reaction leading to increase in moles.

$B \rightarrow 3 A$
"""

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
#begin
function MinDemo(dim; times=nothing, mfluxin = nothing)
	
	if dim == 1
		mygrid=grid1D()
		inb = [Γ_left]
		irrb = []
		outb = [Γ_right]
		sb = []
		times = isnothing(times) ? [0,20.0] : times
	elseif dim == 2
		mygrid=grid2D()
		inb = [Γ_left]
		irrb = []
		outb = [Γ_right]
		sb = []
		times = isnothing(times) ? [0,12.0] : times
	else
		mygrid=grid3D()
		inb = [Γ_left]
		irrb = []
		outb = [Γ_right]
		sb = []
		times = isnothing(times) ? [0,3.0] : times
	end

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
	
	mydata = ReactorData(;
		inlet_boundaries=inb,
		irradiated_boundaries=irrb,
		outlet_boundaries=outb,
		side_boundaries=sb,
		kinpar=MinKin,
		lcat = 1.0,
		Tamb = 273.15,
		m = [2.0,6.0,21.0]*ufac"g/mol",
		X0 = [0.2, 0.8, 0.0],
		mfluxin = isnothing(mfluxin) ? 0.01*ufac"kg/(m^2*s)" : mfluxin,
		solve_T_equation = false,
		constant_properties = true,
		is_reactive = true
	)
	
	(;p,ip,X0)=mydata
	ng=ngas(mydata)
	
	sys=VoronoiFVM.System( 	mygrid;
							data=mydata,
							flux=MultiComponentReactiveMixtureProjectponentReactiveMixtureProject.DMS_flux,
							reaction=MultiComponentReactiveMixtureProject.DMS_reaction,
							storage=MultiComponentReactiveMixtureProject.DMS_storage,
							bcondition=MultiComponentReactiveMixtureProject.PTR_bcond,
							boutflow=MultiComponentReactiveMixtureProject.DMS_boutflow,
							outflowboundaries=outb,
							assembly=:edgewise
							)
	
	enable_species!(sys; species=collect(1:(ng+1))) # gas phase species pi & ptotal	
	
	inival=unknowns(sys)

	inival[ip,:].=p
	for i=1:ng
		inival[i,:] .= X0[i]
	end

	if dim == 1 || dim == 2
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
		control.Δu_opt=1.0e5

	solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="nae")

	return solt,mygrid,sys,mydata
end

# ╔═╡ b9f80e49-bf9f-4a92-afba-867ab1cb1d1f
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
if RunSim
	solt,grid,sys,data=MinDemo(dim);
end;
  ╠═╡ =#

# ╔═╡ 05949759-2bb9-475b-b2f4-900b32c30e00
md"""
## Molar Fractions
"""

# ╔═╡ f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╠═╡ skip_as_script = true
#=╠═╡
@bind t Slider(solt.t,show_value=true,default=solt.t[end])
  ╠═╡ =#

# ╔═╡ 5588790a-73d4-435d-950f-515ae2de923c
# ╠═╡ skip_as_script = true
#=╠═╡
sol = solt(t);
  ╠═╡ =#

# ╔═╡ e29848dd-d787-438e-9c32-e9c2136aec4f
# ╠═╡ skip_as_script = true
#=╠═╡
MultiComponentReactiveMixtureProject._checkinout(sol,sys,data)
  ╠═╡ =#

# ╔═╡ 862bf54f-8700-4956-9024-07fdf809c922
#=╠═╡
MultiComponentReactiveMixtureProject.Print_summary(sol,grid,sys,data)
  ╠═╡ =#

# ╔═╡ 06ec3b97-5532-4a86-9abb-61b91e94b4e7
#=╠═╡
if dim==3
	MultiComponentReactiveMixtureProject.WriteSolution3D(sol,grid,data;desc="MinDemo")
end
  ╠═╡ =#

# ╔═╡ 111b1b1f-51a5-4069-a365-a713c92b79f4
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;ip,p) = data
	if dim == 1
		vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300))
		scalarplot!(vis, grid, sol[1,:], clear=false, label="x1")
		scalarplot!(vis, grid, sol[2,:], clear=false, color=:red, label="x2")
		scalarplot!(vis, grid, sol[3,:], clear=false, color=:blue, label="x3")
	elseif dim == 2
		vis=GridVisualizer(layout=(1,3), resolution=(700,300))
		scalarplot!(vis[1,1], grid, sol[1,:])
		scalarplot!(vis[1,2], grid, sol[2,:])
		scalarplot!(vis[1,3], grid, sol[3,:])
	else
		vis=GridVisualizer(layout=(3,1), resolution=(400,1200), outlinealpha=0.0)
		scalarplot!(vis[1,1], grid, sol[1,:])
		scalarplot!(vis[2,1], grid, sol[2,:])
		scalarplot!(vis[3,1], grid, sol[3,:])
	end
	
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ d4aaf146-e551-444a-a2f1-e70b05104b53
md"""
## Pressure, Density, Flow
1) Total Pressure
2) Density
3) Velocity - X
4) Velocity - Y
"""

# ╔═╡ de69f808-2618-4add-b092-522a1d7e0bb7
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;p,m,ip,Tamb, mfluxin) = data
	ng = ngas(data)
	mmix = []
	for j in 1:length(sol[1,:])
		_mmix=0
		for i=1:ng
			_mmix += sol[i,j]*m[i]
		end
		push!(mmix, _mmix)
	end
	
	w1 = sol[1,:]*m[1] ./mmix
	w2 = sol[2,:]*m[2] ./mmix
	w3 = sol[3,:]*m[3] ./mmix

	ps = sol[ip,:]
	rho = @. ps * mmix /(ph"R"*Tamb)
	
	if dim == 1
		vis=GridVisualizer(legend=:lt, resolution=(600,400), layout = (2, 1))
		scalarplot!(vis[1,1], grid, w1, clear=false, label="w1")
		scalarplot!(vis[1,1], grid, w2, clear=false, color=:red, label="w2")
		scalarplot!(vis[1,1], grid, w3, clear=false, color=:blue, label="w3")

		p0 = sol[ip,1]
		rho0 = @. p0 * mmix[1] /(ph"R"*Tamb)
		scalarplot!(vis[2,1], grid, rho/rho0, clear=false, label="Rho / Rho0")
		scalarplot!(vis[2,1], grid, rho0./rho, clear=false, color=:red, label="v / v0")
		scalarplot!(vis[2,1], grid, ps/p0, clear=false, color=:blue, label="p / p0")
	elseif dim == 2
		vis=GridVisualizer(layout=(2,2), resolution=(660,660))
		scalarplot!(vis[1,1], grid, ps)
		scalarplot!(vis[1,2], grid, rho)
		nf = nodeflux(sys, sol)
		massflux = nf[:,ip,:]
		scalarplot!(vis[2,1], grid, massflux[1,:]./rho,)
		scalarplot!(vis[2,2], grid, massflux[2,:]./rho,)
	else
		vis=GridVisualizer(layout=(1,2), resolution=(660,660), outlinealpha=0.0)
		scalarplot!(vis[1,1], grid, ps, )
		scalarplot!(vis[1,2], grid, rho, )
		
	end
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═2c4b16ab-8f60-467b-b608-2fea9fbc741c
# ╟─a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╟─4e05ab31-7729-4a4b-9c14-145118477715
# ╠═83fa22fa-451d-4c30-a4b7-834974245996
# ╠═4dae4173-0363-40bc-a9ca-ce5b4d5224cd
# ╠═561e96e2-2d48-4eb6-bb9d-ae167a622aeb
# ╠═832f3c15-b75a-4afe-8cc5-75ff3b4704d6
# ╟─0fadb9d2-1ccf-4d44-b748-b76d911784ca
# ╟─b94513c2-c94e-4bcb-9342-47ea48fbfd14
# ╟─c886dd12-a90c-40ab-b9d0-32934c17baee
# ╟─3440d4d8-3e03-4ff3-93f1-9afd7aaf9c41
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═b9f80e49-bf9f-4a92-afba-867ab1cb1d1f
# ╠═5588790a-73d4-435d-950f-515ae2de923c
# ╠═e29848dd-d787-438e-9c32-e9c2136aec4f
# ╠═862bf54f-8700-4956-9024-07fdf809c922
# ╠═06ec3b97-5532-4a86-9abb-61b91e94b4e7
# ╟─05949759-2bb9-475b-b2f4-900b32c30e00
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╟─111b1b1f-51a5-4069-a365-a713c92b79f4
# ╟─d4aaf146-e551-444a-a2f1-e70b05104b53
# ╟─de69f808-2618-4add-b092-522a1d7e0bb7
