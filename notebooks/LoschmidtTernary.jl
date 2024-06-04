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
    \partial_t \rho_i +\nabla\cdot(\rho_i \vec v) + \nabla \cdot \vec J_i  = r_i(\varrho),\qquad i=1,\dots,n.
    \end{align}   
```
"""

# ╔═╡ 43148504-814c-46ec-985a-2d790e1265e4
md"""
To close the model, expressions for the convective velocity $\vec v$, diffusive species mass fluxes $J_i$ and diffusive thermal energy flux $\vec Q$ are introduced corresponding to equations in Section 2.3:

Convective velocity $\vec v$:
```math
    \begin{align}
         \vec v  = -\frac{\kappa}{\nu} \nabla p
    \end{align}
```

Diffusive speceis mass fluxes $J_i$:
```math
\begin{equation}
\begin{split}
	\frac{p}{RT}\frac{1}{M_{\rm mix}}\mathsf{\vec d}_i' &=-\sum_{j:j\not=i}\frac{w_j\vec J_i-w_i\vec J_j}{M_iM_jD_{ij}},\qquad i=1,\dots,n,\\
	\mathsf{\vec d}_i'&=\mathsf{\vec d}_i+x_i\widetilde{\mathcal{X}_i}\nabla\log T\\&=\nabla x_i+(x_i{-}w_i)\nabla \log p+x_i\widetilde{\mathcal{X}_i} \nabla\log T,\qquad i=1,\dots,n. 
\end{split}
\end{equation}
```

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

# ╔═╡ e31b7177-1005-4a33-b876-58cbc9870495
L = 2*40.56ufac"cm"

# ╔═╡ f4286a46-ceb4-40b2-8ce5-7fcd231e3168
md"""
The one-dimensional modeling domain corresponds to Loschmidt diffusion cell of L
= $(L/ufac"cm") cm.

It consists of two chambers of $(L/2/ufac"cm") cm length each and is shown below.
"""

# ╔═╡ 6b3bc830-92b2-4d1a-8e87-ceca93af8c9f
begin
	Ω_bottom = 1
	Ω_top = 2
end;

# ╔═╡ bfebc943-28ef-477b-bc15-6cab9d398d92
function grid_1D(;nref=0,L=L)
	n = 100 * 2^nref
	
 	XLeft = geomspace(0.0, L/2, L*10/n, L/n)
    XRight = geomspace(L/2, L, L/n, L*10/n)
	
    X = glue(XLeft, XRight)
	
	grid=simplexgrid(X)
	cellmask!(grid, [L/2], [L], Ω_top)
	return grid
end

# ╔═╡ 6939978d-9590-407b-80dc-54721c3f672d
md"""
Select problem dimension: $(@bind dim Select([1, 2], default=1))

Select grid refinement level: $(@bind nref Select([0,1,2,3], default=0))

"""

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
CH4ArH2 = KinData{}(;
    ng = 3,
    gnames = [:CH4, :Ar, :H2],
    Fluids = [CH4, Ar, H2],
    gn = Dict(1:3 .=> [:CH4, :Ar, :H2]),
    gni = Dict(value => key for (key, value) in Dict(1:3 .=> [:CH4, :Ar, :H2])),
	nr = 0,
	rnames = [],
)

# ╔═╡ 4c2be37b-f47e-4841-9ad9-a8c222965e06
md"""
In accordance with [1] the following constant values for the binary Maxwell-Stefan diffusivites are applied:
-  $D_{12}$: CH4-Ar = 21.57 mm²/s
-  $D_{13}$: CH4-H2 = 77.16 mm²/s
-  $D_{23}$: Ar-H2 = 83.35 mm²/s
"""

# ╔═╡ 4b80f177-9591-49db-8c51-60e37ce2dab5
constant_binary_diff_coeffs = [21.57, 77.16, 83.35] * ufac"mm^2/s"

# ╔═╡ 7f1d9cf8-7785-48c1-853c-74680188121f
data = ReactorData(
	dim = dim,
	X0 = [1.0,1.0,1.0],
	kinpar = CH4ArH2,
	Tamb = 25.0 + 273.15,
	p = 101.3*ufac"kPa",

	constant_properties = true,
	solve_T_equation = false,
	is_reactive = false,
	include_Soret_Dufour = false,
	
	γ_τ = 1.0,
	poros=1.0,

	perm = 1.23e-10*ufac"m^2" * 1.0e6,
	
	constant_binary_diff_coeffs = constant_binary_diff_coeffs,
)

# ╔═╡ 246c0a65-9504-44cc-91ee-d5d6c186d4fe
md"""
Define simulation parameters e.g. the initial compositions in both chambers of the Loschmidt diffusion cell, temperature and pressure.
"""

# ╔═╡ 9a341a65-b560-48a6-b954-96ea6740c44f
begin
	X0_bottom = [
		0.0,
		0.509,
		0.491
	]
	X0_top = [
		0.515,
		0.485,
		0.0
	]
end;

# ╔═╡ bcfc3138-3cbb-4386-a33b-573b6c39caf9
md"""
## Boundary conditions
"""

# ╔═╡ 0c5e24c0-4aa5-44a2-b2fd-db78795485af
md"""
The modeling domain corresponds to a closed system, so no-flux boundary conditions for species mass apply on all domain boundaries. The homogeneous Neumann boundary conditions correspond to the default boundary conditon and thus nothing needs to be specified in particular.
"""

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

# ╔═╡ c40265d6-a9e0-443f-bfec-ca70418d8361
md"""
$(LocalResource("../data/Goell2012/Loschmidt_Krishna1993.png"))
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

# ╔═╡ 379dc67a-41c8-479d-b86f-afa3b77594b1
Xpm(i, X0m, X0p, grid) = map((x) -> x == L/2 ? (X0m[i]+X0p[i])/2 : (x < L/2 ? X0m[i] : X0p[i]) , grid)

# ╔═╡ 1e51701d-a893-4056-8336-a3772b85abe4
function setup_run_sim(grid, data)
	(;ng, ip, iT, Tamb, p, X0, outlet_boundaries) = data
	sys=VoronoiFVM.System(
		grid;
		data=data,
		flux=DMS_flux,
		reaction=DMS_reaction,
		storage=DMS_storage,
		assembly=:edgewise,
		unknown_storage=:dense
	)
	
	enable_species!(sys; species=collect(1:(ng+1))) # gas phase species xi, ptotal, T
		
	inival=unknowns(sys)
	inival[ip,:].=p


	for i=1:ng
		#inival[i,:] .= X0[i]
		inival[i,:] .= Xpm(i, X0_bottom, X0_top, grid)
	end

	control = SolverControl(nothing, sys;)
	control.handle_exceptions=true
	control.Δt_max=100.0
	#control.Δt_min = 1.0e-4
	#control.Δt = 1.0e-3
	

	times=[0,2ufac"hr"]
	
	sol=VoronoiFVM.solve(sys;inival=inival,times,control,log=true)
	return (sol,sys)
end

# ╔═╡ ff501186-d8a0-4666-8991-fe576f8ff6ad
function grid_2D(;nref=0)
	X = collect(0:(1/2^nref):18)*ufac"cm"
    Y = collect(0:(0.5/2^nref):8)*ufac"cm"
    grid = simplexgrid(X, Y)
    
	rect!(grid, [4, 0]*ufac"cm", [14, 2.5]*ufac"cm"; region = 2)
	rect!(grid, [4, 5.5]*ufac"cm", [14, 8]*ufac"cm"; region = 2, bregions= [3, 3, 3, 3])
	bfacemask!(grid, [0, 0]*ufac"cm", [4, 0]*ufac"cm", 4)
	bfacemask!(grid, [0, 8]*ufac"cm", [4, 8]*ufac"cm", 4)
	
	subgrid(grid, [1])
end

# ╔═╡ 138b4d12-af0f-4a1c-a4aa-4f55e6038664
if dim == 1
	grid = grid_1D(nref=nref)
	const Γ_left = 1
	const Γ_right = 2
elseif dim == 2
	grid = grid_2D(nref=nref)
	const Γ_left = 4
	const Γ_right = 2
end;

# ╔═╡ b55537bf-9982-4997-8a2a-1972127bdd86
let
	vis = GridVisualizer(resolution=(450,225))
	gridplot!(vis, grid, zoom=2, linewidth=1)
	reveal(vis)
end

# ╔═╡ 035d4123-7092-4429-8cfd-1e5926e84493
solt, sys = setup_run_sim(grid, data);

# ╔═╡ 5a0900cc-df10-4176-b903-358b3e00415c
md"""
Move the slider to change $t$: $(if isa(solt, TransientSolution)
	@bind t PlutoUI.Slider(solt.t,show_value=true,default=solt.t[end])	
end)
"""

# ╔═╡ 26bab6eb-7457-4fb7-b8e2-5148769891ff
sol = solt(t);

# ╔═╡ ae6e4bb7-46e8-4f95-b337-0b4589c43cbf
let
	(;ip,gni) = data
	vis =GridVisualizer(layout=(1,2), resolution=(600,300))
	scalarplot!(vis[1,1], grid, sol[gni[:Ar],:], title = "Ar molar fraction",)
	scalarplot!(vis[1,2], grid, sol[gni[:CH4],:], title = "CH4 molar fraction",)

	reveal(vis)
end

# ╔═╡ 076b4a28-be0f-46f0-9857-e6f886c4b118
md"""
Plot the transient solution at time $t=$ $(round(t)):
"""

# ╔═╡ 0aef8cc3-daea-4dd0-98bc-4188b1baffc9
function fcn_identity(f,u,node,data)
	(;ip) = data
	ng=ngas(data)

    for i=1:ng
	    f[i] = u[i]
    end
    f[ip] = u[ip]
end

# ╔═╡ cf0f4137-17a4-4859-8058-ec76b1f28290
function xAvgBottomTop(sys, solt, data; L=L)
	ng = ngas(data)
	l = length(solt.t)

	x_avg_bottom = Array{Float64, 2}(undef, ng,l)
	x_avg_top = Array{Float64, 2}(undef, ng,l)
	

	for (i,t) in enumerate(solt.t)
		x_avg = integrate(sys, fcn_identity, solt(t); boundary=false)
		x_avg ./= L/2
		for j=1:ng
			x_avg_bottom[j,i] = x_avg[j,Ω_bottom]
			x_avg_top[j,i] = x_avg[j,Ω_top]
		end
	end
	
	x_avg_bottom, x_avg_top
end

# ╔═╡ 9ba456c7-f6ed-49c5-9878-d9e0deb96384
let
	(; gni,iT,ip) = data;
		
	f = Figure(size = (600, 300))
	
	xs = solt.t/ufac"hr"
	Bottom = 1
	Top = length(grid[Coordinates])

	x_avg_bottom, x_avg_top = xAvgBottomTop(sys, solt, data)

	ax1 = Axis(f[1, 1], xlabel="time / hr", ylabel="Mole fraction", title="Argon", limits = (nothing, nothing, 0.39, 0.61))

	# Ar
	lines!(xs, x_avg_bottom[gni[:Ar],:], label = "Bottom")
	scatter!(xAr_Bottom.time, xAr_Bottom.xAr_Bottom)
	
	lines!(xs, x_avg_top[gni[:Ar],:], label = "Top")
	scatter!(xAr_Top.time, xAr_Top.xAr_Top)
	axislegend(ax1, merge=true)
	
	
	# CH4
	ax2 = Axis(f[1, 2], xlabel="time / hr", ylabel="Mole fraction", title="Methane", limits = (nothing, nothing, nothing, 0.61))
	
	lines!(xs, x_avg_bottom[gni[:CH4],:], label = "Bottom")
	scatter!(xCH4_Bottom.time, xCH4_Bottom.xCH4_Bottom)
	
	lines!(xs, x_avg_top[gni[:CH4],:], label = "Top")
	scatter!(xCH4_Top.time, xCH4_Top.xCH4_Top)
	axislegend(ax2, merge=true)
	f	
end

# ╔═╡ 3c4ae416-a47b-451d-b870-fa10166d97de
function plotting_movie(; filename = "plotting_video.gif", Plotter = default_plotter())
    vis =GridVisualizer(layout=(2,2), resolution=(1000,600))

	(; gni,iT,ip) = data;
    movie(vis; file = filename) do vis
        for t in solt.t
            scalarplot!(vis[1, 1],
                        grid,
                        solt(t)[gni[:He],:];
                        clear = true,
                        title = "He molar fraction, t=$(Integer(round(t))) s",
                        limits = (0.32,0.345),
                        levels = 7,
                        colormap = :coolwarm)
            scalarplot!(vis[1, 2],
                        grid,
                        solt(t)[gni[:Kr],:];
                        clear = true,
                        title = "Kr molar fraction, t=$(Integer(round(t))) s",
                        limits = (0.325,0.34),
                        levels = 7,
                        colormap = :coolwarm)
			scalarplot!(vis[2, 1],
                        grid,
                        solt(t)[iT,:];
                        clear = true,
                        title = "Temperature (K), t=$(Integer(round(t))) s",
                        limits = (300.0,400.0),
                        levels = 7,
                        colormap = :gist_heat)

			t0, tend = solt.t[1], solt.t[end]
			pmax = round(maximum(solt(t)[ip,:]))
			f = map((x, y) -> pmax, grid)
			scalarplot!(vis[2, 2],
                        grid,
                        f;
                        clear = true,
                        title = "Pressure (Pa), t=$(Integer(round(t))) s",
                        limits = (solt(t0)[ip,1],solt(tend)[ip,1]),
                        levels = 7,
                        colormap = :viridis)
            reveal(vis)
        end
    end
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
# ╠═e31b7177-1005-4a33-b876-58cbc9870495
# ╠═26bab6eb-7457-4fb7-b8e2-5148769891ff
# ╟─f4286a46-ceb4-40b2-8ce5-7fcd231e3168
# ╠═6b3bc830-92b2-4d1a-8e87-ceca93af8c9f
# ╠═bfebc943-28ef-477b-bc15-6cab9d398d92
# ╠═b55537bf-9982-4997-8a2a-1972127bdd86
# ╟─6939978d-9590-407b-80dc-54721c3f672d
# ╠═7f1d9cf8-7785-48c1-853c-74680188121f
# ╠═ae6e4bb7-46e8-4f95-b337-0b4589c43cbf
# ╟─5a0900cc-df10-4176-b903-358b3e00415c
# ╠═138b4d12-af0f-4a1c-a4aa-4f55e6038664
# ╟─e9cb07eb-cfbb-4802-bc7f-6de7a6ad8ac6
# ╟─f008f30f-0137-4c54-8d4e-a6f589c4a952
# ╠═e7ca4902-0e14-48ca-bcc6-96b06c85a39d
# ╟─4c2be37b-f47e-4841-9ad9-a8c222965e06
# ╠═4b80f177-9591-49db-8c51-60e37ce2dab5
# ╟─246c0a65-9504-44cc-91ee-d5d6c186d4fe
# ╠═9a341a65-b560-48a6-b954-96ea6740c44f
# ╟─bcfc3138-3cbb-4386-a33b-573b6c39caf9
# ╟─0c5e24c0-4aa5-44a2-b2fd-db78795485af
# ╟─13cfe122-eea8-4bbb-aaa0-5bcc74a247d1
# ╟─a4fd9977-577e-45c9-a5c3-c0cd8a5dd012
# ╟─7b0a84b5-60d3-4dd9-89e9-29c88282cb25
# ╟─3a063481-d447-4bf5-9c49-ecde37a0fcea
# ╠═035d4123-7092-4429-8cfd-1e5926e84493
# ╟─076b4a28-be0f-46f0-9857-e6f886c4b118
# ╟─9ba456c7-f6ed-49c5-9878-d9e0deb96384
# ╟─c40265d6-a9e0-443f-bfec-ca70418d8361
# ╠═af0a2719-3e0a-420e-9c8c-4cbbcb828cb1
# ╟─65dbb492-4795-44ca-afcb-fb2a2c925d92
# ╟─e8425c71-666a-462e-9c4d-fc480810f922
# ╠═379dc67a-41c8-479d-b86f-afa3b77594b1
# ╠═1e51701d-a893-4056-8336-a3772b85abe4
# ╠═ff501186-d8a0-4666-8991-fe576f8ff6ad
# ╠═0aef8cc3-daea-4dd0-98bc-4188b1baffc9
# ╠═cf0f4137-17a4-4859-8058-ec76b1f28290
# ╟─3c4ae416-a47b-451d-b870-fa10166d97de
