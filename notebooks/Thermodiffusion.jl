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
	using ExtendableGrids, GridVisualize, PlutoVista, ColorSchemes, CairoMakie
	using StaticArrays, LinearAlgebra
	using Revise
	using Printf
	using MultiComponentReactiveMixtureProject
	#GridVisualize.default_plotter!(PlutoVista)
	GridVisualize.default_plotter!(CairoMakie)	
end;

# ╔═╡ 0f102f06-3ff3-4bcc-8892-8d9190a87849
TableOfContents(aside=true)

# ╔═╡ d14462c6-f63b-4a61-a1d9-4bcdb8e30e3d
md"""
Supplementary notebook for:$br
__Transport of heat and mass for reactive gas mixtures in porous media: modeling and application__.

This notebook can be used to reproduce the results shown in Section 4.2 _Thermodiffusion_.
"""

# ╔═╡ 38c3ddb4-b44a-4981-9000-0a1d303bd9ac
md"""
# Introduction
The modeling framework is applied to demonstrate the effect of thermodiffusion. As an example the separation of an equimolar ternary gas mixture of Helium, Argon and Krypton in a closed separation chamber driven by a temperature gradient as presented in [1] is reproduced.
"""

# ╔═╡ 2628cb2d-c1ef-4ad0-8ee4-38e45f864838
md"""
## Model equations
"""

# ╔═╡ 75fddee6-e057-4e91-a239-2033370b00fc
md"""
Reiterating the presented modeling equations from Section 2.3, the species mass balance and thermal energy equations are solved for the non-isothermal ternary gas mixture in the separation chamber:
```math
\begin{align}
    \partial_t \rho_i +\nabla\cdot(\rho_i \vec v) + \nabla \cdot \vec J_i  = r_i(\varrho),\qquad i=1,\dots,n.
    \end{align}   
```


```math
\begin{align}
	\partial_t(\rho h)+\nabla\cdot(\rho h \vec v)+\nabla\cdot\vec Q=\partial_tp.
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

Diffusive thermal energy flux $\vec Q$:
```math
\begin{align}
\vec Q = -\lambda \nabla T + \sum_{i=1}^{n} \left (h_i + RT \widetilde{\mathcal{X}_i}/M_i \right) \vec J_i
\end{align}
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

# ╔═╡ f4286a46-ceb4-40b2-8ce5-7fcd231e3168
md"""
## Modeling domain
The two-dimensional modeling domain corresponding to the separation chamber of 18 cm total length consisting of two chambers of 8 cm height and 2 cm width that are connect by a pipe of 3 cm height and 10 cm length is shown below.
"""

# ╔═╡ 6939978d-9590-407b-80dc-54721c3f672d
md"""
Select problem dimension: $(@bind dim Select([1, 2], default=2))

Select grid refinement level: $(@bind nref Select([0,1,2,3], default=0))

"""

# ╔═╡ e9cb07eb-cfbb-4802-bc7f-6de7a6ad8ac6
md"""
## Simulation parameters
"""

# ╔═╡ f008f30f-0137-4c54-8d4e-a6f589c4a952
md"""
Setup data structure corresponding to the ernary gas mixture consisting of the noble gases:
1) Helium (He)
2) Argon (Ar)
3) Kryton (Kr)
"""

# ╔═╡ e7ca4902-0e14-48ca-bcc6-96b06c85a39d
HeArKr = KinData{}(;
    ng = 3,
    gnames = [:He, :Ar, :Kr],
    Fluids = [He, Ar, Kr],
    gn = Dict(1:3 .=> [:He, :Ar, :Kr]),
    gni = Dict(value => key for (key, value) in Dict(1:3 .=> [:He, :Ar, :Kr])),
	nr = 0,
	rnames = [],
)

# ╔═╡ 4c2be37b-f47e-4841-9ad9-a8c222965e06
md"""
In accordance with [1] the following constant values for the binary Maxwell-Stefan diffusivites are applied:
-  $D_{12}$: He-Ar = 0.756 cm²/s
-  $D_{13}$: He-Kr = 0.6553 cm²/s
-  $D_{23}$: Ar-Kr = 0.1395 cm²/s
"""

# ╔═╡ 4b80f177-9591-49db-8c51-60e37ce2dab5
constant_binary_diff_coeffs = [0.756, 0.6553, 0.1395] * ufac"cm^2/s"

# ╔═╡ 9a7206ae-35f5-430f-b933-58c9494b9f0c
md"""
### Thermodiffusion coefficients
In a multi-component setting the thermal diffusion coefficients always appear in differences of pairs $\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}$, thus a matrix of so called _Newman–Soret diffusivities_ [1] can be defined that is calculated from binary Soret diffusion coefficients:
```math
\begin{align}
\mathcal{A}_{ij} &= \frac{D_i^{\text T}}{\rho_i} - \frac{D_j^{\text T}}{\rho_j} = \mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T} \\
\end{align}
```
"""

# ╔═╡ 3f9eea9d-7c2f-4d83-981a-a243fdf0531a
md"""
Adopting the notation presented in [2] which introduces the _rescaled thermal diffusion ratios_ $\widetilde{\mathcal{X}_i}$:
```math
\begin{align}
\sum_{j=1}^{n} \frac{x_j}{D_{ij}} (\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}) = \sum_{j=1}^{n} x_j \frac{\mathcal{A}_{ij}}{D_{ij}} = \widetilde{\mathcal{X}_i}
\end{align}
```
"""

# ╔═╡ c0478858-4509-4b8c-98c4-919d580f58c5
md"""
Therefore, in accordance with [1] the following constant values for the _Newman–Soret diffusivities_ $\mathcal{A}_{ij}$ are applied:
-  $\mathcal{A}_{12}$: He-Ar = 0.756 cm²/s
-  $\mathcal{A}_{13}$: He-Kr = 0.6553 cm²/s
-  $\mathcal{A}_{23}$: Ar-Kr = 0.1395 cm²/s
"""

# ╔═╡ e9449825-b747-4c89-9cc8-b52ae639e1c0
constant_newman_soret_diff_coeffs = [-0.291, -0.2906, 0.004] * ufac"cm^2/s"

# ╔═╡ 246c0a65-9504-44cc-91ee-d5d6c186d4fe
md"""
Define simulation parameters e.g. the initial composition (equimolar), temperature (300 K) and pressure (1 bara).
"""

# ╔═╡ 7f1d9cf8-7785-48c1-853c-74680188121f
data = ReactorData(
	dim = dim,
	permeable_regions = [1],

	kinpar = HeArKr,
	X0 = [1,1,1] / 3,
	Tamb = 300.0,
	p = 1*ufac"bar",
	
	solve_T_equation = true,
	is_reactive = false,
	include_Soret_Dufour = true,
	#include_Soret_Dufour = false,
	
	γ_τ = [1.0],
	poros= [1.0],

	perm = [1.0]*1.23e-10*ufac"m^2" * 1.0e6,
	
	constant_binary_diff_coeffs = constant_binary_diff_coeffs,
	constant_newman_soret_diff_coeffs = constant_newman_soret_diff_coeffs
)

# ╔═╡ bcfc3138-3cbb-4386-a33b-573b6c39caf9
md"""
## Boundary conditions
"""

# ╔═╡ 0c5e24c0-4aa5-44a2-b2fd-db78795485af
md"""
The modeling domain corresponds to a closed system, so no-flux boundary conditions for species mass apply on all domain boundaries.
A fixed temperature of 300 K is imposed on boundary 4, corresponding to the left hand side chamber while a constant heat flux is imposed on boundary 2, corresponding to the right hand side chamber.
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

# ╔═╡ 78ce253e-49b3-46ee-90d8-54a43369ce9f
md"""
Uncomment the following line to export the time evolution of the solution as an animation.
"""

# ╔═╡ 8b30f68c-9111-4803-b3dc-16e4c440865b
#plotting_movie(filename="Soret_demo_transient.mp4")

# ╔═╡ 8c53810e-330f-4eef-9402-62d31fb5d753
md"""
## Stationary solution
"""

# ╔═╡ a94b7d3e-981a-4e68-9214-471bb854615d
md"""
Using the transient solution at the final time step as initial value for the solution of the stationary system, solve the stationary system and return the stationary solution:
"""

# ╔═╡ dc79ae74-1480-472d-92ff-7e4c36c69215
md"""
Plot the stationary solution corresponding to Figure 8:
"""

# ╔═╡ 65dbb492-4795-44ca-afcb-fb2a2c925d92
md"""
# References
1) Van_Brunt, Alexander; Farrell, Patrick E.; Monroe, Charles W. (2022): Consolidated theory of fluid thermodiffusion. In: AIChE Journal 68 (5), Artikel e17599. DOI: 10.1002/aic.17599                                   .
1) Giovangigli, Vincent (2016): Solutions for Models of Chemically Reacting Mixtures. In: Yoshikazu Giga und Antonin Novotny (Hg.): Handbook of Mathematical Analysis in Mechanics of Viscous Fluids. Cham: Springer International Publishing, S. 1–52.
"""

# ╔═╡ e8425c71-666a-462e-9c4d-fc480810f922
md"""
# Function definitions
"""

# ╔═╡ 56b18561-2d4e-42a8-a363-98c783d0f991
function run_ss(solt,sys)	
	sol_steadystate = VoronoiFVM.solve(
		sys;
		time = solt.t[end],
		inival=solt(solt.t[end]),
		verbose="na"
	)
end

# ╔═╡ bfebc943-28ef-477b-bc15-6cab9d398d92
function grid_1D(;L=18*ufac"cm", n=20)
	X=(0:(L/n):L)
	simplexgrid(X)
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

# ╔═╡ b55537bf-9982-4997-8a2a-1972127bdd86
let
	vis = GridVisualizer(resolution=(450,225))
	gridplot!(vis, grid_2D(nref=0), zoom=2, linewidth=1)
	reveal(vis)
	#fn = "../../../data/out/2024-03-25/Soret_domain_nref0.pdf"
	#GridVisualize.save(fn, vis)
end

# ╔═╡ 138b4d12-af0f-4a1c-a4aa-4f55e6038664
if dim == 1
	grid = grid_1D()
	const Γ_left = 1
	const Γ_right = 2
elseif dim == 2
	grid = grid_2D(nref=nref)
	const Γ_left = 4
	const Γ_right = 2
end;

# ╔═╡ 93970c02-91c6-499a-9318-f7f632604bb5
function bcondition(f,u,bnode,data)
	(;p,ip,iT,Tamb,inflow_boundaries,dt_hf_enth)=data

	eps_=1/1e-4
	boundary_robin!(f,u,bnode, species=iT,region=Γ_left,value=Tamb*eps_,factor=eps_)
	
	if dim==2
		heatflux_right = 11.3*ufac"W/m^2"
	elseif dim==1
		heatflux_right = 25*ufac"W/m^2"
	end
	boundary_neumann!(f,u,bnode, species=iT, region=Γ_right, 
		value=ramp(bnode.time; du=(0,heatflux_right), dt=dt_hf_enth) )
end

# ╔═╡ 1e51701d-a893-4056-8336-a3772b85abe4
function setup_run_sim(grid, data)
	(;ng, ip, iT, Tamb, p, X0, outflow_boundaries) = data
	sys=VoronoiFVM.System(
		grid;
		data=data,
		flux=DMS_flux,
		reaction=DMS_reaction,
		storage=DMS_storage,
		bcondition=bcondition,
		assembly=:edgewise,
		unknown_storage=:dense
	)
	
	enable_species!(sys; species=collect(1:(ng+2))) # gas phase species xi, ptotal, T
		
	inival=unknowns(sys)
	inival[ip,:].=p
	inival[iT,:].=Tamb

	for i=1:ng
		inival[i,:] .= X0[i]
	end

	control = SolverControl(nothing, sys;)
	control.handle_exceptions=true
	if dim == 2
		control.Δu_opt=1000
		control.Δt_max=20.0
		times=[0,3000.0]
	elseif dim == 1
		control.Δu_opt=10
		control.Δt_max=10.0		
		times=[0,3000.0]
	end
		

	
	sol=VoronoiFVM.solve(sys;inival=inival,times,control,log=true, verbose="nae")
	return (sol,sys)
end

# ╔═╡ 035d4123-7092-4429-8cfd-1e5926e84493
solt, sys = setup_run_sim(grid, data)

# ╔═╡ 5a0900cc-df10-4176-b903-358b3e00415c
md"""
Move the slider to change $t$: $(if isa(solt, TransientSolution)
	@bind t PlutoUI.Slider(solt.t,show_value=true,default=solt.t[end])	
end)
"""

# ╔═╡ 076b4a28-be0f-46f0-9857-e6f886c4b118
md"""
Plot the transient solution at time $t=$ $(round(t)):
"""

# ╔═╡ ad68e43e-df7e-4a06-a697-fa5824f54d3e
if isa(solt, TransientSolution)
	sol = solt(t)
end;

# ╔═╡ 4296aa28-9f52-4d40-a968-ee583ffc7d3c
sol_ss = run_ss(solt,sys)

# ╔═╡ 9ba456c7-f6ed-49c5-9878-d9e0deb96384
let
	(; gni,iT,ip) = data;
	
	vis =GridVisualizer(layout=(2,2), resolution=(1000,600))
	scalarplot!(vis[1,1], grid, sol[gni[:He],:], colormap=:coolwarm, title = "He molar fraction",)
	scalarplot!(vis[1,2], grid, sol[gni[:Kr],:], colormap=:coolwarm, title = "Kr molar fraction",)
	scalarplot!(vis[2,1], grid, sol[iT,:], colormap=:gist_heat, title = "Temperature (K)",)

	#f = map((x, y) -> maximum(sol[ip,:]), grid)
	pmax = round(maximum(sol[ip,:]))
	f = map((x, y) -> pmax, grid)
	scalarplot!(vis[2,2], grid, f, title = "Pressure (Pa)",)
	reveal(vis)
end

# ╔═╡ 2dc80ec7-f9ad-4e98-92de-faeac2beb368
let
	(; gni,iT,ip) = data;
	
	#vis =GridVisualizer(layout=(2,2), resolution=(1000,600))
	vis =GridVisualizer(layout=(2,2), resolution=(900,450))
	
	scalarplot!(vis[1,1], grid, sol_ss[gni[:He],:], limits=(0.32, 0.345), levels=51, linewidth=1, colorbarticks=[0.32,0.325,0.33,0.335,0.34,0.345], colormap=:coolwarm, xlabel="X (cm)", ylabel="Y (cm)", title = "He molar fraction")

	scalarplot!(vis[1,2], grid, sol_ss[gni[:Kr],:], limits=(0.325, 0.34), levels=28, linewidth=1, colorbarticks=[0.325,0.33,0.335,0.34], colormap=:coolwarm, xlabel="X (cm)", ylabel="Y (cm)", title = "Kr molar fraction")

	scalarplot!(vis[2,1], grid, sol_ss[iT,:], limits=(300, 400), levels=51, linewidth=1, colorbarticks=[300,320,340,360,380,400], colormap=:gist_heat, xlabel="X (cm)", ylabel="Y (cm)", title = "Temperature (K)")
	
	pmax = round(maximum(sol_ss[ip,:]))
	f = map((x, y) -> pmax, grid)
	scalarplot!(vis[2,2], grid, f, xlabel="X (cm)", ylabel="Y (cm)", title = "Pressure (Pa)",)
	fn = "../../../data/out/2024-03-25/Soret_result_highres_nref1.pdf"
	
	reveal(vis)
	#GridVisualize.save(fn, vis)
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
# ╟─f4286a46-ceb4-40b2-8ce5-7fcd231e3168
# ╟─b55537bf-9982-4997-8a2a-1972127bdd86
# ╟─6939978d-9590-407b-80dc-54721c3f672d
# ╠═138b4d12-af0f-4a1c-a4aa-4f55e6038664
# ╟─e9cb07eb-cfbb-4802-bc7f-6de7a6ad8ac6
# ╟─f008f30f-0137-4c54-8d4e-a6f589c4a952
# ╠═e7ca4902-0e14-48ca-bcc6-96b06c85a39d
# ╟─4c2be37b-f47e-4841-9ad9-a8c222965e06
# ╠═4b80f177-9591-49db-8c51-60e37ce2dab5
# ╟─9a7206ae-35f5-430f-b933-58c9494b9f0c
# ╟─3f9eea9d-7c2f-4d83-981a-a243fdf0531a
# ╟─c0478858-4509-4b8c-98c4-919d580f58c5
# ╠═e9449825-b747-4c89-9cc8-b52ae639e1c0
# ╟─246c0a65-9504-44cc-91ee-d5d6c186d4fe
# ╠═7f1d9cf8-7785-48c1-853c-74680188121f
# ╟─bcfc3138-3cbb-4386-a33b-573b6c39caf9
# ╟─0c5e24c0-4aa5-44a2-b2fd-db78795485af
# ╠═93970c02-91c6-499a-9318-f7f632604bb5
# ╟─13cfe122-eea8-4bbb-aaa0-5bcc74a247d1
# ╟─a4fd9977-577e-45c9-a5c3-c0cd8a5dd012
# ╟─7b0a84b5-60d3-4dd9-89e9-29c88282cb25
# ╟─3a063481-d447-4bf5-9c49-ecde37a0fcea
# ╠═035d4123-7092-4429-8cfd-1e5926e84493
# ╟─076b4a28-be0f-46f0-9857-e6f886c4b118
# ╟─9ba456c7-f6ed-49c5-9878-d9e0deb96384
# ╟─5a0900cc-df10-4176-b903-358b3e00415c
# ╟─78ce253e-49b3-46ee-90d8-54a43369ce9f
# ╠═8b30f68c-9111-4803-b3dc-16e4c440865b
# ╟─ad68e43e-df7e-4a06-a697-fa5824f54d3e
# ╟─8c53810e-330f-4eef-9402-62d31fb5d753
# ╟─a94b7d3e-981a-4e68-9214-471bb854615d
# ╠═4296aa28-9f52-4d40-a968-ee583ffc7d3c
# ╟─dc79ae74-1480-472d-92ff-7e4c36c69215
# ╟─2dc80ec7-f9ad-4e98-92de-faeac2beb368
# ╟─65dbb492-4795-44ca-afcb-fb2a2c925d92
# ╟─e8425c71-666a-462e-9c4d-fc480810f922
# ╠═56b18561-2d4e-42a8-a363-98c783d0f991
# ╠═1e51701d-a893-4056-8336-a3772b85abe4
# ╠═bfebc943-28ef-477b-bc15-6cab9d398d92
# ╠═ff501186-d8a0-4666-8991-fe576f8ff6ad
# ╠═3c4ae416-a47b-451d-b870-fa10166d97de
