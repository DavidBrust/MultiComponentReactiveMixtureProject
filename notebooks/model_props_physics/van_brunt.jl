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
	Pkg.activate(joinpath(@__DIR__,"../.."))
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

# ╔═╡ f4286a46-ceb4-40b2-8ce5-7fcd231e3168
md"""
# Modelling Domain
$(LocalResource("img/vanbrunt_domain.png", :width=> 400))
"""

# ╔═╡ 6939978d-9590-407b-80dc-54721c3f672d
md"""
Select problem dimension: $(@bind dim Select([1, 2], default=2))

Select grid refinement level: $(@bind nref Select([0,1,2,3], default=0))

"""

# ╔═╡ ff501186-d8a0-4666-8991-fe576f8ff6ad
function grid_2D(;nref=0)
    #X = collect(0:0.5:18)*ufac"cm"
	X = collect(0:(1/2^nref):18)*ufac"cm"
    Y = collect(0:(0.5/2^nref):8)*ufac"cm"
    grid = simplexgrid(X, Y)
    
	rect!(grid, [4, 0]*ufac"cm", [14, 2.5]*ufac"cm"; region = 2)
	rect!(grid, [4, 5.5]*ufac"cm", [14, 8]*ufac"cm"; region = 2, bregions= [3, 3, 3, 3])
	bfacemask!(grid, [0, 0]*ufac"cm", [4, 0]*ufac"cm", 4)
	bfacemask!(grid, [0, 8]*ufac"cm", [4, 8]*ufac"cm", 4)
	
	subgrid(grid, [1])
end

# ╔═╡ bfebc943-28ef-477b-bc15-6cab9d398d92
function grid_1D(;L=18*ufac"cm", n=20)
	X=(0:(L/n):L)
	simplexgrid(X)
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

# ╔═╡ b55537bf-9982-4997-8a2a-1972127bdd86
gridplot(grid, resolution=(500,300), zoom=2)

# ╔═╡ 9fde21eb-1112-40ec-8fe2-0c8a12c9926d
md"""
# M-S flux-force relationship
```math
\begin{align}
	\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} + x_i \sum_{j=1}^{n} \frac{x_j}{D_{ij}} (\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}) \frac{\nabla T}{T} \right) &= -\sum_{j=1 \atop j \neq i}^{n} \frac{w_j \vec J_i-w_i \vec \PhJ{D_{ij} M_i M_j} \\
\end{align}
```
"""

# ╔═╡ 9a7206ae-35f5-430f-b933-58c9494b9f0c
md"""
## Thermo-diffusion coefficients
In a multi-component setting the thermal diffusion coefficients always appear in differences of pairs $\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}$, thus a matrix of so called _Newman–Soret diffusivities_ [1] can be defined that is calculated from binary Soret diffusion coefficients:
```math
\begin{align}
\mathcal{A}_{ij} &= \frac{D_i^{\text T}}{\rho_i} - \frac{D_j^{\text T}}{\rho_j} = \mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T} \\
\end{align}
```
 $(LocalResource("img/vanbrunt_newman_soret_coeff_regularized.png", :width=>500))

The _thermal diffusion factors_ $\alpha_{ij}$ are commonly used which are related to the _Newman–Soret diffusivities_ via [1]:
```math
\begin{align}
\alpha_{ij} = \frac{\mathcal{A}_{ij} }{D_{ij}}
\end{align}
```
"""

# ╔═╡ 2fcc681c-1b79-46cc-b96f-20143f8e331e
md"""
$(LocalResource("img/vanbrunt_thermo_diff_factors.png", :width=>500))
"""

# ╔═╡ 3f9eea9d-7c2f-4d83-981a-a243fdf0531a
md"""
Adopting the notation presented in [2] which introduces the _rescaled thermal diffusion ratios_ $\widetilde{\mathcal{X}_i}$:
```math
\begin{align}
\sum_{j=1}^{n} \frac{x_j}{D_{ij}} (\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}) = \sum_{j=1}^{n} x_j \alpha_{ij} = \widetilde{\mathcal{X}_i}
\end{align}
```
"""

# ╔═╡ 94dd7674-751f-4128-b5eb-303fb9693c22
md"""
The Maxwell-Stefan flux-force relationship including the Soret effect thus reads:

```math
\begin{align}
	\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} + x_i \widetilde{\mathcal{X}_i} \frac{\nabla T}{T} \right) &= -\sum_{j=1 \atop j \neq i}^{n} \frac{w_j \vec J_i-w_i \vec \PhJ{D_{ij} M_i M_j} \\
\end{align}
```
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

# ╔═╡ e9cb07eb-cfbb-4802-bc7f-6de7a6ad8ac6
md"""
# Specify Sim Data
"""

# ╔═╡ 7f1d9cf8-7785-48c1-853c-74680188121f
data = ReactorData(
	dim = dim,

	#dt_hf_enth = (0, 200),
	kinpar = HeArKr,
	X0 = [1,1,1] / 3,
	Tamb = 300.0,
	p = 1*ufac"bar",
	
	solve_T_equation = true,
	is_reactive = false,
	#constant_properties = true,
	include_Soret_Dufour = true,
	
	γ_τ = 1.0,
	poros=1.0,

	perm = 1.23e-10*ufac"m^2" * 1.0e6,
	# 1) He, 2) Ar, 3) Kr
	# D12:= He-Ar, D13:= He-Kr, D23:= Ar-Kr
	constant_binary_diff_coeffs = [0.756, 0.6553, 0.1395] * ufac"cm^2/s",
	constant_newman_soret_diff_coeffs = [-0.291, -0.2906, 0.004] * ufac"cm^2/s",
	
	outlet_boundaries=[Γ_left],
	inlet_boundaries=[Γ_right]	
)

# ╔═╡ 8b30f68c-9111-4803-b3dc-16e4c440865b
#plotting_movie(filename="Soret_demo_transient.mp4")

# ╔═╡ 7b0a84b5-60d3-4dd9-89e9-29c88282cb25
md"""
# Transient Solution
"""

# ╔═╡ 93970c02-91c6-499a-9318-f7f632604bb5
function bcondition(f,u,bnode,data)
	(;p,ip,iT,Tamb,inlet_boundaries,dt_hf_enth)=data

	#boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_left,value=Tamb)
	eps_=1/1e-4
	boundary_robin!(f,u,bnode, species=iT,region=Γ_left,value=Tamb*eps_,factor=eps_)
	
	#boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_right,value=ramp(bnode.time; du=(Tamb,Tamb+100), dt=(0,5) ) )

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
	(;ng, ip, iT, Tamb, p, X0, outlet_boundaries) = data
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
		

	
	sol=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="aen",log=true)
	return (sol,sys)
end

# ╔═╡ 035d4123-7092-4429-8cfd-1e5926e84493
solt, sys = setup_run_sim(grid, data);

# ╔═╡ cf1d3089-a0d2-445d-b004-571776f1c9a0
if isa(solt, TransientSolution)
	@bind t PlutoUI.Slider(solt.t,show_value=true,default=solt.t[end])	
end

# ╔═╡ ad68e43e-df7e-4a06-a697-fa5824f54d3e
if isa(solt, TransientSolution)
	sol = solt(t)
end;

# ╔═╡ 07e97ba1-357a-4de8-ad5a-64e7a27b0cb8
md"""
| [1] | This work |
|:----------:|:----------:|
| $(LocalResource("img/vanbrunt_result_xHe_2.png", :width=> 250)) | __He__ molar fraction |
| $(LocalResource("img/vanbrunt_result_xHe.png", :width=> 250)) 	| $((; gni) = data; scalarplot(grid, sol[gni[:He],:], resolution=(300,200), colormap=:coolwarm, zoom=2.5) ) |
| $(LocalResource("img/vanbrunt_result_xKr_2.png", :width=> 250))     | __Kr__ molar fraction |
| $(LocalResource("img/vanbrunt_result_xKr.png", :width=> 250))  	| $((;gni) = data; scalarplot(grid, sol[gni[:Kr],:], resolution=(300,200), colormap=:coolwarm, zoom=2.5) ) |
| $(LocalResource("img/vanbrunt_result_T_2.png", :width=> 250))  	| __Temperature (K)__ |
| $(LocalResource("img/vanbrunt_result_T.png", :width=> 250))  	| $((;iT) = data; scalarplot(grid, sol[iT,:], resolution=(300,200), colormap=:gist_heat, zoom=2.5) ) |
"""

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

# ╔═╡ 5a517383-c597-4fa5-b7dc-441e1952cb97
plothistory(solt, Plotter=PlutoVista)

# ╔═╡ 56b18561-2d4e-42a8-a363-98c783d0f991
function run_ss(solt,sys)	
	sol_steadystate = VoronoiFVM.solve(
		sys;
		time = solt.t[end],
		inival=solt(solt.t[end]),
		verbose="na"
	)
end

# ╔═╡ a31d1583-8b59-4f30-9fae-cc4c3ceea1cd
sol

# ╔═╡ 8c53810e-330f-4eef-9402-62d31fb5d753
md"""
# Steady-state solution
"""

# ╔═╡ 4296aa28-9f52-4d40-a968-ee583ffc7d3c
sol_ss = run_ss(solt,sys)

# ╔═╡ b0007963-cf73-49b5-92a7-5a2ef1bbd2f5
md"""
| [1] | This work |
|:----------:|:----------:|
| $(LocalResource("img/vanbrunt_result_xHe_2.png", :width=> 250)) | __He__ molar fraction |
| $(LocalResource("img/vanbrunt_result_xHe.png", :width=> 250)) 	| $(scalarplot(grid, sol_ss[gni[:He],:], limits=(0.32, 0.345), levels=51, linewidth=1, colorbarticks=[0.32,0.325,0.33,0.335,0.34,0.345], resolution=(400,300), colormap=:coolwarm, zoom=2.5) ) |
| $(LocalResource("img/vanbrunt_result_xKr_2.png", :width=> 250))     | __Kr__ molar fraction |
| $(LocalResource("img/vanbrunt_result_xKr.png", :width=> 250))  	| $( scalarplot(grid, sol_ss[gni[:Kr],:], limits=(0.325, 0.34), levels=28, linewidth=1, colorbarticks=[0.325,0.33,0.335,0.34], resolution=(400,300), colormap=:coolwarm, zoom=2.5) ) |
| $(LocalResource("img/vanbrunt_result_T_2.png", :width=> 250))  	| __Temperature (K)__ |
| $(LocalResource("img/vanbrunt_result_T.png", :width=> 250))  	| $(scalarplot(grid, sol_ss[iT,:], limits=(300, 400), levels=51, linewidth=1, colorbarticks=[300,320,340,360,380,400], resolution=(400,300), colormap=:gist_heat, zoom=2.5) ) |
"""

# ╔═╡ 65dbb492-4795-44ca-afcb-fb2a2c925d92
md"""
# References
1) Van_Brunt, Alexander; Farrell, Patrick E.; Monroe, Charles W. (2022): Consolidated theory of fluid thermodiffusion. In: AIChE Journal 68 (5), Artikel e17599. DOI: 10.1002/aic.17599                         .
1) Giovangigli, Vincent (2016): Solutions for Models of Chemically Reacting Mixtures. In: Yoshikazu Giga und Antonin Novotny (Hg.): Handbook of Mathematical Analysis in Mechanics of Viscous Fluids. Cham: Springer International Publishing, S. 1–52.
"""

# ╔═╡ e8425c71-666a-462e-9c4d-fc480810f922
md"""
# Auxiliary functions
"""

# ╔═╡ 3c4ae416-a47b-451d-b870-fa10166d97de
function plotting_movie(; filename = "plotting_video.gif", Plotter = default_plotter())
    vis =GridVisualizer(layout=(2,2), resolution=(1000,600))

	(; gni,iT,ip) = data;
    movie(vis; file = filename) do vis
        for t in solt.t
            #f = map((x, y) -> sin(x - t) * cos(y - t), grid)
            #g = map((x, y) -> sin(t) * sin(x) * cos(y), grid)
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
# ╟─f4286a46-ceb4-40b2-8ce5-7fcd231e3168
# ╟─6939978d-9590-407b-80dc-54721c3f672d
# ╠═138b4d12-af0f-4a1c-a4aa-4f55e6038664
# ╠═b55537bf-9982-4997-8a2a-1972127bdd86
# ╠═ff501186-d8a0-4666-8991-fe576f8ff6ad
# ╠═bfebc943-28ef-477b-bc15-6cab9d398d92
# ╟─9fde21eb-1112-40ec-8fe2-0c8a12c9926d
# ╟─9a7206ae-35f5-430f-b933-58c9494b9f0c
# ╟─2fcc681c-1b79-46cc-b96f-20143f8e331e
# ╟─3f9eea9d-7c2f-4d83-981a-a243fdf0531a
# ╟─94dd7674-751f-4128-b5eb-303fb9693c22
# ╠═e7ca4902-0e14-48ca-bcc6-96b06c85a39d
# ╟─e9cb07eb-cfbb-4802-bc7f-6de7a6ad8ac6
# ╠═7f1d9cf8-7785-48c1-853c-74680188121f
# ╟─07e97ba1-357a-4de8-ad5a-64e7a27b0cb8
# ╟─9ba456c7-f6ed-49c5-9878-d9e0deb96384
# ╠═8b30f68c-9111-4803-b3dc-16e4c440865b
# ╠═cf1d3089-a0d2-445d-b004-571776f1c9a0
# ╠═ad68e43e-df7e-4a06-a697-fa5824f54d3e
# ╠═5a517383-c597-4fa5-b7dc-441e1952cb97
# ╟─7b0a84b5-60d3-4dd9-89e9-29c88282cb25
# ╠═035d4123-7092-4429-8cfd-1e5926e84493
# ╠═93970c02-91c6-499a-9318-f7f632604bb5
# ╠═1e51701d-a893-4056-8336-a3772b85abe4
# ╠═56b18561-2d4e-42a8-a363-98c783d0f991
# ╠═a31d1583-8b59-4f30-9fae-cc4c3ceea1cd
# ╟─8c53810e-330f-4eef-9402-62d31fb5d753
# ╟─4296aa28-9f52-4d40-a968-ee583ffc7d3c
# ╟─b0007963-cf73-49b5-92a7-5a2ef1bbd2f5
# ╟─65dbb492-4795-44ca-afcb-fb2a2c925d92
# ╟─e8425c71-666a-462e-9c4d-fc480810f922
# ╠═3c4ae416-a47b-451d-b870-fa10166d97de
