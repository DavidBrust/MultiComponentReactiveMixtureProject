### A Pluto.jl notebook ###
# v0.19.38

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
	using ExtendableGrids, GridVisualize, PlutoVista
	using StaticArrays
	using Revise
	using Printf
	using MultiComponentReactiveMixtureProject
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 0f102f06-3ff3-4bcc-8892-8d9190a87849
TableOfContents(aside=true)

# ╔═╡ f4286a46-ceb4-40b2-8ce5-7fcd231e3168
md"""
# Modelling Domain
$(LocalResource("img/vanbrunt_domain.png"))
"""

# ╔═╡ 6939978d-9590-407b-80dc-54721c3f672d
md"""
Select problem dimension: $(@bind dim Select([1, 2], default=1))
"""

# ╔═╡ b55537bf-9982-4997-8a2a-1972127bdd86
#gridplot(grid_2D())

# ╔═╡ ff501186-d8a0-4666-8991-fe576f8ff6ad
function grid_2D()
    #X = collect(0:0.5:18)*ufac"cm"
	X = collect(0:1:18)*ufac"cm"
    Y = collect(0:0.5:8)*ufac"cm"
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
	grid = grid_2D()
	const Γ_left = 4
	const Γ_right = 2
end;

# ╔═╡ 9fde21eb-1112-40ec-8fe2-0c8a12c9926d
md"""
# M-S flux-force relationship
```math
\begin{align}
	\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} + x_i \sum_{j=1}^{n} \frac{x_j}{D_{ij}} (\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}) \frac{\nabla T}{T} \right) &= -\sum_{j=1 \atop j \neq i}^{n} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
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
 $(LocalResource("img/vanbrunt_newman_soret_coeff_regularized.png"))

The _thermal diffusion factors_ $\alpha_{ij}$ are commonly used which are related to the _Newman–Soret diffusivities_ via [1]:
```math
\begin{align}
\alpha_{ij} = \frac{\mathcal{A}_{ij} }{D_{ij}}
\end{align}
```
"""

# ╔═╡ 2fcc681c-1b79-46cc-b96f-20143f8e331e
md"""
$(LocalResource("img/vanbrunt_thermo_diff_factors.png"))
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
	\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} + x_i \widetilde{\mathcal{X}_i} \frac{\nabla T}{T} \right) &= -\sum_{j=1 \atop j \neq i}^{n} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
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

# ╔═╡ 7f1d9cf8-7785-48c1-853c-74680188121f
data = ReactorData(
	#dim = dim,
	dim = 1,
	kinpar = HeArKr,
	X0 = [1,1,1] / 3,
	Tamb = 300,
	p = 1*ufac"bar",
	γ_τ = 1.0,
	#rhos=5.0*ufac"kg/m^3", # low value for solid density -> low thermal inertia
	solve_T_equation = true,
	is_reactive = false,
	constant_properties = true,
	include_Soret_Dufour = true,
	rhos=1.0*ufac"kg/m^3", # low value for solid density -> low thermal inertia
	#poros=0.99,
	#poros=1.0,
	# 1) He, 2) Ar, 3) Kr
	# D12:= He-Ar, D13:= He-Kr, D23:= Ar-Kr
	constant_binary_diff_coeffs = [0.756, 0.6553, 0.1395] * ufac"cm^2/s",
	constant_newman_soret_diff_coeffs = [-0.291, -0.2906, 0.004] * ufac"cm^2/s",
	#constant_species_thermal_conductivities = [
	#	thermcond_gas(He, 300),
	#	thermcond_gas(Ar, 300),
	#	thermcond_gas(Kr, 300)
	#],
	#constant_species_viscosities = [
	#	dynvisc_gas(He, 300),
	#	dynvisc_gas(Ar, 300),
	#	dynvisc_gas(Kr, 300)
	#],
	#constant_newman_soret_diff_coeffs = [-0.3012, -0.2804, -0.0099] * ufac"cm^2/s"
	#outlet_boundaries=[],
	#inlet_boundaries=[]
	outlet_boundaries=[Γ_left],
	inlet_boundaries=[Γ_right]
	
	
)

# ╔═╡ 0cc8d870-e13e-4ca0-99e9-48a374939c6b
md"""
Results from [1]:

$(LocalResource("img/vanbrunt_result_xHe_2.png", :width=> 400))
$(LocalResource("img/vanbrunt_result_xHe.png", :width=> 400))

$(LocalResource("img/vanbrunt_result_xKr_2.png", :width=> 400))
$(LocalResource("img/vanbrunt_result_xKr.png", :width=> 400))
"""

# ╔═╡ 5ae4d7c6-11b2-4f8f-8234-9b658afaa83b
md"""
Results from this work:
1) __He__ molar fraction
2) __Kr__ molar fraction
""";

# ╔═╡ 93970c02-91c6-499a-9318-f7f632604bb5
function bcondition(f,u,bnode,data)
	(;p,ip,iT,Tamb,inlet_boundaries,dt_hf_enth)=data

	boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_left,value=300)
	boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_right,value=ramp(bnode.time; du=(300,400), dt=(0,5) ) )
	
	#heatflux_right = 135.0
	#boundary_neumann!(f,u,bnode, species=iT,region=Γ_right,value=ramp(bnode.time; du=(0,heatflux_right), dt=dt_hf_enth ) )
end

# ╔═╡ 1e51701d-a893-4056-8336-a3772b85abe4
function setup_run_sim(grid, data)
	(;ng, ip, iT, Tamb, p, X0) = data
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
	control.Δu_opt=100
	
	times=[0,200]
	sol=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="aen",log=true)
	#sol=VoronoiFVM.solve(sys;inival=inival,verbose="an",log=true)
	return (sol,sys)
end

# ╔═╡ 035d4123-7092-4429-8cfd-1e5926e84493
solt, sys = setup_run_sim(grid, data);

# ╔═╡ c995a528-fa98-4860-9bf3-648af15693e9
let
	inflow_rate, outflow_rate, reaction_rate, stored_amount, I_in, I_out, I_reac = BoundaryFlows_Integrals(solt, sys, data)
	(;ng, gn, gni, iT, ip) = data

	#k=gni[:H2]
	#k=gni[:CO]
	k=iT
	#k=ip

	if k in 1:ng
		name = gn[k]
	elseif k == ip
		name = "Total Mass"
	elseif k == iT
		name = "Enthalpy Flow"
	end
	
	@printf "%s In: %2.2e \t Out: %2.2e \t React: %2.2e \nIn - Out: %2.4e \nStorage tEnd -t0: %2.4e" name I_in[k] I_out[k] I_reac[k] (I_in+I_out-I_reac)[k] (stored_amount[end]-stored_amount[1])[k]

	vis=GridVisualizer(resolution=(500,600), layout=(4,1), xlabel="Time / s", ylabel="Inflow / Outflow / Reaction Rate")

	function plot_flows!(k,vis)
		scalarplot!(vis, solt.t[2:end], stack(inflow_rate, dims=1)[:,k], label="Inflow rate")
		scalarplot!(vis, solt.t[2:end], stack(outflow_rate, dims=1)[:,k], label="Outflow rate", color=:red, clear=false)	
		scalarplot!(vis, solt.t[2:end], stack(-reaction_rate, dims=1)[:,k], label="Reaction rate",  color=:blue, clear=false)
		#scalarplot!(vis, solt.t[2:end], stack(stored_amount, dims=1)[:,k], label="Stored amount", color=:green, clear=false, )
	end

	plot_flows!(ip,vis[1,1])
	plot_flows!(gni[:He],vis[2,1])
	plot_flows!(gni[:Ar],vis[3,1])
	plot_flows!(gni[:Kr],vis[4,1])
	
	reveal(vis)
end

# ╔═╡ cf1d3089-a0d2-445d-b004-571776f1c9a0
if isa(solt, TransientSolution)
	@bind t Slider(solt.t,show_value=true,default=solt.t[end])	
end

# ╔═╡ ad68e43e-df7e-4a06-a697-fa5824f54d3e
if isa(solt, TransientSolution)
	sol = solt(t)
end;

# ╔═╡ e4486776-8e7a-4590-b10d-1b797396dd39
MultiComponentReactiveMixtureProject._checkinout(sol,sys,data)

# ╔═╡ 728abf6f-3ff1-45c8-8b10-213f08b1b4dd
let
	(;iT, ip, gni, gn) = data
	ng=ngas(data)
	if dim==1
		cs = [:black, :red, :blue]
		vis=GridVisualizer(resolution=(600,400), layout=(3,1))
		for i=1:ng
			scalarplot!(vis[1,1], grid, sol[i,:], clear=false, label=gn[i], color=cs[i])
		end
		
		scalarplot!(vis[2,1], grid, sol[iT,:], clear=false, label="Temperature / K")
		scalarplot!(vis[3,1], grid, sol[ip,:], clear=false, label="Pressure / Pa")
	elseif dim==2
		vis=GridVisualizer(resolution=(600,600), layout=(2,1), zoom=2)
		scalarplot!(vis[1,1], grid, sol[gni[:He],:])
		scalarplot!(vis[2,1], grid, sol[gni[:Kr],:])
		
	end
	reveal(vis)
end

# ╔═╡ 7fdc10ec-2d72-4656-a094-e9b2b1c54ecf
let
	(;iT, ip, gni, gn) = data
	vis=GridVisualizer(resolution=(600,600), layout=(2,1), zoom=2)
	scalarplot!(vis[1,1], grid, sol[iT,:])
	scalarplot!(vis[2,1], grid, sol[ip,:])
	reveal(vis)
end

# ╔═╡ 5a517383-c597-4fa5-b7dc-441e1952cb97
plothistory(solt)

# ╔═╡ 56b18561-2d4e-42a8-a363-98c783d0f991
function run_ss(solt,sys)
	control = SolverControl(nothing, sys)
	
	sol_steadystate = VoronoiFVM.solve(
		sys;
		time = 100.0,
		inival=solt(solt.t[end]),
		control,
		verbose="na"
	)
end

# ╔═╡ 8c53810e-330f-4eef-9402-62d31fb5d753
md"""
# Steady-state solution
"""

# ╔═╡ 4296aa28-9f52-4d40-a968-ee583ffc7d3c
sol_ss = run_ss(solt,sys)

# ╔═╡ 206c143d-0af5-4e48-9fd9-f6e7a15e5083
md"""
1) He
1) Kr
1) Temperature
"""

# ╔═╡ 25a57aa7-7b67-4733-b631-994af5118134
let
	(;iT, ip, gni, gn) = data
	ng=ngas(data)
	if dim==1
		cs = [:black, :red, :blue]
		vis=GridVisualizer(resolution=(600,400), layout=(3,1))
		for i=1:ng
			scalarplot!(vis[1,1], grid, sol_ss[i,:], clear=false, label=gn[i], color=cs[i])
		end
		
		scalarplot!(vis[2,1], grid, sol_ss[iT,:], clear=false, label="Temperature / K")
		scalarplot!(vis[3,1], grid, sol_ss[ip,:], clear=false, label="Pressure / Pa")
	elseif dim==2
		vis=GridVisualizer(resolution=(600,600), layout=(3,1), zoom=2.5)
		scalarplot!(vis[1,1], grid, sol_ss[gni[:He],:])
		scalarplot!(vis[2,1], grid, sol_ss[gni[:Kr],:])
		scalarplot!(vis[3,1], grid, sol_ss[iT,:])
		
	end
	reveal(vis)
end

# ╔═╡ 65dbb492-4795-44ca-afcb-fb2a2c925d92
md"""
# References
1) Van‐Brunt, Alexander; Farrell, Patrick E.; Monroe, Charles W. (2022): Consolidated theory of fluid thermodiffusion. In: AIChE Journal 68 (5), Artikel e17599. DOI: 10.1002/aic.17599             .
1) Giovangigli, Vincent (2016): Solutions for Models of Chemically Reacting Mixtures. In: Yoshikazu Giga und Antonin Novotny (Hg.): Handbook of Mathematical Analysis in Mechanics of Viscous Fluids. Cham: Springer International Publishing, S. 1–52.
"""

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
# ╟─7f1d9cf8-7785-48c1-853c-74680188121f
# ╠═e4486776-8e7a-4590-b10d-1b797396dd39
# ╠═c995a528-fa98-4860-9bf3-648af15693e9
# ╟─0cc8d870-e13e-4ca0-99e9-48a374939c6b
# ╟─5ae4d7c6-11b2-4f8f-8234-9b658afaa83b
# ╠═728abf6f-3ff1-45c8-8b10-213f08b1b4dd
# ╠═7fdc10ec-2d72-4656-a094-e9b2b1c54ecf
# ╠═cf1d3089-a0d2-445d-b004-571776f1c9a0
# ╠═ad68e43e-df7e-4a06-a697-fa5824f54d3e
# ╠═5a517383-c597-4fa5-b7dc-441e1952cb97
# ╠═035d4123-7092-4429-8cfd-1e5926e84493
# ╠═93970c02-91c6-499a-9318-f7f632604bb5
# ╠═1e51701d-a893-4056-8336-a3772b85abe4
# ╠═56b18561-2d4e-42a8-a363-98c783d0f991
# ╟─8c53810e-330f-4eef-9402-62d31fb5d753
# ╠═4296aa28-9f52-4d40-a968-ee583ffc7d3c
# ╟─206c143d-0af5-4e48-9fd9-f6e7a15e5083
# ╟─25a57aa7-7b67-4733-b631-994af5118134
# ╟─65dbb492-4795-44ca-afcb-fb2a2c925d92
