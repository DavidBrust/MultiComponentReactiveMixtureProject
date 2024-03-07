### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

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
	using MultiComponentReactiveMixtureProject
end;

# ╔═╡ 0f102f06-3ff3-4bcc-8892-8d9190a87849
TableOfContents(aside=true)

# ╔═╡ f4286a46-ceb4-40b2-8ce5-7fcd231e3168
md"""
# Modelling Domain
$(LocalResource("img/vanbrunt_domain.png"))
"""

# ╔═╡ bfebc943-28ef-477b-bc15-6cab9d398d92
function grid_1D(;L=18*ufac"cm", n=20)
	X=(0:(L/n):L)
	simplexgrid(X)
end

# ╔═╡ 138b4d12-af0f-4a1c-a4aa-4f55e6038664
grid = grid_1D();

# ╔═╡ 63f0afeb-d0fc-4d58-8e24-7f8ccc587657
begin
	const Γ_left = 1
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
	# 1) He, 2) Ar, 3) Kr
	# D12:= He-Ar, D13:= He-Kr, D23:= Ar-Kr
	constant_binary_diff_coeffs = [0.756, 0.6553, 0.1395] * ufac"cm^2/s",
	constant_newman_soret_diff_coeffs = [-0.291, -0.2906, 0.004] * ufac"cm^2/s"
	#constant_newman_soret_diff_coeffs = [-0.3012, -0.2804, -0.0099] * ufac"cm^2/s"
	
)

# ╔═╡ 175b0dc1-3f13-4c5d-bdd1-d5530855180f
let
	(;Tamb,p) = data
	ng=ngas(data)
	D = MMatrix{ng,ng,Float64}(undef)
	A = MMatrix{ng,ng,Float64}(undef)
	#D_matrix!(D, Tamb, p, data)
	MultiComponentReactiveMixtureProject.D_A_matrices!(D, A, Tamb, p, data)
end

# ╔═╡ 728abf6f-3ff1-45c8-8b10-213f08b1b4dd
let
	(;iT, gni, gn) = data
	vis=GridVisualizer(resolution=(600,400), layout=(2,1), Plotter=PlutoVista)
	scalarplot!(vis[1,1], grid, sol[2,:], clear=false, label=gn[2])
	scalarplot!(vis[1,1], grid, sol[5,:], clear=false, color=:red, label=gn[5])
	scalarplot!(vis[1,1], grid, sol[6,:], clear=false, color=:blue, label=gn[6])
	
	scalarplot!(vis[2,1], grid, sol[iT,:], clear=false, label="T")
	reveal(vis)
end

# ╔═╡ 035d4123-7092-4429-8cfd-1e5926e84493
#sol, sys = setup_run_sim(grid, data);

# ╔═╡ 93970c02-91c6-499a-9318-f7f632604bb5
function bcondition(f,u,bnode,data)
	(;p,ip,iT,Tamb)=data
	
	boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_left,value=p)
	boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_right,value=p)

	boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_left,value=300)
	boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_right,value=400)	
    
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

	#times=[0,20]
	#sol=VoronoiFVM.solve(sys;inival=inival,times,verbose="an",log=true)
	sol=VoronoiFVM.solve(sys;inival=inival,verbose="an",log=true)
	return (sol,sys)
end

# ╔═╡ 65dbb492-4795-44ca-afcb-fb2a2c925d92
md"""
# References
1) Van‐Brunt, Alexander; Farrell, Patrick E.; Monroe, Charles W. (2022): Consolidated theory of fluid thermodiffusion. In: AIChE Journal 68 (5), Artikel e17599. DOI: 10.1002/aic.17599   .
1) Giovangigli, Vincent (2016): Solutions for Models of Chemically Reacting Mixtures. In: Yoshikazu Giga und Antonin Novotny (Hg.): Handbook of Mathematical Analysis in Mechanics of Viscous Fluids. Cham: Springer International Publishing, S. 1–52.
"""

# ╔═╡ Cell order:
# ╠═349e7220-dc69-11ee-13d2-8f95e6ee5c96
# ╠═0f102f06-3ff3-4bcc-8892-8d9190a87849
# ╠═f4286a46-ceb4-40b2-8ce5-7fcd231e3168
# ╠═bfebc943-28ef-477b-bc15-6cab9d398d92
# ╠═138b4d12-af0f-4a1c-a4aa-4f55e6038664
# ╠═63f0afeb-d0fc-4d58-8e24-7f8ccc587657
# ╟─9fde21eb-1112-40ec-8fe2-0c8a12c9926d
# ╟─9a7206ae-35f5-430f-b933-58c9494b9f0c
# ╟─2fcc681c-1b79-46cc-b96f-20143f8e331e
# ╟─3f9eea9d-7c2f-4d83-981a-a243fdf0531a
# ╟─94dd7674-751f-4128-b5eb-303fb9693c22
# ╠═e7ca4902-0e14-48ca-bcc6-96b06c85a39d
# ╠═7f1d9cf8-7785-48c1-853c-74680188121f
# ╠═175b0dc1-3f13-4c5d-bdd1-d5530855180f
# ╠═728abf6f-3ff1-45c8-8b10-213f08b1b4dd
# ╠═035d4123-7092-4429-8cfd-1e5926e84493
# ╠═93970c02-91c6-499a-9318-f7f632604bb5
# ╠═1e51701d-a893-4056-8336-a3772b85abe4
# ╟─65dbb492-4795-44ca-afcb-fb2a2c925d92
