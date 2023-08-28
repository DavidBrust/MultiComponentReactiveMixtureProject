### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	
	using VoronoiFVM
	using ExtendableGrids, GridVisualize
	using LinearAlgebra
	using LessUnitful
	using PlutoVista	
	using PlutoUI
	#using Colors

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="Stefan-Maxwell")

# ╔═╡ fa2002f1-1b35-4fe8-bbd1-f9fff62acbdf
function grid2(;L_CL=5.0*ufac"μm",L_DM=325.0*ufac"μm",nref=0)
	X_CL=range(0,L_CL,length=10*2^nref+1)
	X_DM=range(L_CL,L_DM,length=10*2^nref+1)
	grid=simplexgrid(glue(X_CL,X_DM))

	cellmask!(grid,[L_CL],[L_DM],2)
    #bfacemask!(grid,[L_CL], [L_CL],2)
	#bfacemask!(grid,[L_DM], [L_DM],3)
    grid
end

# ╔═╡ b659879f-3706-4937-b655-116f1be85069
begin
	const Γ_EL_CL=1
	const Γ_DM_GC=2
	#const Γ_CL_DM=3

	# bulk regions
	# CL = 1
	# DM = 2
end;

# ╔═╡ 7068922e-c797-416d-ad8e-197d80e3af1a
gridplot(grid2(L_CL=50.0*ufac"μm"),resolution=(600,200),legend=:rt,show=true)

# ╔═╡ 2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
md"""
# Governing system
"""

# ╔═╡ f4dcde90-6d8f-4b17-b4ec-367d2372637f
md"""
#### Species Mass balances
There are $N$ gas phase species in the system, species 'concentration' is expressed in mole fractions $x_{\alpha}, \alpha = 1 ... N$.
"""

# ╔═╡ 3703afb0-93c4-4664-affe-b723758fb56b
md"""
```math
\begin{align}
\partial_t x_{\alpha} + \partial_i J_{\alpha} &= 0
~,
\alpha = 1 ... N
\end{align}

```
"""

# ╔═╡ c4e08410-3d1a-4921-aa5c-2f232daaa07a
md"""
##### Force-flux relations
Express the Stefan-Maxwell diffusion equations in Matrix-Vector notation, notation taken from __R. Taylor and R. Krishna__, Multicomponent mass transfer.
"""

# ╔═╡ 313cd88a-1497-4d3a-b09a-9b98a6dad9c2
md"""
In a multi-component system of $N$ species, there are only $N-1$ independent equations. The system is completed by the constraint
```math
\sum_{\alpha} x_{\alpha} = 1.
```
The species flux of the last component $J_N$ can be eliminated via
```math
J_N = -\sum_{\alpha=1}^{N-1} J_{\alpha}.
```
"""

# ╔═╡ 8528e15f-cce7-44d7-ac17-432f92cc5f53
md"""
In the simplified case we look at ideal gas species at constant pressure. The driving force for transport of species $\alpha$ relative to the other species then reduces to $\nabla x_{\alpha}$. Using vector notation, the force-flux relations read
"""

# ╔═╡ 19a46c60-d7b5-45da-a3c5-3687316b2a2c
md"""
```math
\begin{align}
c_{\text{t}} \nabla \mathbf{x} &= - \mathbf{M} ~\mathbf{J} \\
c_{\text{t}} \nabla x_{\alpha} &= - M_{\alpha \alpha} J_{\alpha} - \sum_{\beta=1 \atop \beta \neq \alpha}^{N-1} M_{\alpha \beta} J_{\beta}
\end{align}
```
"""

# ╔═╡ d5987ad8-e71b-4d40-b1c1-88142558265d
md"""
In the follwing, equation numbers refer to R. Taylor and R. Krishna, Multicomponent mass transfer.
The coefficients $M_{\alpha \alpha}$ and $M_{\alpha \beta}$ are defined by (2.1.21, 2.1.22):

"""

# ╔═╡ d2adc467-d366-4461-9a00-9832e7d5a737
md"""
```math
	M_{\alpha \alpha} = \frac{x_{\alpha}}{D_{\alpha N}} + \sum_{\gamma=1 \atop \gamma \neq \alpha}^{N} \frac{x_{\gamma}}{D_{\alpha \gamma}}
```
```math
M_{\alpha \beta} = -x_{\alpha} \left ( \frac{1}{D_{\alpha \beta}} - \frac{1}{D_{\alpha N}} \right)
```
"""

# ╔═╡ a6afe118-dcbd-4126-8646-c7268acfacf3
md"""
The numerical fluxes $\textbf{J}^{\textrm{num}}_{kl}$ are computed using the backslash operator (Julia solver for linear systems) 

$\textbf{J}^{\textrm{num}}_{kl} 
= 
M^{\textrm{num}}_{kl} 
\ \backslash_\textrm{Julia}\ 
\textbf{F}^{\textrm{num}}_{kl}$

"""

# ╔═╡ f3dc68cb-f9c3-4e85-892f-34358831e3eb
function Dmatrix(n, D)
	v = zeros(eltype(D),n,n)
	k=1
	for i=1:(n-1)
		for j=(i+1):n
			v[j,i] = D[k]
			k +=1
		end
	end
	Symmetric(v, :L)
end

# ╔═╡ b6381008-0280-404c-a86c-9c9c3c9f82eb
function M_matrix(data,x)
	n=data.ng
	#D=Dmatrix(n, data.Dg)
	D=data.D
	M=zeros(eltype(x), n-1, n-1)
	for i=1:(n-1)
		M[i,i] = x[i]/D[i,n]
		for j=1:(n-1)
			if j != i
				M[i,j] = - x[i]*(1/D[i,j]-1/D[i,n])
			end
		end
		for k=1:n
			if k != i
				M[i,i] += x[k]/D[i,k]
			end
		end
	end
	-1.0*M./data.ct
end

# ╔═╡ 02b76cda-ffae-4243-ab40-8d0fe1325776
md"""
##### Auxiliary functions
"""

# ╔═╡ dd7e659d-2fbf-494e-973e-03d913c60ca0
density_means(u,n) = [0.5*(u[i,1]+u[i,2]) for i in 1:n]

# ╔═╡ ed7941c4-0485-4d84-ad5b-383eb5cae70a
function flux(f,u,edge,data)
	n=data.ng
	F=zeros(eltype(u), n-1)
	for i=1:(n-1)
		F[i] = u[i,1]-u[i,2]
	end
	
	#F3 = u[3,1]-u[3,2]
	#F1 = -1*sedan(0.0, u[1,1], u[1,2])
    #F2 = -1*sedan(0.0, u[2,1], u[2,2])
	#F3 = -1*sedan(0.0, u[3,1], u[3,2])
	
	#x1e, x2e, x3e = density_means(u)
    
	# computation of fluxes J 
    

	J = M_matrix(data,density_means(u,n)) \ F

	
	# 
    #f[1] = J1 
    #f[2] = J2
	f[1:(n-1)] = J
	
	#f[3] = J3
	#f[3] via reaction: ∑xi = 1

end

# ╔═╡ c2cd05b5-dafb-4959-80b0-d5edb5f421aa
function sedan(g, u1, u2)
	# flux: j = - ( u' + u*q)
	# dirichlet BCs u1        u2 
	# line segment  |---------|
	# colloc. pts  x1        x2
	# h = |x2-x1|
	# g = q*h
  	bp,bm=fbernoulli_pm(g) # Bernoulli function for B(+g) and B(-g)
  	return (-1.0)*( u2*bm - u1*bp )
end

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)	
	boundary_dirichlet!(f,u,bnode,1,Γ_DM_GC,0.1) # CO
	boundary_dirichlet!(f,u,bnode,2,Γ_DM_GC,0.1) # H2
	boundary_dirichlet!(f,u,bnode,3,Γ_DM_GC,0.8) # CO2
	#boundary_dirichlet!(f,u,bnode,4,Γ_DM_GC,0.45) # N2
	
	#boundary_dirichlet!(f,u,bnode,1,Γ_EL_CL,0.4) # CO
	
end

# ╔═╡ f9af078f-6262-4bee-9b1d-04cdf2e9adec
grid=grid2(L_CL=50*ufac"μm",nref=0)

# ╔═╡ aa5facd4-ffac-4b24-8868-d94b1e7b615d
md"""

### Charge transfer reaction
In the simplified case only consider a single charge transfer reachtion, that both creates H2 and CO and consumes CO2:
```math
CO_2 \rightarrow CO + H2
```

This is preliminary.
"""

# ╔═╡ 8c639001-b3ef-494c-a3eb-c600b011c159
md"""
### Darcy's law
Calculate gas velocity in porous medium. This gas velocity is responsible for convective transport.
"""

# ╔═╡ bba07772-a923-482e-b24b-acfaa24ac81c
md"""
```math
	\vec u_{\text{g}} = - \frac{\kappa^{\text{eff}}}{\mu_{\text{g}}} \nabla p
```
"""

# ╔═╡ c13b04ae-854e-4d26-93b4-101c9f69d784
function κ_μ_darcy(data)
	ϵ_sat=data.ϵ0*(1-data.S)
	κ_eff=data.κ*ϵ_sat^3/(1-ϵ_sat)^2*(1-data.S)^3
	κ_eff/data.μ_g
end

# ╔═╡ 59bb0503-52d8-4a76-a3f3-ce0fdcfcb170
@Base.kwdef mutable struct ModelData
# diffusion medium thickness
L_DM::Float64 			= 325.0*ufac"μm"	
# catalyst layer thickness
L_CL::Float64 			= 5.0*ufac"μm"
	
# number of gas phase species
ng::Int64		 		= 4
# 1:CO, 2:H2, 3:CO
gn::Dict{Int, String} 	= Dict(1=>"CO", 2=>"H2", 3=>"CO2", 4=>"N2")
# stochiometric coeff for CT reaction
nug::Vector{Int} 		= [1,1,-2,0]
# number of transferred e- per CT reaction turnover
nCT::Vector{Int} 		= [2]
	
# gas density
ρg::Float64 			= 1.0*ufac"kg/(m^3)"
# gas dynamic viscosity (use air at 25 °C)
μ_g::Float64 			= 18.37e-6*ufac"Pa*s"
# total molar concentration
ct::Float64 			= 1.0*ufac"mol/(dm^3)"

# convective gas transport via Darcy's law
# saturation, fraction of flooded pores
S::Float64 				= 0.0
# dry porosity 
ϵ0::Float64 			= 0.8
# gas permeability in porous medium
κ::Float64 				= 1.34e-12*ufac"m^2"
# assumed pressure drop by gas flow through porous media
Δp::Float64 			= 0.01*ufac"bar"
	
# gas velocity, convective transport, different for the 2 regions
vg::Vector{Float64}		= 0*[1.7e-5,0.35]*ufac"m/s"
#vg::Vector{Float64}		= -1e-5*[1.7e-5,0.35]*ufac"m/s"

# gas phase diffusivity coefficients
Dg::Vector{Float64} 	= [1.0,10.0,3.0,8.0,2.0,1.8]*1e-5*ufac"cm^2/s"
D::Symmetric{Float64, Matrix{Float64}} = Dmatrix(ng, Dg)

# gas phase species molar mass
Mg::Vector{Float64} 	= fill(10.0*ufac"g/mol",ng)

# number of aqueous phase species
ni::Int64 				= 2

# ion charge numbers
z::Vector{Int}   		= [-1,1]

# reference voltage
E_ref::Float64   		= 0.0*ufac"V"

# solid phase electrical conductivity
σs::Float64 	= 200.0*ufac"S/m"

# volume specific active surface area
#av::Float64 	= 1.0*ufac"1/m"
av::Float64 	= 1.5e7*ufac"1/m"

e::Float64 = ph"e"
F::Float64 = ph"e"*ph"N_A"

end

# ╔═╡ cc94646d-b51c-4e1e-ba69-3bba8f6a771b
let
	data=ModelData()
	typeof(Dmatrix(data.ng, data.Dg))
end

# ╔═╡ 3f643eb7-ff5f-4e88-8a57-5c559687ad1b
begin
	const ngas=ModelData().ng
	const nion=ModelData().ni
	
	# variable indices	
	const igas=1
end;

# ╔═╡ 71e05017-adb5-430f-8477-10b62eff8ba7
function R_CT(u, data::ModelData)
 
	# current density of CT reaction, mA cm-2, reaction rate limited by availability of reactant CO2
	rr=50*ufac"mA/cm^2"*u[3]	
	
	# species molar flux from CT reaction, eq. (26)
	@. data.nug*data.av*rr/data.nCT/data.F	
end

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	# volumetric charge-transfer reaction in CL region
	if node.region == 1
		rr=R_CT(u,data)
		for i=1:(ngas-1)
			f[i]=rr[i]
		end
	end
	
	# ∑xi = 1
	f[ngas]=log(sum(u[1:ngas]))
	
end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
begin
	sys=VoronoiFVM.System( 	grid;
							data=ModelData(),
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)
	enable_species!(sys; species=collect(igas:(igas+ngas-1)), regions=[1,2])

	inival=unknowns(sys)
	#inival[:,:].=1
	inival[igas:(igas+ngas-1),:].=1/ngas
	
	sol=solve(inival,sys)
end;

# ╔═╡ 4b41d985-8ebc-4cab-a089-756fce0d3060
let
	# visualization
	vis=GridVisualizer(resolution=(600,300))
	
	c1 = colorant"blue"
	c2 = colorant"red"
	cols=range(c1, stop=c2, length=ngas)

	for i in 1:ngas
		scalarplot!(vis, grid, sol[i,:], color=cols[i], label="$(ModelData().gn[i])", clear=false, legend=:ct)
	end

	reveal(vis)

end

# ╔═╡ 60a5bdda-7ebf-4585-b388-cba42f23347b
let
	data=ModelData()
	κ_μ_darcy(data) * data.Δp/(data.L_DM+data.L_CL) # m/s
	#darcy_perm(data)
end

# ╔═╡ Cell order:
# ╠═11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╟─863c9da7-ef45-49ad-80d0-3594eca4a189
# ╠═fa2002f1-1b35-4fe8-bbd1-f9fff62acbdf
# ╠═b659879f-3706-4937-b655-116f1be85069
# ╠═7068922e-c797-416d-ad8e-197d80e3af1a
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─f4dcde90-6d8f-4b17-b4ec-367d2372637f
# ╟─3703afb0-93c4-4664-affe-b723758fb56b
# ╟─c4e08410-3d1a-4921-aa5c-2f232daaa07a
# ╟─313cd88a-1497-4d3a-b09a-9b98a6dad9c2
# ╟─8528e15f-cce7-44d7-ac17-432f92cc5f53
# ╠═19a46c60-d7b5-45da-a3c5-3687316b2a2c
# ╟─d5987ad8-e71b-4d40-b1c1-88142558265d
# ╟─d2adc467-d366-4461-9a00-9832e7d5a737
# ╠═ed7941c4-0485-4d84-ad5b-383eb5cae70a
# ╟─a6afe118-dcbd-4126-8646-c7268acfacf3
# ╠═f3dc68cb-f9c3-4e85-892f-34358831e3eb
# ╠═cc94646d-b51c-4e1e-ba69-3bba8f6a771b
# ╠═b6381008-0280-404c-a86c-9c9c3c9f82eb
# ╟─02b76cda-ffae-4243-ab40-8d0fe1325776
# ╠═dd7e659d-2fbf-494e-973e-03d913c60ca0
# ╠═c2cd05b5-dafb-4959-80b0-d5edb5f421aa
# ╠═78cf4646-c373-4688-b1ac-92ed5f922e3c
# ╠═29d66705-3d9f-40b1-866d-dd3392a1a268
# ╠═3f643eb7-ff5f-4e88-8a57-5c559687ad1b
# ╠═333b5c80-259d-47aa-a441-ee7894d6c407
# ╠═4b41d985-8ebc-4cab-a089-756fce0d3060
# ╠═f9af078f-6262-4bee-9b1d-04cdf2e9adec
# ╟─aa5facd4-ffac-4b24-8868-d94b1e7b615d
# ╠═71e05017-adb5-430f-8477-10b62eff8ba7
# ╟─8c639001-b3ef-494c-a3eb-c600b011c159
# ╟─bba07772-a923-482e-b24b-acfaa24ac81c
# ╠═c13b04ae-854e-4d26-93b4-101c9f69d784
# ╠═60a5bdda-7ebf-4585-b388-cba42f23347b
# ╠═59bb0503-52d8-4a76-a3f3-ce0fdcfcb170
