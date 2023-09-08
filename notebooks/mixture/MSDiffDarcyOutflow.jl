### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,"..\\.."))
	using Revise
	using VoronoiFVM, VoronoiFVM.SolverStrategies
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve
	using LinearAlgebra
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots
	using PlutoUI, HypertextLiteral
	using Colors

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="M-S Transport + Darcy")

# ╔═╡ 57117131-8c83-4788-9f8e-f0b1b9404f0a
function grid1D()
	X=0:0.01:1
	grid=simplexgrid(X)

	# catalyst region
	cellmask!(grid,[.46],[0.54],2)
	#cellmask!(grid,[0.0],[0.2],2)
	grid
end

# ╔═╡ 173364be-0b52-445c-817d-f90b2a77e6e9
function grid2D()
	X=0:0.05:1
	Y=0:0.2:1
	grid=simplexgrid(X,Y)

	# catalyst region
	cellmask!(grid,[0.7,0.0],[0.9,1.0],2)
	bfacemask!(grid, [1,0.2],[1,0.8],5)
	
	grid
end

# ╔═╡ 87630484-de43-46bc-a179-3e00b6e63c2a
gridplot(grid1D())
#gridplot(grid2D())

# ╔═╡ 0a911687-aff4-4c77-8def-084293329f35
begin
	# for 1D domain
	const Γ_left = 1
	const Γ_right = 2

	# for 2D domain
	#const Γ_bottom = 1
	#const Γ_right = 2
	#const Γ_right = 5
	#const Γ_top = 3
	#const Γ_left = 4
	
end;

# ╔═╡ 2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
md"""
# Governing system
"""

# ╔═╡ 9a22118f-cbf0-4412-96dc-66e43a14b0cd
md"""
Following __Weber, A. Z., et al. (2014) Journal of The Electrochemical Society 161(12): F1254-F1299.__ and __Weber, A. Z. and J. Newman (2005) International Communications in Heat and Mass Transfer 32(7): 855-860.__ to formulate a system of equations for multi-component species transport through porous media. This system is proposed as an alternative to the Dusty-Gas-Model but explicitly includes convective transport through D'arcys law.

There are $\nu$ gas phase species in the system, gas phase composition is expressed in partial pressures $p_i$ and in molar fractions for an ideal gas mixture $x_i=p_i/p$,  $i = 1 ... \nu$.
"""

# ╔═╡ 83e4a00a-2834-4a68-87ee-81394adea30f
md"""
```math
\begin{align}
	\frac{\partial \rho}{\partial t} - \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\frac{1}{RT} \frac{\partial p_i}{\partial t} - \nabla \cdot \vec N_i + R_i &= 0
~,
i = 1 ... \nu \\
	\sum_{j=1 \atop j \neq i}^{\nu} \frac{x_j \vec N_i-x_i \vec N_j}{D_{ij}^{\text{eff}}} + \frac{\vec N_i}{D_{i,\text K}^{\text{eff}}} &= -\frac{\nabla p_i}{RT}  \\
	\sum_{i=1}^\nu p_i &= p
\end{align}
```

where $\rho$ is the mixture density, $\vec v$ is the (mass-averaged ?) superficial mixture velocity $\vec N_i$ is the molar flux ($\frac{\text{mol}}{\text{m}^2 \text{s}}$) and $R_i$ is the molar volumetric source/sink ($\frac{\text{mol}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
"""

# ╔═╡ c0f65543-a64c-486a-a936-7dc55d75d5f1
md"""
## Implementation
"""

# ╔═╡ 9b2e8232-c4c2-4c69-acf5-49ce9f6e9270
md"""
Implement the equation system described above in a way that inherently conserves total mass by replacing the species mass flux of the last species of index $N$ via the sum of the remaining $N-1$ species mass fluxes and their difference to total mass flux.
"""

# ╔═╡ b28f3573-751b-4154-8249-d140baa65d0b
md"""
```math
	\vec N_N=\frac{\vec m - \sum_i^{N-1}M_i \vec N_i}{M_N}
```
Here $\vec N_i$ are molar species fluxes, $\vec m$ is total mass flux and $M_i$ are species molar masses.

Substituting above expression for $\vec N_N$ into the M-S equations with Knudsen addition and writing out the resulting linear equation for $N=4$ for the first species $i=1$:
```math
\begin{align}
	-\frac{\nabla p_1}{RT} &= \frac{\vec N_1}{D_1^K} + \frac{x_2 \vec N_1 - x_1 \vec N_2}{D_{12}^B} + \frac{x_3 \vec N_1 - x_1 \vec N_3}{D_{13}^B} \\
	&+ \frac{x_4 \vec N_1 - x_1 \left( \vec m - M_1 \vec N_1 - M_2 \vec N_2 - M_3 \vec N_3 \right)/M_4} {D_{14}^B}
\end{align}
```
"""

# ╔═╡ 0f455649-cd10-41b4-896d-29b85cdbdb89
md"""
To represent the system of equations in matrix form, we can derive the (N-1) x (N-1) matrix ``M=(m_{ij})`` with
```math
\begin{align}
	m_{ii} &= \frac{1}{D^K_i} + \frac{x_i M_i}{M_N D^B_{iN}} + \sum_{j\neq i}^N \frac{x_j}{D^B_{ij}} \quad (i=1\dots N-1) \\ 
	m_{ij} &= -x_i \sum_{j\neq i \atop j \neq N } \left( \frac{1}{D^B_{ij}} - \frac{M_j}{M_N D^B_{iN}} \right)
\end{align}
```
"""

# ╔═╡ e9d0be5a-acbe-44fc-88d8-e1b40c0a46bd
md"""
and right hand side vector ``f = (f_i)``
```math
	f_i= -  \frac{\nabla p_i}{RT} + \frac{x_i \vec m}{M_N D^B_{iN}}
```
"""

# ╔═╡ 45b34b7b-f816-41ca-9f47-bd05f862db59
md"""
such that
"""

# ╔═╡ 8490e4f0-04f9-4cfa-a547-0cb4e4fffc59
md"""
```math
	M\begin{pmatrix}
\vec N_1\\
\vdots\\
\vec N_{N-1}
\end{pmatrix}
=
f
```
"""

# ╔═╡ e8b49b2e-d19e-4d46-a399-9972919cd680
md"""
This function describes the outflow boundary condition.
It is called on edges (including interior ones) which have at least one node
on one of the outflow boundaries. Within this function
`outflownode` can be used to identify
that node.  There is some ambiguity in the case that both nodes are outflow
nodes, in that case it is assumed that the contribution is zero. In the present case this is guaranteed by the constant Dirichlet boundary condition for the pressure.
"""

# ╔═╡ 02b76cda-ffae-4243-ab40-8d0fe1325776
md"""
# Auxiliary functions
"""

# ╔═╡ 44aa5b49-d595-4982-bbc8-100d2f199415
md"""
# System Setup and Solution
"""

# ╔═╡ e25e7b7b-47b3-457c-995b-b2ee4a87710a
md"""
## Model Data
"""

# ╔═╡ 3a35ac76-e1b7-458d-90b7-d59ba4f43367
begin	
	Base.@kwdef mutable struct ModelData{NG} <:AbstractModelData
	
	# catalyst / chemistry data
	kinpar::FixedBed.KinData{nreac(XuFroment)} = XuFroment
	#kinpar::FixedBed.KinData{nreac(ReOrderSpec)} = ReOrderSpec

	# kinetic constants for simple reactions
	kp::Float64=1.0
	km::Float64=0.1
	isreactive::Bool=true
		
	ip::Int64=NG+1 # index of total pressure variable
	
	gn::Dict{Int, Symbol} 	= kinpar.gn # names and fluid indices
	gni::Dict{Symbol, Int}  = kinpar.gni # inverse names and fluid indices
	Fluids::Vector{FluidProps} = kinpar.Fluids # fluids and respective properties in system
	
	X0::Vector{Float64} = let
		x=zeros(Float64, NG)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		#x[gni[:CO]] = 1.0
		#x[gni[:CH4]] = 1.0
		#x[gni[:H2O]] = 1.0
		#x[gni[:N2]] = 1.0		
		x/sum(x)
	end # inlet composition	
	
	## porous filter data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0
	
	# for cylindrical geometry
	#Ac::Float64=pi*wi^2/4*ufac"m^2" # cross-sectional area
	# 1D grid
	Ac::Float64=1.0*ufac"m^2" # cross-sectional area
	
	#ϕ::Float64=0.36 # porosity, exp determined
	ϕ::Float64=0.33 # porosity, VitraPor sintetered filter class 0
	
	# approximation from Wesselingh, J. A., & Krishna, R. (2006). Mass Transfer in Multicomponent Mixtures
	γ_τ::Float64=ϕ^1.5 # constriction/tourtuosity factor

	#k::Float64=1.23e-10*ufac"m^2" # permeability , por class 0
	k::Float64=1.23e-12*ufac"m^2" # permeability , por class 0
	#k::Float64=1.23e-13*ufac"m^2" # permeability , por class 0


	## Flow data
	#norm conditions
	pn::Float64 = 1.0*ufac"bar"
	Tn::Float64 = 273.15*ufac"K"


	#Qflow::Float64=0.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	#Qflow::Float64=148.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	#Qflow::Float64=1480.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	#Qflow::Float64=14800.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	Qflow::Float64=50000.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		

	MWin::Float64 = molarweight_mix(Fluids, X0)
	mdotin::Float64=MWin*Qflow*pn/(ph"R"*Tn)*ufac"kg/s"
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	p::Float64=1.0*ufac"atm" # reactor pressure

	u0::Float64=Qflow/Ac*ufac"m/s" # mean inlet superficial velocity
	nfluxin::Float64=Qflow*pn/(ph"R"*Tn)/Ac
	mfluxin::Float64=mdotin/Ac*ufac"kg/(m^2*s)" # mean inlet mass flux
	#ubot::Float64=u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate
	
end;
	
	# !!!ALLOC Method to be called instead of data.ng
	FixedBed.ngas(::ModelData{NG}) where NG = NG
	
	# !!!ALLOC Additional constructo taking ng as parameter	
	ModelData(;ng=S3P.ng, kwargs...) = ModelData{ng}(;kwargs...)

end;

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	(;ip,kp,km,gni,isreactive)=data

	if node.region == 2 && isreactive # catalyst layer

		# 4 H2 + CO2 <-> 2 H2O + CH4
		r = kp *(u[gni[:H2]]/ufac"bar")^4 *u[gni[:CO2]]/ufac"bar" -km *(u[gni[:H2O]]/ufac"bar")^2	*u[gni[:CH4]]/ufac"bar"		
		r *=embedparam(node)
		
		f[gni[:H2]]=4*r
		f[gni[:CO2]]=r
		f[gni[:H2O]]=-2*r
		f[gni[:CH4]]=-r
	end
		
	ng=ngas(data)
	# last species i=ng via ∑pi = p
	for i=1:ng
		f[ng] += u[i]
	end
	#f[ng] = log(f[ng]/u[ip])
	f[ng] = f[ng]-u[ip]
end

# ╔═╡ a8d57d4c-f3a8-42b6-9681-29a1f0724f15
function darcyvelo(u,data)
	(;ip, k, Tamb) = data

	ng=ngas(data)
	μ = 2.0e-5*ufac"Pa*s"
	k/μ*(u[ip,1]-u[ip,2])	
end

# ╔═╡ d4bc847b-052f-4d40-9211-12dbe7e06ee1
function boutflow(f,u,edge,data)
	(;Tamb)=data
	ng=ngas(data)
		
	for i=1:ng
		# specify flux at boundary
		k=outflownode(edge)		
		f[i] = -darcyvelo(u,data)*u[i,k]/(ph"R"*Tamb)		
	end	
end

# ╔═╡ 7da59e27-62b9-4b89-b315-d88a4fd34f56
function inlet(f,u,bnode,data)
	# left boundary 1D
	if bnode.region==Γ_left 
		(;ip,X0,u0,pn,Tn,mfluxin,nfluxin)=data
		ng=ngas(data)
		
				
		# flow velocity is normal to top boundary
		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative

		for i=1:(ng-1)
			# specify flux at inlet
			f[i] = -X0[i]*u0*pn/(ph"R"*Tn)*embedparam(bnode)		
		end

		# specify inlet mass fluxs
		f[ip] = -mfluxin*embedparam(bnode)
	end
end

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
	(;ip,p,pn,Tn,u0,X0,Tamb,mfluxin,nfluxin)=data
	ng=ngas(data)

	# specifying (N-1) species molar fluxes and total mass flux does not converge...
	#inlet(f,u,bnode,data)

	boundary_neumann!(f,u,bnode, species=ip, region=Γ_left,value=mfluxin*embedparam(bnode))

	# specifying the inlet partial pressures converges
	for i=1:ng	
		boundary_dirichlet!(f,u,bnode, species=i,region=Γ_left,value=u[ip]*X0[i])
	end
	
	boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_right,value=p)

end

# ╔═╡ 2191bece-e186-4d8e-8a21-3830441baf11
#function D_matrix!(data, D, T, p)
function D_matrix!(data, D)
	(;p,Tamb)=data
	ng=ngas(data)
	for i=1:(ng-1)
		for j=(i+1):ng
			#Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
			Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], Tamb, p)
			D[j,i] = Dji
			D[i,j] = Dji
		end
	end
	nothing
end

# ╔═╡ 22447e9e-df5f-4773-a83e-8645ad0c7873
function rho_mix_flux(u,data)
	ng=ngas(data)

	(;Tamb,Fluids)=data
	rho=zero(eltype(u))
	for i=1:ng
		rho += Fluids[i].MW*0.5*(u[i,1]+u[i,2])
	end
	rho /= ph"R"*Tamb	
end

# ╔═╡ ed7941c4-0485-4d84-ad5b-383eb5cae70a
function flux(f,u,edge,data)
	(;ip,gni,X0,Tamb,k,γ_τ) = data

	ng=ngas(data)

	F = MVector{ng-1,eltype(u)}(undef)

	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)
	D = MMatrix{ng,ng,eltype(u)}(undef)

	# compute total mass flux
	@inline f[ip] =rho_mix_flux(u,data)*darcyvelo(u,data)
		
	@inline D_matrix!(data, D)

	xi=zero(eltype(u))
	xj=zero(eltype(u))

	
	# assemble matrix M and right hand side vector f
	for i=1:(ng-1)
		xi= (u[i,1]+u[i,2])/(u[ip,1]+u[ip,2])
		# Matrix M
		M[i,i] =  1.0/DK_eff(data,Tamb,i)		
	
		M[i,i] += xi*data.Fluids[i].MW/data.Fluids[ng].MW/D[i,ng]
	
		for j=1:ng
			if j != i
				xj=(u[j,1]+u[j,2])/(u[ip,1]+u[ip,2])
				M[i,i] += xj/D[i,j]

				if j != ng
					M[i,j] = -xi/D[i,j]
					M[i,j] += xi*data.Fluids[j].MW/data.Fluids[ng].MW/D[i,ng]
				end
			end
		end
		# right hand side vector f
		F[i] = (u[i,1]-u[i,2])/(ph"R"*Tamb) + xi*f[ip]/data.Fluids[ng].MW/D[i,ng]
		
	end
			
	# computation of fluxes N
	@inline inplace_linsolve!(M,F)

	for i=1:(ng-1)
		f[i] = F[i]		
	end	
end

# ╔═╡ dd5e0a30-80ef-4eac-a014-b5a0f6a3c5fe
function rho_mix_storage(u,data)
	ng=ngas(data)

	(;Tamb,Fluids)=data
	rho=zero(eltype(u))
	for i=1:ng
		rho += Fluids[i].MW*u[i]
	end
	rho /= ph"R"*Tamb	
end

# ╔═╡ 6779db7a-3823-4a41-8b2d-5558dcd73943
function storage(f,u,node,data)
	(;Tamb)=data
	ng=ngas(data)

	for i=1:ng
		f[i]=u[i]/(ph"R"*Tamb)
	end
	
	# total pressure
	f[ng+1] = rho_mix_storage(u,data)
end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
function main(;data=ModelData(),nref=0,verbose="aen")

	grid=grid1D()
	#grid=grid2D()
	
	(;ip,p,X0,)=data
	ng=ngas(data)
	
	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=flux,
							reaction=reaction,
							storage=storage,
							bcondition=bcond,					
							boutflow=boutflow,
							outflowboundaries=[Γ_right]
							)
	
	enable_species!(sys; species=collect(1:(ng+1))) # gas phase species pi & ptotal	
	
	inival=unknowns(sys)

	inival[:,:].=1.0*p	
	for i=1:ng
		inival[i,:] .*= X0[i]
	end
	
	# solve problem with parameter embedding (homotopy)
	control = SolverControl(nothing, sys;)
		control.Δp=1.0e-1
		#control.Δp_min=1.0e-6
		control.Δp_grow=2.0
		control.handle_exceptions=true
		control.Δu_opt=1.0e5

	embed=[0.0,1.0]
	sol=solve(sys;inival=inival,control,embed,verbose=verbose)

	sol,grid,sys,data
end;

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
# ╠═╡ skip_as_script = true
#=╠═╡
sol_,grid,sys,data_embed=main(;data=ModelData());
  ╠═╡ =#

# ╔═╡ 2790b550-3105-4fc0-9070-d142c19678db
md"""
## (Partial) Pressure Profiles
"""

# ╔═╡ 36d02d06-6d0d-4b98-8d60-f2df0afaeaad
#=╠═╡
@bind t Slider(sol_.t,show_value=true)
  ╠═╡ =#

# ╔═╡ 55f305b8-47a2-4fdf-b1cb-39f95f3dfa36
#=╠═╡
sol=sol_(t)
  ╠═╡ =#

# ╔═╡ 8fcf636c-2330-4eda-9bb3-298e6a53dfd1
#=╠═╡
let
	(;gni,ip,p)=data_embed
	vis = GridVisualizer(legend = :rt, ylabel="Pressure / Pa")
	scalarplot!(vis, grid, sol[ip,:], label="total")
	#scalarplot!(vis, grid, sol[gni[:N2],:], label="N2")
	#scalarplot!(vis, grid, sol[gni[:H2],:], label="H2")	
	#scalarplot!(vis, grid, sol[gni[:CO2],:], label="CO2")
	#scalarplot!(vis, grid, sol[gni[:CH4],:], label="CH4")
	#scalarplot!(vis, grid, sol[gni[:CO],:], label="CO")
	#scalarplot!(vis, grid, sol[gni[:H2O],:], label="H2O")
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ 9e02a838-2c26-453e-97d2-1de46b4d66ea
md"""
# Flow Variation
"""

# ╔═╡ bfa85ada-9ccb-48e7-b1bf-f3c6562a2cb3
let
	
	vis = GridVisualizer(legend = :rt, ylabel="Pressure / Pa")
	cs=colormap("Blues", 6)

	sol_,grid,_,data_embed=main(;data=ModelData(Qflow=0.0*ufac"l/minute"),verbose=false);
	(;gni)=data_embed

	idx=gni[:CO2]
	
	sol = sol_(sol_.t[end])
	scalarplot!(vis, grid, sol[idx,:], label="CH4 0 L/min", color=cs[2])
	
	sol_,_,_,_=main(;data=ModelData(Qflow=0.15*ufac"l/minute"),verbose=false);
	sol = sol_(sol_.t[end])
	scalarplot!(vis, grid, sol[idx,:], clear=false, label="CH4 0.15 L/min", color=cs[3])

	sol_,_,_,_=main(;data=ModelData(Qflow=1.5*ufac"l/minute"),verbose=false);
	sol = sol_(sol_.t[end])
	scalarplot!(vis, grid, sol[idx,:], clear=false, label="CH4 1.5 L/min", 
	color=cs[4])

	sol_,_,_,_=main(;data=ModelData(Qflow=15*ufac"l/minute"),verbose=false);
	sol = sol_(sol_.t[end])
	scalarplot!(vis, grid, sol[idx,:], clear=false, label="CH4 15 L/min", 
	color=cs[5])

	sol_,_,_,_=main(;data=ModelData(Qflow=50*ufac"l/minute"),verbose=false);
	sol = sol_(sol_.t[end])
	scalarplot!(vis, grid, sol[idx,:], clear=false, label="CH4 50 L/min", 
	color=cs[6])
	
	reveal(vis)

end

# ╔═╡ 8e0fdbc7-8e2b-4ae9-8a31-6f38baf36ef3
md"""
## Flows Over Boundaries
"""

# ╔═╡ 3c4b22b9-4f55-45b9-98d8-6cd929c737c2
md"""
1. CO
1. H2
1. CH4
1. H2O
1. CO2
1. N2
"""

# ╔═╡ 2a4c8d15-168f-4908-b24b-8b65ec3ea494
#=╠═╡
let
	(;pn,Tn,Qflow,X0,mfluxin,nfluxin,Fluids)=data_embed
	nfluxin.* X0, mfluxin
end
  ╠═╡ =#

# ╔═╡ e1602429-74fa-4949-bb0f-ecd681f52e42
function checkinout(sys,sol)
	
	tfact=TestFunctionFactory(sys)
	tf_in=testfunction(tfact,[Γ_right],[Γ_left])
	tf_out=testfunction(tfact,[Γ_left],[Γ_right])
	(;in=integrate(sys,tf_in,sol),out=integrate(sys,tf_out,sol) )
end

# ╔═╡ eb44075f-fbd3-4717-a440-41cb4eda8de1
#=╠═╡
checkinout(sys,sol)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╠═863c9da7-ef45-49ad-80d0-3594eca4a189
# ╠═57117131-8c83-4788-9f8e-f0b1b9404f0a
# ╠═173364be-0b52-445c-817d-f90b2a77e6e9
# ╠═87630484-de43-46bc-a179-3e00b6e63c2a
# ╠═0a911687-aff4-4c77-8def-084293329f35
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─9a22118f-cbf0-4412-96dc-66e43a14b0cd
# ╟─83e4a00a-2834-4a68-87ee-81394adea30f
# ╟─c0f65543-a64c-486a-a936-7dc55d75d5f1
# ╟─9b2e8232-c4c2-4c69-acf5-49ce9f6e9270
# ╟─b28f3573-751b-4154-8249-d140baa65d0b
# ╟─0f455649-cd10-41b4-896d-29b85cdbdb89
# ╟─e9d0be5a-acbe-44fc-88d8-e1b40c0a46bd
# ╟─45b34b7b-f816-41ca-9f47-bd05f862db59
# ╟─8490e4f0-04f9-4cfa-a547-0cb4e4fffc59
# ╠═ed7941c4-0485-4d84-ad5b-383eb5cae70a
# ╠═78cf4646-c373-4688-b1ac-92ed5f922e3c
# ╠═6779db7a-3823-4a41-8b2d-5558dcd73943
# ╠═a8d57d4c-f3a8-42b6-9681-29a1f0724f15
# ╟─e8b49b2e-d19e-4d46-a399-9972919cd680
# ╠═d4bc847b-052f-4d40-9211-12dbe7e06ee1
# ╠═7da59e27-62b9-4b89-b315-d88a4fd34f56
# ╠═29d66705-3d9f-40b1-866d-dd3392a1a268
# ╟─02b76cda-ffae-4243-ab40-8d0fe1325776
# ╠═2191bece-e186-4d8e-8a21-3830441baf11
# ╠═22447e9e-df5f-4773-a83e-8645ad0c7873
# ╠═dd5e0a30-80ef-4eac-a014-b5a0f6a3c5fe
# ╟─44aa5b49-d595-4982-bbc8-100d2f199415
# ╠═333b5c80-259d-47aa-a441-ee7894d6c407
# ╠═aa498412-e970-45f2-8b11-249cc5c2b18d
# ╟─e25e7b7b-47b3-457c-995b-b2ee4a87710a
# ╠═3a35ac76-e1b7-458d-90b7-d59ba4f43367
# ╟─2790b550-3105-4fc0-9070-d142c19678db
# ╠═55f305b8-47a2-4fdf-b1cb-39f95f3dfa36
# ╠═8fcf636c-2330-4eda-9bb3-298e6a53dfd1
# ╠═36d02d06-6d0d-4b98-8d60-f2df0afaeaad
# ╟─9e02a838-2c26-453e-97d2-1de46b4d66ea
# ╠═bfa85ada-9ccb-48e7-b1bf-f3c6562a2cb3
# ╟─8e0fdbc7-8e2b-4ae9-8a31-6f38baf36ef3
# ╟─3c4b22b9-4f55-45b9-98d8-6cd929c737c2
# ╠═eb44075f-fbd3-4717-a440-41cb4eda8de1
# ╠═2a4c8d15-168f-4908-b24b-8b65ec3ea494
# ╠═e1602429-74fa-4949-bb0f-ecd681f52e42
