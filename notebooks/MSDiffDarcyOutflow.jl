### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using VoronoiFVM, VoronoiFVM.SolverStrategies
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve
	using LinearAlgebra
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots
	using PlutoUI, HypertextLiteral

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="M-S Transport + Darcy")

# ╔═╡ 57117131-8c83-4788-9f8e-f0b1b9404f0a
function grid1D()
	X=0:0.005:1
	grid=simplexgrid(X)

	# catalyst region
	cellmask!(grid,[0.4],[0.6],2)
	grid
end

# ╔═╡ 87630484-de43-46bc-a179-3e00b6e63c2a
gridplot(grid1D()) 

# ╔═╡ 0a911687-aff4-4c77-8def-084293329f35
begin
	# for 1D domain
	const Γ_left = 1
	const Γ_right = 2
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

where $\rho$ is the mixture density, $\vec v$ is the mass-averaged superficial mixture velocity $\vec N_i$ is the molar flux ($\frac{\text{mol}}{\text{m}^2 \text{s}}$) and $R_i$ is the molar volumetric source/sink ($\frac{\text{mol}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
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

	# kinetic constants for simple reactions
	kp::Float64=1.0
	km::Float64=0.1
		
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
		x[gni[:N2]] = 1.0
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
	k::Float64=1.23e-13*ufac"m^2" # permeability , por class 0


	## Flow data
	#norm conditions
	pn::Float64 = 1.0*ufac"bar"
	Tn::Float64 = 273.15*ufac"K"


	#Qflow::Float64=0.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	Qflow::Float64=148.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	#Qflow::Float64=1480.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	#Qflow::Float64=14800.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		

	MWin::Float64 = molarweight_mix(Fluids, X0)
	mdotin::Float64=MWin*Qflow*pn/(ph"R"*Tn)*ufac"kg/s"
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	p::Float64=1.0*ufac"atm" # reactor pressure

	u0::Float64=Qflow/Ac*ufac"m/s" # mean inlet superficial velocity
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
	(;ip,kp,km,gni,Tamb)=data

	if node.region == 2 # catalyst layer

		# H2 + CO2 <-> H2O + CO
		#r = kp*u[gni[:H2]]/ufac"bar"*u[gni[:CO2]]/ufac"bar"-km*u[gni[:H2O]]/ufac"bar"*u[gni[:CO]]/ufac"bar"		
		
		#f[gni[:H2]]=r
		#f[gni[:CO2]]=r
		#f[gni[:H2O]]=-r
		#f[gni[:CO]]=-r

		# 4 H2 + CO2 <-> 2 H2O + CH4
		#r = kp *(u[gni[:H2]]/ufac"bar")^4 *u[gni[:CO2]]/ufac"bar" -km *(u[gni[:H2O]]/ufac"bar")^2	*u[gni[:CH4]]/ufac"bar"		
		r=zero(eltype(u))
		
		f[gni[:H2]]=4*r
		f[gni[:CO2]]=r
		f[gni[:H2O]]=-2*r
		f[gni[:CH4]]=-r
	end
	
	ng=ngas(data)
	# last species i=ng via ∑pi = p
	#@views f[ng]=u[ip]-sum(u[1:ng])
	#for i=1:ng
		#f[ng] += u[i]
	#end
	#f[ng] = log(f[ng]/u[ip])
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
		#f[i] = -darcyvelo(u,data)*u[i,k]/(ph"R"*Tamb)
		f[i] = -darcyvelo(u,data)*u[i,k]/(ph"R"*Tamb)*embedparam(edge)
	end
	
end

# ╔═╡ 7da59e27-62b9-4b89-b315-d88a4fd34f56
function inlet(f,u,bnode,data)
	# left boundary 1D
	if bnode.region==Γ_left 
		(;ip,X0,u0,pn,Tn,mfluxin)=data
		ng=ngas(data)
		
				
		# flow velocity is normal to top boundary
		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative

		for i=1:ng
			# specify flux at inlet
			f[i] = -X0[i]*u0*pn/(ph"R"*Tn)*embedparam(bnode)
		
		end

		# inlet flow velocity
		#f[ip] = -u0
		# specify inlet mass fluxs
		#f[ip] = -mfluxin
		
	end
end

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
	(;ip,p,pn,Tn,u0,X0,Tamb,mfluxin)=data
	ng=ngas(data)
	
	#inlet(f,u,bnode,data)
	for i=1:ng
		#boundary_dirichlet!(f,u,bnode, species=i,region=Γ_left,value=u[ip]*X0[i])

		boundary_robin!(f,u,bnode, species=i,region=Γ_left,factor=1.0e4,value=1.0e4*u[ip]*X0[i])

		# interpolate between dirichlet and neumann b.c. (specify influx) via embedding
		ϵ = embedparam(bnode) == 0.0 ? 1.0e-5 : embedparam(bnode)

		val_dbc = p*X0[i]
		#val_nbc = -X0[i]*u0*pn/(ph"R"*Tn)
		val_nbc = X0[i]*u0*pn/(ph"R"*Tn)
		
		
		#boundary_robin!(f,u,bnode,i,Γ_left, factor=(1.0-ϵ)/ϵ, value=(1.0-ϵ)/ϵ*( (1.0-ϵ)*val_dbc + ϵ*val_nbc ))
	end
	
	boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_right,value=p)
	#boundary_neumann!(f,u,bnode, species=ip, region=Γ_left,value=mfluxin)
	
	boundary_neumann!(f,u,bnode, species=ip, region=Γ_left,value=mfluxin*embedparam(bnode))
end

# ╔═╡ c0de2aff-8f7f-439c-b931-8eb8fbfcd45d
function mole_frac!(data,X,u::VoronoiFVM.EdgeUnknowns)
	ng=ngas(data)
	sump=zero(eltype(u))
	#X=zeros(eltype(u), ng)
	for i=1:ng
		X[i] = 0.5*(u[i,1]+u[i,2])
		sump += X[i]
	end
	for i=1:ng
		X[i] = X[i] / sump
	end
	nothing
end

# ╔═╡ 51b03d3b-3c7c-4019-8ed8-bf1aaa0b1ddb
function mole_frac!(data,X,u::VoronoiFVM.BNodeUnknowns)
	ng=ngas(data)
	sump=zero(eltype(u))
	#X=zeros(eltype(u), ng)
	for i=1:ng
		X[i] = u[i]
		sump += X[i]
	end
	for i=1:ng
		X[i] = X[i] / sump
	end
	nothing
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

# ╔═╡ b6381008-0280-404c-a86c-9c9c3c9f82eb
#function M_matrix!(M,D,T,p,x,data)
function M_matrix!(M,D,x,data)
	(;Tamb)=data
	ng=ngas(data)
	# !!!ALLOC all methods to be called with arrays to be stack allocated
	# have to  be inlined - here we use callsite inline from Julia 1.8
	#@inline D_matrix!(data, D, T, p)
	@inline D_matrix!(data, D)
	for i=1:ng
		M[i,i] = 1/DK_eff(data,Tamb,i)
		for j=1:ng
			if j != i
				M[i,i] += x[j]/(data.γ_τ*D[i,j])
				M[i,j] = -x[i]/(data.γ_τ*D[i,j])
			end
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
	(;ip,gni,X0,Tamb,k,Fluids) = data

	ng=ngas(data)

	F = MVector{ng,eltype(u)}(undef)
	X = MVector{ng,eltype(u)}(undef)


	M = MMatrix{ng,ng,eltype(u)}(undef)
	D = MMatrix{ng,ng,eltype(u)}(undef)

	#@inline M_matrix!(M,D,Tm,pm,X,data)
	#@inline M_matrix!(M,D,Tamb,0.5*(u[ip,1]+u[ip,2]),X,data)
	#@inline M_matrix!(M,D,X,data)

	rho = rho_mix_flux(u,data)
	vh=darcyvelo(u,data)

	f[ip]=rho*vh
	
	mu = 2.0e-5*ufac"Pa*s"
	
	@inline D_matrix!(data, D)
	for i=1:ng
		#M[i,i] = 1/DK_eff(data,Tamb,i)
		#M[i,i] = 1/(DK_eff(data,Tamb,i) + ph"R"*Tamb*rho*k/Fluids[i].MW/mu)
		M[i,i] = 1/(ph"R"*Tamb*rho*k/Fluids[i].MW/mu)
	
		for j=1:ng
			if j != i
				M[i,i] += (u[j,1]+u[j,2])/(u[ip,1]+u[ip,2])/(data.γ_τ*D[i,j])
				M[i,j] = -(u[i,1]+u[i,2])/(u[ip,1]+u[ip,2])/(data.γ_τ*D[i,j])
			end
		end
		
	end
	#M .*= (ph"R"*Tamb)
	
	
	
	for i=1:ng	
		F[i] = (u[i,1]-u[i,2])/(ph"R"*Tamb)
		#F[i] = (u[i,1]-u[i,2])
	end
	
		
	# computation of fluxes J
	# inplace_linsolve  is made available from VoronoiFVM 1.3.2
	@inline inplace_linsolve!(M,F)


	# index of last species present in mixture
	# apply closing relation for flux calc for this species
	idls = gni[:N2]
	for i=1:ng
		# total flux = diffusive + convective
		#f[i] = F[i] + vh*(u[i,1]+u[i,2])/2/(ph"R"*Tamb)
		f[i] = F[i]
		
		f[idls] += f[i]
	end
	f[idls] = (f[ip] - f[idls])/Fluids[idls].MW

	#f[ng] = (rho*vh - f[ng])/Fluids[ng].MW
	
		
	#@views f[1:ng] .= F
	#f[ng] = zero(eltype(u))
	
	#@views f[ng] = -sum(f[1:(ng-1)])
	#@views f[1:(ng-1)] .= F[1:(ng-1)]	
	
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
	#for i=1:(ng+1)
	for i=1:ng
		f[i]=u[i]/(ph"R"*Tamb)
		#f[i]=u[i]
	end
	# total pressure
	f[ng+1] = rho_mix_storage(u,data)
	#f[ng+1] = u[ng+1]
end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
function main(;data=ModelData(),nref=0)

	grid=grid1D()
	(;ip,p,X0,Tamb)=data
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
	
	enable_species!(sys; species=collect(1:(ng+1))) # gas phase species + p)	

	lg=length(grid1D()[Coordinates])
	
	inival=unknowns(sys)
	inival[:,:].=1.0*p
	#inival[:,:].=1.5*p	
	#inival[ip,:].=linspace(1.1,1.0,lg)
	
	for i=1:ng
		#inival[i,:] .= inival[ip,:]*X0[i]
		inival[i,:] .*= X0[i]
	end

	# transient solver
	#control = SolverControl(nothing, sys;)
	#	control.handle_exceptions=true
	#	control.damp_initial=1
	#	control.damp_growth=1.2
	#	control.force_first_step=true
	#	control.Δt=1.0e-4
	#	control.Δu_opt=1.0e5		
	#	control.Δt_min=1.0e-6
	#	control.Δt_max=1.0e3
	#sol=solve(sys;inival=inival,times,control,verbose="en")

	# steady-state solver
	control = SolverControl(nothing, sys;)
		control.handle_exceptions=true
		control.Δp=1.0e-3
		control.Δp_min=1.0e-6
		control.Δp_grow=1.2
		control.Δu_opt=1.0e5
	
	#sol=solve(sys;inival=inival,control,verbose="en")
	embed=[0.0,1.0]
	sol=solve(sys;inival=inival,control,embed,verbose="en")
		

	sol,grid,sys,data
end;

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	
	sol_,grid,sys,data_embed=main(;data=ModelData());
	
	#if sol_ isa VoronoiFVM.TransientSolution
	#	#sol = copy(sol_(sol_.t[end]))
	#	sol = sol_(t)
	#else
	#	sol = copy(sol_)
	#end
end;
  ╠═╡ =#

# ╔═╡ 2790b550-3105-4fc0-9070-d142c19678db
md"""
## (Partial) Pressure Profiles
"""

# ╔═╡ dcd68d6a-5661-4efd-bbce-f7de99c99d79
times=[0,100_000]

# ╔═╡ 37cd3558-15cc-4133-a0e0-e501a9b73889
#=╠═╡
sol_.t
  ╠═╡ =#

# ╔═╡ 36d02d06-6d0d-4b98-8d60-f2df0afaeaad
#=╠═╡
@bind t Slider(sol_.t,show_value=true)
  ╠═╡ =#

# ╔═╡ 55f305b8-47a2-4fdf-b1cb-39f95f3dfa36
#=╠═╡
sol=sol_(t)
#sol=sol_[1]
  ╠═╡ =#

# ╔═╡ 8fcf636c-2330-4eda-9bb3-298e6a53dfd1
#=╠═╡
let
	(;gni,ip,p)=data_embed
	vis = GridVisualizer(legend = :rt, ylabel="Pressure / Pa")
	scalarplot!(vis, grid, sol[ip,:], label="total")
	scalarplot!(vis, grid, sol[gni[:N2],:], label="N2", color=:green, clear=false)
	scalarplot!(vis, grid, sol[gni[:H2],:], label="H2", color=:red, clear=false)	
	scalarplot!(vis, grid, sol[gni[:CO2],:], label="CO2", color=:blue, clear=false)
	scalarplot!(vis, grid, sol[gni[:CH4],:], label="CO", color=:purple, clear=false)
	scalarplot!(vis, grid, sol[gni[:H2O],:], label="H2O", color=:cyan, clear=false)
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ 87d34a89-d53f-4f27-b8c4-113a6facaad5
#=╠═╡
R=integrate(sys,reaction,sol)
  ╠═╡ =#

# ╔═╡ 8e0fdbc7-8e2b-4ae9-8a31-6f38baf36ef3
md"""
## Flows Over Boundaries
"""

# ╔═╡ 2a4c8d15-168f-4908-b24b-8b65ec3ea494
#=╠═╡
let
	(;pn,Tn,Qflow,X0,mfluxin)=data_embed
	pn*Qflow/(ph"R"*Tn) .* X0 , mfluxin
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
# ╠═87630484-de43-46bc-a179-3e00b6e63c2a
# ╠═0a911687-aff4-4c77-8def-084293329f35
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─9a22118f-cbf0-4412-96dc-66e43a14b0cd
# ╟─83e4a00a-2834-4a68-87ee-81394adea30f
# ╠═ed7941c4-0485-4d84-ad5b-383eb5cae70a
# ╠═78cf4646-c373-4688-b1ac-92ed5f922e3c
# ╠═6779db7a-3823-4a41-8b2d-5558dcd73943
# ╠═a8d57d4c-f3a8-42b6-9681-29a1f0724f15
# ╟─e8b49b2e-d19e-4d46-a399-9972919cd680
# ╠═d4bc847b-052f-4d40-9211-12dbe7e06ee1
# ╠═7da59e27-62b9-4b89-b315-d88a4fd34f56
# ╠═29d66705-3d9f-40b1-866d-dd3392a1a268
# ╟─02b76cda-ffae-4243-ab40-8d0fe1325776
# ╠═c0de2aff-8f7f-439c-b931-8eb8fbfcd45d
# ╠═51b03d3b-3c7c-4019-8ed8-bf1aaa0b1ddb
# ╠═2191bece-e186-4d8e-8a21-3830441baf11
# ╠═b6381008-0280-404c-a86c-9c9c3c9f82eb
# ╠═22447e9e-df5f-4773-a83e-8645ad0c7873
# ╠═dd5e0a30-80ef-4eac-a014-b5a0f6a3c5fe
# ╟─44aa5b49-d595-4982-bbc8-100d2f199415
# ╠═333b5c80-259d-47aa-a441-ee7894d6c407
# ╠═aa498412-e970-45f2-8b11-249cc5c2b18d
# ╟─e25e7b7b-47b3-457c-995b-b2ee4a87710a
# ╠═3a35ac76-e1b7-458d-90b7-d59ba4f43367
# ╟─2790b550-3105-4fc0-9070-d142c19678db
# ╠═dcd68d6a-5661-4efd-bbce-f7de99c99d79
# ╠═37cd3558-15cc-4133-a0e0-e501a9b73889
# ╠═55f305b8-47a2-4fdf-b1cb-39f95f3dfa36
# ╠═8fcf636c-2330-4eda-9bb3-298e6a53dfd1
# ╠═36d02d06-6d0d-4b98-8d60-f2df0afaeaad
# ╠═87d34a89-d53f-4f27-b8c4-113a6facaad5
# ╟─8e0fdbc7-8e2b-4ae9-8a31-6f38baf36ef3
# ╠═eb44075f-fbd3-4717-a440-41cb4eda8de1
# ╠═2a4c8d15-168f-4908-b24b-8b65ec3ea494
# ╠═e1602429-74fa-4949-bb0f-ecd681f52e42
