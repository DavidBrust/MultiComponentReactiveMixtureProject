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
PlutoUI.TableOfContents(title="PC Reactor 2D")

# ╔═╡ 0a911687-aff4-4c77-8def-084293329f35
begin
	# cylindrical geometry, symmetry
	#const Γ_bottom = 1 # outflow bc	
	const Γ_bottom = 5 # outflow bc	
	const Γ_outer = 2 # wall bc	
	const Γ_top = 3 # inflow bc
	const Γ_axis = 4 # symmetry
	#const Γ_cat = 5 # symmetry
end;

# ╔═╡ d95873cc-ad5a-4581-b8d7-b0147eb2491c
function cyl_sym(data; nref=0)
	(;wi,h,cath)=data
	
	#R = linspace(0.0,wi/2.0,9)
	R = linspace(0.0,wi/2.0,30)
 	Z = linspace(0.0,2*h,50)

	grid = simplexgrid(R, Z)
    circular_symmetric!(grid)

	#bfacemask!(grid,[0,h],[wi/2,h],Γ_cat)
	bfacemask!(grid,[wi/2*0,0],[wi/2*0.5,0],Γ_bottom)
	
	grid
end

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
	\vec v  &= -\frac{\kappa}{\mu} \nabla p\\
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
	-\frac{\nabla p_1}{RT} &= \frac{x_2 \vec N_1 - x_1 \vec N_2}{D_{12}} + \frac{x_3 \vec N_1 - x_1 \vec N_3}{D_{13}} \\
	&+ \frac{x_4 \vec N_1 - x_1 \left( \vec m - M_1 \vec N_1 - M_2 \vec N_2 - M_3 \vec N_3 \right)/M_4} {D_{14}}
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

# ╔═╡ 0f455649-cd10-41b4-896d-29b85cdbdb89
md"""
To represent the system of equations in matrix form, we can derive the matrix ``M=(m_{ij})`` with
```math
\begin{align}
	m_{ii} &= \frac{1}{D^K_i} + \frac{x_i M_i}{M_N D^B_{iN}} + \sum_{j\neq i} \frac{x_j}{D^B_{ij}} \quad (i=1\dots N-1) \\ 
	m_{ij} &= -x_i \sum_{j\neq i \atop j \neq N } \frac{1}{D^B_{ij}} - \frac{M_j}{M_N D^B_{iN}}
\end{align}
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
\vec N_N
\end{pmatrix}
=
f
```
"""

# ╔═╡ 9f751b41-9974-4f08-b607-be22806867c3
function darcyvelo(u,μ,data)
	(;ip,k) = data
	k/μ*(u[ip,1]-u[ip,2])	
end

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

# ╔═╡ dfc9b34b-f85f-46b5-910e-4b2c76a0aa96
SolverControl()

# ╔═╡ 778326a2-4ad6-4b8b-b39a-a93a354c5ec8
md"""
# Peclet Number
"""

# ╔═╡ be76db3b-e442-4208-bced-6ab215156291
md"""
The problem of multi-component mass-transport is focused on inter diffusion of species. The numerics struggle to converge in a high Peclet number or convection dominated case. An alternative problem formulation in form of a convective-diffusive flux for species transport might be more appropriate. 
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
	isreactive::Bool=false
		
	ip::Int64=NG+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable
	
	gn::Dict{Int, Symbol} 	= kinpar.gn # names and fluid indices
	gni::Dict{Symbol, Int}  = kinpar.gni # inverse names and fluid indices
	Fluids::Vector{FluidProps} = kinpar.Fluids # fluids and respective properties in system
	
	X0::Vector{Float64} = let
		x=zeros(Float64, NG)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		x[gni[:CO]] = 1.0
		x[gni[:CH4]] = 1.0
		x[gni[:H2O]] = 1.0
		x[gni[:N2]] = 1.0		
		x/sum(x)
	end # inlet composition	
	
	## porous filter data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0

		
	# domain geometry: cylinder
	wi::Float64=16.0*ufac"cm" # corresponds to Diameter
	h::Float64=5.0*ufac"mm"
	#cath::Float64 = 250.0*ufac"μm"
	cath::Float64 = 500.0*ufac"μm"
	Ac::Float64=pi*wi^2.0/4.0*ufac"m^2" # cross-sectional area
	#Ac::Float64=1.0*ufac"m^2" # cross-sectional area
	
	
	#ϕ::Float64=0.36 # porosity, exp determined
	ϕ::Float64=0.33 # porosity, VitraPor sintetered filter class 0
	
	# approximation from Wesselingh, J. A., & Krishna, R. (2006). Mass Transfer in Multicomponent Mixtures
	γ_τ::Float64=ϕ^1.5 # constriction/tourtuosity factor

	#k::Float64=1.23e-10*ufac"m^2" # permeability , por class 0
	k::Float64=1.23e-12*ufac"m^2" # permeability , por class 0
	#k::Float64=1.23e-13*ufac"m^2" # permeability , por class 0

	# Solid Boro-Silikatglas
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	#λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
	λs::Float64=1.13*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
		
	## Flow data
	#norm conditions
	pn::Float64 = 1.0*ufac"bar"
	Tn::Float64 = 273.15*ufac"K"


	#Qflow::Float64=0.0*ufac"ml/minute" # volumetric feed flow rate (sccm)
	#Qflow::Float64=14.80*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	Qflow::Float64=148.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	#Qflow::Float64=1480.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		
	#Qflow::Float64=14800.0*ufac"ml/minute" # volumetric feed flow rate (sccm)		

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

# ╔═╡ 87630484-de43-46bc-a179-3e00b6e63c2a
gridplot(cyl_sym(ModelData()), aspect=5, resolution=(650,600),zoom=1.2)

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

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
	(;ip,iT,p,pn,Tn,u0,X0,Tamb,mfluxin,nfluxin)=data
	ng=ngas(data)

	## Top boundary ##
	
	boundary_neumann!(f,u,bnode, species=ip, region=Γ_top, value=mfluxin*embedparam(bnode))
	
	# specifying the inlet partial pressures 
	for i=1:ng	
		boundary_dirichlet!(f,u,bnode, species=i,region=Γ_top,value=u[ip]*X0[i])
	end
	#if bnode.region == Γ_cat
	#	f[iT] = -10.0*ufac"kW/m^2"*embedparam(bnode)
	#end
	boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_top,value=Tamb)
	
	boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_bottom,value=Tamb)



	## Bottom boundary ##
	boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_bottom,value=p)
	
	
	# apply ∑pi = p for last species (N2) at outflow boundary
    if bnode.region == Γ_bottom
		for i=1:ng
			f[ng] += u[i]
		end
		f[ng] = f[ng]-u[ip]
    end
	
end

# ╔═╡ 202aefe3-1086-4fd3-8f22-8b8e2cee6cec
function mole_frac!(data,X,u::VoronoiFVM.EdgeUnknowns)
	(;ip)=data
	ng=ngas(data)
	#sump=zero(eltype(u))
	
	#for i=1:ng
	#	X[i] = u[i,1]+u[i,2] > 0 ? 0.5*(u[i,1]+u[i,2]) : zero(eltype(u))
	#	sump += X[i]
	#end
	for i=1:ng
		#X[i] = X[i] / sump
		X[i] = u[i,1]+u[i,2] > 0 ? (u[i,1]+u[i,2])/(u[ip,1]+u[ip,2]) : zero(eltype(u))
	end
	nothing
end

# ╔═╡ 6220ceec-5c29-45c0-8bde-0444c916e820
function mole_frac!(data,X,u::VoronoiFVM.BNodeUnknowns)
	(;ip)=data
	ng=ngas(data)
	#sump=zero(eltype(u))
	#X=zeros(eltype(u), ng)
	#for i=1:ng
		#X[i] = u[i]
		#sump += u[i]
	#end
	for i=1:ng
		#X[i] = u[i] / sump 
		X[i] = u[i] > 0 ? u[i]/u[ip] : zero(eltype(u))
	end
	nothing
end

# ╔═╡ d4bc847b-052f-4d40-9211-12dbe7e06ee1
function boutflow(f,u,edge,data)
	(;iT,Tamb)=data
	ng=ngas(data)

	k=outflownode(edge)

	Tm = 0.5*(u[iT,1]+u[iT,2])
	X = MVector{ngas(data),eltype(u)}(undef)
	#@inline mole_frac!(data,X,u)	
	@inline mole_frac!(data,X,u)	
	@inline μ,_=dynvisc_thermcond_mix(data, Tm, X)

	hconv=zero(eltype(u))
	for i=1:ng
		# specify flux at boundary
		
		#f[i] = -darcyvelo(u,data)*u[i,k]/(ph"R"*Tamb)		
		f[i] = -darcyvelo(u,μ,data)*u[i,k]/(ph"R"*u[iT,k])
		hconv += f[i] * enthalpy_gas(data.Fluids[i], u[iT,k])
	end
	#f[iT] = hconv
end

# ╔═╡ 28d0e651-2969-482f-8da5-8bac1c037d38
function D_matrix!(data, u, D)
	(;ip,iT,γ_τ)=data
	ng=ngas(data)
	
	Tm = 0.5*(u[iT,1]+u[iT,2])
	pm = 0.5*(u[ip,1]+u[ip,2])
	for i=1:(ng-1)
		for j=(i+1):ng
			Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], Tm, pm)
			# effective binary diffusion coefficient: porous media correction
			#D[j,i] = Dji
			#D[i,j] = Dji
			D[j,i] = Dji*γ_τ
			D[i,j] = Dji*γ_τ
		end
	end
	nothing
end

# ╔═╡ 7b856860-5eb9-44ea-8244-d8166e8fdc5f
function rho_mix_flux(u,data)
	(;iT)=data
	ng=ngas(data)

	rho=zero(eltype(u))
	for i=1:ng
		rho += data.Fluids[i].MW*0.5*(u[i,1]+u[i,2])
	end
	rho /= ph"R"*0.5*(u[iT,1]+u[iT,2])
end

# ╔═╡ ed7941c4-0485-4d84-ad5b-383eb5cae70a
function flux(f,u,edge,data)
	(;ip,iT,gni,X0,Tamb,k,γ_τ) = data

	ng=ngas(data)

	F = MVector{ng-1,eltype(u)}(undef)
	X = MVector{ngas(data),eltype(u)}(undef)

	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)
	D = MMatrix{ng,ng,eltype(u)}(undef)

	Tm = 0.5*(u[iT,1]+u[iT,2])
	pm = 0.5*(u[ip,1]+u[ip,2])

	@inline mole_frac!(data,X,u)
	@inline μ,λf=dynvisc_thermcond_mix(data, Tm, X)

	# mass continuity
	@inline f[ip] =rho_mix_flux(u,data)*darcyvelo(u,μ,data)
	
	@inline D_matrix!(data, u, D)

	xi=zero(eltype(u))
	xj=zero(eltype(u))

	
	# assemble matrix M and right hand side vector f
	for i=1:(ng-1)
		xi= (u[i,1]+u[i,2])/(u[ip,1]+u[ip,2])
		# Matrix M
		#M[i,i] =  1.0/DK_eff(data,Tamb,i)		
		M[i,i] =  1.0/DK_eff(data,Tm,i)		
	
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
		#F[i] = (u[i,1]-u[i,2])/(ph"R"*Tamb) + xi*f[ip]/data.Fluids[ng].MW/D[i,ng]
		F[i] = (u[i,1]-u[i,2])/(ph"R"*Tm) + xi*f[ip]/data.Fluids[ng].MW/D[i,ng]
		
	end
			
	# computation of fluxes N
	@inline inplace_linsolve!(M,F)

	hconv=zero(eltype(u))
	sumNiMi=zero(eltype(u))
	for i=1:(ng-1)
		f[i] = F[i]
		hconv += f[i] * enthalpy_gas(data.Fluids[i], Tm)
		sumNiMi += f[i] * data.Fluids[i].MW
	end
	hconv += enthalpy_gas(data.Fluids[ng], Tm) * (f[ip] - sumNiMi)/data.Fluids[ng].MW
	
	λbed=kbed(data,λf)*λf
	Bp,Bm = fbernoulli_pm(hconv/λbed/Tm)
	#f[iT]= λbed*(Bm*u[iT,1]-Bp*u[iT,2])
	f[iT]=λbed*(u[iT,1]-u[iT,2])
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
	(;Tamb,iT,ip)=data
	ng=ngas(data)

	cp_mix = zero(eltype(u))
	for i=1:ng
		#f[i]=u[i]/(ph"R"*Tamb)
		f[i]=u[i]/(ph"R"*u[iT])
		cp_mix += u[i]/u[ip] * heatcap_gas(data.Fluids[i], u[iT])
	end
	
	# total pressure
	f[ng+1] = rho_mix_storage(u,data)

	# heat storage term
	f[iT] = cp_mix*u[ip]/ph"R"
end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
function main(;data=ModelData(),nref=0)


	grid=cyl_sym(data)
	
	(;ip,iT,gni,p,Tamb,X0,wi,h)=data
	ng=ngas(data)



	sys=VoronoiFVM.System( 	grid;
						data=data,
						flux=flux,
						reaction=reaction,
						storage=storage,
						bcondition=bcond,					
						boutflow=boutflow,
						outflowboundaries=[Γ_bottom],
						assembly=:edgewise
						)

	enable_species!(sys; species=collect(1:(ng+2))) # pi, ptotal, T		
	inival=unknowns(sys)
	
		
	inival[:,:].=1.0*p	
	for i=1:ng
		inival[i,:] .*= X0[i]
	end
	inival[iT,:] .= Tamb
	
	control = SolverControl(nothing, sys;)
		#control.Δp=0.1
		#control.Δp_grow=2.0
		control.handle_exceptions=true
		control.Δu_opt=1.0e5
		control.abstol_linear=1.0e-10
		

	embed=[0.0,1.0]

	sol=solve(sys;inival=inival,control,embed,verbose="aen")
	sol,grid,sys,data	
end;

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
# ╠═╡ skip_as_script = true
#=╠═╡
sol_,grid,sys,data_embed=main(;data=ModelData());
  ╠═╡ =#

# ╔═╡ c1c935e1-406f-4536-90e9-f11d02e68894
#=╠═╡
let
	data=data_embed
	(;gni,Tamb,p,Ac,Qflow,h,γ_τ)=data
	Dij_B=binary_diff_coeff_gas(data.Fluids[gni[:H2]], data.Fluids[gni[:CO2]], Tamb, p)
	Dij_eff = Dij_B*γ_τ
	
	u=Qflow/Ac
	Pe=h*u/Dij_eff
end
  ╠═╡ =#

# ╔═╡ 2790b550-3105-4fc0-9070-d142c19678db
md"""
## (Partial) Pressure Profiles
"""

# ╔═╡ b9a9b575-4a26-4c79-a8c8-fc3741ae8ec7
md"""
### Including Inert Inlet
"""

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
	(;gni,ip,iT)=data_embed
	vis = GridVisualizer(title="Pressure / Pa", resolution=(650,600), zoom=1.2)
	#idx = ip
	idx = gni[:N2]
	#idx = gni[:H2]
	scalarplot!(vis, grid, sol[idx,:], show=true)
end
  ╠═╡ =#

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
	(;Ac,X0,mfluxin,nfluxin,)=data_embed
	nfluxin.* X0 * Ac, mfluxin * Ac
end
  ╠═╡ =#

# ╔═╡ e1602429-74fa-4949-bb0f-ecd681f52e42
function checkinout(sys,sol)	
	tfact=TestFunctionFactory(sys)
	tf_in=testfunction(tfact,[Γ_bottom],[Γ_top])
	tf_out=testfunction(tfact,[Γ_top],[Γ_bottom])
	(;in=integrate(sys,tf_in,sol),out=integrate(sys,tf_out,sol) )
end

# ╔═╡ 0738a2cc-516d-4a5f-b594-6105bea488ff
#=╠═╡
checkinout(sys,sol)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╠═863c9da7-ef45-49ad-80d0-3594eca4a189
# ╠═d95873cc-ad5a-4581-b8d7-b0147eb2491c
# ╠═87630484-de43-46bc-a179-3e00b6e63c2a
# ╠═0a911687-aff4-4c77-8def-084293329f35
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─9a22118f-cbf0-4412-96dc-66e43a14b0cd
# ╟─83e4a00a-2834-4a68-87ee-81394adea30f
# ╟─c0f65543-a64c-486a-a936-7dc55d75d5f1
# ╟─9b2e8232-c4c2-4c69-acf5-49ce9f6e9270
# ╟─b28f3573-751b-4154-8249-d140baa65d0b
# ╟─e9d0be5a-acbe-44fc-88d8-e1b40c0a46bd
# ╟─0f455649-cd10-41b4-896d-29b85cdbdb89
# ╟─45b34b7b-f816-41ca-9f47-bd05f862db59
# ╟─8490e4f0-04f9-4cfa-a547-0cb4e4fffc59
# ╠═ed7941c4-0485-4d84-ad5b-383eb5cae70a
# ╠═78cf4646-c373-4688-b1ac-92ed5f922e3c
# ╠═6779db7a-3823-4a41-8b2d-5558dcd73943
# ╠═a8d57d4c-f3a8-42b6-9681-29a1f0724f15
# ╠═9f751b41-9974-4f08-b607-be22806867c3
# ╟─e8b49b2e-d19e-4d46-a399-9972919cd680
# ╠═d4bc847b-052f-4d40-9211-12dbe7e06ee1
# ╠═29d66705-3d9f-40b1-866d-dd3392a1a268
# ╟─02b76cda-ffae-4243-ab40-8d0fe1325776
# ╠═202aefe3-1086-4fd3-8f22-8b8e2cee6cec
# ╠═6220ceec-5c29-45c0-8bde-0444c916e820
# ╠═28d0e651-2969-482f-8da5-8bac1c037d38
# ╠═7b856860-5eb9-44ea-8244-d8166e8fdc5f
# ╠═dd5e0a30-80ef-4eac-a014-b5a0f6a3c5fe
# ╟─44aa5b49-d595-4982-bbc8-100d2f199415
# ╠═333b5c80-259d-47aa-a441-ee7894d6c407
# ╠═dfc9b34b-f85f-46b5-910e-4b2c76a0aa96
# ╠═aa498412-e970-45f2-8b11-249cc5c2b18d
# ╟─778326a2-4ad6-4b8b-b39a-a93a354c5ec8
# ╟─be76db3b-e442-4208-bced-6ab215156291
# ╠═c1c935e1-406f-4536-90e9-f11d02e68894
# ╟─e25e7b7b-47b3-457c-995b-b2ee4a87710a
# ╠═3a35ac76-e1b7-458d-90b7-d59ba4f43367
# ╟─2790b550-3105-4fc0-9070-d142c19678db
# ╠═55f305b8-47a2-4fdf-b1cb-39f95f3dfa36
# ╟─b9a9b575-4a26-4c79-a8c8-fc3741ae8ec7
# ╟─8fcf636c-2330-4eda-9bb3-298e6a53dfd1
# ╠═36d02d06-6d0d-4b98-8d60-f2df0afaeaad
# ╟─8e0fdbc7-8e2b-4ae9-8a31-6f38baf36ef3
# ╟─3c4b22b9-4f55-45b9-98d8-6cd929c737c2
# ╠═0738a2cc-516d-4a5f-b594-6105bea488ff
# ╠═2a4c8d15-168f-4908-b24b-8b65ec3ea494
# ╠═e1602429-74fa-4949-bb0f-ecd681f52e42
