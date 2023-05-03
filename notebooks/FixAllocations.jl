### A Pluto.jl notebook ###
# v0.19.24

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

# ╔═╡ 11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using VoronoiFVM
	using ExtendableGrids, GridVisualize
	using LinearAlgebra, LinearSolve, LinearSolvePardiso
    using StaticArrays
	using LessUnitful
	
	using PlutoVista
	using PlutoUI
	using Colors

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="Photo-Catalytic Reactor",depth=5)

# ╔═╡ 3235c072-ed6e-4eb7-821a-523f2520ddda
md"""
Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=true))
"""

# ╔═╡ 2ed3223e-a604-410e-93d4-016580f49093
md"""
# Domain / Grid
"""

# ╔═╡ 0a911687-aff4-4c77-8def-084293329f35
begin
	const Γ_side_front = 1 # symmetry bc
	const Γ_side_right = 2 # wall bc
	const Γ_side_back = 3 # wall bc
	const Γ_side_left = 4 # symmetry bc
	const Γ_bottom = 5 # inflow bc
	const Γ_top_frit = 6 # outflow bc, uncoated porous frit 
	const Γ_top_cat = 7 # outflow bc, catalyst coated porous frit 
end;

# ╔═╡ ada45d4d-adfa-484d-9d0e-d3e7febeb3ef
function prism_sq(data; nref=0, w=data.wi, h=data.h, cath=data.cath, catwi=data.catwi)
	
	hw=(w/2.0)/5.0*2.0^(-nref)
	W=collect(0:hw:(w/2.0))
	
	zCL=h-cath
	hhfrit=zCL/5.0
	Hfrit=collect(0:hhfrit:zCL)
	
	hhCL=cath/5.0
	HCL=collect(zCL:hhCL:h)
	H=glue(Hfrit,HCL)
	
	grid=simplexgrid(W,W,H)
	
	# catalyst layer region
	cellmask!(grid,[0.0,0.0,h-cath],[catwi/2,catwi/2,h],2)
	# catalyst layer boundary
	bfacemask!(grid,[0.0,0.0,h],[catwi/2,catwi/2,h],Γ_top_cat)	
end

# ╔═╡ 2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
md"""
# Governing system
"""

# ╔═╡ f4dcde90-6d8f-4b17-b4ec-367d2372637f
md"""
## Species Mass balances
There are $\nu$ gas phase species in the system, gas phase composition is expressed in partial pressures $p_i, i = 1 ... \nu$.
"""

# ╔═╡ 3703afb0-93c4-4664-affe-b723758fb56b
md"""
```math
\begin{align}
- \nabla N_i + R_i &= 0
~,
i = 1 ... \nu
\end{align}

```
where $N_i$ is the molar flux ($\frac{\text{mol}}{\text{m}^2 \text{s}}$) and $R_i$ is the molar volumetric source/sink ($\frac{\text{mol}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
"""

# ╔═╡ 21d0195b-b170-460d-989e-f9d00b511237
md"""
### Isothermal Dusty Gas Model
"""

# ╔═╡ 8f4843c6-8d2b-4e24-b6f8-4eaf3dfc9bf0
md"""
In the framework of the DGM a system with ``\nu`` gas phase species has ``\nu+1`` unknowns: the species composition in terms of partial pressures ``p_i`` and the total pressure ``p``. The ``\nu`` independent equations expressing the flux-driving force relationships for each of the gas phase species is supplemented by the constraint that the partial pressures sum to the total pressure.
"""

# ╔═╡ 66b55f6b-1af5-438d-aaa8-fe4745e85426
md"""
```math
\begin{align}
	\sum_{j=1 \atop j \neq i}^{\nu} \frac{x_i N_j-x_j N_i}{D_{ij}^{\text{eff}}} - \frac{N_i}{D_{i,\text K}^{\text{eff}}} &= \frac{\nabla p_i}{RT} + \frac{1}{D_{i,K}^{\text{eff}}} \frac{p_i}{RT} \left(\frac{\kappa}{\mu} \nabla p \right) \\
	\sum_{i=1}^\nu p_i &= p
\end{align}
```
"""

# ╔═╡ 8528e15f-cce7-44d7-ac17-432f92cc5f53
md"""
The above form of the DGM neglects thermal diffusion and is formulated for a mixture of ideal gases. The driving force for transport of species $i$ relative to the other species (diffusive transport) then reduces to $\nabla p_i$. Transport of species due to the viscous flow through the porous material in the direction of a negative pressure gradient is captured by D'arcy law (convective transport).
"""

# ╔═╡ a6afe118-dcbd-4126-8646-c7268acfacf3
md"""
The numerical fluxes $\textbf{J}$ could be computed using the backslash operator (Julia solver for linear systems) 
``\textbf{J} 
= 
M
\ \backslash\ 
\textbf{F}``.
However, this allocates a return vector (and internally a pivoting vector).

For a non-allocating variant, VoronoiFVM allows to perform `J.=F ; inplace_linsolve!(M,J)`

In order to be stack allocated, temporary arrays as M, F, J shall be created as `MArray` (available from StaticArrays.jl) with size information known at compile time.
I addition, the call to `inplace_linsolve` needs to be inlined. Conveniently, with Julia 1.8, "call-site  inlining" is available.

An alternative is to use `StrideArray` and `@gc_preserve`, see 
https://discourse.julialang.org/t/what-is-stridearrays-jl/97146/23

Compared to StrideArrays.jl, StaticArrays.jl is the more mature package. Moreover, `@gc_preserve` dose not work for calls with return values. So it seems to be reasonable to stick to `MArray` and inlining, the more with Julia 1.8, callsite inline is aviablable. 


"""

# ╔═╡ a60ce05e-8d92-4172-b4c1-ac3221c54fe5
md"""
### Knudsen effective Diffusivity
"""

# ╔═╡ 24374b7a-ce77-45f0-a7a0-c47a224a0b06
md"""
In the DGM the solid porous matrix is considered as another species in the ideal gas mix. The Knudsen diffusion coefficients describe the interactions between the gas molecules and the solid porous matrix. Analogous to the gas-gas diffusion coefficients, the Knudsen diffusion coefficients quantify the resulting friction force from momentum transfer upon collision with the walls acting on the gas molecules.
"""

# ╔═╡ 4865804f-d385-4a1a-9953-5ac66ea50057
md"""
Calculation of Knudsen diffusion coefficients according to __Wesselingh, J. A., & Krishna, R. (2006).__ Mass Transfer in Multicomponent Mixtures (ch. 21, Fig. 21.8)
"""

# ╔═╡ 3bf71cea-4f73-47da-b5ed-2cae3ec3d18b
md"""
## Thermal Energy Transport
"""

# ╔═╡ 7f94d703-2759-4fe1-a8c8-ddf26732a6ca
md"""
```math
\begin{align}
- \nabla \left(\lambda_{\text{eff}} \nabla T - c_p \left( T-T_{\text{ref}} \right) \frac{p}{RT} \vec u_{\text D}\right) + R_{\text{th}} &= 0

\end{align}

```

where ``\lambda_{\text{eff}}`` is the effective thermal conductivity in ``\frac{\text{W}}{\text{mK}}``, ``c_p`` is the isobaric heat capacity of an ideal gas in ``\frac{\text J}{\text{molK}}``, ``\vec u_{\text D}`` is the D'arcy velocity in ``\frac{\text m}{\text s}`` and ``R_{\text{th}}`` is the volumetric heat source/sink term ``\frac{\text{W}}{\text{m}^3}`` from chemical reactions.
"""

# ╔═╡ 906ad096-4f0c-4640-ad3e-9632261902e3
md"""
# Boundary Conditions
"""

# ╔═╡ 39e74955-aab6-4bba-a1b8-b2307b45e673
md"""
## Thermal
"""

# ╔═╡ 6798d5e2-b8c7-4f54-aa71-6ea1ccab78fb
md"""
The thermal boundary conditions consist of radiative losses (reflection of incoming light, thermal emission of IR), conductive heat transfer to the reactor walls (described by a wall heat transfer coefficient) as well as enthalpy in and outflows through the top and bottom surfaces of the modeling domain. At the symmetry boundaries there are no-flux (homogeneous Neumann) boundary conditions.
"""

# ╔═╡ ed3609cb-8483-4184-a385-dca307d13f17
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
$(LocalResource("../img/ThermalBC.png")) 
"""
  ╠═╡ =#

# ╔═╡ 8139166e-42f9-41c3-a360-50d3d4e5ee86
md"""
## Gas species
"""

# ╔═╡ 44d91c2e-8082-4a90-89cc-81aba783d5ac
md"""
Boundary conditions for the transport of gas phase species cover in and outflow boundary conditions at the bottom and top surfaces of the modelling domain with no-flux conditins applied elsewhere. In the catalyst layer, volumetric catalytic reactions take place.
"""

# ╔═╡ 560feb9c-bd0f-414c-bb6f-ef5fd7cb3c49
MVector

# ╔═╡ edd9fdd1-a9c4-4f45-8f63-9681717d417f
function side(f,u,bnode,data)
	# side wall boundary condition
	(;iT,α_w,Tamb)=data
	
	boundary_robin!(f,u,bnode;species=iT,region=(Γ_side_back,Γ_side_right), factor=α_w, value=Tamb*α_w)	
end

# ╔═╡ 02b76cda-ffae-4243-ab40-8d0fe1325776
md"""
# Auxiliary functions
"""

# ╔═╡ 44aa5b49-d595-4982-bbc8-100d2f199415
md"""
# System Setup and Solution
"""

# ╔═╡ 7759971e-2f5a-439b-a6cf-8150538a0034
md"""
Reactive: $(@bind reactive PlutoUI.CheckBox(default=true))
"""

# ╔═╡ e25e7b7b-47b3-457c-995b-b2ee4a87710a
md"""
## Model Data
"""

# ╔═╡ 3a35ac76-e1b7-458d-90b7-d59ba4f43367
begin
	# !!!ALLOC We need ng as a "static" parameter (known at compile time)
	# !!!ALLOC so we make it a type parameter instead of a struct entry
    Base.@kwdef mutable struct ModelData{NG} <:AbstractModelData
	
	# catalyst / chemistry data
	kinpar::FixedBed.KinData{nreac(S3P)} = S3P
	#kinpar::AbstractKineticsData = XuFroment1989
	
	# names and fluid indices
	gn::Dict{Int, Symbol} 	= S3P.gn
	# inverse names and fluid indices
	gni::Dict{Symbol, Int}  = S3P.gni
	# fluids and respective properties in system
	Fluids::Vector{FluidProps} = S3P.Fluids
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = let
		x=zeros(Float64, NG)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		x/sum(x)
	end # inlet composition

	# volume specific cat mass loading, UPV lab scale PC reactor
	mcats::Float64 =1234.568*ufac"kg/m^3"
	isreactive::Bool = 1
	#isreactive::Bool = 0

	ip::Int64=NG+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable

	α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	
	## irradiation data
	Tau_quartz::Float64=0.9 # transmission coefficient of quartz window
	G_lamp::Float64=1.0*ufac"kW/m^2" # solar simulator irradiation flux
	Abs_lamp_cat::Float64=0.7 # cat avg absorptivity of irradiation from lamp
	Abs_lamp_frit::Float64=0.7 # frit avg absorptivity of irradiation from lamp
	Eps_ir::Float64=0.7 # avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
		
	## porous filter data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0
	#dp::Float64=100.0*ufac"μm" # average pore size, por class 2

	# frit thickness (applies to 2D & 3D)
	h::Float64=0.5*ufac"cm"
	# catalyst layer thickness (applies to 2D & 3D)
	cath::Float64 = 250.0*ufac"μm"
	

	# prism / 3D
	wi::Float64=12.0*ufac"cm" # prism width/side lenght
	le::Float64=wi # prism width/side lenght
	catwi::Float64=10.0*ufac"cm" # prism width/side lenght
	Ac::Float64=wi*le*ufac"m^2" # cross-sectional area, square

	# Solid Boro-Silikatglas
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2 	
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
	
	#ϕ::Float64=0.36 # porosity, class 2
	ϕ::Float64=0.33 # porosity, class 0
	
	# approximation from Wesselingh, J. A., & Krishna, R. (2006). Mass Transfer in Multicomponent Mixtures
	γ_τ::Float64=ϕ^1.5 # constriction/tourtuosity factor

	#k::Float64=2.9e-11*ufac"m^2" # permeability , por class 2
	k::Float64=1.23e-10*ufac"m^2" # permeability , por class 0
	
	# a_s::Float64=0.13*ufac"m^2/g" # specific surface area, por class 2
	a_s::Float64=0.02*ufac"m^2/g" # specific surface area, por class 0
	
	ρfrit::Float64=(1.0-ϕ)*ρs*ufac"kg/m^3" # density of porous frit
	a_v::Float64=a_s*ρfrit*ufac"m^2/m^3" # volume specific interface area
	## END porous filter data


	## Flow data
	#norm conditions
	pn::Float64 = 1.0*ufac"bar"
	Tn::Float64 = 273.15*ufac"K"
	
	Qflow::Float64=3400.0*ufac"ml/minute" # volumetric feed flow rate (sccm)
	#Qflow::Float64=20000.0*ufac"ml/minute" # volumetric feed flow rate (sccm)

	MWin::Float64 = molarweight_mix(Fluids, X0)
	mdotin::Float64=MWin*Qflow*pn/(ph"R"*Tn)*ufac"kg/s"

	
	Tin::Float64=298.15*ufac"K" # inlet temperature
	Tamb::Float64=Tin # ambient temperature
	p::Float64=1.0*ufac"atm" # reactor pressure

	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	utop::Float64=u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate
	
end;

# !!!ALLOC Method to be called instead of data.ng
FixedBed.ngas(::ModelData{NG}) where NG = NG

# !!!ALLOC Additional constructo taking ng as parameter	
ModelData(;ng=S3P.ng, kwargs...) = ModelData{ng}(;kwargs...)
end;

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	(;ip,isreactive)=data
	ng=ngas(data)
	if  node.region == 2 && isreactive # catalyst layer
		(;iT,mcats,kinpar)=data
		(;rni,nuij,ΔHi)=kinpar
		rni=data.kinpar.rni # reaction indices from names
		pi = MVector{ngas(data),eltype(u)}(undef)
		
		for i=1:ng
			pi[i] = u[i]/ufac"bar"
		end
		# Ki =MVector{nreac(data),eltype(u)}(undef)
		
		# negative sign: sign convention of VoronoiFVM: source term < 0 
		# ri returns reaction rate in mol/(h gcat)

		RR = @inline -mcats*ri(data,u[iT],pi)*ufac"mol/(hr*g)"

		for i=1:ng
			f[i] = zero(eltype(u))
			for j=1:nreac(kinpar)
				f[i] += nuij[i,j] * RR[j]
			end			
		end
				
		# temperature eq. / heat source
		f[iT] = @inbounds -(RR[rni[:R1]]*ΔHi[:R1]+RR[rni[:R2]]*ΔHi[:R2] +RR[rni[:R3]]*ΔHi[:R3])
	end
	
	# ∑xi = 1
	@views f[ip]=u[ip]-sum(u[1:ng])
end

# ╔═╡ 40906795-a4dd-4e4a-a62e-91b4639a48fa
function bottom(f,u,bnode,data)
	if bnode.region==Γ_bottom # bottom boundary
		(;iT,Eps_ir,Tamb,u0,X0,pn,Tn)=data
		ng=ngas(data)
		f[iT] = Eps_ir*ph"σ"*(u[iT]^4 - Tamb^4)
		
		for i=1:ng
			# specify flux at boundary: flow velocity is normal to bot boundary	
			f[i] = -u0*X0[i]*pn/(ph"R"*data.Tn)
			
		end
	end
end

# ╔═╡ 7d705c9c-6653-4050-b32e-ea96ff972657
function mole_frac!(node,data,X,u::VoronoiFVM.EdgeUnknowns)
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

# ╔═╡ e4997de1-672a-4060-97b8-b4aa000f5112
function mole_frac!(node,data,X,u::VoronoiFVM.BNodeUnknowns)
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

# ╔═╡ 7da59e27-62b9-4b89-b315-d88a4fd34f56
function top(f,u,bnode,data)
	# top boundaries (cat layer & frit)
	if bnode.region==Γ_top_frit || bnode.region==Γ_top_cat 
		(;ip,iT,Tamb,Fluids,utop,Eps_ir,Abs_lamp_frit,Abs_lamp_cat,Tau_quartz,G_lamp)=data
		ng=ngas(data)
		flux_rerad = Eps_ir*ph"σ"*(u[iT]^4 - Tamb^4)
		
		X =MVector{ngas(data),eltype(u)}(undef)
		@inline mole_frac!(bnode,data,X,u)
		
		@inline cf=heatcap_mix(Fluids, u[iT], X)		
		flux_convec=utop*u[ip]/(ph"R"*u[iT])*cf*(u[iT]-Tamb)

		abs = 0.0
		if bnode.region==Γ_top_frit
			abs=Abs_lamp_frit
		else # catalyst layer
			abs=Abs_lamp_cat
		end
		f[iT] = -abs*Tau_quartz*G_lamp + flux_rerad + flux_convec
		
		
		# flow velocity is normal to top boundary
		for i=1:ng
			f[i] = utop*u[i]/(ph"R"*u[iT])
		end
	end
end

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
    top(f,u,bnode,data)
    bottom(f,u,bnode,data)
    side(f,u,bnode,data)	
end

# ╔═╡ bdf5e360-cc20-443a-8929-5cd2b64fd749
function D_matrix!(data, D, T, p)
	ng=ngas(data)
	for i=1:(ng-1)
		for j=(i+1):ng
			Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
			D[j,i] = Dji
			D[i,j] = Dji
			#v[j,i] = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
		end
	end
	nothing
end

# ╔═╡ 101b88bc-f382-4c60-a890-6b9e0265d206
function M_matrix!(M,D,T,p,x,data)
	ng=ngas(data)
	# !!!ALLOC all methods to be called with arrays to be stack allocated
	# have to  be inlined - here we use callsite inline from Julia 1.8
	@inline D_matrix!(data, D, T, p)
	for i=1:ng
		M[i,i] = -1/DK_eff(data,T,i)
		for j=1:ng
			if j != i
				M[i,i] -= x[j]/(data.γ_τ*D[i,j])
				M[i,j] = x[i]/(data.γ_τ*D[i,j])
			end
		end	
	end
	nothing
end

# ╔═╡ ed7941c4-0485-4d84-ad5b-383eb5cae70a
    function flux(f,u,edge,data)
    (;ip,iT,k,Fluids,Tamb) = data
        ng=ngas(data)

        # !!!ALLOC Use MVector with static size information instef of Vector
        F = MVector{ngas(data),eltype(u)}(undef)
	X = MVector{ngas(data),eltype(u)}(undef)
        
	pk,pl = u[ip,1],u[ip,2]
	δp = pk-pl
	pm = 0.5*(pk+pl)
        
	Tm=0.5*(u[iT,1]+u[iT,2])
        f[iT]=((u[iT,1]-Tamb)-(u[iT,2]-Tamb))	

	# !!!ALLOC All functions which take MVectors to be stack allocated
        # as argument must be inlined. Here, we use the new "callsite inline"
        # wich is new with Julia 1.8.
        
	@inline mole_frac!(edge,data,X,u)
	@inline cf=heatcap_mix(Fluids, Tm, X)
	@inline μ,λf=dynvisc_thermcond_mix(data, Tm, X)
 	λbed=kbed(data,λf)*λf
        
	# Darcy flow
	ud=-k/μ * δp
	# vh=project(edge,(0,ud)) # 2D
	vh=project(edge,(zero(eltype(u)),zero(eltype(u)),ud)) # 3D
	# convective enthalpy flux
	conv=vh*pm/(ph"R"*Tm)*cf/λbed
	
	Bp,Bm = fbernoulli_pm(conv)
	# temperature flux
	f[iT]= λbed*(Bm*(u[iT,1]-Tamb)-Bp*(u[iT,2]-Tamb))		

        # !!!ALLOC Use MMatriy with static size information instef of Matrix
	M = MMatrix{ngas(data),ngas(data),eltype(u)}(undef)
	D = MMatrix{ngas(data),ngas(data),eltype(u)}(undef)
	@inline M_matrix!(M,D,Tm,pm,X,data)
   	
	for i=1:ng
	    DK = DK_eff(data,Tm,i)
	    bp,bm=fbernoulli_pm(vh/DK)		
	    F[i] = -(bm*u[i,1]-bp*u[i,2])/(ph"R"*Tm)
	end
	
	# computation of fluxes J

	# !!!ALLOC All functions which take MMatriy to be stack allocated
        # as argument must be inlined. Here, we use the new "callsite inline"
        # wich is new with Julia 1.8
        # inplace_linsolve  is made available from VoronoiFVM 1.3.2
       @inline inplace_linsolve!(M,F)

	
	@views f[1:ng] .= F	
    #f[ip] via reaction: ∑pi = p

end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
function main(;data=ModelData())

	grid=prism_sq(data)
	
	iT=data.iT
	
	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)
	enable_species!(sys; species=collect(1:(ngas(data)+2))) # gas phase species + p + T
	
	inival=unknowns(sys)
	inival[:,:].=1.0*data.p
	for i=1:ngas(data)
		inival[i,:] .*= data.X0[i]
	end
	inival[iT,:] .= data.Tamb

	sol=solve(sys;inival,)
	
	 function pre(sol,par)
	 	ng=ngas(data)
	 	iT=data.iT
	 	# iteratively adapt top outflow boundary condition
	 	function Inttop(f,u,bnode,data)
			
	 		X=zeros(eltype(u),ngas(data))
			mole_frac!(bnode,data,X,u)
	 		# top boundary(cat/frit)
	 		if bnode.region==Γ_top_frit || bnode.region==Γ_top_cat  
	 			for i=1:ng
	 				f[i] = data.Fluids[i].MW*X[i]
	 			end
	 			f[iT] = u[iT]
	 		end
	 	end
		
	 	MWavg=sum(integrate(sys,Inttop,sol; boundary=true)[1:ng,[Γ_top_frit,Γ_top_cat]])/(data.Ac/4)
	 	ntop=data.mdotin/MWavg
		 
	 	Tavg=sum(integrate(sys,Inttop,sol; boundary=true)[data.iT,[Γ_top_frit,Γ_top_cat]])/(data.Ac/4)
		
	 	utop_calc=ntop*ph"R"*Tavg/(1.0*ufac"bar")/data.Ac
	 	utops=[data.utop, utop_calc]*ufac"m/s"
	 	data.utop = minimum(utops) + par*(maximum(utops)-minimum(utops))
		
		
	 	# specific catalyst loading
	 	mcats=[10.0, 1300.0]*ufac"kg/m^3"
	 	data.mcats= minimum(mcats) + par*(maximum(mcats)-minimum(mcats))

	 	# irradiation flux density
	 	G_lamp=[1.0, 100.0]*ufac"kW/m^2"
	 	data.G_lamp= minimum(G_lamp) + par*(maximum(G_lamp)-minimum(G_lamp))
	 end
	
	 control=SolverControl( ;
	 				  		handle_exceptions=true,
							Δp_min=1.0e-4,					  
	 				  		Δp=0.1,
	 				  		Δp_grow=1.5,
	 				  		Δu_opt=100000.0, # large value, due to unit Pa of pressure?
							verbose="an", # log allocations and newton convergence
		 					method_linear = KrylovJL_BICGSTAB(),
                     		precon_linear = A -> factorize(A, UMFPACKFactorization())
	 				  		)
	
	embed=[0.0,1.0]
	sol=solve(sys;inival,embed,pre,control)
	
	sol,grid,sys,data
end;

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	if RunSim
		sol_,grid,sys,data_embed=main(data=ModelData(isreactive=reactive));
		if sol_ isa VoronoiFVM.TransientSolution
			sol = copy(sol_(sol_.t[end]))
		else
			sol = copy(sol_)
		end
	end
end;
  ╠═╡ =#

# ╔═╡ 985718e8-7ed7-4c5a-aa13-29462e52d709
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Cutplane at ``z=`` $(@bind zcut Slider(range(0.0,data_embed.h,length=101),default=data_embed.h,show_value=true)) m
"""
  ╠═╡ =#

# ╔═╡ ff58b0b8-2519-430e-8343-af9a5adcb135
# ╠═╡ skip_as_script = true
#=╠═╡
let
	vis=GridVisualizer(resolution=(600,400), zoom=1.9, )
	grid_=prism_sq(ModelData{S3P.ng}(),nref=1)

	gridplot!(vis, grid_, zplane=zcut)
	reveal(vis)
	
	#save("../img/out/domain.svg", vis)
end
  ╠═╡ =#

# ╔═╡ bcaf83fb-f215-428d-9c84-f5b557fe143f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
To calculate the temperature field within the modelling domain, the thermal energy equation is solved. The modelling domain consists of the porous frit material (SiO2, silica-glas) and the gaseous phase in the void volume. The domain is considered a "Quasi-homogenous" phase with ``T_{\text{s}}=T_{\text{f}}`` (i.e. porous solid and gas phase are in thermal equil.) This assumption is based on:
- Very high specific surface area of porous material ``a_{\text V} = ``$(data_embed.a_v) ``\frac{\text m^2}{\text m^3}``
- Large interphase heat transfer coefficient (Kuwahara, F., Shirota, M., & Nakayama, A. (2001). doi:10.1016/s0017-9310(00)00166-6)


Heat is transported within the domain via conduction and convective transport. Because of the treatment of the porous domain as a Quasi-homogenous phase, an effective thermal conductivity is used to describe conduction. The velocity of convective gas transport within the porous medium is the D'arcy velocity.
"""
  ╠═╡ =#

# ╔═╡ b2df1087-6628-4889-8cd6-c5ee7629cd93
md"""
## Temperature Plot
"""

# ╔═╡ 3bd80c19-0b49-43f6-9daa-0c87c2ea8093
#=╠═╡
let
	(;ip,iT)=data_embed
    vis=GridVisualizer()
	scalarplot!(vis, grid, sol[iT,:].- 273.15, levelalpha=0, show=true)
	#scalarplot!(vis, grid, sol[ip,:]./ufac"bar", levelalpha=0.5, show=true)

end
  ╠═╡ =#

# ╔═╡ f8e6847c-c021-4eee-89f1-f0453ed1bf40
FixedBed.XuFroment

# ╔═╡ 2790b550-3105-4fc0-9070-d142c19678db
md"""
## Partial Pressure Plot
"""

# ╔═╡ 71f06a02-d27b-4057-8e82-438fe88bfd21
#=╠═╡
let
	data=data_embed
	grid=prism_sq(data)
	(;wi,h,iT) =data
	y_=0.06
	
	bid = maximum(grid[BFaceRegions])+1
	bfacemask!(grid, [0,y_,0],[wi/2,y_,h],bid)

	# keep x-z coordinates of parent grid
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[3]
	end
	#grid_1D  = subgrid(grid, [bid], boundary=true, transform=_3to1) 
	grid2D  = subgrid(grid, [bid], boundary=true, transform=_3to2) 
	grid2D[Coordinates]
	sol2D = view(sol[iT, :], grid2D)
	#sol
	scalarplot(grid2D, sol2D)
	

	#gridplot(grid,yplane=0.024)

	
	#bedgemask!(    grid::ExtendableGrid,    xa,    xb,    ireg::Int64;    tol) -> ExtendableGrid
end
  ╠═╡ =#

# ╔═╡ bd7552d2-2c31-4834-97d9-ccdb4652242f
function SolAlongLine(data,sol)
		
	grid=prism_sq(data)
	mid_x=argmin(abs.(grid[Coordinates][1,:] .-data.wi/4))
	mid_x=grid[Coordinates][1,mid_x]
	mid_y=argmin(abs.(grid[Coordinates][2,:] .-data.le/4))
	mid_y=grid[Coordinates][2,mid_y]
	

	Nodes = findall(x->x[1] == mid_x && x[2] == mid_y, eachcol(grid[Coordinates]))
	
	grid1D = grid[Coordinates][:,Nodes]
	grid1D = grid1D[3,:] # extract z-coordinate


        sol_p = []
	for i=1:(ngas(data)+2)
		push!(sol_p, sol[i,Nodes])
	end

	sol_p,grid1D,mid_x,mid_y
end

# ╔═╡ bea97fb3-9854-411c-8363-15cbef13d033
#=╠═╡
let
	sol_, grid_,midx,midy = SolAlongLine(data_embed,sol)
	
	(;gn, ip)=data_embed
	c1 = colorant"red"
	c2 = colorant"blue"
	cols=range(c1, stop=c2, length=ngas(data_embed)+1)

	vis=GridVisualizer(title="At position x= $(midx), y= $(midy)")
	#vis=GridVisualizer(Plotter=PyPlot)
	for i in 1:ngas(data_embed)
		scalarplot!(vis, grid_, sol_[i], color=cols[i], label="$(gn[i])", clear=false, legend=:best)
	end

	scalarplot!(vis, grid_, sol_[ip], color=cols[ip], label="total p", clear=false,legend=:best, xlabel="Height / m", ylabel="Pressure / Pa")
	
	reveal(vis)
	#save("../img/out/pi_pt.svg", vis)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╟─863c9da7-ef45-49ad-80d0-3594eca4a189
# ╠═3235c072-ed6e-4eb7-821a-523f2520ddda
# ╟─2ed3223e-a604-410e-93d4-016580f49093
# ╠═ada45d4d-adfa-484d-9d0e-d3e7febeb3ef
# ╠═0a911687-aff4-4c77-8def-084293329f35
# ╠═ff58b0b8-2519-430e-8343-af9a5adcb135
# ╟─985718e8-7ed7-4c5a-aa13-29462e52d709
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─f4dcde90-6d8f-4b17-b4ec-367d2372637f
# ╟─3703afb0-93c4-4664-affe-b723758fb56b
# ╟─21d0195b-b170-460d-989e-f9d00b511237
# ╟─8f4843c6-8d2b-4e24-b6f8-4eaf3dfc9bf0
# ╟─66b55f6b-1af5-438d-aaa8-fe4745e85426
# ╟─8528e15f-cce7-44d7-ac17-432f92cc5f53
# ╠═ed7941c4-0485-4d84-ad5b-383eb5cae70a
# ╟─a6afe118-dcbd-4126-8646-c7268acfacf3
# ╠═78cf4646-c373-4688-b1ac-92ed5f922e3c
# ╟─a60ce05e-8d92-4172-b4c1-ac3221c54fe5
# ╟─24374b7a-ce77-45f0-a7a0-c47a224a0b06
# ╟─4865804f-d385-4a1a-9953-5ac66ea50057
# ╟─3bf71cea-4f73-47da-b5ed-2cae3ec3d18b
# ╟─bcaf83fb-f215-428d-9c84-f5b557fe143f
# ╟─7f94d703-2759-4fe1-a8c8-ddf26732a6ca
# ╟─906ad096-4f0c-4640-ad3e-9632261902e3
# ╟─39e74955-aab6-4bba-a1b8-b2307b45e673
# ╟─6798d5e2-b8c7-4f54-aa71-6ea1ccab78fb
# ╟─ed3609cb-8483-4184-a385-dca307d13f17
# ╟─8139166e-42f9-41c3-a360-50d3d4e5ee86
# ╟─44d91c2e-8082-4a90-89cc-81aba783d5ac
# ╠═560feb9c-bd0f-414c-bb6f-ef5fd7cb3c49
# ╠═7da59e27-62b9-4b89-b315-d88a4fd34f56
# ╠═40906795-a4dd-4e4a-a62e-91b4639a48fa
# ╠═edd9fdd1-a9c4-4f45-8f63-9681717d417f
# ╠═29d66705-3d9f-40b1-866d-dd3392a1a268
# ╟─02b76cda-ffae-4243-ab40-8d0fe1325776
# ╠═7d705c9c-6653-4050-b32e-ea96ff972657
# ╠═e4997de1-672a-4060-97b8-b4aa000f5112
# ╠═bdf5e360-cc20-443a-8929-5cd2b64fd749
# ╠═101b88bc-f382-4c60-a890-6b9e0265d206
# ╟─44aa5b49-d595-4982-bbc8-100d2f199415
# ╠═333b5c80-259d-47aa-a441-ee7894d6c407
# ╟─7759971e-2f5a-439b-a6cf-8150538a0034
# ╠═aa498412-e970-45f2-8b11-249cc5c2b18d
# ╟─e25e7b7b-47b3-457c-995b-b2ee4a87710a
# ╠═3a35ac76-e1b7-458d-90b7-d59ba4f43367
# ╟─b2df1087-6628-4889-8cd6-c5ee7629cd93
# ╠═3bd80c19-0b49-43f6-9daa-0c87c2ea8093
# ╠═f8e6847c-c021-4eee-89f1-f0453ed1bf40
# ╟─2790b550-3105-4fc0-9070-d142c19678db
# ╠═bea97fb3-9854-411c-8363-15cbef13d033
# ╠═71f06a02-d27b-4057-8e82-438fe88bfd21
# ╠═bd7552d2-2c31-4834-97d9-ccdb4652242f
