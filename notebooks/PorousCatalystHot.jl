### A Pluto.jl notebook ###
# v0.19.22

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
	
	using VoronoiFVM
	using ExtendableGrids, GridVisualize
	using LinearAlgebra
	using LessUnitful
	using PlutoVista	
	using PlutoUI
	#using Plots
	using Colors
	#using CSV, DataFrames
	using Revise
	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="Dusty Gas Model")

# ╔═╡ b659879f-3706-4937-b655-116f1be85069
begin
	const Γ_Bot=1
	const Γ_Top=2
	const Γ_Cat=3
end;

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
- \nabla J_{\alpha} + R_{\alpha} &= 0
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
The above form of the DGM neglects thermal diffusion and external forces and is formulated for a mixture of ideal gases. The driving force for transport of species $i$ relative to the other species (diffusive transport) then reduces to $\nabla p_i$. Transport of species due to the viscous flow through the porous material in the direction of a negative pressure gradient is captured by D'arcy law (convective transport).
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

# ╔═╡ a60ce05e-8d92-4172-b4c1-ac3221c54fe5
md"""
## Knudsen effective Diffusivity
"""

# ╔═╡ 24374b7a-ce77-45f0-a7a0-c47a224a0b06
md"""
In the DGM the solid porous matrix is considered as another species in the ideal gas mix. The Knudsen diffusion coefficients describe the interactions between the gas molecules and the solid porous matrix. Analogous to the gas-gas diffusion coefficients, the Knudsen diffusion coefficients quantify the resulting friction force from momentum transfer upon collision with the walls acting on the gas molecules.
"""

# ╔═╡ 4865804f-d385-4a1a-9953-5ac66ea50057
md"""
Calculation of Knudsen diffusion coefficients according to __Wesselingh, J. A., & Krishna, R. (2006).__ Mass Transfer in Multicomponent Mixtures (ch. 21, Fig. 21.8)
"""

# ╔═╡ 2fb1a154-2721-4da6-848e-99d46dfce774
md"""
### Darcy's law
Calculate gas velocity in porous medium. This gas velocity is responsible for convective transport.
"""

# ╔═╡ cc12f23a-4a80-463b-baf8-22f58d341d9d
md"""
```math
	\vec u_{\text{g}} = - \frac{\kappa^{\text{eff}}}{\mu_{\text{g}}} \nabla p
```
"""

# ╔═╡ 91d8119f-501f-49f2-93e0-88da8d996f7a
function vel_darcy(data, ∇p, T, X)
	μ=dynvisc_mix(data, T, X)
	data.k/μ * ∇p # m/s
end

# ╔═╡ 4af2237c-9144-4ffc-8966-2b3bf9d3d720
function mole_frac!(data,X,u)
	n=data.ng
	sump = 0.0
	for i=1:n
		X[i] = 0.5*(u[i,1]+u[i,2])
		sump += X[i]
	end
	X .= X / sump
	nothing
end

# ╔═╡ 2191bece-e186-4d8e-8a21-3830441baf11
function D_matrix(data, T, p)
	n=data.ng
	v = zeros(Float64,n,n)
	
	for i=1:(n-1)
		for j=(i+1):n
			v[j,i] = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
		end
	end
	Symmetric(v, :L)
end

# ╔═╡ 02b76cda-ffae-4243-ab40-8d0fe1325776
md"""
##### Auxiliary functions
"""

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	ngas=data.ng
	n=data.gni
	ip=data.ip
	RRc=data.RRc

	if node.region == 2 && data.isreactive # catalyst layer
		RR=RRc*u[n["H2"]]*u[n["CO2"]]
		f[n["H2"]] = RR
		f[n["CO2"]] = RR
		f[n["CO"]] = -RR
		f[n["H2O"]] = -RR
	end
	
	# ∑xi = 1
	f[ip]=u[ip]-sum(u[1:ngas])
	#f[ip]=log(sum(u[1:ngas])) -log(u[ip])
	
end

# ╔═╡ 906ad096-4f0c-4640-ad3e-9632261902e3
md"""
## Boundary Conditions
"""

# ╔═╡ 7da59e27-62b9-4b89-b315-d88a4fd34f56
function top(f,u,bnode,data)
	
	if bnode.region==3 # top boundary
		
		iT = data.iT
		flux_rerad = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
		ρf=density_idealgas(data.Fluids, u[iT], data.p, data.X0)
		cf=heatcap_mix(data.Fluids, u[iT], data.X0)
		flux_convec=data.utop*ρf*cf*(u[iT]-data.Tamb)
		f[iT] = -data.Abs_lamp*data.G_lamp + flux_rerad + flux_convec
		
		
		ng=data.ng
		# flow velocity is normal to top boundary
		for i=1:data.ng
			f[i] = data.utop*u[i]/(ph"R"*data.Tin)
		end
	end
end

# ╔═╡ 40906795-a4dd-4e4a-a62e-91b4639a48fa
function bottom(f,u,bnode,data)
	if bnode.region==1 # bottom boundary
		iT=data.iT
		f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
	end
end

# ╔═╡ edd9fdd1-a9c4-4f45-8f63-9681717d417f
function side(f,u,bnode,data)
	# side wall boundary condition
	iT=data.iT
	boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)	
end

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
	ip=data.ip
	p0=data.p
	ng=data.ng
	X0=data.X0
	#Δp=data.Δp
	# set partial pressures of inlet comp at bottom boundary
	for i=1:ng
		boundary_dirichlet!(f,u,bnode,i,1,X0[i]*p0)
	end

	top(f,u,bnode,data)
	bottom(f,u,bnode,data)
	side(f,u,bnode,data)	
end

# ╔═╡ d7960eb2-bfd7-4be4-a0ae-552d68892206
md"""
Factor for ``u_{\text{top}}`` $(@bind utop_fac Slider(range(0.0,2.0,length=41),default=1.2,show_value=true))
"""

# ╔═╡ f39dd714-972c-4d29-bfa8-d2c3795d2eef
function massflow(data, bflux)
	mdot=0.0
	for i=1:data.ng
		mdot += bflux[i] * data.Fluids[i].MW
	end
	mdot/ufac"kg/hr"
end

# ╔═╡ 8db106cc-c5f1-498b-bf8b-fddd8e21b444
function Bern(x)
	y=0
	sat=1e-3
	x=2*x/sat
	if x≈0.0
		y=0
	else
		y=1-x/(exp(x)-1)
		
	end
	sat*y
end

# ╔═╡ 722e681c-225a-4484-b0b8-c85d4536e5f9
function DK_eff(data,T,i)
	ϕ=data.ϕ # porosity
	dp=data.dp # avg. particle size (assumed to be = pore size)
	DK=dp*ϕ^1.5/(1.0-ϕ)*sqrt(8.0*ph"R"*T/(9.0*π*data.Fluids[i].MW))
	Bern(DK)
end

# ╔═╡ b6381008-0280-404c-a86c-9c9c3c9f82eb
function M_matrix(data,T,p,x)
	n=data.ng
	D=data.γ_τ*D_matrix(data,T,p)
	#D=data.D_0*data.ϵ_τ
	M=zeros(eltype(x), n, n)
	for i=1:n
		M[i,i] = -1/DK_eff(data,T,i)
		#M[i,i] = -1/1.0e-5
		for j=1:n
			if j != i
				M[i,i] -= x[j]/D[i,j]
				M[i,j] = x[i]/D[i,j]
			end
		end
	
	end
	M
end

# ╔═╡ ed7941c4-0485-4d84-ad5b-383eb5cae70a
function flux(f,u,edge,data)
	#(;ng, ip, k, Tin, pscale) = data
	ng=data.ng
	ip=data.ip
	iT=data.iT
	k=data.k
	Fluids=data.Fluids
	p=data.p
	X0=data.X0
	
	F=zeros(eltype(u), ng)
	X=zeros(eltype(u), ng)
	
	

	Tbar=0.5*(u[iT,1]+u[iT,2])
	ρf=density_idealgas(Fluids, Tbar, p, X0)
	cf=heatcap_mix(Fluids, Tbar, X0)
	_,λf=dynvisc_thermcond_mix(data, Tbar, X0)
	λbed=kbed(data,λf)*λf

	#T=data.Tin
	T=Tbar

	mole_frac!(data,X,u)
	μ=dynvisc_mix(data, T, X)


	pk,pl = u[ip,1],u[ip,2]
	δp = pk-pl


	# Darcy flow
	ud=k/μ * δp
	conv=ud*ρf*cf/λbed	
	#vh=project(edge,(0,ud))
	#conv=vh*ρf*cf/λbed
	
	Bp,Bm = fbernoulli_pm(conv)
	# temperature flux
	f[iT]= λbed*(Bm*(u[iT,1]-data.Tamb)-Bp*(u[iT,2]-data.Tamb))		

	
	for i=1:ng

		DK = DK_eff(data,T,i)
		#DK = data.D_K_eff[i]
		
		bp,bm=fbernoulli_pm(ud/DK)
		
		
		F[i] = -(bm*u[i,1]-bp*u[i,2])/(ph"R"*T)
	end
	
    
	# computation of fluxes J
	#J = M_matrix(data,X) \ F
	pm = 0.5*(pk+pl)
	J = M_matrix(data,T,pm,X) \ F
	
	f[1:ng] = J
	
	#f[ip] via reaction: ∑pi = p
end

# ╔═╡ 2f24ba85-748a-44a6-bde1-8e320106198d
md"""
## Diffusion Coefficients
"""

# ╔═╡ a96be76c-df66-4149-994c-5b14cc0c3a37
md"""
Estimate binary diffusion coefficients for gases at __low pressures__ (up to 10 bar) according to Fuller
__Fuller EN, Ensley K, Giddings JC (1969)__ Diffusion of halogenated hydrocarbons in helium. J Phys Chem 73:3679.

For pressure up to 10 bar, it is inversely proportional to pressure, and almost independent from the gas phase composition.

The binary diffusion coefficients are then used within the DGM (after correcting for the porous media transport).
"""

# ╔═╡ 3a35ac76-e1b7-458d-90b7-d59ba4f43367
Base.@kwdef mutable struct ModelData <:AbstractModelData
	# number of gas phase species
	ng::Int64		 		= 4
	#ng::Int64		 		= 1
	# names and fluid indices
	gn::Dict{Int, String} 	= Dict(1=>"H2", 2=>"CO2", 3=>"CO", 4=>"H2O")
	#gn::Dict{Int, String} 	= Dict(1=>"N2")
	# inverse names and fluid indices
	gni::Dict{String, Int}  = Dict(value => key for (key, value) in gn)
	# fluids and respective properties in system
	Fluids::Vector{AbstractFluidProps} = [H2,CO2,CO,H2O]
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = [0.5, 0.5, 0.0, 0.0]
	isreactive::Bool = 1

	ip::Int64=ng+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable

	
	RRc::Float64 = 1.0e-5 # reaction rate constant	
		
	α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	
	## irradiation data
	G_lamp::Float64=0.0*ufac"kW/m^2" # solar simulator irradiation flux
	Abs_lamp::Float64=0.7 # avg absorptivity of cat. of irradiation coming from lamp
	Eps_ir::Float64=0.7 # avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
		
	
	## porous filter data
	dp::Float64=50.0*ufac"μm" # average pore size
	#dp::Float64=1.0*ufac"μm" # average pore size
	cath::Float64 = 500.0*ufac"μm" # catalyst layer thickness
	# cylindrical disc / 2D
    D::Float64=12.0*ufac"cm" # disc diameter
	

	# prism / 3D
	wi::Float64=12.0*ufac"cm" # prism width/side lenght
	le::Float64=wi # prism width/side lenght
	h::Float64=0.5*ufac"cm" # frit thickness (applies to 2D & 3D)

	Ac::Float64=pi*D^2.0/4.0*ufac"m^2" # cross-sectional area, circular
	#Ac::Float64=wi^2*ufac"m^2" # cross-sectional area, square
	
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2 	
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
	
	ϕ::Float64=0.36 # porosity, class 2
	# approximation from Wesselingh, J. A., & Krishna, R. (2006). Mass Transfer in Multicomponent Mixtures
	γ_τ::Float64=ϕ^1.5 # constriction/tourtuosity factor
	
	k::Float64=2.9e-11*ufac"m^2" # permeability
	a_s::Float64=0.13*ufac"m^2/g" # specific surface area
	ρfrit::Float64=(1.0-ϕ)*ρs*ufac"kg/m^3" # density of porous frit
	a_v::Float64=a_s*ρfrit # volume specific interface area
	## END porous filter data


	## Flow data
	Qflow::Float64=3400.0*ufac"ml/minute" # volumetric feed flow rate (sccm)
	#Qflow::Float64=34000.0*ufac"ml/minute" # volumetric feed flow rate (sccm)
	# inlet composition

	MWin::Float64 = molarweight_mix(Fluids, X0)
	mdotin::Float64=MWin*Qflow*1.0*ufac"bar"/(ph"R"*273.15*ufac"K")*ufac"kg/s"
	# mass flux in kg/(m^2 s)
	mfluxin::Float64 = mdotin/Ac*ufac"kg/(m^2*s)"
	
	Tin::Float64=298.15*ufac"K" # inlet temperature
	Tamb::Float64=Tin # ambient temperature
	#Tin::Float64=600.0*ufac"K" # inlet temperature
	p::Float64=1.0*ufac"atm" # reactor pressure
	Δp::Float64=0.1*p # pressure differential over reactor
	# u0::Float64=Qflow/(Ac*ϕ)*ufac"m/s" # mean superficial velocity
	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	utop::Float64=u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate
	
end;

# ╔═╡ 7530df59-03e7-4bb6-83f2-86369edc13ee
data=ModelData()

# ╔═╡ 077d4ede-9e0f-4f94-afb0-01bd36c584fc
function catcylinder(;nref=1, r=data.D/2, h=data.h, cath=data.cath)
    	
	hr=r/10.0*2.0^(-nref)
	
	hh=h/10.0*2.0^(-nref)
	hhf=h/50.0*2.0^(-nref)
	
	#RLeft=geomspace(0.0,r/2,hrf, hr)
	#RRight=geomspace(r/2,r,hr,hrf)
	#R=glue(RLeft, RRight)
	Z=geomspace(0,h,hh,hhf)
	
    R=collect(0:hr:r)
    #Z=collect(0:hh:h)
    grid=simplexgrid(R,Z)
    circular_symmetric!(grid)
	#cat_th = 500.0*ufac"μm"
	cellmask!(grid,[0.0,h-cath],[5.0*ufac"cm",h],2) # catalyst layer
	grid
end

# ╔═╡ 8d90a6c3-95c5-4076-a3a1-2c01d119edb9
gridplot(catcylinder(nref=0),legend=:rt,show=true)

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
function main(;data=ModelData())
	#grid=grid_(data,nref=nref)
	#grid=cylinder()
	grid=catcylinder()

	ngas=data.ng
	iT=data.iT
	
	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)
	enable_species!(sys; species=collect(1:(ngas+2))) # gas phase species + p + T
	#enable_species!(sys; species=collect(1:ngas))

	inival=unknowns(sys)
	inival[:,:].=1.0*data.p
	for i=1:ngas
		inival[i,:] .*= data.X0[i]
	end
	inival[iT,:] .= map( (r,z)->(data.Tamb+500*z/data.h),grid)

	
	sol=solve(sys;inival,)

	
	sol,grid,sys
end;

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
sol,grid,sys=main(data=ModelData(utop=ModelData().u0*utop_fac));

# ╔═╡ 33851920-ef27-4041-b67f-face63e1eadb
R=integrate(sys,reaction,sol)

# ╔═╡ 3d660986-f6d7-41a6-800b-68ccd920c7ac
begin
	tf=VoronoiFVM.TestFunctionFactory(sys);
	#testfunction(factory::TestFunctionFactory, bc0, bc1) -> Any

	# bottom - inflow
	Tbot=testfunction(tf,[2,3,4],1)
	Ibot=integrate(sys,Tbot,sol)
	# top - outflow
	Ttop=testfunction(tf,[1,2,4],3)
	Itop=integrate(sys,Ttop,sol)
end;

# ╔═╡ bd7552d2-2c31-4834-97d9-ccdb4652242f
function plane(data,sol)
	#grid=cylinder()
	grid=catcylinder()
	
	bfacemask!(grid, [data.D/4,0],[data.D/4,data.h],5)

	# transform y coordinate of parent grid into x coordinate of subgrid
	function _2to1(a,b)
		a[1]=b[2]
	end
	grid_1D  = subgrid(grid, [5], boundary=true, transform=_2to1) 

	sol_p = []
	for i=1:(data.ng+1)
		sol_i = view(sol[i, :], grid_1D)
		push!(sol_p, collect(sol_i))
	end
		
	sol_p, grid_1D	
	#sol_cutplane, grid_2D	
end

# ╔═╡ fbd53732-2fef-490f-aef0-cac0e7d60574
function cylinder(;nref=2, r=data.D/2, h=data.h)
    #step=0.1*ufac"cm"*2.0^(-nref)

	
	hr=r/10.0*2.0^(-nref)
	#hrf=r/100.0*2.0^(-nref)
	hh=h/10.0*2.0^(-nref)

	#RLeft=geomspace(0.0,r/2,hrf, hr)
	#RRight=geomspace(r/2,r,hr,hrf)
	#R=glue(RLeft, RRight)
	
    R=collect(0:hr:r)
    Z=collect(0:hh:h)
    grid=simplexgrid(R,Z)
    circular_symmetric!(grid)
	grid
end

# ╔═╡ 8f2b1d6a-04eb-4a79-9807-99905c8ef398
gridplot(cylinder(nref=0),legend=:rt)

# ╔═╡ 1d8bcb9b-cff4-4868-8dfc-ba2900084645
data.mdotin*3600

# ╔═╡ 8b1a0902-2542-40ed-9f91-447bffa4290f
md"""
Mass flows:

- through __bottom__ boundary __$(round(massflow(data, Ibot),sigdigits=3))__ kg/h
- through __top__ boundary __$(round(massflow(data, Itop),sigdigits=3))__ kg/h
"""

# ╔═╡ a2dd4745-91bf-4f89-aa85-b57ed889f50a
md"""
Flow of $(data.gn[1]):

- through __bottom__ boundary __$(round(Ibot[1],sigdigits=2))__ mol/s
- through __top__ boundary __$(round(Itop[1],sigdigits=2))__ mol/s

Flow of $(data.gn[2]):

- through __bottom__ boundary __$(round(Ibot[2],sigdigits=2))__ mol/s
- through __top__ boundary __$(round(Itop[2],sigdigits=2))__ mol/s

Flow of $(data.gn[3]):

- through __bottom__ boundary __$(round(Ibot[3],sigdigits=2))__ mol/s
- through __top__ boundary __$(round(Itop[3],sigdigits=2))__ mol/s

Flow of $(data.gn[4]):

- through __bottom__ boundary __$(round(Ibot[4],sigdigits=2))__ mol/s
- through __top__ boundary __$(round(Itop[4],sigdigits=2))__ mol/s
"""

# ╔═╡ 3bd80c19-0b49-43f6-9daa-0c87c2ea8093
let
	iT=data.iT
	vis=GridVisualizer()
	scalarplot!(vis, grid, sol[iT,:].- 273.15, resolution=(600,400), show=true)
end

# ╔═╡ 4b41d985-8ebc-4cab-a089-756fce0d3060
let	
	ngas=data.ng
	ip=data.ip
	# visualization
	vis=GridVisualizer(layout=(2,1),resolution=(600,800))

	# 2D plot
	scalarplot!(vis[1,1], grid, sol[1,:], zoom=1.5)
	
	# 1D plot, over reactor height
	c1 = colorant"red"
	c2 = colorant"blue"
	cols=range(c1, stop=c2, length=ngas+1)

	sol_1D, grid_1D=plane(data,sol)

	for i in 1:ngas
		scalarplot!(vis[2,1], grid_1D, sol_1D[i], color=cols[i], label="$(data.gn[i])", clear=false, legend=:ct)
	end
	

	scalarplot!(vis[2,1], grid_1D, sol_1D[ip], color=cols[ip], label="total p", clear=false, legend=:ct)

	reveal(vis)

end

# ╔═╡ Cell order:
# ╠═11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╟─863c9da7-ef45-49ad-80d0-3594eca4a189
# ╠═077d4ede-9e0f-4f94-afb0-01bd36c584fc
# ╠═8d90a6c3-95c5-4076-a3a1-2c01d119edb9
# ╠═fbd53732-2fef-490f-aef0-cac0e7d60574
# ╠═8f2b1d6a-04eb-4a79-9807-99905c8ef398
# ╠═b659879f-3706-4937-b655-116f1be85069
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─f4dcde90-6d8f-4b17-b4ec-367d2372637f
# ╟─3703afb0-93c4-4664-affe-b723758fb56b
# ╟─c4e08410-3d1a-4921-aa5c-2f232daaa07a
# ╟─313cd88a-1497-4d3a-b09a-9b98a6dad9c2
# ╟─21d0195b-b170-460d-989e-f9d00b511237
# ╟─8f4843c6-8d2b-4e24-b6f8-4eaf3dfc9bf0
# ╟─66b55f6b-1af5-438d-aaa8-fe4745e85426
# ╟─8528e15f-cce7-44d7-ac17-432f92cc5f53
# ╠═ed7941c4-0485-4d84-ad5b-383eb5cae70a
# ╟─a6afe118-dcbd-4126-8646-c7268acfacf3
# ╟─a60ce05e-8d92-4172-b4c1-ac3221c54fe5
# ╟─24374b7a-ce77-45f0-a7a0-c47a224a0b06
# ╟─4865804f-d385-4a1a-9953-5ac66ea50057
# ╠═722e681c-225a-4484-b0b8-c85d4536e5f9
# ╟─2fb1a154-2721-4da6-848e-99d46dfce774
# ╟─cc12f23a-4a80-463b-baf8-22f58d341d9d
# ╠═91d8119f-501f-49f2-93e0-88da8d996f7a
# ╠═4af2237c-9144-4ffc-8966-2b3bf9d3d720
# ╠═2191bece-e186-4d8e-8a21-3830441baf11
# ╠═b6381008-0280-404c-a86c-9c9c3c9f82eb
# ╟─02b76cda-ffae-4243-ab40-8d0fe1325776
# ╠═78cf4646-c373-4688-b1ac-92ed5f922e3c
# ╟─906ad096-4f0c-4640-ad3e-9632261902e3
# ╠═7da59e27-62b9-4b89-b315-d88a4fd34f56
# ╠═40906795-a4dd-4e4a-a62e-91b4639a48fa
# ╠═edd9fdd1-a9c4-4f45-8f63-9681717d417f
# ╠═29d66705-3d9f-40b1-866d-dd3392a1a268
# ╠═333b5c80-259d-47aa-a441-ee7894d6c407
# ╠═aa498412-e970-45f2-8b11-249cc5c2b18d
# ╟─d7960eb2-bfd7-4be4-a0ae-552d68892206
# ╠═1d8bcb9b-cff4-4868-8dfc-ba2900084645
# ╠═f39dd714-972c-4d29-bfa8-d2c3795d2eef
# ╟─8b1a0902-2542-40ed-9f91-447bffa4290f
# ╟─a2dd4745-91bf-4f89-aa85-b57ed889f50a
# ╠═33851920-ef27-4041-b67f-face63e1eadb
# ╠═3d660986-f6d7-41a6-800b-68ccd920c7ac
# ╠═3bd80c19-0b49-43f6-9daa-0c87c2ea8093
# ╠═4b41d985-8ebc-4cab-a089-756fce0d3060
# ╟─8db106cc-c5f1-498b-bf8b-fddd8e21b444
# ╟─bd7552d2-2c31-4834-97d9-ccdb4652242f
# ╟─2f24ba85-748a-44a6-bde1-8e320106198d
# ╟─a96be76c-df66-4149-994c-5b14cc0c3a37
# ╠═3a35ac76-e1b7-458d-90b7-d59ba4f43367
# ╠═7530df59-03e7-4bb6-83f2-86369edc13ee
