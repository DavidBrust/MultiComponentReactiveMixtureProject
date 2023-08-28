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
	using PyPlot	
	using Colors

	using ForwardDiff, DiffResults, Preferences

	using Revise
	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
	set_preferences!(ForwardDiff, "nansafe_mode" => true)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="Sensitivity")

# ╔═╡ 2ed3223e-a604-410e-93d4-016580f49093
md"""
# Domain / Grid
"""

# ╔═╡ f9205777-f2af-4192-88f5-a72e193f13df
md"""
## 3D
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

# ╔═╡ 2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
md"""
# Governing system
"""

# ╔═╡ f4dcde90-6d8f-4b17-b4ec-367d2372637f
md"""
#### Species Mass balances
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

# ╔═╡ 722e681c-225a-4484-b0b8-c85d4536e5f9
function DK_eff(data,T,i)
	ϕ=data.ϕ # porosity
	dp=data.dp # avg. particle size (assumed to be = pore size)
	DK=dp*ϕ^1.5/(1.0-ϕ)*sqrt(8.0*ph"R"*T/(9.0*π*data.Fluids[i].MW))
	#Bern(DK)
end

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

# ╔═╡ 2191bece-e186-4d8e-8a21-3830441baf11
function D_matrix(data, T, p)
	n=data.ng
	v = zeros(typeof(T),n,n)
	
	for i=1:(n-1)
		for j=(i+1):n
			v[j,i] = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
		end
	end
	Symmetric(v, :L)
end

# ╔═╡ b6381008-0280-404c-a86c-9c9c3c9f82eb
function M_matrix(data,T,p,x)
	n=data.ng
	D=data.γ_τ*D_matrix(data,T,p)
	M=zeros(eltype(x), n, n)
	for i=1:n
		M[i,i] = -1/DK_eff(data,T,i)
		for j=1:n
			if j != i
				M[i,i] -= x[j]/D[i,j]
				M[i,j] = x[i]/D[i,j]
			end
		end	
	end
	M
end

# ╔═╡ 02b76cda-ffae-4243-ab40-8d0fe1325776
md"""
##### Auxiliary functions
"""

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	ngas=data.ng
	ip=data.ip
	
	if node.region == 2 && data.isreactive # catalyst layer
		#ng=data.gni # gas species indices from names
		nr=data.kinpar.rni # reaction indices from names
		iT=data.iT
		
		pi = u[1:ngas]./ufac"bar"
		# negative sign: sign convention of VoronoiFVM: source term < 0 
		# ri returns reaction rate in mol/(h gcat)
		RR = -data.mcats*ri(data.kinpar,u[iT],pi)*ufac"mol/hr"*ufac"1/g"
		# reactions in S3P kinetics model
		# R1: CO + H2O = CO2 + H2
		# R2: CH4 + 2 H2O = CO2 + 4 H2
		# R3: CH4 + H2O = CO + 3 H2

		for i=1:ngas
			f[i] = sum(data.kinpar.nuij[i,:] .* RR)
		end
		
		# temperature eq. / heat source
		ΔHi=data.kinpar.ΔHi
		f[iT] = -(RR[nr[:R1]]*ΔHi[:R1]+RR[nr[:R2]]*ΔHi[:R2] +RR[nr[:R3]]*ΔHi[:R3])
	end
	
	# ∑xi = 1
	f[ip]=u[ip]-sum(u[1:ngas])
	
end

# ╔═╡ 906ad096-4f0c-4640-ad3e-9632261902e3
md"""
## Boundary Conditions
"""

# ╔═╡ 40906795-a4dd-4e4a-a62e-91b4639a48fa
function bottom(f,u,bnode,data)
	if bnode.region==Γ_bottom # bottom boundary
		iT=data.iT
		f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)

		ng=data.ng
		
		for i=1:data.ng
			# specify flux at boundary: flow velocity is normal to bot boundary	
			f[i] = -data.u0*data.X0[i]*data.pn/(ph"R"*data.Tn)
			
		end
	end
end

# ╔═╡ edd9fdd1-a9c4-4f45-8f63-9681717d417f
function side(f,u,bnode,data)
	# side wall boundary condition
	iT=data.iT
	boundary_robin!(f,u,bnode;species=iT,region=[Γ_side_back,Γ_side_right], factor=data.α_w, value=data.Tamb*data.α_w)	
end

# ╔═╡ c4521a0c-c5af-43cd-97bc-a4a7a42d27b1
md"""
## Metrics of Reactor Operation
1. Catalyst layer average temperature
1. CO Yield (Converions x Selectivity)
1. Solar-to-chemical efficiency
"""

# ╔═╡ a6e61592-7958-4094-8614-e77446eb2223
md"""
### Surface Temperature
"""

# ╔═╡ 808e5a71-572f-4b0c-aeb3-9513d969dada
function CatDims(grid)
	Cells = findall(x->x == 2, grid[CellRegions]) # catalyst layer
	Nodes = grid[CellNodes][:,Cells]
	Nodes = unique(reshape(Nodes,1,:))
	coords = grid[Coordinates][:,Nodes]

	mincatwi, maxcatwi = minimum(coords[1,:]), maximum(coords[1,:])
	mincatle, maxcatle = minimum(coords[2,:]), maximum(coords[2,:])
	mincathe, maxcathe = minimum(coords[3,:]), maximum(coords[3,:])

	catA = (maxcatwi-mincatwi)*(maxcatle-mincatle)
	catV = (maxcatwi-mincatwi)*(maxcatle-mincatle)*(maxcathe-mincathe)
	catA, catV
end

# ╔═╡ 68e2628a-056a-4ec3-827f-2654f49917d9
function Tcatavg(sol,sys,grid,data)

	#sol_1D, grid_1D=planeTop(data,sol)

	function Tcat_(f,u,bnode,data)
		iT = data.iT		
		f[iT] = u[iT]		
	end

	catA, catV = CatDims(grid)
	
	
	Tcat_avg_surf=integrate(sys,Tcat_,sol; boundary=true)[data.iT,Γ_top_cat] / catA - 273.15
	Tcat_avg_vol=integrate(sys,Tcat_,sol; )[data.iT,2] / catV - 273.15
	#Tcat_avg_surf, Tcat_avg_vol
	Tcat_avg_surf
end

# ╔═╡ 4cde8752-bbdf-4b83-869e-46b78bb4adb5
md"""
### Yield of CO
"""

# ╔═╡ 7bf4d925-56e1-4c85-a40a-fef1917f501e
md"""
```math
Y_{\text{CO}} = \frac{\dot n_{\text{CO}}}{\dot n^0_{\text{CO}_2}}
```
"""

# ╔═╡ aebcb161-425c-4a46-aaed-d4b13b3e6654
md"""
### Solar to Chemical Efficiency (STC)
"""

# ╔═╡ 6bbd0496-c275-4cd5-bf48-9ea4e77a9091
md"""
Defined with the reaction enthalpy of the RWGS reaction at standart conditions.
"""

# ╔═╡ 1f3b2a0e-c76f-462c-a85c-e85815fd7e5a
md"""
```math
	\text{STC} = \frac{\dot n_{\text{CO}} \Delta H^0_{\text{RWGS}}}{I_{\text{lamp}} A}
```
"""

# ╔═╡ b34c1d1b-a5b8-4de9-bea9-3f2d0503c1c0
md"""
### Mass and Molar flows
"""

# ╔═╡ f39dd714-972c-4d29-bfa8-d2c3795d2eef
function massflow(data, bflux)
	mdot=0.0
	for i=1:data.ng
		mdot += bflux[i] * data.Fluids[i].MW
	end
	mdot/ufac"kg/hr"
end

# ╔═╡ 91d8119f-501f-49f2-93e0-88da8d996f7a
function vel_darcy(data, ∇p, T, X)
	μ=dynvisc_mix(data, T, X)
	data.k/μ * ∇p # m/s
end

# ╔═╡ b2df1087-6628-4889-8cd6-c5ee7629cd93
md"""
### Temperature Plot
"""

# ╔═╡ 2790b550-3105-4fc0-9070-d142c19678db
md"""
### Partial Pressure Plots
"""

# ╔═╡ e25e7b7b-47b3-457c-995b-b2ee4a87710a
md"""
## Model Data
"""

# ╔═╡ 3a35ac76-e1b7-458d-90b7-d59ba4f43367
#Base.@kwdef mutable struct ModelData{Tv} <:AbstractModelData
Base.@kwdef mutable struct ModelData <:AbstractModelData
	#S::Tv = 1.0
	
	# catalyst / chemistry data
	# kinetic parameters, S3P="simple 3 parameter" kinetics fit to UPV lab scale experimental data
	# kinpar::AbstractKineticsData = S3P
	kinpar::AbstractKineticsData = XuFroment1989
	
	# number of gas phase species
	ng::Int64		 		= S3P.ng
	# names and fluid indices
	gn::Dict{Int, Symbol} 	= S3P.gn
		# inverse names and fluid indices
	gni::Dict{Symbol, Int}  = S3P.gni
	# fluids and respective properties in system
	Fluids::Vector{FluidProps} = S3P.Fluids
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = let
		x=zeros(Float64, ng)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		x/sum(x)
	end # inlet composition
	#X0::Vector{Float64} = [1.0]
	
	# volume specific cat mass loading, UPV lab scale PC reactor
	#mcats::Float64 =1234.568*ufac"kg/m^3"
	mcats::Float64 =20.0*ufac"kg/m^3"
	isreactive::Bool = 1
	#isreactive::Bool = 0

	ip::Int64=ng+1 # index of total pressure variable
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
	#dp::Float64=25.0*ufac"μm" # average pore size
	

	# frit thickness (applies to 2D & 3D)
	h::Float64=0.5*ufac"cm"
	# catalyst layer thickness (applies to 2D & 3D)
	cath::Float64 = 1000.0*ufac"μm"
	
	# cylindrical disc / 2D
    D::Float64=12.0*ufac"cm" # disc diameter
	catD::Float64 = 10.0*ufac"cm" # catalyst layer diameter
	

	# prism / 3D
	wi::Float64=12.0*ufac"cm" # prism width/side lenght
	le::Float64=wi # prism width/side lenght
	catwi::Float64=10.0*ufac"cm" # prism width/side lenght
	
	

	#Ac::Float64=pi*D^2.0/4.0*ufac"m^2" # cross-sectional area, circular
	Ac::Float64=wi*le*ufac"m^2" # cross-sectional area, square
	
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
	a_v::Float64=a_s*ρfrit # volume specific interface area
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
	#Tin::Float64=600.0*ufac"K" # inlet temperature
	p::Float64=1.0*ufac"atm" # reactor pressure

	# u0::Float64=Qflow/(Ac*ϕ)*ufac"m/s" # mean superficial velocity
	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	utop::Float64=u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate
	
end;

# ╔═╡ e73a3dfb-740e-4a20-9101-47cf18bcb9be
data=ModelData()

# ╔═╡ 077d4ede-9e0f-4f94-afb0-01bd36c584fc
function catcylinder(;nref=0, r=data.D/2, h=data.h, catr=data.catD/2,cath=data.cath)
    	
	hr=r/10.0*2.0^(-nref)
	
	hh=h/10.0*2.0^(-nref)
	#hhf=h/20.0*2.0^(-nref)
	
	
	#Z=geomspace(0,h,hh,hhf)
	
    R=collect(0:hr:r)
	Z=collect(0:hh:h)
    
    grid=simplexgrid(R,Z)
    circular_symmetric!(grid)
	cellmask!(grid,[0.0,h-cath],[catr,h],2) # catalyst layer region
	
	bfacemask!(grid,[0.0,h],[catr,h],5) # catalyst layer boundary
	grid
end

# ╔═╡ 8d90a6c3-95c5-4076-a3a1-2c01d119edb9
let
	vis=GridVisualizer(Plotter=PyPlot, resolution=(600,400))
	gridplot!(vis, catcylinder(nref=0,cath=2*ufac"mm"),legend=:none,linewidth=0.1,)
	reveal(vis)
	#save("../img/out/domain.svg", vis)
end

# ╔═╡ 985718e8-7ed7-4c5a-aa13-29462e52d709
md"""
Cutplane at ``z=`` $(@bind zcut Slider(range(0.0,data.h,length=101),default=data.h,show_value=true)) m
"""

# ╔═╡ ada45d4d-adfa-484d-9d0e-d3e7febeb3ef
function prism_sq(;nref=0, w=data.wi, h=data.h, cath=data.cath, catwi=data.catwi)
	
	hw=w/2.0/5.0*2.0^(-nref)
	hh=h/5.0*2.0^(-nref)
	W=collect(0:hw:(w/2.0))
    H=collect(0:hh:h)
	grid=simplexgrid(W,W,H)
	
	# catalyst layer region
	cellmask!(grid,[0.0,0.0,h-cath],[catwi/2,catwi/2,h],2)
	# catalyst layer boundary
	bfacemask!(grid,[0.0,0.0,h],[catwi/2,catwi/2,h],Γ_top_cat)	
end

# ╔═╡ ff58b0b8-2519-430e-8343-af9a5adcb135
let
	vis=GridVisualizer(resolution=(600,400), zoom=1.9, )
	#zcut = data.h-data.cath/2
	gridplot!(vis, prism_sq(nref=0,cath=2*ufac"mm"), zplane=zcut)
	reveal(vis)
	#save("../img/out/domain.svg", vis)
end

# ╔═╡ fec9ca6d-d815-4b50-bec7-f8fb3d8195ba
function TopPlane(data,sol)

	grid=prism_sq()
	wi=data.wi/2
	
	bid = maximum(grid[BFaceRegions])+1
	bfacemask!(grid, [0,0,data.h],[wi,wi,data.h],bid)

	# keep x-y coordinates of parent grid
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[2]
	end
	#grid_1D  = subgrid(grid, [bid], boundary=true, transform=_3to1) 
	grid_2D  = subgrid(grid, [bid], boundary=true, transform=_3to2) 

	sol_p = []
	for i=1:(data.ng+2)
		sol_i = view(sol[i, :], grid_2D)
		push!(sol_p, collect(sol_i))
	end

	sol_p, grid_2D
end

# ╔═╡ bd7552d2-2c31-4834-97d9-ccdb4652242f
function SolAlongLine(data,sol)
		
	grid=prism_sq()
	mid_xy=data.wi/5*2

	Nodes = findall(x->x[1] == mid_xy && x[2] == mid_xy, eachcol(grid[Coordinates]))
	
	grid1D = grid[Coordinates][:,Nodes]
	grid1D = grid1D[3,:] # extract z-coordinate


	sol_p = []
	for i=1:(data.ng+2)
		push!(sol_p, sol[i,Nodes])
	end

	sol_p,grid1D
end

# ╔═╡ a9ec34cb-5c01-423f-a9b1-027cba13dca3
data.mdotin/ufac"kg/hr"

# ╔═╡ 6258c7c7-dbf5-4d5b-a8e3-de2250f5e66a
md"""
# Sensitivity
"""

# ╔═╡ 2abbe557-9869-4861-9c3b-208f1f2edf59
@bind RunSensitivity CheckBox()

# ╔═╡ ec3dce96-a6df-4ff5-aa7b-ee8631555b67
SensPar = :Abs_lamp_cat
#SensPar = :mcats

# ╔═╡ cd2547b4-0e2f-4618-8cca-3e02ab7cb236
md"""
### Parametric Model Data
"""

# ╔═╡ 1046820c-777f-4180-a570-1efe90924ec9
md"""
Use this data structure for sensitivity calculation.
"""

# ╔═╡ 2551d24d-91fe-40cd-8ade-eaf6b549efa3
Base.@kwdef mutable struct ModelDataSens{Tv} <:AbstractModelData
	
	##############################################################################
	# BEGIN SENSITIVITY PARAMETERS
	S::Tv = 1.0
	# cat avg absorptivity of irradiation from lamp
	Abs_lamp_cat::Tv=S
	#Abs_lamp_cat::Float64=0.7 # cat avg absorptivity of irradiation from lamp
	
	# volume specific cat mass loading, UPV lab scale PC reactor
	mcats::Float64 =1234.568*ufac"kg/m^3"
	#mcats::Tv=S
	
	# frit avg absorptivity of irradiation from lamp
	#Abs_lamp_frit::Tv=S
	Abs_lamp_frit::Float64=0.7
	
	# avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
	#Eps_ir::Tv=S
	Eps_ir::Float64=0.7

	λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2 	
	#λs::Tv=S

	# END SENSITIVITY PARAMETERS
	##############################################################################
	
	# catalyst / chemistry data
	# kinetic parameters, S3P="simple 3 parameter" kinetics fit to UPV lab scale experimental data
	kinpar::AbstractKineticsData = XuFroment1989
	
	# number of gas phase species
	ng::Int64		 		= S3P.ng
	# names and fluid indices
	gn::Dict{Int, Symbol} 	= S3P.gn
	# inverse names and fluid indices
	gni::Dict{Symbol, Int}  = S3P.gni
	# fluids and respective properties in system
	Fluids::Vector{FluidProps} = S3P.Fluids
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = let
		x=zeros(Float64, ng)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		x/sum(x)
	end # inlet composition
	#X0::Vector{Float64} = [1.0]
	
	
	
	isreactive::Bool = 1
	#isreactive::Bool = 0

	ip::Int64=ng+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable

		
	α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	
	## irradiation data
	Tau_quartz::Float64=0.9 # transmission coefficient of quartz window
	G_lamp::Float64=1.0*ufac"kW/m^2" # solar simulator irradiation flux

	## porous filter data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0

	

	# frit thickness (applies to 2D & 3D)
	h::Float64=0.5*ufac"cm"
	# catalyst layer thickness (applies to 2D & 3D)
	cath::Float64 = 1000.0*ufac"μm"
	
	# cylindrical disc / 2D
    D::Float64=12.0*ufac"cm" # disc diameter
	catD::Float64 = 10.0*ufac"cm" # catalyst layer diameter
	

	# prism / 3D
	wi::Float64=12.0*ufac"cm" # prism width/side lenght
	le::Float64=wi # prism width/side lenght
	catwi::Float64=10.0*ufac"cm" # prism width/side lenght
	
	

	#Ac::Float64=pi*D^2.0/4.0*ufac"m^2" # cross-sectional area, circular
	Ac::Float64=wi*le*ufac"m^2" # cross-sectional area, square
	
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	
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
	a_v::Float64=a_s*ρfrit # volume specific interface area
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
	#Tin::Float64=600.0*ufac"K" # inlet temperature
	p::Float64=1.0*ufac"atm" # reactor pressure

	# u0::Float64=Qflow/(Ac*ϕ)*ufac"m/s" # mean superficial velocity
	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	utop::Float64=u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate
	
end;

# ╔═╡ 4af2237c-9144-4ffc-8966-2b3bf9d3d720
function mole_frac!(node,data,X,u)
	n=data.ng
	sump = 0.0
	for i=1:n
		if u isa VoronoiFVM.EdgeUnknowns
			X[i] = 0.5*(u[i,1]+u[i,2])
		elseif u isa VoronoiFVM.BNodeUnknowns # boundary node
			if data isa ModelDataSens
				X[i] = u[i].value
			else
				X[i] = u[i]
			end
		end
		sump += X[i]
	end
	X .= X / sump
	nothing
end

# ╔═╡ ed7941c4-0485-4d84-ad5b-383eb5cae70a
function flux(f,u,edge,data)
	#(;ng, ip, k, Tin, pscale) = data
	ng=data.ng
	ip=data.ip
	iT=data.iT
	k=data.k
	Fluids=data.Fluids
	
	F=zeros(eltype(u), ng)
	X=zeros(eltype(u), ng)

	pk,pl = u[ip,1],u[ip,2]
	δp = pk-pl
	pm = 0.5*(pk+pl)

	T=0.5*(u[iT,1]+u[iT,2])
	mole_frac!(edge,data,X,u)
	
	μ=dynvisc_mix(data, T, X)
	cf=heatcap_mix(Fluids, T, X)
	_,λf=dynvisc_thermcond_mix(data, T, X)
	λbed=kbed(data,λf)*λf

	# Darcy flow
	ud=-k/μ * δp
	# vh=project(edge,(0,ud)) # 2D
	vh=project(edge,(0,0,ud)) # 3D
	# convective enthalpy flux
	conv=vh*pm/(ph"R"*T)*cf/λbed
	
	Bp,Bm = fbernoulli_pm(conv)
	# temperature flux
	f[iT]= λbed*(Bm*(u[iT,1]-data.Tamb)-Bp*(u[iT,2]-data.Tamb))		

	
	for i=1:ng
		DK = DK_eff(data,T,i)
		bp,bm=fbernoulli_pm(vh/DK)		
		F[i] = -(bm*u[i,1]-bp*u[i,2])/(ph"R"*T)
	end
	
    
	# computation of fluxes J
	J = M_matrix(data,T,pm,X) \ F
	
	f[1:ng] = J	
	#f[ip] via reaction: ∑pi = p
end

# ╔═╡ 7da59e27-62b9-4b89-b315-d88a4fd34f56
function top(f,u,bnode,data)
	# top boundaries (cat layer & frit)
	if bnode.region==Γ_top_frit || bnode.region==Γ_top_cat 
		ng=data.ng
		ip = data.ip
		iT = data.iT
		flux_rerad = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)

		
		X=zeros(eltype(u), ng)
		mole_frac!(bnode,data,X,u)

		cf=heatcap_mix(data.Fluids, u[iT], X)
		
		flux_convec=data.utop*u[ip]/(ph"R"*u[iT])*cf*(u[iT]-data.Tamb)

		abs = 0.0
		if bnode.region==Γ_top_frit
			abs=data.Abs_lamp_frit
		else # catalyst layer
			abs=data.Abs_lamp_cat
		end
		f[iT] = -abs*data.Tau_quartz*data.G_lamp + flux_rerad + flux_convec
		
		
		# flow velocity is normal to top boundary
		for i=1:data.ng
			f[i] = data.utop*u[i]/(ph"R"*u[iT])
		end
	end
end

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
	ip=data.ip
	p0=data.p
	ng=data.ng
	X0=data.X0
	
	# set partial pressures of inlet comp at bottom boundary
	#for i=1:ng
	#	boundary_dirichlet!(f,u,bnode,i,1,X0[i]*p0)
	#end
	#boundary_dirichlet!(f,u,bnode,ip,1,p0)


	top(f,u,bnode,data)
	bottom(f,u,bnode,data)
	side(f,u,bnode,data)	
end

# ╔═╡ 4042a8de-8df8-4c9d-b1c9-6323afe20d95
function runSensitivity(grid, sol, data, Tv)
	ngas=data.ng
	
	sys=VoronoiFVM.System( 	grid;
							valuetype = Tv,
							data=data,
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)
	enable_species!(sys; species=collect(1:(ngas+2))) # gas phase species + p + T
	inival=unknowns(sys)
	inival[:,:] .= sol
	solSens = solve(sys;inival,)
	solSens,sys
end

# ╔═╡ 6ae6d894-4923-4408-9b77-1067ba9e2aff
function MoleFlows(sol,sys,data)
	# bottom - inflow
	Ibot=integrate(sys,bottom,sol; boundary=true)[:,Γ_bottom]
	# top: outflow
	# uncoated outer frit area, inner cat coated area
	Itop=integrate(sys,top,sol; boundary=true)[:,[Γ_top_frit,Γ_top_cat]] 
	Itop=sum(Itop, dims=2)
	Ibot[1:data.ng],Itop[1:data.ng]
end

# ╔═╡ b3bf7b7d-eb38-4a32-87f7-5aef098ad03e
function Yield_CO(sol,sys,data)
	ndot_bot,ndot_top = MoleFlows(sol,sys,data)
	
	ndot_top[data.gni[:CO]] / abs(ndot_bot[data.gni[:CO2]])
end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
function main(;data=ModelData())
	#grid=grid_(data,nref=nref)
	#grid=cylinder()
	#grid=catcylinder()
	grid=prism_sq()

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
	# inival[iT,:] .= map( (r,z)->(data.Tamb+25.0*z/data.h),grid)
	inival[iT,:] .= data.Tamb

	sol=solve(sys;inival,)
	
	 function pre(sol,par)
		(;ng,iT)=data
	 	
	 	# iteratively adapt top outflow boundary condition
	 	function Inttop(f,u,bnode,data)
			
	 		X=zeros(eltype(u), ng)
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
		
	# 	# embedding parameter: Qflow / ml/minute
	# 	#Qflows = [5000.0,50000.0]*ufac"ml/minute"
	# 	#Qflow = minimum(Qflows) + par*(maximum(Qflows)-minimum(Qflows))
	# 	#data.Qflow=Qflow
	# 	#data.u0=Qflow/(data.Ac)*ufac"m/s" # mean superficial velocity
	# 	#data.utop=data.u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate

		
	 	# specific catalyst loading
	 	#mcats=[10.0, 1300.0]*ufac"kg/m^3"
		 mcats=[10.0, 20.0]*ufac"kg/m^3"
	 	data.mcats= minimum(mcats) + par*(maximum(mcats)-minimum(mcats))

	 	# irradiation flux density
	 	G_lamp=[1.0, 100.0]*ufac"kW/m^2"
	 	data.G_lamp= minimum(G_lamp) + par*(maximum(G_lamp)-minimum(G_lamp))
	 end
	
	 control=SolverControl(;
	 				  handle_exceptions=true,
	 				  Δp_min=1.0e-4,					  
	 				  Δp=0.1,
	 				  Δp_grow=1.2,
	 				  Δu_opt=100000.0, # large value, due to unit Pa of pressure?
	 				  )
	
	# #sol=solve(sys;inival,)
	
	 sol=solve(sys;inival,embed=[0.0,1.0],pre,control)

	
	sol,grid,sys,data
end;

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
begin
	sol_,grid,sys,data_embed=main(data=ModelData());
	if sol_ isa VoronoiFVM.TransientSolution
		sol = copy(sol_(sol_.t[end]))
	else
		sol = copy(sol_)
	end
end;

# ╔═╡ 2739bcff-3fb0-4169-8a1a-2b0a14998cec
md"""
Average Catalyst layer surface temperature: $(round(Tcatavg(sol,sys,grid,data_embed)[1])) °C
"""

# ╔═╡ 30393c90-298c-412d-86ce-e36106613d35
let
	sol_, grid_ =TopPlane(data,sol)

	vis=GridVisualizer(resolution=(600,400), Plotter=PyPlot)
	scalarplot!(vis, grid_, sol_[data.iT].- 273.15, show=true)
	#save("../img/out/Test.svg", vis)
end

# ╔═╡ 0f7cce89-3add-4316-bcc2-a924064af884
md"""
Yield of CO relative to inflow of ``\text{CO}_2``: $(round(Yield_CO(sol,sys,data_embed),sigdigits=2))
"""

# ╔═╡ 47161886-9a5c-41ac-abf5-bbea82096d5a
function STCefficiency(sol,sys,data)
	_,ndot_top = MoleFlows(sol,sys,data)
	# R2 = RWGS in Xu & Froment kinetic model
	STC_Atot = ndot_top[data.gni[:CO]] * -data.kinpar.ΔHi[:R2] / (data.G_lamp*data.Ac/4)
	catA, _ = CatDims(grid)
	STC_Acat = ndot_top[data.gni[:CO]] * -data.kinpar.ΔHi[:R2] / (data.G_lamp*catA)
	STC_Atot, STC_Acat
end

# ╔═╡ 428afb22-55ab-4f64-805f-7e15ad4cf23f
md"""
Solar-to-chemical efficiency as defined above: $(round(STCefficiency(sol,sys,data_embed)[1]*100, sigdigits=2)) %
"""

# ╔═╡ 7ab38bc2-9ca4-4206-a4c3-5fed673557f1
begin
	ndot_bot, ndot_top = MoleFlows(sol,sys,data_embed)
end;

# ╔═╡ 8b1a0902-2542-40ed-9f91-447bffa4290f
md"""
Mass flows:

- through __bottom__ boundary __$(round(massflow(data, ndot_bot),sigdigits=3))__ kg/h
- through __top__ boundary __$(round(massflow(data, ndot_top),sigdigits=3))__ kg/h
"""

# ╔═╡ b06a7955-6c91-444f-9bf3-72cfb4a011ec
md"""

Chemical species flows through porous frit from bottom and top:

|    | Bottom     | Top     |    |
|----|-------|-------|-------|
| $(data.gn[1]) | $(abs(round(ndot_bot[1]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_top[1]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[2]) | $(abs(round(ndot_bot[2]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_top[2]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[3]) | $(abs(round(ndot_bot[3]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_top[3]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[4]) | $(abs(round(ndot_bot[4]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_top[4]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[5]) | $(abs(round(ndot_bot[5]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_top[5]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[6]) | $(abs(round(ndot_bot[6]/ufac"mol/hr")))  |   $(round(ndot_top[6]/ufac"mol/hr"))   |  mol/hr  |


"""

# ╔═╡ 3bd80c19-0b49-43f6-9daa-0c87c2ea8093
let
	iT=data.iT
	#vis=GridVisualizer(Plotter=PyPlot)
	vis=GridVisualizer()
	scalarplot!(vis, grid, sol[iT,:].- 273.15, show=true)

	#reveal(vis)
	#save("../img/out/Temeprature.svg", vis)
end

# ╔═╡ bea97fb3-9854-411c-8363-15cbef13d033
let
	sol_, grid_ = SolAlongLine(data,sol)
	
	ngas=data.ng
	c1 = colorant"red"
	c2 = colorant"blue"
	cols=range(c1, stop=c2, length=ngas+1)

	vis=GridVisualizer()
	#vis=GridVisualizer(Plotter=PyPlot)
	for i in 1:ngas
		scalarplot!(vis, grid_, sol_[i], color=cols[i], label="$(data.gn[i])", clear=false, legend=:best)
	end

	scalarplot!(vis, grid_, sol_[data.ip], color=cols[data.ip], label="total p", clear=false,legend=:best, xlabel="Height / m", ylabel="Pressure / Pa")
	
	reveal(vis)
	
	
	#vis=GridVisualizer(layout=(2,1),resolution=(600,800), Plotter=PyPlot)
	#scalarplot!(vis[1,1], grid_, sol_[1])
	#scalarplot!(vis[2,1], grid_, sol_[data.iT] .- 273.15, show=true)
	#save("../img/out/pi_pt.svg", vis)
end

# ╔═╡ b66d9a7f-3023-4780-9999-5e07862fc583
	function simSens(P)
		Tv = eltype(P)
		dataSens=ModelDataSens{Tv}(S=P[1])
		for n in fieldnames(typeof(data_embed))
	   		if !(getfield(dataSens,n) isa Tv)
				v = getfield(data_embed,n)
				setfield!(dataSens, n, v)
			end
		end
		
		solSens,sysSens=runSensitivity(grid, sol, dataSens, Tv)
		
		Tca=Tcatavg(solSens,sysSens,grid,dataSens)
		YCO=Yield_CO(solSens,sysSens,dataSens)
		STC=STCefficiency(solSens,sysSens,dataSens)
		#[Tca, YCO, STC]
		[Tca, YCO]
	end

# ╔═╡ 7530df59-03e7-4bb6-83f2-86369edc13ee
begin
	if RunSensitivity
		# specify sensitivity parameter (α_cat, ϵ_cat, λ_eff, mcats)
		P = getfield(data_embed, SensPar)
		P *= 1.0 .+ [-0.5,-0.2,-0.1,0.0,0.1,0.2,0.5]
		
		dresult = DiffResults.JacobianResult(ones(2))
	    Tca = zeros(0)
	    DTca = zeros(0)
		YCO = zeros(0)
		DYCO = zeros(0)
	
	    for p ∈ P
			ForwardDiff.jacobian!(dresult, simSens, [p,p])
	        push!(Tca, DiffResults.value(dresult)[1])
	        push!(DTca, DiffResults.jacobian(dresult)[1])
			push!(YCO, DiffResults.value(dresult)[2])
	        push!(DYCO, DiffResults.jacobian(dresult)[2])			
	    end

	end

end

# ╔═╡ db16ecbe-cb69-46b6-8ac8-9ffa51c11dac
let
	vis = GridVisualizer(layout=(2,1), legend = :rb, xlabel=string(SensPar), markershape=:circle, markevery=1, Plotter=PyPlot)
	
    scalarplot!(vis[1,1], P, Tca; color = :red, label = "Tcat")
    scalarplot!(vis[1,1], P, DTca; color = :blue, label = "dTcat/d($(string(SensPar)))", clear = false, show = true)

	scalarplot!(vis[2,1], P, YCO; color = :red, label = "YCO")
    scalarplot!(vis[2,1], P, DYCO; color = :blue, label = "dYCO/d($(string(SensPar)))", clear = false, show = true)

	#save("../img/out/Sens_abs_cat.svg", vis)
end

# ╔═╡ Cell order:
# ╠═11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╟─863c9da7-ef45-49ad-80d0-3594eca4a189
# ╟─2ed3223e-a604-410e-93d4-016580f49093
# ╠═077d4ede-9e0f-4f94-afb0-01bd36c584fc
# ╠═8d90a6c3-95c5-4076-a3a1-2c01d119edb9
# ╟─f9205777-f2af-4192-88f5-a72e193f13df
# ╠═0a911687-aff4-4c77-8def-084293329f35
# ╠═ff58b0b8-2519-430e-8343-af9a5adcb135
# ╟─985718e8-7ed7-4c5a-aa13-29462e52d709
# ╠═ada45d4d-adfa-484d-9d0e-d3e7febeb3ef
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─f4dcde90-6d8f-4b17-b4ec-367d2372637f
# ╟─3703afb0-93c4-4664-affe-b723758fb56b
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
# ╟─c4521a0c-c5af-43cd-97bc-a4a7a42d27b1
# ╟─a6e61592-7958-4094-8614-e77446eb2223
# ╟─2739bcff-3fb0-4169-8a1a-2b0a14998cec
# ╠═30393c90-298c-412d-86ce-e36106613d35
# ╠═68e2628a-056a-4ec3-827f-2654f49917d9
# ╟─fec9ca6d-d815-4b50-bec7-f8fb3d8195ba
# ╟─808e5a71-572f-4b0c-aeb3-9513d969dada
# ╟─4cde8752-bbdf-4b83-869e-46b78bb4adb5
# ╟─0f7cce89-3add-4316-bcc2-a924064af884
# ╟─7bf4d925-56e1-4c85-a40a-fef1917f501e
# ╠═b3bf7b7d-eb38-4a32-87f7-5aef098ad03e
# ╟─aebcb161-425c-4a46-aaed-d4b13b3e6654
# ╟─6bbd0496-c275-4cd5-bf48-9ea4e77a9091
# ╟─1f3b2a0e-c76f-462c-a85c-e85815fd7e5a
# ╟─428afb22-55ab-4f64-805f-7e15ad4cf23f
# ╠═47161886-9a5c-41ac-abf5-bbea82096d5a
# ╟─b34c1d1b-a5b8-4de9-bea9-3f2d0503c1c0
# ╠═a9ec34cb-5c01-423f-a9b1-027cba13dca3
# ╟─8b1a0902-2542-40ed-9f91-447bffa4290f
# ╟─b06a7955-6c91-444f-9bf3-72cfb4a011ec
# ╠═f39dd714-972c-4d29-bfa8-d2c3795d2eef
# ╠═6ae6d894-4923-4408-9b77-1067ba9e2aff
# ╠═7ab38bc2-9ca4-4206-a4c3-5fed673557f1
# ╠═91d8119f-501f-49f2-93e0-88da8d996f7a
# ╟─b2df1087-6628-4889-8cd6-c5ee7629cd93
# ╠═3bd80c19-0b49-43f6-9daa-0c87c2ea8093
# ╟─2790b550-3105-4fc0-9070-d142c19678db
# ╟─bea97fb3-9854-411c-8363-15cbef13d033
# ╠═bd7552d2-2c31-4834-97d9-ccdb4652242f
# ╟─e25e7b7b-47b3-457c-995b-b2ee4a87710a
# ╠═3a35ac76-e1b7-458d-90b7-d59ba4f43367
# ╠═e73a3dfb-740e-4a20-9101-47cf18bcb9be
# ╟─6258c7c7-dbf5-4d5b-a8e3-de2250f5e66a
# ╠═2abbe557-9869-4861-9c3b-208f1f2edf59
# ╠═ec3dce96-a6df-4ff5-aa7b-ee8631555b67
# ╠═7530df59-03e7-4bb6-83f2-86369edc13ee
# ╠═b66d9a7f-3023-4780-9999-5e07862fc583
# ╠═4042a8de-8df8-4c9d-b1c9-6323afe20d95
# ╠═db16ecbe-cb69-46b6-8ac8-9ffa51c11dac
# ╟─cd2547b4-0e2f-4618-8cca-3e02ab7cb236
# ╟─1046820c-777f-4180-a570-1efe90924ec9
# ╠═2551d24d-91fe-40cd-8ade-eaf6b549efa3
