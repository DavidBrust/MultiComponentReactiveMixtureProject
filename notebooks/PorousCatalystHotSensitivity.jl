### A Pluto.jl notebook ###
# v0.19.22

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
	using PyPlot
	#using Plots
	using Colors
	#using CSV, DataFrames
	using Revise
	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="Dusty Gas Model")

# ╔═╡ 077d4ede-9e0f-4f94-afb0-01bd36c584fc
function catcylinder(;nref=0, r=data.D/2, h=data.h, catr=data.catD/2,cath=data.cath)
    	
	hr=r/10.0*2.0^(-nref)
	
	hh=h/10.0*2.0^(-nref)
	hhf=h/50.0*2.0^(-nref)
	
	
	Z=geomspace(0,h,hh,hhf)
	
    R=collect(0:hr:r)
    
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

# ╔═╡ 4af2237c-9144-4ffc-8966-2b3bf9d3d720
function mole_frac!(node,data,X,u)
	n=data.ng
	sump = 0.0
	for i=1:n
		if node isa VoronoiFVM.Edge
			X[i] = 0.5*(u[i,1]+u[i,2])
		else # boundary node
			X[i] = u[i]
		end
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
	vh=project(edge,(0,ud))
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

# ╔═╡ 02b76cda-ffae-4243-ab40-8d0fe1325776
md"""
##### Auxiliary functions
"""

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	ngas=data.ng
	ip=data.ip
	
	if node.region == 2 && data.isreactive # catalyst layer
		RRc=data.RRc
		ng=data.gni # gas species indices from names
		nr=data.kinpar.rni # reaction indices from names
		iT=data.iT
		#RR=RRc*u[n["H2"]]*u[n["CO2"]]
		pi = u[1:ngas]./ufac"bar"
		# negative sign: sign convention of VoronoiFVM: source term < 0 
		RR = -data.mcats*ri(data.kinpar,u[iT],pi)
		# reactions
		# R1: CO + H2O = CO2 + H2
		# R2: CH4 + 2 H2O = CO2 + 4 H2
		# R3: CH4 + H2O = CO + 3 H2
		f[ng["H2"]] = RR[nr["R1"]] + 4*RR[nr["R2"]] + 3*RR[nr["R3"]]
		f[ng["CO2"]] = RR[nr["R1"]] + RR[nr["R2"]]
		f[ng["CO"]] = -RR[nr["R1"]] + RR[nr["R3"]]
		f[ng["H2O"]] = -RR[nr["R1"]] - 2*RR[nr["R2"]] - RR[nr["R3"]]
		f[ng["CH4"]] = -RR[nr["R2"]] - RR[nr["R3"]]
		f[ng["N2"]] = 0.0
		# temperature eq. / heat source
		ΔHi=data.kinpar.ΔHi
		f[iT] = -(RR[nr["R1"]]*ΔHi["R1"]+RR[nr["R2"]]*ΔHi["R2"] +RR[nr["R3"]]*ΔHi["R3"])
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
	
	if bnode.region==3 || bnode.region==5 # top boundaries (cat layer & frit)
		ng=data.ng
		ip = data.ip
		iT = data.iT
		flux_rerad = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
				
		X=zeros(eltype(u), ng)
		mole_frac!(bnode,data,X,u)
		cf=heatcap_mix(data.Fluids, u[iT], X)
		
		flux_convec=data.utop*u[ip]/(ph"R"*u[iT])*cf*(u[iT]-data.Tamb)

		abs = 0.0
		if bnode.region==3 # frit
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

# ╔═╡ 40906795-a4dd-4e4a-a62e-91b4639a48fa
function bottom(f,u,bnode,data)
	if bnode.region==1 # bottom boundary
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
	boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)	
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

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
begin
	# sol_,grid,sys=main(data=ModelData());
	if sol_ isa VoronoiFVM.TransientSolution
		sol = sol_(sol_.t[end])
	else
		sol = sol_
	end
end;

# ╔═╡ f39dd714-972c-4d29-bfa8-d2c3795d2eef
function massflow(data, bflux)
	mdot=0.0
	for i=1:data.ng
		mdot += bflux[i] * data.Fluids[i].MW
	end
	mdot/ufac"kg/hr"
end

# ╔═╡ 229d36fc-b033-4aa5-a889-fba6ce461b29
data.mdotin/ufac"kg/hr"

# ╔═╡ 247e422a-f9f8-4350-9939-e308e0fef6cf
let
	T=sol[data.iT,:]
	sum(T)/length(T) - 273.15
end

# ╔═╡ 7ccd80aa-7be0-48c9-b971-7b54ab13977a
let
	function Tcat_(f,u,bnode,data)
		iT = data.iT
		#if bnode.region==5 # cat layer surface
			f[iT] = u[iT]
		#end
	end

	Tcat=integrate(sys,Tcat_,sol; boundary=true)[data.iT,5] / (π*data.catD^2/4)
	#Tbot=integrate(sys,Tcat,sol; boundary=true)[data.iT,1] / data.Ac
	#Tcat=integrate(sys,Tcat,sol; boundary=true)[data.iT,5] / (π*data.catD^2/4)
	#Tcat = Tcat .- 273.15
end

# ╔═╡ c4521a0c-c5af-43cd-97bc-a4a7a42d27b1
md"""
## Metrics of Reactor Operation
1. Catalyst layer average temperature
1. CO Yield (Converions x Selectivity)
1. Solar-to-chemical efficiency
"""

# ╔═╡ 3d660986-f6d7-41a6-800b-68ccd920c7ac
begin
	# bottom - inflow	
	Ibot=integrate(sys,bottom,sol; boundary=true)[:,1]
	# top - outflow
	Itop=integrate(sys,top,sol; boundary=true)[:,[3,5]]
	Itop=sum(Itop, dims=2)
end;

# ╔═╡ 8b1a0902-2542-40ed-9f91-447bffa4290f
md"""
Mass flows:

- through __bottom__ boundary __$(round(massflow(data, Ibot),sigdigits=3))__ kg/h
- through __top__ boundary __$(round(massflow(data, Itop),sigdigits=3))__ kg/h
"""

# ╔═╡ b06a7955-6c91-444f-9bf3-72cfb4a011ec
md"""

Chemical species flows through porous frit from bottom and top:

|    | Bottom     | Top     |    |
|----|-------|-------|-------|
| $(data.gn[1]) | $(abs(round(Ibot[1]/ufac"mol/hr",sigdigits=2)))  |   $(round(Itop[1]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[2]) | $(abs(round(Ibot[2]/ufac"mol/hr",sigdigits=2)))  |   $(round(Itop[2]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[3]) | $(abs(round(Ibot[3]/ufac"mol/hr",sigdigits=2)))  |   $(round(Itop[3]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[4]) | $(abs(round(Ibot[4]/ufac"mol/hr",sigdigits=2)))  |   $(round(Itop[4]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[5]) | $(abs(round(Ibot[5]/ufac"mol/hr",sigdigits=2)))  |   $(round(Itop[5]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data.gn[6]) | $(abs(round(Ibot[6]/ufac"mol/hr",sigdigits=2)))  |   $(round(Itop[6]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |


"""

# ╔═╡ 91d8119f-501f-49f2-93e0-88da8d996f7a
function vel_darcy(data, ∇p, T, X)
	μ=dynvisc_mix(data, T, X)
	data.k/μ * ∇p # m/s
end

# ╔═╡ b2df1087-6628-4889-8cd6-c5ee7629cd93
md"""
### Temperature Plot
"""

# ╔═╡ 3bd80c19-0b49-43f6-9daa-0c87c2ea8093
let
	iT=data.iT
	vis=GridVisualizer(resolution=(600,400), Plotter=PyPlot)
	scalarplot!(vis, grid, sol[iT,:].- 273.15, show=true)

	#reveal(vis)
	#save("../img/out/Temeprature.svg", vis)
end

# ╔═╡ 2790b550-3105-4fc0-9070-d142c19678db
md"""
### Partial Pressure Plots
"""

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

# ╔═╡ bd7552d2-2c31-4834-97d9-ccdb4652242f
function plane(data,sol)
	#grid=cylinder()
	grid=catcylinder()
	xcut=data.D/4
	#xcut=data.D/2 * 9/10
	bfacemask!(grid, [xcut,0],[xcut,data.h],6)

	# transform y coordinate of parent grid into x coordinate of subgrid
	function _2to1(a,b)
		a[1]=b[2]
	end
	grid_1D  = subgrid(grid, [6], boundary=true, transform=_2to1) 

	sol_p = []
	for i=1:(data.ng+1)
		sol_i = view(sol[i, :], grid_1D)
		push!(sol_p, collect(sol_i))
	end
		
	sol_p, grid_1D	
	#sol_cutplane, grid_2D	
end

# ╔═╡ 99638215-3197-492a-8ac3-dcbf16f41fbe
function pdrop(data, sol)
	ip=data.ip
	Qflows=[5000.0,50000.0]*ufac"ml/minute"

	sol_0, grid_0=plane(data,sol(Qflows[1]))
	sol_1, grid_1=plane(data,sol(Qflows[end]))
	Δp1 = (sol_1[ip][1]-sol_1[ip][end])
end

# ╔═╡ 4b41d985-8ebc-4cab-a089-756fce0d3060
let	
	ngas=data.ng
	ip=data.ip
	#vis=GridVisualizer(layout=(2,1),resolution=(800,1000), Plotter=PyPlot)
	vis=GridVisualizer(layout=(2,1),resolution=(600,800))

	sol=sol_(sol_.t[end])
	
	
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
	#save("../img/out/partial_pressure.svg", vis)

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
Base.@kwdef mutable struct ModelData{Tv} <:AbstractModelData
	S::Tv = 1.0
	
	# catalyst / chemistry data
	# kinetic parameters, S3P="simple 3 parameter" kinetics fit to UPV lab scale experimental data
	kinpar::AbstractKineticsData = S3P
	
	# number of gas phase species
	ng::Int64		 		= S3P.ng
	#ng::Int64		 		= 1
	# names and fluid indices
	gn::Dict{Int, String} 	= S3P.gn
	#gn::Dict{Int, String} 	= Dict(1=>"N2")
	# inverse names and fluid indices
	gni::Dict{String, Int}  = S3P.gni
	# fluids and respective properties in system
	Fluids::Vector{AbstractFluidProps} = S3P.Fluids
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = let
		x=zeros(Float64, ng)
		x[gni["H2"]] = 1.0
		x[gni["CO2"]] = 1.0
		x/sum(x)
	end # inlet composition
	#X0::Vector{Float64} = [1.0]
	

	# reaction enthalpy
	ΔHR::Float64 = 10.0*ufac"kJ/mol"
	# volume specific cat mass loading, UPV lab scale PC reactor
	mcats::Float64 =1234.568*ufac"kg/m^3"
	isreactive::Bool = 1
	#isreactive::Bool = 0

	ip::Int64=ng+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable

	
	RRc::Float64 = 1.0e-5 # reaction rate constant	
		
	α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	
	## irradiation data
	Tau_quartz::Float64=0.9 # transmission coefficient of quartz window
	G_lamp::Float64=0.0*ufac"kW/m^2" # solar simulator irradiation flux
	Abs_lamp_cat::Float64=0.7 # cat avg absorptivity of irradiation from lamp
	Abs_lamp_frit::Float64=0.7 # frit avg absorptivity of irradiation from lamp
	Eps_ir::Float64=0.7 # avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
		
	
	## porous filter data
	dp::Float64=50.0*ufac"μm" # average pore size
	#dp::Float64=25.0*ufac"μm" # average pore size
	cath::Float64 = 500.0*ufac"μm" # catalyst layer thickness
	catD::Float64 = 10.0*ufac"cm" # catalyst layer diameter
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
	#k::Float64=2.9e-14*ufac"m^2" # permeability
	#k::Float64=2.9e-6*ufac"m^2" # permeability
	a_s::Float64=0.13*ufac"m^2/g" # specific surface area
	ρfrit::Float64=(1.0-ϕ)*ρs*ufac"kg/m^3" # density of porous frit
	a_v::Float64=a_s*ρfrit # volume specific interface area
	## END porous filter data


	## Flow data
	#norm conditions
	pn::Float64 = 1.0*ufac"bar"
	Tn::Float64 = 273.15*ufac"K"
	
	Qflow::Float64=3400.0*ufac"ml/minute" # volumetric feed flow rate (sccm)

	MWin::Float64 = molarweight_mix(Fluids, X0)
	mdotin::Float64=MWin*Qflow*pn/(ph"R"*Tn)*ufac"kg/s"
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

	sol_=solve(sys;inival,)
	
	function pre(sol,par)
		ng=data.ng
		iT=data.iT
		# iteratively adapt top outflow boundary condition
		function Inttop(f,u,bnode,data)
			
			X=zeros(eltype(u), ng)
			mole_frac!(bnode,data,X,u)
			if bnode.region==3 # top boundary
				for i=1:ng
					f[i] = data.Fluids[i].MW*X[i]
				end
				f[iT] = u[iT]
			end
		end
		
		MWavg=sum(integrate(sys,Inttop,sol; boundary=true)[1:ng,3])/data.Ac
		ntop=data.mdotin/MWavg
		Tavg=integrate(sys,Inttop,sol; boundary=true)[data.iT,3]/data.Ac
		utop_calc=ntop*ph"R"*Tavg/(1.0*ufac"bar")/data.Ac
		utops=[data.utop, utop_calc]*ufac"m/s"
		data.utop = minimum(utops) + par*(maximum(utops)-minimum(utops))
		
		# embedding parameter: Qflow / ml/minute
		#Qflows = [5000.0,50000.0]*ufac"ml/minute"
		#Qflow = minimum(Qflows) + par*(maximum(Qflows)-minimum(Qflows))
		#data.Qflow=Qflow
		#data.u0=Qflow/(data.Ac)*ufac"m/s" # mean superficial velocity
		#data.utop=data.u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate

		
		# specific catalyst loading
		mcats=[10.0, 1300.0]*ufac"kg/m^3"
		data.mcats= minimum(mcats) + par*(maximum(mcats)-minimum(mcats))

		# irradiation flux density
		G_lamp=[1.0, 125.0]*ufac"kW/m^2"
		data.G_lamp= minimum(G_lamp) + par*(maximum(G_lamp)-minimum(G_lamp))
	end
	
	control=SolverControl(;
					  handle_exceptions=true,
					  Δp_min=5.0e-3,
					  Δp=0.1,
					  Δp_grow=1.2,
					  Δu_opt=10000.0, # large value, due to unit Pa of pressure?
					  )
	
	#sol=solve(sys;inival,)
	
	sol=solve(sys;inival,embed=[0.0,0.999,1.0],pre,control)

	
	sol,grid,sys
end;

# ╔═╡ 7530df59-03e7-4bb6-83f2-86369edc13ee
let
	function g(P)
 		Tv = eltype(P)
		data=ModelData(S=one(Tv))
		[1.0]
	end
	
	dresult = DiffResults.JacobianResult(ones(1))

	P = 0.1:0.05:2
    G = zeros(0)
    DG = zeros(0)
    @time for p ∈ P
        ForwardDiff.jacobian!(dresult, g, [p])
        push!(G, DiffResults.value(dresult)[1])
        push!(DG, DiffResults.jacobian(dresult)[1])
    end



end

# ╔═╡ Cell order:
# ╠═11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╟─863c9da7-ef45-49ad-80d0-3594eca4a189
# ╠═077d4ede-9e0f-4f94-afb0-01bd36c584fc
# ╠═8d90a6c3-95c5-4076-a3a1-2c01d119edb9
# ╠═b659879f-3706-4937-b655-116f1be85069
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
# ╠═f39dd714-972c-4d29-bfa8-d2c3795d2eef
# ╠═229d36fc-b033-4aa5-a889-fba6ce461b29
# ╠═247e422a-f9f8-4350-9939-e308e0fef6cf
# ╟─8b1a0902-2542-40ed-9f91-447bffa4290f
# ╟─b06a7955-6c91-444f-9bf3-72cfb4a011ec
# ╠═7ccd80aa-7be0-48c9-b971-7b54ab13977a
# ╟─c4521a0c-c5af-43cd-97bc-a4a7a42d27b1
# ╠═3d660986-f6d7-41a6-800b-68ccd920c7ac
# ╠═99638215-3197-492a-8ac3-dcbf16f41fbe
# ╠═91d8119f-501f-49f2-93e0-88da8d996f7a
# ╟─b2df1087-6628-4889-8cd6-c5ee7629cd93
# ╠═3bd80c19-0b49-43f6-9daa-0c87c2ea8093
# ╟─2790b550-3105-4fc0-9070-d142c19678db
# ╠═4b41d985-8ebc-4cab-a089-756fce0d3060
# ╟─8db106cc-c5f1-498b-bf8b-fddd8e21b444
# ╠═bd7552d2-2c31-4834-97d9-ccdb4652242f
# ╟─2f24ba85-748a-44a6-bde1-8e320106198d
# ╟─a96be76c-df66-4149-994c-5b14cc0c3a37
# ╠═3a35ac76-e1b7-458d-90b7-d59ba4f43367
# ╠═7530df59-03e7-4bb6-83f2-86369edc13ee
