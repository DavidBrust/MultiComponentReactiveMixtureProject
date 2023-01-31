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
	using Plots
	using Colors
	using CSV, DataFrames

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="Dusty Gas Model")

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

# ╔═╡ 2f2007ae-df64-4001-944e-1e0b6e29ac1d
function grid_(data;nref=0)
	X=range(0,data.R0,length=10*2^nref+1)
	grid=simplexgrid(X)
	spherical_symmetric!(grid)
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
	D=data.D_0*data.ϵ_τ
	D_K_eff=data.D_K_eff
	M=zeros(eltype(x), n, n)
	for i=1:n
		M[i,i] = -1/D_K_eff[i]
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
	ngas=data.ng
	ip=data.ip
	B0=data.B0
	μ=data.μ
	T=data.Tamb
	DK=data.D_K_eff
	F=zeros(eltype(u), ngas)
	X=zeros(eltype(u), ngas)
	δp = u[ip,1]-u[ip,2]
	ud = B0/μ*δp # darcy velocity
	pm = 0.5*(u[ip,1]+u[ip,2])
	
	for i=1:ngas
		δpi = u[i,1]-u[i,2]
		pim = 0.5*(u[i,1]+u[i,2])
		#F[i] = (δpi + pim/DK[i]*B0/μ*δp)/(ph"R"*T)

		bp,bm=fbernoulli_pm(ud/DK[i])
		F[i] = -(bm*u[i,1]-bp*u[i,2])/(ph"R"*T)
		#F[i] = (δpi)/(ph"R"*T)
		X[i] = pim/pm
	end
	
    
	# computation of fluxes J     

	#J = M_matrix(data,density_means(u,n)) \ F
	J = M_matrix(data,X) \ F

	
	f[1:ngas] = J
	#f[3] via reaction: ∑xi = 1
	#F
	
end

# ╔═╡ 02b76cda-ffae-4243-ab40-8d0fe1325776
md"""
##### Auxiliary functions
"""

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	ngas=data.ng
	ip=data.ip
	RRc=data.RRc
	m=data.m
	
	RR=RRc*u[1]
	f[1] = RR
	f[2] = -m*RR
	
	# ∑xi = 1
	f[ip]=u[ip]-sum(u[1:ngas])
	#f[ip]=log(sum(u[1:ngas])) -log(u[ip])
	
end

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
	R0=data.R0
	p0=data.pamb
	boundary_dirichlet!(f,u,bnode,1,2,1.0*p0) # A
	boundary_dirichlet!(f,u,bnode,2,2,0.0*p0) # B
	
end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
function main(;data=ModelData(),nref=1)
	grid=grid_(data,nref=nref)
	ngas=data.ng
	R0=data.R0
	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)
	enable_species!(sys; species=collect(1:(ngas+1)))

	inival=unknowns(sys)
	inival[:,:].=1.0*data.pamb


	
	sol=solve(sys;inival)

	
	sol,grid,data,sys
end;

# ╔═╡ ea20ff65-4ecb-4847-af46-37b71d474efd
md"""
# Model Verification
"""

# ╔═╡ 50e852a3-0aae-4cf3-9950-5548e34545a4
md"""
Compare against data from __Veldsink, J. W., van Damme, R. M. J., Versteeg, G. F., & van Swaaij, W. P. M. (1995).__ The use of the dusty-gas model for the description of mass transport with chemical reaction in porous media. The Chemical Engineering Journal and the Biochemical Engineering Journal, 57(2), 115-125. doi:10.1016/0923-0467(94)02929-6

"""

# ╔═╡ 39452e46-ac88-4900-8f15-d8d67f9ac993
md"""
$(LocalResource("../img/params_DGM.png"))
"""

# ╔═╡ d6178bb7-4d33-4102-8d52-d1c2ef2f70b9
begin
	p_p0_dat = CSV.read("../data/Veldsink1995/p_p0.csv", DataFrame, delim=";", header=["r_R0","p_P0"] )
	xA_dat = CSV.read("../data/Veldsink1995/x_A.csv", DataFrame, delim=";", header=["r_R0","xA"] )
	xB_dat = CSV.read("../data/Veldsink1995/x_B.csv", DataFrame, delim=";", header=["r_R0","xB"] )
end;

# ╔═╡ 1bd2ffdb-c364-4a3f-8955-9de8273d2acd
let
	sol,grid,data,sys=main(nref=5)
	tf=VoronoiFVM.TestFunctionFactory(sys)
	Γ_where_T_equal_1=[2]
	Γ_where_T_equal_0=[1]
	T=testfunction(tf,Γ_where_T_equal_0,Γ_where_T_equal_1)
	I=integrate(sys,T,sol)
end

# ╔═╡ cb908660-7fcb-4e2d-b1c3-f02d60a31221
function thiele_mod(data)
	δ=data.R0/3
	n=1
	kn=data.RRc
	D_Ae=data.Dg[1]*data.ϵ_τ
	δ*sqrt( (n+1)*kn*ph"R"*data.Tamb/(2*D_Ae))
end

# ╔═╡ b94129e3-684e-4b97-94af-2d3c029a0244
md"""
Reaction rate coefficinet pre factor, used to matched the rate of reaction that was used in Veldsink 1995 and described there as "fast kinetics".
"""

# ╔═╡ bc0bfd62-6057-4c07-89ac-cde6e8c840e3
md"""
$(@bind RRc_mod Slider(range(1.0,5.0,length=41),show_value=true,default=2.4))
"""

# ╔═╡ 5abd4e0e-bcb3-4657-a73d-f558dc404e7f
thiele_mod(ModelData(RRc=5.0e-2*RRc_mod))

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
sol,grid=main(data=ModelData(RRc=5.0e-2*RRc_mod),nref=2);

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

# ╔═╡ 62dfbf09-c2dc-4fb9-82b8-07c631a5f826
md"""
From __Veldsink, J. W., van Damme, R. M. J., Versteeg, G. F., & van Swaaij, W. P. M. (1995).__ The use of the dusty-gas model for the description of mass transport with chemical reaction in porous media. The Chemical Engineering Journal and the Biochemical Engineering Journal, 57(2), 115-125. doi:10.1016/0923-0467(94)02929-6
"""

# ╔═╡ 3a35ac76-e1b7-458d-90b7-d59ba4f43367
Base.@kwdef mutable struct ModelData <:AbstractModelData
	# number of gas phase species
	ng::Int64		 		= 2
	
	iT::Int64=1 # index of Temperature variable
	ip::Int64=ng+1 # index of total pressure variable

	m::Float64=3.0 # stochiometry factor A -> mB
	Mi::Vector{Float64} = 20.0e-3*[1.0, 1.0/m]*ufac"kg/mol" # molar masses
	Tamb::Float64=600.0*ufac"K" # ambient temperature
	pamb::Float64=1.0e5*ufac"Pa"

	# Parameters from Veldsink1995 / DGM
	R0::Float64=1.0e-3 # particle radius
	rp::Float64=2.9e-7*ufac"m" # pore radius inside cat. particles
	B0::Float64=4.79e-16*ufac"m^2" # permeability
	ϵ_τ::Float64=0.0456 # porosity/tourtuosity factor
	
	# 1:CO, 2:H2, 3:CO
	gn::Dict{Int, String} 	= Dict(1=>"CO", 2=>"H2", 3=>"CO2", 4=>"N2")
	#Dg::Vector{Float64} 	= [1.0,10.0,3.0,8.0,2.0,1.8]*1e-5*ufac"cm^2/s"
	Dg::Vector{Float64} 	= [1.5]*1e-4*ufac"m^2/s"
	D_0::Symmetric{Float64, Matrix{Float64}} = Dmatrix(ng, Dg)
	D_K_eff::Vector{Float64} = @. 2.0/3.0*ϵ_τ*rp*sqrt((8*ph"R"*Tamb)/(π*Mi))*ufac"m^2/s"
	
	RRc::Float64 = 5.0e-2 # reaction rate constant
	

	
	
	μ::Float64=1.0e-5*ufac"Pa*s" # gas viscosity
	
	#α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	#α_nc::Float64=15.0*ufac"W/(m^2*K)" # natural convection heat transfer coefficient

	## irradiation data
	#G_lamp::Float64=1.0*ufac"kW/m^2" # solar simulator irradiation flux
	#Abs_lamp::Float64=0.7 # avg absorptivity of cat. of irradiation coming from lamp
	#Eps_ir::Float64=0.7 # avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
		
	
	## porous filter data
	d::Float64=100.0*ufac"μm" # average pore size
	# cylindrical disc / 2D
    D::Float64=12.0*ufac"cm" # disc diameter
	

	# prism / 3D
	wi::Float64=12.0*ufac"cm" # prism width/side lenght
	le::Float64=wi # prism width/side lenght
	h::Float64=0.5*ufac"cm" # frit thickness (applies to 2D & 3D)

	#Ac::Float64=pi*D^2.0/4.0*ufac"m^2" # cross-sectional area, circular
	Ac::Float64=wi^2*ufac"m^2" # cross-sectional area, square
	
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2 	
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
	
	ϕ::Float64=0.36 # porosity, class 2
	k::Float64=2.9e-11*ufac"m^2" # permeability
	a_s::Float64=0.13*ufac"m^2/g" # specific surface area
	ρfrit::Float64=(1.0-ϕ)*ρs+ϕ*density_idealgas(Air, 298.15, 1.0*ufac"atm")*ufac"kg/m^3" # density of porous frit
	a_v::Float64=a_s*ρfrit # volume specific interface area
	## END porous filter data

	## fluid data
	
	Qflow::Float64=3400.0*ufac"ml/minute" # volumetric feed flow rate
	Tin::Float64=298.15*ufac"K" # inlet temperature
	p::Float64=1.0*ufac"atm" # reactor pressure		
	# u0::Float64=Qflow/(Ac*ϕ)*ufac"m/s" # mean superficial velocity
	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	# fluid properties: Air
	# values taken from VDI heat atlas 2010 chapter D3.1
	Fluid::FluidProps=Air
	## END fluid data
	
end;

# ╔═╡ 7530df59-03e7-4bb6-83f2-86369edc13ee
data=ModelData()

# ╔═╡ 3a95f282-2347-4dc8-bf0e-b6d55be65323
gridplot(grid_(data),resolution=(600,200),legend=:rt,show=true)

# ╔═╡ 65f8a11b-c7bd-4df4-846b-e0b658e8f570
let
	ngas=data.ng
	par=data.B0*[0.1,1.0,10.0]
	lp=length(par)
	c1 = colorant"red"
	c2 = colorant"blue"
	cols=range(c1, stop=c2, length=ngas*lp)
	lbls=["A","B"]
	vis=GridVisualizer(resolution=(600,300))
	for (j,p) in enumerate(par)
		sol,grid,data=main(data=ModelData(B0=p),nref=3)
		for i in 1:ngas
			scalarplot!(vis, grid, sol[i,:], color=cols[(j-1)*ngas+i], label="$(lbls[i]) $(data.B0)", clear=false, legend=:ct)
		end
	end
	reveal(vis)
end

# ╔═╡ cecfe64e-1caa-447e-b3ba-a12ecf05e805
let
	# Veldsink 1995 data
	Plots.plot(xguide="Non-dim radius / -", yguide="Mole frac or Non-dim pressure / -", legend=:outertopright)
	
	Plots.scatter!(xA_dat.r_R0, xA_dat.xA, label="xA Veldsink95", c=1)
	Plots.scatter!(xB_dat.r_R0, xB_dat.xB, label="xB Veldsink95", c=2)
	Plots.scatter!(p_p0_dat.r_R0, p_p0_dat.p_P0, label="p/p0 Veldsink95", c=3)

	# DGM model
	p_tot = sol[3,:]
	xA = sol[1,:] ./ p_tot
	xB = sol[2,:] ./ p_tot
	p_p0 = p_tot / data.pamb
	r_R0 = grid[Coordinates] / data.R0
	r_R0 = vec(r_R0)

	Plots.plot!(r_R0, xA, label="xA DGM", c=1)
	Plots.plot!(r_R0, xB, label="xB DGM", c=2)
	Plots.plot!(r_R0, p_p0, label="p/p0 DGM", c=3)	

end

# ╔═╡ 4b41d985-8ebc-4cab-a089-756fce0d3060
let	
	ngas=data.ng
	ip=data.ip
	# visualization
	vis=GridVisualizer(resolution=(600,300))
	

	#lt=length(sol.t)
	c1 = colorant"red"
	c2 = colorant"blue"
	cols=range(c1, stop=c2, length=ngas)

	lbls=["A","B"]

	#for (j,t) in enumerate(sol.t)
	#	_sol=sol(t)
		for i in 1:ngas
			scalarplot!(vis, grid, sol[i,:], color=cols[i], label="$(lbls[i])", clear=false, legend=:ct)
		end
	#end

	#scalarplot!(vis, grid, sol[ip,:], color=cols[ip], label="total p", clear=false, legend=:ct)

	reveal(vis)

end

# ╔═╡ 154ab2c8-1956-4258-8c10-f2b0314ed069
data.D_K_eff

# ╔═╡ b4edce69-a693-44c4-b0f1-604ce30ecfad
data.D_0*data.ϵ_τ

# ╔═╡ Cell order:
# ╠═11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╟─863c9da7-ef45-49ad-80d0-3594eca4a189
# ╠═fa2002f1-1b35-4fe8-bbd1-f9fff62acbdf
# ╠═2f2007ae-df64-4001-944e-1e0b6e29ac1d
# ╠═b659879f-3706-4937-b655-116f1be85069
# ╠═3a95f282-2347-4dc8-bf0e-b6d55be65323
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─f4dcde90-6d8f-4b17-b4ec-367d2372637f
# ╟─3703afb0-93c4-4664-affe-b723758fb56b
# ╟─c4e08410-3d1a-4921-aa5c-2f232daaa07a
# ╟─313cd88a-1497-4d3a-b09a-9b98a6dad9c2
# ╟─21d0195b-b170-460d-989e-f9d00b511237
# ╟─8f4843c6-8d2b-4e24-b6f8-4eaf3dfc9bf0
# ╟─66b55f6b-1af5-438d-aaa8-fe4745e85426
# ╟─8528e15f-cce7-44d7-ac17-432f92cc5f53
# ╟─19a46c60-d7b5-45da-a3c5-3687316b2a2c
# ╟─d5987ad8-e71b-4d40-b1c1-88142558265d
# ╟─d2adc467-d366-4461-9a00-9832e7d5a737
# ╠═ed7941c4-0485-4d84-ad5b-383eb5cae70a
# ╟─a6afe118-dcbd-4126-8646-c7268acfacf3
# ╠═f3dc68cb-f9c3-4e85-892f-34358831e3eb
# ╠═b6381008-0280-404c-a86c-9c9c3c9f82eb
# ╟─02b76cda-ffae-4243-ab40-8d0fe1325776
# ╠═78cf4646-c373-4688-b1ac-92ed5f922e3c
# ╠═29d66705-3d9f-40b1-866d-dd3392a1a268
# ╠═333b5c80-259d-47aa-a441-ee7894d6c407
# ╠═65f8a11b-c7bd-4df4-846b-e0b658e8f570
# ╟─ea20ff65-4ecb-4847-af46-37b71d474efd
# ╟─50e852a3-0aae-4cf3-9950-5548e34545a4
# ╟─39452e46-ac88-4900-8f15-d8d67f9ac993
# ╠═d6178bb7-4d33-4102-8d52-d1c2ef2f70b9
# ╠═1bd2ffdb-c364-4a3f-8955-9de8273d2acd
# ╠═cb908660-7fcb-4e2d-b1c3-f02d60a31221
# ╠═5abd4e0e-bcb3-4657-a73d-f558dc404e7f
# ╟─cecfe64e-1caa-447e-b3ba-a12ecf05e805
# ╠═aa498412-e970-45f2-8b11-249cc5c2b18d
# ╟─b94129e3-684e-4b97-94af-2d3c029a0244
# ╠═bc0bfd62-6057-4c07-89ac-cde6e8c840e3
# ╠═4b41d985-8ebc-4cab-a089-756fce0d3060
# ╟─8c639001-b3ef-494c-a3eb-c600b011c159
# ╟─bba07772-a923-482e-b24b-acfaa24ac81c
# ╟─62dfbf09-c2dc-4fb9-82b8-07c631a5f826
# ╠═154ab2c8-1956-4258-8c10-f2b0314ed069
# ╠═b4edce69-a693-44c4-b0f1-604ce30ecfad
# ╠═3a35ac76-e1b7-458d-90b7-d59ba4f43367
# ╠═7530df59-03e7-4bb6-83f2-86369edc13ee
