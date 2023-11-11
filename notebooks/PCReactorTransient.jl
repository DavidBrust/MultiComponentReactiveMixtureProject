### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ c21e1942-628c-11ee-2434-fd4adbdd2b93
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using VoronoiFVM
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve, LinearSolve
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots
	using PlutoUI, Colors
	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 6da83dc0-3b0c-4737-833c-6ee91552ff5c
md"""
Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=false))
"""

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
PlutoUI.TableOfContents(title="M-S Transport + Darcy")

# ╔═╡ 83fa22fa-451d-4c30-a4b7-834974245996
function grid1D()
	X=(0:0.02:1)*ufac"cm"
	grid=simplexgrid(X)
	# catalyst region
	cellmask!(grid,[0.4]*ufac"cm",[0.6]*ufac"cm",2)	
	#cellmask!(grid,[0.0]*ufac"cm",[0.2]*ufac"cm",2)	
	grid
end

# ╔═╡ 4dae4173-0363-40bc-a9ca-ce5b4d5224cd
function grid2D()
	R=(0:1:7)*ufac"cm"
	Z=(0:0.05:0.5)*ufac"cm"
	grid=simplexgrid(R,Z)
	circular_symmetric!(grid)

	cellmask!(grid,[0.0,0.4].*ufac"cm",[6.0,0.5].*ufac"cm",2) # catalyst region
	bfacemask!(grid, [0.0,0.5].*ufac"cm",[6.0,0.5].*ufac"cm",5)
	
	grid
end

# ╔═╡ 561e96e2-2d48-4eb6-bb9d-ae167a622aeb
function grid3D()
	X=(0:1:14)*ufac"cm"
	Y=(0:0.05:0.5)*ufac"cm"
	Z=(0:1:14)*ufac"cm"
	grid=simplexgrid(X,Y,Z)

	# catalyst region
	cellmask!(grid,[2,0.0,2].*ufac"cm",[12,0.1,12].*ufac"cm",2)
	#bfacemask!(grid, [2,1,2].*ufac"cm",[12,1,12].*ufac"cm",7) # mask outflow
	bfacemask!(grid, [2,0,2].*ufac"cm",[12,0,12].*ufac"cm",7) # mask inflow
	
	grid
end

# ╔═╡ 107a6fa3-60cb-43f0-8b21-50cd1eb5065a
const dim = 2

# ╔═╡ 4e05ab31-7729-4a4b-9c14-145118477715
if dim == 3
	@bind xcut Slider(linspace(0,14,21)*ufac"cm",show_value=true,default=6.5*ufac"cm")
end

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
let
	if dim == 1
		gridplot(grid1D(), resolution=(600,200))
	elseif dim == 2
		gridplot(grid2D())
	else
		gridplot(grid3D(); xplane=xcut, show=true, outlinealpha=0.0 )
	end
end

# ╔═╡ 832f3c15-b75a-4afe-8cc5-75ff3b4704d6
begin
	if dim == 1
		const Γ_left = 1
		const Γ_right = 2
	elseif dim == 2
		const Γ_out = 1
		const Γ_sym = 4		
		const Γ_outer = 2
		#const Γ_left = 4
		const Γ_in = 5
	else
		# mask inflow
		const Γ_left = 7 
		const Γ_right = 3
		# mask outflow
		# const Γ_left = 1 
		# const Γ_right = 7		
		const Γ_front = 2
		const Γ_back = 4
		const Γ_bottom = 5
		const Γ_top = 6
	end
end;

# ╔═╡ a078e1e1-c9cd-4d34-86d9-df4a052b6b96
md"""
Multi-component transport of ideal-gas mixture through porous medium.
"""

# ╔═╡ 0fadb9d2-1ccf-4d44-b748-b76d911784ca
md"""
## Overall Mass Continuity
Mixture mass flow (overall mass flow) through the pore space of the porous medium. Mixture mass averaged velocity is calculated from Darcy equation. The void fraction (porosity) is given by $\epsilon$.

```math
\begin{align}
	\frac{\partial \epsilon \rho}{\partial t} + \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\end{align}
```
"""

# ╔═╡ b94513c2-c94e-4bcb-9342-47ea48fbfd14
md"""
## Species Mass Continuity and Transport
```math
\begin{align}
	\frac{\partial \epsilon \rho_i}{\partial t} + \nabla \cdot \left( \vec \Phi_i + \rho_i \vec v \right ) - R_i &= 0 ~, \qquad i = 1 ... \nu \\
		\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= -\sum_{j=1 \atop j \neq i}^{\nu} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
		\sum_{i=1}^\nu x_i &= 1
\end{align}
```
"""

# ╔═╡ c886dd12-a90c-40ab-b9d0-32934c17baee
md"""
where $\rho$ is the (total) mixture density, $\vec v$ is the mass-averaged (barycentric)  mixture velocity calculated with the Darcy equation, $x_i$, $w_i$ and $M_i$ are the molar fraction, mass fraction and molar mass of species $i$ respectively, $\vec \Phi_i$ is the mass flux of species $i$ ($\frac{\text{kg}}{\text{m}^2 \text{s}}$) and $R_i$ is the species mass volumetric source/sink ($\frac{\text{kg}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
"""

# ╔═╡ 8f2549f4-b0a6-440f-af94-6880e0814dc2
md"""
## Thermal Energy Transport
"""

# ╔═╡ 78589a1e-2507-4279-ba42-1aaec90d87d0
md"""
```math
\begin{align}
- \nabla \cdot \left(\lambda_{\text{eff}} \nabla T - \sum_i^{\nu} \vec N_i h_i \right)  &= 0
\end{align}
```

where $\lambda_{\text{eff}}$ is the effective thermal conductivity. The heat release from chemical reactions is considered as part of the species enthalpies. The convective heat transport within the porous medium is expressed via the sum over the product of species (mass) fluxes with species mass-specific enthalpies $\vec H_i= \vec N_i h_i$.

As part of the initialisation strategy (see next section) the convective contribution of heat flux is ramped up at the same time with the heat transport boundary conditions after the flow field has been established.
"""

# ╔═╡ 927dccb1-832b-4e83-a011-0efa1b3e9ffb
md"""
## Initialisation
The simulation is setup as a transient simulation. An initialisation strategy is employed where different physics are enabled step by step once a stationary state is established. Initially, no heat is transported and no chemical reactions take place. 

1. Velocity field (mass flow is ramped up from 1-100 % in T=[0;1] s)
2. Temperature field (heat transport is ramped up from 0-100 % in T=[1;2] s)
3. Catalyst reactivity (ramped up from 0-100 % in T=[2;3] s)

The mass flow boundary condition into the reactor domain is "ramped up" starting from a low value and linearly increasing until the final value is reached. A time delay is given to let the flow stabilize. Once the flow field is established, heat transport is ramped up until a stable temperature field is established. Finally, the reactivity of the catalyst is "ramped up" until its final reactivity value is reached.
"""

# ╔═╡ 68ca72ae-3b24-4c09-ace1-5e340c8be3d4
function len(grid)
	coord = grid[Coordinates]
	L=0.0
	if dim == 1
		L=coord[end]
	elseif dim == 2
		L=coord[1,end]
	else
		L=coord[2,end]
	end
	L*ufac"m"
end

# ╔═╡ db77fca9-4118-4825-b023-262d4073b2dd
md"""
### Peclet Number
```math
\text{Pe}_L= \frac{L \vec v}{D}

```
"""

# ╔═╡ e7497364-75ef-4bd9-87ca-9a8c2d97064c
md"""
### Damköhler Number
```math
\text{Da}= \frac{k_{\text{react}}}{k_{\text{conv}}} = \frac{\dot r}{c_0} \tau = \frac{\dot r L}{c_0 v}

```
"""

# ╔═╡ f6e54602-b63e-4ce5-b8d5-f626fbe5ae7a
md"""
### Reynolds Number
```math
\text{Re} = \frac{\rho v D}{\eta}  
```
where $\rho$ is the density of the fluid, ideal gas in this case, $v$ is the magnitufe of the velocity, $D$ is a representative length scale, here it is the mean pore diameter, and $\eta$ is the dynamic viscosity of the fluid.
"""

# ╔═╡ 0d507c5e-eb0f-4094-a30b-a4a82fd5c302
md"""
### Knudsen Number
```math
\text{Kn} = \frac{\lambda}{l} = \frac{k_{\text B}T}{\sqrt 2 \pi \sigma^2pl}
```
where $\lambda$ is the mean free path length of the fluid, ideal gas in this case and $l$ is the pore diameter of the porous medium. Determine the mean free path length for ideal gas from kinetic gas theory and kinetic collision cross-sections (diameter).

-  $\text{Kn} < 0.01$: Continuum flow
-  $0.01 <\text{Kn} < 0.1$: Slip flow
-  $0.1 < \text{Kn} < 10$: Transitional flow
-  $\text{Kn} > 10$: Free molecular flow
"""

# ╔═╡ 67adda35-6761-4e3c-9d05-81e5908d9dd2
md"""
## Model Data
"""

# ╔═╡ f40a8111-b5cb-40f0-8b12-f57cf59637f1
begin
	
Base.@kwdef mutable struct ModelData{NG}
	kinpar::FixedBed.KinData{nreac(XuFroment)} = XuFroment
	mcat::Float64=2000.0*ufac"mg"
	Vcat::Float64=1.0*ufac"m^2"*0.02*ufac"cm"
	lcat::Float64=mcat/Vcat
	
	#isreactive::Bool = 0
	isreactive::Bool = 1
	
	ip::Int64 = NG+1
	iT::Int64 = ip+1 
	
	p::Float64 = 1.0*ufac"bar"
	T::Float64 = 923.15*ufac"K"
	Tamb::Float64 = 298.15*ufac"K"

	gn::Dict{Int, Symbol} 	= kinpar.gn # names and fluid indices
	gni::Dict{Symbol, Int}  = kinpar.gni # inverse names and fluid indices
	Fluids::Vector{FluidProps} = kinpar.Fluids # fluids and respective properties in system
	
	#m::Vector{Float64} = [2.0,6.0,12.0,21.0]*ufac"g/mol"
	m::Vector{Float64} = let
		m=zeros(Float64, NG)
		for i=1:NG
			m[i] = Fluids[i].MW
		end
		m
	end
	
	X0::Vector{Float64} = let
		x=zeros(Float64, NG)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		#x[gni[:N2]] = 1.0
		x/sum(x)
	end # inlet composition

	
	mmix0::Float64 = sum(X0 .* m)
	W0::Vector{Float64} = @. m*X0/mmix0

	#mfluxin::Float64=0.01*ufac"kg/(m^2*s)"
	mfluxin::Float64 = 0.05*ufac"kg/(m^2*s)"
	nfluxin::Float64 = mfluxin/mmix0
	#mfluxin::Float64=nfluxin*mmix0

	# VitraPor data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0
	ϵ::Float64=0.33 # porosity, VitraPor sintetered filter class 0
	perm::Float64=1.23e-10*ufac"m^2" # perm. of porous medium, use in Darcy Eq.
	γ_τ::Float64=ϵ^1.5 # constriction/tourtuosity factor
	# Solid (non-porous) Borosilica glass (frit material)
	rhos::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	lambdas::Float64=1.13*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
	
end

ModelData(;ng=XuFroment.ng, kwargs...) = ModelData{ng}(;kwargs...)

FixedBed.ngas(::ModelData{NG}) where NG = NG
end;

# ╔═╡ 3bb2deff-7816-4749-9f1e-c1e451372b1e
function reaction(f,u,node,data)
	(;m,ip,T,isreactive)=data
	ng=ngas(data)

	if node.region == 2 && isreactive # catalyst layer
		(;T,lcat,kinpar,Fluids,gni)=data
		(;rni,nuij)=kinpar
		
		pi = MVector{ng,eltype(u)}(undef)
		for i=1:ng
            pi[i] = u[ip]*u[i]
		end

		rf = ramp(node.time; du=(0,1), dt=(2.0,3.0))
		RR = @inline -lcat*ri(data,T,pi)*rf	
        #RR = @inline -lcat*ri(data,T,pi)
		for i=1:ng
			f[i] = zero(eltype(u))
			for j=1:nreac(kinpar)
				f[i] += nuij[(j-1)*ng+i] * RR[j] * m[i]
			end			
		end
	end
	
	for i=1:ng
		f[ng] += u[i]
	end
	f[ng] = f[ng] - 1.0
end

# ╔═╡ 4af1792c-572e-465c-84bf-b67dd6a7bc93
function storage(f,u,node,data)
	(;ip,iT,m,ϵ,rhos,cs)=data
	ng=ngas(data)

	c = u[ip]/(ph"R"*u[iT])
	#c = u[ip]/(ph"R"*T)
	for i=1:ng
		f[i]=c*u[i]*m[i]*ϵ
	end
	
	# total pressure
	mmix = zero(eltype(u))
	cpmix = 0.0
	for i=1:ng
		mmix += u[i]*m[i]
		@inline cpmix += heatcap_gas(data.Fluids[i], u[iT])*u[i]
	end
	
	f[ip] = mmix*c*ϵ

	# solid heat capacity is 4 orders of magnitude larger than gas phase heat cap
	#f[iT] = u[iT] * (rhos*cs*(1-ϵ) + cpmix*c*ϵ)
	f[iT] = u[iT] * (rhos*cs*(1-ϵ) + cpmix*c*ϵ) / 100
end

# ╔═╡ 5f88937b-5802-4a4e-81e2-82737514b9e4
function bcond(f,u,bnode,data)
	(;p,ip,iT,T,Tamb,mfluxin,X0,W0,m,mmix0)=data
	ng=ngas(data)

	r_mfluxin = mfluxin*ramp(bnode.time; du=(0.01,1), dt=(0.0,1.0))
	@inline r_hfluxin = mfluxin/mmix0 * enthalpy_mix(data.Fluids, Tamb+200, X0) * ramp(bnode.time; du=(0.0,1), dt=(1.0,2.0))
	
	
	if dim==2		
		for i=1:(ng-1)		
			boundary_neumann!(f,u,bnode, species=i,region=Γ_in,value=r_mfluxin*W0[i])
		end		

		boundary_neumann!(f,u,bnode, species=iT,region=Γ_in,value=r_hfluxin)
		
		boundary_neumann!(f,u,bnode, species=ip, region=Γ_in, value=r_mfluxin)
		boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_out,value=p)

	else
		boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_left,value=T)
		for i=1:(ng-1)
		boundary_neumann!(f,u,bnode, species=i,region=Γ_left,value=r_mfluxin*W0[i])
		end
		boundary_neumann!(f,u,bnode, species=ip, region=Γ_left, value=r_mfluxin)
		boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_right,value=p)
		boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_right,value=T)
	end

end

# ╔═╡ c29f9187-e79c-4e56-8063-c76c98839523
function darcyvelo(u,data,mu)
	(;ip,perm) = data

	#μ = 2.0e-5*ufac"Pa*s"
	-perm/mu*(u[ip,1]-u[ip,2])	
end

# ╔═╡ 389a4798-a9ee-4e9c-8b44-a06201b4c457
function boutflow(f,u,edge,data)
	(;iT,T,ip,m)=data
	ng=ngas(data)

	k=outflownode(edge)

	pout = u[ip,k]
	#cout = pout/(ph"R"*T)
	cout = pout/(ph"R"*u[iT,k])
	X = MVector{ng,eltype(u)}(undef)
	
	for i=1:ng
		X[i] = u[i,k]
	end
	@inline mumix, _ = dynvisc_thermcond_mix(data, u[iT,k], X)
	v = darcyvelo(u,data,mumix)
	
	for i=1:(ng-1)
		# specify flux at boundary
		# f[i] = darcyvelo(u,data) * cout*u[i,k]*m[i]
		f[i] = v*cout*u[i,k]*m[i]
	end

	#@inline f[iT] = darcyvelo(u,data) *cout * enthalpy_mix(data.Fluids, u[iT,k], X) * ramp(edge.time; du=(0.0,1.0), dt=(1.0,2.0))
	@inline f[iT] = v *cout * enthalpy_mix(data.Fluids, u[iT,k], X) * ramp(edge.time; du=(0.0,1.0), dt=(1.0,2.0))
	
end

# ╔═╡ f4dba346-61bc-420d-b419-4ea2d46be6da
function D_matrix!(data, D, T, p)
	(;m,γ_τ)=data
	ng=ngas(data)
	@inbounds for i=1:(ng-1)
		for j=(i+1):ng
			Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
			#Dji *= m[i]*m[j]
			Dji *= m[i]*m[j]*γ_τ # porosity corrected eff. diffusivity
			D[j,i] = Dji
			D[i,j] = Dji
		end
	end
end

# ╔═╡ 5547d7ad-dd58-4b00-8238-6e1abb32874e
function flux(f,u,edge,data)
	(;m,ip,iT,T)=data
	ng=ngas(data)
		
	F = MVector{ng-1,eltype(u)}(undef)
	X = MVector{ng,eltype(u)}(undef)
	W = MVector{ng,eltype(u)}(undef)
	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)
	D = MMatrix{ng,ng,eltype(u)}(undef)

	pm = 0.5*(u[ip,1]+u[ip,2])
	Tm = 0.5*(u[iT,1]+u[iT,2])
	#c = pm/(ph"R"*T)
	c = pm/(ph"R"*Tm)
	
	δp = u[ip,1]-u[ip,2]
	
	mmix = zero(eltype(u))
	for i=1:ng
		X[i] = 0.5*(u[i,1]+u[i,2])
		mmix += X[i]*m[i]
	end

	@inline D_matrix!(data, D, Tm, pm)
	@inline mumix, lambdamix = dynvisc_thermcond_mix(data, Tm, X)
	
	rho = c*mmix
	#v = darcyvelo(u,data)
	v = darcyvelo(u,data,mumix)
	
	f[ip] = -rho*v

	for i=1:ng
		W[i] = m[i]/mmix*X[i]
	end
	
	@inbounds for i=1:(ng-1)
		for j=1:(ng-1)
			M[i,j] = zero(eltype(u))
		end
		 for j=1:ng
			if i != j
				M[i,i] -= W[j]/D[i,j]
				if j == ng
					 for k=1:(ng-1)
						M[i,k] -= W[i]/D[i,ng]
					end
				else
					M[i,j] += W[i]/D[i,j]
				end
			end
		end
		F[i] = ( u[i,1]-u[i,2] + (X[i]-W[i])*δp/pm )*c/mmix
	end				

	@inline inplace_linsolve!(M,F)

	hf_conv = zero(eltype(u))
	@inbounds for i=1:(ng-1)
		f[i] = -(F[i] + c*X[i]*m[i]*v)
	end
	#@inline hf_conv = f[ip] * enthalpy_mix(data.Fluids, Tm, X) / mmix
	@inline hf_conv = f[ip] * enthalpy_mix(data.Fluids, Tm, X) * ramp(edge.time; du=(0.0,1), dt=(1.0,2.0)) / mmix
	
	
	# TODO: include λbed
	Bp,Bm = fbernoulli_pm(hf_conv/1.0/Tm)
	f[iT]= 1.0*(Bm*u[iT,1]-Bp*u[iT,2])
	#f[iT]= u[iT,1]-u[iT,2]
	
end

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
begin
	
	if dim == 1
		mygrid=grid1D()
		strategy = nothing
		times=[0,10]
	elseif dim == 2
		mygrid=grid2D()
		strategy = nothing
		times=[0,5.0]
	else
		mygrid=grid3D()
		strategy = GMRESIteration(UMFPACKFactorization())
		#times=[0,1.2]
		times=[0,20.0]
	end
	mydata=ModelData()
	(;p,ip,Tamb,iT,X0)=mydata
	ng=ngas(mydata)
	
	sys=VoronoiFVM.System( 	mygrid;
							data=mydata,
							flux=flux,
							reaction=reaction,
							storage=storage,
							bcondition=bcond,
							boutflow=boutflow,
							outflowboundaries=
								[dim == 2 ? Γ_out : Γ_right],
							assembly=:edgewise
							)
	
	#enable_species!(sys; species=collect(1:(ng+1))) # gas phase species pi & ptotal	
	enable_species!(sys; species=collect(1:(ng+2))) # gas phase species xi, ptotal & T
	
	inival=unknowns(sys)

	inival[ip,:].=p
	inival[iT,:].=Tamb
	for i=1:ng
		inival[i,:] .= 1.0/ng
		#inival[i,:] .= X0[i]
	end

	control = SolverControl(strategy, sys;)
		#control.Δt=1.0e-3
		control.Δt_min=1.0e-6
		control.Δt_max=1.0
		control.handle_exceptions=true
		control.Δu_opt=1000.0
	function post(sol,oldsol, t, Δt)
		@info "t= "*string(round(t,sigdigits=2))*"\t Δt= "*string(round(Δt,sigdigits=2))
	end

	if RunSim
		solt=solve(sys;inival=inival,times,control,post)
	end
end;

# ╔═╡ f798e27a-1d7f-40d0-9a36-e8f0f26899b6
@bind t Slider(solt.t,show_value=true,default=solt.t[end])

# ╔═╡ 5588790a-73d4-435d-950f-515ae2de923c
sol = solt(t);

# ╔═╡ 99b59260-7651-45d0-b364-4f86db9927f8
let
	#vis=GridVisualizer(layout=(1,3), resolution=(700,300))
	(;iT)=mydata
	scalarplot(mygrid, sol[iT,:] .- 273.15)
end

# ╔═╡ 111b1b1f-51a5-4069-a365-a713c92b79f4
let
	(;ip,p,gn) = mydata
	ng=ngas(mydata)
	if dim == 1
		cols = distinguishable_colors(ng, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
		pcols = map(col -> (red(col), green(col), blue(col)), cols)
		vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300))
		for i=1:(ng-1)
			scalarplot!(vis, mygrid, sol[i,:], clear=false, color=pcols[i],label=gn[i])
		end
	elseif dim == 2
		vis=GridVisualizer(layout=(1,3), resolution=(700,300))
		scalarplot!(vis[1,1], mygrid, sol[1,:], clear=false, title="x1")
		scalarplot!(vis[1,2], mygrid, sol[2,:], clear=false, title="x2")
		scalarplot!(vis[1,3], mygrid, sol[3,:], clear=false, title="x3")
	else
		vis=GridVisualizer(layout=(3,1), resolution=(400,1200), outlinealpha=0.0)
		scalarplot!(vis[1,1], mygrid, sol[1,:], clear=false, label="x1")
		scalarplot!(vis[2,1], mygrid, sol[2,:], clear=false, color=:red, label="x2")
		scalarplot!(vis[3,1], mygrid, sol[3,:], clear=false, color=:blue, label="x3")
	end
	
	reveal(vis)
end

# ╔═╡ de69f808-2618-4add-b092-522a1d7e0bb7
let
	mydata = ModelData()
	(;p,m,ip,iT,Tamb,mfluxin) = mydata
	ng = ngas(mydata)
	mmix = []
	for j in 1:length(sol[1,:])
		_mmix=0
		for i=1:ng
			_mmix += sol[i,j]*m[i]
		end
		push!(mmix, _mmix)
	end
	
	ps = sol[ip,:]
	Ts = sol[iT,:]
	#rho = @. ps * mmix /(ph"R"*T)
	rho = @. ps * mmix /(ph"R"*Ts)
	
	if dim == 1
		vis=GridVisualizer(legend=:lt, resolution=(600,400), layout = (2, 1))

		scalarplot!(vis[1,1], mygrid, sol[iT,:]/Tamb, clear=false, label="T / Tamb")

		p0 = sol[ip,1]
		rho0 = @. p0 * mmix[1] /(ph"R"*T)
		scalarplot!(vis[2,1], mygrid, rho/rho0, clear=false, label="Rho / Rho0")
		scalarplot!(vis[2,1], mygrid, rho0./rho, clear=false, color=:red, label="v / v0")
		scalarplot!(vis[2,1], mygrid, ps/p0, clear=false, color=:blue, label="p / p0")
	elseif dim == 2
		vis=GridVisualizer(layout=(2,2), resolution=(660,660))
		scalarplot!(vis[1,1], mygrid, ps, title="Total Pressure")
		scalarplot!(vis[1,2], mygrid, rho, title="Total Density")
		nf = nodeflux(sys, sol)
		massflux = nf[:,ip,:]
		scalarplot!(vis[2,1], mygrid, massflux[1,:]./rho, title="Velocity - X")
		scalarplot!(vis[2,2], mygrid, massflux[2,:]./rho, title="Velocity - Y")
	else
		vis=GridVisualizer(layout=(2,1), resolution=(400,800), outlinealpha=0.0)
		scalarplot!(vis[1,1], mygrid, ps, title="Total Pressure")
		scalarplot!(vis[2,1], mygrid, rho, title="Total Density")
		
	end
	reveal(vis)
end

# ╔═╡ e000c100-ee46-454e-b049-c1c29daa9a56
let
	L=len(mygrid)
	(;mfluxin,mmix0,lcat,T,p,X0,gni) = mydata
	c0 = p*X0/(ph"R"*T)
	rho0 = p*mmix0/(ph"R"*T)
	v0 = mfluxin / rho0	
	RR = -lcat*ri(mydata,T,p*X0)
	tau = L/v0
	Da = maximum(RR)/c0[gni[:H2]] * tau
end

# ╔═╡ f690c19f-22e3-4428-bc1c-3ed7d1646e71
let	
	(;mfluxin,mmix0,X0,p,T,m,dp,Fluids) = mydata
	ng = ngas(mydata)
	rho0 = p*mmix0/(ph"R"*T)
	v0 = mfluxin / rho0

	ηf, λf = dynvisc_thermcond_mix(mydata, T, X0)
    
    cf = heatcap_mix(Fluids, T, X0)
	
	Re = v0*rho0*dp/ηf # Reynolds number
	#Pr = cf*ηf/λf # Prandtl number
	#Pe = u0*ρf*cf*d/λf # Peclet
	#Re,Pr,Pe
end

# ╔═╡ 5bbe72b2-2f80-4dae-9706-7ddb0b8b6dbe
let
	(;dp,T,p) = mydata
	σs = Dict(:H2 => 289*ufac"pm", :CH4 => 380*ufac"pm", :H2O => 265*ufac"pm", :N2 => 364*ufac"pm", :CO => 376*ufac"pm", :CO2 => 330*ufac"pm")

	Kn = ph"k_B"*T/(sqrt(2)*pi*σs[:H2O]^2*p*dp)	
end

# ╔═╡ 37b5908c-dd4e-4fb8-9d5b-68402493e10d
function DiffCoeffsMass(ng,m)
	k=0
	D=zeros(Float64, ng,ng)
	for i=1:(ng-1)
		for j=(i+1):ng
			k +=1
			#Dji = k*1.0e-5*ufac"m^2/s"
			Dji = k*1.0e-5*ufac"m^2/s"
			Dji *= m[i]*m[j]
			D[j,i] = Dji
			D[i,j] = Dji
		end			
	end
	D
end

# ╔═╡ 2cbd3e87-289c-47a2-b837-10133974ae82
function DiffCoeffs(data)
	(;m,γ_τ,T,p)=data
	ng=ngas(data)
	D=zeros(Float64, ng,ng)
	for i=1:(ng-1)
		for j=(i+1):ng

			Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
			Dji *= γ_τ # porosity corrected eff. diffusivity
			D[j,i] = Dji
			D[i,j] = Dji
		end			
	end
	D
end

# ╔═╡ ae8c7993-a89f-438a-a72a-d4a0c9a8ce57
let
	L=len(mygrid)
	(;mfluxin,mmix0,p,T) = mydata
	ng = ngas(mydata)
	rho0 = p*mmix0/(ph"R"*T)
	v0 = mfluxin / rho0
	D = DiffCoeffs(mydata)
	Pe = L*v0 / minimum(D[D.>0])
end

# ╔═╡ 1224970e-8a59-48a9-b0ef-76ed776ca15d
function checkinout(sys,sol)	
	tfact=TestFunctionFactory(sys)
	if dim == 2
		tf_in=testfunction(tfact,[Γ_out],[Γ_in])
		tf_out=testfunction(tfact,[Γ_in],[Γ_out])
	else
		tf_in=testfunction(tfact,[Γ_right],[Γ_left])
		tf_out=testfunction(tfact,[Γ_left],[Γ_right])
	end

	(;in=integrate(sys,tf_in,sol),out=integrate(sys,tf_out,sol) )
end

# ╔═╡ e29848dd-d787-438e-9c32-e9c2136aec4f
checkinout(sys,sol)

# ╔═╡ 7c7d2f10-d8d2-447e-874b-7be365e0b00c
let
	(;gn,m,nfluxin,X0) = mydata
	in_,out_=checkinout(sys,sol)
	println("Molar species in/outlfows obtained from Testfunction integration:")
	for i = 1:ngas(mydata)
		println(string(gn[i])*" IN:\t"* string(round(in_[i]/m[i]/ufac"mol/hr",sigdigits=4)) * "\t OUT: " * string(round(out_[i]/m[i]/ufac"mol/hr",sigdigits=4)) *"\t mol/hr")
	end
	println("\nMolar species in/outlfows as specified")
	for i = 1:ngas(mydata)
		println(string(gn[i])*" IN:\t"* string(round(nfluxin*X0[i]/ufac"mol/hr",sigdigits=5)) )
	end	
end

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─6da83dc0-3b0c-4737-833c-6ee91552ff5c
# ╠═d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═83fa22fa-451d-4c30-a4b7-834974245996
# ╠═4dae4173-0363-40bc-a9ca-ce5b4d5224cd
# ╠═561e96e2-2d48-4eb6-bb9d-ae167a622aeb
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═4e05ab31-7729-4a4b-9c14-145118477715
# ╠═107a6fa3-60cb-43f0-8b21-50cd1eb5065a
# ╠═832f3c15-b75a-4afe-8cc5-75ff3b4704d6
# ╟─a078e1e1-c9cd-4d34-86d9-df4a052b6b96
# ╟─0fadb9d2-1ccf-4d44-b748-b76d911784ca
# ╟─b94513c2-c94e-4bcb-9342-47ea48fbfd14
# ╟─c886dd12-a90c-40ab-b9d0-32934c17baee
# ╟─8f2549f4-b0a6-440f-af94-6880e0814dc2
# ╟─78589a1e-2507-4279-ba42-1aaec90d87d0
# ╠═5547d7ad-dd58-4b00-8238-6e1abb32874e
# ╠═3bb2deff-7816-4749-9f1e-c1e451372b1e
# ╠═4af1792c-572e-465c-84bf-b67dd6a7bc93
# ╠═5f88937b-5802-4a4e-81e2-82737514b9e4
# ╠═389a4798-a9ee-4e9c-8b44-a06201b4c457
# ╠═927dccb1-832b-4e83-a011-0efa1b3e9ffb
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═5588790a-73d4-435d-950f-515ae2de923c
# ╠═e29848dd-d787-438e-9c32-e9c2136aec4f
# ╠═7c7d2f10-d8d2-447e-874b-7be365e0b00c
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╠═99b59260-7651-45d0-b364-4f86db9927f8
# ╠═111b1b1f-51a5-4069-a365-a713c92b79f4
# ╠═de69f808-2618-4add-b092-522a1d7e0bb7
# ╟─68ca72ae-3b24-4c09-ace1-5e340c8be3d4
# ╟─db77fca9-4118-4825-b023-262d4073b2dd
# ╠═ae8c7993-a89f-438a-a72a-d4a0c9a8ce57
# ╟─e7497364-75ef-4bd9-87ca-9a8c2d97064c
# ╠═e000c100-ee46-454e-b049-c1c29daa9a56
# ╟─f6e54602-b63e-4ce5-b8d5-f626fbe5ae7a
# ╠═f690c19f-22e3-4428-bc1c-3ed7d1646e71
# ╟─0d507c5e-eb0f-4094-a30b-a4a82fd5c302
# ╠═5bbe72b2-2f80-4dae-9706-7ddb0b8b6dbe
# ╟─67adda35-6761-4e3c-9d05-81e5908d9dd2
# ╠═f40a8111-b5cb-40f0-8b12-f57cf59637f1
# ╠═c29f9187-e79c-4e56-8063-c76c98839523
# ╠═f4dba346-61bc-420d-b419-4ea2d46be6da
# ╠═37b5908c-dd4e-4fb8-9d5b-68402493e10d
# ╠═2cbd3e87-289c-47a2-b837-10133974ae82
# ╠═1224970e-8a59-48a9-b0ef-76ed776ca15d
