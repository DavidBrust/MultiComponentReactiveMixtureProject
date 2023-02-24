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

# ╔═╡ 404e4e80-b28d-11ed-3de7-d77a058a4923
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	
	using VoronoiFVM
	using ExtendableGrids, GridVisualize
	using LinearAlgebra
	using StaticArrays
	using PreallocationTools
	using LessUnitful
	
	using PlutoVista	
	using PlutoUI
	using PyPlot	
	using Colors

	using Revise
	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ d7bb0f6d-48b2-4a51-bbcf-151081cd8da9
begin
	const Γ_side_front = 1 # symmetry bc
	const Γ_side_right = 2 # wall bc
	const Γ_side_back = 3 # wall bc
	const Γ_side_left = 4 # symmetry bc
	const Γ_bottom = 5 # inflow bc
	const Γ_top_frit = 6 # outflow bc, uncoated porous frit 
	const Γ_top_cat = 7 # outflow bc, catalyst coated porous frit 
end;

# ╔═╡ 1f1bdae4-b7c4-4b8f-bf59-8030c67ce686
function Dmatrix!(Dmatrix, n, Fluids, T, p)
	#v = zeros(typeof(T),n,n)
	#Dmatrix[:] .= zero(eltype(Dmatrix))
	for i=1:(n-1)
		for j=(i+1):n
			Dji = binary_diff_coeff_gas(Fluids[j], Fluids[i], T, p)
			Dmatrix[j,i] = Dji
			Dmatrix[i,j] = Dji # symmetric
		end
	end
	nothing
	#Symmetric(v, :L)
end

# ╔═╡ 253dd8ea-cc65-438b-992e-9f40ae10b138
function Mmatrix!(Mmatrix, Dmatrix, DK_eff, n, γ_τ, x)
	#n=data.ng
	#D=data.γ_τ*D_matrix(data,T,p)
	#M=zeros(eltype(x), n, n)
	#Mmatrix[:] .= zero(eltype(Mmatrix))
	for i=1:n
		Mmatrix[i,i] = -1/DK_eff[i]
		for j=1:n
			if j != i
				Mmatrix[i,i] -= x[j]/(γ_τ*Dmatrix[i,j])
				Mmatrix[i,j] = x[i]/(γ_τ*Dmatrix[i,j])
			end
		end	
	end
	nothing
end

# ╔═╡ 13350e83-5f03-4b6c-b7e8-c911b085ff1b
function DK_eff(data,T,i)
	ϕ=data.ϕ # porosity
	dp=data.dp # avg. particle size (assumed to be = pore size)
	DK=dp*ϕ^1.5/(1.0-ϕ)*sqrt(8.0*ph"R"*T/(9.0*π*data.Fluids[i].MW))
	#Bern(DK)
end

# ╔═╡ 4c1749e1-84d2-42d6-ab60-5f9deffe54f2
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

# ╔═╡ 2a4d8f57-c159-4507-a61e-55dcd235c437
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

# ╔═╡ 684c526d-25bc-4dcf-85a1-2434e5d185dc
function project!(vh, edge::VoronoiFVM.Edge, vec)
    vh[1]=zero(eltype(vec))
    for i=1:size(edge.coord)[1]
        vh[1]+=(edge.coord[i,edge.node[2]]-edge.coord[i,edge.node[1]])*vec[i]
    end
    nothing
end

# ╔═╡ caf7503f-eb32-4132-957e-0d024c14fbf5
function DK_eff!(DK_eff,i,ϕ,dp,Fluid,T,)
	DK_eff[i]=dp*ϕ^1.5/(1.0-ϕ)*sqrt(8.0*ph"R"*T/(9.0*π*Fluid.MW))
	nothing
end

# ╔═╡ 6185cf12-879d-455e-8462-54bf5f37c34a
function reaction(f,u,node,data)
	(;ng,ip,pt) = data	
	
	if node.region == 2 && data.isreactive # catalyst layer
		#ng=data.gni # gas species indices from names
		nr=data.kinpar.rni # reaction indices from names
		iT=data.iT
		
		pi = u[1:ng]./ufac"bar"
		# negative sign: sign convention of VoronoiFVM: source term < 0 
		# ri returns reaction rate in mol/(h gcat)
		RR = -data.mcats*ri(data.kinpar,u[iT],pi)*ufac"mol/hr"*ufac"1/g"
		# reactions in S3P kinetics model
		# R1: CO + H2O = CO2 + H2
		# R2: CH4 + 2 H2O = CO2 + 4 H2
		# R3: CH4 + H2O = CO + 3 H2

		for i=1:ng
			f[i] = sum(data.kinpar.nuij[i,:] .* RR)
		end
		
		# temperature eq. / heat source
		ΔHi=data.kinpar.ΔHi
		f[iT] = -(RR[nr["R1"]]*ΔHi["R1"]+RR[nr["R2"]]*ΔHi["R2"] +RR[nr["R3"]]*ΔHi["R3"])
	end
	
	# ∑xi = 1
	pt[1]=zero(eltype(u))
	for i=1:ng
		pt[1] += u[i]
	end
	f[ip]=u[ip]-pt[1]
	
end

# ╔═╡ 1d07dbba-ed06-418b-a12f-6697120b5a87
function bottom(f,u,bnode,data)
	if bnode.region==Γ_bottom # bottom boundary
		(;ng,iT,Eps_ir,Tamb,u0,X0,pn,Tn) = data

		f[iT] = Eps_ir*ph"σ"*(u[iT]^4 - Tamb^4)
		
		for i=1:ng
			# specify flux at boundary: flow velocity is normal to bot boundary	
			f[i] = -u0*X0[i]*pn/(ph"R"*Tn)
			
		end
	end
end

# ╔═╡ 6cf7645b-9e1a-4b42-87ba-8d46f3a4456f
function side(f,u,bnode,data)
	# side wall boundary condition
	(;iT,α_w,Tamb) = data
	boundary_robin!(f,u,bnode;species=iT,region=[Γ_side_back,Γ_side_right], factor=α_w, value=Tamb*α_w)	
end

# ╔═╡ 825e6eac-7f14-4dd8-aa10-03afc65f1d77
md"""
## Activate / Deactivate Chemistry
"""

# ╔═╡ 714addba-0c07-4ee3-8f14-bed0ef2ec7e1
md"""
Upon activating the chemical reaction, the solver fails to solve the system with "Converge Error".
"""

# ╔═╡ 0f1771f7-cafc-4d7d-ba35-b9c172daf001
md"Turn on chemical reactions $(@bind isreactive CheckBox())"

# ╔═╡ f9906982-cbba-494a-9c42-f74ab0f9db6e
isreactive

# ╔═╡ 2ec69cd6-2614-473a-b7ad-3840085a430b
Base.@kwdef mutable struct ModelData <:AbstractModelData
	#S::Tv = 1.0
	
	# catalyst / chemistry data
	# kinetic parameters, S3P="simple 3 parameter" kinetics fit to UPV lab scale experimental data
	# kinpar::AbstractKineticsData = S3P
	kinpar::AbstractKineticsData = XuFroment1989
	
	# number of gas phase species
	ng::Int64		 		= kinpar.ng
	#ng::Int64		 		= 1
	# names and fluid indices
	gn::Dict{Int, String} 	= kinpar.gn
	#gn::Dict{Int, String} 	= Dict(1=>"N2")
	# inverse names and fluid indices
	gni::Dict{String, Int}  = kinpar.gni
	# fluids and respective properties in system
	Fluids::Vector{FluidProps} = kinpar.Fluids
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = let
		x=zeros(Float64, ng)
		x[gni["H2"]] = 1.0
		x[gni["CO2"]] = 1.0
		x/sum(x)
	end # inlet composition
	#X0::Vector{Float64} = [1.0]
	
	# volume specific cat mass loading, UPV lab scale PC reactor
	mcats::Float64 =1234.568*ufac"kg/m^3"
	isreactive::Bool = 1
	#isreactive::Bool = 0

	ip::Int64=ng+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable
	#iT::Int64=1 # index of Temperature variable

		
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

	MWin::Float64 = molarweight_mix(Fluids, X0)
	mdotin::Float64=MWin*Qflow*pn/(ph"R"*Tn)*ufac"kg/s"

	
	Tin::Float64=298.15*ufac"K" # inlet temperature
	Tamb::Float64=Tin # ambient temperature
	#Tin::Float64=600.0*ufac"K" # inlet temperature
	p::Float64=1.0*ufac"atm" # reactor pressure

	# u0::Float64=Qflow/(Ac*ϕ)*ufac"m/s" # mean superficial velocity
	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	utop::Float64=u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate

	# intermediate variables to avoid allocations
	X::Vector{Float64} = ones(Float64, ng)
	cmix::Vector{Float64} = ones(Float64, 1)*ufac"J/mol"
	cf::Vector{Float64} = ones(Float64, ng)*ufac"J/mol"
	μmix::Vector{Float64} = ones(Float64, 1)*ufac"Pa*s"
	μf::Vector{Float64} = ones(Float64, ng)*ufac"Pa*s"
	λmix::Vector{Float64} = ones(Float64, 1)*ufac"W/(m*K)"
	λf::Vector{Float64} = ones(Float64, ng)*ufac"W/(m*K)"
	kbed::Vector{Float64} = ones(Float64, 1)
	uvec::Vector{Float64} = zeros(Float64,3)
	vh::Vector{Float64} = ones(Float64, 1)*ufac"m/s"
	DK_eff::Vector{Float64} = ones(Float64, ng)*ufac"m^2/s"
	Dmatrix::Matrix{Float64}=ones(Float64,ng,ng)*ufac"m^2/s"
	Mmatrix::Matrix{Float64}=ones(Float64,ng,ng)*ufac"m^2/s"
	
	pt::Vector{Float64} = ones(Float64, 1)*ufac"Pa"
end;

# ╔═╡ 273ea122-3adb-42db-89c0-b85ca5c16756
function mole_frac!(node,data::ModelData,X,u::VoronoiFVM.EdgeUnknowns)
	n=data.ng
	sump = zero(eltype(u))
	for i=1:n
		sump += 0.5*(u[i,1]+u[i,2])
	end
	for i=1:n
		X[i] = 0.5*(u[i,1]+u[i,2])/sump
	end
	#X .= X / sump
	nothing
end

# ╔═╡ 73b5ad9b-49ed-48b5-9cac-be949593229d
function mole_frac!(node,data::ModelData,X,u::VoronoiFVM.BNodeUnknowns)
	n=data.ng
	sump = zero(eltype(u))
	for i=1:n
		sump += u[i]
	end
	for i=1:n
		X[i] = u[i]/sump
	end
	#X .= X / sump
	nothing
end

# ╔═╡ 292ea1d4-2e2f-4c58-aa5f-aeff0e385d57
function flux_diffcache(f,u,edge,data)
	(;ng,ip,iT,Fluids,k,ϕ,dp,γ_τ,λs,X,cf,cmix,μf,μmix,λf,λmix,kbed,uvec,vh,Mmatrix,Dmatrix,DK_eff,F) = data
	
	#F= @MArray zeros(eltype(u), ng)
	F = get_tmp(F, u)

	pk,pl = u[ip,1],u[ip,2]
	δp = pk-pl
	p = 0.5*(pk+pl)

	T=0.5*(u[iT,1]+u[iT,2])
	mole_frac!(edge,data,X,u)	
	
	FixedBed.heatcap_mix!(cmix, cf, Fluids, T, X)
	FixedBed.dynvisc_thermcond_mix!(μmix, μf, λmix, λf, Fluids, T, X)
	FixedBed.kbed!(kbed,ϕ,λs,λmix)
	λbed=kbed[1]*λmix[1]

	# Darcy flow
	ud=-k/μmix[1] * δp	
	uvec[3]=ud # 3D
	
	project!(vh, edge, uvec)
	
	# convective enthalpy flux
	conv=vh[1]*p/(ph"R"*T)*cmix[1]/λbed
	Bp,Bm = fbernoulli_pm(conv)
	# thermal energy flux
	f[iT]= λbed*(Bm*(u[iT,1]-data.Tamb)-Bp*(u[iT,2]-data.Tamb))

	
	for i=1:ng
		DK_eff!(DK_eff,i,ϕ,dp,Fluids[i],T)
		bp,bm=fbernoulli_pm(vh[1]/DK_eff[i])
		
		F[i] = -(bm*u[i,1]-bp*u[i,2])/(ph"R"*T)
	end
	
	
	Dmatrix!(Dmatrix, ng, Fluids, T, p)
    Mmatrix!(Mmatrix, Dmatrix, DK_eff, ng, γ_τ, X)

	f[1:ng] = Mmatrix \ F

		
	#f[ip] via reaction: ∑pi = p
end

# ╔═╡ b88a1970-462d-4d48-b4b4-c7e83323b8c8
function flux(f,u,edge,data)
	(;ng,ip,iT,Fluids,k,ϕ,dp,γ_τ,λs,X,cf,cmix,μf,μmix,λf,λmix,kbed,uvec,vh,Mmatrix,Dmatrix,DK_eff,) = data
	
	F= zeros(eltype(u), ng)

	pk,pl = u[ip,1],u[ip,2]
	δp = pk-pl
	p = 0.5*(pk+pl)

	T=0.5*(u[iT,1]+u[iT,2])
	mole_frac!(edge,data,X,u)	
	
	FixedBed.heatcap_mix!(cmix, cf, Fluids, T, X)
	FixedBed.dynvisc_thermcond_mix!(μmix, μf, λmix, λf, Fluids, T, X)
	FixedBed.kbed!(kbed,ϕ,λs,λmix)
	λbed=kbed[1]*λmix[1]

	# Darcy flow
	ud=-k/μmix[1] * δp	
	uvec[3]=ud # 3D
	
	project!(vh, edge, uvec)
	
	# convective enthalpy flux
	conv=vh[1]*p/(ph"R"*T)*cmix[1]/λbed
	Bp,Bm = fbernoulli_pm(conv)
	# thermal energy flux
	f[iT]= λbed*(Bm*(u[iT,1]-data.Tamb)-Bp*(u[iT,2]-data.Tamb))

	
	for i=1:ng
		DK_eff!(DK_eff,i,ϕ,dp,Fluids[i],T)
		bp,bm=fbernoulli_pm(vh[1]/DK_eff[i])
		
		F[i] = -(bm*u[i,1]-bp*u[i,2])/(ph"R"*T)
	end
	
	
	Dmatrix!(Dmatrix, ng, Fluids, T, p)
    Mmatrix!(Mmatrix, Dmatrix, DK_eff, ng, γ_τ, X)

	f[1:ng] = Mmatrix \ F

		
	#f[ip] via reaction: ∑pi = p
end

# ╔═╡ 17ae8bd5-b948-46a9-8f88-d7a016f6f4ba
function top(f,u,bnode,data)
	# top boundaries (cat layer & frit)
	if bnode.region==Γ_top_frit || bnode.region==Γ_top_cat 
		(;ng,ip,iT,G_lamp,Eps_ir,Tau_quartz,Tamb,X,cmix,cf,Fluids,utop,) = data
		
		flux_rerad = Eps_ir*ph"σ"*(u[iT]^4 - Tamb^4)

		
		#X=zeros(eltype(u), ng)
		mole_frac!(bnode,data,X,u)
		FixedBed.heatcap_mix!(cmix, cf, Fluids, u[iT], X)
		#cf=heatcap_mix(data.Fluids, u[iT], X)
		
		flux_convec=utop*u[ip]/(ph"R"*u[iT])*cmix[1]*(u[iT]-Tamb)

		abs = 0.0
		if bnode.region==Γ_top_frit
			abs=data.Abs_lamp_frit
		else # catalyst layer
			abs=data.Abs_lamp_cat
		end
		f[iT] = -abs*Tau_quartz*G_lamp + flux_rerad + flux_convec
		
		
		# flow velocity is normal to top boundary
		for i=1:ng
			f[i] = utop*u[i]/(ph"R"*u[iT])
		end
	end
end

# ╔═╡ 98175d54-b6db-4eb5-bcad-be67a060ff61
function bcond(f,u,bnode,data)
	top(f,u,bnode,data)
	bottom(f,u,bnode,data)
	side(f,u,bnode,data)
end

# ╔═╡ b7b0361e-cfc3-43c3-b6ba-e429a8aee40d
data=ModelData();

# ╔═╡ 53d678aa-ecc4-41b1-9ba9-26fc0c006603
md"""
Cutplane at ``z=`` $(@bind zcut Slider(range(0.0,data.h,length=101),default=data.h,show_value=true)) m
"""

# ╔═╡ fcd5cd1b-8b10-4326-b177-ee703dd030dd
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

# ╔═╡ 7b4aab28-0cd6-4615-9e29-e9110766db98
let
	vis=GridVisualizer()
	gridplot!(vis, prism_sq(), zplane=zcut,show=true)
end

# ╔═╡ 0e23757f-81ca-4145-bd81-44bedca2cd4c
function main()
	data = ModelData()
	data.isreactive = isreactive
	grid = prism_sq()
	iT = data.iT
	ng = data.ng

	sys=VoronoiFVM.System( grid;
						data=data,
						flux=flux,
						reaction=reaction,
						bcondition=bcond,
						)
	enable_species!(sys; species=collect(1:(ng+2)))

	
	inival=unknowns(sys)
	inival[:,:].=1.0*data.p
	for i=1:ng
		inival[i,:] .*= data.X0[i]
	end
	inival[iT,:] .= data.Tamb

	sol=solve(sys;inival,)

	function pre(sol,par)
		ng=data.ng
		iT=data.iT
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

		utop_calc=ntop*ph"R"*Tavg/(1.0*ufac"bar")/(data.Ac)
		utops=[data.utop, utop_calc]*ufac"m/s"
		data.utop = minimum(utops) + par*(maximum(utops)-minimum(utops))
		
		# specific catalyst loading
		mcats=[10.0, 1300.0]*ufac"kg/m^3"
		data.mcats= minimum(mcats) + par*(maximum(mcats)-minimum(mcats))

		# irradiation flux density
		G_lamp=[1.0, 125.0]*ufac"kW/m^2"
		data.G_lamp= minimum(G_lamp) + par*(maximum(G_lamp)-minimum(G_lamp))
	end
	
	control=SolverControl(;
					  handle_exceptions=true,
					  Δp_min=5.0e-4,					  
					  Δp=0.1,
					  Δp_grow=1.2,
					  Δu_opt=100000.0, # large value, due to unit Pa of pressure?
					  )
	

	
	sol=solve(sys;inival,embed=[0.0,1.0],pre,control)


	

	sol,sys,grid,data
end

# ╔═╡ c6552a92-c4eb-4d3a-84fc-5518ea0ad48f
begin
	sol_,sys,grid,data_embed=main();
	if sol_ isa VoronoiFVM.TransientSolution
		sol = copy(sol_(sol_.t[end]))
	else
		sol = copy(sol_)
	end
end;

# ╔═╡ 4dc8b256-3725-4dce-ae6c-e251f0c35d37
let
	vis=GridVisualizer(Plotter=PyPlot)
	scalarplot!(vis, grid, sol[data.iT,:].- 273.15, show=true,)
end

# ╔═╡ Cell order:
# ╠═404e4e80-b28d-11ed-3de7-d77a058a4923
# ╠═d7bb0f6d-48b2-4a51-bbcf-151081cd8da9
# ╠═7b4aab28-0cd6-4615-9e29-e9110766db98
# ╟─53d678aa-ecc4-41b1-9ba9-26fc0c006603
# ╠═fcd5cd1b-8b10-4326-b177-ee703dd030dd
# ╠═1f1bdae4-b7c4-4b8f-bf59-8030c67ce686
# ╠═253dd8ea-cc65-438b-992e-9f40ae10b138
# ╠═13350e83-5f03-4b6c-b7e8-c911b085ff1b
# ╠═4c1749e1-84d2-42d6-ab60-5f9deffe54f2
# ╠═2a4d8f57-c159-4507-a61e-55dcd235c437
# ╠═273ea122-3adb-42db-89c0-b85ca5c16756
# ╠═73b5ad9b-49ed-48b5-9cac-be949593229d
# ╠═684c526d-25bc-4dcf-85a1-2434e5d185dc
# ╠═caf7503f-eb32-4132-957e-0d024c14fbf5
# ╟─292ea1d4-2e2f-4c58-aa5f-aeff0e385d57
# ╠═b88a1970-462d-4d48-b4b4-c7e83323b8c8
# ╠═6185cf12-879d-455e-8462-54bf5f37c34a
# ╠═17ae8bd5-b948-46a9-8f88-d7a016f6f4ba
# ╠═1d07dbba-ed06-418b-a12f-6697120b5a87
# ╠═6cf7645b-9e1a-4b42-87ba-8d46f3a4456f
# ╠═98175d54-b6db-4eb5-bcad-be67a060ff61
# ╟─825e6eac-7f14-4dd8-aa10-03afc65f1d77
# ╟─714addba-0c07-4ee3-8f14-bed0ef2ec7e1
# ╠═0f1771f7-cafc-4d7d-ba35-b9c172daf001
# ╠═f9906982-cbba-494a-9c42-f74ab0f9db6e
# ╠═0e23757f-81ca-4145-bd81-44bedca2cd4c
# ╠═c6552a92-c4eb-4d3a-84fc-5518ea0ad48f
# ╠═4dc8b256-3725-4dce-ae6c-e251f0c35d37
# ╠═b7b0361e-cfc3-43c3-b6ba-e429a8aee40d
# ╠═2ec69cd6-2614-473a-b7ad-3840085a430b
