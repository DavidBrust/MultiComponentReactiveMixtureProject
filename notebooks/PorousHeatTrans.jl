### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 5c3adaa0-9285-11ed-3ef8-1b57dd870d6f
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	
	using VoronoiFVM
	using ExtendableGrids
	using GridVisualize
	using LessUnitful
	using PlutoVista	
	using PlutoUI

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 7d8eb6f5-3ba6-46ef-8058-1f24a0938ed1
PlutoUI.TableOfContents(title="Heat Transfer in Fixed Beds")

# ╔═╡ f353e09a-4a61-4def-ab8a-1bd6ce4ed58f
md"""
# Porous filter disc
"""

# ╔═╡ 2015c8e8-36cd-478b-88fb-94605283ac29
md"""
Specifications of porous filter disc from sintered silica glas (SiO₂): __VitraPOR P2__ (40-100 μm)
$(LocalResource("../img/filter1.png", :width => 1000))
$(LocalResource("../img/filter2.png", :width => 1000))
$(LocalResource("../img/filter3.png", :width => 1000))
"""

# ╔═╡ 03d0c88a-462b-43c4-a589-616a8870be64
md"""
# Experimental Conditions
"""

# ╔═╡ b4cee9ac-4d14-4169-90b5-25d12ac9c003
md"""
## Porous media effective thermal conductivity
"""

# ╔═╡ bbd0b076-bcc1-43a1-91cb-d72bb17d3c88
md"""
Effective thermal conductivity of porous filter frit according to:

__Zehner, P., & Schlünder, E. U. (1970).__ Wärmeleitfähigkeit von Schüttungen bei mäßigen Temperaturen. Chemie Ingenieur Technik, 42(14), 933-941. doi:10.1002/cite.330421408

Implementation follows the notation in __VDI Heat Atlas 2010, ch. D6.3 eqs. (5a-5e).__
"""

# ╔═╡ d7317b2d-e2c7-4114-8985-51979f2205ba
md"""
# Grid
"""

# ╔═╡ 2fe11550-683d-4c4b-b940-3e63a4f8a87d
function cylinder(;nref=0, r=5.0*ufac"cm", h=0.5*ufac"cm")
    step=0.1*ufac"cm"*2.0^(-nref)
    R=collect(0:step:r)
    Z=collect(0:step:h)
    grid=simplexgrid(R,Z)
    circular_symmetric!(grid)
	grid
end

# ╔═╡ 9d8c6ddc-2662-4055-b636-649565c36287
md"""
# Simulation
"""

# ╔═╡ e148212c-bc7c-4553-8ffa-26e6c20c5f47
md"""
# Thermo-physical properties
"""

# ╔═╡ 9be8cb32-88ae-43e2-b7b6-6ff8428c8f3c
#Dynamic viscosity of gases at low pressures, Pa*s
function dynvisc_gas(data, T)
	(;A,B,C,D,E) = data.DynVisc
	# VDI heat atlas 2010 D3.1 Equation (3)
	A+B*T+C*T^2+D*T^3+E*T^4 * ufac"Pa*s"
end

# ╔═╡ 534f29f2-2563-465b-9a37-bf1791ca22d9
#Thermal conductivity of gases at low pressures, W/(m*K)
function thermcond_gas(data, T)
	(;A,B,C,D,E) = data.ThermCond
	# VDI heat atlas 2010 D3.1 Equation (5)
	A+B*T+C*T^2+D*T^3+E*T^4 * ufac"W/(m*K)"
end

# ╔═╡ 7c4d4083-33e2-4b17-8575-214ef458d75d
function kbed(data)
	(;ϕ,λs,Fluid,Tin) = data
	λf=thermcond_gas(Fluid, Tin)
	B=1.25*((1.0-ϕ)/ϕ)^(10.0/9.0)
	kp=λs/λf
	N=1.0-(B/kp)
	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	1.0-sqrt(1.0-ϕ)+sqrt(1.0-ϕ)*kc
end

# ╔═╡ bd0fd520-217c-471e-89af-fd568ea70be5
#Molar heat capacity of ideal gases, J/(kg*K)
function heatcap_gas(data, T)
	(;A,B,C,D,E,F,G) = data.HeatCap
	# VDI heat atlas 2010 D3.1 Equation (10)
	T_ApT = (T/(A+T))
	(B+(C-B)*T_ApT^2*(1- (A/(A+T))*(D+E*T_ApT+F*T_ApT^2+G*T_ApT^3) ) ) * ph"R" / data.MW  * ufac"J/(kg*K)"
end

# ╔═╡ d4bd9a44-f96b-4c94-af47-fdb2c09eda82
function density_idealgas(data, T, p)
	p/(ph"R"*T)*data.MW*ufac"kg/m^3"
end

# ╔═╡ 4f9e10ae-8c37-4c89-9658-fb47095a1c1e
abstract type AbstractFluidProps end

# ╔═╡ aaee8aa5-a41c-4879-bea5-485f1188f9d0
abstract type AbstractPropsCoeffs end

# ╔═╡ 37cf0697-0091-49f8-bda0-9864314b2dff
Base.@kwdef mutable struct PropsCoeffs <: AbstractPropsCoeffs
    A::Float64=1.0
	B::Float64=1.0
	C::Float64=1.0
	D::Float64=1.0
	E::Float64=1.0
	F::Float64=1.0
	G::Float64=1.0
end

# ╔═╡ 29db6427-f1a0-4a2e-8a3b-82d2ad48b42b
Base.@kwdef mutable struct FluidProps <: AbstractFluidProps
	name::String="Air"
	MW::Float64=28.96*ufac"g/mol"
	HeatCap::AbstractPropsCoeffs=PropsCoeffs(
	A=2548.9320,
	B=3.5248,
	C=-0.6366,
	D=-3.4281,
	E=49.8238,
	F=-120.3466,
	G=98.8658
	)
	ThermCond::AbstractPropsCoeffs=PropsCoeffs( 
	A=-0.908e-3,
	B=0.112e-3,
	C=-0.084333e-6,
	D=0.056964e-9,
	E=-0.015631e-12
	)
	DynVisc::AbstractPropsCoeffs=PropsCoeffs(
	A=-0.01702e-5,
	B=0.79965e-7,
	C=-0.72183e-10,
	D=0.04960e-12,
	E=-0.01388e-15
	)	
end

# ╔═╡ 98063329-31e1-4d87-ba85-70419beb07e9
Base.@kwdef mutable struct ModelData
	iT::Int64=1 # index of Temperature variable
	
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	α_w::Float64=20*ufac"W/(m^2*K)" # wall heat transfer coefficient
	α_nc::Float64=10*ufac"W/(m^2*K)" # natural convection heat transfer coefficient
	
	## porous filter data
    D::Float64=10.0*ufac"cm" # disc diameter
	h::Float64=0.5*ufac"cm" # disc thickness
	d::Float64=100.0*ufac"μm" # average pore size
	Ac::Float64=pi*D^2.0/4.0*ufac"m^2" # cross-sectional area
	
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous SiO2
	λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2 	
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
	
	ϕ::Float64=0.36 # porosity, class 2
	k::Float64=2.9e-11*ufac"m^2" # permeability
	a_s::Float64=0.13*ufac"m^2/g" # specific surface area
	## END porous filter data

	## fluid data
	Qflow::Float64=2000.0*ufac"ml/minute" # volumetric feed flow rate
	Tin::Float64=298.15*ufac"K" # inlet temperature
	p::Float64=1.0*ufac"atm" # reactor pressure		
	u0::Float64=Qflow/(Ac*ϕ)*ufac"m/s" # mean superficial velocity
	# fluid properties: Air
	# values taken from VDI heat atlas 2010 chapter D3.1
	Fluid::FluidProps=FluidProps(
		name="Air",
		MW=28.96*ufac"g/mol",
		HeatCap=PropsCoeffs( # table 6
			A=2548.9320,
			B=3.5248,
			C=-0.6366,
			D=-3.4281,
			E=49.8238,
			F=-120.3466,
			G=98.8658
		),
		ThermCond=PropsCoeffs( # table 6
			A=-0.908e-3,
			B=0.112e-3,
			C=-0.084333e-6,
			D=0.056964e-9,
			E=-0.015631e-12
		),	
		DynVisc=PropsCoeffs( # table 6
			A=-0.01702e-5,
			B=0.79965e-7,
			C=-0.72183e-10,
			D=0.04960e-12,
			E=-0.01388e-15
		)
	)
	## END fluid data
	
end;

# ╔═╡ 3b3595c4-f53d-4827-918e-edcb74dd81f8
data = ModelData(;p=1.0*ufac"atm",Qflow=3400*ufac"ml/minute")

# ╔═╡ 8cd85a0e-3d11-4bcc-8a7d-f30313b31363
gridplot(cylinder(;r=data.D/2,h=data.h))

# ╔═╡ 6a92c1c7-fb51-4363-b8d0-18eeb24087a8
let
	(;Fluid)=data
	Ts=(273.15:1:573.15)
	λf = zeros(Float64, length(Ts))
	for (i,T) in enumerate(Ts)
		λf[i]=thermcond_gas(Fluid, T)
	end
	plot(Ts.-273.15, λf,xlabel="Temperature / °C",ylabel="Thermal Conductivity / W m⁻¹ K⁻¹", title=Fluid.name )
end

# ╔═╡ d725f9b9-61c4-4724-a1d9-6a04ba42499d
function main(;p=1.0*ufac"atm",Qflow=3400*ufac"ml/minute")
	data=ModelData(Qflow=Qflow,	p=p,)
	iT=data.iT

	# function return 2D velocity vector: flow upward in z-direction
    function fup(r,z)
        return 0,-data.u0
    end    
	
	function flux(f,u,edge,data)
		(;Fluid,u0,p)=data
		Tbar=0.5*(u[iT,1]+u[iT,2])
		ρf=density_idealgas(Fluid, Tbar, p)
		cf=heatcap_gas(Fluid, Tbar)
		λf=thermcond_gas(Fluid, Tbar)
		#λf=thermcond_gas(Fluid, Tin)
		λbed=kbed(data)*λf
		
		#f[iT] = λbed*(u[iT,1]-u[iT,2])
		conv=evelo[edge.index]*ρf*cf/λbed
		Bp,Bm = fbernoulli_pm(conv)
		f[iT]= λbed*(Bm*u[iT,1]-Bp*u[iT,2])
		
		#f[iT] = λbed*(u[iT,1]-u[iT,2])
		
	end

	function bcondition(f,u,bnode,data)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=1,value=data.Tamb)
		boundary_robin!(f,u,bnode;species=iT,region=1, factor=data.α_nc, value=data.Tamb*data.α_nc)
		boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
		boundary_dirichlet!(f,u,bnode;species=iT,region=3,value=data.Tamb+300.0)
	end
	

	
	grid=cylinder(;r=data.D/2,h=data.h)
	evelo=edgevelocities(grid,fup)
	
	sys=VoronoiFVM.System(grid;
                          data=data,
                          flux=flux,
    #                      reaction=pnpreaction,
    #                      #storage=pnpstorage,
                          bcondition,
                          species=[iT],
	#					  regions=[1,2],
    #                      kwargs...
                          )
	inival=unknowns(sys)
	inival[iT,:] .= map( (r,z)->(data.Tamb+300*z/data.h),grid)
	sol=solve(inival,sys)
	sys,sol,data
end

# ╔═╡ f15fd785-010c-4fda-ab4f-7947642556dd
let
	sys,sol,data=main()
	iT=data.iT
	vis=GridVisualizer()
	@. sol[iT,:] -= 273.15
	scalarplot!(vis,sys,sol;species=iT,title="T",colormap=:summer,show=true)

end

# ╔═╡ b4f7d4d5-c277-4a1d-9dff-0a35853d80ef
let
	Air=FluidProps()
	typeof(Air)
end

# ╔═╡ 453db496-4688-48ce-8e2a-1a960540f69f
# all values taken from VDI heat atlas 2010 chapter D3.1
Air_=(
	# Table 1
	MW=28.96*ufac"g/mol",
	# Table 6
	HeatCap=PropsCoeffs(
	A=2548.9320,
	B=3.5248,
	C=-0.6366,
	D=-3.4281,
	E=49.8238,
	F=-120.3466,
	G=98.8658
	),
	# Table 10
	ThermCond = PropsCoeffs( 
	A=-0.908e-3,
	B=0.112e-3,
	C=-0.084333e-6,
	D=0.056964e-9,
	E=-0.015631e-12
	),
	# Table 8
	DynVisc= PropsCoeffs(
	A=-0.01702e-5,
	B=0.79965e-7,
	C=-0.72183e-10,
	D=0.04960e-12,
	E=-0.01388e-15
)
);

# ╔═╡ d912a1ca-1b69-4ea1-baa5-69794e004693
begin
	Ac=pi*data.D^2/4
	u0 = data.Qflow / (Ac*data.ϕ)
	# fluid properties: Air
	ρf = density_idealgas(Air_, data.Tin, data.p)
	cf = heatcap_gas(Air_, data.Tin)
	λf = thermcond_gas(Air_, data.Tin)
	# Peclet number
	Pe0 = u0*ρf*cf*data.d/λf
	kbed(data)
end

# ╔═╡ 6d5a7d83-53f9-43f3-9ccd-dadab08f62c1
md"""
Atmospheric pressure operation: __p = $(data.p/ufac"bar") bar__

Reactor to be fed with 1/1 mixture of CO₂/H₂

Max. volumetric feed flow rate for each: Q = $(0.5*data.Qflow/ufac"ml/minute") ml/min

Total feed volumetric flow rate: __$(data.Qflow/ufac"ml/minute") ml/min__

For frit diameter of __$(data.D/ufac"cm") cm__, porosity of __$(data.ϕ)__ the mean superficial velocity is __$(round(u0/ufac"cm/s",sigdigits=2)) cm/s__.
"""

# ╔═╡ 0727ad12-2c8f-47ba-9599-2843d2dac016
# all values taken from VDI heat atlas 2010 chapter D3.1
N2_=(
	# Table 1
	MW=28.01*ufac"g/mol",
	# Table 6
	HeatCap=PropsCoeffs(
	A=432.2027,
	B=3.5160,
	C=2.8021,
	D=-4.1924,
	E=42.0153,
	F=-114.2500,
	G=111.1019
	),
	# Table 10
	ThermCond = PropsCoeffs(
	A=-0.133e-3,
	B=0.101e-3,
	C=-0.060650e-6,
	D=0.033610e-9,
	E=-0.0071e-12
	),
	# Table 8
	DynVisc= PropsCoeffs(
	A=-0.01020e-5,
	B=0.74785e-7,
	C=-0.59037e-10,
	D=0.03230e-12,
	E=-0.00673e-15
)
);

# ╔═╡ Cell order:
# ╠═7d8eb6f5-3ba6-46ef-8058-1f24a0938ed1
# ╠═5c3adaa0-9285-11ed-3ef8-1b57dd870d6f
# ╟─f353e09a-4a61-4def-ab8a-1bd6ce4ed58f
# ╟─2015c8e8-36cd-478b-88fb-94605283ac29
# ╠═98063329-31e1-4d87-ba85-70419beb07e9
# ╟─03d0c88a-462b-43c4-a589-616a8870be64
# ╟─6d5a7d83-53f9-43f3-9ccd-dadab08f62c1
# ╠═3b3595c4-f53d-4827-918e-edcb74dd81f8
# ╠═d912a1ca-1b69-4ea1-baa5-69794e004693
# ╟─b4cee9ac-4d14-4169-90b5-25d12ac9c003
# ╟─bbd0b076-bcc1-43a1-91cb-d72bb17d3c88
# ╠═7c4d4083-33e2-4b17-8575-214ef458d75d
# ╟─d7317b2d-e2c7-4114-8985-51979f2205ba
# ╠═2fe11550-683d-4c4b-b940-3e63a4f8a87d
# ╠═8cd85a0e-3d11-4bcc-8a7d-f30313b31363
# ╟─9d8c6ddc-2662-4055-b636-649565c36287
# ╠═f15fd785-010c-4fda-ab4f-7947642556dd
# ╠═d725f9b9-61c4-4724-a1d9-6a04ba42499d
# ╟─e148212c-bc7c-4553-8ffa-26e6c20c5f47
# ╠═6a92c1c7-fb51-4363-b8d0-18eeb24087a8
# ╠═9be8cb32-88ae-43e2-b7b6-6ff8428c8f3c
# ╠═534f29f2-2563-465b-9a37-bf1791ca22d9
# ╠═bd0fd520-217c-471e-89af-fd568ea70be5
# ╠═d4bd9a44-f96b-4c94-af47-fdb2c09eda82
# ╠═4f9e10ae-8c37-4c89-9658-fb47095a1c1e
# ╠═29db6427-f1a0-4a2e-8a3b-82d2ad48b42b
# ╠═453db496-4688-48ce-8e2a-1a960540f69f
# ╠═b4f7d4d5-c277-4a1d-9dff-0a35853d80ef
# ╠═0727ad12-2c8f-47ba-9599-2843d2dac016
# ╠═aaee8aa5-a41c-4879-bea5-485f1188f9d0
# ╠═37cf0697-0091-49f8-bda0-9864314b2dff
