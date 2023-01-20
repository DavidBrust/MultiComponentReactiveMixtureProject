### A Pluto.jl notebook ###
# v0.19.19

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

# ╔═╡ 5c3adaa0-9285-11ed-3ef8-1b57dd870d6f
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	
	using VoronoiFVM
	using ExtendableGrids, SimplexGridFactory, TetGen
	using GridVisualize
	using LessUnitful
	using PlutoVista
	using PlutoUI
	using PyPlot

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
	#GridVisualize.default_plotter!(PyPlot)
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

# ╔═╡ fd52fbb0-afca-4361-8604-2cb9185af275
abstract type AbstractModelData end

# ╔═╡ 98063329-31e1-4d87-ba85-70419beb07e9
Base.@kwdef mutable struct ModelData <:AbstractModelData
	#iT::Int64=1 # index of Temperature variable
	iTs::Int64=1
	iTf::Int64=2
	
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	α_nc::Float64=15.0*ufac"W/(m^2*K)" # natural convection heat transfer coefficient

	## irradiation data
	G_lamp::Float64=50.0*ufac"kW/m^2" # solar simulator irradiation flux
	Abs_lamp::Float64=0.7 # avg absorptivity of cat. of irradiation coming from lamp
	Eps_ir::Float64=0.7 # avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
	
	
	## porous filter data
	d::Float64=100.0*ufac"μm" # average pore size
	# cylindrical disc / 2D
    D::Float64=10.0*ufac"cm" # disc diameter
	Ac::Float64=pi*D^2.0/4.0*ufac"m^2" # cross-sectional area

	# prism / 3D
	wi::Float64=10.0*ufac"cm" # prism width/side lenght
	le::Float64=wi # prism width/side lenght
	h::Float64=0.5*ufac"cm" # frit thickness (applies to 2D & 3D)
	
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2 	
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
	
	ϕ::Float64=0.36 # porosity, class 2
	k::Float64=2.9e-11*ufac"m^2" # permeability
	a_s::Float64=0.13*ufac"m^2/g" # specific surface area
	ρfrit::Float64=(1.0-ϕ)*ρs+ϕ*density_idealgas(Air, 298.15, 1.0*ufac"atm")*ufac"kg/m^3" # density of porous frit
	a_v::Float64=a_s*ρfrit # volume specific interface area
	#a_v::Float64=a_s*ρfrit*1.0e-3 # volume specific interface area
	## END porous filter data

	## fluid data
	Qflow::Float64=2000.0*ufac"ml/minute" # volumetric feed flow rate
	Tin::Float64=298.15*ufac"K" # inlet temperature
	p::Float64=1.0*ufac"atm" # reactor pressure		
	u0::Float64=Qflow/(Ac*ϕ)*ufac"m/s" # mean superficial velocity
	# fluid properties: Air
	# values taken from VDI heat atlas 2010 chapter D3.1
	Fluid::FluidProps=Air
	## END fluid data
	
end;

# ╔═╡ 03d0c88a-462b-43c4-a589-616a8870be64
md"""
# Experimental Conditions
"""

# ╔═╡ 3b3595c4-f53d-4827-918e-edcb74dd81f8
data = ModelData(;p=1.0*ufac"atm",Qflow=3400*ufac"ml/minute")

# ╔═╡ 6d5a7d83-53f9-43f3-9ccd-dadab08f62c1
md"""
Atmospheric pressure operation: __p = $(data.p/ufac"bar") bar__

Reactor to be fed with 1/1 mixture of CO₂/H₂

Max. volumetric feed flow rate for each: Q = $(0.5*data.Qflow/ufac"ml/minute") ml/min

Total feed volumetric flow rate: __$(data.Qflow/ufac"ml/minute") ml/min__

For frit diameter of __$(data.D/ufac"cm") cm__, porosity of __$(data.ϕ)__ the mean superficial velocity is __$(round(data.u0/ufac"cm/s",sigdigits=2)) cm/s__.
"""

# ╔═╡ 4bcdb950-ed22-496c-ad70-e0c0fa4d7f52
md"""
## Dimensionless numbers
"""

# ╔═╡ 73b71898-0268-42bd-b2e6-d0c5118700dd
function RePrPe(data::AbstractModelData,T,p)
	d=data.d
	u0=data.u0
	ρf = density_idealgas(data.Fluid, T, p)
	ηf = dynvisc_gas(data.Fluid, T)
	cf = heatcap_gas(data.Fluid, T)
	λf = thermcond_gas(data.Fluid, T)
	
	Re = u0*ρf*d/ηf # Reynolds number
	Pr = cf*ηf/λf # Prandtl number
	Pe = u0*ρf*cf*d/λf # Peclet
	Re,Pr,Pe
end

# ╔═╡ 7e83918e-3ba4-4bbb-be8c-839eb32def13
Re,Pr,Pe = RePrPe(data,data.Tin,data.p)

# ╔═╡ 13e66a6a-b329-40e8-9098-05f4077d1242
md"""
At given experimental conditions the Reynolds, Prandtl and Peclet numbers assume the following values:
- Re = $(round(Re,sigdigits=2))
- Pr = $(round(Pr,sigdigits=2))
- Pe = $(round(Pe,sigdigits=2))
"""

# ╔═╡ cb6a357f-e244-4725-a04a-3e006dd4b53d
md"""
## Irradiation Boundary Condition
"""

# ╔═╡ 463a9a2b-8437-407f-b31a-dde3165f49ad
md"""
### Irradiation conditions

Catalyst surface temperature critically depends on the __irradition__ concentration and the __optical properties__ of the catalyst and components involved.

- Solar simulator (lamp) mean irradiance: $(data.G_lamp) W m⁻²
- absorptivity of catalyst material of irradiation coming from lamp: $(data.Abs_lamp)
- absorptivity/emissivity of catalyst material of IR irradiation from surroundings/ emission: $(data.Eps_ir)
- Temperature of surroundings: $(data.Tamb) K

There can be __different values for optical properties__ of the catalyst depending on __wavelength__: e.g. a different absorption coefficients for radiation coming from the solar simulator (lamp) and from the surroundings.

This needs to be investigated further.
"""

# ╔═╡ 387b5b8e-a466-4a65-a360-fa2cf08092e3
md"""
$(LocalResource("../img/IrradBC.png", :width => 800))
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

# ╔═╡ 14d0d6e9-adde-4e51-a3c5-18a326e6faf8
md"""
# Heterogeneous phase model
"""

# ╔═╡ 9198330d-7c05-440a-916e-fca0fe796b7b
md"""
Treat the fluid and the solid phase (porous material) separately: one energy balance for each phase.
"""

# ╔═╡ 16f5e0bc-8e3d-40cd-b67b-694eda6b67d9
md"""
## Interfacial heat transfer coefficient
"""

# ╔═╡ 1459c3db-5ffc-46bd-9c94-8c8964519f39
md"""
When working with a heterogeneous phase model (separate energy balances for both fluid and porous solid material), the exchange of energe between the phases can be described by an interfacial heat transfer coefficient. It can be calculated according to:

__Kuwahara, F., Shirota, M., & Nakayama, A. (2001).__ A numerical study of interfacial convective heat transfer coefficient in two-energy equation model for convection in porous media. International Journal of Heat and Mass Transfer, 44(6), 1153-1159. doi:10.1016/s0017-9310(00)00166-6

```math
\frac{h_{\text{sf}} \text D}{k_{\text f}}= \left( 1+ \frac{4(1- \phi)}{\phi} \right) + \frac{1}{2} (1-\phi)^{\frac{1}{2}} \text{Re}^{0.6}_D\text{Pr}^{\frac{1}{3}}
```

"""

# ╔═╡ 4eb00d7b-d10e-478d-a1df-9eea3362ef5f
function hsf(data::AbstractModelData,T,p)
	Re,Pr,_ = RePrPe(data,T,p)
	λf = thermcond_gas(data.Fluid, T)
	ϕ = data.ϕ
	d = data.d
	λf/d*((1.0 + 4*(1.0-ϕ)/ϕ) + 0.5*(1.0-ϕ)^0.5*Re^0.6*Pr^(1.0/3.0))*ufac"W/(m^2*K)"
end

# ╔═╡ 641d61ca-448b-4240-95fb-486e83b6b768
md"""
For the given porous medium with very fine pore and particle sizes (__$(round(data.d/ufac"μm",sigdigits=2)) μm__), the volume specific interfacial area ``A_{\text v} = `` $(round(data.a_v,sigdigits=4)) `` \text{m}^2`` and interfacial heat transfer coefficient ``h_{\text{sf}} = `` $(round(hsf(data,data.Tin,data.p),sigdigits=4)) ``\text W/ \text m^2 \text K`` take on very large values. Therefore the solid and gas phases are in thermal equilibrium. This justifies the use of a quasi-homogeneous model, describing both phases by a single temperature.
"""

# ╔═╡ d7317b2d-e2c7-4114-8985-51979f2205ba
md"""
# Grid 
"""

# ╔═╡ 3c75c762-a44c-4328-ae41-a5016ce181f1
md"""
## 2D
"""

# ╔═╡ 2c31e63a-cf42-45cd-b367-112438a02a97
md"""
Assume axysymmetric geometry (thin cylindrical disk) to allow a 2-dimensional (rotationally symmetric) representation of the domain. 
"""

# ╔═╡ 2fe11550-683d-4c4b-b940-3e63a4f8a87d
function cylinder(;nref=0, r=5.0*ufac"cm", h=0.5*ufac"cm")
    #step=0.1*ufac"cm"*2.0^(-nref)
	hr=r/10.0*2.0^(-nref)
	hh=h/10.0*2.0^(-nref)
    R=collect(0:hr:r)
    Z=collect(0:hh:h)
    grid=simplexgrid(R,Z)
    circular_symmetric!(grid)
	grid
end

# ╔═╡ 8cd85a0e-3d11-4bcc-8a7d-f30313b31363
gridplot(cylinder(;r=data.D/2,h=data.h))

# ╔═╡ a190862c-2251-4110-8274-9960c495a2c4
md"""
## 3D
"""

# ╔═╡ 4d9145f3-06aa-4a7c-82d0-feee0fa01865
function prism_sq(;nref=0, l=10.0*ufac"cm", w=10.0*ufac"cm", h=0.5*ufac"cm")
	
	hw=w/2.0/10.0*2.0^(-nref)
	hl=l/2.0/10.0*2.0^(-nref)
	hh=h/10.0*2.0^(-nref)
	W=collect(0:hw:(w/2.0))
    L=collect(0:hl:(l/2.0))
    H=collect(0:hh:h)
	
	simplexgrid(W,L,H)	
end

# ╔═╡ 61a67079-cb15-4283-ac15-96b49c461b6e
gridplot(prism_sq(nref=0))

# ╔═╡ ed50c2d4-25e9-4159-84c7-e0c70ffa63a1
md"""
Using symmetry, it is sufficient to include 1/4 of the square prismatic domain. (in principle 1/8 would be sufficient, but that would preclude the application of a simple grid)

The 3D geomtry has 6 outer facets, whose boundary conditions need to be specified:
- facet 1: symmetry (no flux)
- facet 4: symmetry (no flux)
- facet 2: convective heat transfer to the wall, air gab between porous frit and Al reactor wall (robin bc.)
- facet 3: convective heat transfer to the wall, air gab between porous frit and Al reactor wall (robin bc.)
- facet 5: convective heat transfer to inflowing gas stream (robin bc.)
- facet 6: radiation + convection (outflow of gas stream) 
"""

# ╔═╡ 9d8c6ddc-2662-4055-b636-649565c36287
md"""
# Simulation
"""

# ╔═╡ e9003129-31db-4f9d-b289-638511c7ec26
# ╠═╡ disabled = true
#=╠═╡
Sim2D=main(nref=1);
  ╠═╡ =#

# ╔═╡ ba5c2095-4858-444a-99b5-ae6cf40374f9
md"""
## 2D
"""

# ╔═╡ f15fd785-010c-4fda-ab4f-7947642556dd
#=╠═╡
let
	sys,sol,data=Sim2D
	iT=data.iT
	vis=GridVisualizer()
	solC = copy(sol)
	@. solC[iT,:] -= 273.15

	scalarplot!(vis,sys,solC;species=iT,title="Temperature / °C",xlabel="Radial coordinate / m", ylabel="Axial coordinate / m",legend=:best,colormap=:summer,show=true)

end
  ╠═╡ =#

# ╔═╡ 667c095d-f7c7-4244-806f-a70f1250146e
md"""
Fluid phase:
```math
	-\nabla \cdot \left ( k^{\text f} \nabla T_{\text f} - \phi \rho_{\text f} c_{p, \text f} \vec{u} T_{\text f} \right) =  h_{\text{sf}}A_{\text V} (T_{\text s}-T_{\text f}) 
```
"""

# ╔═╡ f0054c8a-921c-4603-9567-fb98beab4b69
md"""
Solid phase:
```math
	-\nabla \cdot \left ( k_{\text{eff}}^{\text s} \nabla T_{\text s} \right) = - h_{\text{sf}}A_{\text V} (T_{\text s}-T_{\text f}) + \dot{q}_{\text{chem}} + \dot{q}_{\text{rad}}
```
"""

# ╔═╡ d725f9b9-61c4-4724-a1d9-6a04ba42499d
function main2phase(;nref=0,p=1.0*ufac"atm",Qflow=3400*ufac"ml/minute")
	data=ModelData(Qflow=Qflow,	p=p,)
	
	iTs=data.iTs
	iTf=data.iTf

	# function return 2D velocity vector: flow upward in z-direction
    function fup(r,z)
        return 0,-data.u0
    end    
	
	function flux(f,u,edge,data)
		(;Fluid,p,ϕ)=data
		# Fluid phase
		Tbar=0.5*(u[iTf,1]+u[iTf,2])
		ρf=density_idealgas(Fluid, Tbar, p)
		cf=heatcap_gas(Fluid, Tbar)
		λf=thermcond_gas(Fluid, Tbar)
		λbed=kbed(data)*λf
		#conv=ϕ*evelo[edge.index]*ρf*cf/λf
		#Bp,Bm = fbernoulli_pm(conv)
		#f[iTf]= λf*(Bm*u[iTf,1]-Bp*u[iTf,2])
		conv=evelo[edge.index]*ρf*cf/λbed
		Bp,Bm = fbernoulli_pm(conv)
		f[iTf]= λbed*(Bm*u[iTf,1]-Bp*u[iTf,2])
		
		# Solid phase		
		λf0=thermcond_gas(Fluid, data.Tin)
		#λbed=kbed(data)*λf0
		
		#f[iTs]= λbed*(u[iTs,1]-u[iTs,2])
		f[iTs]= λbed*(Bm*u[iTs,1]-Bp*u[iTs,2])

		
	end

	function reaction(f,u,edge,data)
		(;Fluid,p,ϕ,a_v)=data
		hsf_=hsf(data,u[iTf],data.p)
		ip_htx = a_v*hsf_*(u[iTs]-u[iTf]) # heat exchange between solid and gas phase
		f[iTs] = ip_htx
		f[iTf] = -ip_htx
		
	end

	function irrad_bc(f,u,bnode,data)
		if bnode.region==3 # top boundary
			flux_rerad = data.Eps_ir*ph"σ"*(u[iTs]^4 - data.Tamb^4)
			flux_convec = data.α_nc*(u[iTs]-data.Tamb)
			f[iTs] = -(data.Abs_lamp*data.G_lamp - flux_rerad - flux_convec)
			f[iTf] = -(data.Abs_lamp*data.G_lamp - flux_rerad - flux_convec)
		end
	end

	function bcondition(f,u,bnode,data)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=1,value=data.Tamb)
		boundary_robin!(f,u,bnode;species=iTs,region=1, factor=data.α_nc, value=data.Tamb*data.α_nc)
		boundary_robin!(f,u,bnode;species=iTf,region=1, factor=data.α_nc, value=data.Tamb*data.α_nc)
		
		#boundary_dirichlet!(f,u,bnode;species=iTf,region=1, value=data.Tamb)
		
		boundary_robin!(f,u,bnode;species=iTs,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
		boundary_robin!(f,u,bnode;species=iTf,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
		
		#boundary_dirichlet!(f,u,bnode;species=iT,region=3,value=data.Tamb+300.0)
		# irradiation boundary condition
		irrad_bc(f,u,bnode,data)
	end
	

	
	grid=cylinder(;nref=nref,r=data.D/2,h=data.h)
	evelo=edgevelocities(grid,fup)
	
	sys=VoronoiFVM.System(grid;
                          data=data,
                          flux=flux,
                          reaction=reaction,
    #                      #storage=pnpstorage,
                          bcondition,
                          species=[iTs,iTf],
	#					  regions=[1,2],
    #                      kwargs...
                          )
	inival=unknowns(sys)
	inival[iTs,:] .= map( (r,z)->(data.Tamb+500*z/data.h),grid)
	inival[iTf,:] .= map( (r,z)->(data.Tamb+500*z/data.h),grid)
	#inival[iT,:] .= data.Tamb
	sol=solve(inival,sys)
	sys,sol,data
end

# ╔═╡ 186c0b6f-a049-4841-a69c-34b982c3d17c
Sim2Dhet=main2phase(nref=1);

# ╔═╡ d0ed3983-5118-479f-855a-1cd3c4778771
let
	sys,sol,data=Sim2Dhet
	iTs=data.iTs
	iTf=data.iTf
	vis=GridVisualizer(layout=(2,1))
	solC = copy(sol)
	@. solC[iTs,:] -= 273.15
	@. solC[iTf,:] -= 273.15

	scalarplot!(vis[1,1],sys,solC;species=iTs,title="Solid Temp / °C",xlabel="Radial coordinate / m", ylabel="Axial coordinate / m",legend=:best,colormap=:summer,show=true)

	scalarplot!(vis[2,1],sys,solC;species=iTf,title="Fluid Temp / °C",xlabel="Radial coordinate / m", ylabel="Axial coordinate / m",legend=:best,colormap=:summer,show=true)

end

# ╔═╡ 28a2230f-5a59-4034-86af-e3d58dcceb6c
md"""
## 3D
"""

# ╔═╡ f435b4df-9162-42bd-8154-5f3434ea0e2a
# ╠═╡ disabled = true
#=╠═╡
function main3D(;nref=0,p=1.0*ufac"atm",Qflow=3400*ufac"ml/minute")
	data=ModelData(Qflow=Qflow,	p=p,)
	iT=data.iT

	# function return 3D velocity vector: flow upward in z-direction
    function fup(x,y,z)
        return 0,0,-data.u0
    end    
	
	function flux(f,u,edge,data)
		(;Fluid,u0,p)=data
		Tbar=0.5*(u[iT,1]+u[iT,2])
		ρf=density_idealgas(Fluid, Tbar, p)
		cf=heatcap_gas(Fluid, Tbar)
		λf=thermcond_gas(Fluid, Tbar)
		#λf=thermcond_gas(Fluid, Tin)
		λbed=kbed(data)*λf
		
		vh=project(edge,(0,0,data.u0))
		conv=vh*ρf*cf/λbed
		#conv=evelo[edge.index]*ρf*cf/λbed
		Bp,Bm = fbernoulli_pm(conv)
		f[iT]= λbed*(Bm*u[iT,1]-Bp*u[iT,2])
		
		#f[iT] = λbed*(u[iT,1]-u[iT,2])
		
	end

	function irrad_bc(f,u,bnode,data)
		if bnode.region==6 # top boundary
			flux_rerad = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
			flux_convec = data.α_nc*(u[iT]-data.Tamb)
			f[iT] = -(data.Abs_lamp*data.G_lamp - flux_rerad - flux_convec)
		end
	end

	function bcondition(f,u,bnode,data)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=1,value=data.Tamb)
		boundary_robin!(f,u,bnode;species=iT,region=5, factor=data.α_nc, value=data.Tamb*data.α_nc) # bottom
		boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
		# for prism of square
		boundary_robin!(f,u,bnode;species=iT,region=3, factor=data.α_w, value=data.Tamb*data.α_w)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=2,value=data.Tamb+300.0)
		# irradiation boundary condition
		irrad_bc(f,u,bnode,data)
	end
	

	
	#grid=prism(;nref=nref, w=data.wi, l=data.le, h=data.h)
	grid=prism_sq(;nref=nref, w=data.wi, l=data.le, h=data.h)
	
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
	#inival[iT,:] .= map( (x,y,z)->(data.Tamb+300*z/data.h),grid)
	inival[iT,:] .= data.Tamb
	sol=solve(inival,sys)
	sys,sol,data,nref
end
  ╠═╡ =#

# ╔═╡ ea81db48-d02c-4438-a0ba-9c54a3ea0e52
#=╠═╡
Sim3D=main3D(nref=1);
  ╠═╡ =#

# ╔═╡ fb72aede-8997-48ef-9a2d-98acbf372747
#=╠═╡
md"""
f=$(@bind flevel Slider(range(([minimum(Sim3D[2]),maximum(Sim3D[2])].-273.15)... , length=20),default=(minimum(Sim3D[2])+maximum(Sim3D[2]))/2-273.15,show_value=true))

x=$(@bind xplane Slider(range(0.0,data.wi/2,length=20),default=0.0,show_value=true))

y=$(@bind yplane Slider(range(0.0,data.le/2,length=20),default=0.0,show_value=true))

z=$(@bind zplane Slider(range(0.0,data.h,length=20),default=0.0,show_value=true))

"""
  ╠═╡ =#

# ╔═╡ 8107492e-fdaf-4a86-96ef-ea34b45c0e67
#=╠═╡
let
	sys,sol,data,nref=Sim3D
	iT=data.iT
	vis=GridVisualizer(resolution=(800,800),)
	solC=copy(sol)
	@. solC[iT,:] -= 273.15

	scalarplot!(vis,sys,solC;species=iT,xlabel="X / m", ylabel="Y / m", zlabel="Z / m", levels=[flevel], xplanes=[xplane], yplanes=[yplane], zplanes=[zplane], levelalpha=1.0,title="Temperature / °C", outlinealpha=0.05)
	scene=reveal(vis)
	#save("plot_3D.svg",scene)
end
  ╠═╡ =#

# ╔═╡ 64dd5097-16aa-4c44-b000-6177cd4be226
md"""
### Cut planes
"""

# ╔═╡ bd179765-f996-4f91-ac4e-57d5817a2ed6
md"""
Analyse 2D slices of the 3D solution to be able to compare with the 2D axisymmetric solution. The goal is to show, that the geometry can be treated as 2D with sufficient accuracy.
Below the difference between 3D and 2D __(3D - 2D)__ calculation is shown for different cross-sections:
"""

# ╔═╡ 641988da-8888-4e0d-b720-4f21a9900aca
md"""
Y - cutplane $(@bind ycut Slider(range(0.0,data.wi/2,length=21),default=0.0,show_value=true))
"""

# ╔═╡ b4175b62-198d-4605-9cf0-04b0be52c9c0
function plane(ypos,sol,data,nref)
	grid=prism_sq(;nref=nref, w=data.wi, l=data.le, h=data.h)
	
	bfacemask!(grid, [0,ypos,0],[data.wi/2.0,ypos,data.h],7)

	# transform z coordinate of parent grid into y coordinate of subgrid
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[3]
	end
	grid_2D  = subgrid(grid, [7], boundary=true, transform=_3to2) 
	
	sol_cutplane = view(sol[data.iT, :], grid_2D)
		
	collect(sol_cutplane), grid_2D	
	#sol_cutplane, grid_2D	
end

# ╔═╡ 3612d83d-a8bc-4c5f-8ea0-a4c975929500
#=╠═╡
let
	sys,sol,data=Sim2D


	sys3,sol3,data,nref=Sim3D
	sol_cutplane, grid_2D = plane(ycut,sol3,data,nref)

	sol_reorder=zeros(length(sol))

	k=1
	for coord in eachcol(grid_2D[Coordinates])
		id = findfirst(x->x==coord, eachcol(sys.grid[Coordinates]))
		sol_reorder[k] = sol[id]
		k +=1
	end

	vis=GridVisualizer()
	#scalarplot!(vis,grid_2D,sol_reorder;colormap=:summer,show=true)
	scalarplot!(vis,grid_2D,sol_cutplane.-sol_reorder;colormap=:bwr,show=true)

end
  ╠═╡ =#

# ╔═╡ e148212c-bc7c-4553-8ffa-26e6c20c5f47
md"""
# Thermo-physical properties
"""

# ╔═╡ 6a92c1c7-fb51-4363-b8d0-18eeb24087a8
let
	(;Fluid)=data
	Ts=(273.15:1:573.15)
	λf = zeros(Float64, length(Ts))
	for (i,T) in enumerate(Ts)
		λf[i]=thermcond_gas(Fluid, T)
	end
	PlutoVista.plot(Ts.-273.15, λf,xlabel="Temperature / °C",ylabel="Thermal Conductivity / W m⁻¹ K⁻¹", title=Fluid.name )
end

# ╔═╡ 8310812b-97af-45d9-aff8-234ffbb169af
md"""
# Analysis
"""

# ╔═╡ 909044b9-d40d-4129-bfef-6cd131a932f2
md"""
Calculate heat flow ``(\text W)`` into the domain from surface integral over boundaries.
"""

# ╔═╡ 30f0e723-a5a6-4b49-80d2-a29efb0f417b
data.Abs_lamp*data.G_lamp*data.Ac

# ╔═╡ a1f5c0c2-290e-418f-b381-4626d61c2bfd
let
	sys,sol,data=main(nref=1)
	tf=VoronoiFVM.TestFunctionFactory(sys)
	
	# testfunction arguments:
		# 2nd arg: boundaries, where the test func. value is 0 -> these boundaries are excludes from the surface integral
		# 3rd arg: boundaries, where the test func. value is 1 -> these boundaries are included in the surface integral
	
	# flux integral over bottom boundary: diffusive in/out flux
	tf_irr=testfunction(tf,[1,2,4],3)
	integrate(sys,tf_irr,sol)	
end

# ╔═╡ Cell order:
# ╠═7d8eb6f5-3ba6-46ef-8058-1f24a0938ed1
# ╠═5c3adaa0-9285-11ed-3ef8-1b57dd870d6f
# ╟─f353e09a-4a61-4def-ab8a-1bd6ce4ed58f
# ╟─2015c8e8-36cd-478b-88fb-94605283ac29
# ╠═fd52fbb0-afca-4361-8604-2cb9185af275
# ╠═98063329-31e1-4d87-ba85-70419beb07e9
# ╟─03d0c88a-462b-43c4-a589-616a8870be64
# ╟─6d5a7d83-53f9-43f3-9ccd-dadab08f62c1
# ╠═3b3595c4-f53d-4827-918e-edcb74dd81f8
# ╟─4bcdb950-ed22-496c-ad70-e0c0fa4d7f52
# ╠═73b71898-0268-42bd-b2e6-d0c5118700dd
# ╠═7e83918e-3ba4-4bbb-be8c-839eb32def13
# ╟─13e66a6a-b329-40e8-9098-05f4077d1242
# ╟─cb6a357f-e244-4725-a04a-3e006dd4b53d
# ╟─463a9a2b-8437-407f-b31a-dde3165f49ad
# ╟─387b5b8e-a466-4a65-a360-fa2cf08092e3
# ╟─b4cee9ac-4d14-4169-90b5-25d12ac9c003
# ╟─bbd0b076-bcc1-43a1-91cb-d72bb17d3c88
# ╠═7c4d4083-33e2-4b17-8575-214ef458d75d
# ╟─14d0d6e9-adde-4e51-a3c5-18a326e6faf8
# ╟─9198330d-7c05-440a-916e-fca0fe796b7b
# ╟─16f5e0bc-8e3d-40cd-b67b-694eda6b67d9
# ╟─1459c3db-5ffc-46bd-9c94-8c8964519f39
# ╠═4eb00d7b-d10e-478d-a1df-9eea3362ef5f
# ╟─641d61ca-448b-4240-95fb-486e83b6b768
# ╟─d7317b2d-e2c7-4114-8985-51979f2205ba
# ╟─3c75c762-a44c-4328-ae41-a5016ce181f1
# ╟─2c31e63a-cf42-45cd-b367-112438a02a97
# ╠═2fe11550-683d-4c4b-b940-3e63a4f8a87d
# ╠═8cd85a0e-3d11-4bcc-8a7d-f30313b31363
# ╟─a190862c-2251-4110-8274-9960c495a2c4
# ╠═4d9145f3-06aa-4a7c-82d0-feee0fa01865
# ╠═61a67079-cb15-4283-ac15-96b49c461b6e
# ╟─ed50c2d4-25e9-4159-84c7-e0c70ffa63a1
# ╟─9d8c6ddc-2662-4055-b636-649565c36287
# ╠═e9003129-31db-4f9d-b289-638511c7ec26
# ╟─ba5c2095-4858-444a-99b5-ae6cf40374f9
# ╠═f15fd785-010c-4fda-ab4f-7947642556dd
# ╟─667c095d-f7c7-4244-806f-a70f1250146e
# ╟─f0054c8a-921c-4603-9567-fb98beab4b69
# ╠═186c0b6f-a049-4841-a69c-34b982c3d17c
# ╠═d0ed3983-5118-479f-855a-1cd3c4778771
# ╠═d725f9b9-61c4-4724-a1d9-6a04ba42499d
# ╟─28a2230f-5a59-4034-86af-e3d58dcceb6c
# ╠═ea81db48-d02c-4438-a0ba-9c54a3ea0e52
# ╠═8107492e-fdaf-4a86-96ef-ea34b45c0e67
# ╟─fb72aede-8997-48ef-9a2d-98acbf372747
# ╠═f435b4df-9162-42bd-8154-5f3434ea0e2a
# ╟─64dd5097-16aa-4c44-b000-6177cd4be226
# ╟─bd179765-f996-4f91-ac4e-57d5817a2ed6
# ╟─3612d83d-a8bc-4c5f-8ea0-a4c975929500
# ╟─641988da-8888-4e0d-b720-4f21a9900aca
# ╠═b4175b62-198d-4605-9cf0-04b0be52c9c0
# ╟─e148212c-bc7c-4553-8ffa-26e6c20c5f47
# ╠═6a92c1c7-fb51-4363-b8d0-18eeb24087a8
# ╟─8310812b-97af-45d9-aff8-234ffbb169af
# ╟─909044b9-d40d-4129-bfef-6cd131a932f2
# ╠═30f0e723-a5a6-4b49-80d2-a29efb0f417b
# ╠═a1f5c0c2-290e-418f-b381-4626d61c2bfd
