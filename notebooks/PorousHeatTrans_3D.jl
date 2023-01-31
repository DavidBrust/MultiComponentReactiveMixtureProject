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

# ╔═╡ c258e0e6-72cc-4b3e-8f6d-e629e4a0d7fd
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 1600px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

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

# ╔═╡ 5136828c-6e88-40fc-b9a5-3c86db2f0879
abstract type AbstractModelData end

# ╔═╡ 98063329-31e1-4d87-ba85-70419beb07e9
Base.@kwdef mutable struct ModelData <:AbstractModelData
	iT::Int64=1 # index of Temperature variable
	
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	#α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	α_w::Float64=0.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	α_nc::Float64=15.0*ufac"W/(m^2*K)" # natural convection heat transfer coefficient

	## irradiation data
	G_lamp::Float64=1.0*ufac"kW/m^2" # solar simulator irradiation flux
	Abs_lamp::Float64=0.7 # avg absorptivity of cat. of irradiation coming from lamp
	Eps_ir::Float64=0.7 # avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
	#Eps_ir::Float64=0.0 # avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
	
	
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
	#Qflow::Float64=100000.0*ufac"ml/minute" # volumetric feed flow rate
	#Qflow::Float64=3400.0*ufac"ml/minute" # volumetric feed flow rate
	Qflow::Float64=0.0*ufac"ml/minute" # volumetric feed flow rate
	Tin::Float64=298.15*ufac"K" # inlet temperature
	p::Float64=1.0*ufac"atm" # reactor pressure		
	# u0::Float64=Qflow/(Ac*ϕ)*ufac"m/s" # mean superficial velocity
	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	# fluid properties: Air
	# values taken from VDI heat atlas 2010 chapter D3.1
	Fluid::FluidProps=Air
	## END fluid data
	
end;

# ╔═╡ 03d0c88a-462b-43c4-a589-616a8870be64
md"""
# Experimental Conditions
"""

# ╔═╡ 821aecc3-cf5d-4db2-bace-dab82c8cd711
md"""
## Average Superficial Velocity
"""

# ╔═╡ 561a9857-d29e-40cd-809f-cc667613251c
md"""
Following the definition provided in VDI heat atlas, the average superficiel veloctiy ``u_0`` is defined based on the volumetric flow rate and the cross-sectional area of the flow path according to:
```math
\begin{align}
	u_0 &= \frac{\dot V}{A_{\text c}} \\
	A_{\text c} &= \frac{\pi D^2}{4}, \quad \text{circular cross section} \\
	A_{\text c} &= W^2, \quad \text{square cross section}
\end{align}
```
"""

# ╔═╡ 3b3595c4-f53d-4827-918e-edcb74dd81f8
data = ModelData()

# ╔═╡ 6d5a7d83-53f9-43f3-9ccd-dadab08f62c1
md"""
Atmospheric pressure operation: __p = $(data.p/ufac"bar") bar__

Reactor to be fed with 1/1 mixture of CO₂/H₂

Max. volumetric feed flow rate for each: Q = $(0.5*data.Qflow/ufac"ml/minute") ml/min

Total feed volumetric flow rate: __$(data.Qflow/ufac"ml/minute") ml/min__

"""

# ╔═╡ fdddb318-831c-4277-a5a3-b4be8f63c430
md"""
For a (sqare) frit with side length of __$(data.wi/ufac"cm") cm__ the mean superficial velocity is __$(round(data.u0/ufac"cm/s",sigdigits=2)) cm/s__.
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

# ╔═╡ c9f326c1-abc0-412d-8348-8e23a7027661
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
gridplot(cylinder(;r=data.D/2,h=data.h),resolution=(1000,1000))

# ╔═╡ a190862c-2251-4110-8274-9960c495a2c4
md"""
## 3D
"""

# ╔═╡ 4d9145f3-06aa-4a7c-82d0-feee0fa01865
function prism_sq(data;nref=0, w=data.wi, h=data.h)
	
	hw=w/2.0/10.0*2.0^(-nref)
	#hl=l/2.0/10.0*2.0^(-nref)
	hh=h/10.0*2.0^(-nref)
	W=collect(0:hw:(w/2.0))
    #L=collect(0:hl:(l/2.0))
    H=collect(0:hh:h)
	
	simplexgrid(W,W,H)	
end

# ╔═╡ 61a67079-cb15-4283-ac15-96b49c461b6e
gridplot(prism_sq(data,nref=0),Plotter=PlutoVista,resolution=(1000,1000))

# ╔═╡ ebe2ed24-e2d6-4652-ae69-1a59747b0c4c
function prism(;nref=0, l=10.0*ufac"cm", w=10.0*ufac"cm", h=0.5*ufac"cm")
	builder=SimplexGridBuilder(Generator=TetGen)
	W=w/2
	L=l/2
	H=h
	p1=point!(builder,0,0,0)
    p2=point!(builder,W,0,0)
    p3=point!(builder,W,L,0)
    p4=point!(builder,0,0,H)
    p5=point!(builder,W,0,H)
    p6=point!(builder,W,L,H)

	facetregion!(builder,1)
    facet!(builder,p1 ,p2 ,p3)
    facetregion!(builder,2)
    facet!(builder,p4 ,p5 ,p6)
    facetregion!(builder,3)
    facet!(builder,p1 ,p2 ,p5 ,p4)
    facetregion!(builder,4)
    facet!(builder,p2 ,p3 ,p6, p5)
    facetregion!(builder,5)
    facet!(builder, p3, p1 ,p4 ,p6)

	vol=0.5*W*L*H
	simplexgrid(builder, maxvolume=vol/(100.0*2.0^nref))
	
end

# ╔═╡ 7c6c81db-1920-49af-a101-462228614f95
#gridplot(prism(nref=1),azim=20,elev=20,linewidth=0.5,outlinealpha=0.3)
gridplot(prism(nref=1),Plotter=PlutoVista)

# ╔═╡ ed50c2d4-25e9-4159-84c7-e0c70ffa63a1
md"""
The prismatic geometry represents 1/8 of the square plate, using symmetry to make the modelling domain smaller. 

The 3D geomtry has 5 outer facets, whose boundary conditions need to be specified:
- facet 3: symmetry (no flux)
- facet 5: symmetry (no flux)
- facet 4: convective heat transfer to the wall, air gab between porous frit and Al reactor wall (robin bc.)
- facet 1: convective heat transfer to inflowing gas stream (robin bc.)
- facet 2: radiation + convection (custom bc.)
"""

# ╔═╡ 9d8c6ddc-2662-4055-b636-649565c36287
md"""
# Simulation
"""

# ╔═╡ ba5c2095-4858-444a-99b5-ae6cf40374f9
md"""
## 2D
"""

# ╔═╡ d725f9b9-61c4-4724-a1d9-6a04ba42499d
function main(;nref=0)
	data=ModelData()
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
		λbed=kbed(data)*λf

		conv=evelo[edge.index]*ρf*cf/λbed
		Bp,Bm = fbernoulli_pm(conv)
		
		#f[iT]= λbed*(Bm*u[iT,1]-Bp*u[iT,2])
		f[iT]= λbed*(Bm*(u[iT,1]-data.Tamb)-Bp*(u[iT,2]-data.Tamb))
		
		
	end

	function top(f,u,bnode,data)
		if bnode.region==3 # top boundary
			flux_rerad = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
			
			#flux_convec = data.α_nc*(u[iT]-data.Tamb)
			ρf=density_idealgas(data.Fluid, u[iT], data.p)
			cf=heatcap_gas(data.Fluid, u[iT])
            flux_convec = -bfvelo[bnode.ibnode,bnode.ibface]*ρf*cf*(u[iT]-data.Tamb)
			
			f[iT] = -data.Abs_lamp*data.G_lamp + flux_rerad + flux_convec
		end
	end

	function bottom(f,u,bnode,data)
		if bnode.region==1 # bottom boundary
            f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
		end
	end


	function bcondition(f,u,bnode,data)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=1,value=data.Tamb)
		#boundary_robin!(f,u,bnode;species=iT,region=1, factor=data.α_nc, value=data.Tamb*data.α_nc)
		bottom(f,u,bnode,data) # bottom boundary
		
		boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=3,value=data.Tamb+300.0)
		# irradiation boundary condition
		top(f,u,bnode,data)
	end
	

	
	grid=cylinder(;nref=nref,r=data.D/2,h=data.h)
	evelo=edgevelocities(grid,fup)
    bfvelo=bfacevelocities(grid,fup)
	
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
	inival[iT,:] .= map( (r,z)->(data.Tamb+500*z/data.h),grid)
	#inival[iT,:] .= data.Tamb
	sol=solve(inival,sys)
	sys,sol,data
end

# ╔═╡ e9003129-31db-4f9d-b289-638511c7ec26
Sim2D=main(nref=1);

# ╔═╡ f15fd785-010c-4fda-ab4f-7947642556dd
let
	sys,sol,data=Sim2D
	iT=data.iT
	vis=GridVisualizer(resolution=(800,800))
	solC = copy(sol)
	@. solC[iT,:] -= 273.15

	scalarplot!(vis,sys,solC;species=iT,title="Temperature / °C",xlabel="Radial coordinate / m", ylabel="Axial coordinate / m",legend=:best,colormap=:summer,show=true)

end

# ╔═╡ 28a2230f-5a59-4034-86af-e3d58dcceb6c
md"""
## 3D
"""

# ╔═╡ f435b4df-9162-42bd-8154-5f3434ea0e2a
function main3D(;nref=0)
	data=ModelData()
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
		Bp,Bm = fbernoulli_pm(conv)
		#f[iT]= λbed*(Bm*u[iT,1]-Bp*u[iT,2])
		f[iT]= λbed*(Bm*(u[iT,1]-data.Tamb)-Bp*(u[iT,2]-data.Tamb))
		
		
		
	end

	function irrad_bc(f,u,bnode,data)
		if bnode.region==6 # top boundary
			flux_rerad = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
			flux_convec = data.α_nc*(u[iT]-data.Tamb)
			f[iT] = -(data.Abs_lamp*data.G_lamp - flux_rerad - flux_convec)
		end
	end

	function top(f,u,bnode,data)
		if bnode.region==6 # top boundary
			#flux_rerad = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
			flux_rerad = 0
			
			ρf=density_idealgas(data.Fluid, u[iT], data.p)
			cf=heatcap_gas(data.Fluid, u[iT])
            flux_convec = data.u0*ρf*cf*(u[iT]-data.Tamb) # flow velocity is normal to top boundary
			
			f[iT] = -data.Abs_lamp*data.G_lamp + flux_rerad + flux_convec
		end
	end

	function bottom(f,u,bnode,data)
		if bnode.region==5 # bottom boundary
            f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
		end
	end

	function bcondition(f,u,bnode,data)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=1,value=data.Tamb)
		#boundary_robin!(f,u,bnode;species=iT,region=5, factor=data.α_nc, value=data.Tamb*data.α_nc) # bottom
		bottom(f,u,bnode,data)
		boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
		# for prism of square
		boundary_robin!(f,u,bnode;species=iT,region=3, factor=data.α_w, value=data.Tamb*data.α_w)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=2,value=data.Tamb+300.0)
		# irradiation boundary condition
		#irrad_bc(f,u,bnode,data)
		top(f,u,bnode,data)
	end
	

	
	#grid=prism(;nref=nref, w=data.wi, l=data.le, h=data.h)
	grid=prism_sq(data;nref=nref,)
	
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

# ╔═╡ ea81db48-d02c-4438-a0ba-9c54a3ea0e52
Sim3D=main3D(nref=0);

# ╔═╡ fb72aede-8997-48ef-9a2d-98acbf372747
md"""
f=$(@bind flevel Slider(range(([minimum(Sim3D[2]),maximum(Sim3D[2])].-273.15)... , length=40),default=(minimum(Sim3D[2])+maximum(Sim3D[2]))/2-273.15,show_value=true))

x=$(@bind xplane Slider(range(0.0,data.wi/2,length=20),default=0.0,show_value=true))

y=$(@bind yplane Slider(range(0.0,data.le/2,length=20),default=0.0,show_value=true))

z=$(@bind zplane Slider(range(0.0,data.h,length=20),default=0.0,show_value=true))

"""

# ╔═╡ 8107492e-fdaf-4a86-96ef-ea34b45c0e67
let
	sys,sol,data,nref=Sim3D
	iT=data.iT
	vis=GridVisualizer(resolution=(1000,1000),Plotter=PlutoVista)
	solC=copy(sol)
	@. solC[iT,:] -= 273.15

	scalarplot!(vis,sys,solC;species=iT,xlabel="X / m", ylabel="Y / m", zlabel="Z / m", levels=[flevel], xplanes=[xplane], yplanes=[yplane], zplanes=[zplane], levelalpha=1.0,title="Temperature / °C", outlinealpha=0.05)
	scene=reveal(vis)
	#save("plot_3D.svg",scene)
end

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
	grid=prism_sq(data;nref=nref)
	
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

	vis=GridVisualizer(resolution=(800,800))
	#scalarplot!(vis,grid_2D,sol_reorder;colormap=:summer,show=true)
	scalarplot!(vis,grid_2D,sol_cutplane.-sol_reorder;colormap=:bwr,show=true)

end

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
	# sys,sol,data=main(nref=1)
	sys,sol,data=Sim3D
	
	tff=VoronoiFVM.TestFunctionFactory(sys)
	
	# testfunction arguments:
		# 2nd arg: boundaries, where the test func. value is 0 -> these boundaries are excludes from the surface integral
		# 3rd arg: boundaries, where the test func. value is 1 -> these boundaries are included in the surface integral
	
	# flux integral over top boundary: radiation (abs+emit), convective trans
	tf=testfunction(tff,[1,2,3,4,5],6)
	integrate(sys,tf,sol)*4
end

# ╔═╡ 09632d67-a2cd-4c6f-b390-688e9900fc1e
gridplot(prism_sq(data,nref=0),)

# ╔═╡ Cell order:
# ╠═c258e0e6-72cc-4b3e-8f6d-e629e4a0d7fd
# ╠═7d8eb6f5-3ba6-46ef-8058-1f24a0938ed1
# ╠═5c3adaa0-9285-11ed-3ef8-1b57dd870d6f
# ╟─f353e09a-4a61-4def-ab8a-1bd6ce4ed58f
# ╠═2015c8e8-36cd-478b-88fb-94605283ac29
# ╠═5136828c-6e88-40fc-b9a5-3c86db2f0879
# ╠═98063329-31e1-4d87-ba85-70419beb07e9
# ╟─03d0c88a-462b-43c4-a589-616a8870be64
# ╟─6d5a7d83-53f9-43f3-9ccd-dadab08f62c1
# ╟─821aecc3-cf5d-4db2-bace-dab82c8cd711
# ╟─561a9857-d29e-40cd-809f-cc667613251c
# ╟─fdddb318-831c-4277-a5a3-b4be8f63c430
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
# ╟─16f5e0bc-8e3d-40cd-b67b-694eda6b67d9
# ╠═1459c3db-5ffc-46bd-9c94-8c8964519f39
# ╠═4eb00d7b-d10e-478d-a1df-9eea3362ef5f
# ╟─c9f326c1-abc0-412d-8348-8e23a7027661
# ╟─d7317b2d-e2c7-4114-8985-51979f2205ba
# ╟─3c75c762-a44c-4328-ae41-a5016ce181f1
# ╟─2c31e63a-cf42-45cd-b367-112438a02a97
# ╠═2fe11550-683d-4c4b-b940-3e63a4f8a87d
# ╠═8cd85a0e-3d11-4bcc-8a7d-f30313b31363
# ╟─a190862c-2251-4110-8274-9960c495a2c4
# ╠═4d9145f3-06aa-4a7c-82d0-feee0fa01865
# ╠═61a67079-cb15-4283-ac15-96b49c461b6e
# ╠═ebe2ed24-e2d6-4652-ae69-1a59747b0c4c
# ╠═7c6c81db-1920-49af-a101-462228614f95
# ╟─ed50c2d4-25e9-4159-84c7-e0c70ffa63a1
# ╟─9d8c6ddc-2662-4055-b636-649565c36287
# ╠═e9003129-31db-4f9d-b289-638511c7ec26
# ╟─ba5c2095-4858-444a-99b5-ae6cf40374f9
# ╟─f15fd785-010c-4fda-ab4f-7947642556dd
# ╠═d725f9b9-61c4-4724-a1d9-6a04ba42499d
# ╟─28a2230f-5a59-4034-86af-e3d58dcceb6c
# ╠═ea81db48-d02c-4438-a0ba-9c54a3ea0e52
# ╠═8107492e-fdaf-4a86-96ef-ea34b45c0e67
# ╠═fb72aede-8997-48ef-9a2d-98acbf372747
# ╠═f435b4df-9162-42bd-8154-5f3434ea0e2a
# ╟─64dd5097-16aa-4c44-b000-6177cd4be226
# ╟─bd179765-f996-4f91-ac4e-57d5817a2ed6
# ╠═3612d83d-a8bc-4c5f-8ea0-a4c975929500
# ╟─641988da-8888-4e0d-b720-4f21a9900aca
# ╠═b4175b62-198d-4605-9cf0-04b0be52c9c0
# ╟─e148212c-bc7c-4553-8ffa-26e6c20c5f47
# ╠═6a92c1c7-fb51-4363-b8d0-18eeb24087a8
# ╟─8310812b-97af-45d9-aff8-234ffbb169af
# ╟─909044b9-d40d-4129-bfef-6cd131a932f2
# ╠═30f0e723-a5a6-4b49-80d2-a29efb0f417b
# ╠═a1f5c0c2-290e-418f-b381-4626d61c2bfd
# ╠═09632d67-a2cd-4c6f-b390-688e9900fc1e
