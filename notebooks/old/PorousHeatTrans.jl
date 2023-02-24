### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

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

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 7d8eb6f5-3ba6-46ef-8058-1f24a0938ed1
PlutoUI.TableOfContents(title="Heat Transfer in Fixed Beds")

# ╔═╡ f353e09a-4a61-4def-ab8a-1bd6ce4ed58f
md"""
# Porous filter disc
"""

# ╔═╡ 98063329-31e1-4d87-ba85-70419beb07e9
Base.@kwdef mutable struct ModelData
	iT::Int64=1 # index of Temperature variable
	
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	α_nc::Float64=15.0*ufac"W/(m^2*K)" # natural convection heat transfer coefficient

	## irradiation data
	G_lamp::Float64=50.0*ufac"kW/m^2" # solar simulator irradiation flux
	Abs_lamp::Float64=0.7 # avg absorptivity of cat. of irradiation coming from lamp
	Eps_ir::Float64=0.7 # avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
	
	
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
	Fluid::FluidProps=Air
	## END fluid data
	
end;

# ╔═╡ 03d0c88a-462b-43c4-a589-616a8870be64
md"""
# Experimental Conditions
"""

# ╔═╡ 3b3595c4-f53d-4827-918e-edcb74dd81f8
data = ModelData(;p=1.0*ufac"atm",Qflow=3400*ufac"ml/minute")

# ╔═╡ 2015c8e8-36cd-478b-88fb-94605283ac29
md"""
Specifications of porous filter disc from sintered silica glas (SiO₂): __VitraPOR P2__ (40-100 μm)

Since no value of average particle size is given, the particle size is in a first approximation assumed to equal the pore size of $(round(data.d/ufac"μm",sigdigits=1)) μm.
$(LocalResource("../img/filter1.png", :width => 1000))
$(LocalResource("../img/filter2.png", :width => 1000))
$(LocalResource("../img/filter3.png", :width => 1000))
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
$(LocalResource("../img/IrradBC.png", :width => 1000))
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

# ╔═╡ d912a1ca-1b69-4ea1-baa5-69794e004693
begin
	Ac=pi*data.D^2/4
	u0 = data.Qflow / (Ac*data.ϕ)
	# fluid properties: Air
	Fluid=data.Fluid
	ρf = density_idealgas(Fluid, data.Tin, data.p)
	cf = heatcap_gas(Fluid, data.Tin)
	λf = thermcond_gas(Fluid, data.Tin)
	# Peclet number
	Pe0 = u0*ρf*cf*data.d/λf

	# Reynolds number
	ηf=dynvisc_gas(Fluid, data.Tin)
	Re0 = u0*ρf*data.d/ηf
	# Prandtl number
	Pr = cf*ηf/λf
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

# ╔═╡ 0cff03d0-00ca-4026-b041-ad9de2908d87
md"""
At given experimental conditions the Reynolds, Prandtl and Peclet numbers assume the following values:
- Re = $(round(Re0,sigdigits=2))
- Pr = $(round(Pr,sigdigits=2))
- Pe = $(round(Pe0,sigdigits=2))
"""

# ╔═╡ 0510e643-8eab-4423-a0e2-aacd0f83a0c4
md"""
## Interfacial heat transfer coefficient
"""

# ╔═╡ 239212da-7911-4a51-9b22-08cd07699d4f
md"""
When working with a heterogeneous phase model (separate energy balances for both fluid and porous solid material), the exchange of energe between the phases can be described by an interfacial heat transfer coefficient. It can be calculated according to:

__Kuwahara, F., Shirota, M., & Nakayama, A. (2001).__ A numerical study of interfacial convective heat transfer coefficient in two-energy equation model for convection in porous media. International Journal of Heat and Mass Transfer, 44(6), 1153-1159. doi:10.1016/s0017-9310(00)00166-6

```math
\frac{h_{\text{sf}} \text D}{k_{\text f}}= \left( 1+ \frac{4(1- \phi)}{\phi} \right) + \frac{1}{2} (1-\phi)^{\frac{1}{2}} \text{Re}^{0.6}_D\text{Pr}^{\frac{1}{3}}
```

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
    step=0.1*ufac"cm"*2.0^(-nref)
    R=collect(0:step:r)
    Z=collect(0:step:h)
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
	simplexgrid(builder, maxvolume=vol/100.0)
	
end

# ╔═╡ 7c6c81db-1920-49af-a101-462228614f95
gridplot(prism(),azim=20,elev=20,linewidth=0.5,outlinealpha=0.3)

# ╔═╡ 9d8c6ddc-2662-4055-b636-649565c36287
md"""
# Simulation
"""

# ╔═╡ d725f9b9-61c4-4724-a1d9-6a04ba42499d
function main(;nref=0,p=1.0*ufac"atm",Qflow=3400*ufac"ml/minute")
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

	function irrad_bc(f,u,bnode,data)
		if bnode.region==3 # top boundary
			flux_rerad = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
			flux_convec = data.α_nc*(u[iT]-data.Tamb)
			f[iT] = -(data.Abs_lamp*data.G_lamp - flux_rerad - flux_convec)
		end
	end

	function bcondition(f,u,bnode,data)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=1,value=data.Tamb)
		boundary_robin!(f,u,bnode;species=iT,region=1, factor=data.α_nc, value=data.Tamb*data.α_nc)
		boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
		#boundary_dirichlet!(f,u,bnode;species=iT,region=3,value=data.Tamb+300.0)
		# irradiation boundary condition
		irrad_bc(f,u,bnode,data)
	end
	

	
	grid=cylinder(;nref=nref,r=data.D/2,h=data.h)
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
	inival[iT,:] .= map( (r,z)->(data.Tamb+500*z/data.h),grid)
	#inival[iT,:] .= data.Tamb
	sol=solve(inival,sys)
	sys,sol,data
end

# ╔═╡ f15fd785-010c-4fda-ab4f-7947642556dd
let
	sys,sol,data=main()
	iT=data.iT
	vis=GridVisualizer()
	@. sol[iT,:] -= 273.15
	scalarplot!(vis,sys,sol;species=iT,title="Temperature / °C",xlabel="Radial coordinate / m", ylabel="Axial coordinate / m",legend=:best,colormap=:summer,show=true)

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
	plot(Ts.-273.15, λf,xlabel="Temperature / °C",ylabel="Thermal Conductivity / W m⁻¹ K⁻¹", title=Fluid.name )
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
# ╠═98063329-31e1-4d87-ba85-70419beb07e9
# ╟─03d0c88a-462b-43c4-a589-616a8870be64
# ╟─6d5a7d83-53f9-43f3-9ccd-dadab08f62c1
# ╟─3b3595c4-f53d-4827-918e-edcb74dd81f8
# ╠═d912a1ca-1b69-4ea1-baa5-69794e004693
# ╟─0cff03d0-00ca-4026-b041-ad9de2908d87
# ╟─cb6a357f-e244-4725-a04a-3e006dd4b53d
# ╟─463a9a2b-8437-407f-b31a-dde3165f49ad
# ╟─387b5b8e-a466-4a65-a360-fa2cf08092e3
# ╟─b4cee9ac-4d14-4169-90b5-25d12ac9c003
# ╟─bbd0b076-bcc1-43a1-91cb-d72bb17d3c88
# ╠═7c4d4083-33e2-4b17-8575-214ef458d75d
# ╟─0510e643-8eab-4423-a0e2-aacd0f83a0c4
# ╟─239212da-7911-4a51-9b22-08cd07699d4f
# ╟─d7317b2d-e2c7-4114-8985-51979f2205ba
# ╟─3c75c762-a44c-4328-ae41-a5016ce181f1
# ╟─2c31e63a-cf42-45cd-b367-112438a02a97
# ╠═2fe11550-683d-4c4b-b940-3e63a4f8a87d
# ╠═8cd85a0e-3d11-4bcc-8a7d-f30313b31363
# ╟─a190862c-2251-4110-8274-9960c495a2c4
# ╠═ebe2ed24-e2d6-4652-ae69-1a59747b0c4c
# ╠═7c6c81db-1920-49af-a101-462228614f95
# ╟─9d8c6ddc-2662-4055-b636-649565c36287
# ╠═f15fd785-010c-4fda-ab4f-7947642556dd
# ╠═d725f9b9-61c4-4724-a1d9-6a04ba42499d
# ╟─e148212c-bc7c-4553-8ffa-26e6c20c5f47
# ╠═6a92c1c7-fb51-4363-b8d0-18eeb24087a8
# ╟─8310812b-97af-45d9-aff8-234ffbb169af
# ╟─909044b9-d40d-4129-bfef-6cd131a932f2
# ╠═30f0e723-a5a6-4b49-80d2-a29efb0f417b
# ╠═a1f5c0c2-290e-418f-b381-4626d61c2bfd
