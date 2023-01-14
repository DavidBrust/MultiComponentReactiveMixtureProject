### A Pluto.jl notebook ###
# v0.19.19

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
	Fluid::FluidProps=Air
	## END fluid data
	
end;

# ╔═╡ 03d0c88a-462b-43c4-a589-616a8870be64
md"""
# Experimental Conditions
"""

# ╔═╡ 3b3595c4-f53d-4827-918e-edcb74dd81f8
data = ModelData(;p=1.0*ufac"atm",Qflow=3400*ufac"ml/minute")

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

# ╔═╡ 8cd85a0e-3d11-4bcc-8a7d-f30313b31363
gridplot(cylinder(;r=data.D/2,h=data.h))

# ╔═╡ 9d8c6ddc-2662-4055-b636-649565c36287
md"""
# Simulation
"""

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
