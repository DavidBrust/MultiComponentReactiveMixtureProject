### A Pluto.jl notebook ###
# v0.19.24

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

# ╔═╡ 21925980-d94a-11ed-0a02-0fb6d3fd6888
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	
	using LessUnitful
	using Plots, LaTeXStrings
	using PlutoUI
	using Colors
	using Statistics

	using Revise
	using FixedBed
	
end;

# ╔═╡ 0fa3775d-4957-4644-8c9f-039a27b895b1
begin
	scalefontsizes()
	scalefontsizes(1.3)
	plot_font = "Computer Modern"
	default(fontfamily=plot_font,
		linewidth=2, framestyle=:box, label=nothing, grid=false)
end

# ╔═╡ f80fa3a7-a0f6-40b0-95eb-505502f770a6
begin
	const λeff_meas = [
		0.392*ufac"W/(m*K)",
		0.392*ufac"W/(m*K)",
		0.390*ufac"W/(m*K)",
		0.385*ufac"W/(m*K)",
		0.385*ufac"W/(m*K)",
		0.385*ufac"W/(m*K)"
	]
	λeff_meas_mean=mean(λeff_meas)
	λeff_meas_std=std(λeff_meas)
	ϕ_meas=0.36
end;

# ╔═╡ 5779b6f9-0d08-4111-8131-5a95c8290da2
md"""
__Schuetz, M.A. and L.R. Glicksman__, A Basic Study of Heat-Transfer through Foam Insulation. Journal of Cellular Plastics, 1984. 20(2): p. 114-121
"""

# ╔═╡ cb9ae049-3cff-4844-b85b-e465aaf58385
function lambda_eff_schuetz(data,λf)
	(;ϕ,λs) = data
	λf+2/3*(1-ϕ)*λs
end

# ╔═╡ 67690f27-6def-480e-8b54-3ca1db966dfc
md"""
Equation for tight backfill:

__Chudnovsky, A.F.__, thermal-physical characteristics of materials. 1962.
"""

# ╔═╡ 4927c598-997e-4237-91cc-08f39a7e896d
function lambda_eff_tight_backfill(data,λf)
	(;ϕ) = data
	phi=100*ϕ
	if phi-26 > 1e-3
		3*π*λf*log((43.0+0.31*phi)/(phi-26))
	else
		NaN
	end
end

# ╔═╡ 97814644-a5e9-4a7f-a4dd-b4f49286c220
md"""
The flattening coefficient ψ takes into account a deviation from point contacts between the particles that would result from packings of spheres. In a sintered material produced from partially melting individual particles, flat area contacts are to be expected.

__ψ = 0: point contacts, ψ large: large flat areas__

A large value of ψ would not be surprising in case of sintered materials.

The inclusion of the flattening coefficient ψ introduced a small small change compared to the standard model:
```math
	k_{\text{bed,flattening}} = 1-\sqrt{1-\phi} + \sqrt{1-\phi}(\psi k_{\text p} + (1-\psi) k_{\text c})
```
"""

# ╔═╡ 15e11af7-1f81-4d06-83b9-f7acc6c2ad97
function kbed_VDI_flattening(data,λf,ψ)
	(;ϕ,λs) = data
	B=1.25*((1.0-ϕ)/ϕ)^(10.0/9.0)
	kp=λs/λf
	N=1.0-(B/kp)
	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	1 - sqrt(1-ϕ) + sqrt(1-ϕ)*(ψ*kp+(1-ψ)*kc)
end

# ╔═╡ 81d9a193-3683-44d1-899a-810e02bdb1c7
md"""
The standard model:
```math
	k_{\text{bed}} = 1-\sqrt{1-\phi} + \sqrt{1-\phi}k_{\text c}
```
"""

# ╔═╡ 81bf69c2-7ef1-4964-a852-c13ea51b3785
function kbed_VDI(data,λf)
	(;ϕ,λs) = data
	B=1.25*((1.0-ϕ)/ϕ)^(10.0/9.0)
	kp=λs/λf
	N=1.0-(B/kp)
	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	1.0-sqrt(1.0-ϕ)+sqrt(1.0-ϕ)*kc
end

# ╔═╡ 8c13c9a2-3859-4eac-9732-adae7da2fab3
md"""
# Effective thermal conductivity AC 2
Andrii proposed a new calculation approcah for $\lambda_{\text{eff}}$ that was originally developed for open volumetric receivers with open channel structure.
"""

# ╔═╡ 12fffc7d-ae4d-4988-8051-a75a605b8ad2
md"""
```math
\lambda_{\text{eff}} = \frac{\frac{\lambda_f \lambda_s}{\lambda_s \phi + \lambda_f(1- \phi)}}{\frac{\lambda_s}{\lambda_f}\phi + 2 \psi + \frac{\lambda_f}{\lambda_s}(1-\phi)+1} \left[ \frac{\lambda_s}{\lambda_f} \phi +(1-\phi) + \psi \right] \left[ \phi + \frac{\lambda_f}{\lambda_s} (1-\phi) + \psi \right]

```
"""

# ╔═╡ 294e2e7c-6926-4ed6-a2a5-27ca024170b6
md"""
with the geometry dependent parameter
```math
\psi=\frac{1}{k_y}(\phi-1)\frac{\lambda_s-\lambda_f}{\lambda_s \lambda_f}[\lambda_f(\phi-1) -\lambda_s\phi]
```
"""

# ╔═╡ 26f7b04c-4e3a-42e5-ae27-c6ede74a77cb
function lambda_eff_AC(data,λf)
	(;ϕ,λs) = data
    # contact angle between particles in packing
	ky=1/sin(45*π/180)
    # initial version, developed for unordered foams
	# Ψ=ky*(1-ϕ)*(1-λf/λs+sqrt((1-λf/λs)^2+4))/2+ky*ϕ*(1-λs/λf+sqrt((1-λs/λf)^2+4))/2
    # version developed for open volumetric receivers with open channel structure
    Ψ=1/ky*(ϕ-1)*(λs-λf)/(λs*λf)*(λf*(ϕ-1)-λs*ϕ)
	λeff = (λs*λf/(λs*ϕ+λf*(1-ϕ)))*(ϕ*λs/λf+(1-ϕ)+Ψ)*(ϕ+λf/λs*(1-ϕ)+Ψ)/(ϕ*λs/λf+2*Ψ+(1-ϕ)*λf/λs+1)
end

# ╔═╡ 643b3f87-9283-4174-be00-12fcea4e9984
@bind ψ Slider(range(0,1,length=101), show_value=true, default=0.26)

# ╔═╡ b8835dc4-a71c-4e50-b090-6ece3bd8cc71
Base.@kwdef mutable struct ModelData <:AbstractModelData
	
	# catalyst / chemistry data
	kinpar::AbstractKineticsData = XuFroment1989
	
	# number of gas phase species
	ng::Int64		 		= kinpar.ng
	ip::Int64=ng+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable

	# names and fluid indices
	gn::Dict{Int, Symbol} 	= kinpar.gn

	# inverse names and fluid indices
	gni::Dict{Symbol, Int}  = kinpar.gni
	# fluids and respective properties in system
	Fluids::Vector{FluidProps} = kinpar.Fluids
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = let
		x=zeros(Float64, ng)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		x/sum(x)
	end # inlet composition
	
	## porous filter data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0
	#dp::Float64=100.0*ufac"μm" # average pore size, por class 2

	
	# Solid Boro-Silikatglas
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	#λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
	λs::Float64=1.13*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2

	ϕ::Float64=0.36 # porosity, exp determined
	#ϕ::Float64=0.33 # porosity, VitraPor sintetered filter class 0
	
	
	# approximation from Wesselingh, J. A., & Krishna, R. (2006). Mass Transfer in Multicomponent Mixtures
	γ_τ::Float64=ϕ^1.5 # constriction/tourtuosity factor

	#k::Float64=2.9e-11*ufac"m^2" # permeability , por class 2
	k::Float64=1.23e-10*ufac"m^2" # permeability , por class 0
	
	# a_s::Float64=0.13*ufac"m^2/g" # specific surface area, por class 2
	a_s::Float64=0.02*ufac"m^2/g" # specific surface area, por class 0

	
	ρfrit::Float64=(1.0-ϕ)*ρs*ufac"kg/m^3" # density of porous frit
	a_v::Float64=a_s*ρfrit*ufac"m^2/m^3" # volume specific interface area
	## END porous filter data


	## Flow data
	#norm conditions
	pn::Float64 = 1.0*ufac"bar"
	Tn::Float64 = 273.15*ufac"K"
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	p::Float64=1.0*ufac"atm" # reactor pressure
	
end;

# ╔═╡ 15d388f1-0258-481d-a902-3292f40b11b0
let
	data=ModelData()
	(;Tamb)=data
	#λf=thermcond_gas(Air, Tamb)
	λf=0.021
	lambda_eff_AC(data,λf)
	#kbed_VDI_flattening(data,λf,ψ)*λf
end

# ╔═╡ 3b0dd4c8-1a6c-48e0-9864-e3010c6f8a51
let
	lp=101
	phis=range(1.0e-6,1.0-1.0e-6,length=lp)
	λbed_VDI=zeros(lp)
	λbed_VDI_flattening=zeros(lp)
	λeff_AC=zeros(lp)
	#λeff_tight_backfill=zeros(lp)
	λeff_schuetz=zeros(lp)
	
	(;Tamb,λs,ϕ)=ModelData()
	λf=thermcond_gas(Air, Tamb)
	for (i,phi) in enumerate(phis)
		data=ModelData(ϕ=phi)
		λbed_VDI[i] = kbed_VDI(data,λf)*λf
		λbed_VDI_flattening[i] = kbed_VDI_flattening(data,λf,ψ)*λf
		λeff_AC[i] = lambda_eff_AC(data,λf)
		#λeff_tight_backfill[i] = lambda_eff_tight_backfill(data,λf)
		λeff_schuetz[i] = lambda_eff_schuetz(data,λf)
	end
	

	p=plot(xguide="Porosity / %", yguide=L"\lambda_\textrm{eff}~/~\mathrm{W}(\mathrm{m~K})^{-1}", legend=:outertopright)
	plot!(p,100*phis,λbed_VDI, label="VDI heat atlas", ls=:auto)
	plot!(p,100*phis,λbed_VDI_flattening, label="VDI heat atlas\nflattening", ls=:auto)
	plot!(p,100*phis,λeff_AC, label="Cheilytko", ls=:auto)
	#plot!(p,100*phis,λeff_tight_backfill, label="Chudnovsky\nTight Backfill", ls=:auto)
	plot!(p,100*phis,λeff_schuetz, label="Schuetz", ls=:auto)

	plot!(p,[100ϕ_meas],[λeff_meas_mean],yerror=[λeff_meas_std], m=:circle, label="λeff exp")

	lens!(p,[33, 39], [0.37, 0.41], inset = (1, bbox(0.3, 0.15, 0.3, 0.3)))

	hline!(p,100*phis,[λs,λs],ls=:dash,c=:black,label=:none)
	annotate!(p,88, 0.95*λs, text(L"\lambda_{\textrm{s}}="*string(λs),plot_font,12))

	hline!(p,100*phis,[λf,λf],ls=:dash,c=:black,label=:none)
	annotate!(p,13, 0.07, text(L"\lambda_{\textrm{f}}="*string(round(λf, sigdigits=2)),plot_font,12) )

	#savefig(p, "../img/out/lambda_eff_comp.pdf")
end

# ╔═╡ 2374882b-ef7e-47f3-911e-1d929febc851
ModelData()

# ╔═╡ Cell order:
# ╠═21925980-d94a-11ed-0a02-0fb6d3fd6888
# ╠═0fa3775d-4957-4644-8c9f-039a27b895b1
# ╠═f80fa3a7-a0f6-40b0-95eb-505502f770a6
# ╟─5779b6f9-0d08-4111-8131-5a95c8290da2
# ╠═cb9ae049-3cff-4844-b85b-e465aaf58385
# ╟─67690f27-6def-480e-8b54-3ca1db966dfc
# ╠═4927c598-997e-4237-91cc-08f39a7e896d
# ╟─97814644-a5e9-4a7f-a4dd-b4f49286c220
# ╠═15e11af7-1f81-4d06-83b9-f7acc6c2ad97
# ╠═81d9a193-3683-44d1-899a-810e02bdb1c7
# ╠═81bf69c2-7ef1-4964-a852-c13ea51b3785
# ╟─8c13c9a2-3859-4eac-9732-adae7da2fab3
# ╟─12fffc7d-ae4d-4988-8051-a75a605b8ad2
# ╠═294e2e7c-6926-4ed6-a2a5-27ca024170b6
# ╠═26f7b04c-4e3a-42e5-ae27-c6ede74a77cb
# ╠═15d388f1-0258-481d-a902-3292f40b11b0
# ╠═643b3f87-9283-4174-be00-12fcea4e9984
# ╠═3b0dd4c8-1a6c-48e0-9864-e3010c6f8a51
# ╠═2374882b-ef7e-47f3-911e-1d929febc851
# ╠═b8835dc4-a71c-4e50-b090-6ece3bd8cc71
