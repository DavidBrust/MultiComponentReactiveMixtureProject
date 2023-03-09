### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 8d82851a-4313-4cb9-aa5a-8ee15ac40657
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))	
	
	using LinearAlgebra
	using LessUnitful
	
	using Plots
	using PlutoUI

	using Revise
	using FixedBed
	
end;

# ╔═╡ 50237926-dd49-45db-8ebe-63d9ba98dc9a
md"""
Irradiation exchange in upper chamber of PC reactor. The inside glass surface facing towards the catalyst and the catalyst itself are designated as surfaces 1 and 2 respectively. The two surfaces are coupled through irradiation exchange. The two unknows are the surface temperatures $T_1$ and $T_2$.

Only the highlighted parameters must be specified explicitly, the remaining parameters can be determined from those.
The following table lists the necessary optical properties. 
"""

# ╔═╡ b17ef005-cf0f-4a7d-86b1-0c10ce60ff3d
begin
	# index 1 = glass
	const alpha1_IR = 0.2 # obtain from datasheet, absorptance
	const tau1_IR = 0.5 # obtain from datasheet, transmittance
	const eps1 = alpha1_IR
	const rho1_IR = 1.0 - tau1_IR - alpha1_IR
	const alpha1_vis = 0.0
	const tau1_vis = 0.9 # obtain from datasheet, transmittance
	# the following is probably wrong since it depends on incidence angle, which differs between the irradiation coming from the lamp (relevant for tau1_vis) and the diffuse reflected light coming from the catalyst (surface 2)
	const rho1_vis = 1.0 - tau1_vis - alpha1_vis 

	# index 2 = catalyst surface
	const eps2 = 0.7 # obtain from measurement
	const alpha2_IR = eps2
	const tau2_IR = 0.0
	const rho2_IR = 1.0 - alpha2_IR - tau2_IR
	const alpha2_vis = 0.7 # obtain from measurement
	const tau2_vis = 0.0
	const rho2_vis = 1.0 -alpha2_vis - tau2_vis
	
end;

# ╔═╡ a4caa320-bd20-11ed-0e42-693f0a5328c9
md"""
|  glass (1) |   |   catalyst (2)|    |
|-------|----|---|-------|
| $\epsilon_1$ | $(eps1) | $*$ $\epsilon_2$  | $(eps2)     |
| $*$ $\alpha_1^{\text{IR}}$   | $(alpha1_IR) | $\alpha_2^{\text{IR}}$ |   $(alpha2_IR)  |
| $*$ $\tau_1^{\text{IR}}$   | $(tau1_IR) | $\tau_2^{\text{IR}}$ |   $(tau2_IR)  |
| $\rho_1^{\text{IR}}$   | $(rho1_IR) | $\rho_2^{\text{IR}}$ |   $(round(rho2_IR,sigdigits=1))  |
|||||
| $\alpha_1^{\text{vis}}$   | $(alpha1_vis) | $*$ $\alpha_2^{\text{vis}}$ |   $(alpha2_vis)  |
| $*$ $\tau_1^{\text{vis}}$   | $(round(tau1_vis,sigdigits=1)) | $\tau_2^{\text{vis}}$ |   $(tau2_vis)  |
| $\rho_1^{\text{vis}}$   | $(round(rho1_vis,sigdigits=1)) | $\rho_2^{\text{vis}}$ |   $(round(rho2_vis,sigdigits=1))  |
"""

# ╔═╡ 07a18632-4c57-4820-815d-24a3afa71347
Base.@kwdef struct IrradiationExchange
	# surface 1
	# effective optical parameters intergrated over infrared spectrum (IR)
	alpha1_IR::Float64 = 0.2 # absorptance: obtain from datasheet/measurement
	tau1_IR::Float64 = 0.5 # transmittance: obtain from datasheet/measurement
	rho1_IR::Float64 = 1.0 - tau1_IR - alpha1_IR
	eps1::Float64 = alpha1_IR
	# effective opt. par. integrated over visible spectrum (vis)
	alpha1_vis::Float64 = 0.0
	tau1_vis::Float64 = 0.9 # obtain from datasheet, transmittance
	# the following is probably inaccurate since it depends on incidence angle, which differs between the irradiation coming from the lamp (relevant for tau1_vis) and the diffuse reflected light coming from the catalyst (surface 2)
	rho1_vis::Float64 = 1.0 - tau1_vis - alpha1_vis 

	# surface 2
	# effective optical parameters intergrated over infrared spectrum (IR)
	alpha2_IR::Float64 = 0.7 # obtain from datasheet, absorptance
	tau2_IR::Float64 = 0.0 # obtain from datasheet, transmittance
	rho2_IR::Float64 = 1.0 - tau2_IR - alpha2_IR
	eps2::Float64 = alpha2_IR
	# effective opt. par. integrated over visible spectrum (vis)
	alpha2_vis::Float64 = 0.7
	tau2_vis::Float64 = 0.0 
	rho2_vis::Float64 = 1.0 - tau2_vis - alpha2_vis 
	
end

# ╔═╡ 67f418db-5c38-4e30-86d1-ef2a6258b9c9
# surface 1: glass
upper_chamber = IrradiationExchange(
	
)

# ╔═╡ 6c246116-d73b-4c2d-bb4b-577d2578f2ab
const Glamp = 50.0*ufac"kW/m^2";

# ╔═╡ e5b9e84d-880a-485a-8429-12e5584aa8e4
begin
	denum = 1-rho1_IR*rho2_IR
	F = zeros(2)
	F[1] = 0.0
	F[2] = alpha2_vis*tau1_vis/(1.0-rho1_vis*rho2_vis)*Glamp
	M = zeros(2,2)
	
	M[1,1] = eps1*ph"σ" - alpha1_IR*rho2_IR*eps1*ph"σ"/denum
	M[1,2] = -alpha1_IR*eps2*ph"σ"*(1+rho1_IR*rho2_IR/denum)
	M[2,1] = -alpha2_IR*eps1*ph"σ"/denum
	M[2,2] = eps2*ph"σ"*(1-alpha2_IR*rho1_IR/denum)

	T = M \ F
	T .^(1/4) .- 273.15
end

# ╔═╡ Cell order:
# ╠═8d82851a-4313-4cb9-aa5a-8ee15ac40657
# ╟─50237926-dd49-45db-8ebe-63d9ba98dc9a
# ╟─a4caa320-bd20-11ed-0e42-693f0a5328c9
# ╠═b17ef005-cf0f-4a7d-86b1-0c10ce60ff3d
# ╠═07a18632-4c57-4820-815d-24a3afa71347
# ╠═67f418db-5c38-4e30-86d1-ef2a6258b9c9
# ╠═6c246116-d73b-4c2d-bb4b-577d2578f2ab
# ╠═e5b9e84d-880a-485a-8429-12e5584aa8e4
