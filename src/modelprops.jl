abstract type AbstractModelData end

# Base.@kwdef mutable struct ModelData <:AbstractModelData
# 	iT::Int64=1 # index of Temperature variable
#     ng::Int64=2 # number of gas phase species
	
	
	
# 	Tamb::Float64=298.15*ufac"K" # ambient temperature
# 	α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
# 	α_nc::Float64=15.0*ufac"W/(m^2*K)" # natural convection heat transfer coefficient

# 	## irradiation data
# 	G_lamp::Float64=1.0*ufac"kW/m^2" # solar simulator irradiation flux
# 	Abs_lamp::Float64=0.7 # avg absorptivity of cat. of irradiation coming from lamp
# 	Eps_ir::Float64=0.7 # avg absorptivity/emissivity of cat. of IR irradiation coming from surroundings / emitted
		
	
# 	## porous filter data
# 	d::Float64=100.0*ufac"μm" # average pore size
# 	# cylindrical disc / 2D
#     D::Float64=12.0*ufac"cm" # disc diameter
	

# 	# prism / 3D
# 	wi::Float64=12.0*ufac"cm" # prism width/side lenght
# 	le::Float64=wi # prism width/side lenght
# 	h::Float64=0.5*ufac"cm" # frit thickness (applies to 2D & 3D)

# 	#Ac::Float64=pi*D^2.0/4.0*ufac"m^2" # cross-sectional area, circular
# 	Ac::Float64=wi^2*ufac"m^2" # cross-sectional area, square
	
# 	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
# 	λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2 	
# 	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
	
# 	ϕ::Float64=0.36 # porosity, class 2
# 	k::Float64=2.9e-11*ufac"m^2" # permeability
# 	a_s::Float64=0.13*ufac"m^2/g" # specific surface area
# 	ρfrit::Float64=(1.0-ϕ)*ρs*ufac"kg/m^3" # density of porous frit
# 	a_v::Float64=a_s*ρfrit # volume specific interface area
# 	## END porous filter data

# 	## fluid data
	
# 	Qflow::Float64=3400.0*ufac"ml/minute" # volumetric feed flow rate
# 	Tin::Float64=298.15*ufac"K" # inlet temperature
# 	p::Float64=1.0*ufac"atm" # reactor pressure		
# 	# u0::Float64=Qflow/(Ac*ϕ)*ufac"m/s" # mean superficial velocity
# 	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
# 	# fluid properties: Air
# 	# values taken from VDI heat atlas 2010 chapter D3.1
# 	Fluid::FluidProps=Air
# 	## END fluid data
	
# end;

# calculate non-dimensional numbers: Reynolds, Prandtl, Peclet
function RePrPe(data::AbstractModelData,T,p,x)
	d=data.d
	u0=data.u0
    Fluids=data.Fluids
	ρf = density_idealgas(Fluids, T, p, x)
	ηf, λf = dynvisc_thermcond_mix(data, T, x)
    
    cf = heatcap_mix(Fluids, T, x)
	
	Re = u0*ρf*d/ηf # Reynolds number
	Pr = cf*ηf/λf # Prandtl number
	Pe = u0*ρf*cf*d/λf # Peclet
	Re,Pr,Pe
end

#"""
#Effective thermal conductivity of porous filter frit according to:
#
#__Zehner, P., & Schlünder, E. U. (1970).__ Wärmeleitfähigkeit von Schüttungen bei mäßigen Temperaturen. Chemie Ingenieur Technik, 42(14), 933-941. doi:10.1002/cite.330421408
#
#Implementation follows the notation in __VDI Heat Atlas 2010, ch. D6.3 eqs. (5a-5e).__
#"""
# calculate fixed bed (porous material) effective thermal conductivity
function kbed(data::AbstractModelData)
	(;ϕ,λs,Fluid,Tin) = data
	λf=thermcond_gas(Fluid, Tin)
	B=1.25*((1.0-ϕ)/ϕ)^(10.0/9.0)
	kp=λs/λf
	N=1.0-(B/kp)
	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	1.0-sqrt(1.0-ϕ)+sqrt(1.0-ϕ)*kc
end

function kbed(data::AbstractModelData,λf)
	(;ϕ,λs) = data
	B=1.25*((1.0-ϕ)/ϕ)^(10.0/9.0)
	kp=λs/λf
	N=1.0-(B/kp)
	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	1.0-sqrt(1.0-ϕ)+sqrt(1.0-ϕ)*kc
end

function kbed!(kbed,ϕ,λs,λf)
	#(;ϕ,λs) = data
	B=1.25*((1.0-ϕ)/ϕ)^(10.0/9.0)
	kp=λs/λf[1]
	N=1.0-(B/kp)
	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	kbed[1] = 1.0-sqrt(1.0-ϕ)+sqrt(1.0-ϕ)*kc
end

# mostly equivalent to VDI Heat Atlas 2010, ch. D6.3 eqs. (5a-5e). with addition
# of flattening coefficient ψ, that might be relevant for sintered materials
function kbed_VDI_flattening(data,λf)
	(;ϕ,ψ,λs) = data
	B=1.25*((1.0-ϕ)/ϕ)^(10.0/9.0)
	kp=λs/λf
	N=1.0-(B/kp)
	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	1 - sqrt(1-ϕ) + sqrt(1-ϕ)*(ψ*kp+(1-ψ)*kc)
end

# effective thermal conductivity for porous materials with random structure
# as proposed by Andrii Cheilytko (Andrii.Cheilytko@dlr.de)
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

#"""
#When working with a heterogeneous phase model (separate energy balances for both fluid and porous solid material), the exchange of energe between the phases can be described by an interfacial heat transfer coefficient. It can be calculated according to:
#
#__Kuwahara, F., Shirota, M., & Nakayama, A. (2001).__ A numerical study of interfacial convective heat transfer coefficient in two-energy equation model for convection in porous media. International Journal of Heat and Mass Transfer, 44(6), 1153-1159. doi:10.#1016/s0017-9310(00)00166-6
#
#```math
#\frac{h_{\text{sf}} \text D}{k_{\text f}}= \left( 1+ \frac{4(1- \phi)}{\phi} \right) + \frac{1}{2} (1-\phi)^{\frac{1}{2}} \text{Re}^{0.6}_D\text{Pr}^{\frac{1}{3}}
#```
#
#"""
function hsf(data::AbstractModelData,T,p,x)
	Re,Pr,_ = RePrPe(data,T,p,x)
	_,λf = dynvisc_thermcond_mix(data, T, x)
	ϕ = data.ϕ
	d = data.d
	λf/d*((1.0 + 4*(1.0-ϕ)/ϕ) + 0.5*(1.0-ϕ)^0.5*Re^0.6*Pr^(1.0/3.0))*ufac"W/(m^2*K)"
end

# Knudsen effective Diffusivity

# In the DGM the solid porous matrix is considered as another species in the ideal gas mix. The Knudsen diffusion coefficients describe the interactions between the gas molecules and the solid porous matrix. Analogous to the gas-gas diffusion coefficients, the Knudsen diffusion coefficients quantify the resulting friction force from momentum transfer upon collision with the walls acting on the gas molecules.

# Calculation of Knudsen diffusion coefficients according to __Wesselingh, J. A., & Krishna, R. (2006).__ Mass Transfer in Multicomponent Mixtures (ch. 21, Fig. 21.8)


function DK_eff(data,T,i)
	ϕ=data.ϕ # porosity
	dp=data.dp # avg. particle size (assumed to be = pore size)
	DK=dp*ϕ^1.5/(1.0-ϕ)*sqrt(8.0*ph"R"*T/(9.0*π*data.Fluids[i].MW))
	#Bern(DK)
end