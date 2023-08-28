abstract type AbstractModelData end


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
	kbed=1.0-sqrt(1.0-ϕ)+sqrt(1.0-ϕ)*kc
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

#
# ## Irradiation
# Irradiation flux coming from the solar simulator enters the aperture of the reactor with the specified flux profile as determined via photo-metric measurement. The measured profile is imported from a csv-datafile, handled by DataFrames.jl and interpolated by methodes provided by Interpolations.jl to be used as boundary condition in the simulation.
#

# domain width (porous frit = 15.7 cm, ~ 16.0)
function sel12by12(M;wi=16.0*ufac"cm",wi_apt=12.0*ufac"cm")
	
	
	
	# starting coordinates for optimum 10 cm x 10 cm selection
	sr=33
	sc=33

	# (inverse) resolution of flux measurements, distance between data points
	Dx = 0.031497*ufac"cm"
	Dy = 0.031497*ufac"cm"
	
	
	# calculate coordinate offsets, when deviating (expanding/shrinking from the center) from 10cm x 10cm selection 	
	Dsr=Integer(round((wi_apt-10.0*ufac"cm")/(2*Dy)))
	Dsc=Integer(round((wi_apt-10.0*ufac"cm")/(2*Dx)))

	sr-=Dsr
	sc-=Dsc

	nx=Integer(round(wi_apt/Dx))
	ny=Integer(round(wi_apt/Dy))
	
	M=Matrix(FluxMap)
	# coordinate system of measurement data matrix has its origin at bottom left corner, the excel data matrix (and its coordiinates) start in top left corner
	reverse!(M; dims=1) 
	@views M_ = M[sr:(sr+ny-1),sc:(sc+nx-1)]*ufac"kW/m^2"

	# pad Flux Map with zeros outside of aperture area
	nx_dom = Integer(round(wi/Dx))
	ny_dom = Integer(round(wi/Dy))
	M__ = zeros(nx_dom,ny_dom)
		
	Dsr_dom_apt=Integer(round((wi - wi_apt)/(2*Dy)))
	Dsc_dom_apt=Integer(round((wi - wi_apt)/(2*Dx)))

	M__[(Dsr_dom_apt+1):(Dsr_dom_apt+ny), (Dsc_dom_apt+1):(Dsc_dom_apt+nx)] = M[sr:(sr+ny-1),sc:(sc+nx-1)]*ufac"kW/m^2"
	#Dsc=Integer(round((wi_apt-10.0*ufac"cm")/(2*Dy)))
	
	# origin of coordinate system in the center of the plane
	#x = range(-wi/2,wi/2,length=nx)
	#y = range(-wi/2,wi/2,length=ny)

	#itp = Interpolations.interpolate((x,y), M_, Gridded(Linear()))
	x = range(-wi/2,wi/2,length=nx_dom)
	y = range(-wi/2,wi/2,length=ny_dom)

	itp = Interpolations.interpolate((x,y), M__, Gridded(Linear()))

	#M_,itp	
	M__,itp	
end


