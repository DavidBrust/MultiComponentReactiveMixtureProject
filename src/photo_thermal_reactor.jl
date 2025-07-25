# abstract type AbstractReactorData end

@doc raw"""
Helper function to calculate radiation emitted from the window towards the surface
    of the catalyst layer in the upper chamber in the photo thermal catalytic
    reactor (PTR) model.
"""
function PTR_radiosity_window(f,u,bnode,data)
    # (;dim,iTw,G_lamp,nom_flux,FluxIntp,Tamb,uc_window,irradiated_boundaries)=data
	(;dim,iTw,ibf,nom_flux,Tamb,uc_window,irradiated_boundaries)=data
	
    # irrad. exchange between quartz window (1), cat surface (2), masked sruface (3) 
    tau1_vis=uc_window.tau_vis
    rho1_vis=uc_window.rho_vis
    tau1_IR=uc_window.tau_IR
    rho1_IR=uc_window.rho_IR
    eps1=uc_window.eps
	
	# G_lamp = u[ibf]
	G_lamp = zero(eltype(u))
	if dim == 2
		G_lamp += nom_flux
	elseif dim == 3
		# obtain local irradiation flux value from boundary species
		# !!! DEBUG !!!
		# G_lamp += nom_flux
		G_lamp += u[ibf]
		# !!! DEBUG !!!
	end
	

    Tglass = u[iTw] # local tempererature of quartz window
	G1_bot_IR = eps1*σ*(Tglass^4-Tamb^4)+σ*Tamb^4
	#G1_bot_IR = σ*Tamb^4
	G1_bot_vis = zero(eltype(u))
    if bnode.region in irradiated_boundaries		
		G1_bot_vis += G_lamp # flux profile measured behind quarz in plane of cat layer
    end
    return G1_bot_vis,G1_bot_IR
end

@doc raw"""
Function defining the top/inlet boundary condition in the photo thermal catalytic
    reactor (PTR) model.
"""
function PTR_top(f,u,bnode,data)
	(;ip,iT,iTw,Tamb,T_gas_in,mfluxin,X0,W0,mmix0,uc_window,uc_cat,uc_h,Nu,k_nat_conv,dt_mf,dt_hf_enth,dt_hf_irrad,inlet_boundaries,irradiated_boundaries,solve_T_equation,Fluids)=data
	ng=ngas(data)

	# irrad. exchange between quartz window (1), cat surface (2), masked surface (3)
	alpha1_vis=uc_window.alpha_vis
	alpha1_IR=uc_window.alpha_IR
	
	rho2_vis=uc_cat.rho_vis
	rho2_IR=uc_cat.rho_IR
	alpha2_vis=uc_cat.alpha_vis
	alpha2_IR=uc_cat.alpha_IR
	eps2=uc_cat.eps
	
	if bnode.region in inlet_boundaries
		r_mfluxin = mfluxin*ramp(bnode.time; du=(0.1,1), dt=dt_mf)
		
		f[ip] = -r_mfluxin # Neumann bc. for total mass flux
		for i=1:(ng-1)
			f[i] = -r_mfluxin*W0[i] # Neumann bc. for species mass flux
		end

	end

	if solve_T_equation
		if bnode.region in inlet_boundaries
			# heatflux from enthalpy inflow
			# !!! DEBUG !!!
            # @inline r_hf_enth = mfluxin/mmix0 * enthalpy_mix(data, T_gas_in, X0) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_enth)
			# Tfilm = (T_gas_in+u[iT])
			# @inline hin = heatcap_mix(data, Tfilm, X0) * (Tfilm - 298.15)
			hin = zero(eltype(u))
			for i=1:ng
				hin += enthalpy_gas_thermal(Fluids[i], T_gas_in)
			end
			# @inline hin = heatcap_mix(data, T_gas_in, X0) * (T_gas_in - Tref)
			r_hf_enth = mfluxin/mmix0 * hin
			f[iT] = -r_hf_enth * ramp(bnode.time; du=(0.0,1), dt=dt_hf_enth)
			# !!! DEBUG !!!
			
			if bnode.region in irradiated_boundaries
				@inline G1_vis, G1_IR = PTR_radiosity_window(f,u,bnode,data)

				hflux_irrad = (-eps2*σ*u[iT]^4 + alpha2_vis*G1_vis + alpha2_IR*G1_IR) # irradiation heatflux 
				#hflux_irrad = (-eps2*σ*(u[iT]^4-Tamb^4) + alpha2_vis*G1_vis) # irradiation heatflux 
				
				##################### NEGLECT INTERNAL CONVECTION #####################
				# Tm=0.5*(u[iT] + u[iTw]) # mean temperature: catalyst layer and window
		        # @inline _,λf=dynvisc_thermcond_mix(data, Tm, X0) # thermal conductivity at Tm and inlet composition X0
				# dh=2*uc_h
				# kconv=Nu*λf/dh*ufac"W/(m^2*K)"
				# hflux_conv = kconv*(u[iT]-Tm) # convection heatflux through top chamber
				##################### NEGLECT INTERNAL CONVECTION #####################
				
				# !!! DEBUG !!!
				# f[iT] += (-hflux_irrad + hflux_conv) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)
				f[iT] -= hflux_irrad * ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)
				# f[iT] -= hflux_irrad 
				# !!! DEBUG !!!
				
				# f[iT] += zero(eltype(u))
				# f[iT] += (-hflux_irrad) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)

				
				
				# calculate local window temperature from (local) flux balance
				hflux_conv_top_w = k_nat_conv*(u[iTw]-Tamb)
				G2_vis = rho2_vis*G1_vis		
				G2_IR = eps2*σ*u[iT]^4 + rho2_IR*G1_IR
				hflux_abs_w = alpha1_vis*G2_vis + alpha1_IR*G2_IR + alpha1_IR*σ*Tamb^4
				hflux_emit_w = uc_window.eps*σ*u[iTw]^4
				# !!! DEBUG !!!
				# f[iTw] = -hflux_conv -hflux_abs_w +hflux_conv_top_w +2*hflux_emit_w
				# f[iTw] *= ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)
				f[iTw] = -hflux_abs_w +hflux_conv_top_w +2*hflux_emit_w
                # f[iTw] = zero(eltype(u))
				# !!! DEBUG !!!
			end
		end
		# !!! DEBUG !!!
		# f[iT] = u[iT] - (273.15 + 650.0)
		# !!! DEBUG !!!
	end
	
end

@doc raw"""
Function defining the side boundary condition in the photo thermal catalytic
    reactor (PTR) model.
"""
function PTR_side(f,u,bnode,data)
	(;iT,k_nat_conv,delta_gap,Tamb,side_boundaries,solve_T_equation,dt_hf_enth)=data
	ng=ngas(data)

	if solve_T_equation && bnode.region in side_boundaries
	    # X=MVector{ng,eltype(u)}(undef)
		X = get_tmp(data.X, u); X = @view X[1:ng]
        MoleFrac!(X,u,data)
        # @inline _,λf=dynvisc_thermcond_mix(data, u[iT], X)
		mu = get_tmp(data.dynvisc, u); mu = @view mu[1:ng]
		lambda = get_tmp(data.thermcond, u); lambda = @view lambda[1:ng]

		_, λf = dynvisc_thermcond_mix(data, mu, lambda, u[iT], X)


		# !!! DEBUG !!!
        # w/o shell height
        f[iT] = (u[iT]-Tamb)/(delta_gap/λf+1/k_nat_conv) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_enth)
        # f[iT] = zero(eltype(u))		
		# !!! DEBUG !!!
    end
end

@doc raw"""
Function defining the bottom/outlet boundary condition in the photo thermal catalytic
    reactor (PTR) model.
"""
function PTR_bottom(f,u,bnode,data)
	(;iT,iTp,k_nat_conv,Tamb,lc_h,Nu,dt_hf_irrad,lc_frit,lc_plate,solve_T_equation,outlet_boundaries) = data
	ng=ngas(data)
	
	if solve_T_equation && bnode.region in outlet_boundaries

		# irradiation heat flux
		rho1_IR=lc_frit.rho_IR
		alpha1_IR=lc_frit.alpha_IR
		eps1=lc_frit.eps
		rho2_IR=lc_plate.rho_IR
		alpha2_IR=lc_plate.alpha_IR
		eps2=lc_plate.eps

		# heatflux from irradiation exchange with bottom plate
		Tplate = u[iTp]	
		hflux_irrad = -eps1*σ*u[iT]^4 + alpha1_IR/(1-rho1_IR*rho2_IR)*(eps2*σ*Tplate^4+rho2_IR*eps1*σ*u[iT]^4)
	

		##################### NEGLECT INTERNAL CONVECTION #####################
		# heatflux from convective/conductive heat exchange with bottom plate
		# Tm=0.5*(u[iT] + Tplate)

		# X=MVector{ng,eltype(u)}(undef)
		# @inline MoleFrac!(X,u,data)
        
        # @inline _,λf=dynvisc_thermcond_mix(data, Tm, X) # thermal conductivity at Tm and outlet composition X
        # dh=2*lc_h
		# kconv=Nu*λf/dh*ufac"W/(m^2*K)"
		# hflux_conv = kconv*(u[iT]-Tm) # positive flux in negative z coord. (towards plate)
		##################### NEGLECT INTERNAL CONVECTION #####################


		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative

		# !!! DEBUG
        # f[iT] = (-hflux_irrad + hflux_conv) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)
		f[iT] = -hflux_irrad * ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)
		# f[iT] = -hflux_irrad
		# f[iT] = zero(eltype(u))
		# !!! DEBUG

		# f[iT] = (-hflux_irrad) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)
		# f[iT] = zero(eltype(u))

		# calculate (local) plate temperature from (local) flux balance
		hflux_conv_bot_p = k_nat_conv*(u[iTp]-Tamb)
		hflux_emit_p = lc_plate.eps*σ*u[iTp]^4
		G1_IR = (eps1*σ*u[iT]^4 + rho1_IR*eps2*σ*u[iTp]^4)/(1-rho1_IR*rho2_IR)
		hflux_abs_p = alpha2_IR*G1_IR + alpha2_IR*σ*Tamb^4

		# !!! DEBUG
		# f[iTp] = -hflux_conv -hflux_abs_p +2*hflux_emit_p +hflux_conv_bot_p
		# f[iTp] *= ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)
		f[iTp] = -hflux_abs_p + 2*hflux_emit_p +hflux_conv_bot_p
        # f[iTp] = zero(eltype(u))
		# !!! DEBUG

	end
end

@doc raw"""
Wrapper function defining the boundary conditions in the photo thermal catalytic
    reactor (PTR) model. Used with the VoronoiFVM.jl package
"""
function PTR_bcond(f,u,bnode,data)
	(;p,ip,outlet_boundaries)=data
	
    PTR_top(f,u,bnode,data)
    PTR_side(f,u,bnode,data)
    PTR_bottom(f,u,bnode,data)
    for boundary in outlet_boundaries
        boundary_dirichlet!(f,u,bnode, species=ip,region=boundary,value=p)
    end	
end

@doc raw"""
Function defining the flux function in the boundary species in the photo thermal catalytic
    reactor (PTR) model. Calculates the heat flux occuring in the window and 
    bottom plate respectively. Only relevant if temperature equation is solved.
"""
function PTR_bflux(f,u,bedge,data)
	(;irradiated_boundaries,outlet_boundaries)=data

	if bedge.region in irradiated_boundaries # window, upper chamber
		(;iTw,lambda_window) = data
		f[iTw] = lambda_window * (u[iTw, 1] - u[iTw, 2])
	elseif bedge.region in outlet_boundaries # bottom plate, lower chamber
		(;iTp,lambda_Al) = data
		f[iTp] = lambda_Al * (u[iTp, 1] - u[iTp, 2])
	end
end

@doc raw"""
Function defining the storage function of the boundary species in the photo
    thermal catalytic reactor (PTR) model. Calculates the heat storage capacity
    for transient simulations in the window and bottom plate respectively.
    Only relevant if temperature equation is solved.
"""
function PTR_bstorage(f,u,bnode,data)
	(;irradiated_boundaries,outlet_boundaries,solve_T_equation)=data
	if solve_T_equation
		if bnode.region in irradiated_boundaries # window, upper chamber
			(;iTw) = data
			f[iTw] = u[iTw]
		elseif bnode.region in outlet_boundaries # bottom plate, lower chamber
			(;iTp) = data
			f[iTp] = u[iTp]
		end
	end
end


#"""
#Effective thermal conductivity of porous filter frit according to:
#
#__Zehner, P., & Schlünder, E. U. (1970).__ Wärmeleitfähigkeit von Schüttungen bei mäßigen Temperaturen. Chemie Ingenieur Technik, 42(14), 933-941. doi:10.1002/cite.330421408
#
#Implementation follows the notation in __VDI Heat Atlas 2010, ch. D6.3 eqs. (5a-5e).__
#"""
# calculate fixed bed (porous material) effective thermal conductivity


# function kbed(data, lambdaf)
# 	(;poros,lambdas) = data
# 	# !!! regularize function to catch poros=1 (pure gas phase) !!!
# 	ϵ_reg = 1.0e-12
# 	B=1.25*((1.0-poros)/poros)^(10.0/9.0)
# 	kp=lambdas/lambdaf
# 	N=1.0-(B/kp)
# 	#kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
# 	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/(B+ϵ_reg)) - (B+1.0)/2.0 - (B-1.0)/N)
# 	1.0-sqrt(1.0-poros)+sqrt(1.0-poros)*kc
# end

# mostly equivalent to VDI Heat Atlas 2010, ch. D6.3 eqs. (5a-5e). with addition
# of flattening coefficient ψ, that might be relevant for sintered materials
function kbed_VDI_flattening(data,lambdaf)
	(;poros,ψ,lambdas) = data
	B=1.25*((1.0-poros)/poros)^(10.0/9.0)
	kp=lambdas/lambdaf
	N=1.0-(B/kp)
	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	1 - sqrt(1-poros) + sqrt(1-poros)*(ψ*kp+(1-ψ)*kc)
end

# effective thermal conductivity for porous materials with random structure
# as proposed by Andrii Cheilytko (Andrii.Cheilytko@dlr.de)
function lambda_eff_AC(data,lambdaf)
	(;poros,lambdas) = data
    # contact angle between particles in packing
	ky=1/sin(45*π/180)
    # initial version, developed for unordered foams
	# Ψ=ky*(1-poros)*(1-lambdaf/lambdas+sqrt((1-lambdaf/lambdas)^2+4))/2+ky*poros*(1-lambdas/lambdaf+sqrt((1-lambdas/lambdaf)^2+4))/2
    # version developed for open volumetric receivers with open channel structure
    Ψ=1/ky*(poros-1)*(lambdas-lambdaf)/(lambdas*lambdaf)*(lambdaf*(poros-1)-lambdas*poros)
	λeff = (lambdas*lambdaf/(lambdas*poros+lambdaf*(1-poros)))*(poros*lambdas/lambdaf+(1-poros)+Ψ)*(poros+lambdaf/lambdas*(1-poros)+Ψ)/(poros*lambdas/lambdaf+2*Ψ+(1-poros)*lambdaf/lambdas+1)
end

Base.@kwdef struct SurfaceOpticalProps
	alpha_IR::Float64 = 0.2 # absorptance: obtain from datasheet/measurement
	tau_IR::Float64 = 0.5 # transmittance: obtain from datasheet/measurement
	rho_IR::Float64 = 1.0 - tau_IR - alpha_IR
	eps::Float64 = alpha_IR
	# effective opt. par. integrated over visible spectrum (vis)
	alpha_vis::Float64 = 0.0
	tau_vis::Float64 = 0.9 # obtain from datasheet, transmittance
	rho_vis::Float64 = 1.0 - tau_vis - alpha_vis 
end;

# optical parameters for quartz window in upper chamber (uc)
const uc_window = SurfaceOpticalProps(	
	alpha_IR=0.67,
	tau_IR=0.33,	
	alpha_vis=0.0,
	tau_vis=0.93,
)

#  optical parameters for catalyst layer in upper chamber (uc)
const uc_cat = SurfaceOpticalProps(
	alpha_IR=0.56, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.39, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# optical parameters for uncoated frit in upper chamber (uc)
const uc_mask = SurfaceOpticalProps(
	alpha_IR=0.55, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.15, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# optical parameters for uncoated frit in lower chamber (lc) = frit in upper chamber
const lc_frit = uc_mask

# optical parameters for Al bottom plate in lower chamber (lc)
const lc_plate = SurfaceOpticalProps(
	# aluminium bottom plate (anodized Al)
	alpha_IR=0.8, # https://www.osti.gov/servlets/purl/1483461
	tau_IR=0.0, # opaque surface
	alpha_vis=0.1, # see in lit, assume large reflectivity
	tau_vis=0.0 # opaque surface	
)

#
# ## Irradiation
# Irradiation flux coming from the solar simulator enters the aperture of the reactor with the specified flux profile as determined via photo-metric measurement. The measured profile is imported from a csv-datafile, handled by DataFrames.jl and interpolated by methodes provided by Interpolations.jl to be used as boundary condition in the simulation.
#

function readFlux(flux)

	path = pkgdir(MultiComponentReactiveMixtureProject)*"/data/IrradiationFluxProfiles/"
	fluxes = [40.0,60.0,80.0,100.0]
	@assert flux in fluxes
	d_flux_n = Dict(
	40.0=>("FluxMap_Import_20230523_110038_40.csv","Target_coords_20230523_110038_40.csv"),
	60.0=>	("FluxMap_Import_20230523_105203_60.csv","Target_coords_20230523_105203_60.csv"),
	80.0=>("FluxMap_Import_20230523_105025_80.csv","Target_coords_20230523_105025_80.csv"),
	100.0=>("FluxMap_Import_20230523_104820_100.csv","Target_coords_20230523_104820_100.csv")
	)
	fn_flux, fn_coord = d_flux_n[flux]
	FluxMap = CSV.read(path*fn_flux, DataFrame, header=false,delim=";")
	coords = CSV.read(path*fn_coord, DataFrame, header=1,delim=";")
	M=Matrix(FluxMap)
	reverse!(M; dims=1)
	return M, coords
end

function flux_interpol(flux)
	flux /= ufac"kW/m^2"
	nom_flux = 80.0
	if flux >= 30.0 && flux < 50.0
		nom_flux = 40.0
	elseif flux >= 50.0 && flux < 70.0
		nom_flux = 60.0
	elseif flux >= 70.0 && flux < 90.0
		nom_flux = 80.0
	elseif flux >= 90.0 && flux <= 100.0
		nom_flux = 100.0
	end
	M, coords = readFlux(nom_flux);
	M .*= flux/nom_flux*ufac"kW/m^2"
	

	# Interpolations.interpolate((coords.X*ufac"cm",coords.Y*ufac"cm"), M, Gridded(Linear()) )
	extrapolate(Interpolations.interpolate((coords.X*ufac"cm",coords.Y*ufac"cm"), M, Gridded(Linear()) ), Flat())
	#Interpolations.linear_interpolation((coords.X*ufac"cm",coords.Y*ufac"cm"), M,extrapolation_bc=Flat())	
	# return_intp(M)
end


"""
$(TYPEDEF)
Mutable data structure to hold modeling parameters of photo-thermal reactor
$(TYPEDFIELDS)
"""

##################################### BAK ########################################
# @kwdef mutable struct ReactorData{NG, KP}
#     "Spatial dimension, default=2"
# 	dim::Int64 = 2
#     # time constants for ramp functions
#     dt_mf::Tuple{Float64, Float64}=(0.0,1.0)
#     dt_hf_enth::Tuple{Float64, Float64}=(2.0,10.0)
#     dt_hf_irrad::Tuple{Float64, Float64}=(3.0,10.0)
    
#     # vectors holding boundary information
#     inlet_boundaries::Vector{Int64} = Int64[]
#     irradiated_boundaries::Vector{Int64} = Int64[]
#     outlet_boundaries::Vector{Int64} = Int64[]
#     side_boundaries::Vector{Int64} = Int64[]
#     catalyst_regions::Vector{Int64} =[2]
#     impermeable_regions::Vector{Int64} =[3]

#     kinpar::KP = XuFroment
#     ng::Int64 = kinpar.ng

#     # switches to control the simulation
# 	is_reactive::Bool = true
#     solve_T_equation::Bool = true
#     constant_properties::Bool = false
# 	include_Soret_Dufour::Bool = false
# 	include_dpdt::Bool = false

#     ip::Int64 = ng+1
# 	iT::Int64 = ip+1
# 	# inlcude window & plate temperatures as boundary species
# 	iTw::Int64=iT+1 # index of window Temperature (upper chamber)
# 	iTp::Int64=iTw+1 # index of plate Temperature (lower chamber)
	
# 	ibf::Int64 = dim == 3 ? iTp+1 : iTp # index of boundary flux species, workaround to include spatially varying boundary flux

# 	idpdt::Int64=ibf+1 # index of auxiliary variable holding dp/dt for use in enthalpy equation
	
# 	p::Float64 = 1.0*ufac"bar"
# 	Tamb::Float64 = 298.15*ufac"K"
# 	Treac::Float64 = 298.15*ufac"K"

# 	gn::Dict{Int, Symbol} = kinpar.gn # names and fluid indices
# 	gni::Dict{Symbol, Int} = kinpar.gni # inverse names and fluid indices
# 	Fluids::Vector{FluidProps} = kinpar.Fluids # fluids and respective properties in system
	
# 	m::Vector{Float64} = let
#         m=zeros(Float64, ng)
#         for i=1:ng
# 			m[i] = Fluids[i].MW
# 		end
# 		m
# 	end
	
# 	X0::Vector{Float64} = let
#         x=zeros(Float64, ng)
# 		x[gni[:H2]] = 1.0
# 		x[gni[:CO2]] = 1.0
# 		#x[gni[:N2]] = 1.0
# 		x/sum(x)
# 	end # inlet / initial composition
	
# 	mmix0::Float64 = sum(X0 .* m)
# 	W0::Vector{Float64} = @. m*X0/mmix0

# 	# Constant property storage
# 	constant_binary_diff_coeffs::Vector{Float64} = ones(Float64, ng*(ng-1)÷2) * 0.2*ufac"cm^2/s"
# 	constant_newman_soret_diff_coeffs::Vector{Float64} = ones(Float64, ng*(ng-1)÷2) * 0.1
# 	constant_species_viscosities::Vector{Float64} = ones(Float64, ng) * 2.0e-5*ufac"Pa*s"
# 	constant_species_thermal_conductivities::Vector{Float64} = ones(Float64, ng) * 0.02*ufac"W/(m*K)"

# 	# Inflow properties
# 	nflowin::Float64 = 7.4*ufac"mol/hr"
# 	mflowin::Float64 = nflowin*mmix0
# 	mfluxin::Float64 = mflowin/(100*ufac"cm^2")*ufac"kg/(m^2*s)"
# 	T_gas_in::Float64 = Tamb

#     # catalyst mass, volume and loading
# 	mcat::Float64=500.0*ufac"mg"
# 	Vcat::Float64=1.0*ufac"m^2"*0.02*ufac"cm"
# 	lcat::Float64=mcat/Vcat

# 	# upper and lower chamber heights for calculation of conduction/convection b.c.
# 	uc_h::Float64=17.0*ufac"mm"
# 	lc_h::Float64=18.0*ufac"mm"
# 	Nu::Float64=4.861	

# 	# optical properties of surfaces in IR/visible ranges
# 	uc_window::SurfaceOpticalProps = uc_window
# 	uc_cat::SurfaceOpticalProps = uc_cat
# 	uc_mask::SurfaceOpticalProps = uc_mask
# 	lc_frit::SurfaceOpticalProps = lc_frit
# 	lc_plate::SurfaceOpticalProps = lc_plate

# 	# reactor data
# 	delta_gap::Float64=1.5*ufac"mm" # gas gap between frit and reactor wall
# 	delta_wall::Float64=15*ufac"mm" # reactor wall thickness
# 	shell_h::Float64=20*ufac"mm" # height of heat transfer area adjancent to domain boundary 
# 	k_nat_conv::Float64=20.0*ufac"W/(m^2*K)" #17.5*ufac"W/(m^2*K)" # natural+forced convection heat transfer coeff.
# 	nom_flux::Float64 = 70.0*ufac"kW/m^2" # nominal irradiation flux density / kW/m^2

# 	# VitraPor data
# 	dp::Float64=200.0*ufac"μm" # average pore size, por class 0
# 	poros::Float64=0.33 # porosity, VitraPor sintetered filter class 0
# 	perm::Float64=1.23e-10*ufac"m^2" # perm. of porous medium, use in Darcy Eq.
# 	γ_τ::Float64=poros^1.5 # constriction/tourtuosity factor
	
# 	# Solid (non-porous) Borosilica glass (frit material)
# 	rhos::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
# 	lambdas::Float64=1.13*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
# 	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2

# 	# quartz window / Al bottom plate thermal conductivities
# 	lambda_window::Float64=1.38*ufac"W/(m*K)"
# 	lambda_Al::Float64=235.0*ufac"W/(m*K)"

# 	logreg::Float64=1.0e-20

#     function ReactorData(dim,dt_mf,dt_hf_enth,dt_hf_irrad,inlet_boundaries,irradiated_boundaries,outlet_boundaries,side_boundaries,catalyst_regions,impermeable_regions,kinpar,ng,is_reactive,solve_T_equation,constant_properties,include_Soret_Dufour,include_dpdt,ip,iT,iTw,iTp,ibf,idpdt,p,Tamb,Treac,gn,gni,Fluids,m,X0,mmix0,W0,constant_binary_diff_coeffs,constant_newman_soret_diff_coeffs,constant_species_viscosities,constant_species_thermal_conductivities,nflowin,mflowin,mfluxin,T_gas_in,mcat,Vcat,lcat,uc_h,lc_h,Nu,uc_window,uc_cat,uc_mask,lc_frit,lc_plate,delta_gap,delta_wall,shell_h,k_nat_conv,nom_flux,dp,poros,perm,γ_τ,rhos,lambdas,cs,lambda_window,lambda_Al,logreg)
#         KP = MultiComponentReactiveMixtureProject.KinData{nreac(kinpar)}
# 		# FluxIntp = flux_interpol(nom_flux)
# 		# flux_inner, flux_outer = flux_inner_outer(nom_flux)
#         new{ng,KP}(dim,dt_mf,dt_hf_enth,dt_hf_irrad,inlet_boundaries,irradiated_boundaries,outlet_boundaries,side_boundaries,catalyst_regions,impermeable_regions,kinpar,ng,is_reactive,solve_T_equation,constant_properties,include_Soret_Dufour,include_dpdt,ip,iT,iTw,iTp,ibf,idpdt,p,Tamb,Treac,gn,gni,Fluids,m,X0,mmix0,W0,constant_binary_diff_coeffs,constant_newman_soret_diff_coeffs,constant_species_viscosities,constant_species_thermal_conductivities,nflowin,mflowin,mfluxin,T_gas_in,mcat,Vcat,lcat,uc_h,lc_h,Nu,uc_window,uc_cat,uc_mask,lc_frit,lc_plate,delta_gap,delta_wall,shell_h,k_nat_conv,nom_flux,dp,poros,perm,γ_τ,rhos,lambdas,cs,lambda_window,lambda_Al,logreg)

#     end
# end
##################################### BAK ########################################


@kwdef mutable struct ReactorData
    "Spatial dimension, default=2"
	dim::Int64 = 2
    # time constants for ramp functions
    dt_mf::Tuple{Float64, Float64}=(0.0,1.0)
    dt_hf_enth::Tuple{Float64, Float64}=(2.0,10.0)
    dt_hf_irrad::Tuple{Float64, Float64}=(3.0,10.0)
    
    # vectors holding boundary information
    inlet_boundaries::Vector{Int64} = Int64[]
    irradiated_boundaries::Vector{Int64} = Int64[]
    outlet_boundaries::Vector{Int64} = Int64[]
    side_boundaries::Vector{Int64} = Int64[]
    catalyst_regions::Vector{Int64} =[2]
    impermeable_regions::Vector{Int64} =[3]

    # kinpar::KP = XuFroment
	# kinpar::KinData{3} = XuFroment
	kinpar::KinData = XuFroment
    ng::Int64 = kinpar.ng

	

    # switches to control the simulation
	is_reactive::Bool = true
    solve_T_equation::Bool = true
    constant_properties::Bool = false
	include_Soret_Dufour::Bool = false
	include_dpdt::Bool = false

    ip::Int64 = ng+1
	iT::Int64 = solve_T_equation ? ip+1 : 0
	# inlcude window & plate temperatures as boundary species
	iTw::Int64 = solve_T_equation ? iT+1 : 0 # index of window Temperature (upper chamber)
	iTp::Int64 = solve_T_equation ? iTw+1 : 0 # index of plate Temperature (lower chamber)
	
	ibf::Int64 = (solve_T_equation && dim == 3) ? iTp+1 : iTp # index of boundary flux species, workaround to include spatially varying boundary flux

	idpdt::Int64 = (solve_T_equation && include_dpdt) ? (dim == 3 ? ibf+1 : iTp+1) : 0 # index of auxiliary variable holding dp/dt for use in enthalpy equation
	
	# number of unknows / dofs / variables per node
	nvar::Int64 = maximum([ip, iT, iTw, iTp, ibf, idpdt])

	# pre-allocated buffers via PreallocationTools.jl
	X::DiffCache{Vector{Float64}, Vector{Float64}} = DiffCache(ones(nvar), 2 * nvar)	
	W::DiffCache{Vector{Float64}, Vector{Float64}} = DiffCache(ones(nvar), 2 * nvar)
    
	F::DiffCache{Vector{Float64}, Vector{Float64}} = DiffCache(ones(nvar), 2 * (nvar))

	D::DiffCache{Matrix{Float64}, Vector{Float64}} = DiffCache(ones(nvar, nvar), 2 * nvar)
	M::DiffCache{Matrix{Float64}, Vector{Float64}} = DiffCache(ones(nvar, nvar), 2 * (nvar))

	A::DiffCache{Matrix{Float64}, Vector{Float64}} = DiffCache(ones(nvar, nvar), 2 * (nvar))
	TDR::DiffCache{Vector{Float64}, Vector{Float64}} = DiffCache(ones(nvar), 2 * (nvar))

	dynvisc::DiffCache{Vector{Float64}, Vector{Float64}} = DiffCache(ones(nvar), 2 * nvar)
	thermcond::DiffCache{Vector{Float64}, Vector{Float64}} = DiffCache(ones(nvar), 2 * nvar)

	# vector holding pivoting information for linear system solution via RecursiveFactorizations.jl (see vfvm_functions.jl)
	ipiv::Vector{Int} = zeros(Int, ng-1)
	
	
	p::Float64 = 1.0*ufac"bar"
	Tamb::Float64 = 298.15*ufac"K"
	Treac::Float64 = 298.15*ufac"K"

	gn::Dict{Int, Symbol} = kinpar.gn # names and fluid indices
	gni::Dict{Symbol, Int} = kinpar.gni # inverse names and fluid indices
	Fluids::Vector{FluidProps} = kinpar.Fluids # fluids and respective properties in system
	
	m::Vector{Float64} = let
        m=zeros(Float64, ng)
        for i=1:ng
			m[i] = Fluids[i].MW
		end
		m
	end
	
	X0::Vector{Float64} = let
        x=zeros(Float64, ng)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		#x[gni[:N2]] = 1.0
		x/sum(x)
	end # inlet / initial composition
	
	mmix0::Float64 = sum(X0 .* m)
	W0::Vector{Float64} = @. m*X0/mmix0

	# Constant property storage
	constant_binary_diff_coeffs::Vector{Float64} = ones(Float64, ng*(ng-1)÷2) * 0.2*ufac"cm^2/s"
	constant_newman_soret_diff_coeffs::Vector{Float64} = ones(Float64, ng*(ng-1)÷2) * 0.1
	constant_species_viscosities::Vector{Float64} = ones(Float64, ng) * 2.0e-5*ufac"Pa*s"
	constant_species_thermal_conductivities::Vector{Float64} = ones(Float64, ng) * 0.02*ufac"W/(m*K)"

	# Inflow properties
	nflowin::Float64 = 7.4*ufac"mol/hr"
	mflowin::Float64 = nflowin*mmix0
	mfluxin::Float64 = mflowin/(100*ufac"cm^2")*ufac"kg/(m^2*s)"
	T_gas_in::Float64 = Tamb

    # catalyst mass, volume and loading
	mcat::Float64=500.0*ufac"mg"
	Vcat::Float64=1.0*ufac"m^2"*0.02*ufac"cm"
	lcat::Float64=mcat/Vcat

	# upper and lower chamber heights for calculation of conduction/convection b.c.
	uc_h::Float64=17.0*ufac"mm"
	lc_h::Float64=18.0*ufac"mm"
	Nu::Float64=4.861	

	# optical properties of surfaces in IR/visible ranges
	uc_window::SurfaceOpticalProps = uc_window
	uc_cat::SurfaceOpticalProps = uc_cat
	uc_mask::SurfaceOpticalProps = uc_mask
	lc_frit::SurfaceOpticalProps = lc_frit
	lc_plate::SurfaceOpticalProps = lc_plate

	# reactor data
	delta_gap::Float64=1.5*ufac"mm" # gas gap between frit and reactor wall
	delta_wall::Float64=15*ufac"mm" # reactor wall thickness
	shell_h::Float64=20*ufac"mm" # height of heat transfer area adjancent to domain boundary 
	k_nat_conv::Float64=20.0*ufac"W/(m^2*K)" #17.5*ufac"W/(m^2*K)" # natural+forced convection heat transfer coeff.
	nom_flux::Float64 = 70.0*ufac"kW/m^2" # nominal irradiation flux density / kW/m^2

	# VitraPor data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0
	poros::Float64=0.33 # porosity, VitraPor sintetered filter class 0
	perm::Float64=1.23e-10*ufac"m^2" # perm. of porous medium, use in Darcy Eq.
	γ_τ::Float64=poros^1.5 # constriction/tourtuosity factor
	
	# Solid (non-porous) Borosilica glass (frit material)
	rhos::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	lambdas::Float64=1.13*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2

	# quartz window / Al bottom plate thermal conductivities
	lambda_window::Float64=1.38*ufac"W/(m*K)"
	lambda_Al::Float64=235.0*ufac"W/(m*K)"

	logreg::Float64=1.0e-20

end


"""
Regularized logarithm:
```
   rlog(u,data)= log(u+electrolyte.logreg)
```
"""
rlog(x,data::ReactorData)=rlog(x,data.logreg)

function rlog(x,eps)
    if x<eps
        return log(eps)+(x-eps)/eps
    else
        return log(x)
    end
end

"""
    rexp(x;trunc=500.0)

Regularized exponential. Linear continuation for `x>trunc`,  
returns 1/rexp(-x) for `x<-trunc`.
"""
function rexp(x; trunc = 20.0)
    if x < -trunc
        1.0 / rexp(-x; trunc)
    elseif x <= trunc
        exp(x)
    else
        exp(trunc) * (x - trunc + 1)
    end
end

function kbed(data, lambdaf)
	(;poros,lambdas) = data
	B=1.25*((1.0-poros)/poros)^(10/9)
	kp=lambdas/lambdaf
	N=1-(B/kp)
	# !!! regularize function to catch poros=1 (pure gas phase) !!!
	#kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	kc=2/N* (B/N^2 * (kp-1)/kp * (log(kp)-rlog(B,data)) - (B+1)/2 - (B-1)/N)
	# !!! reglog
	1-sqrt(1-poros)+sqrt(1-poros)*kc
end

# ngas(::ReactorData{NG,KP}) where {NG,KP} = NG
ngas(data::ReactorData) = data.ng
nvar(data::ReactorData) = data.nvar

#function PTR_grid_boundaries_regions(dim;nref=0)
function PTR_grid_boundaries_regions(dim;uf=ufac"cm",nref=0,H=0.5,W=16)
	Ω_catalyst = 2
	# W=16
	# H=0.5

	nz = 10*2^(nref)
	Z=(0:(H/nz):H)*uf

	if dim == 2
		Γ_bottom = 1
		#Γ_bottom_insulating = 7
		Γ_side = 2
		Γ_sym = 4		
		Γ_top_permeable = 5
		Γ_top_irradiated = 6

		W=W/2 # axisymmetry, half domain is sufficient
		nr=W*2^(nref)
		R=(0:(W/nr):W)*uf
		
		grid=simplexgrid(R,Z)
		circular_symmetric!(grid)
	
		cellmask!(grid,[0,9/10*H].*uf,[W,H].*uf,Ω_catalyst) # catalyst layer	
		bfacemask!(grid, [0,H].*uf,[W-1,H].*uf,Γ_top_permeable; allow_new=false)
		bfacemask!(grid, [0,H].*uf,[W-2,0.5].*uf,Γ_top_irradiated; allow_new=false) 
				
		inb = [Γ_top_permeable,Γ_top_irradiated]
		irrb = [Γ_top_irradiated]
		outb = [Γ_bottom]
		sb = [Γ_side]
	else
		Γ_side_1 = 1 
		Γ_side_2 = 2
		Γ_side_3 = 3
		Γ_side_4 = 4		
		Γ_bottom = 5
		Γ_top_permeable = 7
		Γ_top_irradiated = 8
		
		nxy=W*2^(nref)
		X=(0:(W/nxy):W)*uf
		Y=X
		# X=(0:1:W)*uf
		# Y=(0:1:W)*uf
		# Z=(0:H/10:H)*uf	
		grid=simplexgrid(X,Y,Z)


	
		# catalyst region
		cellmask!(grid,[0,0,9/10*H].*uf,[W,W,H].*uf,Ω_catalyst) # catalyst layer

		# mark gas permeable region in inlet plane
		# bfacemask!(grid, [1,1,H].*uf,[W-1,W-1,H].*uf,Γ_top_permeable; allow_new=false)
		w_top_perm = (W/2)/8
		bfacemask!(grid, [w_top_perm,w_top_perm,H].*uf,[W-w_top_perm,W-w_top_perm,H].*uf,Γ_top_permeable; allow_new=false)

		# mark irradiated region in inlet plane
		w_top_irrad = (W/2)/4
		# bfacemask!(grid, [2,2,H].*uf,[W-2,W-2,H].*uf,Γ_top_irradiated; allow_new=false)
		bfacemask!(grid, [w_top_irrad,w_top_irrad,H].*uf,[W-w_top_irrad,W-w_top_irrad,H].*uf,Γ_top_irradiated; allow_new=false)

		inb = [Γ_top_permeable,Γ_top_irradiated]
		irrb = [Γ_top_irradiated]
		outb = [Γ_bottom]
		sb = [Γ_side_1,Γ_side_2,Γ_side_3,Γ_side_4]
	end

	return grid, inb, irrb, outb, sb, [Ω_catalyst]
end

function PTR_init_system(dim, grid, data)

	(;p,ip,Tamb,iT,iTw,iTp,ibf,idpdt,inlet_boundaries,irradiated_boundaries,outlet_boundaries,catalyst_regions,X0,solve_T_equation,include_dpdt,nom_flux) = data
	ng=ngas(data)

	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=DMS_flux,
							reaction=DMS_reaction,
							storage=DMS_storage,
							bcondition=PTR_bcond,
							bflux=PTR_bflux,
							bstorage=PTR_bstorage,
							boutflow=DMS_boutflow,
							outflowboundaries=outlet_boundaries,
							assembly=:edgewise,
							unknown_storage=:dense
							)

	enable_species!(sys; species=collect(1:ng)) # gas phase species xi
	enable_species!(sys; species=ip) # gas phase pressure
	if solve_T_equation
		enable_species!(sys; species=iT) # temperature T
		enable_boundary_species!(sys, iTw, irradiated_boundaries) # window temperature as boundary species in upper chamber
		enable_boundary_species!(sys, iTp, outlet_boundaries) # plate temperature as boundary species in lower chamber

		# for 3 dimensional domain, apply measured irradiation flux density as boundary condition
		if dim == 3
			# boundary flux species, workaround to implement spatially varying irradiation
			enable_boundary_species!(sys, ibf, irradiated_boundaries)
		end
		if include_dpdt
			enable_species!(sys; species=idpdt) # auxiliary variable holding dpdt term
		end		
	end

	inival=unknowns(sys)
	inival[ip,:].=p

	for i=1:ng
		inival[i,:] .= X0[i]
	end

	if solve_T_equation
		inival[[iT,iTw,iTp],:] .= Tamb
		if dim == 3

            FluxIntp = flux_interpol(nom_flux)
			function d3tod2(a,b)
				a[1]=b[1]
				a[2]=b[2]
			end
			inival[ibf,:] .= zero(eltype(inival))
			# sub=ExtendableGrids.subgrid(grid,irradiated_boundaries,boundary=true, transform=d3tod2 )
			subg_window = ExtendableGrids.subgrid(grid, irradiated_boundaries, boundary=true, transform=d3tod2)
			W_window = maximum(subg_window[Coordinates][1,:]) - minimum(subg_window[Coordinates][1,:]) # square window
			W_domain = maximum(grid[Coordinates][1,:]) - minimum(grid[Coordinates][1,:]) # square prismatic domain
			
			coord_offset = 0.5*(W_domain-W_window)
				
			for inode in subg_window[CellNodes]
				c = subg_window[Coordinates][:,inode]
				inodeip = subg_window[ExtendableGrids.NodeParents][inode]
				# inival[ibf,inodeip] = FluxIntp(c[1]-0.02, c[2]-0.02)
				inival[ibf,inodeip] = FluxIntp(c[1]-coord_offset, c[2]-coord_offset)
			end
		end
		if include_dpdt
			inival[idpdt,:] .= zero(eltype(inival))
		end
	end

	if dim > 1
		catalyst_nodes = []
		for reg in catalyst_regions
			catalyst_nodes = vcat(catalyst_nodes, unique(grid[CellNodes][:,grid[CellRegions] .== reg]) )
		end
			
		cat_vol = sum(nodevolumes(sys)[unique(catalyst_nodes)])

		data.lcat = data.mcat/cat_vol
		local Ain = 0.0
		for boundary in inlet_boundaries
			Ain += bareas(boundary,sys,grid)
		end
		data.mfluxin = data.mflowin / Ain
	end
	
	return inival,sys
end




# Function to calculate surface areas of boundaries. Used for calculation of 
# flows / fluxes through boundaries
function bareas(bfaceregion,sys,grid)
	area = 0.0
	for ibface =1:VoronoiFVM.num_bfaces(grid)
        for inode=1:VoronoiFVM.num_nodes(grid[VoronoiFVM.ExtendableGrids.BFaceGeometries][1])
			if grid[VoronoiFVM.ExtendableGrids.BFaceRegions][ibface] == bfaceregion
                area +=bfacenodefactors(sys)[inode,ibface]
			end
		end
	end
	area
end

# calculate non-dimensional numbers: Reynolds, Prandtl, Peclet
function RePrPe(data,T,p,x)
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
