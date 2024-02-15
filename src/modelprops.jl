abstract type AbstractModelData end

@doc raw"""
# The model
Model equations for multi-component transport of reacting (ideal) gas phase mixtures in porous media.
## Overall Mass Continuity
Mixture mass flow (overall mass flow) through the pore space of the porous medium.
Mixture mass averaged velocity is calculated from Darcy equation. The void fraction (porosity) is given by $\epsilon$.

```math
\begin{align}
	\frac{\partial \epsilon \rho}{\partial t} + \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\end{align}
```

## Species Mass Continuity and Transport
```math
\begin{align}
	\frac{\partial \epsilon \rho_i}{\partial t} + \nabla \cdot \left( \vec \Phi_i + \rho_i \vec v \right ) - r_i(\varrho) &= 0 ~, \qquad i = 1 ... n \\
		\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= -\sum_{j=1 \atop j \neq i}^{n} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
		\sum_{i=1}^n x_i &= 1
\end{align}
```

where $\rho$ is the (total) mixture density, $\vec v$ is the mass-averaged (barycentric)  mixture velocity calculated with the Darcy equation,
$x_i$, $w_i$ and $M_i$ are the molar fraction, mass fraction and molar mass of species $i$ respectively,
$\vec \Phi_i$ is the __diffusive__ mass flux of species $i$ ($\frac{\text{kg}}{\text{m}^2 \text{s}}$)
and $r_i$ is the species mass volumetric source/sink ($\frac{\text{kg}}{\text{m}^3 \text{s}}$) of gas phase species $i$.

The __convective__ mass flux of species $i$ is the product of the superficial mean velocity with the partial mass density $\rho_i \vec v$.
The combined __convective-diffusive__ species mass flux $\vec \Psi_i = \vec \Phi_i + \rho_i \vec v$ can be introduced as an auxiliary variable.


"""
function DMS_Info_isothermal() end

@doc raw"""
## Thermal Energy Transport

The thermal energy equation considers convective-diffusive transport of thermal energy and is formulated with effective transport parameters that are a consequence of the quasi-homogeneous phase approach:

```math
\begin{align}
\frac{\partial \overline h}{\partial t} - \nabla \cdot \left(\lambda_{\text{eff}} \nabla T - \sum_i^n \vec \Psi_i h_i / M_i \right) &= 0 \\
\overline h &= (\varepsilon \sum_i^n \rho_i h_i / M_i + [1-\varepsilon] \rho_{\text s} h_{\text s}) \\
h_i(T) &= h_i^0(T_{\text{ref}}) + \int_{T_{\text{ref}}}^{T} c_{p,i} \,dT 
\end{align}
```
where $\overline h$ is the volumetric enthalpy [$\text{J/m}^3$] of the quasi-homogeneous domain comprised of gas and solid volume fractions, the gas phase species enthalpy $h_i$ and $\rho_{\text s}$ and $h_{\text s}$ are the mass density and enthalpy of the solid matrix respectively.

Here $\lambda_{\text{eff}}$ is the effective thermal conductivity in the quasi-homogeneous approach, which is a (complicated) function of the thermal conductivities of the solid and gas phases.

Another process contributing to heat transport within the porous medium is "thermal drift", resulting from gas phase mass flow carrying enthalpy.
It is expressed as the sum  of species mass fluxes and species specific enthalpies $\sum_i^n \vec \Psi_i h_i / M_i$ (division by $M_i$ as $h_i$ is given on a molar basis).
The gas phase species enthalpies consist of a chemical contribution originating from the potential energy stored within the chemical bonds of the species, and a thermal contribution depending on the species' heat capacity and temperature.
During chemical reactions the chemical bonds in the species involved in the reaction are rearranged. If in the process of rearranging of chemical bonds the products have a different chemical potential energy than the reactants, the difference is made up from thermal energy.
In case of an exothermal reaction, the products have less chemical energy thus the difference is released as heat. In contrast for endothermal reactions the products hava more chemical energy than the reactants and thus consume heat during the reaction.
"""
function DMS_Info_thermal() end

@doc raw"""
Assemble symmetric Maxwell-Stefan matrix M expressing the flux-force relations in the Maxwell-Stefan
framework for multicomponent diffusion.
"""
function M_matrix!(M, W, D, data)
	ng=ngas(data)
	@inbounds for i=1:(ng-1)
		for j=1:(ng-1)
			M[i,j] = zero(eltype(W))
		end
		 for j=1:ng
			if i != j
				M[i,i] -= W[j]/D[i,j]
				if j == ng
					 for k=1:(ng-1)
						M[i,k] -= W[i]/D[i,ng]
					end
				else
					M[i,j] += W[i]/D[i,j]
				end
			end
		end
	end			
end
@doc raw"""
Assemble symmetric Maxwell-Stefan diffusivity matrix D. Account for porous material
through the constriction and tourtuosity factor γ_τ which lowers the diffusivities
 ~1 order of magnitude compared to diffusion through free space.
"""
function D_matrix!(D, T, p, data)
	(;m,γ_τ,constant_properties)=data
	ng=ngas(data)
	@inbounds for i=1:(ng-1)
		for j=(i+1):ng
            if !constant_properties
			    Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
			    Dji *= m[i]*m[j]*γ_τ # eff. diffusivity
            else
                Dji = 2.0e-5*m[i]*m[j]
            end
			D[j,i] = Dji
			D[i,j] = Dji
		end
	end
end


@doc raw"""
Mixture mass flow (bulk convective mass flow) through the pore space of the porous medium. Mixture mass averaged (barycentric) velocity for convective bulk mass flow through the porous medium is calculated from Darcy equation. 

```math
	\vec v  = -\frac{\kappa}{\mu} \vec \nabla p
```
"""

function DarcyVelo(u,data,mu)
	(;ip,perm) = data
	-perm/mu*(u[ip,1]-u[ip,2])	
end

function MoleFrac!(X,u::VoronoiFVM.EdgeUnknowns,data)
	@inbounds for i=1:ngas(data)
		X[i] = 0.5*(u[i,1]+u[i,2])
	end
	nothing
end

function MoleFrac!(X,u::VoronoiFVM.BNodeUnknowns,data)
	@inbounds for i=1:ngas(data)
		X[i] = u[i]
	end
	nothing
end

function MoleFrac!(X,u::VoronoiFVM.NodeUnknowns,data)
	@inbounds for i=1:ngas(data)
		X[i] = u[i]
	end
	nothing
end

function MassFrac!(X,W,data)
	(;m) = data
	@inline mmix = molarweight_mix(X,data)
	@inbounds for i=1:ngas(data)
		W[i] = X[i]*m[i]/mmix
	end
	nothing
end

@doc raw"""
Flux function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_flux(f,u,edge,data)
	(;m,ip,iT,dt_hf_enth,solve_T_equation,Tamb,Fluids)=data
	ng=ngas(data)
		
	F = MVector{ng-1,eltype(u)}(undef)
	X = MVector{ng,eltype(u)}(undef)
	W = MVector{ng,eltype(u)}(undef)
	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)
	D = MMatrix{ng,ng,eltype(u)}(undef)

	pm = 0.5*(u[ip,1]+u[ip,2])
    Tm = solve_T_equation ? 0.5*(u[iT,1]+u[iT,2]) : one(eltype(u))*Tamb
	c = pm/(ph"R"*Tm)
	
	δp = u[ip,1]-u[ip,2]
	
	@inline MoleFrac!(X,u,data)
	@inline mmix = molarweight_mix(X,data)
	@inline MassFrac!(X,W,data)
	
	
	@inline D_matrix!(D, Tm, pm, data)
	@inline mumix, lambdamix = dynvisc_thermcond_mix(data, Tm, X)
		
	rho = c*mmix
	v = DarcyVelo(u,data,mumix)
	
	f[ip] = -rho*v

	@inline M_matrix!(M, W, D, data)
	
	@inbounds for i=1:(ng-1)
		F[i] = ( u[i,1]-u[i,2] + (X[i]-W[i])*δp/pm )*c/mmix
	end				

	@inline inplace_linsolve!(M,F)

	enthalpy_flux = zero(eltype(u))
	mass_flux = zero(eltype(u))
	@inbounds for i=1:(ng-1)
		f[i] = -(F[i] + c*X[i]*m[i]*v)
		if solve_T_equation
			mass_flux += f[i]
			# !!! DEBUG !!!
			# enthalpy_flux += f[i] * enthalpy_gas(Fluids[i], Tm) / m[i]
			# enthalpy_flux += f[i] * heatcap_gas(Fluids[i], Tm) / m[i] * (Tm - Tref)
			enthalpy_flux += f[i] * enthalpy_gas_thermal(Fluids[i], Tm) / m[i]
			# !!! DEBUG !!!
		end
	end	
	
    if solve_T_equation
		# !!! DEBUG !!!
		# enthalpy_flux += (f[ip] - mass_flux) * enthalpy_gas(Fluids[ng], Tm) / m[ng]# species n
		# enthalpy_flux += (f[ip] - mass_flux) * heatcap_gas(Fluids[ng], Tm) / m[ng] * (Tm - Tref)# species n
		enthalpy_flux += (f[ip] - mass_flux) * enthalpy_gas_thermal(Fluids[ng], Tm) / m[ng] # species n
		# !!! DEBUG !!!
        lambda_bed=kbed(data,lambdamix)*lambdamix
        # @inline hf_conv = f[ip] * enthalpy_mix(data, Tm, X) / mmix * ramp(edge.time; du=(0.0,1), dt=dt_hf_enth) 
        # Bp,Bm = fbernoulli_pm(hf_conv/lambda_bed/Tm)
        # f[iT] = lambda_bed*(Bm*u[iT,1]-Bp*u[iT,2])

		# !!! DEBUG !!!
		f[iT] = lambda_bed*(u[iT,1]-u[iT,2]) + enthalpy_flux * ramp(edge.time; du=(0.0,1), dt=dt_hf_enth) 
		# f[iT] = lambda_bed*(u[iT,1]-u[iT,2]) + enthalpy_flux
		# f[iT] = lambda_bed*(u[iT,1]-u[iT,2])
		# !!! DEBUG !!!
    end
end

@doc raw"""
Reaction function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_reaction(f,u,node,data)
	(;m,ip,iT,is_reactive,catalyst_regions)=data
	ng=ngas(data)

	if node.region in catalyst_regions && is_reactive # catalyst layer
		(;lcat,kinpar,iT,Treac,solve_T_equation)=data
		(;nuij)=kinpar
		
		pi = MVector{ng,eltype(u)}(undef)
		for i=1:ng
            pi[i] = u[ip]*u[i]
		end

        T = solve_T_equation ? u[iT] : one(eltype(u))*Treac
        RR = @inline -lcat*ri(data,T,pi)
        #RR = @inline -lcat*ri(data,u[iT],pi)
		
		for i=1:ng
			f[i] = zero(eltype(u))
			for j=1:nreac(kinpar)
				f[i] += nuij[(j-1)*ng+i] * RR[j] * m[i]
				if solve_T_equation
					f[iT] -= nuij[(j-1)*ng+i] * RR[j] * enthalpy_gas(kinpar.Fluids[i], T)
				end
			end			
		end
	end

	for i=1:ng
		f[ng] += u[i]
	end
	f[ng] = f[ng] - 1.0
end

@doc raw"""
Storage function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_storage(f,u,node,data)
	(;ip,iT,Tamb,m,poros,rhos,cs,solve_T_equation)=data
	ng=ngas(data)

    T = solve_T_equation ? u[iT] : one(eltype(u))*Tamb
    c = u[ip]/(ph"R"*T)
    mmix = zero(eltype(u))
	for i=1:ng
		f[i]=c*u[i]*m[i]*poros
        mmix += u[i]*m[i]
	end
	
	# total pressure
	f[ip] = mmix*c*poros

    if solve_T_equation
        X=MVector{ng,eltype(u)}(undef)
        @inline MoleFrac!(X,u,data)
        @inline cpmix = heatcap_mix(data, T, X)
		# !!! DEBUG !!!
		f[iT] = u[iT] * (rhos*cs*(1-poros) + cpmix*c*poros)
		# f[iT] = (u[iT]-Tref) * (rhos*cs*(1-poros) + cpmix*c*poros)
		# !!! DEBUG !!!
    end
	
end

@doc raw"""
Outflow boundary function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_boutflow(f,u,edge,data)
	(;iT,ip,m,dt_hf_enth,Tamb,solve_T_equation,Fluids)=data
	ng=ngas(data)

	k=outflownode(edge)
	pout = u[ip,k]
    Tout = solve_T_equation ? u[iT,k] : one(eltype(u))*Tamb
    cout = pout/(ph"R"*Tout)

	X = MVector{ng,eltype(u)}(undef)
	
	for i=1:ng
		X[i] = u[i,k]
	end
	
    @inline mumix, _ = dynvisc_thermcond_mix(data, Tout, X)
	v = DarcyVelo(u,data,mumix)
	
	for i=1:(ng-1)
		f[i] = v*cout*u[i,k]*m[i] # species mass flux at outflow
	end

    if solve_T_equation
        # @inline r_hf_enth = v *cout * enthalpy_mix(data, Tout, X) * ramp(edge.time; du=(0.0,1.0), dt=dt_hf_enth) # enthalpy heat flux
		# !!! DEBUG !!!
		
		hout = zero(eltype(u))
        @inbounds for i=1:ng
            hout += X[i] * enthalpy_gas_thermal(Fluids[i], Tout)
        end

		#@inline hout = heatcap_mix(data, Tout, X) * (Tout - Tref)
		r_hf_enth = v * cout * hout
        # f[iT] = r_hf_enth
		f[iT] = r_hf_enth * ramp(edge.time; du=(0.0,1.0), dt=dt_hf_enth)
		# f[iT] = zero(eltype(u))
		# !!! DEBUG !!!
    end
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

#"""
#Effective thermal conductivity of porous filter frit according to:
#
#__Zehner, P., & Schlünder, E. U. (1970).__ Wärmeleitfähigkeit von Schüttungen bei mäßigen Temperaturen. Chemie Ingenieur Technik, 42(14), 933-941. doi:10.1002/cite.330421408
#
#Implementation follows the notation in __VDI Heat Atlas 2010, ch. D6.3 eqs. (5a-5e).__
#"""
# calculate fixed bed (porous material) effective thermal conductivity


function kbed(data, lambdaf)
	(;poros,lambdas) = data
	B=1.25*((1.0-poros)/poros)^(10.0/9.0)
	kp=lambdas/lambdaf
	N=1.0-(B/kp)
	kc=2.0/N* (B/N^2.0*(kp-1.0)/kp*log(kp/B) - (B+1.0)/2.0 - (B-1.0)/N)
	1.0-sqrt(1.0-poros)+sqrt(1.0-poros)*kc
end

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

	path = pkgdir(FixedBed)*"/data/IrradiationFluxProfiles/"
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

# function flux_inner_outer(flux)
# 	flux /= ufac"kW/m^2"
# 	nom_flux = 0.0
# 	if flux >= 30.0 && flux < 50.0
# 		nom_flux = 40.0
# 	elseif flux >= 50.0 && flux < 70.0
# 		nom_flux = 60.0
# 	elseif flux >= 70.0 && flux < 90.0
# 		nom_flux = 80.0
# 	elseif flux >= 90.0 && flux <= 100.0
# 		nom_flux = 100.0
# 	end
# 	M, coords = readFlux(nom_flux);
# 	M .*= flux/nom_flux*ufac"kW/m^2"
# 	sc, ec, sr, er = 32, 32, 33, 31
# 	width10 =  coords.X[end-ec] - coords.X[sc]
# 	height10 =  coords.Y[end-er] - coords.X[sr]
# 	M10 = M[sc:end-ec,sr:end-er]
# 	pwr_10 = sum(M10)*width10*height10*ufac"cm"^2/length(M10)
# 	flux_inner = sum(M10)/length(M10)

# 	width12 =  coords.X[end] - coords.X[1]
# 	height12 =  coords.Y[end] - coords.X[1]
# 	pwr_12 =  sum(M)*width12*height12*ufac"cm"^2/length(M)

# 	area_12 = (width12*height12-width10*height10)*ufac"cm"^2
# 	flux_outer = (pwr_12-pwr_10)/area_12

# 	#extrapolate(Interpolations.interpolate((coords.X*ufac"cm",coords.Y*ufac"cm"), M, Gridded(Linear()) ), Line())
# 	#Interpolations.linear_interpolation((coords.X*ufac"cm",coords.Y*ufac"cm"), M,extrapolation_bc=Flat())	
# 	#return_intp(M)
# 	return (flux_inner, flux_outer)
# end

# function return_intp(M;wi=16.0*ufac"cm",wi_apt=12.0*ufac"cm")
		
# 	# starting coordinates for optimum 10 cm x 10 cm selection
# 	sr=33
# 	sc=33

# 	# (inverse) resolution of flux measurements, distance between data points
# 	Dx = 0.031497*ufac"cm"
# 	Dy = 0.031497*ufac"cm"
	
	
# 	# calculate coordinate offsets, when deviating (expanding/shrinking from the center) from 10cm x 10cm selection 	
# 	Dsr=Integer(round((wi_apt-10.0*ufac"cm")/(2*Dy)))
# 	Dsc=Integer(round((wi_apt-10.0*ufac"cm")/(2*Dx)))

# 	sr-=Dsr
# 	sc-=Dsc

# 	nx=Integer(round(wi_apt/Dx))
# 	ny=Integer(round(wi_apt/Dy))
	
# 	@views M_ = M[sr:(sr+ny-1),sc:(sc+nx-1)]

# 	# pad Flux Map with zeros outside of aperture area
# 	nx_dom = Integer(round(wi/Dx))
# 	ny_dom = Integer(round(wi/Dy))
# 	M__ = zeros(nx_dom,ny_dom)
		
# 	Dsr_dom_apt=Integer(round((wi - wi_apt)/(2*Dy)))
# 	Dsc_dom_apt=Integer(round((wi - wi_apt)/(2*Dx)))

# 	M__[(Dsr_dom_apt+1):(Dsr_dom_apt+ny), (Dsc_dom_apt+1):(Dsc_dom_apt+nx)] = M[sr:(sr+ny-1),sc:(sc+nx-1)]

# 	x = range(0.0,wi,length=nx_dom)
# 	y = range(0.0,wi,length=nx_dom)
	
# 	Interpolations.interpolate((x,y), M__, Gridded(Linear()))
# end

const itp12 = flux_interpol(70.0*ufac"kW/m^2")

@kwdef mutable struct ReactorData{NG, KP}
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

    kinpar::KP = XuFroment
    ng::Int64 = kinpar.ng

    # switches to control the simulation
	is_reactive::Bool = true
    solve_T_equation::Bool = true
    constant_properties::Bool = false

	#ip::Int64 = NG+1
    ip::Int64 = ng+1
	iT::Int64 = ip+1
	# inlcude window & plate temperatures as boundary species
	iTw::Int64=iT+1 # index of window Temperature (upper chamber)
	iTp::Int64=iTw+1 # index of plate Temperature (lower chamber)
	
	ibf::Int64=iTp+1 # index of boundary flux species, workaround to include spatially varying boundary flux
	
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
	end # inlet composition
	
	mmix0::Float64 = sum(X0 .* m)
	W0::Vector{Float64} = @. m*X0/mmix0

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
	FluxIntp::typeof(itp12)=itp12
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

    function ReactorData(dim,dt_mf,dt_hf_enth,dt_hf_irrad,inlet_boundaries,irradiated_boundaries,outlet_boundaries,side_boundaries,catalyst_regions,impermeable_regions,kinpar,ng,is_reactive,solve_T_equation,constant_properties,ip,iT,iTw,iTp,ibf,p,Tamb,Treac,gn,gni,Fluids,m,X0,mmix0,W0,nflowin,mflowin,mfluxin,T_gas_in,mcat,Vcat,lcat,uc_h,lc_h,Nu,uc_window,uc_cat,uc_mask,lc_frit,lc_plate,delta_gap,delta_wall,shell_h,k_nat_conv,nom_flux,FluxIntp,dp,poros,perm,γ_τ,rhos,lambdas,cs,lambda_window,lambda_Al)
        KP = FixedBed.KinData{nreac(kinpar)}
		FluxIntp = flux_interpol(nom_flux)
		# flux_inner, flux_outer = flux_inner_outer(nom_flux)
        new{ng,KP}(dim,dt_mf,dt_hf_enth,dt_hf_irrad,inlet_boundaries,irradiated_boundaries,outlet_boundaries,side_boundaries,catalyst_regions,impermeable_regions,kinpar,ng,is_reactive,solve_T_equation,constant_properties,ip,iT,iTw,iTp,ibf,p,Tamb,Treac,gn,gni,Fluids,m,X0,mmix0,W0,nflowin,mflowin,mfluxin,T_gas_in,mcat,Vcat,lcat,uc_h,lc_h,Nu,uc_window,uc_cat,uc_mask,lc_frit,lc_plate,delta_gap,delta_wall,shell_h,k_nat_conv,nom_flux,FluxIntp,dp,poros,perm,γ_τ,rhos,lambdas,cs,lambda_window,lambda_Al)

    end
end

ngas(::ReactorData{NG,KP}) where {NG,KP} = NG

@doc raw"""
Helper function to calculate flux integrals over in- and outflow boundaries in 
    the Darcy-Maxwell-Stefan (DMS) model.
"""
function DMS_checkinout(sol,sys,data)	
	(;inlet_boundaries,outlet_boundaries)=data
	tfact=TestFunctionFactory(sys)

	tf_in=testfunction(tfact,outlet_boundaries,inlet_boundaries)
	tf_out=testfunction(tfact,inlet_boundaries,outlet_boundaries)
	
	(;in=integrate(sys,tf_in,sol),out=integrate(sys,tf_out,sol) )
end

@doc raw"""
Helper function to print a summary based on calculated flux integrals over in- 
    and outflow boundaries in the Darcy-Maxwell-Stefan (DMS) model.
"""
function DMS_print_summary(sol,grid,sys,data)
    (;ip,m,mfluxin,mmix0,X0,ng,inlet_boundaries) = data
    in_,out_=DMS_checkinout(sol,sys,data)

    nout(i) = out_[i]/m[i]
    local Ain = 0.0
	for boundary in inlet_boundaries
		Ain += bareas(boundary,sys,grid)
	end
    nin(i) = mfluxin/mmix0*Ain*X0[i]

    RI=sum(integrate(sys,sys.physics.reaction,sol),dims=2) # reaction integral

    println("Total mass inflows and outflows:")
    @printf "IN: %2.6e \t OUT: %2.6e \t REACT: %2.6e kg/hr \nSUM: %2.6e kg/hr\n\n" in_[ip]/ufac"kg/hr" out_[ip]/ufac"kg/hr" RI[ip]/ufac"kg/hr" (in_[ip]+out_[ip]+RI[ip])/ufac"kg/hr"

    println("Molar species inflows, outflows and reaction integrals:")
    for i = 1:ng
        @printf "%i\tIN: %2.6e \t OUT: %2.6e \t REACT: %2.6e mol/hr \n\tSUM: %2.6e mol/hr\n" i nin(i)/ufac"mol/hr" nout(i)/ufac"mol/hr" -RI[i]/m[i]/ufac"mol/hr" (nin(i)+nout(i)-RI[i]/m[i])/ufac"mol/hr"
    end
end
@doc raw"""
Helper function to print an extended summary based on calculated flux integrals over in- 
    and outflow boundaries in the Darcy-Maxwell-Stefan (DMS) model.
"""
function DMS_print_summary_ext(solt,grid,sys,data)

	t = solt.t[end]
	sol = solt(t)

	(;gn,gni,m,nflowin,X0) = data
	ng=ngas(data)
	in_,out_=FixedBed.DMS_checkinout(sol,sys,data)

	nout(i) = -out_[i]/m[i]
	nin(i) = nflowin*X0[i]
	nout_dry = 0.0
	
	println("Molar species in- & outflows:")
	for i = 1:ng
		@printf "%s\tIN: %2.2f\t OUT: %2.2f mol/hr\n" gn[i] nin(i)/ufac"mol/hr" nout(i)/ufac"mol/hr"
		if i != gni[:H2O] 
			nout_dry += nout(i)
		end
	end


	println("\nDry Product Molar Fractions:")
	for i=1:ng
		if i != gni[:H2O] && i != gni[:N2] 
		@printf "%3s: %2.1f%%\n" gn[i] nout(i)/nout_dry*100
		end
	end
	
	println("\nConversion:")
	for i in(gni[:CO2],gni[:H2])
		@printf "X%4s: %2.2f\n" gn[i] (nin(i)-nout(i))/nin(i)
	end
	println("\nYield & Selectivity (CO2 based):")
	for i in(gni[:CO],gni[:CH4])
		@printf "Y%4s: %2.2f \tS%4s: %2.2f\n" gn[i] nout(i)/nin(gni[:CO2]) gn[i] nout(i)/(nin(gni[:CO2])-nout(gni[:CO2]))
	end

	fluxes = HeatFluxes_EB_I(t,solt,grid,sys,data)
	ns = keys(fluxes)
	vs = values(fluxes)
	#QG_01,	QG_10, H_thermal, H_reaction, Qconv_10,	Qconv_34, QG_34, QG_43,	Qsides = fluxes
	println("\nEnergy Balancing [W]:")
	for i=1:length(fluxes)
		@printf "%s:\t%.2f \n" String(ns[i]) vs[i]
	end
	@printf "Sum:\t\t%.6f \n" sum(fluxes)

end

@doc raw"""
Helper function to export to VTK format for visualization 3D solutions of the
    photo thermal catalytic reactor (PCR) model. Exported are the pressure, 
    species molar fractions and temperature fields in the 3D domain.
"""
function PCR_writeSol3D(sol,grid,data;desc="")
    (;dim,ip,iT,ibf,gn,nflowin,solve_T_equation) = data
    ng=ngas(data)
    _t = now()
    tm = "$(hour(_t))_$(minute(_t))_$(second(_t))"
    desc = isempty(desc) ? desc : "_"*desc
    path = "../data/out/$(Date(_t))/$(tm)$(desc)"
    try
        mkpath(path)
    catch e
        println("Directory " * path * " already exists.")
    end
    #mkdir(string(Date(_t)))
    
    VoronoiFVM.writeVTK("$(path)/$(tm)_3D_ptot_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").vtu", grid; point_data = sol[ip,:])
	if solve_T_equation
        VoronoiFVM.writeVTK("$(path)/$(tm)_3D_T_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").vtu", grid; point_data = sol[iT,:] .-273.15)
    end
    for i=1:ng
        VoronoiFVM.writeVTK("$(path)/$(tm)_3D_x$(gn[i])_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").vtu", grid; point_data = sol[i,:])
    end
	if dim == 3
		VoronoiFVM.writeVTK("$(path)/$(tm)_3D_irradiation_flux_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").vtu", grid; point_data = sol[ibf,:])
	end
end

@doc raw"""
Helper function to calculate radiation emitted from the window towards the surface
    of the catalyst layer in the upper chamber in the photo thermal catalytic
    reactor (PCR) model.
"""
function radiosity_window(f,u,bnode,data)
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
		G_lamp += nom_flux
		# G_lamp += u[ibf]
		# !!! DEBUG !!!
	end
	

    Tglass = u[iTw] # local tempererature of quartz window
    G1_bot_IR = eps1*ph"σ"*(Tglass^4-Tamb^4)+ph"σ"*Tamb^4
	#G1_bot_IR = ph"σ"*Tamb^4
	G1_bot_vis = zero(eltype(u))
    if bnode.region in irradiated_boundaries		
		G1_bot_vis += G_lamp # flux profile measured behind quarz in plane of cat layer
    end
    return G1_bot_vis,G1_bot_IR
end

@doc raw"""
Function defining the top/inlet boundary condition in the photo thermal catalytic
    reactor (PCR) model.
"""
function PCR_top(f,u,bnode,data)
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

#	G1_vis, G1_IR = radiosity_window(f,u,bnode,data)
	
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
			@inbounds for i=1:ng
				hin += enthalpy_gas_thermal(Fluids[i],T_gas_in)
			end
			# @inline hin = heatcap_mix(data, T_gas_in, X0) * (T_gas_in - Tref)
			r_hf_enth = mfluxin/mmix0 * hin
			f[iT] = -r_hf_enth * ramp(bnode.time; du=(0.0,1), dt=dt_hf_enth)
			# !!! DEBUG !!!
			
			if bnode.region in irradiated_boundaries
				@inline G1_vis, G1_IR = radiosity_window(f,u,bnode,data)

				hflux_irrad = (-eps2*ph"σ"*u[iT]^4 + alpha2_vis*G1_vis + alpha2_IR*G1_IR) # irradiation heatflux 
				#hflux_irrad = (-eps2*ph"σ"*(u[iT]^4-Tamb^4) + alpha2_vis*G1_vis) # irradiation heatflux 
				
				Tm=0.5*(u[iT] + u[iTw]) # mean temperature: catalyst layer and window
		        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X0) # thermal conductivity at Tm and inlet composition X0
				dh=2*uc_h
				kconv=Nu*λf/dh*ufac"W/(m^2*K)"
				hflux_conv = kconv*(u[iT]-Tm) # convection heatflux through top chamber
				
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
				G2_IR = eps2*ph"σ"*u[iT]^4 + rho2_IR*G1_IR
				hflux_abs_w = alpha1_vis*G2_vis + alpha1_IR*G2_IR + alpha1_IR*ph"σ"*Tamb^4
				hflux_emit_w = uc_window.eps*ph"σ"*u[iTw]^4
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
    reactor (PCR) model.
"""
function PCR_side(f,u,bnode,data)
	(;iT,k_nat_conv,delta_gap,Tamb,side_boundaries,solve_T_equation,dt_hf_enth)=data
	ng=ngas(data)

	if solve_T_equation && bnode.region in side_boundaries
	    X=MVector{ng,eltype(u)}(undef)
        @inline MoleFrac!(X,u,data)
        @inline _,λf=dynvisc_thermcond_mix(data, u[iT], X)

		# !!! DEBUG !!!
        # w/o shell height
        f[iT] = (u[iT]-Tamb)/(delta_gap/λf+1/k_nat_conv) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_enth)
        # f[iT] = zero(eltype(u))		
		# !!! DEBUG !!!
    end
end

@doc raw"""
Function defining the bottom/outlet boundary condition in the photo thermal catalytic
    reactor (PCR) model.
"""
function PCR_bottom(f,u,bnode,data)
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
		hflux_irrad = -eps1*ph"σ"*u[iT]^4 + alpha1_IR/(1-rho1_IR*rho2_IR)*(eps2*ph"σ"*Tplate^4+rho2_IR*eps1*ph"σ"*u[iT]^4)
	
		# heatflux from convective/conductive heat exchange with bottom plate
		Tm=0.5*(u[iT] + Tplate)

		X=MVector{ng,eltype(u)}(undef)
		@inline MoleFrac!(X,u,data)
        
        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X) # thermal conductivity at Tm and outlet composition X
        dh=2*lc_h
		kconv=Nu*λf/dh*ufac"W/(m^2*K)"
		hflux_conv = kconv*(u[iT]-Tm) # positive flux in negative z coord. (towards plate)

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
		hflux_emit_p = lc_plate.eps*ph"σ"*u[iTp]^4
		G1_IR = (eps1*ph"σ"*u[iT]^4 + rho1_IR*eps2*ph"σ"*u[iTp]^4)/(1-rho1_IR*rho2_IR)
		hflux_abs_p = alpha2_IR*G1_IR + alpha2_IR*ph"σ"*Tamb^4

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
    reactor (PCR) model. Used with the VoronoiFVM.jl package
"""
function PCR_bcond(f,u,bnode,data)
	(;p,ip,outlet_boundaries)=data
	
    FixedBed.PCR_top(f,u,bnode,data)
    FixedBed.PCR_side(f,u,bnode,data)
    FixedBed.PCR_bottom(f,u,bnode,data)
    for boundary in outlet_boundaries
        boundary_dirichlet!(f,u,bnode, species=ip,region=boundary,value=p)
    end	
end

@doc raw"""
Function defining the flux function in the boundary species in the photo thermal catalytic
    reactor (PCR) model. Calculates the heat flux occuring in the window and 
    bottom plate respectively. Only relevant if temperature equation is solved.
"""
function PCR_bflux(f,u,bedge,data)
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
    thermal catalytic reactor (PCR) model. Calculates the heat storage capacity
    for transient simulations in the window and bottom plate respectively.
    Only relevant if temperature equation is solved.
"""
function PCR_bstorage(f,u,bnode,data)
	(;irradiated_boundaries,outlet_boundaries)=data
	if bnode.region in irradiated_boundaries # window, upper chamber
		(;iTw) = data
		f[iTw] = u[iTw]
	elseif bnode.region in outlet_boundaries # bottom plate, lower chamber
		(;iTp) = data
		f[iTp] = u[iTp]
	end
end


function grid_boundaries_regions(dim;nref=0)
	Ω_catalyst = 2
	W=16
	H=0.5

	nz = 10*2^(nref)
	Z=(0:(H/nz):H)*ufac"cm"

	if dim == 2
		Γ_bottom = 1
		#Γ_bottom_insulating = 7
		Γ_side = 2
		Γ_sym = 4		
		Γ_top_permeable = 5
		Γ_top_irradiated = 6

		W=W/2 # axisymmetry, half domain is sufficient
		nr=W*2^(nref)
		R=(0:(W/nr):W)*ufac"cm"
		
		grid=simplexgrid(R,Z)
		circular_symmetric!(grid)
	
		cellmask!(grid,[0,9/10*H].*ufac"cm",[W,H].*ufac"cm",Ω_catalyst) # catalyst layer	
		bfacemask!(grid, [0,H].*ufac"cm",[W-1,H].*ufac"cm",Γ_top_permeable)
		bfacemask!(grid, [0,H].*ufac"cm",[W-2,0.5].*ufac"cm",Γ_top_irradiated) 
				
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
		X=(0:(W/nxy):W)*ufac"cm"
		Y=X
		# X=(0:1:W)*ufac"cm"
		# Y=(0:1:W)*ufac"cm"
		# Z=(0:H/10:H)*ufac"cm"	
		grid=simplexgrid(X,Y,Z)
	
		# catalyst region
		cellmask!(grid,[0,0,9/10*H].*ufac"cm",[W,W,H].*ufac"cm",Ω_catalyst) # catalyst layer	
		bfacemask!(grid, [1,1,H].*ufac"cm",[W-1,W-1,H].*ufac"cm",Γ_top_permeable)
		bfacemask!(grid, [2,2,H].*ufac"cm",[W-2,W-2,H].*ufac"cm",Γ_top_irradiated)

		inb = [Γ_top_permeable,Γ_top_irradiated]
		irrb = [Γ_top_irradiated]
		outb = [Γ_bottom]
		sb = [Γ_side_1,Γ_side_2,Γ_side_3,Γ_side_4]
	end

	return grid, inb, irrb, outb, sb, [Ω_catalyst]
end;

function init_system(dim, grid, data::ReactorData)

	(;p,ip,Tamb,iT,iTw,iTp,ibf,inlet_boundaries,irradiated_boundaries,outlet_boundaries,catalyst_regions,FluxIntp,ng,X0,solve_T_equation)=data

	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=FixedBed.DMS_flux,
							reaction=FixedBed.DMS_reaction,
							storage=FixedBed.DMS_storage,
							bcondition=FixedBed.PCR_bcond,
							bflux=FixedBed.PCR_bflux,
							bstorage=FixedBed.PCR_bstorage,
							boutflow=FixedBed.DMS_boutflow,
							outflowboundaries=outlet_boundaries,
							assembly=:edgewise,
							)

	if solve_T_equation
		enable_species!(sys; species=collect(1:(ng+2))) # gas phase species xi, ptotal & T
		enable_boundary_species!(sys, iTw, irradiated_boundaries) # window temperature as boundary species in upper chamber

		# for 3 dimensional domain, apply measured irradiation flux density as boundary condition
		if dim == 3
			# boundary flux species, workaround to implement spatially varying irradiation
			enable_boundary_species!(sys, ibf, irradiated_boundaries)
		end
		enable_boundary_species!(sys, iTp, outlet_boundaries) # plate temperature as boundary species in lower chamber
	else
		enable_species!(sys; species=collect(1:(ng+1))) # gas phase species xi, ptotal
	end

	inival=unknowns(sys)
	inival[ip,:].=p

	for i=1:ng
		inival[i,:] .= X0[i]
	end

	if solve_T_equation
		inival[[iT,iTw,iTp],:] .= Tamb
		if dim == 3
			function d3tod2(a,b)
				a[1]=b[1]
				a[2]=b[2]
			end
			inival[ibf,:] .= 0.0
			sub=ExtendableGrids.subgrid(grid,irradiated_boundaries,boundary=true, transform=d3tod2 )
				
			for inode in sub[CellNodes]
				c = sub[Coordinates][:,inode]
				inodeip = sub[ExtendableGrids.NodeInParent][inode]
				inival[ibf,inodeip] = FluxIntp(c[1]-0.02, c[2]-0.02)
			end
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
function hsf(data,T,p,x)
	Re,Pr,_ = RePrPe(data,T,p,x)
	_,lambdaf = dynvisc_thermcond_mix(data, T, x)
    (;poros,d) = data
	lambdaf/d*((1.0 + 4*(1.0-poros)/poros) + 0.5*(1.0-poros)^0.5*Re^0.6*Pr^(1.0/3.0))*ufac"W/(m^2*K)"
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


