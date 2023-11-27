abstract type AbstractModelData end

@doc raw"""
## Overall Mass Continuity
Mixture mass flow (overall mass flow) through the pore space of the porous medium. Mixture mass averaged velocity is calculated from Darcy equation. The void fraction (porosity) is given by $\epsilon$.

```math
\begin{align}
	\frac{\partial \epsilon \rho}{\partial t} + \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\end{align}
```


## Species Mass Continuity and Transport
```math
\begin{align}
	\frac{\partial \epsilon \rho_i}{\partial t} + \nabla \cdot \left( \vec \Phi_i + \rho_i \vec v \right ) - R_i &= 0 ~, \qquad i = 1 ... \nu \\
		\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= -\sum_{j=1 \atop j \neq i}^{\nu} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
		\sum_{i=1}^\nu x_i &= 1
\end{align}
```


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


function D_matrix!(D,data)
	ng=ngas(data)
	(;m) = data
	@inbounds for i=1:(ng-1)
		for j=(i+1):ng
            Dji = one(eltype(D)) *1.0e-5*ufac"m^2/s" * m[i]*m[j]
			@views D[j,i] = Dji
			@views D[i,j] = Dji
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


function kbed(data, thermal_cond_fluid)
	(;poros,thermal_cond_solid) = data
	B=1.25*((1.0-poros)/poros)^(10.0/9.0)
	kp=thermal_cond_solid/thermal_cond_fluid
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
	alpha_IR=0.3,
	tau_IR=0.7,	
	alpha_vis=0.0,
	tau_vis=0.9,
)

#  optical parameters for catalyst layer in upper chamber (uc)
const uc_cat = SurfaceOpticalProps(
	alpha_IR=0.45, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.45, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# optical parameters for uncoated frit in upper chamber (uc)
const uc_mask = SurfaceOpticalProps(
	alpha_IR=0.2, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.2, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# optical parameters for uncoated frit in lower chamber (lc) = frit in upper chamber
const lc_frit = uc_mask

# optical parameters for Al bottom plate in lower chamber (lc)
const lc_plate = SurfaceOpticalProps(
	# aluminium bottom plate (machined surface)
	alpha_IR=0.1, # see in lit, assume large reflectivity
	tau_IR=0.0, # opaque surface
	alpha_vis=0.1, # see in lit, assume large reflectivity
	tau_vis=0.0 # opaque surface	
)

Base.@kwdef mutable struct ModelData{NG}
    
	dt_mf::Tuple{Float64, Float64}=(0.0,1.0)
	dt_hf_enth::Tuple{Float64, Float64}=(2.0,10.0)
	dt_hf_irrad::Tuple{Float64, Float64}=(3.0,10.0)
    kinetics::Symbol = :XuFroment
	#kinpar::FixedBed.KinData{nreac(FixedBed.kin_obj(kinetics))} = FixedBed.kin_obj(kinetics)
    kinpar::FixedBed.KinData{} = XuFroment
    ng::Int64 = kinpar.ng
	mcat::Float64=500.0*ufac"mg"
	Vcat::Float64=1.0*ufac"m^2"*0.02*ufac"cm"
	lcat::Float64=mcat/Vcat

	# upper and lower chamber heights for calculation of conduction/convection b.c.
	uc_h::Float64=17.0*ufac"mm"
	lc_h::Float64=18.0*ufac"mm"
	Nu::Float64=4.861
	
	#isreactive::Bool = 0
	isreactive::Bool = 1
	
	#ip::Int64 = NG+1
    ip::Int64 = ng+1
	iT::Int64 = ip+1
	# inlcude window & plate temperatures as boundary species
	iTw::Int64=iT+1 # index of window Temperature (upper chamber)
	iTp::Int64=iTw+1 # index of plate Temperature (lower chamber)
	
	p::Float64 = 1.0*ufac"bar"
	Tamb::Float64 = 298.15*ufac"K"

	gn::Dict{Int, Symbol} = kinpar.gn # names and fluid indices
	gni::Dict{Symbol, Int} = kinpar.gni # inverse names and fluid indices
	Fluids::Vector{FluidProps} = kinpar.Fluids # fluids and respective properties in system
	
	m::Vector{Float64} = let
		#m=zeros(Float64, NG)
        m=zeros(Float64, ng)
		#for i=1:NG
        for i=1:ng
			m[i] = Fluids[i].MW
		end
		m
	end
	
	X0::Vector{Float64} = let
		#x=zeros(Float64, NG)
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
	#mfluxin::Float64 = 0.007*ufac"kg/(m^2*s)"
	#nfluxin::Float64 = mfluxin/mmix0

	G_lamp::Float64 = 70.0*ufac"kW/m^2"
	
	# optical properties of surfaces in IR/visible ranges
	uc_window::SurfaceOpticalProps = uc_window
	uc_cat::SurfaceOpticalProps = uc_cat
	uc_mask::SurfaceOpticalProps = uc_mask
	lc_frit::SurfaceOpticalProps = lc_frit
	lc_plate::SurfaceOpticalProps = lc_plate

	# reactor data
	delta_gap::Float64=1.5*ufac"mm" # gas gap between frit and reactor wall
	k_nat_conv::Float64=17.5*ufac"W/(m^2*K)" # natural convection heat transfer coeff.

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

    ModelData(dt_mf,dt_hf_enth,dt_hf_irrad,kinetics,kinpar,ng,mcat,Vcat,lcat,uc_h,lc_h,Nu,isreactive,ip,iT,iTw,iTp,p,Tamb,gn,gni,Fluids,m,X0,mmix0,W0,nflowin,mflowin,mfluxin,G_lamp,uc_window,uc_cat,uc_mask,lc_frit,lc_plate,delta_gap,k_nat_conv,dp,poros,perm,γ_τ,rhos,lambdas,cs,lambda_window,lambda_Al) = 
    new{ng}(dt_mf,dt_hf_enth,dt_hf_irrad,kinetics,kinpar,ng,mcat,Vcat,lcat,uc_h,lc_h,Nu,isreactive,ip,iT,iTw,iTp,p,Tamb,gn,gni,Fluids,m,X0,mmix0,W0,nflowin,mflowin,mfluxin,G_lamp,uc_window,uc_cat,uc_mask,lc_frit,lc_plate,delta_gap,k_nat_conv,dp,poros,perm,γ_τ,rhos,lambdas,cs,lambda_window,lambda_Al)
end

ngas(::ModelData{NG}) where NG = NG


function DiffCoeffsMass(ng,m)
	k=0
	D=zeros(Float64, ng,ng)
	for i=1:(ng-1)
		for j=(i+1):ng
			k +=1
			#Dji = k*1.0e-5*ufac"m^2/s"
			Dji = k*1.0e-5*ufac"m^2/s"
			Dji *= m[i]*m[j]
			D[j,i] = Dji
			D[i,j] = Dji
		end			
	end
	D
end

mutable struct MD_test{NG, KP}
    
    kinetics::Symbol
    isreactive::Int64
    kinpar::KP
    ng::Int64
    dt_mf::Tuple{Float64, Float64}
	dt_hf_enth::Tuple{Float64, Float64}
	dt_hf_irrad::Tuple{Float64, Float64}
    mcat::Float64
	Vcat::Float64
	lcat::Float64
	uc_h::Float64
	lc_h::Float64
	Nu::Float64
    ip::Int64
    iT::Int64
	iTw::Int64
	iTp::Int64
	p::Float64
	Tamb::Float64
    gn::Dict{Int, Symbol}
	gni::Dict{Symbol, Int}
	Fluids::Vector{FluidProps}
	m::Vector{Float64}    
	X0::Vector{Float64}	
	mmix0::Float64
	W0::Vector{Float64}	
	nflowin::Float64
	mflowin::Float64
	mfluxin::Float64
    G_lamp::Float64
    uc_window::SurfaceOpticalProps
	uc_cat::SurfaceOpticalProps
	uc_mask::SurfaceOpticalProps
	lc_frit::SurfaceOpticalProps
	lc_plate::SurfaceOpticalProps
    delta_gap::Float64
	k_nat_conv::Float64
    dp::Float64
	poros::Float64
	perm::Float64
	constriction_tourtuosity_fac::Float64
    density_solid::Float64
	thermal_cond_solid::Float64
	heatcap_solid::Float64
    thermal_cond_window::Float64
	thermal_cond_reactorwall::Float64	

    function MD_test(;
        kinetics=:XuFroment,
        isreactive=1,
        ng=nothing,
        dt_mf=(0.0,1.0),
        dt_hf_enth=(2.0,10.0),
        dt_hf_irrad=(3.0,10.0),
        mcat=500.0*ufac"mg",
        Vcat=100.0*ufac"cm^2"*1.0*ufac"mm",
        lcat=mcat/Vcat,
        uc_h=17.0*ufac"mm",
        lc_h=18.0*ufac"mm",
        Nu=4.861,
        p=1.0*ufac"bar",
        Tamb=298.15*ufac"K",
        m=nothing,
        X0=nothing,        
        nflowin=7.4*ufac"mol/hr",
        G_lamp=70.0*ufac"kW/m^2",
        uc_window=uc_window,
        uc_cat=uc_cat,
        lc_frit=lc_frit,
        lc_plate=lc_plate,
        delta_gap=1.5*ufac"mm",
        k_nat_conv=17.5*ufac"W/(m^2*K)",
        dp=200.0*ufac"μm",
        poros=0.33,
        perm=1.23e-11*ufac"m^2",
        constriction_tourtuosity_fac=poros^1.5,
        density_solid=2.23e3*ufac"kg/m^3",
        thermal_cond_solid=1.13*ufac"W/(m*K)",
        heatcap_solid=0.8e3*ufac"J/(kg*K)",
        thermal_cond_window=1.38*ufac"W/(m*K)",
        thermal_cond_reactorwall=235.0*ufac"W/(m*K)"
        )

        kinpar = FixedBed.kin_obj(kinetics)
        KP = FixedBed.KinData{nreac(kinpar)}
        ng = isnothing(ng) ? kinpar.ng : ng
        ip = ng+1
        iT = ip+1
        iTw = iT+1
        iTp = iTw+1
        (;gn,gni,Fluids) = kinpar
        if isnothing(m)
            m=Vector{Float64}(undef,ng)
            for i=1:ng
                m[i] = Fluids[i].MW
            end
        end
        if isnothing(X0)
            X0=zeros(Float64, ng)
            X0[gni[:H2]] = 1.0
            X0[gni[:CO2]] = 1.0
            #X0[gni[:N2]] = 1.0
            X0 /= sum(X0)
        end
        @assert length(X0) == length(m) "X0 and m must have same length"
        mmix0 = sum(X0 .* m)
        W0 = @. m*X0/mmix0

        mflowin = nflowin*mmix0
        mfluxin = mflowin/(100*ufac"cm^2")*ufac"kg/(m^2*s)"

        new{ng,KP}(
            kinetics,
            isreactive,
            kinpar,
            ng,
            dt_mf,
            dt_hf_enth,
            dt_hf_irrad,
            mcat,
            Vcat,
            lcat,
            uc_h,
            lc_h,
            Nu,
            ip,
            iT,
            iTw,
            iTp,
            p,
            Tamb,
            gn,
            gni,
            Fluids,
            m,
            X0,
            mmix0,
            W0,
            nflowin,
            mflowin,
            mfluxin,
            G_lamp,
            uc_window,
            uc_cat,
            uc_mask,
            lc_frit,
            lc_plate,
            delta_gap,
            k_nat_conv,
            dp,
            poros,
            perm,
            constriction_tourtuosity_fac,
            density_solid,
            thermal_cond_solid,
            heatcap_solid,
            thermal_cond_window,
            thermal_cond_reactorwall)
    end

end

ngas(::MD_test{NG,KP}) where {NG,KP} = NG

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


