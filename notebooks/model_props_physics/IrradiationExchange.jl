### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 8d82851a-4313-4cb9-aa5a-8ee15ac40657
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))	

	using NLsolve, ForwardDiff
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

# ╔═╡ df21e384-2691-4feb-a076-9288dcb9c8f5
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

# ╔═╡ c0111aa7-6f49-4b6b-845c-9d6c089383d3
# optical parameters for quartz window in upper chamber (uc)
const uc_window = SurfaceOpticalProps(
	# quartz glass windows
	alpha_IR=0.2, # datasheet integrated over spectrum of HFSS
	tau_IR=0.5,# datasheet integrated over spectrum of HFSS
	alpha_vis=0.0, # quartz glass is transparent in visible spectrum
	tau_vis=0.9, # datasheet, take reflection losses at adouble interface into account
	# the calculated value for rho_vis is probably inaccurate since it depends on incidence angle, which differs between the irradiation coming from the lamp (relevant for tau_vis) and the diffuse reflected light coming from the catalyst 
)

# ╔═╡ aeca609b-e6da-4e81-85bd-5fcd6d1e694f
#  optical parameters for catalyst layer in upper chamber (uc)
const uc_cat = SurfaceOpticalProps(
	# quartz glass windows
	alpha_IR=0.7, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.7, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# ╔═╡ b8e086b6-7def-452f-a625-5b0c498acf1d
# optical parameters for uncoated frit in upper chamber (uc)
const uc_frit = SurfaceOpticalProps(
	# quartz glass windows
	alpha_IR=0.2, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.2, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# ╔═╡ df2a7787-fa68-4f82-b8c4-b7078de80a97
function EquilWindowTemperature(data, T2, T3)

	(;ng,iT,Glamp,X0,α_nat_conv,Tamb,uc_window,uc_cat,uc_frit,vf_uc_window_cat,vf_uc_window_frit,uc_h)=data

	σ=ph"σ"

	
	# irradiation exchange between quartz window (1), cat surface (2) & frit surface (3)
	# window properties (1)
	tau1_vis=uc_window.tau_vis
	rho1_vis=uc_window.rho_vis
	rho1_IR=uc_window.rho_IR
	alpha1_IR=uc_window.alpha_IR
	alpha1_vis=uc_window.alpha_vis
	eps1=uc_window.eps
	# catalyst layer properties (2)
	rho2_vis=uc_cat.rho_vis
	rho2_IR=uc_cat.rho_IR
	alpha2_vis=uc_cat.alpha_vis
	alpha2_IR=uc_cat.alpha_IR
	eps2=uc_cat.eps
	# uncoated frit properties (3)
	rho3_vis=uc_frit.rho_vis
	rho3_IR=uc_frit.rho_IR
	alpha3_vis=uc_frit.alpha_vis
	alpha3_IR=uc_frit.alpha_IR
	eps3=uc_frit.eps
	#view factors
	ϕ12=vf_uc_window_cat
	ϕ13=vf_uc_window_frit
	
	function f!(F, x)		
		
		# surface brigthness of quartz window (1) in vis & IR	
		G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))
		# here the simplification is applied: only local value of T (u[iT]) is available, it is used for both surfaces
		G1_bot_IR = (eps1*σ*x[1]^4 + rho1_IR*(ϕ12*eps2*σ*T2+ϕ13*eps3*σ*T3^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))		

        G2_vis = rho2_vis*G1_bot_vis
        G2_IR = eps2*σ*T2^4 + rho2_IR*G1_bot_IR

        G3_vis = rho3_vis*G1_bot_vis
        G3_IR = eps3*σ*T3^4 + rho3_IR*G1_bot_IR

		q_abs_23_IR = alpha1_IR*(ϕ12*G2_IR+ϕ13*G3_IR)
		q_abs_23_vis = alpha1_vis*(ϕ12*G2_vis+ϕ13*G3_vis)

		# conductive heat flux through top chamber
		# mean temperature
        Tm2=0.5*(T2 + x[1])
		Tm3=0.5*(T3 + x[1])
        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        _,λf2=dynvisc_thermcond_mix(data, Tm2, X0)
		_,λf3=dynvisc_thermcond_mix(data, Tm3, X0)

        q_cond_bot = -λf2*(x[1]-T2)/uc_h*ϕ12 -λf3*(x[1]-T3)/uc_h*ϕ13
		
		q_conv_top = α_nat_conv*(x[1]-Tamb)

		q_emit = eps1*σ*x[1]^4		
		
    	F[1] = -q_conv_top + q_cond_bot - 2*q_emit + q_abs_23_IR + q_abs_23_vis
	end

	nlsolve(f!, [573.15], autodiff=:forward)
end

# ╔═╡ db854ec4-8bb0-41a2-a870-3b9827ee3579
md"""
## Energy Balance: Window
"""

# ╔═╡ 52e8545f-9efa-4cd3-aefd-165eadeec95a
md"""
$(LocalResource("../img/WindowEB.png")) 
"""

# ╔═╡ 4f7e75ca-285f-45c7-b214-10406a321ccd
md"""
At steady state, the net transfer of energy to/from the balance region around the quartz window must vanish. This leads to:

```math
0=-\dot q_{\text{conv,top}} - 2 \dot q_{\text{emit}} + \dot q_{\text{cond,bot,2}} + \dot q_{\text{cond,bot,3}} + \dot q_{\text{abs,bot,2}} + \dot q_{\text{abs,bot,3}}
```

with
```math
	\begin{align}
	\dot q_{\text{conv,top}} &= \alpha_{\text{conv}} (T_1-T_{\text{amb}}) \\
	\dot q_{\text{emit}} &= \epsilon_1\sigma T_1^4 \quad(\text{once for top \& bottom surfaces}) \\
	\dot q_{\text{cond,bot,2}} &= -\lambda(\frac{T_1+T_2}{2})\frac{T_1-T_2}{\Delta z} \\
	\dot q_{\text{cond,bot,3}} &= -\lambda(\frac{T_1+T_3}{2})\frac{T_1-T_3}{\Delta z} \\
	\dot q_{\text{abs,bot,2}} &= \alpha_1^{\text{IR}} \phi_{12} G_2^{\text{IR}} + \alpha_1^{\text{vis}} \phi_{12} G_2^{\text{vis}} \\
	\dot q_{\text{abs,bot,3}} &= \alpha_1^{\text{IR}} \phi_{13} G_3^{\text{IR}} + \alpha_1^{\text{vis}} \phi_{13} G_3^{\text{vis}}
	\end{align}
```

From this balance, the steady-state window temperature $T_1$ can be obtained.
"""

# ╔═╡ 6c903b77-f899-4796-af3e-6cd66ba6561f
md"""
## Energy Balance: Bottom Chamber
"""

# ╔═╡ d33366ed-4a43-4274-b8a7-b589b0300f2f
md"""
$(LocalResource("../img/BotChamberEB.png")) 
"""

# ╔═╡ 9ce02c32-6812-4dfc-bf4f-43def0b2c374
md"""
At steady state, the net transfer of energy to/from the balance region around the bottom plate must vanish. This leads to:

```math
0=-\dot q_{\text{conv}} - 2 \dot q_{\text{emit}} - \dot q_{\text{cond}} + \dot q_{\text{abs,1}} 
```

with
```math
	\begin{align}
	\dot q_{\text{conv}} &= \alpha_{\text{conv}} (T_2-T_{\text{amb}}) \\

	\dot q_{\text{emit}} &= \epsilon_2\sigma T_2^4 \quad(\text{once for top \& bottom surfaces}) \\
	\dot q_{\text{cond}} &= -\lambda(\frac{T_1+T_2}{2})\frac{T_1-T_2}{\Delta z} \\
	
	\dot q_{\text{abs,1}} &= \alpha_2^{\text{IR}} G_1^{\text{IR}} \\
	
	\end{align}
```

The surface brigthness of the frit $G_1$ is computed via:

```math
	G_1^{\text{IR}} = \frac{\epsilon_1 \sigma T_1^4 + \rho_1^{\text{IR}} \epsilon_2 \sigma T_2^4}{1-\rho_1^{\text{IR}} \rho_2^{\text{IR}}}

```

We obtain the mean temperature of the frit $T_1$ from the simulation. From this balance, the steady-state bottom plate temperature $T_2$ can be obtained.
"""

# ╔═╡ 2c5003b1-1ca1-4749-bdf6-21d546f1da50
function runSolveBC(data)
	
	function f!(F, x)

		# x[1] = T2 = temperature of bottom plate
		(;ng,iT,gni,α_nat_conv,Tamb,lc_frit,lc_plate,lc_h)=data

		σ=ph"σ"
		T1=773.15

		# frit properties
		eps1=lc_frit.eps
		rho1_IR=lc_frit.rho_IR
		# plate properties
		alpha2_IR=lc_plate.alpha_IR
		eps2=lc_plate.eps
		rho2_IR=lc_plate.rho_IR
		
		# conductive heat flux through bottom chamber
		# mean temperature
        Tm=0.5*(T1 + x[1])
		
        # thermal conductivity at Tm and outlet composition X (CO/H2/CO2/H2O = 1/1/1/1)
		MoleFrac = let
			mf=zeros(Float64, ng)
			mf[gni[:H2]] = 1.0
			mf[gni[:CO2]] = 1.0
			mf[gni[:CO]] = 1.0
			mf[gni[:H2O]] = 1.0
			mf/sum(mf)
		end
        _,λf=dynvisc_thermcond_mix(data, Tm, MoleFrac)
		
        q_cond = -λf*(T1-x[1])/lc_h
		
		q_conv = α_nat_conv*(x[1]-Tamb)

		q_emit = eps2*σ*x[1]^4

		G1_IR = (eps1*σ*T1^4 + rho1_IR*eps2*σ*x[1]^4)/(1-rho1_IR*rho2_IR)

		q_abs_1 = alpha2_IR * G1_IR


		F[1] = -q_conv -q_cond - 2*q_emit + q_abs_1
	end

	nlsolve(f!, [573.15], autodiff=:forward)
end

# ╔═╡ f80199ae-8902-4cb4-b87d-731a51859c6c
# optical parameters for uncoated frit in lower chamber (lc) = frit in upper chamber
const lc_frit = uc_frit

# ╔═╡ 809361da-8e71-4028-b0a1-e25463c2d0aa
# optical parameters for Al bottom plate in lower chamber (lc)
const lc_plate = SurfaceOpticalProps(
	# aluminium bottom plate (machined surface)
	alpha_IR=0.1, # see in lit, assume large reflectivity
	tau_IR=0.0, # opaque surface
	alpha_vis=0.1, # see in lit, assume large reflectivity
	tau_vis=0.0 # opaque surface	
)

# ╔═╡ 8d0f463b-897b-491e-9876-322af9dc1614
Base.@kwdef mutable struct ModelData <:AbstractModelData
	
	# catalyst / chemistry data
	# kinetic parameters, S3P="simple 3 parameter" kinetics fit to UPV lab scale experimental data
	# kinpar::AbstractKineticsData = S3P
	kinpar::AbstractKineticsData = XuFroment1989
	
	# number of gas phase species
	ng::Int64		 		= S3P.ng
	ip::Int64=ng+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable

	# names and fluid indices
	gn::Dict{Int, Symbol} 	= S3P.gn

	# inverse names and fluid indices
	gni::Dict{Symbol, Int}  = S3P.gni
	# fluids and respective properties in system
	Fluids::Vector{FluidProps} = S3P.Fluids
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = let
		x=zeros(Float64, ng)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		x/sum(x)
	end # inlet composition

	# volume specific cat mass loading, UPV lab scale PC reactor
	mcats::Float64 =1234.568*ufac"kg/m^3"
	isreactive::Bool = 1
	#isreactive::Bool = 0
		
	α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	α_nat_conv::Float64=15.0*ufac"W/(m^2*K)" # natural convection heat transfer coefficient	
	
	## porous filter data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0
	#dp::Float64=100.0*ufac"μm" # average pore size, por class 2

	# frit thickness (applies to 2D & 3D)
	h::Float64=0.5*ufac"cm"
	# catalyst layer thickness (applies to 2D & 3D)
	cath::Float64 = 500.0*ufac"μm"

	# upper and lower chamber heights for calculation of conduction b.c.
	uc_h::Float64=17.0*ufac"mm"
	lc_h::Float64=18.0*ufac"mm"

	# prism / 3D
	wi::Float64=12.0*ufac"cm" # prism width/side lenght
	le::Float64=wi # prism width/side lenght
	catwi::Float64=10.0*ufac"cm" # prism width/side lenght	
	
	Ac::Float64=wi*le*ufac"m^2" # cross-sectional area, square

	## irradiation data
	Glamp::Float64=1.0*ufac"kW/m^2" # solar simulator irradiation flux

	# upper chamber: quartz window eff. optical properties
	uc_window::SurfaceOpticalProps = uc_window
	# upper chamber: catalyst layer eff. optical properties
	uc_cat::SurfaceOpticalProps = uc_cat
	# upper chamber: (uncoated) porous frir eff. optical properties
	uc_frit::SurfaceOpticalProps = uc_frit

	#view factors
	vf_uc_window_cat::Float64 = catwi^2/Ac # A2/A1
	vf_uc_window_frit::Float64 = 1-vf_uc_window_cat # A3/A1 = 1-A2/A1

	# lower chamber: (uncoated) porous frir eff. optical properties
	lc_frit::SurfaceOpticalProps = lc_frit
	# lower chamber:  Al bottom plate eff. optical properties
	lc_plate::SurfaceOpticalProps = lc_plate
	
	Tglass::Float64 = 373.15*ufac"K" # glass temperature on surface facing catalyst
	Tplate::Float64 = 373.15*ufac"K" # bottom Al plate temperature facing frit
	

	# Solid Boro-Silikatglas
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
	
	#ϕ::Float64=0.36 # porosity, class 2
	ϕ::Float64=0.33 # porosity, class 0
	
	
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

	Qflow::Float64=3400.0*ufac"ml/minute" # volumetric feed flow rate (sccm)
	#Qflow::Float64=50000.0*ufac"ml/minute" # volumetric feed flow rate (sccm)

	MWin::Float64 = molarweight_mix(Fluids, X0)
	mdotin::Float64=MWin*Qflow*pn/(ph"R"*Tn)*ufac"kg/s"
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	p::Float64=1.0*ufac"atm" # reactor pressure

	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	ubot::Float64=u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate
	
end;

# ╔═╡ f6c41087-cdb2-405f-8b31-936d8d5b7a50
let
	EquilWindowTemperature(ModelData(),	873.15, 773.15)
end

# ╔═╡ f04ea704-6a18-47e8-b2a5-b673f7f0862c
let
	runSolveBC(ModelData())
end

# ╔═╡ df701006-d6c1-428a-95b7-bfa0a1ccf11e
ModelData()

# ╔═╡ Cell order:
# ╠═8d82851a-4313-4cb9-aa5a-8ee15ac40657
# ╟─50237926-dd49-45db-8ebe-63d9ba98dc9a
# ╠═df21e384-2691-4feb-a076-9288dcb9c8f5
# ╠═c0111aa7-6f49-4b6b-845c-9d6c089383d3
# ╠═aeca609b-e6da-4e81-85bd-5fcd6d1e694f
# ╠═b8e086b6-7def-452f-a625-5b0c498acf1d
# ╠═df2a7787-fa68-4f82-b8c4-b7078de80a97
# ╠═f6c41087-cdb2-405f-8b31-936d8d5b7a50
# ╟─db854ec4-8bb0-41a2-a870-3b9827ee3579
# ╠═52e8545f-9efa-4cd3-aefd-165eadeec95a
# ╠═4f7e75ca-285f-45c7-b214-10406a321ccd
# ╟─6c903b77-f899-4796-af3e-6cd66ba6561f
# ╠═d33366ed-4a43-4274-b8a7-b589b0300f2f
# ╟─9ce02c32-6812-4dfc-bf4f-43def0b2c374
# ╠═2c5003b1-1ca1-4749-bdf6-21d546f1da50
# ╠═f04ea704-6a18-47e8-b2a5-b673f7f0862c
# ╠═df701006-d6c1-428a-95b7-bfa0a1ccf11e
# ╠═f80199ae-8902-4cb4-b87d-731a51859c6c
# ╠═809361da-8e71-4028-b0a1-e25463c2d0aa
# ╠═8d0f463b-897b-491e-9876-322af9dc1614