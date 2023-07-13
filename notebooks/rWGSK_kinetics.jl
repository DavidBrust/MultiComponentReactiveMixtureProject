### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 5f58cde0-0b62-11ee-3544-1524261952f1
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using VoronoiFVM, VoronoiFVM.SolverStrategies
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve
	using LinearAlgebra
	using StaticArrays
	
	using OrdinaryDiffEq

	using LessUnitful
	
	using PlutoVista, Plots
	using PlutoUI
	using CSV,DataFrames
	using Interpolations

	using FixedBed
	
	#GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 7ad57b26-f249-4397-8229-dfc2273c45e2
PlutoUI.TableOfContents(title="Kinetics for rWGS",depth=6)

# ╔═╡ 2945aa66-180e-46e8-ad5c-ac8bc4cdd3d9
md"""
# Riedel 2001
"""

# ╔═╡ 64ae775d-1e4c-47b4-8b7f-bfaf10e10199
md"""
Implement kinetics of promoted iron based catalyst for the reverse water gas shift reaction.
Utilize kinetics published in  __Riedel, T., et al. (2001).__ "Kinetics of CO2 Hydrogenation on a K-Promoted Fe Catalyst." Industrial & Engineering Chemistry Research 40(5): 1355-1363.
"""

# ╔═╡ be27eabe-abc7-406c-bbda-e81e46809ab6
md"""
## Rate Equation
$(LocalResource("../img/RWGS_RR_eq.png")) 
"""

# ╔═╡ 044ac89d-3b61-4d5d-ad29-ddb0127b52a4
md"""
## Temperature Dependence / Arrhenius
The temperature effect on the kinetic rate coefficient kj was assumed to follow the Arrhenius law:
$(LocalResource("../img/Arrhenius.png")) 
"""

# ╔═╡ edad4a40-969d-4122-83fe-bee6af9985fa
md"""
## Kinetic Constants
"""

# ╔═╡ ad6d2019-082c-44e9-826e-bd101a95647c
md"""
$(LocalResource("../img/kinetics_par.png")) 
"""

# ╔═╡ 8a012774-4d41-4c13-9de2-0790b2fd67de
md"""
## Equilibrium constant for WGS (forward)
Equatino from __Zimmerman, W. H.; Bukur, D. B.__ Reaction Kinetics Over
Iron Catalysts Used for the Fischer-Tropsch Synthesis. Can. J.
Chem. Eng. 1990, 68, 292.

$(LocalResource("../img/RWGS_Kequil.png"))
"""

# ╔═╡ 36aecfd9-dd75-4ab2-b192-f54f69310125
function K_WGS(T)
	logK = 2073.0/T-2.029
	exp(logK)
end

# ╔═╡ a41d4ab2-0911-41be-96c7-5efda6664706
function rr_Riedel2001(T,p)
	aH2O = 65.0 # H2O inhibition
	bCO2 = 7.4 # CO2 inhibition
		
	k0 = 1.51e7
	Ea = 55.0*ufac"kJ/mol"
	
	k = k0 * exp(-Ea/(ph"R"*T))
	Keq = K_WGS(T)

	pCO,pH2O,pCO2,pH2,pN2 = p

	# K_WGS = 1/K_rWGS
	k * (pCO2*pH2 - pCO*pH2O*Keq) / (pCO + aH2O*pH2O + bCO2*pCO2 ) * ufac"mol/(s*g)"
end

# ╔═╡ f8fbd14b-ca42-42da-8ac4-a28952336375
let
	T = 800 + 273.15
	p = [0.0, 0.0, 50.0, 50.0, 0.0]*ufac"kPa"
	rr_Riedel2001(T,p./ufac"MPa")
end

# ╔═╡ dcbf5972-ed70-4c65-85a7-a8db773d569b
md"""
# Wolf 2016
"""

# ╔═╡ 32bcf50e-b99d-4612-a325-027575110bc0
md"""
Implement kinetics of Al2O3 supported __nickel__ based catalyst (normally used for methane steam reforming) for the reverse water gas shift reaction.
Utilize kinetics published in  __Wolf, A., et al. (2016).__ "Syngas Production via Reverse Water-Gas Shift Reaction over a Ni-Al2O3 Catalyst: Catalyst Stability, Reaction Kinetics, and Modeling." Chemical Engineering & Technology 39(6): 1040-1048.
"""

# ╔═╡ 54f5deb8-bc11-4870-8a09-765ebb1ce7d4
md"""
## Rate Equation
In wolf 2016, intrinsic kinetics, obtained from experiments with crushed catalysts, were recorded.
"""

# ╔═╡ a3830cd2-5c29-4bd6-98f2-80840281811d
md"""
$(LocalResource("../img/rWGS_kinetics/Wolf2016/intrinsic_forwar_RR.png")) 
"""

# ╔═╡ 82234a76-903f-49fd-ba91-4bd69f970b39
md"""
### Net rate
Obtain the net rate by substracting the reverse rate (WGS) from the forward rate (rWGS)
"""

# ╔═╡ 52757d44-a2c3-41a6-abab-4afbe290c98a
md"""
$(LocalResource("../img/rWGS_kinetics/Wolf2016/net_rate.png")) 
"""

# ╔═╡ 1716cfca-0b51-41c0-943f-17bffcb6fd05
md"""
## Kinetic constants
"""

# ╔═╡ 687c5159-c91c-4f2d-aec7-0117e20f8d40
md"""
$(LocalResource("../img/rWGS_kinetics/Wolf2016/kinetic_constants.png")) 
"""

# ╔═╡ 384b1551-52ee-4c33-8ae2-50809f7923db
md"""
## Equilibrium Constant
Calculation according to formulation proposed by Twigg
"""

# ╔═╡ 39459b80-d652-4703-a2fc-2cb351a9c6ac
md"""
$(LocalResource("../img/rWGS_kinetics/Wolf2016/k_equil.png")) 
"""

# ╔═╡ 51ba6385-1ceb-4978-84cf-10247bb3978f
function K_rWGS_Twigg(T)
    Z = 1000.0/T -1
	1.0/exp(Z*(Z*(0.63508-0.29353*Z)+4.1778)+0.31688)	
end

# ╔═╡ 23147811-c812-4039-8a64-31b46023336b
md"""
Species indices:
1) CO
2) H2
3) CH4
4) H2O
5) CO2
6) N2
"""

# ╔═╡ 86273521-619c-4769-b3fd-60e56b7c5f62
#function rr_Wolf2016(T,c)	
function rr_Wolf2016(T,p)	
	
	k0 = 3100.0
	Ea = 82.0*ufac"kJ/mol"
	
	k = k0 * exp(-Ea/(ph"R"*T))
	#Keq = K_rWGS_Twigg(T)

	c=p/(ph"R"*T)

	#cCO2,cH2,cCO,cH2O,cN2 = c
	cCO,cH2,cCH4,cH2O,cCO2,cN2 = c

	cH2_in = 0.45*ufac"bar"/(ph"R"*T)
	
	k * cCO2*cH2^0.3 # mol/(s*kgcat)
	#k * cCO2*cH2_in^0.3 # mol/(s*kgcat)
	
end

# ╔═╡ ad2a01af-df95-4c7d-9756-2a5e0461b84c
#function rr_Wolf2016_fr(T,c)
function rr_Wolf2016_fr(T,p)	
		
	k0 = 3100.0
	Ea = 82.0*ufac"kJ/mol"
	
	k = k0 * exp(-Ea/(ph"R"*T))
	Keq = K_rWGS_Twigg(T)

	c=p/(ph"R"*T)

	#cCO2,cH2,cCO,cH2O,cN2 = c
	cCO,cH2,cCH4,cH2O,cCO2,cN2 = c
	
	k * (cCO2*cH2^0.3 - cCO*cH2O/cH2^0.7/Keq) # mol/(s*kgcat)


end

# ╔═╡ 3534a150-4af1-4ee1-b18c-94f523b66ae5
let
	T = 800 + 273.15
	p = [0.0, 50.0, 0.0, 0.0, 50.0, 0.0]*ufac"kPa"
	c = p/(ph"R"*T)
	@show rr_Wolf2016_fr(T,p)/c[2]
	
	# Wolf
	#data1=ModelData(kinpar=Wolf_rWGS)
	#ri(data1,T,p)

	# Xu /Froment
	#data2=ModelData(kinpar=XuFroment)
	#-ri(data2,T,p)[2]
	
	
end;

# ╔═╡ 9975c842-d660-474d-b6be-627141e609c2
md"""
# Xu/Froment
"""

# ╔═╡ 7a09f721-3445-4001-b96a-2b00bc36859e
begin
	
	Base.@kwdef mutable struct ModelData{NG} <:AbstractModelData
	
	# catalyst / chemistry data
	# kinetic parameters, S3P="simple 3 parameter" kinetics fit to UPV lab scale experimental data

	#kinpar::FixedBed.KinData{nreac(S3P)} = S3P
	#kinpar::FixedBed.KinData{nreac(XuFroment)} = XuFroment
	kinpar::FixedBed.KinData{T} where {T}= Riedel_rWGS
	#kinpar::FixedBed.KinData{nreac(Wolf_rWGS)} = Wolf_rWGS
	
    
	
	
	# number of gas phase species
	#ng::Int64		 		= kinpar.ng
	ip::Int64=NG+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable
	# register window & plate temperatures as boundary species
	iTw::Int64=iT+1 # index of window Temperature (upper chamber)
	iTp::Int64=iTw+1 # index of plate Temperature (lower chamber)

	# names and fluid indices
	gn::Dict{Int, Symbol} 	= kinpar.gn

	# inverse names and fluid indices
	gni::Dict{Symbol, Int}  = kinpar.gni
	# fluids and respective properties in system
	Fluids::Vector{FluidProps} = kinpar.Fluids
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = let
		x=zeros(Float64, NG)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		#x[gni[:N2]] = 1.0
		x/sum(x)
	end # inlet composition

	#mcat::Float64=3.3*ufac"g" # total catalyst loading
	mcat::Float64=200.0*ufac"mg"
	# volume specific cat mass loading, UPV lab scale PC reactor
	lcats::Float64 =1000.0*ufac"kg/m^3"
	#mcats::Float64=80.0*ufac"kg/m^3" # 200 mg cat total loading, 250μm CL 
	isreactive::Bool = 1
	#isreactive::Bool = 0
		
	#α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	k_nat_conv::Float64=10.0*ufac"W/(m^2*K)" # natural convection heat transfer coefficient	
	
	## porous filter data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0
	#dp::Float64=100.0*ufac"μm" # average pore size, por class 2

	# frit thickness (applies to 2D & 3D)
	h::Float64=0.5*ufac"cm"
	# catalyst layer thickness (applies to 2D & 3D)
	#cath::Float64 = 500.0*ufac"μm"
	cath::Float64 = 250.0*ufac"μm"

	# upper and lower chamber heights for calculation of conduction/convection b.c.
	uc_h::Float64=17.0*ufac"mm"
	lc_h::Float64=18.0*ufac"mm"
	Nu::Float64=4.861

	# prism / 3D
	wi::Float64=12.0*ufac"cm" # prism width/side lenght
	le::Float64=wi # prism width/side lenght
	catwi::Float64=10.0*ufac"cm" # prism width/side lenght	
	shellh::Float64=3.6*ufac"cm" # height of reactor shell contacting domain
	Ac::Float64=wi*le*ufac"m^2" # cross-sectional area, square
	
end;
	
	# !!!ALLOC Method to be called instead of data.ng
	FixedBed.ngas(::ModelData{NG}) where NG = NG
	
	# !!!ALLOC Additional constructo taking ng as parameter	
	ModelData(;ng=S3P.ng, kwargs...) = ModelData{ng}(;kwargs...)

end;

# ╔═╡ b409e381-7a86-4244-837e-cd145248bfb9
let
	T = 800 + 273.15
	p = [0.0, 50.0, 0.0, 0.0, 50.0, 0.0]*ufac"kPa"/ufac"bar"

	data=ModelData(kinpar=XuFroment)
	
	ri(data,T,p)
	
end

# ╔═╡ db283cc7-a2a2-4cc3-9412-862347141696
md"""
# Fixed Bed Reactor Model
Reformulate fixed bed reactor equation provided in __Wolf, A., et al. (2016).__ "Syngas Production via Reverse Water-Gas Shift Reaction over a Ni-Al2O3 Catalyst: Catalyst Stability, Reaction Kinetics, and Modeling." Chemical Engineering & Technology 39(6): 1040-1048.
"""

# ╔═╡ 61bb5aab-3c2b-4a45-99da-4acc420bc97d
md"""
```math
- \frac{1}{\rho_{\text{bed}}} \frac{dc_i}{dt}=r_i
```
"""

# ╔═╡ e8113023-5c2c-4db6-bbe2-153b07ee6181
md"""
$(LocalResource("../img/rWGS_kinetics/Wolf2016/InletComp.png")) 
"""

# ╔═╡ 8f8af545-0de2-449c-aa98-669f364b611e
md"""
Extracted data from article
"""

# ╔═╡ dd111a2f-2bf6-4294-b9fe-f27dc6c7f122
begin
	Wolf2016 = CSV.read("../data/Wolf2016/IntrinsicKin.csv", DataFrame, delim=";")
end;

# ╔═╡ 3033f460-9bab-4030-9ce0-3526e5385100
function Da(T)
	rho_bed = 1020.0*ufac"kg/m^3"

	k0 = 3100.0
	Ea = 82.0*ufac"kJ/mol"
	
	k = k0 * exp(-Ea/(ph"R"*T))

	cH2 = 0.45*ufac"bar"/(ph"R"*T) * ufac"mol/m^3"

	tau_550 = 1.45*ufac"ms"
	tau=tau_550*(550+273.15)/T
	k*cH2^0.3*tau*rho_bed
end

# ╔═╡ d80e8a7b-b6bb-43e1-846a-7c6f2cbddd4b
function X_Da(T)
	1-exp(-Da(T))
end

# ╔═╡ 6c07f849-1c42-4306-a09d-574189fa20b1
begin
	Jess_fr = CSV.read("../data/Wolf2016/IntrinsicKinFixedBed/Jess_fr.csv", DataFrame, delim=";", header=["T", "X"])
	Jess_f = CSV.read("../data/Wolf2016/IntrinsicKinFixedBed/Jess_f.csv", DataFrame, delim=";", header=["T", "X"])
end;

# ╔═╡ f7e53314-72bb-49fa-9a59-1f4114483518
let
	Ts = range(550,950, length=11)
	
	p=Plots.plot(xguide="Temperature / °C", yguide="CO2 Conversion", ylim=(0,1))

	#Plots.plot!(p,Wolf2016.T_C, Wolf2016.X, label="Wolf 2016 Fig. 6a) Fwd")
	
	Plots.plot!(p,Jess_fr.T, Jess_fr.X, label="Jess pers. comm. Fwd Rev")
	Plots.plot!(p,Jess_f.T, Jess_f.X, label="Jess pers. comm. Fwd")

	Plots.plot!(p,Ts, X_Da.(Ts.+273.13), label="X from Da", marker=:circle)

end

# ╔═╡ 6f1c20c1-471a-4e95-adfd-e71f05e22ac0
md"""
Species indices:
1) CO
2) H2
3) CH4
4) H2O
5) CO2
6) N2
"""

# ╔═╡ cd654166-54b0-40c4-9b13-23dc08035d60
function PBR!(du,u,par,t)

	T=par[1]
	rho_bed = 1020 *ufac"kg/m^3"
	
	# RR = rho_bed*rr_Wolf2016(T,u)
	RR = rho_bed*rr_Wolf2016(T,u)*ph"R"*T
	
	du[1] = RR # CO
	du[2] = -RR # H2
	du[3] = 0 # CH4
	du[4] = RR # H2O
	du[5] = -RR # CO2
	du[6] = 0 # N2
	
end

# ╔═╡ a2ba7bfd-13e8-4939-a95a-15a4df8a1035
function PBR_fr!(du,u,par,t)

	T=par[1]
	rho_bed = 1020 *ufac"kg/m^3"

	#RR = rho_bed*rr_Wolf2016_fr(T,u)
	RR = rho_bed*rr_Wolf2016_fr(T,u)*ph"R"*T
	
	du[1] = RR # CO
	du[2] = -RR # H2
	du[3] = 0 # CH4
	du[4] = RR # H2O
	du[5] = -RR # CO2
	du[6] = 0 # N2
	
end

# ╔═╡ 75ea3673-bf3d-4193-814e-fb85e842d5ec
let
	
    # CO2 , H2 , CO , H2O , N2
	#p0 = [0.15, 0.45, 0.0, 0.0, 0.4]*ufac"bar"
	
	# CO , H2 , CH4 , H2O , CO2, N2
	p0 = [0.0, 0.45, 0.0, 0.0, 0.15, 0.4]*ufac"bar"
	

	tau_550 = 1.45*ufac"ms"	
	
	XCO2_f = []
	XCO2_fr = []
	
	Ts = range(550,950, length=11)
	
	for T in Ts
		par=[T+273.15]

		c0 = p0/(ph"R" * (T+273.15))
		tau = tau_550 * (550+273.15)/(T+273.15)
		tspan = (0.0,tau)
		
		# prob_f = ODEProblem(PBR!,c0,tspan,par)
		prob_f = ODEProblem(PBR!,p0,tspan,par)
		sol_f = solve(prob_f,Tsit5(),reltol=1e-8,abstol=1e-8)
		

		# prob_fr = ODEProblem(PBR_fr!,c0,tspan,par)
		prob_fr = ODEProblem(PBR_fr!,p0,tspan,par)					
		sol_fr = solve(prob_fr,Tsit5(),reltol=1e-8,abstol=1e-8)
		
		push!(XCO2_f, (sol_f[1][5]-sol_f[end][5])/sol_f[1][5])
		push!(XCO2_fr, (sol_fr[1][5]-sol_fr[end][5])/sol_fr[1][5])
	end
	
	p=Plots.plot(xguide="Temperature / °C", yguide="CO2 Conversion", ylim=(0,1), size=(400,300))
	
	#Plots.plot!(p,Wolf2016.T_C, Wolf2016.X, label="Wolf 2016 Fig. 6a)")
	
	Plots.plot!(p,Ts, XCO2_f, label="reproduction Eq. 6)")
	Plots.plot!(p,Ts, XCO2_fr, label="reproduction Eq. 9)")


	Plots.plot!(p,Jess_f.T, Jess_f.X, label="Jess pers. comm. Fwd")
	Plots.plot!(p,Jess_fr.T, Jess_fr.X, label="Jess pers. comm. Fwd Rev")

	#Plots.savefig(p,"..\\data\\Wolf2016\\IntrinsicKinFixedBed\\Compare.png")
end

# ╔═╡ Cell order:
# ╠═5f58cde0-0b62-11ee-3544-1524261952f1
# ╟─7ad57b26-f249-4397-8229-dfc2273c45e2
# ╟─2945aa66-180e-46e8-ad5c-ac8bc4cdd3d9
# ╟─64ae775d-1e4c-47b4-8b7f-bfaf10e10199
# ╟─be27eabe-abc7-406c-bbda-e81e46809ab6
# ╟─044ac89d-3b61-4d5d-ad29-ddb0127b52a4
# ╟─edad4a40-969d-4122-83fe-bee6af9985fa
# ╟─ad6d2019-082c-44e9-826e-bd101a95647c
# ╟─8a012774-4d41-4c13-9de2-0790b2fd67de
# ╠═36aecfd9-dd75-4ab2-b192-f54f69310125
# ╠═a41d4ab2-0911-41be-96c7-5efda6664706
# ╠═f8fbd14b-ca42-42da-8ac4-a28952336375
# ╟─dcbf5972-ed70-4c65-85a7-a8db773d569b
# ╟─32bcf50e-b99d-4612-a325-027575110bc0
# ╟─54f5deb8-bc11-4870-8a09-765ebb1ce7d4
# ╟─a3830cd2-5c29-4bd6-98f2-80840281811d
# ╟─82234a76-903f-49fd-ba91-4bd69f970b39
# ╟─52757d44-a2c3-41a6-abab-4afbe290c98a
# ╟─1716cfca-0b51-41c0-943f-17bffcb6fd05
# ╟─687c5159-c91c-4f2d-aec7-0117e20f8d40
# ╟─384b1551-52ee-4c33-8ae2-50809f7923db
# ╟─39459b80-d652-4703-a2fc-2cb351a9c6ac
# ╠═51ba6385-1ceb-4978-84cf-10247bb3978f
# ╠═23147811-c812-4039-8a64-31b46023336b
# ╠═86273521-619c-4769-b3fd-60e56b7c5f62
# ╠═ad2a01af-df95-4c7d-9756-2a5e0461b84c
# ╠═3534a150-4af1-4ee1-b18c-94f523b66ae5
# ╟─9975c842-d660-474d-b6be-627141e609c2
# ╠═b409e381-7a86-4244-837e-cd145248bfb9
# ╟─7a09f721-3445-4001-b96a-2b00bc36859e
# ╟─db283cc7-a2a2-4cc3-9412-862347141696
# ╟─61bb5aab-3c2b-4a45-99da-4acc420bc97d
# ╟─e8113023-5c2c-4db6-bbe2-153b07ee6181
# ╟─8f8af545-0de2-449c-aa98-669f364b611e
# ╠═dd111a2f-2bf6-4294-b9fe-f27dc6c7f122
# ╠═d80e8a7b-b6bb-43e1-846a-7c6f2cbddd4b
# ╠═3033f460-9bab-4030-9ce0-3526e5385100
# ╠═75ea3673-bf3d-4193-814e-fb85e842d5ec
# ╠═f7e53314-72bb-49fa-9a59-1f4114483518
# ╠═6c07f849-1c42-4306-a09d-574189fa20b1
# ╟─6f1c20c1-471a-4e95-adfd-e71f05e22ac0
# ╠═cd654166-54b0-40c4-9b13-23dc08035d60
# ╠═a2ba7bfd-13e8-4939-a95a-15a4df8a1035
