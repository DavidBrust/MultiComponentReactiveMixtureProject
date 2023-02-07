# reaction kinetics for RWGS and Sabatier reactions
# based on Vidal Vázquez, F., et al. (2017).
# "Catalyst Screening and Kinetic Modeling for CO Production by High Pressure and Temperature Reverse Water Gas Shift for Fischer–Tropsch Applications." Industrial & Engineering Chemistry Research 56(45): 13262-13272.

# paratemters fit to experimental data from FPC project partner UPV

#Species indices:
#1) CO
#2) H2
#3) CH4
#4) H2O
#5) CO2
#6) N2

# reactions
# R1: CO + H2O = CO2 + H2
# R2: CH4 + 2 H2O = CO2 + 4 H2
# R3: CH4 + H2O = CO + 3 H2

#Parameter values obtained for reaction R1 & R2, fit to experimental data from UPV

#    |    | k     | E     |
#    |----|-------|-------|
#    | R1 | -1.5  |   -   |
#    | R2 | -11.2 | 175.9 |
#    |    |       |       |
#    

# the simplified model cosists of 3 fitted model parameters, therefore the name: S3P
S3P_old = [
    [-1.51715, -11.2379, -Inf], # ki_ref
    [0.0, 175.904, 0.0], # Ei
    [0.0, 0.0296, 0.0, 0.0, 0.0, 0.0] # Kj_ref
]

abstract type AbstractKineticsData end

Base.@kwdef mutable struct KineticsData <:AbstractKineticsData
    ng::Int64=6 # number of gas phase species
    gnames::Vector{String} = ["CO", "H2", "CH4", "H2O", "CO2", "N2"]
    Fluids::Vector{AbstractFluidProps} = [CO,H2,CH4,H2O,CO2,N2]
    gn::Dict{Int, String} 	= Dict(1:ng .=> gnames)
	# inverse names and fluid indices
	gni::Dict{String, Int}  = Dict(value => key for (key, value) in gn)
    
    nr::Int64=3 # number of reactions
    rnames::Vector{String} = ["R1", "R2", "R3"]
    rn::Dict{Int, String} 	= Dict(1:nr .=> rnames)
	rni::Dict{String, Int}  = Dict(value => key for (key, value) in rn)

    # values of reaction rate constants @ Tref
	ki_ref::Dict{String,Float64} = Dict("R1"=>-1.51715, "R2"=>-11.2379, "R3"=>-Inf)
    Tki_ref::Dict{String,Float64} = Dict( rnames .=> [648.0, 648.0, 648.0]*ufac"K")
    # Activation energies
    Ei::Dict{String,Float64} = Dict( rnames .=> [0.0, 175.904, 0.0]*ufac"kJ/mol")

    # equilibrium constants Ki
    Ki_ref::Dict{String,Float64} = Dict( rnames .=> [10.18, 1.947e-3, 1.913e-4])
    TKi_ref::Dict{String,Float64} = Dict( rnames .=> [693.0, 693.0, 693.0]*ufac"K")
    # reaction enthalpies
    ΔHi::Dict{String,Float64} = Dict( rnames .=> [-37.92, 182.09, 220.01]*ufac"kJ/mol")
	

    # values of gas phase species adsorption coefficients @ Tref
	Kj_ref::Dict{String,Float64} = Dict( gnames .=> [0.0, 0.0296, 0.0, 0.0, 0.0, 0.0])
    TKj_ref::Dict{String,Float64} = Dict( gnames .=> [648.0, 648.0, 823.0, 823.0, 823.0, 823.0]*ufac"K" )
    ΔHj::Dict{String,Float64} = Dict( gnames .=> [-70.65, -82.9, -38.28, 88.68, 0, 0]*ufac"kJ/mol")
end


S3P = KineticsData()

function pi(u,vn,u0,T,pt)
	#pt*u./sum(u)
	v=volFlowRate(u,vn,u0,T,pt)
	#pt*u./sum(u, dims=1) # allow operation on matrix
	u./(v/60/1000)*R*T*1e-5 #bar 
end

# equilibrium constants Ki
function Ki(T)
	Ki_ref=[10.18, 1.947e-3, 1.913e-4]
	ΔHi=[-37.92, 182.09, 220.01]	
	TKi_ref=[693, 693, 693]	
	@. Ki_ref*exp(-1000*ΔHi/ph"R"*(1/T-1/TKi_ref))
end

function Ki(kindata::AbstractKineticsData, T)
    nr=kindata.nr
    Ki=zeros(Float64, nr)
    Ki_ref = kindata.Ki_ref
    TKi_ref = kindata.TKi_ref
    ΔHi = kindata.ΔHi
    rn = kindata.rn
    for i=1:nr
        Ki[i] = Ki_ref[rn[i]] * exp(-ΔHi[rn[i]]/ph"R"*(1.0/T - 1.0/TKi_ref[rn[i]]))
    end
	Ki
end

# kinetic pre-factors ki
function ki(T,par)
	ki_ref=exp.(par[1])
	Ei=par[2]
	#ki_ref = [exp(par[1]), exp(par[2]), 0.0] # params 1 & 2 in log space, disable R3
	#Ei=[0.0, par[3], 0.0] # no activation energy for r1: not kinetically limited
	Tki_ref=[648, 648, 648]
	@. ki_ref*exp(-1000*Ei/ph"R"*(1/T-1/Tki_ref))
end

function ki(kindata::AbstractKineticsData, T)
    nr=kindata.nr
    ki=zeros(Float64, nr)
    ki_ref = kindata.ki_ref
    Tki_ref = kindata.Tki_ref
    Ei = kindata.Ei
    rn = kindata.rn
    for i=1:nr
        ki[i] = exp(ki_ref[rn[i]]) * exp(-Ei[rn[i]]/ph"R"*(1.0/T - 1.0/Tki_ref[rn[i]]))
        # @. ki_ref*exp(-1000*Ei/ph"R"*(1/T-1/Tki_ref))
    end
	ki
end

#Species indices:
#1) CO
#2) H2
#3) CH4
#4) H2O
#5) CO2
#6) N2


# adsorption constants Kj
function Kj(T,par)
	Kj_ref=par[3]
	#Kj_ref=[0, 0.0296, 0, 0, 0, 0] # 1/bar turn DEN ~ 1
	ΔHj=[-70.65, -82.9, -38.28, 88.68, 0, 0] # kJ/mol
	TKj_ref=[648, 648, 823, 823, 823, 823] # K
	@. Kj_ref*exp(-1000*ΔHj/ph"R"*(1/T-1/TKj_ref))
end

function Kj(kindata::AbstractKineticsData, T)
    ng=kindata.ng
    Kj=zeros(Float64, ng)
    Kj_ref = kindata.Kj_ref
    TKj_ref = kindata.TKj_ref
    ΔHj = kindata.ΔHj
    gn = kindata.gn
    for j=1:ng
        Kj[j] = Kj_ref[gn[j]] * exp(-ΔHj[gn[j]]/ph"R"*(1.0/T - 1.0/TKj_ref[gn[j]]))
        
    end
	Kj
end

function DEN(T,p,par)
	p_=p
	p_[4]=p_[4]/p_[2] # p_H2O/p_H2
	1+sum(Kj(T,par).*p_)
end

function DEN(kindata::AbstractKineticsData,T,p)
	n=kindata.gni
    p_=p
	p_[n["H2O"]]=p_[n["H2O"]]/p_[n["H2"]] # p_H2O/p_H2
	1+sum(Kj(kindata, T).*p_)
end

function ri(T,p,par)
	K=Ki(T)
	ki(T,par) .* @views [
		1/p[2]*(p[1]*p[4]-p[2]*p[5]/K[1]),
		1/p[2]^3.5*(p[3]*p[4]^2-p[2]^4*p[5]/K[2]),
		1/p[2]^2.5*(p[3]*p[4]-p[2]^3*p[1]/K[3])
	] / DEN(T,p,par)^2
end




function ri(kindata::AbstractKineticsData,T,p)
    n=kindata.gni
    r=kindata.rni
	K=Ki(kindata, T)
	ki(kindata,T) .* @views [
		1/p[n["H2"]]*(p[n["CO"]]*p[n["H2O"]]-p[n["H2"]]*p[n["CO2"]]/K[r["R1"]]),
		1/p[n["H2"]]^3.5*(p[n["CH4"]]*p[n["H2O"]]^2-p[n["H2"]]^4*p[n["CO2"]]/K[r["R2"]]),
		1/p[n["H2"]]^2.5*(p[n["CH4"]]*p[n["H2O"]]-p[n["H2"]]^3*p[n["CO"]]/K[r["R3"]])
	] / DEN(kindata,T,p)^2
end