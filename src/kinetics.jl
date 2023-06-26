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
# S3P_old = [
#     [-1.51715, -11.2379, -Inf], # ki_ref
#     [0.0, 175.904, 0.0], # Ei
#     [0.0, 0.0296, 0.0, 0.0, 0.0, 0.0] # Kj_ref
# ]



abstract type AbstractKineticsData end



# Base.@kwdef mutable struct KineticsData{F<:Function} <:AbstractKineticsData
#     ng::Int64=6 # number of gas phase species
#     gnames::Vector{Symbol} = [:CO, :H2, :CH4, :H2O, :CO2, :N2]
#     Fluids::Vector{FluidProps} = [CO,H2,CH4,H2O,CO2,N2]
#     gn::Dict{Int, Symbol} 	= Dict(1:ng .=> gnames)
# 	# inverse names and fluid indices
# 	gni::Dict{Symbol, Int}  = Dict(value => key for (key, value) in gn)
    
#     nr::Int64=3 # number of reactions
#     rnames::Vector{Symbol} = [:R1, :R2, :R3]
#     rn::Dict{Int, Symbol} 	= Dict(1:nr .=> rnames)
# 	rni::Dict{Symbol, Int}  = Dict(value => key for (key, value) in rn)
#     nuij::Array{Int, 2} = vcat(
#         [-1 0 1], # CO
#         [1 4 3], # H2
#         [0 -1 -1], # CH4
#         [-1 -2 -1], # H2O
#         [1 1 0], # CO2
#         [0 0 0], # N2
#     )
#     RR::F = RRS3P # RR function
#     #RR::F = RRS3P_ # RR function
#     # values of reaction rate constants @ Tref
# 	ki_ref::Dict{Symbol,Float64} = Dict(:R1=>-1.51715, :R2=>-11.2379, :R3=>-Inf)
#     Tki_ref::Dict{Symbol,Float64} = Dict( rnames .=> [648.0, 648.0, 648.0]*ufac"K")
#     # Activation energies
#     Ei::Dict{Symbol,Float64} = Dict( rnames .=> [0.0, 175.904, 0.0]*ufac"kJ/mol")

#     # equilibrium constants Ki
#     Ki_ref::Dict{Symbol,Float64} = Dict( rnames .=> [10.18, 1.947e-3, 1.913e-4])
#     TKi_ref::Dict{Symbol,Float64} = Dict( rnames .=> [693.0, 693.0, 693.0]*ufac"K")
#     # reaction enthalpies
#     ΔHi::Dict{Symbol,Float64} = Dict( rnames .=> [-37.92, 182.09, 220.01]*ufac"kJ/mol")
	

#     # values of gas phase species adsorption coefficients @ Tref
# 	Kj_ref::Dict{Symbol,Float64} = Dict( gnames .=> [0.0, 0.0296, 0.0, 0.0, 0.0, 0.0])
#     TKj_ref::Dict{Symbol,Float64} = Dict( gnames .=> [648.0, 648.0, 823.0, 823.0, 823.0, 823.0]*ufac"K" )
#     ΔHj::Dict{Symbol,Float64} = Dict( gnames .=> [-70.65, -82.9, -38.28, 88.68, 0, 0]*ufac"kJ/mol")
# end

# Base.@kwdef struct Foo{NR}
#     nr::Int64=3
#     gnames::Vector{Symbol} = [:CO, :H2, :CH4, :H2O, :CO2, :N2]
#     nuij::Vector{Int64} = vcat(
#         [-1, 1, 0, -1, 1, 0], #R1 : CO + H2O -> CO2 + H2
#         [0, 4, -1, -2, 1, 0], #R2 : CH4 + 2 H2O -> 4 H2 + CO2
#         [1, 3, -1, -1, 0, 0], #R3 : CH4 + H2O -> 3 H2 + CO
#     )
#     Foo(nr,gnames,nuij) = new{nr}(nr,gnames,nuij)
# end

Base.@kwdef struct KinData{NR} <:AbstractKineticsData
    ng::Int64=6 # number of gas phase species
    gnames::Vector{Symbol} = [:CO, :H2, :CH4, :H2O, :CO2, :N2]
    Fluids::Vector{FluidProps} = [CO,H2,CH4,H2O,CO2,N2]
    gn::Dict{Int, Symbol} 	= Dict(1:ng .=> gnames)
	# inverse names and fluid indices
	gni::Dict{Symbol, Int}  = Dict(value => key for (key, value) in gn)
    
    nr::Int64=3 # number of reactions
    rnames::Vector{Symbol} = [:R1, :R2, :R3]
    rn::Dict{Int, Symbol} 	= Dict(1:nr .=> rnames)
	rni::Dict{Symbol, Int}  = Dict(value => key for (key, value) in rn)

    nuij::Vector{Int64} = vcat(
        [-1, 1, 0, -1, 1, 0], #R1 : CO + H2O -> CO2 + H2
        [0, 4, -1, -2, 1, 0], #R2 : CH4 + 2 H2O -> 4 H2 + CO2
        [1, 3, -1, -1, 0, 0], #R2 : CH4 + H2O -> 3 H2 + CO
    )

    # values of reaction rate constants @ Tref
	ki_ref::Dict{Symbol,Float64} = Dict(:R1=>-1.51715, :R2=>-11.2379, :R3=>-Inf)
    Tki_ref::Dict{Symbol,Float64} = Dict( rnames .=> [648.0, 648.0, 648.0]*ufac"K")
    # Activation energies
    Ei::Dict{Symbol,Float64} = Dict( rnames .=> [0.0, 175.904, 0.0]*ufac"kJ/mol")

    # equilibrium constants Ki
    Ki_ref::Dict{Symbol,Float64} = Dict( rnames .=> [10.18, 1.947e-3, 1.913e-4])
    TKi_ref::Dict{Symbol,Float64} = Dict( rnames .=> [693.0, 693.0, 693.0]*ufac"K")
    # reaction enthalpies
    ΔHi::Dict{Symbol,Float64} = Dict( rnames .=> [-37.92, 182.09, 220.01]*ufac"kJ/mol")
	

    # values of gas phase species adsorption coefficients @ Tref
	Kj_ref::Dict{Symbol,Float64} = Dict( gnames .=> [0.0, 0.0296, 0.0, 0.0, 0.0, 0.0])
    TKj_ref::Dict{Symbol,Float64} = Dict( gnames .=> [648.0, 648.0, 823.0, 823.0, 823.0, 823.0]*ufac"K" )
    ΔHj::Dict{Symbol,Float64} = Dict( gnames .=> [-70.65, -82.9, -38.28, 88.68, 0, 0]*ufac"kJ/mol")

    KinData(ng,gnames,Fluids,gn,gni,nr,rnames,rn,rni,nuij,ki_ref,Tki_ref,Ei,Ki_ref,TKi_ref,ΔHi,Kj_ref,TKj_ref,ΔHj) = 
    new{nr}(ng,gnames,Fluids,gn,gni,nr,rnames,rn,rni,nuij,ki_ref,Tki_ref,Ei,Ki_ref,TKi_ref,ΔHi,Kj_ref,TKj_ref,ΔHj)
end

# !!!ALLOC Method to be called instead of data.ng
nreac(::KinData{NR}) where NR = NR


# !!!ALLOC Additional constructo taking nr as parameter	
# KinData(;nr=3, kwargs...) = KinData{nr}(;kwargs...)
# KinData(;kwargs...,) = begin
#     nr=(;kwargs...).nr
#     KinData{nr}(;kwargs...)
# end


#S3P = KinData{:S3P}()
const S3P = KinData()




## Xu & Froment 1989 kinetics
# different ordering of the 3 reactions compared to S3P kinetics:
# R1: CH4 + H2O = CO + 3 H2
# R2: CO + H2O = CO2 + H2
# R3: CH4 + 2 H2O = CO2 + 4 H2

# temporary variables
# ng = 6; gnames = [:CO, :H2, :CH4, :H2O, :CO2, :N2]; gn = Dict(1:ng .=> gnames); nr = 3; rnames = [:R1, :R2, :R3]; rn = Dict(1:nr .=> rnames);

const XuFroment = KinData{}(;
    # ng = ng,
    # gnames = [:CO, :H2, :CH4, :H2O, :CO2, :N2],
    # Fluids = [CO,H2,CH4,H2O,CO2,N2],
    # gn = gn,
    # gni = Dict(value => key for (key, value) in gn),    
    # nr = nr,
    # rnames = [:R1, :R2, :R3],
    # rn = rn,
    # rni = Dict(value => key for (key, value) in rn),
    nuij = vcat(
        [1, 3, -1, -1, 0, 0], #R1 : CH4 + H2O -> 3 H2 + CO
        [-1, 1, 0, -1, 1, 0], #R2 : CO + H2O -> CO2 + H2
        [0, 4, -1, -2, 1, 0], #R3 : CH4 + 2 H2O -> 4 H2 + CO2
        
    ),
    # # values of reaction rate constants @ Tref
    ki_ref = Dict( [:R1, :R2, :R3] .=>  log.(1.225*[1.842e-4, 7.558, 2.193e-5]) ),
    Tki_ref = Dict( [:R1, :R2, :R3] .=> [648.0, 648.0, 648.0]*ufac"K"),
    # # Activation energies
    Ei = Dict( [:R1, :R2, :R3] .=> [240.1, 67.13, 243.9]*ufac"kJ/mol"),
    # # equilibrium constants Ki
    Ki_ref = Dict( [:R1, :R2, :R3] .=> [1.913e-4, 10.18, 1.947e-3]),
    TKi_ref = Dict( [:R1, :R2, :R3] .=> [693.0, 693.0, 693.0]*ufac"K"),
    # # reaction enthalpies
    ΔHi = Dict( [:R1, :R2, :R3] .=> [220.01, -37.92, 182.09]*ufac"kJ/mol"),
    # values of gas phase species adsorption coefficients @ Tref
    Kj_ref = Dict( [:CO, :H2, :CH4, :H2O, :CO2, :N2] .=> [40.91, 0.0296, 0.1791, 0.4152, 0.0, 0.0]),
    TKj_ref = Dict( [:CO, :H2, :CH4, :H2O, :CO2, :N2] .=> [648.0, 648.0, 823.0, 823.0, 823.0, 823.0]*ufac"K" ),
    ΔHj = Dict( [:CO, :H2, :CH4, :H2O, :CO2, :N2] .=> [-70.65, -82.9, -38.28, 88.68, 0.0, 0.0]*ufac"kJ/mol")

)

const Hla_WGS = KinData{}(;
    # ng = ng,
    # gnames = gnames,
    # Fluids = [CO,H2,CH4,H2O,CO2,N2],
    # gn = gn,
    # gni = Dict(value => key for (key, value) in gn),    
    nr = 1,
    rnames = [:R1],
    rn = Dict(1:1 .=> [:R1]),
    rni = Dict(value => key for (key, value) in Dict(1:1 .=> [:R1])),
    nuij = vcat(
        [-1, 1, 0, -1, 1, 0], #R1 : CO + H2O -> CO2 + H2
    ),
    # # values of reaction rate constants
    ki_ref = Dict( [:R1] .=>  log.([700.0]) ),
    Tki_ref = Dict( [:R1] .=> [Inf]*ufac"K"),

    # # Activation energies
    Ei = Dict( [:R1] .=> [111.0]*ufac"kJ/mol"),
    # # reaction enthalpies
    ΔHi = Dict( [:R1] .=> [-37.92]*ufac"kJ/mol"),
    # # equilibrium constants Ki
    Ki_ref = Dict( [:R1] .=> [10.18]),
    TKi_ref = Dict( [:R1] .=> [693.0]*ufac"K"),
)

const Riedel_rWGS = KinData{}(;
    # ng = ng,
    # gnames = gnames,
    # Fluids = [CO,H2,CH4,H2O,CO2,N2],
    # gn = gn,
    # gni = Dict(value => key for (key, value) in gn),    
    nr = 1,
    rnames = [:R1],
    rn = Dict(1:1 .=> [:R1]),
    rni = Dict(value => key for (key, value) in Dict(1:1 .=> [:R1])),
    nuij = vcat(
        [1, -1, 0, 1, -1, 0], #R1 : CO2 + H2 -> CO + H2O
    ),
    # # values of reaction rate constants
    ki_ref = Dict( [:R1] .=>  log.([1.51e7]) ),
    Tki_ref = Dict( [:R1] .=> [Inf]*ufac"K"),

    # # Activation energies
    Ei = Dict( [:R1] .=> [55.0]*ufac"kJ/mol"),
    # # reaction enthalpies
    ΔHi = Dict( [:R1] .=> [-37.92]*ufac"kJ/mol"),
    # # equilibrium constants Ki
    Ki_ref = Dict( [:R1] .=> [10.18]),
    TKi_ref = Dict( [:R1] .=> [693.0]*ufac"K"),
)

const Wolf_rWGS = KinData{}(;
    # ng = ng,
    # gnames = gnames,
    # Fluids = [CO,H2,CH4,H2O,CO2,N2],
    # gn = gn,
    # gni = Dict(value => key for (key, value) in gn),    
    nr = 1,
    rnames = [:R1],
    rn = Dict(1:1 .=> [:R1]),
    rni = Dict(value => key for (key, value) in Dict(1:1 .=> [:R1])),
    nuij = vcat(
        [1, -1, 0, 1, -1, 0], #R1 : CO2 + H2 -> CO + H2O
    ),
    # # values of reaction rate constants
    ki_ref = Dict( [:R1] .=> log.([3100.0]) ),
    Tki_ref = Dict( [:R1] .=> [Inf]*ufac"K"),

    # # Activation energies
    Ei = Dict( [:R1] .=> [82.0]*ufac"kJ/mol"),
    # # reaction enthalpies
    ΔHi = Dict( [:R1] .=> [-37.92]*ufac"kJ/mol"),
    # # equilibrium constants Ki
    Ki_ref = Dict( [:R1] .=> [10.18]),
    TKi_ref = Dict( [:R1] .=> [693.0]*ufac"K"),
)


#
# By default, nreac(data) returns data.kinpar.nr. 
function nreac(data::Any)
    data.kinpar.nr
end



# # equilibrium constants Ki
# function Ki(kindata::AbstractKineticsData, T)
#     nr=kindata.nr
#     #Ki=zeros(Float64, nr)
#     Ki=zeros(typeof(T), nr)
#     Ki_ref = kindata.Ki_ref
#     TKi_ref = kindata.TKi_ref
#     ΔHi = kindata.ΔHi
#     rn = kindata.rn
#     for i=1:nr
#         Ki[i] = Ki_ref[rn[i]] * exp(-ΔHi[rn[i]]/ph"R"*(1.0/T - 1.0/TKi_ref[rn[i]]))
#     end
# 	Ki
# end
# equilibrium constants Ki
#function Ki(kindata, T)
function Ki(data, T)
    (;kinpar)=data
    #nr=kindata.nr
    #Ki=zeros(Float64, nr)
    # Ki=zeros(typeof(T), nr)
    Ki_ref = kinpar.Ki_ref
    TKi_ref = kinpar.TKi_ref
    ΔHi = kinpar.ΔHi
    rn = kinpar.rn
   
    # !!!ALLOC Use MVectors with static size information instef of Vector
    Ki=MVector{nreac(kinpar),eltype(T)}(undef)


    for i=1:nreac(kinpar)
        Ki[i] = Ki_ref[rn[i]] * exp(-ΔHi[rn[i]]/ph"R"*(1.0/T - 1.0/TKi_ref[rn[i]]))
    end
	Ki
end



# kinetic pre-factors ki
# function ki(kindata::AbstractKineticsData, T)
#     nr=kindata.nr
#     #ki=zeros(Float64, nr)
#     ki=zeros(typeof(T), nr)
#     ki_ref = kindata.ki_ref
#     Tki_ref = kindata.Tki_ref
#     Ei = kindata.Ei
#     rn = kindata.rn
#     for i=1:nr
#         ki[i] = exp(ki_ref[rn[i]]) * exp(-Ei[rn[i]]/ph"R"*(1.0/T - 1.0/Tki_ref[rn[i]]))
#         # @. ki_ref*exp(-1000*Ei/ph"R"*(1/T-1/Tki_ref))
#     end
# 	ki
# end
# kinetic pre-factors ki
#function ki(kindata::AbstractKineticsData, T)
function ki(data, T)

    (;kinpar)=data
    # !!!ALLOC Use MVectors with static size information instef of Vector
    ki=MVector{nreac(kinpar),eltype(T)}(undef)

    (;ki_ref,Tki_ref,Ei,rn)=kinpar
    
    for i=1:nreac(kinpar)
        ki[i] = exp(ki_ref[rn[i]]) * exp(-Ei[rn[i]]/ph"R"*(1.0/T - 1.0/Tki_ref[rn[i]]))
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




# # adsorption constants Kj
# function Kj(kindata::AbstractKineticsData, T)
#     ng=kindata.ng
#     # Kj=zeros(Float64, ng)
#     Kj=zeros(typeof(T), ng)

#     Kj_ref = kindata.Kj_ref
#     TKj_ref = kindata.TKj_ref
#     ΔHj = kindata.ΔHj
#     gn = kindata.gn

#     for j=1:ng
#         Kj[j] = Kj_ref[gn[j]] * exp(-ΔHj[gn[j]]/ph"R"*(1.0/T - 1.0/TKj_ref[gn[j]]))
        
#     end
# 	Kj
# end

# adsorption constants Kj
function Kj(data, T)

    (;kinpar)=data
    # !!!ALLOC Use MVectors with static size information instef of Vector
    Kj=MVector{ngas(data),eltype(T)}(undef)

    (;Kj_ref,TKj_ref,ΔHj,gn)=kinpar
    
    for j=1:ngas(data)
        Kj[j] = Kj_ref[gn[j]] * exp(-ΔHj[gn[j]]/ph"R"*(1.0/T - 1.0/TKj_ref[gn[j]]))
        
    end
	Kj
end



# function DEN(kindata::AbstractKineticsData,T,p)
# 	n=kindata.gni
#     # p_=p
#     p_=copy(p)
# 	p_[n[:H2O]]=p_[n[:H2O]]/p_[n[:H2]] # p_H2O/p_H2
# 	1+sum(Kj(kindata, T).*p_)
# end

function DEN(data,T,p)

    #  !!!ALLOC for types stubility & correctness
    #  !!!ALLOC initialize with zero(eltype) instead of 0.0
    DEN=zero(eltype(p))
    (;kinpar)=data
    (;gni)=kinpar
	
    p_=MVector{ngas(data),eltype(T)}(undef)
    p_ .= p

    @inbounds p_[gni[:H2O]]=p_[gni[:H2O]]/p_[gni[:H2]]    

    DEN = @inline 1 + sum(  Kj(data, T).*p_)
end

# equilibrium constant for WGS reaction
# function Kequil_WGS_Twigg(T)
#     K=zero(eltype(T))
#     Z=zero(eltype(T))

# 	Z = 1000.0/T-1.0
# 	K = exp(Z*(Z*(0.63508-0.29353*Z)+4.1778)+0.31688)	
# end

function K_rWGS_Twigg(T)
    Z=zero(eltype(T))
    Z = 1000.0/T
	1.0/exp(Z*(Z*(0.63508-0.29353*Z)+4.1778)+0.31688)	
end

function Kequil_WGS_Zimmermann(T)
	exp(2073.0/T-2.029)
end

function ri(data,T,p)
    
    (;kinpar)=data
    n=kinpar.gni
    r=kinpar.rni
    
    ri_ = MVector{nreac(kinpar),eltype(p)}(undef)
    Ki_ = @inline Ki(data, T)
    pexp = MVector{nreac(kinpar),eltype(p)}(undef)
    unitc=one(eltype(p))

    #p=MVector{ngas(data),eltype(T)}(undef)

    if kinpar == XuFroment
        # supply partial pressures in bar
        #@show p .= p_/ufac"bar"
        # returns (forward)WGS reation rate in mol g−1 hr−1
        unitc *=ufac"mol/(hr*g)"

        pexp[1] = @inbounds @inline 1/p[n[:H2]]^2.5*(p[n[:CH4]]*p[n[:H2O]]-p[n[:H2]]^3*p[n[:CO]]/Ki_[r[:R1]]) / DEN(data,T,p)^2 *unitc
        pexp[2] = @inbounds @inline 1/p[n[:H2]]*(p[n[:CO]]*p[n[:H2O]]-p[n[:H2]]*p[n[:CO2]]/Ki_[r[:R2]]) / DEN(data,T,p)^2 *unitc
        pexp[3] = @inbounds @inline 1/p[n[:H2]]^3.5*(p[n[:CH4]]*p[n[:H2O]]^2-p[n[:H2]]^4*p[n[:CO2]]/Ki_[r[:R3]]) / DEN(data,T,p)^2 *unitc

    elseif kinpar == S3P
        # supply partial pressures in bar
        #p .= p_/ufac"bar"

        # returns (forward)WGS reation rate in mol kg−1 s−1
        pexp[1] = @inbounds @inline 1/p[n[:H2]]*(p[n[:CO]]*p[n[:H2O]]-p[n[:H2]]*p[n[:CO2]]/Ki_[r[:R1]]) / DEN(data,T,p)^2
        pexp[2] = @inbounds @inline 1/p[n[:H2]]^3.5*(p[n[:CH4]]*p[n[:H2O]]^2-p[n[:H2]]^4*p[n[:CO2]]/Ki_[r[:R2]]) / DEN(data,T,p)^2
        pexp[3] = @inbounds @inline 1/p[n[:H2]]^2.5*(p[n[:CH4]]*p[n[:H2O]]-p[n[:H2]]^3*p[n[:CO]]/Ki_[r[:R3]]) / DEN(data,T,p)^2

    elseif kinpar == Riedel_rWGS
        #p .= p_/ufac"MPa"

        unitc *=ufac"mol/(s*g)"
        pexp[1] = @inbounds @inline (p[n[:CO2]]*p[n[:H2]] - p[n[:CO]]*p[n[:H2O]]*Kequil_WGS_Zimmermann(T)) / (p[n[:CO]] + 65.0*p[n[:H2O]] + 7.4*p[n[:CO2]]) *unitc

    elseif kinpar == Wolf_rWGS
        # convert partial pressures (Pa) into molar concentrations for kinetics
        c=MVector{ngas(data),eltype(T)}(undef)
        c .= p ./ (ph"R"*T)

        pexp[1] = @inbounds @inline (c[n[:CO2]]*c[n[:H2]]^0.3 - c[n[:CO]]*c[n[:H2O]]/c[n[:H2]]^0.7/K_rWGS_Twigg(T))
        
    end
    ri_ = @inline ki(data,T) .* pexp
 
end

# function ri(kindata::AbstractKineticsData,T,p)
#     kindata.RR(kindata,T,p)
# end

# function rr(i,kindata::AbstractKineticsData,T,p)
#     kindata.RR(i,kindata,T,p)
# end