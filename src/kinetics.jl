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



abstract type AbstractKineticsData end


#Base.@kwdef struct KinData{NR} <:AbstractKineticsData
Base.@kwdef struct KinData{NR}
    ng::Int64=6 # number of gas phase species
    gnames::Vector{Symbol} = [:CO, :H2, :CH4, :H2O, :CO2, :N2]
    Fluids::Vector{FluidProps} = [CO,H2,CH4,H2O,CO2,N2]
    gn::Dict{Int, Symbol} 	= Dict(1:ng .=> gnames)
	# inverse names and fluid indices
	gni::Dict{Symbol, Int}  = Dict(value => key for (key, value) in gn)
    
    nr::Int64=3 # number of reactions
    rnames::Vector{Symbol} = [:R1, :R2, :R3]
    rn::Dict{Int, Symbol} 	= Dict(1:nr .=> rnames[1:nr])
	rni::Dict{Symbol, Int}  = Dict(value => key for (key, value) in rn)

    nuij::Vector{Int64} = vcat(
        [-1, 1, 0, -1, 1, 0], #R1 : CO + H2O -> CO2 + H2
        [0, 4, -1, -2, 1, 0], #R2 : CH4 + 2 H2O -> 4 H2 + CO2
        [1, 3, -1, -1, 0, 0], #R2 : CH4 + H2O -> 3 H2 + CO
    )
    # values of reaction rate constants @ Tref
	# ki_ref::Dict{Symbol,Float64} = Dict(:R1=>-1.51715, :R2=>-11.2379, :R3=>-Inf)
    ki_ref::Dict{Symbol,Float64} = Dict(rnames .=> [-1.51715, -11.2379, -Inf][1:nr])
    Tki_ref::Dict{Symbol,Float64} = Dict( rnames .=> [648.0, 648.0, 648.0][1:nr]*ufac"K")
    # Activation energies
    Ei::Dict{Symbol,Float64} = Dict( rnames .=> [0.0, 175.904, 0.0][1:nr]*ufac"kJ/mol")

    # equilibrium constants Ki
    Ki_ref::Dict{Symbol,Float64} = Dict( rnames .=> [10.18, 1.947e-3, 1.913e-4][1:nr])
    TKi_ref::Dict{Symbol,Float64} = Dict( rnames .=> [693.0, 693.0, 693.0][1:nr]*ufac"K")
    # reaction enthalpies
    ΔHi::Dict{Symbol,Float64} = Dict( rnames .=> [-37.92, 182.09, 220.01][1:nr]*ufac"kJ/mol")
	
    # values of gas phase species adsorption coefficients @ Tref
	Kj_ref::Dict{Symbol,Float64} = Dict( gnames .=> [0.0, 0.0296, 0.0, 0.0, 0.0, 0.0][1:ng])
    TKj_ref::Dict{Symbol,Float64} = Dict( gnames .=> [648.0, 648.0, 823.0, 823.0, 823.0, 823.0][1:ng]*ufac"K" )
    ΔHj::Dict{Symbol,Float64} = Dict( gnames .=> [-70.65, -82.9, -38.28, 88.68, 0, 0][1:ng]*ufac"kJ/mol")

    KinData(ng,gnames,Fluids,gn,gni,nr,rnames,rn,rni,nuij,ki_ref,Tki_ref,Ei,Ki_ref,TKi_ref,ΔHi,Kj_ref,TKj_ref,ΔHj) = 
    new{nr}(ng,gnames,Fluids,gn,gni,nr,rnames,rn,rni,nuij,ki_ref,Tki_ref,Ei,Ki_ref,TKi_ref,ΔHi,Kj_ref,TKj_ref,ΔHj)
end

# !!!ALLOC Method to be called instead of data.ng
nreac(::KinData{NR}) where NR = NR

# !!!ALLOC Additional constructo taking nr as parameter	


#S3P = KinData{:S3P}()
const S3P = KinData{}()


#Species indices:
#1) CO
#2) H2
#3) CH4
#4) H2O
#5) CO2
#6) N2

## Xu & Froment 1989 kinetics
# different ordering of the 3 reactions compared to S3P kinetics:
# R1: CH4 + H2O = CO + 3 H2
# R2: CO + H2O = CO2 + H2
# R3: CH4 + 2 H2O = CO2 + 4 H2

# temporary variables
ng = 6; gnames = [:CO, :H2, :CH4, :H2O, :CO2, :N2]; gn = Dict(1:ng .=> gnames); nr = 3; rnames = [:R1, :R2, :R3]; rn = Dict(1:nr .=> rnames);

const XuFroment = KinData{}(;
    ng = ng,
    gnames = gnames,
    Fluids = [CO,H2,CH4,H2O,CO2,N2],
    gn = gn,
    gni = Dict(value => key for (key, value) in gn),    
    nr = nr,
    rnames = rnames,
    rn = rn,
    rni = Dict(value => key for (key, value) in rn),

    nuij = vcat(
        [1, 3, -1, -1, 0, 0], #R1 : CH4 + H2O -> 3 H2 + CO
        [-1, 1, 0, -1, 1, 0], #R2 : CO + H2O -> CO2 + H2
        [0, 4, -1, -2, 1, 0], #R3 : CH4 + 2 H2O -> 4 H2 + CO2
        
    ),
    # # values of reaction rate constants @ Tref
    ki_ref = Dict( rnames .=>  log.(1.225*[1.842e-4, 7.558, 2.193e-5]) ),
    Tki_ref = Dict( rnames .=> [648.0, 648.0, 648.0]*ufac"K"),
    # # Activation energies
    Ei = Dict( rnames .=> [240.1, 67.13, 243.9]*ufac"kJ/mol"),
    # # equilibrium constants Ki
    Ki_ref = Dict( rnames .=> [1.913e-4, 10.18, 1.947e-3]),
    TKi_ref = Dict( rnames .=> [693.0, 693.0, 693.0]*ufac"K"),
    # # reaction enthalpies
    ΔHi = Dict( rnames .=> [220.01, -37.92, 182.09]*ufac"kJ/mol"),
    # values of gas phase species adsorption coefficients @ Tref
    Kj_ref = Dict( gnames .=> [40.91, 0.0296, 0.1791, 0.4152, 0.0, 0.0]),
    TKj_ref = Dict( gnames .=> [648.0, 648.0, 823.0, 823.0, 823.0, 823.0]*ufac"K" ),
    ΔHj = Dict( gnames .=> [-70.65, -82.9, -38.28, 88.68, 0.0, 0.0]*ufac"kJ/mol")
)


const Wolf_rWGS = KinData{}(;
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

const Riedel_rWGS = KinData{}(;
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


# By default, nreac(data) returns data.kinpar.nr. 
function nreac(data::Any)
    data.kinpar.nr
end



# equilibrium constants Ki
function Ki(data, T)
    (;kinpar)=data
    (;Ki_ref,TKi_ref,ΔHi,rn)=kinpar
   
    # !!!ALLOC Use MVectors with static size information instef of Vector
    Ki=MVector{nreac(kinpar),eltype(T)}(undef)

    for i=1:nreac(kinpar)
        Ki[i] = Ki_ref[rn[i]] * exp(-ΔHi[rn[i]]/ph"R"*(1.0/T - 1.0/TKi_ref[rn[i]]))
    end
	Ki
end



# kinetic pre-factors ki
function ki(data, T)
    (;kinpar)=data
    (;ki_ref,Tki_ref,Ei,rn)=kinpar

    # !!!ALLOC Use MVectors with static size information instef of Vector
    ki=MVector{nreac(kinpar),eltype(T)}(undef)

    for i=1:nreac(kinpar)
        ki[i] = exp(ki_ref[rn[i]]) * exp(-Ei[rn[i]]/ph"R"*(1.0/T - 1.0/Tki_ref[rn[i]]))
    end
	ki
end

# adsorption constants Kj
function Kj(data, T)
    (;kinpar)=data
    (;Kj_ref,TKj_ref,ΔHj,gn)=kinpar

    # !!!ALLOC Use MVectors with static size information instef of Vector
    Kj=MVector{ngas(data),eltype(T)}(undef)
    
    for j=1:ngas(data)
        Kj[j] = Kj_ref[gn[j]] * exp(-ΔHj[gn[j]]/ph"R"*(1.0/T - 1.0/TKj_ref[gn[j]]))        
    end
	Kj
end


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

function K_rWGS_Twigg(T)
    Z=zero(eltype(T))
    Z = 1000.0/T -1
	1.0/exp(Z*(Z*(0.63508-0.29353*Z)+4.1778)+0.31688)	
end

function Kequil_WGS_Zimmermann(T)
	exp(2073.0/T-2.029)
end

function ri(data,T,p_)
    
    (;kinpar)=data
    n=kinpar.gni
    r=kinpar.rni
    p = MVector{ngas(data),eltype(p_)}(undef)
    p .= p_

    ri_ = MVector{nreac(kinpar),eltype(p_)}(undef)
    
    Ki_ = @inline Ki(data, T)
    pexp = MVector{nreac(kinpar),eltype(p_)}(undef)
    pexp .= zero(eltype(p_))
    unitc=  one(eltype(p_))
    atol = one(eltype(p_))*1.0e-12
    if kinpar == XuFroment
        p ./= ufac"bar"
        unitc *=ufac"mol/(hr*g)"
        if @inbounds @views abs(p[n[:H2]]) > atol
            pexp[1] = @inbounds 1/p[n[:H2]]^2.5*(p[n[:CH4]]*p[n[:H2O]]-p[n[:H2]]^3*p[n[:CO]]/Ki_[r[:R1]])    
            pexp[2] = @inbounds 1/p[n[:H2]]*(p[n[:CO]]*p[n[:H2O]]-p[n[:H2]]*p[n[:CO2]]/Ki_[r[:R2]])
            pexp[3] = @inbounds 1/p[n[:H2]]^3.5*(p[n[:CH4]]*p[n[:H2O]]^2-p[n[:H2]]^4*p[n[:CO2]]/Ki_[r[:R3]])
            pexp ./= @inline DEN(data,T,p)^2
        end
        
    elseif kinpar == S3P
        p ./= ufac"bar"
        unitc *=ufac"mol/(s*kg)"
        pexp[1] = @inbounds 1/p[n[:H2]]*(p[n[:CO]]*p[n[:H2O]]-p[n[:H2]]*p[n[:CO2]]/Ki_[r[:R1]])    
        pexp[2] = @inbounds 1/p[n[:H2]]^3.5*(p[n[:CH4]]*p[n[:H2O]]^2-p[n[:H2]]^4*p[n[:CO2]]/Ki_[r[:R2]])
        pexp[3] = @inbounds 1/p[n[:H2]]^2.5*(p[n[:CH4]]*p[n[:H2O]]-p[n[:H2]]^3*p[n[:CO]]/Ki_[r[:R3]])
        pexp ./= @inline DEN(data,T,p)^2
    
    elseif kinpar == Wolf_rWGS
        # convert partial pressures (Pa) into molar concentrations (mol/m3) for Wolf 2016 kinetics
        c=MVector{ngas(data),eltype(p)}(undef)
        c .= p ./ (ph"R"*T)
        unitc *=ufac"mol/(s*kg)"

        pexp[1] = @inbounds @inline (c[n[:CO2]]*c[n[:H2]]^0.3 - c[n[:CO]]*c[n[:H2O]]/c[n[:H2]]^0.7/K_rWGS_Twigg(T))

    elseif kinpar == Riedel_rWGS
        # convert partial pressures (Pa) into (MPa) for Riedel 2001 kinetics
        p ./= ufac"MPa"
        unitc *=ufac"mol/(s*g)"
        pexp[1] = @inbounds @inline (p[n[:CO2]]*p[n[:H2]] - p[n[:CO]]*p[n[:H2O]]*Kequil_WGS_Zimmermann(T)) / (p[n[:CO]] + 65.0*p[n[:H2O]] + 7.4*p[n[:CO2]])
    else # ad-hoc defined kinetics, apply mass action law
        (;nuij) = kinpar
        c=MVector{ngas(data),eltype(p)}(undef)
        c .= p ./ (ph"R"*T)

        @inbounds for i=1:nreac(kinpar)
            pexp[i] = zero(eltype(pexp))
            for j=1:ngas(data) 
                if nuij[j,i] > 0 # product, contrib of backward reaction
                    pexp[i] -= c[j]^nuij[j,i]/Ki_[j]
                elseif nuij[j,i] < 0 # reactant, contrib of forward reaction
                    pexp[i] += c[j]^(-nuij[j,i])
                end
            end
        end        
    end
    
    ri_ = @inline ki(data,T) .* pexp * unitc
end
