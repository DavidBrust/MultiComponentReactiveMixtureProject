module Sensitivity

using DataFrames, CSV, LinearAlgebra
include("../../notebooks/model_props_physics/PCReactor_sens.jl")

function run(nom_fluxs, nflowins; SensPar=SensPar)

    function getIndices(grid)
        d1 = 0.025*sqrt(2)/2
        d2 = 0.05*sqrt(2)/2
        
        # uc
        T_03 = [0.08,0.08,0.005] # y,x,z
        T_04 = T_03 .+ [d2, -d2, 0]
        T_05 = T_03 .+ [d1, d1, 0]
        T_06 = T_03 .+ [-d1, -d1, 0]
        T_07 = T_03 .+ [-d2, d2, 0]
        # lc
        T_12 = [0.08,0.08,0.0] # y,x,z
        T_13 = T_12 .+ [d1, d1, 0]
        T_14 = T_12 .+ [-d2, d2, 0]	
        
        
        function index(c, coords)
            cols = eachcol(coords)
            findmin(norm.([c - col for col in cols]))[2]
        end

        [index(c,grid[Coordinates]) for c in [T_03,T_04,T_05,T_06,T_07,T_12,T_13,T_14]]
    end

    function probe_Temps(par;dim=2)
        solt,grid,sys,data=PCR_base(dim,par;times=[0.0,10.0]);
        sol = solt(solt.t[end])
        (;iT) = data

        if dim == 2
            return sol[iT,91] - 273.15
        elseif dim == 3
            # return sol[iT,3035] .- 273.15
            return sol[iT,getIndices(grid)] .- 273.15
        end
        
    end

    function SensitivityUncertainty(par;dim=2)
        dT2_dp = T2_part_deriv(dim, par, ReactorData())
        return calc_combined_uncertainty_T2(par,dT2_dp)
    end

    T_03s = Float64[]
    T_04s = Float64[]
    T_05s = Float64[]
    T_06s = Float64[]
    T_07s = Float64[]

    T_12s = Float64[]
    T_13s = Float64[]
    T_14s = Float64[]
    T_03_Ucs = Float64[]
    for nom_flux in nom_fluxs
        par = transformSensPar([nom_flux, nflowins[1], SensPar[3:end]...])
        T_03,T_04,T_05,T_06,T_07,T_12,T_13,T_14 = probe_Temps(par;dim=3)
        push!(T_03s, T_03)
        push!(T_04s, T_04)
        push!(T_05s, T_05)
        push!(T_06s, T_06)
        push!(T_07s, T_07)
        push!(T_12s, T_12)
        push!(T_13s, T_13)
        push!(T_14s, T_14)
        push!(T_03_Ucs, SensitivityUncertainty(par)        )
    end
    
    df = DataFrame()
    df[!, :nom_flux]  = nom_fluxs
    df[!, :nflowin]  .= repeat([nflowins[1]], length(nom_fluxs))
    df[!, :T_03] .= T_03s
    df[!, :T_04] .= T_04s
    df[!, :T_05] .= T_05s
    df[!, :T_06] .= T_06s
    df[!, :T_07] .= T_07s
    df[!, :T_12] .= T_12s
    df[!, :T_13] .= T_13s
    df[!, :T_14] .= T_14s
    df[!, :T_03_Uc] .= T_03_Ucs

    CSV.write("data/out/2024-01-24/Tc_Uc.csv", df)

end


end