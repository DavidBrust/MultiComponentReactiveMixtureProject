module Sensitivity

using DataFrames, CSV
include("../../notebooks/model_props_physics/PCReactor_sens.jl")

function run(G_lamps, nflowins; SensPar=SensPar)

    function Tcenter(par;dim=2)
        solt,grid,sys,data=PCR_base(dim,par);
        sol = solt(solt.t[end])
        (;iT) = data
        if dim == 2
            Tc = sol[iT,91]
        elseif dim == 3
            Tc = sol[iT,3035]
        end
        Tc -= 273.15
        return Tc
    end

    function SensitivityUncertainty(par;dim=2)
        dT2_dp = T2_part_deriv(dim, par, ReactorData())
        return calc_combined_uncertainty_T2(par,dT2_dp)
    end

    Tcs = Float64[]
    Ucs = Float64[]
    for G_lamp in G_lamps
        par = transformSensPar([G_lamp, nflowins[1], SensPar[3:end]...])
        push!(Tcs, Tcenter(par;dim=3))
        push!(Ucs, SensitivityUncertainty(par)        )
    end
    
    df = DataFrame()
    df[!, :Glamp]  = G_lamps
    df[!, :nflowin]  .= repeat([nflowins[1]], length(G_lamps))
    df[!, :Tcenter] .= Tcs
    df[!, :Ucombined] .= Ucs

    CSV.write("data/out/2024-01-24/Tc_Uc.csv", df)

end


end