module FixedBed

using VoronoiFVM, ExtendableGrids
using LessUnitful, Printf, Dates, CSV, DataFrames, Interpolations
using StaticArrays
using Documenter

include("physprops.jl")
export AbstractFluidProps, FluidProps, AbstractPropsCoeffs, PropsCoeffs
export dynvisc_gas, thermcond_gas, heatcap_gas, density_idealgas, binary_diff_coeff_gas, enthalpy_gas, enthalpy_gas_thermal
export dynvisc_mix, heatcap_mix, molarweight_mix, dynvisc_thermcond_mix, enthalpy_mix
export Air, N2, Ar, H2, CO2, CO, H2O, CH4
# export ngas


include("kinetics.jl")
export KinData, ri, S3P, XuFroment, Riedel_rWGS, Wolf_rWGS
export nreac


include("modelprops.jl")
export AbstractModelData, RePrPe, kbed, kbed_VDI_flattening, lambda_eff_AC, hsf, DK_eff
export M_matrix!, D_matrix!, DarcyVelo, MoleFrac!, MassFrac!, sel12by12, SurfaceOpticalProps, ReactorData, ngas, bareas
export grid_boundaries_regions, init_system

include("postprocess.jl")
export HeatFluxes_EB_I, BoundaryFlows_Integrals
end
