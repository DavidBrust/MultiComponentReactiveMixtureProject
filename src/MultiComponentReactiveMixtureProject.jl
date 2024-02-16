module MultiComponentReactiveMixtureProject

using VoronoiFVM, ExtendableGrids
using LessUnitful, Printf, Dates, CSV, DataFrames, Interpolations
using StaticArrays
using Documenter

include("physprops.jl")
export AbstractFluidProps, FluidProps, AbstractPropsCoeffs, PropsCoeffs
export dynvisc_gas, thermcond_gas, heatcap_gas, density_idealgas, binary_diff_coeff_gas, enthalpy_gas, enthalpy_gas_thermal
export dynvisc_mix, heatcap_mix, molarweight_mix, dynvisc_thermcond_mix, enthalpy_mix
export Air, N2, Ar, H2, CO2, CO, H2O, CH4

include("kinetics.jl")
export KinData, ri, S3P, XuFroment, Riedel_rWGS, Wolf_rWGS
export nreac

include("darcy_maxwell_stefan_model.jl")
export M_matrix!, D_matrix!, DarcyVelo, MoleFrac!, MassFrac!
export DMS_flux, DMS_reaction, DMS_storage, DMS_boutflow

include("photo_thermal_reactor.jl")
export PTR_top, PTR_side, PTR_bottom, PTR_bcond, PTR_bflux, PTR_bstorage, PTR_grid_boundaries_regions, PTR_init_system
export RePrPe, kbed, kbed_VDI_flattening, lambda_eff_AC, SurfaceOpticalProps, ReactorData, ngas, bareas


include("postprocess.jl")
export HeatFluxes_EB_I, BoundaryFlows_Integrals, Print_summary, Print_summary_ext, WriteSolution3D

end
