module FixedBed

using LessUnitful
using StaticArrays


include("physprops.jl")
export AbstractFluidProps, FluidProps, AbstractPropsCoeffs, PropsCoeffs
export dynvisc_gas, thermcond_gas, heatcap_gas, density_idealgas, binary_diff_coeff_gas, enthalpy_gas
export dynvisc_mix, heatcap_mix, molarweight_mix, dynvisc_thermcond_mix, enthalpy_mix
export Air, N2, Ar, H2, CO2, CO, H2O, CH4
export ngas


include("modelprops.jl")
#export AbstractModelData, ModelData, RePrPe, kbed, hsf, DK_eff
export AbstractModelData, RePrPe, kbed, kbed_VDI_flattening, lambda_eff_AC, hsf, DK_eff


include("kinetics.jl")
export AbstractKineticsData, ri, rr, S3P, S3P_, XuFroment1989, XuFroment
export nreac
end
