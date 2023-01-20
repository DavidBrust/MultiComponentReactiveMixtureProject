module FixedBed

using LessUnitful

include("physprops.jl")
export AbstractFluidProps, FluidProps, AbstractPropsCoeffs, PropsCoeffs
export dynvisc_gas, thermcond_gas, heatcap_gas, density_idealgas
export Air, N2

include("kinetics.jl")
export ri
export S3P

end
