#Dynamic viscosity of gases at low pressures, Pa*s
function dynvisc_gas(data, T)
	(;A,B,C,D,E) = data.DynVisc
	# VDI heat atlas 2010 D3.1 Equation (3)
	A+B*T+C*T^2+D*T^3+E*T^4 * ufac"Pa*s"
end

#Thermal conductivity of gases at low pressures, W/(m*K)
function thermcond_gas(data, T)
	(;A,B,C,D,E) = data.ThermCond
	# VDI heat atlas 2010 D3.1 Equation (5)
	A+B*T+C*T^2+D*T^3+E*T^4 * ufac"W/(m*K)"
end

#Molar heat capacity of ideal gases, J/(kg*K)
function heatcap_gas(data, T)
	(;A,B,C,D,E,F,G) = data.HeatCap
	# VDI heat atlas 2010 D3.1 Equation (10)
	T_ApT = (T/(A+T))
	(B+(C-B)*T_ApT^2*(1- (A/(A+T))*(D+E*T_ApT+F*T_ApT^2+G*T_ApT^3) ) ) * ph"R" / data.MW  * ufac"J/(kg*K)"
end

function density_idealgas(data, T, p)
	p/(ph"R"*T)*data.MW*ufac"kg/m^3"
end


abstract type AbstractPropsCoeffs end

Base.@kwdef mutable struct PropsCoeffs <: AbstractPropsCoeffs
    A::Float64=1.0
	B::Float64=1.0
	C::Float64=1.0
	D::Float64=1.0
	E::Float64=1.0
	F::Float64=1.0
	G::Float64=1.0
end



abstract type AbstractFluidProps end

Base.@kwdef mutable struct FluidProps <: AbstractFluidProps
	name::String="Air"
	MW::Float64=28.96*ufac"g/mol"
	HeatCap::AbstractPropsCoeffs=PropsCoeffs(
	A=2548.9320,
	B=3.5248,
	C=-0.6366,
	D=-3.4281,
	E=49.8238,
	F=-120.3466,
	G=98.8658
	)
	ThermCond::AbstractPropsCoeffs=PropsCoeffs( 
	A=-0.908e-3,
	B=0.112e-3,
	C=-0.084333e-6,
	D=0.056964e-9,
	E=-0.015631e-12
	)
	DynVisc::AbstractPropsCoeffs=PropsCoeffs(
	A=-0.01702e-5,
	B=0.79965e-7,
	C=-0.72183e-10,
	D=0.04960e-12,
	E=-0.01388e-15
	)	
end

# all values taken from VDI heat atlas 2010 chapter D3.1
N2=FluidProps(
    name="N2",
	# Table 1
	MW=28.01*ufac"g/mol",
	# Table 6
	HeatCap=PropsCoeffs(
	A=432.2027,
	B=3.5160,
	C=2.8021,
	D=-4.1924,
	E=42.0153,
	F=-114.2500,
	G=111.1019
	),
	# Table 10
	ThermCond=PropsCoeffs(
	A=-0.133e-3,
	B=0.101e-3,
	C=-0.060650e-6,
	D=0.033610e-9,
	E=-0.0071e-12
	),
	# Table 8
	DynVisc= PropsCoeffs(
	A=-0.01020e-5,
	B=0.74785e-7,
	C=-0.59037e-10,
	D=0.03230e-12,
	E=-0.00673e-15
)
);

# all values taken from VDI heat atlas 2010 chapter D3.1
Air=FluidProps(
    name="N2",
	# Table 1
	MW=28.96*ufac"g/mol",
	# Table 6
	HeatCap=PropsCoeffs(
	A=2548.9320,
	B=3.5248,
	C=-0.6366,
	D=-3.4281,
	E=49.8238,
	F=-120.3466,
	G=98.8658
	),
	# Table 10
	ThermCond=PropsCoeffs( 
	A=-0.908e-3,
	B=0.112e-3,
	C=-0.084333e-6,
	D=0.056964e-9,
	E=-0.015631e-12
	),
	# Table 8
	DynVisc=PropsCoeffs(
	A=-0.01702e-5,
	B=0.79965e-7,
	C=-0.72183e-10,
	D=0.04960e-12,
	E=-0.01388e-15
)
);