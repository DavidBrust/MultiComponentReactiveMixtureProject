### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ bcf12f50-b3a6-11ee-0f2b-1b63d8c56346
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,"../.."))
	using Revise
	using LessUnitful
	using DataFrames, CSV
	using Statistics, StatsPlots

	using VoronoiFVM
	using GridVisualize
	using ExtendableGrids
	using PlutoVista, Plots, Colors, PlutoUI
	import PlutoUI: combine
	using Roots
	using ForwardDiff, DiffResults
	using FixedBed

	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 640a4e80-b06e-40b2-a659-13feed18f55a
TableOfContents()

# ╔═╡ 07cf6a69-c1eb-49a6-adbb-fee87efb9519
md"""
# Optical Properties
## Thermocouple
$(LocalResource(
	"./img/Opt_props_TC.png",
)
)
https://www.osti.gov/servlets/purl/923075
## Catalyst surface
$(
LocalResource(
	"./img/Opt_props.png",
)
)
"""

# ╔═╡ 31d90ef7-e398-49ab-9338-45eb505b720e
md"""
# Temperature Calculation
"""

# ╔═╡ 6355e676-013a-441b-959f-089ffc3010a4
phi02 = phi20 = 1.0;

# ╔═╡ ed68f853-f81a-4b24-a25b-4ecf41fd493a
md"""
### Sensitivities
"""

# ╔═╡ 9b3a6003-7263-40ed-a11b-3a55c3a7a510
rel_deltas=[0.1,0.2,0.3]

# ╔═╡ caae0418-636e-4a68-8cea-0c3a70704fee
eps1_mo(T) = 0.52 + 2.83*1.0e-4 *(T-25.0) # hemispherical emissivity, moderate oxidation

# ╔═╡ 86b550c8-212b-4458-b701-2f150b87bcfd
eps1_ho(T) = 0.66 + 2.06*1.0e-4 *(T-25.0) # hemispherical emissivity, heavy oxidation

# ╔═╡ 3a0c71bd-339f-437f-9705-ddcc543c3e8d
function G2(T2,abs2,eps2,G_lamp; phi02=phi02)
	rho2 = 1-abs2
	G2_IR = eps2*ph"σ"*T2^4 
	G2_VIS = rho2*phi02*G_lamp
	G2_IR, G2_VIS
end

# ╔═╡ 8bd04d1e-08d6-495e-ad64-435a9b504c17
phi10 = phi12 = 0.5;

# ╔═╡ 850a86cf-97db-40dd-994b-e9dd903a0522
function EB_TC_conv(T1, T2, data; phi10=phi10, phi12=phi12,)
	(;G_lamp, T_gas_in, k_nat_conv) = data
	G2_IR, G2_VIS = G2(T2,data.uc_cat.alpha_vis,data.uc_cat.eps,G_lamp)
	abs1_VIS = data.uc_mask.alpha_vis
	eps1 = data.uc_mask.eps
	abs1_IR = eps1
	Tgas_avg = 0.5*(T2+T_gas_in)
	
	phi12*(abs1_VIS*G2_VIS+abs1_IR*G2_IR) + abs1_VIS*phi10*G_lamp - eps1*ph"σ"*T1^4 - k_nat_conv*(T1-Tgas_avg)
end

# ╔═╡ 1aefe230-88fc-4397-91d6-7899499cb3ea
T1_conv(T2,data) = find_zero(x -> EB_TC_conv(x, T2+273.15, data), 800) -273.15

# ╔═╡ 06fb7b1b-3a36-4cd9-8485-8c2c55eab86c
function transformSensPar(SensPar) 
	return [SensPar[1] + 273.15, SensPar[2] * ufac"kW/m^2", SensPar[3:end]...]
end

# ╔═╡ 4027a2c5-9534-4b71-afe5-b83e0cc809d6
function df_par_unc(SensPar)
	
	df = DataFrame(
		par_name=[:T2, :G0, :abs1, :eps1, :abs2, :eps2, :h_conv],
		mean_value=SensPar,
		uncertainty_rel=[0.0, 0.05, 0.05, 0.05, 0.10, 0.10, 0.05],
	)
	
	df[!, :uncertainty_interval] .= df.uncertainty_rel .* df.mean_value
	replace!(df.uncertainty_interval, 0.0 => 30.0)	
	df[!, :variance] .= df.uncertainty_interval.^2 ./ 3 # rectangular distribution
	df
end

# ╔═╡ 382425a8-93c6-4fff-be01-3d4c196cdb7a
md"""
### Combined Standard Uncertainty
Following calculation method in ISO/IEC Guide 98-3 (Guide to the Expression of Uncertainty in Measurement, GUM):
```math
u^2_{\text c} = \sum_{i=1}^N \left( \frac{\partial f}{\partial x_i} \right)^2 u^2(x_i)
```
Following procedure __Type B__ in determining the measurement uncertainty and assuming __rectangular distribution__ of estimate of measurand aroung the mean value, the variance $u^2(x_i)$ of an estimate $x_i$ is computed as $u^2(x_i)=a^2 /3$. Here $a$ is the symmetric interval of uncertainty of the measurement.
"""

# ╔═╡ 82315b1c-bfe3-48f3-a537-d5ef0cc2312c
function EB1_sens(T1, par, data;  phi10=phi10, phi12=phi12)
	T2, G_lamp, abs1, eps1, abs2, eps2, h_conv_D = par
	(;T_gas_in, X0,) = data
	
	G2_IR, G2_VIS = G2(T2,abs2,eps2,G_lamp)
	Tgas_avg = 0.5*(T2+T_gas_in)

	h_in = enthalpy_mix(data, T_gas_in, X0)
	h_out = enthalpy_mix(data, T2, X0)
	
	#abs1*phi12*G2_ + abs1*phi10*G_lamp - eps1*ph"σ"*T1^4 - h_conv_D*(T1-Tgas_avg)
	phi12*(abs1*G2_VIS+eps1*G2_IR) + abs1*phi10*G_lamp - eps1*ph"σ"*T1^4 - h_conv_D*(T1-Tgas_avg) # abs1_IR = eps1
end

# ╔═╡ c3c54b35-3121-4f33-9c09-8dd659bd4f84
md"""
# Parameter Variation
"""

# ╔═╡ 36a032e5-d0f4-49a2-98fb-e362189b0901
md"""
## Absorbtance (solar) of Inconel
Absorptance of different Ni/Cr materials after oxidation. 
__Y.S. Touloukian and D.P.DeWitt__, "Thermal Radiative Properties: Metallic Elements and Alloys," ThermophysicalProperties of Matter, vol. 7, 1970
Data and graph taken from pp. 1392-1394.

$(
LocalResource(
	"./data/Inconel_solar_reflectance_data_edit.png"
)
)

"""

# ╔═╡ 8dc9a71b-759d-4440-a8c4-030f6a9e57e4
abs1_mean = let 
	rho1_meas = [0.432, 0.361, 0.443, 0.390, 0.559, 0.370]
	1-sum(rho1_meas)/length(rho1_meas) # 1 = abs + rho, abs = 1 - rho
end

# ╔═╡ 9073c7da-cc86-4191-b8a8-b49379803c89
md"""
Mean values of main parameters. Calculate thermocouple hemispherical emissivity according to the proposed correlation for heavily oxidized Inconel 600 surfaces.

Operating Conditions:
-  $G_0$ = $(@bind G_lamp NumberField(40.0:1.0:100.0, default=70.0)) kW/m^2
-  $\dot n_{\text{in,total}}$ = $(@bind nflowin NumberField(1.0:0.2:20.0, default=7.4)) mol/hr
-  $h_{\text{conv}}$ = $(@bind h_conv_D NumberField(0.0:0.5:25.0, default=15)) W/m^2/K
-  $T_2$ = $(@bind T2 NumberField(400.0:1.0:800.0, default=590.0)) °C
Thermocouple optical properties:
-  $\alpha_1$ (vis) = $(@bind abs1 NumberField(0.1:0.01:0.9, default=abs1_mean))
-  $\epsilon_1$ (IR) = $(@bind eps1 NumberField(0.1:0.01:0.9, default=eps1_ho(650.0)))
Catalyst surface optical properties:
-  $\alpha_2$ (vis) = $(@bind abs2 NumberField(0.1:0.01:0.9, default=0.39))
-  $\epsilon_2$ (IR) = $(@bind eps2 NumberField(0.1:0.01:0.9, default=0.56))
"""

# ╔═╡ c5c9bc3b-d460-4730-9637-64fd0443860e
SensPar = [T2, G_lamp, abs1, eps1, abs2, eps2, h_conv_D]
#SensPar = (T2=T2, G_lamp=G_lamp, abs1=abs1, eps1=eps1, abs2=abs2, eps2=eps2, h_conv_D=h_conv_D)

# ╔═╡ 045f3f13-c606-4bf7-8416-2a95248ca35b
df_par_unc(transformSensPar(SensPar))

# ╔═╡ 5229e919-38eb-41e2-acc4-0355bb933528
md"""
```math
\Delta T = T_1 - T_2
```
"""

# ╔═╡ f7a6f0af-b34a-4fa1-ad1b-579eca0862c7
mr=maximum(rel_deltas);

# ╔═╡ 72ed6c5e-a7cc-4d58-aeff-4209621b6fb5
md"""
### One-at-a-time Sensitivity
Variation of one of the sensitivity parameters from - $(Integer(mr*100))% ... mean ... + $(Integer(mr*100))% while keeping the other parameters at their mean values.
"""

# ╔═╡ b7db1ccc-7b9d-4332-bff7-04d9e97b7138
md"""
### One-at-a-time Sensitivity
Variation of one of the sensitivity parameters from - $(Integer(mr*100))% ... mean ... + $(Integer(mr*100))% while keeping the other parameters at their mean values.
"""

# ╔═╡ e7ecc92c-226f-4b00-9a29-72abda41225a
md"""

Parameter ranges (min, mean, max): - $(Integer(mr*100))% ... mean ... + $(Integer(mr*100))%

|                        | Min | Mean | Max | Unit      |
|------------------------|-----|------|-----|-----------|
| Glamp                  | $(round(G_lamp*(1-mr),sigdigits=2))    | $(G_lamp)   |  $(round(G_lamp*(1+mr),sigdigits=2))   | kw/m^2    |
| Feed Flow [nflowin]             |  $(round(nflowin*(1-mr),sigdigits=2))   |  $(nflowin)    |  $(round(nflowin*(1+mr),sigdigits=2))   | mol/hr    |
| Convection Coeff. [h\_conv_D]     |  $(round(h_conv_D*(1-mr),sigdigits=2))   |   $(h_conv_D)   |   $(round(h_conv_D*(1+mr),sigdigits=2))  | W/(m^2*K) |
| TC Absorptance (Vis) [abs1] |  $(round(abs1*(1-mr),sigdigits=2))   |  $(abs1)    |  $(round(abs1*(1+mr),sigdigits=2))   | -         |
| TC Emissivity (IR) [eps1]  |  $(round(eps1*(1-mr),sigdigits=2))   | $(eps1) |   $(round(eps1*(1+mr),sigdigits=2))  | -         |
| Cat. Absorptance (Vis) [abs2] |  $(round(abs2*(1-mr),sigdigits=2))   |  $(abs2)    |  $(round(abs2*(1+mr),sigdigits=2))   | -         |
| Cat. Emissivity (IR) [eps2]  |  $(round(eps2*(1-mr),sigdigits=2))   | $(eps2) |   $(round(eps2*(1+mr),sigdigits=2))  | -         |
"""

# ╔═╡ 8fc0ccbd-8108-4957-8ecd-4a0559fc6eed
D_TC = 1*ufac"mm"

# ╔═╡ d7197aab-a552-47d0-8772-b647ff047d91
md"""
# Convection
Obtain convective heat transfer coefficient for convective heat transfer from the cylindrical thermocouple in the cross-flow to the surrounding gas from Nusselt correlation:

```math
\text{Nu}_D = \frac{\text{convective heat transfer}}{\text{conductive heat transfer}} = \frac{h D}{\lambda}
```
where $D$ is the thermocouple diameter with D= $(D_TC/ufac"mm") mm.
"""

# ╔═╡ 63ea26d4-8e31-4378-ba38-1d3f4132e4e1
md"""
## Nusselt
Nusselt correlation for convective heat transfer coefficient for a cylinder in cross-flow:
$(
LocalResource(
	"./img/Nussel_cylinder_crossflow.png"
)
)

(from WSUE Skript, RWTH)
"""

# ╔═╡ b9da0e38-6309-4469-9818-836f5bb77f52
Nu(Re_d,Pr) = 0.989 * Re_d^0.330 * Pr^0.4

# ╔═╡ 48d04afd-346f-423c-ad0c-94a0cabcdfc6
function h_conv(Nu_D,D_TC,data;T=data.T_gas_in,x=data.X0)
	ηf, λf = dynvisc_thermcond_mix(data, T, x)
	Nu_D * λf / D_TC
end

# ╔═╡ 7a968276-bd16-4b0f-86a8-edfb5c0bfb1a
md"""
## Reynolds
```math
\text{Re}_d = \frac{\rho u d}{\eta}
```
"""

# ╔═╡ e922cfd8-ea0f-416a-b5ef-3325a4905c3f
Across = 140^2*ufac"mm^2";

# ╔═╡ 19a9466c-60a7-4b10-9ce7-f3f8f167712c
function EB_catsurf_conv(T2, data; phi02=phi02, phi20=phi20, Across=Across)
	(;G_lamp, T_gas_in, nflowin, X0) = data
	
	abs2 = data.uc_cat.alpha_vis
	eps2 = data.uc_cat.eps

	h_in = enthalpy_mix(data, T_gas_in, X0)
	h_out = enthalpy_mix(data, T2, X0)
	
	abs2*phi02*G_lamp - nflowin/Across*(h_out-h_in) - eps2*ph"σ"*T2^4
end

# ╔═╡ e1ba6625-dda1-40cd-90f7-6d0604152edb
T2_enth(data::ReactorData) = find_zero(x -> EB_catsurf_conv(x+273.15, data), 800)

# ╔═╡ 854a6121-40a6-4d62-8846-02963ae3c136
function EB2_sens(T2, par, data; phi02=phi02, phi20=phi20, Across=Across)
	G0, abs2, eps2, nflowin = par
	(;T_gas_in, X0,) = data
	

	h_in = enthalpy_mix(data, T_gas_in, X0)
	h_out = enthalpy_mix(data, T2, X0)
	
	abs2*phi02*G0 - nflowin/Across*(h_out-h_in) - eps2*ph"σ"*T2^4
end

# ╔═╡ 2b98d5e6-1da8-46e1-8521-1ae8137e07a3
function flow_velo(Across, data; T=data.T_gas_in, p=data.p, x=data.X0)
	(;mflowin) = data
	ρf = density_idealgas(T, p, x, data)
	Vdot = mflowin/ρf
	u = Vdot / Across
end

# ╔═╡ 91cb28c4-ed83-4427-9494-65ea9e436ca9
function Re(d,u,data; p=data.p, x=data.X0)
	(;mflowin,T_gas_in) = data
	T = 0.5*(T_gas_in+T2_enth(data)+273.15)
	ρf = density_idealgas(T, p, x, data)
	ηf, λf = dynvisc_thermcond_mix(data, T, x)
	
	Re = u*ρf*d/ηf # Reynolds number
end

# ╔═╡ 654c6438-6e0e-4da6-81c8-d34acf96acb4
function Pr(data; p=data.p, x=data.X0)
	(;mflowin, T_gas_in, mmix0) = data
	T = 0.5*(T_gas_in+T2_enth(data)+273.15)
	ηf, λf = dynvisc_thermcond_mix(data, T, x)
	cf = heatcap_mix(data, T, x)
	
	Pr = cf/mmix0*ηf/λf # Prandtl number
end

# ╔═╡ c3a8a9c9-113e-4bcc-9ad7-0ceb045a123c
function newData(;abs1=abs1,eps1=eps1,G_lamp=G_lamp*ufac"kW/m^2",abs2=abs2,eps2=eps2,nflowin=nflowin*ufac"mol/hr",h_conv_D=h_conv_D,D_TC=1.0*ufac"mm",Across=140^2*ufac"mm^2")

	data = ReactorData(
		#X0 = X0,
		G_lamp = G_lamp,
		nflowin = nflowin,
		
		uc_mask = SurfaceOpticalProps( # thermocouple surface 1
			alpha_IR=eps1, 
			tau_IR=0.0, 
			alpha_vis=abs1, 
			tau_vis=0.0 
		),
			
		uc_cat = SurfaceOpticalProps( # catalyst surface 2
			alpha_IR=eps2, 
			tau_IR=0.0, 
			alpha_vis=abs2,
			tau_vis=0.0 
		)
	)

	if isnothing(h_conv_D)
		Re_d = Re(D_TC,flow_velo(Across, data),data)

		Nu_D = Nu(Re_d,Pr(data))
		
		data.k_nat_conv = h_conv(Nu_D,D_TC,data)
	else
		data.k_nat_conv = h_conv_D
	end
	data
end

# ╔═╡ 9cc7316f-2654-425d-bc9d-1e26197f081c
data = newData()

# ╔═╡ a56e991c-f6d7-4435-aaaa-6a1a80fe7fb6
md"""
## Catalyst Surface
Calculate equilibrium temperature of the catalyst surface resulting from lamp irradiation by solving steady-state energy balance around the catalyst surface (marked in red):
$(LocalResource(
	"./img/E_balance_cat_surf.png"
)
)

With the following parameters:
- Glamp = $(data.G_lamp/ufac"kW/m^2") kW/m^2
-  $\alpha_2$ = $(data.uc_cat.alpha_vis)
-  $\epsilon_2$ = $(data.uc_cat.eps)

Energy balance around catalyst surface including enthalpy change of gases that are heated up:
```math
\begin{align}
0 &= \dot Q_{\text{abs,0}} - \Delta \dot H_{\text{gas}} - \dot Q_{\text{emit}} \\
0 &= \alpha_2 \phi_{02} A_0 G_{\text{lamp}} - \dot m_{\text{gas}} [h_{\text{gas,out}}(T_2) - h_{\text{gas,in}}(T_{\text{gas,in}})] - A_2 \epsilon_2 \sigma T_2^4\\
\end{align}
```
Solve non-linear equation for $T_2$.
"""

# ╔═╡ 13fdff73-cb2b-4285-a43e-c81b3aca5c4a
md"""
- Catalyst surface temperature __with__ enthalpy change of gas: $\mathbf{T_2}$ = $(Integer(round(T2_enth(data)))) °C

"""

# ╔═╡ 4f24e591-73d4-4b34-b812-24d80b531879
f2(x,p) = EB2_sens(x+273.15, p, data)

# ╔═╡ 1fcdb9b3-71c6-4b05-8460-d3c86504aad3
begin
	# sensitivity parameter
	p2 = [G_lamp*ufac"kW/m^2", abs2, eps2, nflowin*ufac"mol/hr"] 
	x2_ = find_zero(f2, 800.0, Order1(), p2)
	fx2 = ForwardDiff.derivative(x -> f2(x, p2), x2_)
	fp2 = ForwardDiff.gradient(p -> f2(x2_, p), p2)
	x2_p = -fp2 / fx2
end

# ╔═╡ 533362a7-7e32-4df7-98db-07030004fef2
md"""
Formulation for sensitivity analysis: Calculate partial derivatives with ForwardDiff.jl.
```math
T_2 = T_2(G_0,\alpha_2,\epsilon_2,\dot n_{\text{feed,in}})
```
T2 = $(round(x2_,sigdigits=3)) °C

Sensitivities:
1.  $\partial T_2 / \partial G_0 =$  $(round(x2_p[1]*ufac"kW/m^2",sigdigits=2)) K/(kW/m^2)
1.  $\partial T_2 / \partial \alpha_2 =$  $(round(x2_p[2],sigdigits=2)) K
1.  $\partial T_2 / \partial \epsilon_2 =$  $(round(x2_p[3],sigdigits=2)) K
1.  $\partial T_2 / \partial \dot n_{\text{feed,in}} =$  $(round(x2_p[4]*ufac"mol/hr",sigdigits=2)) K/(mol/hr)
"""

# ╔═╡ 71e4bc55-058d-4c2d-9432-baae23a67ccb
T2_sens(par) = find_zero(x -> EB2_sens(x+273.15, par, data), 800)

# ╔═╡ bd4e7816-7790-4657-aa3c-78ce6bd578c7
md"""
## Thermocouple
Steady-state energy balance around the thermocouple including convective heat transfer (marked in red):
$(
LocalResource(
	"./img/E_balance_sketch_conv.png"
)
)

With the following parameters:
-  $G_0$ = $(data.G_lamp/ufac"kW/m^2") kW/m^2
-  $\alpha_2$ = $(data.uc_cat.alpha_vis)
-  $\epsilon_2$ = $(data.uc_cat.eps)
-  $\mathbf{T_2=}$ __$(Integer(round(T2))) °C__

-  $\alpha_1$ = $(data.uc_mask.alpha_vis)
-  $\epsilon_1$ = $(data.uc_mask.eps)
-  $\rho_2$ = $(data.uc_cat.rho_vis)

-  $\phi_{10}$ = $(phi10)
-  $\phi_{12}$ = $(phi12)
-  $\phi_{02}$ = $(phi02)

Convective heat transfer is calculated for the mean gas temperature in the vicinity of the catalyst surface $\bar T_{\text{gas}} = (T_2+T_{\text{gas,in}})/2$.

```math
\begin{align}
0 &= \dot Q_{\text{abs,2}} + \dot Q_{\text{abs,0}} - \dot Q_{\text{emit}} - \dot Q_{\text{conv}}\\
\end{align}
```
The  radiosity of the catalyst surface (2) consists of contributions in the visible (VIS) and infrared (IR) spectra. The VIS component results from the reflected light coming from the lamp, the IR component results from the thermal emission of the catalyst surface. The absorptivity of the thermocouple might behave differently in the VIS and IR spectral regions, accounted by $\alpha_1^{\text{IR}}$ and $\alpha_1^{\text{VIS}}$.

```math
\begin{align}

G_2 &= G_2^{\text{IR}} + G_2^{\text{VIS}} \\
&= \epsilon_2 \sigma T_2^4 + \rho_2 \phi_{02} G_0 \\
\end{align}
```

```math
\begin{align}
0 &= \phi_{12} A_1 (\alpha_1^{\text{IR}} G_2^{\text{IR}} + \alpha_1^{\text{VIS}} G_2^{\text{VIS}}) + \alpha_1 \phi_{10} A_1 G_0- \epsilon_1 A_1 \sigma T_1^4 - A_1 h_{\text{conv}}(T_1 - \bar T_{\text{gas}})\\
\end{align}
```
"""

# ╔═╡ 28249727-2511-4ec1-9067-e2b6fb0651f1
md"""
Temperature of Thermocouple with convection: $\mathbf{T_1}$ = $(Integer(round(T1_conv(T2,data)))) °C
"""

# ╔═╡ c6d4d255-29e1-4d68-9f7e-77afec7186d6
f1(x,p) = EB1_sens(x+273.15, p, data)

# ╔═╡ 1d6a9361-0473-488d-a56a-914e908b7a83
function calc_partial_derivs_T1(par)
	# sensitivity parameter
	#p1 = transformSensPar(SensPar)
	x1_ = find_zero(f1, 800.0, Order1(), par)
	fx1 = ForwardDiff.derivative(x -> f1(x, par), x1_)
	fp1 = ForwardDiff.gradient(p -> f1(x1_, p), par)
	x1_p = -fp1 / fx1
	x1_, x1_p
end

# ╔═╡ 8aca7493-e1d7-488f-a78b-e16517918fe4
x1_, x1_p = calc_partial_derivs_T1(transformSensPar(SensPar))

# ╔═╡ 687d9260-74dc-471a-a3bd-eb485db47d45
md"""
### Sensitivities
Formulation for sensitivity analysis: Calculate partial derivatives with ForwardDiff.jl.
```math
T_1 = T_1(T_2,G_0,\alpha_1,\epsilon_1,\alpha_2,\epsilon_2,h_{\text{conv}})
```
T1 = $(round(x1_,sigdigits=3)) °C
"""

# ╔═╡ 51d8997b-1169-467e-9da4-50cd4cf2ef8d
md"""
Sensitivities:
1.  $\partial T_1 / \partial T_2 =$  $(round(x1_p[1],sigdigits=2)) K/K
1.  $\partial T_1 / \partial G_0 =$  $(round(x1_p[2]*ufac"kW/m^2",sigdigits=2)) K/(kW/m^2)
1.  $\partial T_1 / \partial \alpha_1 =$  $(round(x1_p[3],sigdigits=2)) K
1.  $\partial T_1 / \partial \epsilon_1 =$  $(round(x1_p[4],sigdigits=2)) K
1.  $\partial T_1 / \partial \alpha_2 =$  $(round(x1_p[5],sigdigits=2)) K
1.  $\partial T_1 / \partial \epsilon_2 =$  $(round(x1_p[6],sigdigits=2)) K
1.  $\partial T_1 / \partial h_{\text{conv}} =$  $(round(x1_p[7],sigdigits=2)) K/(W/m^2/K)
"""

# ╔═╡ d19d7b31-d966-4a3b-a040-31667d8c4648
function calc_combined_uncertainty_T1(SensPar)
	x_, df_dx = calc_partial_derivs_T1(SensPar)
	uc = df_par_unc(SensPar).variance
	sqrt(sum(df_dx.^2 .* uc))
end

# ╔═╡ 379c89bb-07f6-4340-adc1-ca90e06eca90
calc_combined_uncertainty_T1(transformSensPar(SensPar))

# ╔═╡ 1a74b5da-438b-4e2d-9899-a1081d59d638
let
	p1 = Plots.plot(xlabel="Catalyst surface T (T2) / °C", ylabel="Thermocouple T (T1) / °C" )
	#T2_ = T2_enth(data)
	T2_ = T2
	T2r = floor(T2_/100)*100
	T2s = (T2r-100):20:(T2r+150)
	T1s = map((x) -> T1_conv(x,data), T2s)
	Plots.plot!(p1, T2s, T1s, label="T1(T2)")
	#Plots.plot!(p, T2s, T2s, label="T1=T2")
	
	T1_ = T1_conv(T2_,data)
	Plots.scatter!(p1, [T2_], [T1_], label="Parameter Point")
	Plots.annotate!(p1, T2_, T1_-5, text("T1= $(Integer(round(T1_))) °C\nT2= $(Integer(round(T2_))) °C",10))

	p2 = Plots.plot(xlabel="Catalyst surface T (T2) / °C", ylabel="ΔT (T1-T2) / °C", )
	Plots.plot!(p2, T2s, T1s.-T2s, label="ΔT")
	Plots.scatter!(p2, [T2_], [T1_-T2_], label="Parameter Point")
	Plots.annotate!(p2, T2_, T1_-T2_-10, text("ΔT= $(Integer(round(T1_-T2_))) °C",10))

	
	Plots.plot(p1,p2,size=(600,300))
end

# ╔═╡ 4791481a-777f-476b-b6e9-be747317a277
function runSensT2(SensPar;rel_deltas=rel_deltas)
	df = DataFrame(
		Parname=String[],
		Rel_delta=Float64[],
		G_lamp=Float64[],
		abs2=Float64[],
		eps2=Float64[],
		nflowin=Float64[],
		T_cat=Float64[],
		Delta_T_cat_rel=Float64[],
	)

	function calc_T2(SensPar)
		G_lamp, abs2, eps2, nflowin = SensPar	
		data = 
		newData(G_lamp=G_lamp*ufac"kW/m^2",abs2=abs2,eps2=eps2,nflowin=nflowin*ufac"mol/hr")
		T2_enth(data) 
	end

	T2mean = calc_T2(SensPar)
	for i in eachindex(SensPar)
		# mean
		
		push!(df, (names(df)[i+2], 1.0,SensPar..., T2mean, 0.0) )

		for rel_delta in rel_deltas
			# min
			SensPar[i] *= (1-rel_delta)
			push!(df, (names(df)[i+2],(1-rel_delta),SensPar..., calc_T2(SensPar), (calc_T2(SensPar)-T2mean)/T2mean*100) )
			
			# max
			SensPar[i] *= (1+rel_delta)/(1-rel_delta)
			push!(df, (names(df)[i+2],(1+rel_delta), SensPar..., calc_T2(SensPar),(calc_T2(SensPar)-T2mean)/T2mean*100) ) 
			
			SensPar[i] *= 1/(1+rel_delta)
		end
	end
	df
end

# ╔═╡ 61af4717-2a0e-4667-aec2-1e13d20296e2
sensT2 = runSensT2([G_lamp, abs2, eps2, nflowin])

# ╔═╡ 866e05f3-607d-46f5-8b4b-1424d66db3be
let	
	rel_deltas = filter(:Parname => ==("G_lamp"), sensT2).:Rel_delta
	parnames = unique(sensT2.:Parname)
	ctg = repeat(string.(rel_deltas), inner = length(parnames))
	nam = repeat(parnames, outer = length(rel_deltas))

	function get_records(col)
		arec = Float64[]
		for par in parnames
			rec = filter(:Parname => ==(par), sensT2)
			rec = rec[!,col]
			push!(arec, rec...)
		end
		arec = reshape(arec, length(rel_deltas), :)'
	end
	rec_abs = get_records(:T_cat)
	rec_rel = get_records(:Delta_T_cat_rel)

	T2_ = T2_enth(data)
	cpal = palette(:bam, length(rel_deltas))

	p1 = Plots.plot(ylabel="T Cat. Surface (T2) / °C", legend=:outertopright, ylim=(500.0,Inf), title = "T2 base = $(Integer(round(T2_))) °C")
	groupedbar!(p1, nam, rec_abs, group = ctg, bar_position = :dogde, bar_width=0.5, color_palette=cpal)
	Plots.plot!(p1, [0.0,4.0], [T2_,T2_] ,color=:black, label=:none)
	

	p2 = Plots.plot(ylabel="Rel. Δ / %", legend=:outertopright)
	groupedbar!(p2, nam, rec_rel, group = ctg, bar_position = :dogde, bar_width=0.5, color_palette=cpal)
	Plots.plot(p1,p2, layout=(2,1))
end

# ╔═╡ a8114301-8e45-4cd6-93e1-e7596532324f
function runSensT1(SensPar;rel_deltas=rel_deltas)
	df = DataFrame(
		Parname=String[],
		Rel_delta=Float64[],
		T_cat_surf=Float64[],
		G_lamp=Float64[],
		abs1=Float64[],
		eps1=Float64[],
		abs2=Float64[],
		eps2=Float64[],
		h_conv_D=Float64[],
		T_TC=Float64[],
		Delta_T_TC_rel=Float64[],
	)

	function calc_T1(SensPar)
		T2, G_lamp, abs1, eps1, abs2, eps2, h_conv_D = SensPar
		data = 
		newData(G_lamp=G_lamp*ufac"kW/m^2",abs1=abs1,eps1=eps1,abs2=abs2,eps2=eps2,h_conv_D=h_conv_D)
		T1_conv(T2,data)
	end

	T1mean = calc_T1(SensPar)
	for i in eachindex(SensPar)
		# mean
		push!(df, (names(df)[i+2], 1.0,SensPar..., T1mean, 0.0) )

		for rel_delta in rel_deltas
			# min
			SensPar[i] *= (1-rel_delta)
			push!(df, (names(df)[i+2],(1-rel_delta),SensPar..., calc_T1(SensPar), (calc_T1(SensPar)-T1mean)/T1mean*100) )
			
			# max
			SensPar[i] *= (1+rel_delta)/(1-rel_delta)
			push!(df, (names(df)[i+2],(1+rel_delta), SensPar..., calc_T1(SensPar), (calc_T1(SensPar)-T1mean)/T1mean*100 ) ) 
			
			SensPar[i] *= 1/(1+rel_delta)
		end
	end
	df
end

# ╔═╡ 9e282b41-7870-4416-9b5f-e873b25cf032
#sensT1 = runSensT1( [T2_enth(data), G_lamp, abs1, eps1, abs2, eps2, h_conv_D])
sensT1 = runSensT1( [T2, G_lamp, abs1, eps1, abs2, eps2, h_conv_D])

# ╔═╡ 964a2633-9c83-4fb7-8ac2-692018a536de
let	
	rel_deltas = filter(:Parname => ==("G_lamp"), sensT1).:Rel_delta
	parnames = unique(sensT1.:Parname)
	ctg = repeat(string.(rel_deltas), inner = length(parnames))
	nam = repeat(parnames, outer = length(rel_deltas))

	function get_records(col)
		adt = Float64[]
		for par in parnames
			dt = filter(:Parname => ==(par), sensT1)
			dt = dt[!,col]
			push!(adt, dt...)
		end
		adt = reshape(adt, length(rel_deltas), :)'
	end
	rec_abs = get_records(:T_TC)
	rec_rel = get_records(:Delta_T_TC_rel)

	#T1_ = T1_conv(T2_enth(data),data)
	T1_ = T1_conv(T2,data)

	p1 = Plots.plot(ylabel="Thermocouple (T1) / °C", legend=:outertopright, ylim=(0.85*T1_,Inf), title = "base T1= $(Integer(round(T1_))) °C")

	cpal = palette(:bam, length(rel_deltas))
	groupedbar!(p1, nam, rec_abs, group = ctg, bar_position = :dogde, bar_width=0.5, color_palette= cpal)
	Plots.plot!(p1, [0.0,7.0], [T1_,T1_] ,color=:black, label=:none)

	p2 = Plots.plot(ylabel="Rel. Δ / %", legend=:outertopright)
	groupedbar!(p2, nam, rec_rel, group = ctg, bar_position = :dogde, bar_width=0.5,color_palette= cpal)
	Plots.plot(p1,p2, layout=(2,1))
	

end

# ╔═╡ 2d2786e6-2851-4bc1-89b7-81911bfb6621
begin	
	Re_d = Re(D_TC,flow_velo(Across, data),data)
	Pr_ = Pr(data)
end;

# ╔═╡ 2faf0d97-0b79-41dc-86bc-6df1626661e7
Nu_D = Nu(Re_d,Pr_)

# ╔═╡ dc7d7f1d-e3a6-4f0e-bd10-8584854beeb1
h_conv(Nu_D,D_TC,data)

# ╔═╡ d1bd10aa-7593-4f21-805c-a7252dee5cde
md"""
-  $\bf{\text{Nu}}$ = $(round(Nu_D,sigdigits=2))
-  Convective heat transfer coeff. $\bf{h}$ = $(round(h_conv_D,sigdigits=3)) W/m^2/K
"""

# ╔═╡ 8f0c4275-8fe6-4e6f-a79e-e655d6dc104a
md"""
Obtain average flow velocity from mass feed flow and the cross sectional area that is permeable for gas flow:
$(
LocalResource(
	"./img/Topview_domain.png"
)
)
-  $A_{\text{cross}}$ = $(Across) m^2
-  $\dot m_{\text{in}}$ = $(round(data.mflowin, sigdigits=2)) kg/s
-  $\dot n_{\text{in}}$ = $(round(data.nflowin/ufac"mol/hr", sigdigits=2)) mol/h
-  $u_{\text{in}}$ = $(round(flow_velo(Across, data),sigdigits=2)) m/s
-  $T_{\text{gas,in}}$ = $(data.T_gas_in-273.15) °C
-  $\bf{\text{Re}_{\text d}}$ = $(round(Re_d,sigdigits=2))
"""

# ╔═╡ Cell order:
# ╠═bcf12f50-b3a6-11ee-0f2b-1b63d8c56346
# ╟─640a4e80-b06e-40b2-a659-13feed18f55a
# ╟─07cf6a69-c1eb-49a6-adbb-fee87efb9519
# ╠═9cc7316f-2654-425d-bc9d-1e26197f081c
# ╠═c3a8a9c9-113e-4bcc-9ad7-0ceb045a123c
# ╟─31d90ef7-e398-49ab-9338-45eb505b720e
# ╟─a56e991c-f6d7-4435-aaaa-6a1a80fe7fb6
# ╟─13fdff73-cb2b-4285-a43e-c81b3aca5c4a
# ╠═6355e676-013a-441b-959f-089ffc3010a4
# ╠═19a9466c-60a7-4b10-9ce7-f3f8f167712c
# ╠═e1ba6625-dda1-40cd-90f7-6d0604152edb
# ╟─ed68f853-f81a-4b24-a25b-4ecf41fd493a
# ╟─533362a7-7e32-4df7-98db-07030004fef2
# ╟─72ed6c5e-a7cc-4d58-aeff-4209621b6fb5
# ╠═1fcdb9b3-71c6-4b05-8460-d3c86504aad3
# ╠═9b3a6003-7263-40ed-a11b-3a55c3a7a510
# ╟─866e05f3-607d-46f5-8b4b-1424d66db3be
# ╠═61af4717-2a0e-4667-aec2-1e13d20296e2
# ╠═4f24e591-73d4-4b34-b812-24d80b531879
# ╠═854a6121-40a6-4d62-8846-02963ae3c136
# ╠═71e4bc55-058d-4c2d-9432-baae23a67ccb
# ╠═4791481a-777f-476b-b6e9-be747317a277
# ╟─bd4e7816-7790-4657-aa3c-78ce6bd578c7
# ╠═caae0418-636e-4a68-8cea-0c3a70704fee
# ╠═86b550c8-212b-4458-b701-2f150b87bcfd
# ╠═850a86cf-97db-40dd-994b-e9dd903a0522
# ╠═3a0c71bd-339f-437f-9705-ddcc543c3e8d
# ╠═8bd04d1e-08d6-495e-ad64-435a9b504c17
# ╠═1aefe230-88fc-4397-91d6-7899499cb3ea
# ╟─28249727-2511-4ec1-9067-e2b6fb0651f1
# ╠═687d9260-74dc-471a-a3bd-eb485db47d45
# ╟─51d8997b-1169-467e-9da4-50cd4cf2ef8d
# ╠═8aca7493-e1d7-488f-a78b-e16517918fe4
# ╟─b7db1ccc-7b9d-4332-bff7-04d9e97b7138
# ╠═964a2633-9c83-4fb7-8ac2-692018a536de
# ╠═9e282b41-7870-4416-9b5f-e873b25cf032
# ╠═1d6a9361-0473-488d-a56a-914e908b7a83
# ╠═06fb7b1b-3a36-4cd9-8485-8c2c55eab86c
# ╠═4027a2c5-9534-4b71-afe5-b83e0cc809d6
# ╟─382425a8-93c6-4fff-be01-3d4c196cdb7a
# ╠═379c89bb-07f6-4340-adc1-ca90e06eca90
# ╠═045f3f13-c606-4bf7-8416-2a95248ca35b
# ╠═d19d7b31-d966-4a3b-a040-31667d8c4648
# ╠═c6d4d255-29e1-4d68-9f7e-77afec7186d6
# ╠═82315b1c-bfe3-48f3-a537-d5ef0cc2312c
# ╠═a8114301-8e45-4cd6-93e1-e7596532324f
# ╟─c3c54b35-3121-4f33-9c09-8dd659bd4f84
# ╟─36a032e5-d0f4-49a2-98fb-e362189b0901
# ╟─8dc9a71b-759d-4440-a8c4-030f6a9e57e4
# ╠═c5c9bc3b-d460-4730-9637-64fd0443860e
# ╟─9073c7da-cc86-4191-b8a8-b49379803c89
# ╟─5229e919-38eb-41e2-acc4-0355bb933528
# ╠═1a74b5da-438b-4e2d-9899-a1081d59d638
# ╟─e7ecc92c-226f-4b00-9a29-72abda41225a
# ╟─f7a6f0af-b34a-4fa1-ad1b-579eca0862c7
# ╟─d7197aab-a552-47d0-8772-b647ff047d91
# ╠═8fc0ccbd-8108-4957-8ecd-4a0559fc6eed
# ╟─63ea26d4-8e31-4378-ba38-1d3f4132e4e1
# ╠═dc7d7f1d-e3a6-4f0e-bd10-8584854beeb1
# ╠═b9da0e38-6309-4469-9818-836f5bb77f52
# ╟─d1bd10aa-7593-4f21-805c-a7252dee5cde
# ╠═2faf0d97-0b79-41dc-86bc-6df1626661e7
# ╠═48d04afd-346f-423c-ad0c-94a0cabcdfc6
# ╟─7a968276-bd16-4b0f-86a8-edfb5c0bfb1a
# ╟─8f0c4275-8fe6-4e6f-a79e-e655d6dc104a
# ╠═e922cfd8-ea0f-416a-b5ef-3325a4905c3f
# ╠═2d2786e6-2851-4bc1-89b7-81911bfb6621
# ╠═2b98d5e6-1da8-46e1-8521-1ae8137e07a3
# ╠═91cb28c4-ed83-4427-9494-65ea9e436ca9
# ╠═654c6438-6e0e-4da6-81c8-d34acf96acb4
