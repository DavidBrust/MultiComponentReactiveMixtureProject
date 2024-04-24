### A Pluto.jl notebook ###
# v0.19.40

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
	using Plots, Colors, PlutoUI, LaTeXStrings, Measures
	using Roots
	using ForwardDiff, DiffResults
	using MultiComponentReactiveMixtureProject

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

# ╔═╡ c7632b59-18c6-4992-9b81-1665e766a839
eps1_no(T) = 0.18 + 3.38*1.0e-4 *(T-25.0) # hemispherical emissivity, polished, pristine

# ╔═╡ 56e603cd-e6de-47bf-8488-b36f6e1a99df
eps1_no(650.0)

# ╔═╡ caae0418-636e-4a68-8cea-0c3a70704fee
eps1_mo(T) = 0.52 + 2.83*1.0e-4 *(T-25.0) # hemispherical emissivity, moderate oxidation

# ╔═╡ 5f4acf3d-55dd-4194-ad39-d455289658cb
eps1_mo(650.0)

# ╔═╡ 86b550c8-212b-4458-b701-2f150b87bcfd
eps1_ho(T) = 0.66 + 2.06*1.0e-4 *(T-25.0) # hemispherical emissivity, heavy oxidation

# ╔═╡ 983fdd0a-fe11-42d8-9f55-1952ac50f45a
eps1_ho(650.0)

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
function EB_TC_conv(T1, T2, nTC, data; phi10=phi10, phi12=phi12,)
	(;nom_flux, T_gas_in, k_nat_conv) = data
	FluxIntp = MultiComponentReactiveMixtureProject.flux_interpol(nom_flux)
	#G_lamp = nom_flux

   	d1 = 0.025*sqrt(2)/2
	d2 = 0.05*sqrt(2)/2
	#T_03 = [0.08,0.08] # y,x,z
	T_03 = [0.06,0.06] # y,x,z
	dTC_coord = Dict(
		3 => T_03,
		4 => T_03 .+ [d2, -d2],
		5 => T_03 .+ [-d1, d1],
		6 => T_03 .+ [d1, -d1],
		7 => T_03 .+ [-d2, -d2]
	)
	# uc
	
	
	G_lamp = FluxIntp(dTC_coord[nTC]...)
	G2_IR, G2_VIS = G2(T2,data.uc_cat.alpha_vis,data.uc_cat.eps,G_lamp)
	abs1_VIS = data.uc_mask.alpha_vis
	eps1 = data.uc_mask.eps
	abs1_IR = eps1
	Tgas_avg = 0.5*(T2+T_gas_in)
	
	phi12*(abs1_VIS*G2_VIS+abs1_IR*G2_IR) + abs1_VIS*phi10*G_lamp - eps1*ph"σ"*T1^4 - k_nat_conv*(T1-Tgas_avg)
end

# ╔═╡ 1aefe230-88fc-4397-91d6-7899499cb3ea
T1_conv(T2,data;nTC=3) = find_zero(x -> EB_TC_conv(x, T2+273.15, nTC, data), 800) -273.15

# ╔═╡ fd262db7-015a-4102-a909-7a392d75c3fe
md"""
### Process Data
Calculate Thermocouple temperatures for catalyst surface temperatures resulting from simulation at different locations. Data from:
- Simulated temperatures on catalysst surface: ```Tc_Uc_230804.csv```
- Time-averaged experimental TC temperatures: ```20230804_Chemical_Analysis.xlsx```
"""

# ╔═╡ 3086a90e-0a2e-423a-afa4-43ac1a94928c
md"""
```20230804_Chemical_Analysis.xlsx```
"""

# ╔═╡ c6740d22-2a2f-43de-8418-55e1030789c7
begin
	df_230804_Texp = DataFrame()
	df_230804_Texp[!, :exp_T_03_TC]  = [499.66, 562.08, 590.16, 615.60, 636.26, 658.71]
	df_230804_Texp[!, :exp_T_03_TC_Uc]  = [2.164, 2.434, 2.555, 2.666, 2.755, 2.852]

	df_230804_Texp[!, :exp_T_04_TC]  = [517.41, 576.10, 601.06, 623.86, 642.19, 661.62]	
	df_230804_Texp[!, :exp_T_04_TC_Uc]  = [2.240, 2.495, 2.603, 2.701, 2.781, 2.865]

	df_230804_Texp[!, :exp_T_05_TC]  = [515.59, 576.67, 604.24, 628.53, 648.27, 668.83]	
	df_230804_Texp[!, :exp_T_05_TC_Uc]  = [2.233, 2.497, 2.616, 2.722, 2.807, 2.896]
	
	df_230804_Texp[!, :exp_T_06_TC]  = [474.78, 537.28, 564.66, 589.61, 609.78, 631.80]	
	df_230804_Texp[!, :exp_T_06_TC_Uc]  = [2.056, 2.326, 2.445, 2.553, 2.640, 2.736]
	
	df_230804_Texp[!, :exp_T_07_TC]  = [444.53, 508.44, 538.06, 564.18, 584.92, 607.13]	
	df_230804_Texp[!, :exp_T_07_TC_Uc]  = [1.925, 2.202, 2.330, 2.443, 2.533, 2.629]
	
end;

# ╔═╡ 264c3b82-576a-430c-a61b-17091f9d1a99
df_230804_Texp

# ╔═╡ 3be5c8e0-d621-4670-aea3-3c928b38668a
md"""
```20230807_Chemical_Analysis.xlsx```
"""

# ╔═╡ 585e99c9-04ba-4e86-8e4a-434345f6c778
begin
	df_230807_Texp = DataFrame()
	df_230807_Texp[!, :exp_T_03_TC]  = [484.29, 594.97, 641.61, 667.63, 688.49]
	df_230807_Texp[!, :exp_T_03_TC_Uc]  = [2.097, 2.576, 2.7784, 2.897, 2.981]

	df_230807_Texp[!, :exp_T_04_TC]  = [503.09, 607.80, 649.97, 671.39, 689.58]
	df_230807_Texp[!, :exp_T_04_TC_Uc]  = [2.178, 2.632, 2.814, 2.907, 2.986]	

	df_230807_Texp[!, :exp_T_05_TC]  = [487.35, 605.35, 649.97, 674.35, 694.97]
	df_230807_Texp[!, :exp_T_05_TC_Uc]  = [2.110, 2.621, 2.814, 2.920, 3.009]	

	df_230807_Texp[!, :exp_T_06_TC]  = [445.33, 564.51, 608.05, 631.81, 651.56]
	df_230807_Texp[!, :exp_T_06_TC_Uc]  = [1.928, 2.444, 2.633, 2.736, 2.821]	

	df_230807_Texp[!, :exp_T_07_TC]  = [437.65, 552.90, 596.86, 621.45, 641.29]
	df_230807_Texp[!, :exp_T_07_TC_Uc]  = [1.895, 2.394, 2.584, 2.691, 2.777]	
end;

# ╔═╡ bbb6d613-1592-49ae-8929-12d10fa7c39c
function plot_Tcalc_exp(df, fflow)
	p=Plots.plot(xguide="Nominal flux / kW m⁻²", yguide="Temperature / °C", title="$(fflow) lₛ/min", legend=:bottomright)
	
	Plots.plot!(p,df.nom_flux, df.T_03, ribbon=df.T_03_Uc, fillalpha=.1, label="Cat. (calc.)", ls=:dash)
	Plots.plot!(p,df.nom_flux, df.calc_T_03_TC, ribbon=df.T_03_TC_Uc, fillalpha=.1, label="TC (calc.)", ls=:dashdot)
	Plots.plot!(p,df.nom_flux, df.exp_T_03_TC, marker=:utriangle, yerr=df.exp_T_03_TC_Uc, fillalpha=.2, label="TC (exp.)")
end

# ╔═╡ 8afe6a43-f0e0-4849-8405-9c4f72bafcd4
function parity_plot_Tcalc(df_calc,df_exp, fflow)
	df = hcat(df_calc,df_exp)
	p=Plots.plot(xguide="TC Temp. Exp. / °C", yguide="TC Temp. Calc. / °C", title="$(fflow) Lₛ/min", aspect_ratio=1,)
	x=400:1:700
	Plots.plot!(x,x,xlim=(400,700),ylim=(400,700), label=:none, color=:black)

	Plots.scatter!(p,df.exp_T_03_TC, df.calc_T_03_TC, label="T_03",c=1)
	Plots.scatter!(p,df.exp_T_04_TC, df.calc_T_04_TC, label="T_04",c=2)
	Plots.scatter!(p,df.exp_T_05_TC, df.calc_T_05_TC, label="T_05",c=3)
	Plots.scatter!(p,df.exp_T_06_TC, df.calc_T_06_TC, label="T_06",c=4)
	Plots.scatter!(p,df.exp_T_07_TC, df.calc_T_07_TC, label="T_07",c=5)
	
end

# ╔═╡ 14d41f40-059c-4826-92c1-d3a42260385f
md"""
### Comparison with chemical equilibrium
"""

# ╔═╡ 9b287a3a-bce9-483e-a788-ada05c5b46bd
md"""
Import experimental dry mole fractions for CO and CH₄ in the product stream from: ```20230804_Chemical_Analysis.xlsx```
"""

# ╔═╡ eafef0d2-e11d-49a8-9222-f889e6fc4beb
begin
	df_230804_xexp = DataFrame()
	# H2
	df_230804_xexp[!, :exp_xH2_dry]  = [0.453401886, 0.43099713, 0.423702601, 0.417574258, 0.413072117, 0.40916896]
	df_230804_xexp[!, :exp_xH2_dry_uc]  = [0.02044895, 0.014774924, 0.014522456, 0.014275272, 0.014133382, 0.01403622]
	# CO
	df_230804_xexp[!, :exp_xCO_dry] = [0.055340852, 0.098991119, 0.116675375, 0.131620628, 0.14273081, 0.15340121]
	df_230804_xexp[!, :exp_xCO_dry_uc]  = [0.008012769, 0.006369281, 0.006556156, 0.006724424, 0.006863756, 0.007013268]
	# CH4
	df_230804_xexp[!, :exp_xCH4_dry] = [0.006911911, 0.008137549, 0.007866462, 0.007600071, 0.007350697, 0.007018087]
	df_230804_xexp[!, :exp_xCH4_dry_uc] = [0.001243559, 0.00094606, 0.000948585, 0.000950745, 
	0.000952182, 0.000953957]
end;

# ╔═╡ 5ef6c167-2696-4c10-8c1f-8a8959db48a2
md"""
Import experimental dry mole fractions for CO and CH₄ in the product stream from: 
```20230807_Chemical_Analysis.xlsx```
"""

# ╔═╡ fe98be11-5ec2-420e-aa4a-46abcaac5496
begin
	df_230807_xexp = DataFrame()
	# CO
	df_230807_xexp[!, :exp_xCO_dry] = [0.025896312, 0.087946288, 0.111456893, 0.123443755, 0.132070297]	
	df_230807_xexp[!, :exp_xCO_dry_uc]  = [0.00717838, 0.00834811, 0.0077081, 0.007172525, 0.006305082]
	
	# CH4
	df_230807_xexp[!, :exp_xCH4_dry] = [0.005072607, 0.005928999, 0.005040328, 0.004338456, 0.003722562]	
	df_230807_xexp[!, :exp_xCH4_dry_uc] = [0.001037677, 0.001165188, 0.001045072, 0.000957133, 0.000832008]
end;

# ╔═╡ 9b3a6003-7263-40ed-a11b-3a55c3a7a510
rel_deltas=[0.1,0.2,0.3]

# ╔═╡ 06fb7b1b-3a36-4cd9-8485-8c2c55eab86c
function transformSensPar(SensPar)
	#return [SensPar[1] + 273.15, SensPar[2] * ufac"kW/m^2", SensPar[3:end]...]
	return [SensPar[1] + 273.15, SensPar[2] * ufac"kW/m^2", SensPar[3:end-1]..., SensPar[end]*ufac"mol/hr"]
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

# ╔═╡ ccc9f43e-1c1e-43c0-8fd5-6f19e79990be
md"""
## Experimental Thermocouple Uncertainty
A __rectangular distribution__ of the estimate of measurand aroung the mean value is assumed. The variance $u^2(x_i)$ of an estimate $x_i$ can thus be computed as $u^2(x_i)=a^2 /3$ where $a_+-a_-=2a$.
"""

# ╔═╡ 0f6953e7-3a37-4d41-bb15-7d2872770921
function TC_meas_uc(T_meas)
	# thernmocouple type N, class 2
	# in temperature range +333°C - +1200°C  measurement ucertainty of ±0.0075 * t
	a = 0.0075*T_meas
	a^2/3
end

# ╔═╡ 6ee8b930-cbc1-4cc2-8da8-86782c4addcf
md"""
## Experimental Flowmeter Uncertainty
"""

# ╔═╡ 3fc86d8d-9b26-464d-b526-c391183ac438
begin
	const ps = 1.01325*ufac"bar"
	const Ts = 293.15*ufac"K"
end

# ╔═╡ bc840826-828e-4441-98ed-659ddba29bb2
function u_FlowMeter(nflowin; FS=10_000*ufac"ml/minute", u_Rd=0.5/100.0, u_FS=0.1/100.0, u_repeat=0.3/100.0)

	ndot_FS = FS*ps/(ph"R"*Ts)
	a_stated = u_Rd*nflowin + u_FS*ndot_FS
	u_stated = (a_stated^2)/3

	a_repeat = u_repeat*nflowin
	u_repeat = (a_repeat^2)/3

	uc = sqrt(u_stated+u_repeat)
end

# ╔═╡ 4027a2c5-9534-4b71-afe5-b83e0cc809d6
function df_par_unc(SensPar,T2_uc)
	
	df = DataFrame(
		#par_name=[:T2, :G0, :abs1, :eps1, :abs2, :eps2, :h_conv],
		par_name=[:T2, :G0, :abs1, :eps1, :abs2, :eps2, :nflowin],
		mean_value=SensPar,
		#uncertainty_rel=[0.0, 0.05, 0.20, 0.20, 0.20, 0.20, 0.20],
		uncertainty_rel=[0.0, 0.05, 0.20, 0.20, 0.20, 0.20, 0.0],
	)
	
	df[!, :uncertainty_interval] .= df.uncertainty_rel .* df.mean_value
	@. df.uncertainty_interval = ifelse(df.par_name == :T2, T2_uc, df.uncertainty_interval)

	#replace!(df.uncertainty_interval, 0.0 => T2_uc)	
	
	df[!, :variance] .= df.uncertainty_interval.^2 ./ 3 # rectangular distribution
	@. df.variance = ifelse(df.par_name == :nflowin, u_FlowMeter(df.mean_value)^2, df.variance)	
	df[!, :standard_uncertainty] .= sqrt.(df.variance) # rectangular distribution
	df
end

# ╔═╡ 324da705-c27f-4b62-8859-df98146cd181
function u_FlowMeter_Vflow(nflowin; FS=10_000*ufac"ml/minute", u_Rd=0.5/100.0, u_FS=0.1/100.0, u_repeat=0.3/100.0)

	#ndot_FS = FS*ps/(ph"R"*Ts)
	Vflowin = nflowin*ph"R"*Ts/ps
	a_stated = u_Rd*Vflowin + u_FS*FS
	u_stated = (a_stated^2)/3

	a_repeat = u_repeat*Vflowin
	u_repeat = (a_repeat^2)/3

	uc = sqrt(u_stated+u_repeat)
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
-  $G_0$ = $(@bind G_lamp confirm(NumberField(40.0:1.0:100.0, default=70.0))) kW/m^2
-  $\dot n_{\text{in,total}}$ = $(@bind nflowin confirm(NumberField(1.0:0.2:20.0, default=7.4))) mol/hr
-  $h_{\text{conv}}$ = $(@bind h_conv_D confirm(NumberField(0.0:0.5:25.0, default=15))) W/m^2/K
-  $T_2$ = $(@bind T2 confirm(NumberField(400.0:1.0:800.0, default=590.0))) °C
-  $u_{\text c}(T_2)$ = $(@bind T2_uc confirm(NumberField(0.0:1.0:100.0, default=30.0))) °C
Thermocouple optical properties:
-  $\alpha_1$ (vis) = $(@bind abs1 confirm(NumberField(0.1:0.01:0.9, default=abs1_mean)))
-  $\epsilon_1$ (IR) = $(@bind eps1 confirm(NumberField(0.1:0.01:0.9, default=0.54)))
Catalyst surface optical properties:
-  $\alpha_2$ (vis) = $(@bind abs2 confirm(NumberField(0.1:0.01:0.9, default=0.39)))
-  $\epsilon_2$ (IR) = $(@bind eps2 confirm(NumberField(0.1:0.01:0.9, default=0.56)))
"""

# ╔═╡ c5c9bc3b-d460-4730-9637-64fd0443860e
#SensPar = [T2, G_lamp, abs1, eps1, abs2, eps2, h_conv_D]
SensPar = [T2, G_lamp, abs1, eps1, abs2, eps2, nflowin]

# ╔═╡ 045f3f13-c606-4bf7-8416-2a95248ca35b
df_par_unc(transformSensPar(SensPar),T2_uc)

# ╔═╡ 0c0d6d13-dd5e-4c9c-becc-3d41fbe627e1
let
	par_uc = df_par_unc(transformSensPar(SensPar),T2_uc)
	G0_uc = par_uc[par_uc.par_name.==:G0, :standard_uncertainty] / ufac"kW/m^2"
	nflowin_uc = par_uc[par_uc.par_name.==:nflowin, :standard_uncertainty] / ufac"mol/hr"
	G0_uc, nflowin_uc
end

# ╔═╡ 5229e919-38eb-41e2-acc4-0355bb933528
md"""
```math
\Delta T = T_1 - T_2
```
"""

# ╔═╡ f7a6f0af-b34a-4fa1-ad1b-579eca0862c7
mr=maximum(rel_deltas);

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
function h_conv(Nu_D,T2,data;D_TC=D_TC)
	(;X0, T_gas_in) = data
	T2 += 273.15
	T_film = 0.5*(T2 + T_gas_in)
	ηf, λf = dynvisc_thermcond_mix(data, T_film, X0)
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

# ╔═╡ 2b98d5e6-1da8-46e1-8521-1ae8137e07a3
function flow_velo(nflowin, Tcat, data; Across=Across, p=data.p, x=data.X0)
	(;T_gas_in) = data
	T_film = 0.5*(T_gas_in+Tcat+273.15)
	c = p / (ph"R"*T_film)
	Vdot = nflowin / c
	u = Vdot / Across
end

# ╔═╡ 91cb28c4-ed83-4427-9494-65ea9e436ca9
function Re(u,Tcat,data;d=D_TC, p=data.p, x=data.X0)
	(;T_gas_in) = data
	#T = 0.5*(T_gas_in+T2_enth(data)+273.15)
	T_film = 0.5*(T_gas_in+Tcat+273.15)
	ρf = density_idealgas(T_film, p, x, data)
	ηf, λf = dynvisc_thermcond_mix(data, T_film, x)
	
	Re = u*ρf*d/ηf # Reynolds number
end

# ╔═╡ 654c6438-6e0e-4da6-81c8-d34acf96acb4
function Pr(Tcat,data; p=data.p, x=data.X0)
	(;mflowin, T_gas_in, mmix0) = data
	T_film = 0.5*(T_gas_in+Tcat+273.15)
	ηf, λf = dynvisc_thermcond_mix(data, T_film, x)
	cf = heatcap_mix(data, T_film, x)
	
	Pr = cf/mmix0*ηf/λf # Prandtl number
end

# ╔═╡ c3a8a9c9-113e-4bcc-9ad7-0ceb045a123c
function newData(Tcat;abs1=abs1,eps1=eps1,G_lamp=G_lamp*ufac"kW/m^2",abs2=abs2,eps2=eps2,nflowin=nflowin*ufac"mol/hr",h_conv_D=nothing,D_TC=1.0*ufac"mm",Across=140^2*ufac"mm^2")

	data = ReactorData(
		#X0 = X0,
		nom_flux = G_lamp,
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
		Re_d = Re(flow_velo(nflowin,Tcat,data),Tcat,data)
		Nu_D = Nu(Re_d,Pr(Tcat,data))		
		data.k_nat_conv = h_conv(Nu_D,Tcat,data)
	else
		data.k_nat_conv = h_conv_D
	end
	data
end

# ╔═╡ 9cc7316f-2654-425d-bc9d-1e26197f081c
data = newData(T2)

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
-  $G_0$ = $(data.nom_flux/ufac"kW/m^2") kW/m^2
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
		nflowin=Float64[],
		T_TC=Float64[],
		Delta_T_TC_rel=Float64[],
	)

	function calc_T1(SensPar)
		#T2, G_lamp, abs1, eps1, abs2, eps2, h_conv_D = SensPar
		T2, G_lamp, abs1, eps1, abs2, eps2, nflowin = SensPar
		data = 
		newData(T2, G_lamp=G_lamp*ufac"kW/m^2",nflowin=nflowin*ufac"mol/hr",abs1=abs1,eps1=eps1,abs2=abs2,eps2=eps2)
		#newData(T2, G_lamp=G_lamp*ufac"kW/m^2",abs1=abs1,eps1=eps1,abs2=abs2,eps2=eps2,h_conv_D=h_conv_D)
		#newData(T2, G_lamp=G_lamp*ufac"kW/m^2",abs1=abs1,eps1=eps1,abs2=abs2,eps2=eps2)
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
#sensT1 = runSensT1( [T2, G_lamp, abs1, eps1, abs2, eps2, h_conv_D]);
sensT1 = runSensT1( [T2, G_lamp, abs1, eps1, abs2, eps2, nflowin])

# ╔═╡ 964a2633-9c83-4fb7-8ac2-692018a536de
let	
	rel_deltas = filter(:Parname => ==("G_lamp"), sensT1).Rel_delta
	parnames = unique(sensT1.Parname)
	ctg = repeat(string.(rel_deltas), inner = length(parnames))
	
	
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

	parnames = [L"T_2",L"G_0",L"\alpha_1^{\textrm{vis}}",L"\alpha_1^{\textrm{IR}}",L"\alpha_2^{\textrm{vis}}",L"\alpha_2^{\textrm{IR}}",L"\dot V_{\textrm{flow,in}}" ]
	nam = repeat(parnames, outer = length(rel_deltas))

	#T1_ = T1_conv(T2_enth(data),data)
	T1_ = T1_conv(T2,data)
	plot_font = "Computer Modern"
	default(fontfamily=plot_font)

	lbs=["70 %" "80 %" "90 %" "100 %" "110 %" "120 %" "130 %"]
	p1 = Plots.plot(ylabel=L"T_1\,\,/\,\, \textrm{°C}", legend=:outertopright, ylim=(0.85*T1_,Inf))

	cpal = palette(:bam, length(rel_deltas))
	groupedbar!(p1, nam, rec_abs, group = ctg, bar_position = :dogde, bar_width=0.5, color_palette= cpal, labels= lbs)
	#Plots.plot!(p1, [0.0,7.0], [T1_,T1_] ,color=:black, label=:none, xtickfontsize=12)

	p2 = Plots.plot(ylabel="Rel. change / %", legend=:outertopright)
	groupedbar!(p2, nam, rec_rel, group = ctg, bar_position = :dogde, bar_width=0.5,color_palette= cpal, labels=lbs)
	
	p = Plots.plot(p1,p2, layout=(2,1), xtickfontsize=12)
	Plots.reset_defaults()
	p
	#savefig(p, "../../data/out/2024-01-26/T_TC_sens.svg")	
end

# ╔═╡ 82315b1c-bfe3-48f3-a537-d5ef0cc2312c
function EB1_sens(T1, par, data;  phi10=phi10, phi12=phi12)
	#T2, G_lamp, abs1, eps1, abs2, eps2, h_conv_D = par
	T2, G_lamp, abs1, eps1, abs2, eps2, nflowin = par
	(;T_gas_in, X0) = data
	
	#(;FluxIntp) = ReactorData(nom_flux=G_lamp)
	#G_lamp = FluxIntp(0.08,0.08)
	
	G2_IR, G2_VIS = G2(T2,abs2,eps2,G_lamp)
	Tgas_avg = 0.5*(T2+T_gas_in)

	T2_ = T2 -273.15
	Re_d = Re(flow_velo(nflowin, T2_, data),T2_,data)
	Pr_ = Pr(T2_,data)
	Nu_D = Nu(Re_d,Pr_)
	h_conv_D=h_conv(Nu_D,T2_,data)
	
	#abs1*phi12*G2_ + abs1*phi10*G_lamp - eps1*ph"σ"*T1^4 - h_conv_D*(T1-Tgas_avg)
	phi12*(abs1*G2_VIS+eps1*G2_IR) + abs1*phi10*G_lamp - eps1*ph"σ"*T1^4 - h_conv_D*(T1-Tgas_avg) # abs1_IR = eps1
end

# ╔═╡ c6d4d255-29e1-4d68-9f7e-77afec7186d6
f1(x,p) = EB1_sens(x+273.15, p, data)

# ╔═╡ 1d6a9361-0473-488d-a56a-914e908b7a83
function calc_partial_derivs_T1(par)
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
T1 = $(x1_) °C
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
1.  $\partial T_1 / \partial \dot n_{\text{flow,in}} =$  $(round(x1_p[7]*ufac"mol/hr",sigdigits=2)) K/(mol/hr)
"""

# ╔═╡ d19d7b31-d966-4a3b-a040-31667d8c4648
function calc_combined_uncertainty_T1(SensPar; Tcat_uc=30.0)
	x_, df_dx = calc_partial_derivs_T1(SensPar)
	uc = df_par_unc(SensPar,Tcat_uc).variance
	sqrt(sum(df_dx.^2 .* uc))
end

# ╔═╡ 9c55c919-e373-40d3-aab5-8c6989b04e2f
function df_Tcalc(f)
	df = DataFrame(CSV.File(f))
	transform!(df, [:T_03, :nom_flux, :nflowin] => ByRow((T_03, nom_flux, nflowin)->T1_conv(T_03,newData(T_03,G_lamp=nom_flux*ufac"kW/m^2", nflowin=nflowin*ufac"mol/hr"),nTC=3)) => :calc_T_03_TC)

	transform!(df, [:T_04, :nom_flux, :nflowin] => ByRow((T_04, nom_flux, nflowin)->T1_conv(T_04,newData(T_04,G_lamp=nom_flux*ufac"kW/m^2", nflowin=nflowin*ufac"mol/hr"),nTC=4)) => :calc_T_04_TC)

	transform!(df, [:T_05, :nom_flux, :nflowin] => ByRow((T_05, nom_flux, nflowin)->T1_conv(T_05,newData(T_05,G_lamp=nom_flux*ufac"kW/m^2", nflowin=nflowin*ufac"mol/hr"),nTC=5)) => :calc_T_05_TC)

	transform!(df, [:T_06, :nom_flux, :nflowin] => ByRow((T_06, nom_flux, nflowin)->T1_conv(T_06,newData(T_06,G_lamp=nom_flux*ufac"kW/m^2", nflowin=nflowin*ufac"mol/hr"),nTC=6)) => :calc_T_06_TC)

	transform!(df, [:T_07, :nom_flux, :nflowin] => ByRow((T_07, nom_flux, nflowin)->T1_conv(T_07,newData(T_07,G_lamp=nom_flux*ufac"kW/m^2", nflowin=nflowin*ufac"mol/hr"),nTC=7)) => :calc_T_07_TC)
	
	# uncertainty of calculated TC temperature at position T_03
	transform!(df, [:T_03, :nom_flux, :T_03_Uc, :nflowin] => ByRow((T_03, nom_flux, T_03_Uc,nflowin)->calc_combined_uncertainty_T1(transformSensPar([T_03,nom_flux,SensPar[3:end-1]...,nflowin] ), Tcat_uc=T_03_Uc)) => :T_03_TC_Uc)

	transform!(df, [:T_04, :nom_flux, :T_03_Uc, :nflowin] => ByRow((T_04, nom_flux, T_03_Uc,nflowin)->calc_combined_uncertainty_T1(transformSensPar([T_04,nom_flux,SensPar[3:end-1]...,nflowin] ), Tcat_uc=T_03_Uc)) => :T_04_TC_Uc)	

	transform!(df, [:T_05, :nom_flux, :T_03_Uc, :nflowin] => ByRow((T_05, nom_flux, T_03_Uc,nflowin)->calc_combined_uncertainty_T1(transformSensPar([T_05,nom_flux,SensPar[3:end-1]...,nflowin] ), Tcat_uc=T_03_Uc)) => :T_05_TC_Uc)	

	transform!(df, [:T_06, :nom_flux, :T_03_Uc, :nflowin] => ByRow((T_06, nom_flux, T_03_Uc,nflowin)->calc_combined_uncertainty_T1(transformSensPar([T_06,nom_flux,SensPar[3:end-1]...,nflowin] ), Tcat_uc=T_03_Uc)) => :T_06_TC_Uc)	

	transform!(df, [:T_07, :nom_flux, :T_03_Uc, :nflowin] => ByRow((T_07, nom_flux, T_03_Uc,nflowin)->calc_combined_uncertainty_T1(transformSensPar([T_07,nom_flux,SensPar[3:end-1]...,nflowin] ), Tcat_uc=T_03_Uc)) => :T_07_TC_Uc)	
	return df
end

# ╔═╡ 673e65b6-63c3-462f-96c5-e658429ebf46
df_230804_Tcalc = df_Tcalc("../../../data/out/JECE_TemperaturesTC_Uncertainty/Tc_Uc_230804.csv" )

# ╔═╡ 8dbdd15a-7912-47c6-95ea-7b47fcb38c71
let
	path = "../../../data/out/2024-03-20/22_38_55/"
	df_230804_Tcalc_recalc = df_Tcalc(path * "Sim_T_probe_70.0suns_7.4.csv" )
	
	CSV.write(path * "Thermocouple_temps_calc.csv", df_230804_Tcalc_recalc)
	CSV.write(path * "Thermocouple_temps_exp.csv", df_230804_Texp)
end

# ╔═╡ b61ce7e8-8b66-4b42-9453-db6e37193423
df_230807_Tcalc = df_Tcalc("../../../data/out/JECE_TemperaturesTC_Uncertainty/Tc_Uc_230807.csv")

# ╔═╡ fc594076-cb09-4cc1-a447-ed7b24dc38d7
let
	p1 = parity_plot_Tcalc(df_230804_Tcalc,df_230804_Texp,3)
	p2 = parity_plot_Tcalc(df_230807_Tcalc,df_230807_Texp,5.6)
	p=Plots.plot(p1,p2,layout=(1,2), margin = 5mm, size=(550, 270))	
	#Plots.savefig(p, "../../../data/out/2024-04-24_JECE_revisions/ParityPlots.svg")
end

# ╔═╡ 4f67dc46-b5b6-4648-a526-77992c22020f
let
	# 04.08.2023, 40-70 suns, 3NL/min
	p1 = plot_Tcalc_exp(hcat(df_230804_Tcalc, df_230804_Texp), 3)
	# 07.08.2023, 40-80 suns, 5.6NL/min
	p2 = plot_Tcalc_exp(hcat(df_230807_Tcalc, df_230807_Texp), 5.6)
	p=Plots.plot(p1,p2,layout=(1,2), margin = 5mm, size=(550, 270))	
	#Plots.savefig(p, "../../../data/out/2024-04-24_JECE_revisions/Tcat_surf_TC_uc.svg")
end

# ╔═╡ 4389b52e-f286-41d1-93e6-e8090743e50a
begin
	df_TCE = CSV.read("../../data/RWGS_Methanation_equil/Thermo_Chem_equil_dry_composition_conv.csv", DataFrame)
	
	df_230804_T_x = hcat(df_230804_Tcalc, df_230804_xexp)
	df_230807_T_x = hcat(df_230807_Tcalc, df_230807_xexp)
end;

# ╔═╡ a745c8ad-ab53-4924-86e2-7bdeadee7293
function plot_mole_frac(df,fflow)
	p=Plots.plot(xguide="Temperature / °C", yguide="Dry molar fraction / %",
		legend=:left, xlim=(300,800), ylim=(-Inf, 20), title="$(fflow) lₛ/min")

	plot!(p, df_TCE.T,df_TCE.Mole_Frac_Dry_CO*100, label="CO TCE")
	plot!(p, df_TCE.T,df_TCE.Mole_Frac_Dry_CH4*100, label="CH₄ TCE",ls=:auto, c=2)

	# min T
	plot!(p, df.T_07,df.exp_xCO_dry*100, c=1, label=:none)
	plot!(p, df.T_03,df.exp_xCO_dry*100, c=1, label=:none)
	
	scatter!(p,df.T_03, df.exp_xCO_dry*100, xerr=df.T_03_Uc, yerr = df.exp_xCO_dry_uc*100, label="CO exp", c=1)
	scatter!(p, df.T_03, df.exp_xCH4_dry*100, xerr=df.T_03_Uc, yerr = df.exp_xCH4_dry_uc*100, label="CH₄ exp", c=2,)
	p
	
end

# ╔═╡ bf851b81-967f-4827-94dc-cf11d94b6bcb
let
	p1 = plot_mole_frac(df_230804_T_x,3)
	p2 = plot_mole_frac(df_230807_T_x,5.6)

	p=Plots.plot(p1,p2,layout=(1,2), size=(550, 280))	
	#Plots.savefig(p, "../../data/out/2024-01-26/Molefrac_TCE.svg")
end

# ╔═╡ 379c89bb-07f6-4340-adc1-ca90e06eca90
calc_combined_uncertainty_T1(transformSensPar(SensPar))

# ╔═╡ 2d2786e6-2851-4bc1-89b7-81911bfb6621
begin	
	Re_d = Re(flow_velo(nflowin* ufac"mol/hr", T2, data),T2,data)
	Pr_ = Pr(T2,data)
end;

# ╔═╡ 2faf0d97-0b79-41dc-86bc-6df1626661e7
Nu_D = Nu(Re_d,Pr_)

# ╔═╡ dc7d7f1d-e3a6-4f0e-bd10-8584854beeb1
h_conv(Nu_D,T2,data)

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
-  $\dot n_{\text{in}}$ = $(round(nflowin, sigdigits=2)) mol/h
-  $u_{\text{in}}$ = $(round(flow_velo(nflowin*ufac"mol/hr",T2,data),sigdigits=2)) m/s
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
# ╠═6355e676-013a-441b-959f-089ffc3010a4
# ╟─bd4e7816-7790-4657-aa3c-78ce6bd578c7
# ╠═c7632b59-18c6-4992-9b81-1665e766a839
# ╠═56e603cd-e6de-47bf-8488-b36f6e1a99df
# ╠═caae0418-636e-4a68-8cea-0c3a70704fee
# ╠═5f4acf3d-55dd-4194-ad39-d455289658cb
# ╠═86b550c8-212b-4458-b701-2f150b87bcfd
# ╠═983fdd0a-fe11-42d8-9f55-1952ac50f45a
# ╠═850a86cf-97db-40dd-994b-e9dd903a0522
# ╠═3a0c71bd-339f-437f-9705-ddcc543c3e8d
# ╠═8bd04d1e-08d6-495e-ad64-435a9b504c17
# ╠═1aefe230-88fc-4397-91d6-7899499cb3ea
# ╟─28249727-2511-4ec1-9067-e2b6fb0651f1
# ╟─fd262db7-015a-4102-a909-7a392d75c3fe
# ╠═9c55c919-e373-40d3-aab5-8c6989b04e2f
# ╟─3086a90e-0a2e-423a-afa4-43ac1a94928c
# ╠═c6740d22-2a2f-43de-8418-55e1030789c7
# ╠═264c3b82-576a-430c-a61b-17091f9d1a99
# ╟─3be5c8e0-d621-4670-aea3-3c928b38668a
# ╠═585e99c9-04ba-4e86-8e4a-434345f6c778
# ╠═673e65b6-63c3-462f-96c5-e658429ebf46
# ╠═8dbdd15a-7912-47c6-95ea-7b47fcb38c71
# ╠═b61ce7e8-8b66-4b42-9453-db6e37193423
# ╠═bbb6d613-1592-49ae-8929-12d10fa7c39c
# ╠═8afe6a43-f0e0-4849-8405-9c4f72bafcd4
# ╠═fc594076-cb09-4cc1-a447-ed7b24dc38d7
# ╠═4f67dc46-b5b6-4648-a526-77992c22020f
# ╟─14d41f40-059c-4826-92c1-d3a42260385f
# ╟─9b287a3a-bce9-483e-a788-ada05c5b46bd
# ╠═eafef0d2-e11d-49a8-9222-f889e6fc4beb
# ╟─5ef6c167-2696-4c10-8c1f-8a8959db48a2
# ╠═fe98be11-5ec2-420e-aa4a-46abcaac5496
# ╠═bf851b81-967f-4827-94dc-cf11d94b6bcb
# ╠═a745c8ad-ab53-4924-86e2-7bdeadee7293
# ╠═4389b52e-f286-41d1-93e6-e8090743e50a
# ╟─687d9260-74dc-471a-a3bd-eb485db47d45
# ╟─51d8997b-1169-467e-9da4-50cd4cf2ef8d
# ╠═8aca7493-e1d7-488f-a78b-e16517918fe4
# ╟─b7db1ccc-7b9d-4332-bff7-04d9e97b7138
# ╠═9b3a6003-7263-40ed-a11b-3a55c3a7a510
# ╠═964a2633-9c83-4fb7-8ac2-692018a536de
# ╠═9e282b41-7870-4416-9b5f-e873b25cf032
# ╠═1d6a9361-0473-488d-a56a-914e908b7a83
# ╠═06fb7b1b-3a36-4cd9-8485-8c2c55eab86c
# ╠═4027a2c5-9534-4b71-afe5-b83e0cc809d6
# ╟─382425a8-93c6-4fff-be01-3d4c196cdb7a
# ╠═379c89bb-07f6-4340-adc1-ca90e06eca90
# ╠═045f3f13-c606-4bf7-8416-2a95248ca35b
# ╠═0c0d6d13-dd5e-4c9c-becc-3d41fbe627e1
# ╠═d19d7b31-d966-4a3b-a040-31667d8c4648
# ╠═82315b1c-bfe3-48f3-a537-d5ef0cc2312c
# ╠═c6d4d255-29e1-4d68-9f7e-77afec7186d6
# ╠═a8114301-8e45-4cd6-93e1-e7596532324f
# ╟─ccc9f43e-1c1e-43c0-8fd5-6f19e79990be
# ╠═0f6953e7-3a37-4d41-bb15-7d2872770921
# ╟─6ee8b930-cbc1-4cc2-8da8-86782c4addcf
# ╠═3fc86d8d-9b26-464d-b526-c391183ac438
# ╠═bc840826-828e-4441-98ed-659ddba29bb2
# ╠═324da705-c27f-4b62-8859-df98146cd181
# ╟─c3c54b35-3121-4f33-9c09-8dd659bd4f84
# ╟─36a032e5-d0f4-49a2-98fb-e362189b0901
# ╠═8dc9a71b-759d-4440-a8c4-030f6a9e57e4
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
