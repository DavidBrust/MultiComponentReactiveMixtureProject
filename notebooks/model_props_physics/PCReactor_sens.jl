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

# ╔═╡ c21e1942-628c-11ee-2434-fd4adbdd2b93
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,"../.."))
	using Revise
	using VoronoiFVM
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve, LinearSolve,Pardiso, ILUZero, LinearAlgebra
	using StaticArrays

	using FiniteDifferences
	using LessUnitful
	using DataFrames, CSV
	using PlutoVista, Plots
	using PlutoUI, Colors
	using FixedBed
	using Printf
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═╡ skip_as_script = true
#=╠═╡
PlutoUI.TableOfContents(title="Photo Catalytic (PC) Reactor")
  ╠═╡ =#

# ╔═╡ b2791860-ae0f-412d-9082-bb2e27f990bc
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Introduction
Demonstration notebook for the photo thermal catalytic reactor (PCR) model. Solve energy equation alongside multicomponent species transport. Include reactive gas mixture (CO2,H2,CO,CH4,H2O,N2) with variable physical properties and a Ni based catalyst described with kinetics from published literature.

Select problem dimension: $(@bind dim Select([2, 3], default=2))

Check the box to __start the sensitivity calculation__: $(@bind RunSens PlutoUI.CheckBox(default=false))
"""
  ╠═╡ =#

# ╔═╡ 6e82c62c-4d3b-46ec-9d91-bdcbeb8693e2
@doc FixedBed.DMS_Info_isothermal

# ╔═╡ 670ea109-33be-42c0-a1fc-8959bbda5990
@doc FixedBed.DMS_Info_thermal()

# ╔═╡ 4e05ab31-7729-4a4b-9c14-145118477715
# ╠═╡ skip_as_script = true
#=╠═╡
if dim == 3
	@bind xcut Slider(linspace(0,16,17)*ufac"cm",show_value=true,default=8*ufac"cm")
end
  ╠═╡ =#

# ╔═╡ bc811695-8394-4c35-8ad6-25856fa29183
# function grid_boundaries_regions(dim;nref=0)
# 	Ω_catalyst = 2
# 	W=16
# 	H=0.5

# 	nz = 10*2^(nref)
# 	Z=(0:(H/nz):H)*ufac"cm"

# 	if dim == 2
# 		Γ_bottom = 1
# 		#Γ_bottom_insulating = 7
# 		Γ_side = 2
# 		Γ_sym = 4		
# 		Γ_top_permeable = 5
# 		Γ_top_irradiated = 6

# 		W=W/2 # axisymmetry, half domain is sufficient
# 		nr=W*2^(nref)
# 		R=(0:(W/nr):W)*ufac"cm"
		
# 		grid=simplexgrid(R,Z)
# 		circular_symmetric!(grid)
	
# 		cellmask!(grid,[0,9/10*H].*ufac"cm",[W,H].*ufac"cm",Ω_catalyst) # catalyst layer	
# 		bfacemask!(grid, [0,H].*ufac"cm",[W-1,H].*ufac"cm",Γ_top_permeable)
# 		bfacemask!(grid, [0,H].*ufac"cm",[W-2,0.5].*ufac"cm",Γ_top_irradiated) 
				
# 		inb = [Γ_top_permeable,Γ_top_irradiated]
# 		irrb = [Γ_top_irradiated]
# 		outb = [Γ_bottom]
# 		sb = [Γ_side]
# 	else
# 		Γ_side_1 = 1 
# 		Γ_side_2 = 2
# 		Γ_side_3 = 3
# 		Γ_side_4 = 4		
# 		Γ_bottom = 5
# 		Γ_top_permeable = 7
# 		Γ_top_irradiated = 8
		
# 		nxy=W*2^(nref)
# 		X=(0:(W/nxy):W)*ufac"cm"
# 		Y=X
# 		# X=(0:1:W)*ufac"cm"
# 		# Y=(0:1:W)*ufac"cm"
# 		# Z=(0:H/10:H)*ufac"cm"	
# 		grid=simplexgrid(X,Y,Z)
	
# 		# catalyst region
# 		cellmask!(grid,[0,0,9/10*H].*ufac"cm",[W,W,H].*ufac"cm",Ω_catalyst) # catalyst layer	
# 		bfacemask!(grid, [1,1,H].*ufac"cm",[W-1,W-1,H].*ufac"cm",Γ_top_permeable)
# 		bfacemask!(grid, [2,2,H].*ufac"cm",[W-2,W-2,H].*ufac"cm",Γ_top_irradiated)

# 		inb = [Γ_top_permeable,Γ_top_irradiated]
# 		irrb = [Γ_top_irradiated]
# 		outb = [Γ_bottom]
# 		sb = [Γ_side_1,Γ_side_2,Γ_side_3,Γ_side_4]
# 	end

# 	return grid, inb, irrb, outb, sb, [Ω_catalyst]
# end;

# ╔═╡ 289753d9-08a7-4447-94ac-efabdee99fea
md"""
# Sensitivities
Indices for surfaces and the involved properties:
- 1: Window surface (upper chamber, uc)
- 2: Catalyst surface (upper chamber)
- 3: Porous frit surface (lower chamber, lc)
- 4: Al plate surface (lower chamber)

Formulation for sensitivity analysis, $T_2$ corresponds to the temperature in the center of the catalyst sheet.

```math
\begin{align}
T_2 = T_2( &G_0, \dot n_{\text{feed,in}}, T_{\text{gas,in}}, X_0, \alpha^{\text{VIS}}_1,\alpha^{\text{IR}}_1, \alpha^{\text{VIS}}_2, \alpha^{\text{IR}}_2, \alpha^{\text{IR}}_3, \alpha^{\text{IR}}_4, \delta_{\text{uc}}, \delta_{\text{lc}}, \delta_{\text{gap,side}}, \\
& \text{Nu}, h_{\text{conv}}, \lambda_{\text{window}}, \lambda_{\text{Al,plate}})
\end{align}
```

"""

# ╔═╡ d522205e-d7fb-4261-a43c-b746339d4071
md"""
Mean values of main parameters. Calculate thermocouple hemispherical emissivity according to the proposed correlation for heavily oxidized Inconel 600 surfaces.

Operating Conditions:
-  $G_0$ = $(@bind G_lamp NumberField(40.0:1.0:100.0, default=70.0)) kW/m^2
-  $\dot n_{\text{in,total}}$ = $(@bind nflowin NumberField(1.0:0.2:20.0, default=7.4)) mol/hr
-  $T_{\text{gas,in}}$ = $(@bind T_gas_in NumberField(25.0:1.0:300.0, default=25)) °C
Window optical properties (upper chamber):
-  $\alpha_1^{\text{VIS}}$ = $(@bind abs1_vis NumberField(0.0:0.01:0.9, default=0.0))
-  $\alpha_1^{\text{IR}}$ = $(@bind abs1_IR NumberField(0.1:0.01:0.9, default=0.67))
Catalyst surface optical properties (upper chamber):
-  $\alpha_2^{\text{VIS}}$ = $(@bind abs2_vis NumberField(0.1:0.01:0.9, default=0.39))
-  $\alpha_2^{\text{IR}}$ = $(@bind abs2_IR NumberField(0.1:0.01:0.9, default=0.56))
Frit surface optical properties (lower chamber):
-  $\alpha_3^{\text{IR}}$ = $(@bind abs3_IR NumberField(0.1:0.01:0.9, default=0.55))
Al plate surface optical properties (lower chamber):
-  $\alpha_4^{\text{IR}}$ = $(@bind abs4_IR NumberField(0.1:0.01:0.9, default=0.8))
Catalyst loading and porous frit properties:
-  $m_{\text{cat}}$ = $(@bind mcat NumberField(100.0:1.0:750.0, default=500.0)) mg
-  $\text{porosity}$ = $(@bind poros NumberField(0.01:0.01:0.99, default=0.33))
-  $d_{\text p}$ (avg. pore diameter) = $(@bind dp NumberField(20.0:1.0:500.0, default=200.0)) μm
-  $\text{permeability}$ (for use in Darcy eq.) = $(@bind perm NumberField(5.0e-11:1.0e-11:6.0e-10, default=1.2e-10)) m²
-  $\lambda_{\text{frit,solid}}$ = $(@bind lambdas NumberField(0.6:0.01:3.0, default=1.13)) W/m/K
Geometric parameters (internal height in upper and lower chambers, gap betwwen frit and reactor wall):
-  $\delta_{\text{uc}}$ = $(@bind delta_uc NumberField(10.0:0.5:20.0, default=17.0)) mm
-  $\delta_{\text{lc}}$ = $(@bind delta_lc NumberField(10.0:0.5:20.0, default=18.0)) mm
-  $\delta_{\text{gap,side}}$ = $(@bind delta_gap NumberField(0.1:0.1:3.0, default=1.5)) mm
Parameters determining convective heat transport:
-  $\text{Nu}$ = $(@bind Nu NumberField(1.0:0.1:10.0, default=1.0)) 
-  $k_{\text{conv,outer}}$ = $(@bind k_conv NumberField(1.0:1.0:30.0, default=30.0)) W/m^2/K

"""

# ╔═╡ c7901b7f-27b9-469c-86f8-080df0cd3fa4
SensPar = [G_lamp, nflowin, T_gas_in, abs1_vis, abs1_IR, abs2_vis, abs2_IR, abs3_IR, abs4_IR, mcat, poros, dp, perm, lambdas, delta_uc, delta_lc, delta_gap, Nu, k_conv]

# ╔═╡ 2f03c691-ca0e-43d2-9577-c1eb245634bd
function transformSensPar(SensPar) 
	return [SensPar[1]*ufac"kW/m^2", SensPar[2]*ufac"mol/hr", SensPar[3]+273.15, SensPar[4:9]..., SensPar[10]*ufac"mg", SensPar[11], SensPar[12]*ufac"μm", SensPar[13], SensPar[14], SensPar[15:17].*ufac"mm"..., SensPar[18:end]... ]
end

# ╔═╡ 115aebe4-459b-4113-b689-38e381cdf64f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Combined Standard Uncertainty
Following calculation method in ISO/IEC Guide 98-3 (Guide to the Expression of Uncertainty in Measurement, GUM):
```math
u^2_{\text c} = \sum_{i=1}^N \left( \frac{\partial f}{\partial x_i} \right)^2 u^2(x_i)
```
Following procedure __Type B__ in determining the measurement uncertainty and assuming __rectangular distribution__ of estimate of measurand aroung the mean value, the variance $u^2(x_i)$ of an estimate $x_i$ is computed as $u^2(x_i)=a^2 /3$. Here $a$ is the symmetric interval of uncertainty of the measurement.
"""
  ╠═╡ =#

# ╔═╡ 9ed1cea4-e032-4f29-a746-e91519ff1d86
function par_unc(SensPar)
df = DataFrame(
		par_name=[:G0, :nflowin, :T_gas_in, :abs1_vis, :abs1_IR, :abs2_vis, :abs2_IR, :abs3_IR, :abs4_IR, :mcat, :poros, :dp, :perm, :lambdas, :delta_uc, :delta_lc, :delta_gap, :Nu, :k_conv],
		mean_value=SensPar,
		uncertainty_rel=[0.05, 0.02, 0.02, 0.05, 0.05, 0.20, 0.20, 0.20, 0.10, 0.05, 0.05, 0.05, 0.05, 0.05, 0.01, 0.01, 0.01, 1.0, 0.50]
)
	df[!, :uncertainty_interval] .= df.uncertainty_rel .* df.mean_value
	df[!, :variance] .= df.uncertainty_interval.^2 ./ 3 # rectangular distribution
	df
end

# ╔═╡ 6ed95794-b1f1-47ad-b655-5ef71e52776b
function par_unc(SensPar,dT2_dp)
	
	df = par_unc(SensPar)
	#df[!, :uncertainty_interval] .= df.uncertainty_rel .* SensPar
	#df[!, :uncertainty_interval] .= df.uncertainty_rel .* df.mean_value
	#df[!, :variance] .= df.uncertainty_interval.^2 ./ 3 # rectangular distribution
	df[!, :var_times_sens] .= df.variance .* dT2_dp.^2
	df
end

# ╔═╡ 8e710a52-937f-4a3a-9940-bcb5d1268281
function calc_combined_uncertainty_T2(SensPar, dT2_dp)	
	uc = par_unc(SensPar,dT2_dp).variance
	sqrt(sum(dT2_dp.^2 .* uc))
end

# ╔═╡ debdf189-7269-4d96-aa12-ff1c2c162f02
# ╠═╡ skip_as_script = true
#=╠═╡
par_unc(transformSensPar(SensPar))
  ╠═╡ =#

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	grid, inb,irrb,outb,sb,catr =  grid_boundaries_regions(dim)
	if dim == 2
		gridplot(grid, resolution=(660,300), aspect=4.0, zoom=2.8)
	else
		gridplot(grid,  xplane=xcut, resolution=(660,460), zoom=1.8, )
	end
end
  ╠═╡ =#

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
function PCR_base(dim, par; times=nothing)	
	times = isnothing(times) ? [0,20.0] : times
	
	G_lamp, nflowin, T_gas_in, abs1_vis, abs1_IR, abs2_vis, abs2_IR, abs3_IR, abs4_IR, mcat, poros, dp, perm, lambdas, delta_uc, delta_lc, delta_gap, Nu, k_conv= par

	uc_window = SurfaceOpticalProps( # upper chamber: window surface, 1
		alpha_IR=abs1_IR, 
		tau_IR=1-abs1_IR,
		alpha_vis=abs1_vis, 
		tau_vis=0.93 
	)
	
	uc_cat = SurfaceOpticalProps( # upper chamber: catalyst surface, 2
		alpha_IR=abs2_IR, 
		tau_IR=0.0,
		alpha_vis=abs2_vis, 
		tau_vis=0.0 
	)

	lc_frit = SurfaceOpticalProps( # lower chamber: porous frit surface, 3
		alpha_IR=abs3_IR, 
		tau_IR=0.0,
		alpha_vis=0.15, 
		tau_vis=0.0 
	)
	
	lc_plate = SurfaceOpticalProps( # lower chamber: Al plate surface, 4
		alpha_IR=abs4_IR, 
		tau_IR=0.0,
		alpha_vis=0.1, 
		tau_vis=0.0 
	)

	# grid, inb,irrb,outb,sb,catr =  grid_boundaries_regions(dim)
	grid, inb,irrb,outb,sb,catr =  FixedBed.grid_boundaries_regions(dim)
	
	data=ReactorData(
		dim=dim,
		nom_flux=G_lamp,
		nflowin=nflowin,
		T_gas_in=T_gas_in,
		uc_window=uc_window,
		uc_cat=uc_cat,
		lc_frit=lc_frit,
		lc_plate=lc_plate,
		uc_h=delta_uc,
		lc_h=delta_lc,
		Nu=Nu,
		k_nat_conv=k_conv,
		mcat=mcat,
		poros=poros,
		dp=dp,
		perm=perm,
		lambdas=lambdas,
		delta_gap=delta_gap,
		inlet_boundaries=inb,
		irradiated_boundaries=irrb,
		outlet_boundaries=outb,
		side_boundaries=sb,
		catalyst_regions=catr,
		rhos=5.0*ufac"kg/m^3" # set solid density to low value to reduce thermal inertia of system
		)

	# ##########################################################################
	# # BEGIN SYS INIT
	# (;p,ip,Tamb,iT,iTw,iTp,ibf,irradiated_boundaries,FluxIntp,ng,X0)=data
	# ng=ngas(data)

	# sys=VoronoiFVM.System( 	grid;
	# 						data=data,
	# 						flux=FixedBed.DMS_flux,
	# 						reaction=FixedBed.DMS_reaction,
	# 						storage=FixedBed.DMS_storage,
	# 						bcondition=FixedBed.PCR_bcond,
	# 						bflux=FixedBed.PCR_bflux,
	# 						bstorage=FixedBed.PCR_bstorage,
	# 						boutflow=FixedBed.DMS_boutflow,
	# 						outflowboundaries=outb,
	# 						assembly=:edgewise,
	# 						)

	# enable_species!(sys; species=collect(1:(ng+2))) # gas phase species xi, ptotal & T
	# enable_boundary_species!(sys, iTw, irrb) # window temperature as boundary species in upper chamber
	# enable_boundary_species!(sys, ibf, irrb) # boundary flux species, workaround to implement spatially varying irradiation
	# enable_boundary_species!(sys, iTp, outb) # plate temperature as boundary species in lower chamber

	# # END SYS INIT
	# ##########################################################################
	# # BEGIN INIVAL
	
	# inival=unknowns(sys)

	# inival[ip,:].=p
	# inival[[iT,iTw,iTp],:] .= Tamb

	# for i=1:ng
	# 	inival[i,:] .= X0[i]
	# end

	# # Only for variable irradiation flux density bc
	# function d3tod2(a,b)
	# 	a[1]=b[1]
	# 	a[2]=b[2]
	# end
	# inival[ibf,:] .= 0.0
	# sub=subgrid(grid,irradiated_boundaries,boundary=true, transform=d3tod2 )
		
	# for inode in sub[CellNodes]
	# 	c = sub[Coordinates][:,inode]
	# 	inodeip = sub[ExtendableGrids.NodeInParent][inode]
	# 	inival[ibf,inodeip] = FluxIntp(c[1]-0.02, c[2]-0.02)
	# end		

	# # END INIVAL
	# ##########################################################################
	# # BEGIN DATA INIT
	# catalyst_nodes = []
	# for reg in catr
	# 	catalyst_nodes = vcat(catalyst_nodes, unique(grid[CellNodes][:,grid[CellRegions] .== reg]) )
	# end
		
	# cat_vol = sum(nodevolumes(sys)[unique(catalyst_nodes)])

	# data.lcat = data.mcat/cat_vol
	# local Ain = 0.0
	# for boundary in inb
	# 	Ain += bareas(boundary,sys,grid)
	# end
	# data.mfluxin = data.mflowin / Ain
	# # END DATA INIT
	# ##########################################################################
	# BEGIN Solver control
	if dim == 2
		control = SolverControl(nothing, sys;)
	else
		control = SolverControl(;
        method_linear = KrylovJL_GMRES(
        ),
        precon_linear = VoronoiFVM.factorizationstrategy(
			MKLPardisoLU(), NoBlock(), sys),
   		)
	end
	control.handle_exceptions=true
	control.Δu_opt=100
	# END Solver control
	##########################################################################
		
	inival,sys = FixedBed.init_system!(dim, grid, data)

	solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="nae",log=true)

	return solt,grid,sys,data
end

# ╔═╡ 12cd7938-1f1d-417f-861b-b340edbd668d
function T2_part_deriv(dim,p,data)
	(;iT) = data

	function f(p) 
		solt,grid,sys,data=PCR_base(dim,p);
		sol = solt(solt.t[end])
		sol[iT,91] # corresponds to node at (r,z) = (0,5 mmm), center at cat surface
	end
	grad(forward_fdm(3, 1), f, p)[1]
end 

# ╔═╡ 8c8e91bd-7abb-453a-b10c-3b10203b44a0
# ╠═╡ skip_as_script = true
#=╠═╡
if RunSens
	dT2_dp = T2_part_deriv(2,transformSensPar(SensPar), ReactorData())
end
  ╠═╡ =#

# ╔═╡ 2801c5bd-991f-4f85-aa49-43927369605f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Sensitivities:
1.  $\partial T_2 / \partial G_0 =$  $(round(dT2_dp[1]*ufac"kW/m^2",sigdigits=2)) K/(kW/m^2)
1.  $\partial T_2 / \partial \dot n_{\text{in,total}} =$  $(round(dT2_dp[2]*ufac"mol/hr",sigdigits=2)) K/(mol/hr)
1.  $\partial T_2 / \partial T_{\text{gas,in}} =$  $(round(dT2_dp[3],sigdigits=2)) K/K
1.  $\partial T_2 / \partial \alpha_1^{\text{VIS}} =$  $(round(dT2_dp[4],sigdigits=2)) K
1.  $\partial T_2 / \partial \alpha_1^{\text{IR}} =$  $(round(dT2_dp[5],sigdigits=2)) K
1.  $\partial T_2 / \partial \alpha_2^{\text{VIS}} =$  $(round(dT2_dp[6],sigdigits=2)) K
1.  $\partial T_2 / \partial \alpha_2^{\text{IR}} =$  $(round(dT2_dp[7],sigdigits=2)) K
1.  $\partial T_2 / \partial \alpha_3^{\text{IR}} =$  $(round(dT2_dp[8],sigdigits=2)) K
1.  $\partial T_2 / \partial \alpha_4^{\text{IR}} =$  $(round(dT2_dp[9],sigdigits=2)) K
1.  $\partial T_2 / \partial m_{\text{cat}} =$  $(round(dT2_dp[10]*ufac"mg",sigdigits=2)) K/mg
1.  $\partial T_2 / \partial \text{porosity} =$  $(round(dT2_dp[11],sigdigits=2)) K
1.  $\partial T_2 / \partial d_{\text p} =$  $(round(dT2_dp[12]*ufac"μm",sigdigits=2)) K/μm
1.  $\partial T_2 / \partial \text{permeability} =$  $(round(dT2_dp[13],sigdigits=2)) K
1.  $\partial T_2 / \partial \lambda_{\text{frit,solid}} =$  $(round(dT2_dp[14],sigdigits=2)) K/(W/m/K)
1.  $\partial T_2 / \partial \delta_{\text{uc}} =$  $(round(dT2_dp[15]*ufac"mm",sigdigits=2)) K/mm
1.  $\partial T_2 / \partial \delta_{\text{lc}} =$  $(round(dT2_dp[16]*ufac"mm",sigdigits=2)) K/mm
1.  $\partial T_2 / \partial \delta_{\text{gap,side}} =$  $(round(dT2_dp[17]*ufac"mm",sigdigits=2)) K/mm
1.  $\partial T_2 / \partial \text{Nu}$ = $(round(dT2_dp[18],sigdigits=2)) K
1.  $\partial T_2 / \partial k_{\text{conv,outer}}$ =  $(round(dT2_dp[19],sigdigits=2)) K/(W/m^2/K)
"""
  ╠═╡ =#

# ╔═╡ e75c1f26-ba4b-417d-95da-e9f4203298e7
# ╠═╡ skip_as_script = true
#=╠═╡
calc_combined_uncertainty_T2(transformSensPar(SensPar),dT2_dp)
  ╠═╡ =#

# ╔═╡ 5afd82af-4911-4342-8091-e663ea0b4f16
# ╠═╡ skip_as_script = true
#=╠═╡
let
	df_par_unc = par_unc(transformSensPar(SensPar),dT2_dp)
	#CSV.write("./data/T2_parameter_uncertainty.csv", df_par_unc)
end
  ╠═╡ =#

# ╔═╡ fac7a69d-5d65-43ca-9bf3-7d9d0c9f2583
# ╠═╡ skip_as_script = true
#=╠═╡
solt,grid,sys,data=PCR_base(dim,transformSensPar(SensPar));
  ╠═╡ =#

# ╔═╡ 927dccb1-832b-4e83-a011-0efa1b3e9ffb
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Initialisation and Solve
The simulation is setup as a transient simulation. An initialisation strategy is employed where different physics are enabled step by step once a stationary state is established. Initially, no heat is transported and no chemical reactions take place. 

1. Velocity field (mass flow is ramped up from 1-100 % in T=$(data.dt_mf) s)
2. Temperature field through enthalpy flux (ramped up from 0-100 % in T=$(data.dt_hf_enth) s)
3. Temperature field through irradiation b.c. (ramped up from 0-100 % in T=$(data.dt_hf_irrad) s)

The mass flow boundary condition into the reactor domain is "ramped up" starting from a low value and linearly increasing until the final value is reached. A time delay is given to let the flow stabilize. Once the flow field is established, heat transport is ramped up until a stable temperature field is established. Finally, the reactivity of the catalyst is "ramped up" until its final reactivity value is reached.
"""
  ╠═╡ =#

# ╔═╡ 5d5ac33c-f738-4f9e-bcd2-efc43b638109
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;m,ip,iT,gn,gni,poros,mflowin,nflowin,W0,T_gas_in,Tamb,X0,outlet_boundaries,inlet_boundaries,dt_mf,dt_hf_enth)=data
	ng=ngas(data)
	vis=GridVisualizer(resolution=(600,300), xlabel="Time / s", ylabel="Molar flow / Total Moles")
	
	tfact=TestFunctionFactory(sys)	
	tf_out=testfunction(tfact,inlet_boundaries,outlet_boundaries)
	tf_in=testfunction(tfact,outlet_boundaries,inlet_boundaries)
		
	inflow_rate=Float64[]
	#inflow_rate_manual=Float64[]
	outflow_rate=Float64[]
	reaction_rate=Float64[]
	stored_amount=Float64[]

	#k=gni[:N2]
	k=iT
	for i=2:length(solt)
		m_ = 1
		W_ = 1
		#fac = k in 1:ng ? ufac"mol/hr" : ufac"kg/hr"
		if k in 1:ng
			m_ = m[k]
			W_ = W0[k]
			ifr=mflowin*W_*ramp(solt.t[i]; du=(0.1,1), dt=dt_mf)			
		elseif k == iT
			ifr=integrate(sys,tf_in,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])[iT]
			
			#ifr_manual=nflowin*(enthalpy_mix(data, T_gas_in, X0)-enthalpy_mix(data, Tamb, X0)) * ramp(solt.t[i]; du=(0.0,1), dt=dt_hf_enth)
			ifr_manual=nflowin*enthalpy_mix(data, T_gas_in, X0) * ramp(solt.t[i]; du=(0.0,1), dt=dt_hf_enth)
			#push!(inflow_rate_manual,ifr_manual)		
		end
		ofr=integrate(sys,tf_out,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])
		push!(inflow_rate,ifr/m_)
		
		push!(outflow_rate,ofr[k]/m_)		
		rr = integrate(sys,sys.physics.reaction,solt[i])[k,2]
		amount = sum(integrate(sys,sys.physics.storage,solt[i]), dims=2)[k]
		push!(reaction_rate, rr/m_)
		push!(stored_amount, amount/m_)
   	end

	
	# integrals
	I_in=0.0
	I_out=0.0
	I_reac=0.0
	for i=1:length(solt)-1
		I_in+=inflow_rate[i]*(solt.t[i+1]-solt.t[i])
		I_out+=outflow_rate[i]*(solt.t[i+1]-solt.t[i])
		I_reac+=reaction_rate[i]*(solt.t[i+1]-solt.t[i])
	end
	if k in 1:ng
		name = gn[k]
	elseif k == ip
		name = "Total Mass"
	elseif k == iT
		name = "Enthalpy Flow"
	end
	
	@printf "%s In: %2.2e \t Out: %2.2e \t React: %2.2e \nIn - Out: %2.4e \nStorage tEnd -t0: %2.4e" name I_in I_out I_reac I_in+I_out-I_reac stored_amount[end]-stored_amount[1]

	scalarplot!(vis, solt.t[2:end], inflow_rate, label="Inflow rate")
	#scalarplot!(vis, solt.t[2:end], inflow_rate_manual, label="Inflow MANUAL", color=:pink, clear=false)
	scalarplot!(vis, solt.t[2:end], outflow_rate, label="Outflow rate", color=:red, clear=false)	
	scalarplot!(vis, solt.t[2:end], -reaction_rate, label="Reaction rate",  color=:blue, clear=false)
	scalarplot!(vis, solt.t[2:end], stored_amount, label="Stored amount", color=:green, clear=false, )
	reveal(vis)

end
  ╠═╡ =#

# ╔═╡ 98468f9e-6dee-4b0b-8421-d77ac33012cc
md"""
### Temperature
1) Porous frit + catalyst layer domain
2) Window inner surface
3) Bottom plate
"""

# ╔═╡ f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╠═╡ skip_as_script = true
#=╠═╡
@bind t Slider(solt.t,show_value=true,default=solt.t[end])
  ╠═╡ =#

# ╔═╡ 5588790a-73d4-435d-950f-515ae2de923c
# ╠═╡ skip_as_script = true
#=╠═╡
sol = solt(t);
  ╠═╡ =#

# ╔═╡ 994d4a87-3f27-4a51-b061-6111c3346d60
# ╠═╡ skip_as_script = true
#=╠═╡
FixedBed.DMS_print_summary(sol,grid,sys,data)
  ╠═╡ =#

# ╔═╡ 3207839f-48a9-49b6-9861-e5e74bc593a4
# ╠═╡ skip_as_script = true
#=╠═╡
FixedBed.DMS_print_summary_ext(sol,sys,data)
  ╠═╡ =#

# ╔═╡ ec21bd68-27f5-4595-9f2c-ed99b06f503e
#=╠═╡
let
	(;iT) = data
	Tc = 0.0
	if dim == 2
		Tc = sol[iT,91]
	elseif dim == 3
		Tc = sol[iT,3035]
	end
	Tc -= 273.15
	@printf "Temperature on catalyst surface: %2.f °C" Tc 
end
  ╠═╡ =#

# ╔═╡ 99b59260-7651-45d0-b364-4f86db9927f8
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;iT,iTw,iTp,irradiated_boundaries,outlet_boundaries)=data
	#vis=GridVisualizer(layout=(3,1), resolution=(680,900))
	vis=GridVisualizer(layout=(1,1), resolution=(680,300))
	scalarplot!(vis[1,1],grid, sol[iT,:] .- 273.15, zoom = 2.8, aspect=4.0, levelalpha=0.0, show=true)
end
  ╠═╡ =#

# ╔═╡ 58c0b05d-bb0e-4a3f-af05-71782040c8b9
#=╠═╡
if dim == 2
md"""
- (1,1): T-profile at r=0
- (2,1): T-profile at z=0
- (1,2): Window T-profile
- (2,2): Bottom Plate T-profile
"""
end
  ╠═╡ =#

# ╔═╡ 8de4b22d-080c-486f-a6a9-41e8a5489966
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	if dim == 2
		(;iT,iTw,iTp,irradiated_boundaries,outlet_boundaries) = data
		vis=GridVisualizer(layout=(2,2), resolution=(680,600))
		function _2to1(a,b)
			a[1]=b[2]
		end
		_grid,_,_,_,_,_ = grid_boundaries_regions(dim)
		bfacemask!(_grid, [0.0,0.0].*ufac"cm",[0.0,0.5].*ufac"cm",8)
		grid1D = subgrid(_grid, [8]; boundary = true, transform = _2to1)
		sol1D=view(sol[iT, :], grid1D)
		scalarplot!(vis[1,1],grid1D, sol1D .-273.15, label="Temperature along Y-axis", clear=false)
	
		function __2to1(a,b)
			a[1]=b[1]
		end
		grid1D = subgrid(grid, outlet_boundaries; boundary = true, transform = __2to1)
		sol1D=view(sol[iT, :], grid1D)
		scalarplot!(vis[2,1],grid1D, sol1D .-273.15, label="Temperature along X-axis", clear=false)
		
	    # window
		bgridw = subgrid(grid, irradiated_boundaries; boundary = true, transform = __2to1)
		bsolw=view(sol[iTw, :], bgridw)
		scalarplot!(vis[1,2],bgridw, bsolw.-273.15,)
		# bottom plate
		bgridp = subgrid(grid, outlet_boundaries; boundary = true, transform = __2to1)
		bsolp=view(sol[iTp, :], bgridp)
		scalarplot!(vis[2,2],bgridp, bsolp.-273.15,show=true)
	end
end
  ╠═╡ =#

# ╔═╡ c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
md"""
### Molar fractions
"""

# ╔═╡ a4165336-17ae-42a7-823e-d75b58983a34
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;ip,p,gn,gni) = data
	ng=ngas(data)
	grid,_,_,_,_,_ = grid_boundaries_regions(dim)
	max_bfr = maximum(grid[BFaceRegions])
	bfacemask!(grid, [3.0,0.0].*ufac"cm",[3.0,0.5].*ufac"cm",max_bfr+1)
	function _2to1(a,b)
			a[1]=b[2]
	end
	grid1D = subgrid(grid, [max_bfr+1]; boundary = true, transform = _2to1)
	idc = length(grid1D[Coordinates])÷2
	
	xi_center_dry = Float64[]
	for i=1:ng
		sol1D=view(sol[i, :], grid1D)
		if i == gni[:H2O] || i == gni[:N2] 
			push!(xi_center_dry, 0.0)
		else
			push!(xi_center_dry, sol1D[idc])
		end
		
	end

	xi_center_dry /= sum(xi_center_dry)

	println("Dry Product Molar Fractions in center of reactor:")
	for i=1:ng
		#if i != gni[:H2O] && i != gni[:N2] 
		@printf "%3s: %2.1f%%\n" gn[i] xi_center_dry[i]*100
		#end
	end
end
  ╠═╡ =#

# ╔═╡ 2a8d25e4-e5e0-4eba-b25c-931a76e4b248
md"""
1) CO
1) CO2
1) CH4
"""

# ╔═╡ 111b1b1f-51a5-4069-a365-a713c92b79f4
# ╠═╡ show_logs = false
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;ip,p,gn,gni) = data
	ng=ngas(data)
	if dim == 2
		vis=GridVisualizer(layout=(4,1), resolution=(680,900))
		scalarplot!(vis[1,1], grid, sol[gni[:CO],:], aspect = 4.0,zoom = 2.8) # CO
		scalarplot!(vis[2,1], grid, sol[gni[:CO2],:], aspect = 4.0,zoom = 2.8) # CO2
		scalarplot!(vis[3,1], grid, sol[gni[:CH4],:], aspect = 4.0,zoom = 2.8) # CH4

		cols = distinguishable_colors(ng)
		# plot species molar fractions along frit thickness (along y axis)
		function _2to1(a,b)
			a[1]=b[2]
		end
		_grid,_,_,_,_,_ = grid_boundaries_regions(dim)
		max_bfr = maximum(grid[BFaceRegions])
		bfacemask!(_grid, [3.0,0.0].*ufac"cm",[3.0,0.5].*ufac"cm",max_bfr+1)
	    grid1D = subgrid(_grid, [max_bfr+1]; boundary = true, transform = _2to1)
		for i=1:ng
			sol1D=view(sol[i, :], grid1D)
			scalarplot!(vis[4,1],grid1D, sol1D, label=gn[i], color=cols[i],clear=false)
		end
		reveal(vis)
	
	else
		vis=GridVisualizer(layout=(3,1), resolution=(400,1200), outlinealpha=0.0)
		scalarplot!(vis[1,1], grid, sol[gni[:CO],:])
		scalarplot!(vis[2,1], grid, sol[gni[:CO2],:])
		scalarplot!(vis[3,1], grid, sol[gni[:CH4],:])
	end	
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─d3278ac7-db94-4119-8efd-4dd18107e248
# ╟─b2791860-ae0f-412d-9082-bb2e27f990bc
# ╠═6e82c62c-4d3b-46ec-9d91-bdcbeb8693e2
# ╠═670ea109-33be-42c0-a1fc-8959bbda5990
# ╠═4e05ab31-7729-4a4b-9c14-145118477715
# ╠═bc811695-8394-4c35-8ad6-25856fa29183
# ╟─289753d9-08a7-4447-94ac-efabdee99fea
# ╟─d522205e-d7fb-4261-a43c-b746339d4071
# ╠═12cd7938-1f1d-417f-861b-b340edbd668d
# ╟─2801c5bd-991f-4f85-aa49-43927369605f
# ╠═8c8e91bd-7abb-453a-b10c-3b10203b44a0
# ╠═c7901b7f-27b9-469c-86f8-080df0cd3fa4
# ╠═2f03c691-ca0e-43d2-9577-c1eb245634bd
# ╟─115aebe4-459b-4113-b689-38e381cdf64f
# ╠═e75c1f26-ba4b-417d-95da-e9f4203298e7
# ╠═8e710a52-937f-4a3a-9940-bcb5d1268281
# ╠═5afd82af-4911-4342-8091-e663ea0b4f16
# ╠═debdf189-7269-4d96-aa12-ff1c2c162f02
# ╠═9ed1cea4-e032-4f29-a746-e91519ff1d86
# ╠═6ed95794-b1f1-47ad-b655-5ef71e52776b
# ╟─927dccb1-832b-4e83-a011-0efa1b3e9ffb
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═fac7a69d-5d65-43ca-9bf3-7d9d0c9f2583
# ╠═5588790a-73d4-435d-950f-515ae2de923c
# ╠═994d4a87-3f27-4a51-b061-6111c3346d60
# ╟─3207839f-48a9-49b6-9861-e5e74bc593a4
# ╟─5d5ac33c-f738-4f9e-bcd2-efc43b638109
# ╟─98468f9e-6dee-4b0b-8421-d77ac33012cc
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╠═ec21bd68-27f5-4595-9f2c-ed99b06f503e
# ╠═99b59260-7651-45d0-b364-4f86db9927f8
# ╟─58c0b05d-bb0e-4a3f-af05-71782040c8b9
# ╟─8de4b22d-080c-486f-a6a9-41e8a5489966
# ╟─c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
# ╟─a4165336-17ae-42a7-823e-d75b58983a34
# ╟─2a8d25e4-e5e0-4eba-b25c-931a76e4b248
# ╠═111b1b1f-51a5-4069-a365-a713c92b79f4
