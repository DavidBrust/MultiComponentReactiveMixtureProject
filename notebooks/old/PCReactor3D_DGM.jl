### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ 11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using VoronoiFVM, VoronoiFVM.SolverStrategies
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve
	using LinearAlgebra
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots
	using PlutoUI, HypertextLiteral
	using CSV,DataFrames
	using Interpolations

	using FixedBed
	
	#GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="Photo-Catalytic Reactor",depth=6)

# ╔═╡ 94ed0f8b-13c1-4460-a391-5057cff401e0
md"""
Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=false))
"""

# ╔═╡ 2ed3223e-a604-410e-93d4-016580f49093
md"""
# Domain / Grid
"""

# ╔═╡ 390c7839-618d-4ade-b9be-ee9ed09a77aa
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The modelling domain covers a prism of square shape, which represents the porous foam that is supporting the catalyst layer in the photo-catalytic reactor.
$(LocalResource("../img/Domain.png")) 

"""
  ╠═╡ =#

# ╔═╡ 0a911687-aff4-4c77-8def-084293329f35
begin
	const Γ_side_front = 1 # symmetry bc
	const Γ_side_right = 2 # wall bc
	const Γ_side_back = 3 # wall bc
	const Γ_side_left = 4 # symmetry bc
	const Γ_bottom_rim = 5 # outer rim, contacting reactor housing
	#const Γ_top_frit = 6 # inflow bc, uncoated porous frit 
	#const Γ_top_cat = 7 # inflow bc, catalyst coated porous frit
	const Γ_top_cat = 6 # inflow bc, catalyst coated porous frit 
	#const Γ_bottom = 8 # outflow bc
	const Γ_bottom = 7 # outflow bc
end;

# ╔═╡ 37adb8da-3ad5-4b41-8f08-85da19e15a53
function prism_sq_full(data; nref=0, w=data.wi, h=data.h, cath=data.cath, catwi=data.catwi)
	
	hw=(w/2.0)/8.0*2.0^(-nref) # width of ~16 cm, divide into 8 segments
	W=collect(-(w/2.0):hw:(w/2.0))
	
	hhfrit=h/5.0
	#Hfrit=collect(0:hhfrit:h)
	#Hfrit=collect(0:hhfrit:(h-cath))
	Hfrit=linspace(0,h-cath,6)
	
	hhCL=cath/5.0
	#HCL=collect(h:hhCL:h+cath)
	#HCL=collect((h-cath):hhCL:h)
	HCL=linspace((h-cath),h,6)
	H=glue(Hfrit,HCL)
	
	grid=simplexgrid(W,W,H)
	
	# catalyst layer region
	#cellmask!(grid,[-catwi/2,-catwi/2,h],[catwi/2,catwi/2,h+cath],2)
	cellmask!(grid,[-w/2,-w/2,h-cath],[w/2,w/2,h],2)
	# catalyst layer boundary
	bfacemask!(grid,[-w/2,-w/2,h],[w/2,w/2,h],Γ_top_cat)
	# outer rim boundary
	bfacemask!(grid,[-(w/2-1.0*ufac"cm"),-(w/2-1.0*ufac"cm"),0],[w/2-1.0*ufac"cm",w/2-1.0*ufac"cm",0],Γ_bottom)
end

# ╔═╡ e21b9d37-941c-4f2c-9bdf-956964428f90
const grid_fun = prism_sq_full

# ╔═╡ 2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
md"""
# Governing system
"""

# ╔═╡ f4dcde90-6d8f-4b17-b4ec-367d2372637f
md"""
## Species Mass balances
There are $\nu$ gas phase species in the system, gas phase composition is expressed in partial pressures $p_i, i = 1 ... \nu$.
"""

# ╔═╡ 3703afb0-93c4-4664-affe-b723758fb56b
md"""
```math
\begin{align}
- \nabla \cdot N_i + R_i &= 0
~,
i = 1 ... \nu
\end{align}

```
where $N_i$ is the molar flux ($\frac{\text{mol}}{\text{m}^2 \text{s}}$) and $R_i$ is the molar volumetric source/sink ($\frac{\text{mol}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
"""

# ╔═╡ 21d0195b-b170-460d-989e-f9d00b511237
md"""
### Isothermal Dusty Gas Model
"""

# ╔═╡ 8f4843c6-8d2b-4e24-b6f8-4eaf3dfc9bf0
md"""
In the framework of the DGM a system with ``\nu`` gas phase species has ``\nu+1`` unknowns: the species composition in terms of partial pressures ``p_i`` and the total pressure ``p``. The ``\nu`` independent equations expressing the flux-driving force relationships for each of the gas phase species is supplemented by the constraint that the partial pressures sum to the total pressure.
"""

# ╔═╡ 66b55f6b-1af5-438d-aaa8-fe4745e85426
md"""
```math
\begin{align}
	\sum_{j=1 \atop j \neq i}^{\nu} \frac{x_j N_i-x_i N_j}{D_{ij}^{\text{eff}}} + \frac{N_i}{D_{i,\text K}^{\text{eff}}} &= -\frac{\nabla p_i}{RT} - \frac{1}{D_{i,K}^{\text{eff}}} \frac{p_i}{RT} \left(\frac{\kappa}{\mu} \nabla p \right) \\
	\sum_{i=1}^\nu p_i &= p
\end{align}
```
"""

# ╔═╡ 8528e15f-cce7-44d7-ac17-432f92cc5f53
md"""
The above form of the DGM neglects thermal diffusion and is formulated for a mixture of ideal gases. The driving force for transport of species $i$ relative to the other species (diffusive transport) then reduces to $\nabla p_i$. Transport of species due to the viscous flow through the porous material in the direction of a negative pressure gradient is captured by D'arcy law (convective transport).
"""

# ╔═╡ a6afe118-dcbd-4126-8646-c7268acfacf3
md"""
The numerical fluxes $\textbf{J}$ could be computed using the backslash operator (Julia solver for linear systems) 
``\textbf{J} 
= 
M
\ \backslash\ 
\textbf{F}``.
However, this allocates a return vector (and internally a pivoting vector).

For a non-allocating variant, VoronoiFVM allows to perform `J.=F ; inplace_linsolve!(M,J)`

In order to be stack allocated, temporary arrays as M, F, J shall be created as `MArray` (available from StaticArrays.jl) with size information known at compile time.
I addition, the call to `inplace_linsolve` needs to be inlined. Conveniently, with Julia 1.8, "call-site  inlining" is available.

An alternative is to use `StrideArray` and `@gc_preserve`, see 
https://discourse.julialang.org/t/what-is-stridearrays-jl/97146/23

Compared to StrideArrays.jl, StaticArrays.jl is the more mature package. Moreover, `@gc_preserve` dose not work for calls with return values. So it seems to be reasonable to stick to `MArray` and inlining, the more with Julia 1.8, callsite inline is aviablable. 


"""

# ╔═╡ a60ce05e-8d92-4172-b4c1-ac3221c54fe5
md"""
### Knudsen effective Diffusivity
"""

# ╔═╡ 24374b7a-ce77-45f0-a7a0-c47a224a0b06
md"""
In the DGM the solid porous matrix is considered as another species in the ideal gas mix. The Knudsen diffusion coefficients describe the interactions between the gas molecules and the solid porous matrix. Analogous to the gas-gas diffusion coefficients, the Knudsen diffusion coefficients quantify the resulting friction force from momentum transfer upon collision with the walls acting on the gas molecules.
"""

# ╔═╡ 4865804f-d385-4a1a-9953-5ac66ea50057
md"""
Calculation of Knudsen diffusion coefficients according to __Wesselingh, J. A., & Krishna, R. (2006).__ Mass Transfer in Multicomponent Mixtures (ch. 21, Fig. 21.8)
"""

# ╔═╡ 3bf71cea-4f73-47da-b5ed-2cae3ec3d18b
md"""
## Thermal Energy Transport
"""

# ╔═╡ 7f94d703-2759-4fe1-a8c8-ddf26732a6ca
md"""
```math
\begin{align}
- \nabla \cdot \left(\lambda_{\text{eff}} \nabla T - \sum_i^{\nu} \vec N_i h_i \right)  &= 0
\end{align}
```

where $\lambda_{\text{eff}}$ is the effective thermal conductivity. The heat  release from chemical reactions is already captured by the species enthalpies. The convective heat transport within the porous medium is expressed via the sum over the product of molar species fluxes with species molar enthalpies $\vec H_i= \vec N_i h_i$.
"""

# ╔═╡ 906ad096-4f0c-4640-ad3e-9632261902e3
md"""
# Boundary Conditions
"""

# ╔═╡ 8139166e-42f9-41c3-a360-50d3d4e5ee86
md"""
## Gas species
"""

# ╔═╡ 44d91c2e-8082-4a90-89cc-81aba783d5ac
md"""
Boundary conditions for the transport of gas phase species cover in and outflow boundary conditions at the bottom and top surfaces of the modelling domain with no-flux conditins applied elsewhere. In the catalyst layer, volumetric catalytic reactions take place.
"""

# ╔═╡ 39e74955-aab6-4bba-a1b8-b2307b45e673
md"""
## Thermal
"""

# ╔═╡ 6798d5e2-b8c7-4f54-aa71-6ea1ccab78fb
md"""
The thermal boundary conditions consist of heat transport processes over the modelling domain boundaries. The processes consist of radiative losses on the top (reflection of incoming light, transmission through window, thermal emission) and bottom (thermal emission), convective heat losses on all outer surfaces, as well as enthalpy in and outflows through the top and bottom surfaces of the modeling domain. Heat is transported through the gas chambers of the reactor via convective and radiative heat transfer.
"""

# ╔═╡ ed3609cb-8483-4184-a385-dca307d13f17
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
$(LocalResource("../img/ThermalBC_new.png")) 
"""
  ╠═╡ =#

# ╔═╡ 58d0610b-1739-4260-8d16-5a31ba362d69
md"""
## Irradiation
Irradiation flux coming from the solar simulator enters the aperture of the reactor with the specified flux profile as determined via photo-metric measurement. The measured profile is imported from a csv-datafile, handled by DataFrames.jl and interpolated by methodes provided by Interpolations.jl to be used as boundary condition in the simulation.
"""

# ╔═╡ d0993435-ed0b-4a43-b848-26f5266017a1
begin
	
	FluxMap = let
		#ff = "FlowPhotoChem_Messdaten_20230523_110038.csv" # 40 suns
		ff = "FlowPhotoChem_Messdaten_20230523_105025.csv" # 80 suns
		#ff = "FlowPhotoChem_Messdaten_20230523_104820.csv" # 100 suns
			if splitdir(pwd())[2]=="FixedBed"
				filepath="data/IrradiationFluxProfiles/"*ff
			else
				filepath="../data/IrradiationFluxProfiles/"*ff
			end
		CSV.read(filepath, DataFrame, delim=";",header=false)
	end
end;

# ╔═╡ 634d1042-b110-45ef-bfbe-51b827fc922f
md"""
In the following, only the section of the total flux profile entering the __12 cm x 12 cm aperture__ of the reactor is needed for the simulation.
"""

# ╔═╡ 29f34e55-91ab-4b6d-adb2-a58412af95f6
# domain width (porous frit = 15.7 cm, ~ 16.0)
function sel12by12(M;wi=16.0*ufac"cm",wi_apt=12.0*ufac"cm")
	
	
	
	# starting coordinates for optimum 10 cm x 10 cm selection
	sr=33
	sc=33

	# (inverse) resolution of flux measurements, distance between data points
	Dx = 0.031497*ufac"cm"
	Dy = 0.031497*ufac"cm"
	
	
	# calculate coordinate offsets, when deviating (expanding/shrinking from the center) from 10cm x 10cm selection 	
	Dsr=Integer(round((wi_apt-10.0*ufac"cm")/(2*Dy)))
	Dsc=Integer(round((wi_apt-10.0*ufac"cm")/(2*Dx)))

	sr-=Dsr
	sc-=Dsc

	nx=Integer(round(wi_apt/Dx))
	ny=Integer(round(wi_apt/Dy))
	
	M=Matrix(FluxMap)
	# coordinate system of measurement data matrix has its origin at bottom left corner, the excel data matrix (and its coordiinates) start in top left corner
	reverse!(M; dims=1) 
	@views M_ = M[sr:(sr+ny-1),sc:(sc+nx-1)]*ufac"kW/m^2"

	# pad Flux Map with zeros outside of aperture area
	nx_dom = Integer(round(wi/Dx))
	ny_dom = Integer(round(wi/Dy))
	M__ = zeros(nx_dom,ny_dom)
		
	Dsr_dom_apt=Integer(round((wi - wi_apt)/(2*Dy)))
	Dsc_dom_apt=Integer(round((wi - wi_apt)/(2*Dx)))

	M__[(Dsr_dom_apt+1):(Dsr_dom_apt+ny), (Dsc_dom_apt+1):(Dsc_dom_apt+nx)] = M[sr:(sr+ny-1),sc:(sc+nx-1)]*ufac"kW/m^2"
	#Dsc=Integer(round((wi_apt-10.0*ufac"cm")/(2*Dy)))
	
	# origin of coordinate system in the center of the plane
	#x = range(-wi/2,wi/2,length=nx)
	#y = range(-wi/2,wi/2,length=ny)

	#itp = Interpolations.interpolate((x,y), M_, Gridded(Linear()))
	x = range(-wi/2,wi/2,length=nx_dom)
	y = range(-wi/2,wi/2,length=ny_dom)

	itp = Interpolations.interpolate((x,y), M__, Gridded(Linear()))

	#M_,itp	
	M__,itp	
end

# ╔═╡ f9ba467a-cefd-4d7d-829d-0889fc6d0f5e
begin
	M12, itp12 = sel12by12(FluxMap)
end;

# ╔═╡ f4ebb596-824a-4124-afb7-c368c1cabb00
md"""
Mean irradiation Flux on __catalyst surface__ (12 cm x 12 cm) : $(7.0/8.0*round(sum(M12)*(0.031497*ufac"cm")^2 / (12*12*ufac"cm^2") /ufac"kW/m^2"))
"""

# ╔═╡ 3aef8203-ce28-4197-a43e-784840f7bc1e
md"""
### Top chamber
"""

# ╔═╡ 2bf61cf4-2f4a-4c48-ad96-5de82403f3a0
md"""
__Edit 05.09.23__:
The most recent experiments involved catalyst, that was directly deposited onto the porous frit, covering the complete area. Thus there is no "uncoated frit" region exposed to the irradiation. The complete top surface is covered with deposited catalyst and can be treated as homogeneous.
"""

# ╔═╡ d9dee38e-6036-46b8-bc06-e545baa06789
md"""
Account for the exchange of irradiation between the surfaces in the top and bottom chambers of the reactor. 
In upper chamber, the surfaces are designated as (see figure below):
1. quartz window
2. catalyst layer


In the following, the net irradiation fluxes through the surfaces of interest marked by the __red__ dashed line (__catalyst layer__, __2__) are stated. Per sign convention applied here, __fluxes entering__ the control surfaces are counted as __positive__ while __fluxes exiting__ the control surfaces are counted as __negative__. 
These irradiation fluxes through the surfaces will be implemented in the code as boundary conditions for the thermal energy transport equation.

"""

# ╔═╡ e58ec04f-023a-4e00-98b8-f9ae85ca506f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""

$(LocalResource("../img/irrad_top_chamber_0923.png"))
```math
\begin{align}
F_2 &=-\epsilon_2\sigma T_2^4 + \alpha_2^{\text{vis}} G_1^{\text{vis}} + \alpha_2^{\text{IR}} G_1^{\text{IR}} \\
\end{align}
```
"""
  ╠═╡ =#

# ╔═╡ 80d2b5db-792a-42f9-b9c2-91d0e18cfcfb
md"""
For catalyst boundary region (2):
```math
G_1^{\text{vis}} = \frac{\tau_1^{\text{vis}} G_{\text{HFSS}}}{1-\rho_1^{\text{vis}}  \rho_2^{\text{vis}}  } \\
```
"""

# ╔═╡ 6de46bcb-9aff-4a21-b8d6-f14e64aac95c
md"""
```math
G_1^{\text{IR}} = \epsilon_1\sigma T_1^4
```
"""

# ╔═╡ b1ef2a89-27db-4f21-a2e3-fd6356c394da
md"""
where $G_1^{\text{vis}}$ and $G_1^{\text{IR}}$ are the surface brightnesses (cumulative irradiation flux resulting from emission, transmission and reflection) of surface 1 in the visible and IR spectral range respectively. The temperature of the quartz window $T_1$ is calculated as part of the problem as a "boundary species".
"""

# ╔═╡ 79e6cb48-53e2-4655-9d77-51f3deb67b94
function radiosity_window(f,u,bnode,data)
    (;iT,iTw,FluxIntp,FluxEmbed,uc_window,uc_cat,uc_frit)=data
    # irradiation exchange between quartz window (1), cat surface (2) 
    # window properties (1)
    tau1_vis=uc_window.tau_vis
    rho1_vis=uc_window.rho_vis
    tau1_IR=uc_window.tau_IR
    rho1_IR=uc_window.rho_IR
    eps1=uc_window.eps

    # obtain local irradiation flux value from interpolation + embedding	
    @views y,x,_ = bnode.coord[:,bnode.index] 
    Glamp =FluxEmbed*FluxIntp(x,y)

	# local tempererature of quartz window
    Tglass = u[iTw]
    G1_bot_IR = eps1*ph"σ"*Tglass^4
    if bnode.region==Γ_top_cat
        # catalyst layer (2)
        rho2_vis=uc_cat.rho_vis
		
        # vis radiosity of quartz window inwards / towards catalyst 
        #G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*rho2_vis)
		# flux profile measured behind quarz in plane of cat layer
		G1_bot_vis = Glamp

    end
    return G1_bot_vis,G1_bot_IR
end

# ╔═╡ 4dae93b2-be63-4ee9-bc1e-871a31ade811
md"""
The optical properties $\tau, \rho, \alpha, \epsilon$ are effective values, integrated over the respective spectral distributions of irradiative power from the lamp (visible, vis) and from emitted thermal radition (infrared, IR) for the surfaces 1, 2 and 3 respectively.
"""

# ╔═╡ 9547ed7c-3304-4c63-a6c1-5f84e0001c54
Base.@kwdef struct SurfaceOpticalProps
	alpha_IR::Float64 = 0.2 # absorptance: obtain from datasheet/measurement
	tau_IR::Float64 = 0.5 # transmittance: obtain from datasheet/measurement
	rho_IR::Float64 = 1.0 - tau_IR - alpha_IR
	eps::Float64 = alpha_IR
	# effective opt. par. integrated over visible spectrum (vis)
	alpha_vis::Float64 = 0.0
	tau_vis::Float64 = 0.9 # obtain from datasheet, transmittance
	rho_vis::Float64 = 1.0 - tau_vis - alpha_vis 
end;

# ╔═╡ 9d191a3a-e096-4ad7-aae6-bbd63d478fa2
# optical parameters for quartz window in upper chamber (uc)
const uc_window = SurfaceOpticalProps(
	# quartz glass windows
	#alpha_IR=0.2, # datasheet integrated over spectrum of HFSS
	#tau_IR=0.5,# datasheet integrated over spectrum of HFSS
	#########!!! introduced change here!!!##########
	alpha_IR=0.3, # datasheet integrated over spectrum of HFSS
	tau_IR=0.7,# datasheet integrated over spectrum of HFSS
	#########!!! introduced change here!!!##########
	alpha_vis=0.0, # quartz glass is transparent in visible spectrum
	tau_vis=0.9, # datasheet, take reflection losses at adouble interface into account
	# the calculated value for rho_vis is probably inaccurate since it depends on incidence angle, which differs between the irradiation coming from the lamp (relevant for tau_vis) and the diffuse reflected light coming from the catalyst 
)

# ╔═╡ 2d51f54b-4cff-4253-b17f-217e3261f36d
#  optical parameters for catalyst layer in upper chamber (uc)
const uc_cat = SurfaceOpticalProps(
	alpha_IR=0.45, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.45, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# ╔═╡ fdaeafb6-e29f-41f8-8468-9c2b03b9eed7
# optical parameters for uncoated frit in upper chamber (uc)
const uc_frit = SurfaceOpticalProps(
	alpha_IR=0.2, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.2, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# ╔═╡ 2e631f58-6a4b-4c6d-86ad-b748fd2d463a
md"""
#### Window Temperature
Calculate window temperature that is used in the temperature boundary condition of the top chamber. It can be obtained from an energy balance in thermal equilibrium. The value of window temperature is obtained as part of the problem via a "boundary species".
"""

# ╔═╡ a395374d-5897-4764-8499-0ebe7f2b4239
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
$(LocalResource("../img/WindowEB_0923.png")) 
"""
  ╠═╡ =#

# ╔═╡ 7db215d7-420e-435b-8889-43c4fff150e0
md"""
At steady state, the net transfer of energy to/from the balance region around the quartz window must vanish. This leads to:

```math
0=-\dot q_{\text{conv,top}} - 2 \dot q_{\text{emit}} + \dot q_{\text{conv,bot,2}} + \dot q_{\text{abs,bot,2}}
```

with
```math
	\begin{align}
	\dot q_{\text{conv,top}} &= \alpha_{\text{conv}} (T_1-T_{\text{amb}}) \\
	\dot q_{\text{emit}} &= \epsilon_1\sigma T_1^4 \quad(\text{once for top \& bottom surfaces}) \\
	\dot q_{\text{conv,bot,2}} &= \alpha_{\text{inner. conv}}(T_{\text m}-T_1), \quad T_{\text m}=\frac{T_1+T_2}{2} \\
	\dot q_{\text{abs,bot,2}} &= \alpha_1^{\text{IR}} G_2^{\text{IR}} \\
	\end{align}
```

From this balance, the steady-state window temperature $T_1$ can be obtained.
"""

# ╔═╡ 8b485112-6b1a-4b38-af91-deb9b79527e0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Bottom chamber
$(LocalResource("../img/irrad_bottom_chamber.png"))
```math
F = -\epsilon_1 \sigma T_1^4 + \frac{\alpha_1^{\text{IR}}}{1-\rho_1^{\text{IR}} \rho_2^{\text{IR}}} \left( \epsilon_2 \sigma T_2^4 + \rho_2^{\text{IR}} \epsilon_1 \sigma T_1^4 \right)
```
"""
  ╠═╡ =#

# ╔═╡ 25edaf7b-4051-4934-b4ad-a4655698a6c7
md"""
where the temperature of the bottom Aluminium plate $T_2$ is calculated as a boundary species as part of the equation system.
"""

# ╔═╡ 3fe2135d-9866-4367-8faa-56cdb42af7ed
md"""

For opaque surfaces $1=\alpha + \rho$ and with Kirchhoff's law $\alpha = \epsilon$, we can simplify the term for the net radiative transfer rate to (positive sign convention) surface __1__:
```math
F = \frac{ \sigma (T_2^4 -T_1^4)}{\frac{1}{\epsilon_1}+\frac{1}{\epsilon_2}-1}
```
"""

# ╔═╡ 652497ee-d07b-45e2-aeaf-87ad5bcc23ad
md"""
#### Bottom plate temperature
Calculate bottom plate temperature $T_2$ that is used in the temperature boundary condition of the bottom chamber. It can be obtained from an energy balance in thermal equilibrium. The value of $T_2$ is stored in the data struct and is updated after every simulation run.
"""

# ╔═╡ 6bd59a54-f059-4646-b053-0fa41ead87fd
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
$(LocalResource("../img/BotChamberEB.png")) 
"""
  ╠═╡ =#

# ╔═╡ f5d78670-a98b-46a1-8bf3-3d2599cfdd88
md"""
At steady state, the net transfer of energy to/from the balance region around the bottom plate must vanish. This leads to:

```math
0=-\dot q^{\text{nat}}_{\text{conv}} - 2 \dot q_{\text{emit}} - \dot q^{\text{int}}_{\text{conv}} + \dot q_{\text{abs,1}} 
```

with
```math
	\begin{align}
	\dot q^{\text{nat}}_{\text{conv}} &= k^{\text{nat}}_{\text{conv}} (T_2-T_{\text{amb}}) \\

	\dot q_{\text{emit}} &= \epsilon_2\sigma T_2^4 \quad(\text{once for top \& bottom surfaces}) \\
	\dot q^{\text{int}}_{\text{conv}} &= k^{\text{int}}_{\text{conv}}(T_2-T_{\text m}), \quad T_{\text m}=\frac{T_1+T_2}{2} \\
	
	\dot q_{\text{abs,1}} &= \alpha_2^{\text{IR}} G_1^{\text{IR}} \\
	
	\end{align}
```

The surface brigthness of the frit, that solely lies in the IR spectral range $G_1^{\text{IR}}$ is computed via:

```math
	G_1^{\text{IR}} = \frac{\epsilon_1 \sigma T_1^4 + \rho_1^{\text{IR}} \epsilon_2 \sigma T_2^4}{1-\rho_1^{\text{IR}} \rho_2^{\text{IR}}}

```

We obtain the mean temperature of the frit $T_1$ from the simulation. From this balance, the steady-state bottom plate temperature $T_2$ can be obtained.
"""

# ╔═╡ 162122bc-12ae-4a81-8df6-86498041be40
# optical parameters for uncoated frit in lower chamber (lc) = frit in upper chamber
const lc_frit = uc_frit

# ╔═╡ 221a1ee4-f7e9-4233-b45c-a715c9edae5f
# optical parameters for Al bottom plate in lower chamber (lc)
const lc_plate = SurfaceOpticalProps(
	# aluminium bottom plate (machined surface)
	alpha_IR=0.1, # see in lit, assume large reflectivity
	tau_IR=0.0, # opaque surface
	alpha_vis=0.1, # see in lit, assume large reflectivity
	tau_vis=0.0 # opaque surface	
)

# ╔═╡ 2228bbd4-bc84-4617-a837-2bf9bba76793
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Convective Transport
Apply convective heat transport through the top and bottom chambers: the heat transport occurs perpendicular to the main direction of flow according to methodology for __Parallel Plate Ducts__ from VDI heat atlas __ch. G2, sec. 8, p. 707__. For this case the dimensionless numbers are defined as:

```math
\text{Nu}=\frac{k^{\text{int}}_{\text{conv}} d_{\text h}}{\lambda} \qquad \text{Re}=\frac{w d_{\text h}}{\nu}
```
The characteristic dimension for flow through the gap, corresponding to the
hydraulic diameter, is twice the gap width:

```math
d_{\text h}=2 s
```

Assume the asymptotic situation of Hydrodynamically and thermally fully developed laminar flow for a fluid flowing between a heated and an isolated surface:

$(LocalResource("../img/ParallelPlateDuct.png", :width => 200)) 
This considers the convective heat transfer between the fluid in the chamber and a single surface.
This situation is addressed by eq. (42) from VDI heat atlas ch. G2, sec. 8, p. 708:

```math
\text{Nu} = 4.861
```
"""
  ╠═╡ =#

# ╔═╡ 09d976a7-f6c6-465b-86f4-9bc654ae158c
md"""
The convective heat flux from the bottom wall into the fluid is expressed via:
```math
\dot q_{\text{bot}} = k^{\text{int}}_{\text{conv}}(T_{\text{bot}} - T_{\text{m}})
```
Similarly, the heat flux from the fluid towards the top wall, which is at a lower temperature is obtained:
```math
\dot q_{\text{top}} = k^{\text{int}}_{\text{conv}}(T_{\text{m}}-T_{\text{top}})
```
where $T_{\text{m}}=(T_{\text{top}}+T_{\text{bot}})/2$.
"""

# ╔═╡ 5c9aad13-914a-4f5e-af32-a9c6403c52d0
md"""
## Side Walls
The heat flux through the side walls of the modelling domain $\dot q_{\text{side}}$ is limited by the convective heat flux on the outside of the reactor shell. In absence of active cooling e.g. via a fan, natural convection prevails with heat transfer coefficient $k^{\text{nat}}_{\text{conv}}$. At the contacting interface between the porous frit and the reactor shell made of Aluminium at the domain boundary a negligable contact resistance is assumed.
"""

# ╔═╡ ff764eea-e282-4fed-90b4-5f418ae426f0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
$(LocalResource("../img/SideBCConv_gap.png", :width=>400)) 
"""
  ╠═╡ =#

# ╔═╡ dfed01d4-8d88-4ce7-afe5-1fd27e2cc746
md"""
```math
\begin{align}
\dot q_{\text{side}}=\frac{\dot Q}{A_{\text{frit}}}&=\frac{1}{A_{\text{frit}}}\frac{1}{W_{\text{cond}}^{\text{gap}}+W^{\text{wall}}_{\text{cond}}+W_{\text{conv}}}(T-T_{\text{amb}})\\
&=\frac{1}{A_{\text{frit}}}\frac{1}{\frac{\delta_{\text{gap}}}{A_{\text{frit}} \lambda_{\text{gap}}}+\frac{\delta_{\text{wall}}}{A_{\text{frit}} \lambda_{\text{wall}}}+\frac{1}{A_{\text{frit}} k^{\text{nat}}_{\text{conv}}}}(T-T_{\text{amb}}), \quad \text{with} \frac{\delta_{\text{wall}}}{\lambda_{\text{wall}}} \ll \frac{\delta_{\text{gap}}}{\lambda_{\text{gap}}} \approx \frac{1}{k^{\text{nat}}_{\text{conv}}}\\
&\approx \frac{T-T_{\text{amb}}}{\frac{\delta_{\text{gap}}}{\lambda_{\text{gap}}}+  \frac{1}{k^{\text{nat}}_{\text{conv}}}}\\
\end{align}
```
"""

# ╔═╡ 4ebbe06f-0993-4c5c-9af3-76b2b645e592
md"""
## Implementation
"""

# ╔═╡ dec81e61-d1eb-4a6c-8648-a07b98a8a7a9
function bflux(f, u, bedge, data)
	# window temperature distribution
	#if bedge.region == Γ_top_cat || bedge.region == Γ_top_frit
	if bedge.region == Γ_top_cat
		(;iTw,λ_window) = data
		f[iTw] = λ_window * (u[iTw, 1] - u[iTw, 2])
	end
	
	if bedge.region == Γ_bottom
		(;iTp,λ_Al) = data
		f[iTp] = λ_Al * (u[iTp, 1] - u[iTp, 2])
	end
end

# ╔═╡ 1a8410fa-caac-4023-8c96-a0042f22d88c
function bottom_rim(f,u,bnode,data)
	(;iT,iTp,k_nat_conv,Tamb)=data

	if bnode.region==Γ_bottom_rim
		f[iT] = k_nat_conv*(u[iTp]-Tamb)		
	end
end

# ╔═╡ 02b76cda-ffae-4243-ab40-8d0fe1325776
md"""
# Auxiliary functions
"""

# ╔═╡ 44aa5b49-d595-4982-bbc8-100d2f199415
md"""
# System Setup and Solution
"""

# ╔═╡ e25e7b7b-47b3-457c-995b-b2ee4a87710a
md"""
## Model Data
"""

# ╔═╡ 3a35ac76-e1b7-458d-90b7-d59ba4f43367
begin
	
	Base.@kwdef mutable struct ModelData{NG} <:AbstractModelData
	
	# catalyst / chemistry data
	# kinetic parameters, S3P="simple 3 parameter" kinetics fit to UPV lab scale experimental data

	#kinpar::FixedBed.KinData{nreac(S3P)} = S3P
	kinpar::FixedBed.KinData{nreac(XuFroment)} = XuFroment
	#kinpar::FixedBed.KinData{nreac(Wolf_rWGS)} = Wolf_rWGS
		
	#ng::Int64		 		= kinpar.ng # number of gas phase species
	ip::Int64=NG+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable
	# register window & plate temperatures as boundary species
	iTw::Int64=iT+1 # index of window Temperature (upper chamber)
	iTp::Int64=iTw+1 # index of plate Temperature (lower chamber)

	
	gn::Dict{Int, Symbol} 	= kinpar.gn # names and fluid indices
	gni::Dict{Symbol, Int}  = kinpar.gni # inverse names and fluid indices
	Fluids::Vector{FluidProps} = kinpar.Fluids # fluids and respective properties in system
	
	X0::Vector{Float64} = let
		x=zeros(Float64, NG)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		x/sum(x)
	end # inlet composition

	mcat::Float64=500.0*ufac"mg"
	# fit factor: increase cat mass to match exp reactivity
	#mcat::Float64=625.0*ufac"mg"
	    
	# volume specific cat mass loading, UPV lab scale PC reactor
	lcats::Float64 =1000.0*ufac"kg/m^3"
	#mcats::Float64=80.0*ufac"kg/m^3" # 200 mg cat total loading, 250μm CL 
	isreactive::Bool = 1
	#isreactive::Bool = 0
		
	# natural convection heat transfer coefficient	
	k_nat_conv::Float64=17.5*ufac"W/(m^2*K)" 
	
	## porous filter data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0

	# frit thickness (applies to 2D & 3D)
	h::Float64=0.5*ufac"cm"
	# catalyst layer thickness (applies to 2D & 3D)
	#cath::Float64 = 500.0*ufac"μm"
	cath::Float64 = 250.0*ufac"μm"

	# upper and lower chamber heights for calculation of conduction/convection b.c.
	uc_h::Float64=17.0*ufac"mm"
	lc_h::Float64=18.0*ufac"mm"
	Nu::Float64=4.861

	# prism / 3D
	#wi::Float64=15.7*ufac"cm" # porous frit width/side lenght
	wi::Float64=16.0*ufac"cm" # simplified porous frit width/side lenght
	le::Float64=wi # prism width/side lenght
	catwi::Float64=10.0*ufac"cm" # prism width/side lenght	
		
	shellh::Float64=3.6*ufac"cm" # height of reactor shell contacting domain
	#shellh::Float64=h+cath # height of reactor shell contacting domain
	delta_gap::Float64=1.5*ufac"mm" # gas gap between frit and reactor wall
		
	Ac::Float64=wi*le*ufac"m^2" # cross-sectional area, square

	## irradiation data
	#Glamp_target::Float64=100.0*ufac"kW/m^2" # solar simulator irradiation flux
	#Glamp::Float64=1.0*ufac"kW/m^2" # solar simulator irradiation flux
	FluxIntp::typeof(itp12)=itp12 # interpolator for irradiation flux
	FluxEmbed::Float64=0.0 # "persistent" embedding parameter avail. outside of solve call w/ embedding
		
    #Flux_target::Float64=1.0
	Flux_target::Float64=7.0/8.0 # interpolate 80 sun flux map down to 70 suns
	
	# upper chamber: quartz window eff. optical properties
	uc_window::SurfaceOpticalProps = uc_window
	# upper chamber: catalyst layer eff. optical properties
	uc_cat::SurfaceOpticalProps = uc_cat
	# upper chamber: (uncoated) porous frir eff. optical properties
	uc_frit::SurfaceOpticalProps = uc_frit


	# lower chamber: (uncoated) porous frir eff. optical properties
	lc_frit::SurfaceOpticalProps = lc_frit
	# lower chamber:  Al bottom plate eff. optical properties
	lc_plate::SurfaceOpticalProps = lc_plate
	
	# Solid Boro-Silikatglas
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	λs::Float64=1.13*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
	cs::Float64=0.8e3*ufac"J/(kg*K)" # heat capacity of non-porous SiO2
	
	ϕ::Float64=0.36 # porosity, exp determined
	#ϕ::Float64=0.33 # porosity, VitraPor sintetered filter class 0
	ψ::Float64=0.26 # flattening coefficient in 
	
	# approximation from Wesselingh, J. A., & Krishna, R. (2006). Mass Transfer in Multicomponent Mixtures
	γ_τ::Float64=ϕ^1.5 # constriction/tourtuosity factor

	#k::Float64=2.9e-11*ufac"m^2" # permeability , por class 2
	k::Float64=1.23e-10*ufac"m^2" # permeability , por class 0
	
	# a_s::Float64=0.13*ufac"m^2/g" # specific surface area, por class 2
	a_s::Float64=0.02*ufac"m^2/g" # specific surface area, por class 0

	
	ρfrit::Float64=(1.0-ϕ)*ρs*ufac"kg/m^3" # density of porous frit
	a_v::Float64=a_s*ρfrit*ufac"m^2/m^3" # volume specific interface area
	## END porous filter data

	# quartz window / Al bottom plate thermal conductivities
	λ_window::Float64=1.38*ufac"W/(m*K)"
	λ_Al::Float64=235.0*ufac"W/(m*K)"

	## Flow data
	#norm conditions
	#pn::Float64 = 1.0*ufac"bar"
	pn::Float64 = 1.0*ufac"atm"
	#Tn::Float64 = 273.15*ufac"K"
	Tn::Float64 = 293.15*ufac"K"

	# volumetric feed flow rate, at norm conditions defined above
	#Qflow::Float64=0.0*ufac"ml/minute" 
	#Qflow::Float64=1480.0*ufac"ml/minute" # volumetric feed flow rate
	Qflow::Float64=3000.0*ufac"ml/minute" # volumetric feed flow rate
	#Qflow::Float64=14800.0*ufac"ml/minute" # volumetric feed flow rate
		

	MWin::Float64 = molarweight_mix(Fluids, X0)
	mdotin::Float64=MWin*Qflow*pn/(ph"R"*Tn)*ufac"kg/s"
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	#p::Float64=1.0*ufac"atm" # reactor pressure
	p::Float64=3.0*ufac"bar" # reactor pressure

	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	ubot::Float64=u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate
	
end;
	
	# !!!ALLOC Method to be called instead of data.ng
	FixedBed.ngas(::ModelData{NG}) where NG = NG
	
	# !!!ALLOC Additional constructo taking ng as parameter	
	ModelData(;ng=S3P.ng, kwargs...) = ModelData{ng}(;kwargs...)

end;

# ╔═╡ 985718e8-7ed7-4c5a-aa13-29462e52d709
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Cutplane at ``z=`` $(d=ModelData{S3P.ng}(); @bind zcut PlutoUI.Slider(range(0.0,(d.h)/ufac"mm",length=101),default=(d.h)/ufac"mm",show_value=true)) mm
"""
  ╠═╡ =#

# ╔═╡ ff58b0b8-2519-430e-8343-af9a5adcb135
# ╠═╡ skip_as_script = true
#=╠═╡
let
	gridplot(grid_fun(ModelData{S3P.ng}(),nref=0),zoom=1.2,resolution=(650,600), zplane=zcut*ufac"mm",Plotter=PlutoVista)
end
  ╠═╡ =#

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	(;ip,isreactive)=data
	ng=ngas(data)
	
	if node.region == 2 && isreactive # catalyst layer
		(;iT,lcats,kinpar,Fluids,gni)=data
		(;rni,nuij)=kinpar
		pi = MVector{ngas(data),eltype(u)}(undef)		

		for i=1:ng
            pi[i] = u[i]
		end

        RR = @inline -lcats*ri(data,u[iT],pi)
		for i=1:ng
			f[i] = zero(eltype(u))
			for j=1:nreac(kinpar)
				#f[i] += nuij[i,j] * RR[j]
                f[i] += nuij[(j-1)*ng+i] * RR[j]
			end			
		end

	end
	# ∑xi = 1
	@views f[ip]=u[ip]-sum(u[1:ng])
end

# ╔═╡ e459bbff-e065-49e1-b91b-119add7d4a71
let	
	(;wi,Flux_target)=ModelData()
	nx,ny=size(M12)
	x_ = range(-wi/2,wi/2,length=nx)
	y_ = range(-wi/2,wi/2,length=ny)
	
	Plots.plot(xguide="X Coordinate / cm", yguide="Y Coordinate / cm", title="Measured Flux Profile", colorbar_title="Irradiation Flux / kW m-2", xlim=(-wi/2/ufac"cm", wi/2/ufac"cm"), ylim=(-wi/2/ufac"cm", wi/2/ufac"cm"),size=(400,400), xticks=(-wi/2/ufac"cm"):2.0:(wi/2/ufac"cm"), yticks=(-wi/2/ufac"cm"):2.0:(wi/2/ufac"cm") )
	#p=Plots.plot(xguide="X Coordinate / cm", yguide="Y Coordinate / cm", xlim=(-wi/2/ufac"cm", wi/2/ufac"cm"), ylim=(-wi/2/ufac"cm", wi/2/ufac"cm"), size=(300,300))
	Plots.contour!(x_./ufac"cm",y_./ufac"cm",Flux_target*M12 / ufac"kW/m^2", lw=0, aspect_ratio = 1, fill = true)
	#Plots.plot!([-5,5],[-5,-5],label=:none,c=:black,lw=2)
	#Plots.plot!([-5,5],[5,5],label=:none,c=:black,lw=2)
	#Plots.plot!([-5,-5],[-5,5],label=:none,c=:black,lw=2)
	#Plots.plot!([5,5],[-5,5],label=:none,c=:black,lw=2)
	#Plots.savefig("../img/out/Flux_meas_70suns.svg")
end

# ╔═╡ 7da59e27-62b9-4b89-b315-d88a4fd34f56
function top(f,u,bnode,data)
	# top boundaries (cat layer & frit)
	if bnode.region==Γ_top_cat
		(;iT,iTw,FluxIntp,FluxEmbed,X0,Fluids,Tamb,uc_window,uc_cat,uc_frit,uc_h,Nu,k_nat_conv)=data
		ng=ngas(data)
		# convective therm. energy flux
		flux_entahlpy=zero(eltype(u))
		
		for i=1:ng

			f[i] = -data.u0*data.X0[i]*data.pn/(ph"R"*data.Tn)
			
			# use species enthalpy (incl. Δh_formation) for conv. therm. eng. flux 
			# reactant gases enter control volume at Tamb
			@inline flux_entahlpy += f[i] * enthalpy_gas(Fluids[i], Tamb)
		end

		
		flux_irrad=zero(eltype(u))
		# calculate local absortion in window
		flux_abs_w=zero(eltype(u))
		σ=ph"σ"

		# irradiation exchange between quartz window (1), cat surface (2)
		# window properties (1)
		alpha1_vis=uc_window.alpha_vis
        alpha1_IR=uc_window.alpha_IR
		# catalyst layer properties (2)
		rho2_vis=uc_cat.rho_vis
        rho2_IR=uc_cat.rho_IR
		alpha2_vis=uc_cat.alpha_vis
		alpha2_IR=uc_cat.alpha_IR
		eps2=uc_cat.eps		

		G1_vis, G1_IR = radiosity_window(f,u,bnode,data)
		
		if bnode.region==Γ_top_cat
			flux_irrad = -eps2*σ*u[iT]^4 + alpha2_vis*G1_vis + alpha2_IR*G1_IR
			G2_vis = rho2_vis*G1_vis		
            G2_IR = eps2*ph"σ"*u[iT]^4 + rho2_IR*G1_IR
			flux_abs_w = alpha1_vis*G2_vis + alpha1_IR*G2_IR
		end

		Tglass = u[iTw]
		# mean temperature
        Tm=0.5*(u[iT] + Tglass)
        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X0)

		# convective heat flux through top chamber
		dh=2*uc_h
		kconv=Nu*λf/dh*ufac"W/(m^2*K)"
		flux_conv = kconv*(u[iT]-Tm)

		f[iT] = -flux_irrad + flux_conv + flux_entahlpy

		# calculate local window temperature from (local) flux balance
		flux_conv_top_w = k_nat_conv*(u[iTw]-Tamb)
		flux_emit_w = uc_window.eps*ph"σ"*u[iTw]^4
		f[iTw] = -flux_conv -flux_abs_w +flux_conv_top_w +2*flux_emit_w
	end
end

# ╔═╡ c0de2aff-8f7f-439c-b931-8eb8fbfcd45d
function mole_frac!(node,data,X,u::VoronoiFVM.EdgeUnknowns)
	ng=ngas(data)
	sump=zero(eltype(u))
	#X=zeros(eltype(u), ng)
	for i=1:ng
		X[i] = 0.5*(u[i,1]+u[i,2])
		sump += X[i]
	end
	for i=1:ng
		X[i] = X[i] / sump
	end
	nothing
end

# ╔═╡ 51b03d3b-3c7c-4019-8ed8-bf1aaa0b1ddb
function mole_frac!(node,data,X,u::VoronoiFVM.BNodeUnknowns)
	ng=ngas(data)
	sump=zero(eltype(u))
	#X=zeros(eltype(u), ng)
	for i=1:ng
		X[i] = u[i]
		sump += X[i]
	end
	for i=1:ng
		X[i] = X[i] / sump
	end
	nothing
end

# ╔═╡ 40906795-a4dd-4e4a-a62e-91b4639a48fa
function bottom(f,u,bnode,data)
	if bnode.region==Γ_bottom # bottom boundary
		(;iT,ip,iTp,Fluids,ubot,Tamb,k_nat_conv,lc_frit,lc_plate,lc_h,Nu) = data
		ng=ngas(data)

		X=MVector{ngas(data),eltype(u)}(undef)
		@inline mole_frac!(bnode,data,X,u)
		

		# use species enthalpy (incl. Δh_formation) for conv. therm. eng. flux 1/2
		@inline hmix=enthalpy_mix(Fluids, u[iT], X)
		flux_entahlpy=ubot*u[ip]/(ph"R"*u[iT])*hmix

		# irradiation exchange between porous frit (1) and Al bottom plate (2)
		# porous frit properties (1)
		rho1_IR=lc_frit.rho_IR
		alpha1_IR=lc_frit.alpha_IR
		eps1=lc_frit.eps
		# Al bottom plate properties (2)
		rho2_IR=lc_plate.rho_IR
		alpha2_IR=lc_plate.alpha_IR
		eps2=lc_plate.eps
		#(;eps1,alpha1_IR,rho1_IR,rho2_IR,eps2,rho2_IR) = lower_chamber
		σ=ph"σ"

		Tplate = u[iTp]
		flux_irrad = -eps1*σ*u[iT]^4 + alpha1_IR/(1-rho1_IR*rho2_IR)*(eps2*σ*Tplate^4+rho2_IR*eps1*σ*u[iT]^4)

		# conductive heat flux through top chamber
		# mean temperature
        Tm=0.5*(u[iT] + Tplate)
        # thermal conductivity at Tm and outlet composition X
        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X)

        # flux_cond = -λf*(u[iT]-Tplate)/lc_h # positive flux in positive z coord.
		dh=2*lc_h
		kconv=Nu*λf/dh*ufac"W/(m^2*K)"
		flux_conv = kconv*(u[iT]-Tm) # positive flux in negative z coord. (towards plate)
		
		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative
		# f[iT] = -flux_irrad + flux_entahlpy - flux_cond
		f[iT] = -flux_irrad + flux_entahlpy + flux_conv		
		
		for i=1:ng
			# specify flux at boundary: flow velocity is normal to bot boundary
			f[i] = ubot*u[i]/(ph"R"*u[iT])			
		end

		# calculate (local) plate temperature from (local) flux balance
		flux_conv_bot_p = k_nat_conv*(u[iTp]-Tamb)
		flux_emit_p = lc_plate.eps*ph"σ"*u[iTp]^4
		G1_IR = (eps1*ph"σ"*u[iT]^4 + rho1_IR*eps2*ph"σ"*u[iTp]^4)/(1-rho1_IR*rho2_IR)
		flux_abs_p = alpha2_IR*G1_IR
		f[iTp] =  -flux_conv -flux_abs_p +flux_conv_bot_p +2*flux_emit_p
	end
end

# ╔═╡ edd9fdd1-a9c4-4f45-8f63-9681717d417f
function side(f,u,bnode,data)
	# side wall boundary condition
	(;iT,h,cath,shellh,k_nat_conv,delta_gap,Tamb)=data

	# all sides for complete domain
	if bnode.region==Γ_side_front || bnode.region==Γ_side_right || bnode.region==Γ_side_back || bnode.region==Γ_side_left

		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative

		X=MVector{ngas(data),eltype(u)}(undef)
		@inline mole_frac!(bnode,data,X,u)
		@inline _,λf=dynvisc_thermcond_mix(data, u[iT], X)

		# w/o shell height
		f[iT] = (u[iT]-Tamb)/(delta_gap/λf+1/k_nat_conv)
		
	end
end

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
	top(f,u,bnode,data)
	bottom(f,u,bnode,data)
	side(f,u,bnode,data)
	bottom_rim(f,u,bnode,data)
	# bottom rim
	#(;iT)=data
 	#boundary_dirichlet!(f,u,bnode,iT,Γ_bottom_rim,50+273.15)
end

# ╔═╡ 2191bece-e186-4d8e-8a21-3830441baf11
function D_matrix!(data, D, T, p)
	ng=ngas(data)
	for i=1:(ng-1)
		for j=(i+1):ng
			Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
			D[j,i] = Dji
			D[i,j] = Dji
			#v[j,i] = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
		end
	end
	nothing
end

# ╔═╡ b6381008-0280-404c-a86c-9c9c3c9f82eb
function M_matrix!(M,D,T,p,x,data)
	ng=ngas(data)
	# !!!ALLOC all methods to be called with arrays to be stack allocated
	# have to  be inlined - here we use callsite inline from Julia 1.8
	@inline D_matrix!(data, D, T, p)
	for i=1:ng
		M[i,i] = 1/DK_eff(data,T,i)
		for j=1:ng
			if j != i
				M[i,i] += x[j]/(data.γ_τ*D[i,j])
				M[i,j] = -x[i]/(data.γ_τ*D[i,j])
			end
		end	
	end
	nothing
end

# ╔═╡ ed7941c4-0485-4d84-ad5b-383eb5cae70a
function flux(f,u,edge,data)
	(;ip, iT, k, Fluids) = data

	ng=ngas(data)
	
	F = MVector{ngas(data),eltype(u)}(undef)
	X = MVector{ngas(data),eltype(u)}(undef)

	pk,pl = u[ip,1],u[ip,2]
	δp = pk-pl
	pm = 0.5*(pk+pl)
	Tm=0.5*(u[iT,1]+u[iT,2])
        
	@inline mole_frac!(edge,data,X,u)
	@inline cf=heatcap_mix(Fluids, Tm, X)
	@inline μ,λf=dynvisc_thermcond_mix(data, Tm, X)
 	λbed=kbed(data,λf)*λf

	# Darcy flow
	ud=-k/μ * δp

	# !!!ALLOC Use MMatriy with static size information instef of Matrix
	M = MMatrix{ngas(data),ngas(data),eltype(u)}(undef)
	D = MMatrix{ngas(data),ngas(data),eltype(u)}(undef)
	@inline M_matrix!(M,D,Tm,pm,X,data)
   	
		
	for i=1:ng
		DK = DK_eff(data,Tm,i)
		bp,bm=fbernoulli_pm(ud/DK)		
		F[i] = (bm*u[i,1]-bp*u[i,2])/(ph"R"*Tm)
	end
	
	# computation of fluxes J

	@inline inplace_linsolve!(M,F)
	
	@views f[1:ng] .= F	

	# use species enthalpy (incl. Δh_formation) for conv. therm. eng. flux 1/3
	conv=zero(eltype(u))
	for i=1:ng
		conv += F[i] * enthalpy_gas(Fluids[i], Tm)
	end	
	Bp,Bm = fbernoulli_pm(conv/λbed/Tm)

	f[iT]= λbed*(Bm*u[iT,1]-Bp*u[iT,2])
end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
function main(;data=ModelData(),nref=0,control = sys->SolverControl(),assembly=:cellwise )

	
	grid=grid_fun(data)
	(;iT,iTw,iTp,p,X0,Tamb)=data
	ng=ngas(data)
	
	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=flux,
							reaction=reaction,
							bcondition=bcond,
							bflux=bflux,
		                    assembly
							)
	enable_species!(sys; species=collect(1:(ng+2))) # gas phase species + p + T
	#enable_boundary_species!(sys, iTw, [Γ_top_cat,Γ_top_frit]) # window temperature as boundary species in upper chamber
	enable_boundary_species!(sys, iTw, [Γ_top_cat]) # window temperature as boundary species in upper chamber
	enable_boundary_species!(sys, iTp, [Γ_bottom]) # plate temperature as boundary species in lower chamber
	

	# this is not  good...
	#precon_linear=BlockPreconditioner(partitioning=partitioning(sys),factorization=UMFPACKFactorization())
	
	inival=unknowns(sys)
	inival[:,:].=1.0*p
	for i=1:ng
		inival[i,:] .*= X0[i]
	end
	inival[[iT,iTw,iTp],:] .= Tamb


	sol=solve(sys;inival,control=control(sys))
	
	function pre(sol,par)
		@info par
		ng=ngas(data)

		
		# iteratively adapt bottom outflow boundary condition
	 	function Intbot(f,u,bnode,data)			
	 		X=zeros(eltype(u), ng)
	 		mole_frac!(bnode,data,X,u)
	 		# top boundary(cat/frit)
	 		if bnode.region==Γ_bottom
	 			for i=1:ng
	 				f[i] = data.Fluids[i].MW*X[i]
	 			end
	 			f[iT] = u[iT]
	 		end
	 	end
		
	 	MWavg=sum(integrate(sys,Intbot,sol; boundary=true)[1:ng,Γ_bottom])/(data.Ac/4)
	 	ndotbot=data.mdotin/MWavg
		 
	 	Tavg=sum(integrate(sys,Intbot,sol; boundary=true)[iT,Γ_bottom])/(data.Ac/4)
	 	
		data.ubot = ndotbot*ph"R"*Tavg/(1.0*ufac"bar")/data.Ac
	 
				
	 	# calculate volume specific catalyst loading lcat in kg/ m^3
		(;catwi,cath,mcat,Flux_target)=data
		Vcat = catwi^2*cath # m^3
		lcat = mcat/Vcat # kg/m^3"
		lcats=[1.0, lcat]*ufac"kg/m^3"
	 	data.lcats = minimum(lcats) + par*(maximum(lcats)-minimum(lcats))

	 	# irradiation flux density
        FluxEmbeds=[0.05, Flux_target]
		data.FluxEmbed = minimum(FluxEmbeds) + par*(maximum(FluxEmbeds)-minimum(FluxEmbeds))
	 end


	
	 mycontrol=control(sys ;
	 				  		handle_exceptions=true,
							#Δp_min=1.0e-4,
	 						Δp=0.05,
	 				  		#Δp=0.25,
	 				  		#Δp_grow=2.0,
	 						Δp_grow=1.2,
	 				  		Δu_opt=1.0e5 )
	
	embed=[0.0,1.0]
	
	sol=solve(sys;control=mycontrol,inival,embed,pre,)

	
	sol(sol.t[end]),grid,sys,data
	#sol,grid,sys,data
end;

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	if RunSim

		function control(sys;kwargs...)
			SolverControl(
			gmres_umfpack(),
			sys;
			verbose="na",
			log=true,
			kwargs...)
		end
		
		sol_,grid,sys,data_embed=main(;data=ModelData(),nref=0,control,
		assembly=:edgewise);
		
		if sol_ isa VoronoiFVM.TransientSolution
			sol = copy(sol_(sol_.t[end]))
		else
			sol = copy(sol_)
		end
	end
end;
  ╠═╡ =#

# ╔═╡ bcaf83fb-f215-428d-9c84-f5b557fe143f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
To calculate the temperature field within the modelling domain, the thermal energy equation is solved. The modelling domain consists of the porous frit material (SiO2, silica-glas) and the gaseous phase in the void volume. The domain is considered a "Quasi-homogenous" phase with ``T_{\text{s}}=T_{\text{f}}`` (i.e. porous solid and gas phase are in thermal equil.) This assumption is based on:
- Very high specific surface area of porous material ``a_{\text V} = ``$(data_embed.a_v) ``\frac{\text m^2}{\text m^3}``
- Large interphase heat transfer coefficient (Kuwahara, F., Shirota, M., & Nakayama, A. (2001). doi:10.1016/s0017-9310(00)00166-6)


Heat is transported within the domain via conduction and convective transport. Because of the treatment of the porous domain as a Quasi-homogenous phase, an effective thermal conductivity is used to describe conduction. The convective heat transport considers both advective flow through the porous medium described by D'arcy's law and the multi-component diffusion flux combined in the total species flux $\vec N_i$.
"""
  ╠═╡ =#

# ╔═╡ ecb8df4c-fb61-472b-ad0b-e55b49905442
#=╠═╡
let
	(;FluxEmbed)=data_embed
	Dx = 0.031497*ufac"cm"
	Dy = 0.031497*ufac"cm"

	wi = 5.0*ufac"cm"
	
	nx = Integer(ceil(2*wi / Dx))
	ny = Integer(ceil(2*wi / Dy))

	
	
	xs=range(-wi, wi, length=nx)
	ys=range(-wi, wi, length=ny)

	int=0.0
	for x in xs
		for y in ys
			int += FluxEmbed*itp12(x,y)*(0.031497*ufac"cm")^2
		end
	end

	int/(2*wi)^2/ufac"kW" # kW/m^2
	#int/ufac"kW" # kW
end
  ╠═╡ =#

# ╔═╡ b2df1087-6628-4889-8cd6-c5ee7629cd93
md"""
## Temperature Plot
"""

# ╔═╡ 3bd80c19-0b49-43f6-9daa-0c87c2ea8093
#=╠═╡
let
	(;iT,ip,gni)=data_embed

	scalarplot(grid, sol[iT,:].- 273.15, levelalpha=0.8, colormap=:inferno, resolution=(650,600), zoom=1.2, Plotter=PlutoVista)	
end
  ╠═╡ =#

# ╔═╡ b99610da-ca0c-46e7-b0a7-bc7858a45a35
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;gni,ip)=data_embed
	#pCO = [p >= 0 ? p : 0.0 for p in sol[gni[:CO],:]]
	#writeVTK("../data/out/3D_ptot_70suns_3000sccm.vtu", grid; point_data = sol[ip,:])
	#writeVTK("../data/out/3D_pCO_70suns_3000sccm.vtu", grid; point_data = sol[gni[:CO],:])
	#writeVTK("../data/out/3D_pCO_pos_70suns_3000sccm.vtu", grid; point_data = pCO)
	#writeVTK("../data/out/3D_T_70suns_3000sccm.vtu", grid; point_data = sol[data_embed.iT,:] .-273.15)
end
  ╠═╡ =#

# ╔═╡ 2790b550-3105-4fc0-9070-d142c19678db
md"""
## Partial Pressures
"""

# ╔═╡ 31add356-6854-43c5-9ebd-ef10add6cc3d
md"""
### CO molar fraction in Cross-section
"""

# ╔═╡ bd7552d2-2c31-4834-97d9-ccdb4652242f
function SolAlongLine(data,sol)
	ng=ngas(data)		
	grid=grid_fun(data)
	
	mid_x=argmin(abs.(grid[Coordinates][1,:] .-0))
	mid_x=grid[Coordinates][1,mid_x]
	mid_y=argmin(abs.(grid[Coordinates][2,:] .-0))
	mid_y=grid[Coordinates][2,mid_y]
	

	Nodes = findall(x->x[1] == mid_x && x[2] == mid_y, eachcol(grid[Coordinates]))
	
	grid1D = grid[Coordinates][:,Nodes]
	grid1D = grid1D[3,:] # extract z-coordinate


	sol_p = []
	for i=1:(ng+2)
		push!(sol_p, sol[i,Nodes])
	end

	sol_p,grid1D,mid_x,mid_y
end

# ╔═╡ bea97fb3-9854-411c-8363-15cbef13d033
#=╠═╡
let
	data=data_embed
	sol_, grid_,midx,midy = SolAlongLine(data,sol)
	ng=ngas(data)
	(;gn, ip)=data

	cols = distinguishable_colors(ng+1, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
	pcols = map(col -> (red(col), green(col), blue(col)), cols)

	#p1=Plots.plot(title="Partial Pressures", size=(450,450), xguide="Height / mm", yguide="Pressure / bar",legend=:outertopright)
	p1=Plots.plot(title="Molar Fractions", size=(450,450), xguide="Height / mm", yguide="Molar Fraction / -",legend=:outertopright)
	for i in 1:ng
		Plots.plot!(p1, grid_./ufac"mm", vec(sol_[i])./vec(sol_[ip]), label="$(gn[i])", lw=2, ls=:auto, color=cols[i])
	end
	#Plots.plot!(p1, grid_./ufac"mm", vec(sol_[ip])./ufac"bar", color=cols[ip], label="total p", lw=2)
	#lens!(10*[0.45, .525], [0.092, 0.1], inset = (1, bbox(0.1, 0.15, 0.25, 0.25)))
	#lens!(10*[0.45, .525], [0.49, 0.5], inset = (1, bbox(0.5, 0.15, 0.25, 0.25)))

	#lens!(10*[0.45, .525], [0.068, 0.078], inset = (1, bbox(0.15, 0.56, 0.20, 0.15)))
	#lens!(10*[0.45, .525], [0.36, 0.38], inset = (1, bbox(0.55, 0.53, 0.25, 0.15)))

	p=Plots.plot(p1, size=(500,300))
	#Plots.plot(p1,p2, layout=(1,2), size=(700,450))
	#Plots.savefig(p, "../img/out/xi_1480sccm_40suns.svg")
	
end
  ╠═╡ =#

# ╔═╡ e81e803a-d831-4d62-939c-1e4a4fdec74f
md"""
# Post-Processing
"""

# ╔═╡ c4521a0c-c5af-43cd-97bc-a4a7a42d27b1
md"""
The following reactor performance metrics are considered in post processing:
1. Mass and molar flow rates
1. Catalyst layer average temperature
1. CO Yield (Converion x Selectivity)
1. Solar-to-chemical efficiency
"""

# ╔═╡ b34c1d1b-a5b8-4de9-bea9-3f2d0503c1c0
md"""
## Mass and Molar flows
"""

# ╔═╡ 84373c2e-f673-42fd-9c39-9d1f8a87c262
#=╠═╡
data_embed.mdotin/ufac"kg/hr"
  ╠═╡ =#

# ╔═╡ f39dd714-972c-4d29-bfa8-d2c3795d2eef
function massflow(data, bflux)
	ng=ngas(data)
	mdot=0.0
	for i=1:ng
		mdot += bflux[i] * data.Fluids[i].MW
	end
	mdot/ufac"kg/hr"
end

# ╔═╡ 6ae6d894-4923-4408-9b77-1067ba9e2aff
function MoleFlows(sol,sys,data)
	ng=ngas(data)
	# top: inflow	
	Itop=integrate(sys,top,sol; boundary=true)[:,Γ_top_cat] 
	# bottom: outflow
	Ibot=integrate(sys,bottom,sol; boundary=true)[:,Γ_bottom]
	
	#Itop=sum(Itop, dims=2)
	Ibot[1:ng],Itop[1:ng]
end

# ╔═╡ 7ab38bc2-9ca4-4206-a4c3-5fed673557f1
#=╠═╡
begin
	ndot_bot, ndot_top = MoleFlows(sol,sys,data_embed)
end;
  ╠═╡ =#

# ╔═╡ 8b1a0902-2542-40ed-9f91-447bffa4290f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Mass flows:

- through __bottom__ boundary __$(round(massflow(data_embed, ndot_bot),sigdigits=3))__ kg/h
- through __top__ boundary __$(round(massflow(data_embed, ndot_top),sigdigits=3))__ kg/h
"""
  ╠═╡ =#

# ╔═╡ b06a7955-6c91-444f-9bf3-72cfb4a011ec
#=╠═╡
md"""

Chemical species flows through porous frit from bottom and top:

|    | Top   | Bottom    |    |
|----|-------|-------|-------|
| $(data_embed.gn[1]) | $(abs(round(ndot_top[1]/ufac"mol/hr",digits=2)))  |   $(round(ndot_bot[1]/ufac"mol/hr",digits=2))   |  mol/hr  |
| $(data_embed.gn[2]) | $(abs(round(ndot_top[2]/ufac"mol/hr",digits=2)))  |   $(round(ndot_bot[2]/ufac"mol/hr",digits=2))   |  mol/hr  |
| $(data_embed.gn[3]) | $(abs(round(ndot_top[3]/ufac"mol/hr",digits=2)))  |   $(round(ndot_bot[3]/ufac"mol/hr",digits=2))   |  mol/hr  |
| $(data_embed.gn[4]) | $(abs(round(ndot_top[4]/ufac"mol/hr",digits=2)))  |   $(round(ndot_bot[4]/ufac"mol/hr",digits=2))   |  mol/hr  |
| $(data_embed.gn[5]) | $(abs(round(ndot_top[5]/ufac"mol/hr",digits=2)))  |   $(round(ndot_bot[5]/ufac"mol/hr",digits=2))   |  mol/hr  |
| $(data_embed.gn[6]) | $(abs(round(ndot_top[6]/ufac"mol/hr",digits=2)))  |   $(round(ndot_bot[6]/ufac"mol/hr",digits=2))   |  mol/hr  |


"""
  ╠═╡ =#

# ╔═╡ d784528f-1c3f-4522-9a6e-20ce648e2ef7
#=╠═╡
function ndo()
	(;gni) = data_embed
	dry_prod = ndot_bot[gni[:CO]]+ndot_bot[gni[:H2]]+ndot_bot[gni[:CH4]]+ndot_bot[gni[:CO2]]
end
  ╠═╡ =#

# ╔═╡ 9f377ad9-3e79-42d2-b340-9020ffce88ca
#=╠═╡
md"""
|    | Dry Molar Fraction    |    |
|----|-------|-------|
| $(data_embed.gn[1]) |   $(round(ndot_bot[1]/ndo(),digits=3)*100)   | % |
| $(data_embed.gn[2]) |   $(round(ndot_bot[2]/ndo(),digits=3)*100)   | % |
| $(data_embed.gn[3]) |   $(round(ndot_bot[3]/ndo(),digits=3)*100)   | % |
| $(data_embed.gn[5]) |   $(round(ndot_bot[5]/ndo(),digits=3)*100)   | % |
"""
  ╠═╡ =#

# ╔═╡ ed0c6737-1237-40d0-81df-32d6a842cbb0
function cat_mass(data)
	data.mcat
end

# ╔═╡ 7f2276d5-bb78-487c-868a-0486595d0113
function avg_RR(sol,sys,data)
	(;mcat)=data
	# sign convention: source term <0
	PRs=-integrate(sys,reaction,sol;)[:,2] # mol/s
	# catalyst mass specific prodcution rate of CO

	PRs[1] / mcat # mol/ (s kgcat)
end

# ╔═╡ 84093807-8ea9-4290-8bf7-bd94c5c28691
#=╠═╡
md"""
### Total reaction turnover
We obtain the total CO production rate (mol/hr) in the flow reactor from integration of the reaction rate over the domain (volume) of the catalyst layer. When dividing by the total catalyst mass in the reactor, we obtain the average catalyst mass specific reaction rate for rWGS.

__Total catalyst mass__ in the catalyst layer:
- __$(round(cat_mass(data_embed)/ufac"g",sigdigits=2))__ g

Average catalyst mass specific __reaction rate__:
- __$(round(avg_RR(sol,sys,data_embed)/ufac"mol/hr/g",sigdigits=2))__ mol/(__hr__ $\cdot$ __g__cat)
- __$(round(avg_RR(sol,sys,data_embed),sigdigits=2))__ mol/(__s $\cdot$ kg__cat)

The obtained values for average catalyst mass specific reaction rates for the rWGS reaction look reasonable when comparing with other photo-thermal catalysts operated at high light concentrations:
__Lou, D., et al. (2021).__ "A core-shell catalyst design boosts the performance of photothermal reverse water gas shift catalysis." Science China Materials 64(9): 2212-2220.
	

"""
  ╠═╡ =#

# ╔═╡ a6e61592-7958-4094-8614-e77446eb2223
md"""
##  Temperature in Domain
"""

# ╔═╡ abbe9d3a-a72b-4d8b-9893-a5e17cc817da
md"""
### Window inner surface
"""

# ╔═╡ 6737a651-0bd3-44f8-b735-21b8d6c9ed90
#=╠═╡
let
	(;iTp,iTw)=data_embed
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[2]
	end

    #bgrid = subgrid(grid, [Γ_top_cat]; boundary = true, transform = _3to2)
	#bsol=view(sol[iTw, :], bgrid)
	bgrid = subgrid(grid, [Γ_bottom]; boundary = true, transform = _3to2)
	bsol=view(sol[iTp, :], bgrid)

	scalarplot(bgrid, bsol.-273.15,colormap=:inferno,Plotter=PlutoVista,show=true)
end
  ╠═╡ =#

# ╔═╡ 94292273-ac9a-4894-909f-d078eb61d81e
md"""
### Top plane
"""

# ╔═╡ 9853479e-cfa9-4078-950e-2539c1e05961
md"""
### Bottom plane
"""

# ╔═╡ c55a6e07-c848-45ea-8412-5aff5495accb
md"""
### Cut plane
"""

# ╔═╡ 808aed68-7077-4079-be75-1bea962c716d
function TmaxSide(sol,sys,grid,data)

    (;iT)=data
	if grid_fun == prism_sq_full 
		sub=subgrid(grid,[Γ_side_front,Γ_side_back,Γ_side_right,Γ_side_left],boundary=true)
	else
		sub=subgrid(grid,[Γ_side_back,Γ_side_right],boundary=true)
	end

    Tsides = view(sol[iT, :], sub) .- 273.15
    maximum(Tsides)
		
	
end

# ╔═╡ b22387c4-cda5-424d-99b7-d7deda24c678
#=╠═╡
md"""
__Frit side maximum__ temperature: $(round(TmaxSide(sol,sys,grid,data_embed))) °C
"""
  ╠═╡ =#

# ╔═╡ 4d9c50ef-3756-492c-90a1-eeb5e51b5515
function PowerIn(sol,sys,grid,data)
	(;iT)=data
	function FluxIn(f,u,bnode,data)
	    if bnode.region==Γ_top_cat
	        (;FluxIntp,FluxEmbed,uc_window,uc_cat)=data
			
			@views x,y,z = bnode.coord[:,bnode.index]		
	        #f[iT] = FluxEmbed*FluxIntp(x,y)
			tau1_vis=uc_window.tau_vis
	        rho1_vis=uc_window.rho_vis
	        rho2_vis=uc_cat.rho_vis
	
	        f[iT] = (1-rho1_vis*rho2_vis)/tau1_vis*FluxEmbed*FluxIntp(x,y)
	    end
	end
	Pin=integrate(sys,FluxIn,sol; boundary=true)[[iT],:]
	#Pin=sum(Pin)
end

# ╔═╡ cfa366bc-f8b9-4219-b210-51b9fc5ff3f6
#=╠═╡
md"""
Total irradiated power on __catalyst surface__ (from measurement, scaled to 70 kW/m2 nominal): __$(Integer(round(sum(M12)*(0.031497*ufac"cm")^2,digits=0))) W__

Total irradiated power on __aperture__ (for interpolation on computational grid, interpolated to nom. avg. irradiation 70 kW/m2): __$(Pwr=PowerIn(sol,sys,grid,data_embed);round(sum(Pwr[Γ_top_cat]))) W__
"""
  ╠═╡ =#

# ╔═╡ 10ded7e6-3b74-43ca-b692-52934c6a95b3
md"""
## Heat flux side walls
Calculate temperature at inner reactor wall from fluxes through the side.
"""

# ╔═╡ 68fec8ae-a3ae-4d9e-9a2c-17618d999115
function TempInsideWall(sol,data)
	# get Temperature at side wall boundary + subgrid

	grid=grid_fun(data)
	(;X0,delta_gap,iT,iTp,Tamb,h,cath,shellh,k_nat_conv)=data

	# keep x-z coordinates of parent grid
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[3]
	end
	grid_xz  = subgrid(grid, [Γ_side_front], boundary=true, transform=_3to2) 
	coords=grid_xz[Coordinates]
	
	sol_xz = []
	for i=1:iTp
		sol_i = view(sol[i, :], grid_xz)
		push!(sol_xz, collect(sol_i))
	end
	sol_xz

	xs=unique(grid_xz[Coordinates][1,:])
	zs=unique(grid_xz[Coordinates][2,:])

	x_=repeat(xs,outer=length(zs))
	z_=repeat(zs,inner=length(xs))
	xz_=hcat(x_,z_)
	#xz_'[:,1] == coords[:,1]
	
	pos=[]
	for (i,xz) in enumerate(eachcol(xz_'))
		push!(pos,findall(all(coords .== xz, dims=1))[1][2])
	end
	pos
	Ts=sol_xz[iT,:][1]
	
	Ts_=reshape(Ts[pos], length(xs),length(zs))

	# calculate heat fluxes through side wall boundary
	#X=MVector{ngas(data),eltype(u)}(undef)
	#@inline mole_frac!(bnode,data,X,u)
	#@inline _,λf=dynvisc_thermcond_mix(data, u[iT], X)

	Twall_ins = []
	for T in Ts
		_,λf = dynvisc_thermcond_mix(data, T, X0)
		# include shell height
		# qside = (T-Tamb)/(delta_gap/λf+(h+cath)/(shellh*k_nat_conv))
		# w/o shell height
		qside = (T-Tamb)/(delta_gap/λf+1/k_nat_conv)
		push!(Twall_ins, T-qside*delta_gap/λf)
	end
	vec(Twall_ins)
	Twall_ins_=reshape(Twall_ins[pos], length(xs),length(zs)) .-273.15
	#Ts_ .-273.15
end

# ╔═╡ 58d5952c-375f-4744-8c3e-093afec65f5d
#=╠═╡
let
	TempInsideWall(sol,data_embed)
end
  ╠═╡ =#

# ╔═╡ a55c5ee7-2274-4447-b0b2-58052f064bc9
function areas(sol,sys,grid,data)
	iT = data.iT
	function area(f,u,bnode,data)
		# repurpose temperature index to hold area information
		f[iT] = one(eltype(u))
	end

	integrate(sys,area,sol; boundary=true)[iT,:]
end

# ╔═╡ 15604034-91fd-4fd4-b09e-e3c5cfe7a265
function T_avg(sol,sys,grid,data)

	iT = data.iT	
	function T_avg_(f,u,bnode,data)			
		f[iT] = u[iT]		
	end

	areas_=areas(sol,sys,grid,data)
	
	T_int=integrate(sys,T_avg_,sol; boundary=true)[iT,:]
	T_int./areas_	
end

# ╔═╡ fec9ca6d-d815-4b50-bec7-f8fb3d8195ba
function TopPlane(data,sol)
	ng=ngas(data)
	grid=grid_fun(data)
	(;iTp,h,cath,wi)=data
	wi=wi/2
	
	bid = maximum(grid[BFaceRegions])+1
	if grid_fun == prism_sq_full
		#bfacemask!(grid, [-wi,-wi,h+cath],[wi,wi,h+cath],bid)
		bfacemask!(grid, [-wi,-wi,h],[wi,wi,h],bid)
	else
		bfacemask!(grid, [0,0,h+cath],[wi,wi,h+cath],bid)
	end

	# keep x-y coordinates of parent grid
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[2]
	end
	#grid_1D  = subgrid(grid, [bid], boundary=true, transform=_3to1) 
	grid_2D  = subgrid(grid, [bid], boundary=true, transform=_3to2) 

	sol_p = []
	for i=1:iTp
		sol_i = view(sol[i, :], grid_2D)
		push!(sol_p, collect(sol_i))
	end

	sol_p, grid_2D
end

# ╔═╡ 48b99de4-682d-439b-935e-408510c44e37
#=╠═╡
let
	data=ModelData()
	(;wi,FluxIntp)=data
	#x=range(-wi/2,wi/2, length=21)
	#y=range(-wi/2,wi/2, length=21)
	_, grid_xy = TopPlane(data_embed,sol)
	#(;wi)=data_embed
	x=@views unique(grid_xy[Coordinates][1,:])
	y=@views unique(grid_xy[Coordinates][2,:])
	Plots.plot(xguide="X Coordinate / cm", yguide="Y Coordinate / cm",title="Interpolated Flux Profile", colorbar_title="Irradiation Flux / kW m-2", xlim=(-wi/2/ufac"cm", wi/2/ufac"cm"), ylim=(-wi/2/ufac"cm", wi/2/ufac"cm"),size=(400,400))

	Plots.contour!(x/ufac"cm",y/ufac"cm",FluxIntp(x,y) / ufac"kW/m^2", lw=0, aspect_ratio = 1, fill = true)
	#Plots.plot!([-5,5],[-5,-5],label=:none,c=:black,lw=2)
	#Plots.plot!([-5,5],[5,5],label=:none,c=:black,lw=2)
	#Plots.plot!([-5,-5],[-5,5],label=:none,c=:black,lw=2)
	#Plots.plot!([5,5],[-5,5],label=:none,c=:black,lw=2)
	#Plots.savefig("../img/out/Flux_interpol.svg")
end
  ╠═╡ =#

# ╔═╡ 06ad2f54-6ae9-4c32-8391-e6cc8c5cfe68
#=╠═╡
let
	(;iT,wi)=data_embed
	sol_xy, grid_xy = TopPlane(data_embed,sol)
	coords=grid_xy[Coordinates]

	xs=unique(grid_xy[Coordinates][1,:])
	ys=unique(grid_xy[Coordinates][2,:])

	x_=repeat(xs,outer=length(ys))
	y_=repeat(ys,inner=length(xs))
	xy_=hcat(x_,y_)
	xy_'[:,1] == coords[:,1]
	
	pos=[]
	for (i,xy) in enumerate(eachcol(xy_'))
		push!(pos,findall(all(coords .== xy, dims=1))[1][2])
	end
	pos
	Ts=sol_xy[iT,:][1] .- 273.15

	p=Plots.plot(xguide="X Coordinate / cm", yguide="Y Coordinate / cm", xlim=(-wi/2/ufac"cm", wi/2/ufac"cm"), ylim=(-wi/2/ufac"cm", wi/2/ufac"cm"), xticks=(-wi/2/ufac"cm"):2.0:(wi/2/ufac"cm"), yticks=(-wi/2/ufac"cm"):2.0:(wi/2/ufac"cm"), size=(300,300))
	Plots.contourf!(p,xs/ufac"cm",ys/ufac"cm", lw=0,Ts[pos],aspect_ratio=1, clim=(0,700), levels=6, colormap=:inferno)
	#Plots.savefig(p,"../img/out/Ttop_70suns.svg")
	
	#@views x,y,z = bnode.coord[:,bnode.index]		
	#Glamp =FluxEmbed*FluxIntp(x,y)
end
  ╠═╡ =#

# ╔═╡ 77e8f837-bacb-428a-ad66-bb6dd878a438
#=╠═╡
function T_uc_lc(sol,data)
	ng=ngas(data)
	grid_=grid_fun(data)
	(;iT,h,cath)=data

	wi_uc=5*ufac"cm"
	#wi_lc=5*ufac"cm"
	wi_lc=7*ufac"cm" # incorporate increased mixing of gas in cross flow direction
	
	function Tcat_(f,u,bnode,data)
		@views x,y,z = bnode.coord[:,bnode.index]	
		if bnode.region == Γ_top_cat
			if x >= -wi_uc && x <= wi_uc && y >= -wi_uc && y <= wi_uc
				f[iT] = u[iT]			
			end
		elseif bnode.region == Γ_bottom
			if x >= -wi_lc && x <= wi_lc && y >= -wi_lc && y <= wi_lc
				f[iT] = u[iT]			
			end
		end
	end
		
	function area_(f,u,bnode,data)
		@views x,y,z = bnode.coord[:,bnode.index]	
		if bnode.region == Γ_top_cat
			if x >= -wi_uc && x <= wi_uc && y >= -wi_uc && y <= wi_uc
				f[iT] = 1		
			end
		elseif bnode.region == Γ_bottom
			if x >= -wi_lc && x <= wi_lc && y >= -wi_lc && y <= wi_lc
				f[iT] = 1
			end
		end
	end
	
	Auc,Alc = integrate(sys,area_,sol; boundary=true)[iT,[Γ_top_cat,Γ_bottom]]
	Tuc,Tlc = integrate(sys,Tcat_,sol; boundary=true)[iT,[Γ_top_cat,Γ_bottom]] ./ (Auc,Alc) .- 273.15

end
  ╠═╡ =#

# ╔═╡ 2739bcff-3fb0-4169-8a1a-2b0a14998cec
#=╠═╡
md"""
For comparison with experiments, compute the average temperatures in the upper and lower chambers.

__Upper Surface__ average temperature: 
$(Integer(round(T_uc_lc(sol,data_embed)[1]))) °C

__Lower Surface__ average temperature: 
$(Integer(round(T_uc_lc(sol,data_embed)[2]))) °C
"""
  ╠═╡ =#

# ╔═╡ ffbd1d36-3bd3-4965-9690-db299d9347f0
#=╠═╡
function Flux_avg(sol,data,wi)
	ng=ngas(data)
	grid_=grid_fun(data)
	(;iT,FluxIntp,FluxEmbed,uc_window,uc_cat)=data

	#wi_uc=5*ufac"cm"
	#wi_lc=7*ufac"cm" # incorporate increased mixing of gas in cross flow direction
	
	function Favg_(f,u,bnode,data)
		@views x,y,z = bnode.coord[:,bnode.index]	
		if bnode.region == Γ_top_cat
			if x >= -wi && x <= wi && y >= -wi && y <= wi			
				tau1_vis=uc_window.tau_vis
		        rho1_vis=uc_window.rho_vis
		        rho2_vis=uc_cat.rho_vis
		
		        f[iT] = (1-rho1_vis*rho2_vis)/tau1_vis*FluxEmbed*FluxIntp(x,y)
			end
		end
	end
		
	function area_(f,u,bnode,data)
		@views x,y,z = bnode.coord[:,bnode.index]	
		if bnode.region == Γ_top_cat
			if x >= -wi && x <= wi && y >= -wi && y <= wi
				f[iT] = 1		
			end
		end
	end
	
	A = integrate(sys,area_,sol; boundary=true)[iT,Γ_top_cat]
	Favg = integrate(sys,Favg_,sol; boundary=true)[iT,Γ_top_cat] / A

end
  ╠═╡ =#

# ╔═╡ d009a262-fe30-49b1-8f43-463758c23093
function BotPlane(data,sol)
	ng=ngas(data)
	grid=grid_fun(data)
	(;iTp,h,cath,wi)=data
	wi=wi/2
	
	bid = maximum(grid[BFaceRegions])+1
	if grid_fun == prism_sq_full
		bfacemask!(grid, [-wi,-wi,0],[wi,wi,0],bid)
	else
		bfacemask!(grid, [0,0,h+cath],[wi,wi,h+cath],bid)
	end

	# keep x-y coordinates of parent grid
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[2]
	end
	#grid_1D  = subgrid(grid, [bid], boundary=true, transform=_3to1) 
	grid_2D  = subgrid(grid, [bid], boundary=true, transform=_3to2) 

	sol_p = []
	for i=1:iTp
		sol_i = view(sol[i, :], grid_2D)
		push!(sol_p, collect(sol_i))
	end

	sol_p, grid_2D
end

# ╔═╡ d2bc148c-badf-4378-ad61-cbeb52c9ab3f
#=╠═╡
let
	(;iT,wi)=data_embed
	sol_xy, grid_xy = BotPlane(data_embed,sol)
	coords=grid_xy[Coordinates]

	xs=unique(grid_xy[Coordinates][1,:])
	ys=unique(grid_xy[Coordinates][2,:])

	x_=repeat(xs,outer=length(ys))
	y_=repeat(ys,inner=length(xs))
	xy_=hcat(x_,y_)
	xy_'[:,1] == coords[:,1]
	
	pos=[]
	for (i,xy) in enumerate(eachcol(xy_'))
		push!(pos,findall(all(coords .== xy, dims=1))[1][2])
	end
	pos
	Ts=sol_xy[iT,:][1] .- 273.15

	p=Plots.plot(xguide="X Coordinate / cm", yguide="Y Coordinate / cm", xlim=(-wi/2/ufac"cm", wi/2/ufac"cm"), ylim=(-wi/2/ufac"cm", wi/2/ufac"cm"), size=(300,300))
	Plots.contourf!(p,xs/ufac"cm",ys/ufac"cm", lw=0,Ts[pos],aspect_ratio=1)
	#Plots.savefig(p,"../img/out/Ttop.svg")
	
	#@views x,y,z = bnode.coord[:,bnode.index]		
	#Glamp =FluxEmbed*FluxIntp(x,y)
end
  ╠═╡ =#

# ╔═╡ 746e1a39-3c9b-478f-b371-3cb8333e93b1
function CutPlane(data,sol)
	ng=ngas(data)
	grid=grid_fun(data)
	(;iTp,h,cath,wi)=data
	wi=data.wi/2
	
	bid = maximum(grid[BFaceRegions])+1
	if grid_fun == prism_sq_full
		bfacemask!(grid, [-wi,0.0,0.0],[wi,0.0,h+cath],bid)
	else
		bfacemask!(grid, [0,0.024,0],[wi,0.024,h+cath],bid)
	end

	# keep x-z coordinates of parent grid
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[3]
	end
	#grid_1D  = subgrid(grid, [bid], boundary=true, transform=_3to1) 
	grid_2D  = subgrid(grid, [bid], boundary=true, transform=_3to2) 

	sol_p = []
	for i=1:iTp
		sol_i = view(sol[i, :], grid_2D)
		push!(sol_p, collect(sol_i))
	end

	sol_p, grid_2D
end

# ╔═╡ 60b475d1-219e-4037-bd36-e0671cfa1893
# ╠═╡ skip_as_script = true
#=╠═╡
let
	data=data_embed
	(;gni,ip)=data
	sol_xz, grid_xz = CutPlane(data,sol)
	#idx = gni[:CO]
	idx = gni[:CH4]
	#idx = gni[:CO2]
	#idx = gni[:H2]
	scalarplot(grid_xz, sol_xz[idx] ./sol_xz[ip], resolution=(650,600), zoom=1.2, aspect=5.0, Plotter=PlutoVista)
end
  ╠═╡ =#

# ╔═╡ 7198a578-1c74-41b7-aec7-568637d7891b
#=╠═╡
let
	(;gni,ip,wi,h,cath)=data_embed
	sol_xz, grid_xz = CutPlane(data_embed,sol)
	coords=grid_xz[Coordinates]

	xs=unique(grid_xz[Coordinates][1,:])
	zs=unique(grid_xz[Coordinates][2,:])

	x_=repeat(xs,outer=length(zs))
	z_=repeat(zs,inner=length(xs))
	xz_=hcat(x_,z_)
	xz_'[:,1] == coords[:,1]
	
	pos=[]
	for (i,xz) in enumerate(eachcol(xz_'))
		push!(pos,findall(all(coords .== xz, dims=1))[1][2])
	end
	pos
	pCOs=sol_xz[gni[:CO],:][1]
	pCOs[pCOs.<0] .= 0
	pts=sol_xz[ip,:][1]
	xCOs = pCOs./pts

	p2=Plots.plot(xguide="X Coordinate / cm", yguide="Z Coordinate / cm", xlim=(-wi/2/ufac"cm", wi/2/ufac"cm"), ylim=(0, h/ufac"cm"), xticks=(-wi/2/ufac"cm"):2.0:(wi/2/ufac"cm"), size=(300,300))
	Plots.contourf!(p2,xs/ufac"cm",zs/ufac"cm", aspect_ratio=10, lw=0, xCOs[pos], colormap=:viridis, clim=(0,0.2), levels=4)
	#Plots.savefig(p2,"../img/out/xCOcut_70suns.svg")
end
  ╠═╡ =#

# ╔═╡ 3dfc1ac2-93dc-4b5b-adbc-e66e92bf76b0
#=╠═╡
let
	(;iT,wi,h,cath)=data_embed
	sol_xz, grid_xz = CutPlane(data_embed,sol)
	coords=grid_xz[Coordinates]

	xs=unique(grid_xz[Coordinates][1,:])
	zs=unique(grid_xz[Coordinates][2,:])

	x_=repeat(xs,outer=length(zs))
	z_=repeat(zs,inner=length(xs))
	xz_=hcat(x_,z_)
	xz_'[:,1] == coords[:,1]
	
	pos=[]
	for (i,xz) in enumerate(eachcol(xz_'))
		push!(pos,findall(all(coords .== xz, dims=1))[1][2])
	end
	pos
	Ts=sol_xz[iT,:][1] .- 273.15

	p2=Plots.plot(xguide="X Coordinate / cm", yguide="Z Coordinate / cm", xlim=(-wi/2/ufac"cm", wi/2/ufac"cm"), ylim=(0, (h+cath)/ufac"cm"), xticks=(-wi/2/ufac"cm"):2.0:(wi/2/ufac"cm"), size=(300,300))
	Plots.contourf!(p2,xs/ufac"cm",zs/ufac"cm", aspect_ratio=10, lw=0,Ts[pos],)
	#Plots.savefig(p2,"../img/out/Tcut_70suns.svg")
end
  ╠═╡ =#

# ╔═╡ 808e5a71-572f-4b0c-aeb3-9513d969dada
function CatDims(grid)
	Cells = findall(x->x == 2, grid[CellRegions]) # catalyst layer
	Nodes = grid[CellNodes][:,Cells]
	Nodes = unique(reshape(Nodes,1,:))
	coords = grid[Coordinates][:,Nodes]

	mincatwi, maxcatwi = minimum(coords[1,:]), maximum(coords[1,:])
	mincatle, maxcatle = minimum(coords[2,:]), maximum(coords[2,:])
	mincathe, maxcathe = minimum(coords[3,:]), maximum(coords[3,:])

	catA = (maxcatwi-mincatwi)*(maxcatle-mincatle)
	catV = (maxcatwi-mincatwi)*(maxcatle-mincatle)*(maxcathe-mincathe)
	catA, catV
end

# ╔═╡ 68e2628a-056a-4ec3-827f-2654f49917d9
function Tcatavg(sol,sys,grid,data)

	function Tcat_(f,u,bnode,data)
		iT = data.iT		
		f[iT] = u[iT]		
	end

	catA, catV = CatDims(grid)
		
	Tcat_avg_surf=integrate(sys,Tcat_,sol; boundary=true)[data.iT,Γ_top_cat] / catA - 273.15
	Tcat_avg_vol=integrate(sys,Tcat_,sol; )[data.iT,2] / catV - 273.15
	Tcat_avg_surf, Tcat_avg_vol
end

# ╔═╡ 9952c815-5459-44ff-b1f8-07ab24ce0c53
#=╠═╡
md"""
__Volumetric average__ Catalyst layer temperature: $(round(Tcatavg(sol,sys,grid,data_embed)[2])) °C
"""
  ╠═╡ =#

# ╔═╡ 4cde8752-bbdf-4b83-869e-46b78bb4adb5
md"""
## Yield of CO
"""

# ╔═╡ 7bf4d925-56e1-4c85-a40a-fef1917f501e
md"""
```math
Y_{\text{CO}} = \frac{\dot n_{\text{CO}}}{\dot n^0_{\text{CO}_2}}
```
"""

# ╔═╡ b3bf7b7d-eb38-4a32-87f7-5aef098ad03e
#=╠═╡
function Yield_CO(sol,sys,data)
	ndot_top,ndot_top = MoleFlows(sol,sys,data)
	
	ndot_bot[data.gni[:CO]] / abs(ndot_top[data.gni[:CO2]])
end
  ╠═╡ =#

# ╔═╡ 0f7cce89-3add-4316-bcc2-a924064af884
#=╠═╡
md"""
Yield of CO relative to inflow of ``\text{CO}_2``: $(round(Yield_CO(sol,sys,data_embed),sigdigits=2))
"""
  ╠═╡ =#

# ╔═╡ aebcb161-425c-4a46-aaed-d4b13b3e6654
md"""
## Solar to Chemical Efficiency (STC)
"""

# ╔═╡ 6bbd0496-c275-4cd5-bf48-9ea4e77a9091
md"""
Defined with the reaction enthalpy of the RWGS reaction at standart conditions.
"""

# ╔═╡ 1f3b2a0e-c76f-462c-a85c-e85815fd7e5a
md"""
```math
	\text{STC} = \frac{\dot n_{\text{CO}} \Delta H^0_{\text{RWGS}}}{I_{\text{lamp}} A}
```
"""

# ╔═╡ 47161886-9a5c-41ac-abf5-bbea82096d5a
#=╠═╡
function STCefficiency(sol,sys,data,grid)
	ndot_bot,ndot_top = MoleFlows(sol,sys,data)
	# R1 = RWGS in S3P kinetic model
	# R2 = RWGS in Xu & Froment kinetic model
	(;gni,kinpar,Glamp,Ac) = data
	DHi = 0.0
	if kinpar == XuFroment
		DHi = kinpar.ΔHi[:R2]
	elseif kinpar == S3P
		DHi = kinpar.ΔHi[:R1]
	end
	if grid_fun == prism_sq_full
		Pwr=PowerIn(sol,sys,grid,data_embed)
		STC_Atot = ndot_bot[gni[:CO]] * -DHi / sum(Pwr[[Γ_top_cat,Γ_top_frit]])
		STC_Acat = ndot_bot[gni[:CO]] * -DHi / sum(Pwr[Γ_top_cat])
	else
		STC_Atot = ndot_bot[gni[:CO]] * -DHi / (Glamp*Ac/4)
		catA, _ = CatDims(grid)
		STC_Acat = ndot_bot[gni[:CO]] * -DHi / (Glamp*catA)
	end
	
	#STC_Acat = ndot_bot[data.gni[:CO]] * -data.kinpar.ΔHi[:R2] / (data.Glamp*catA)
	
	STC_Atot, STC_Acat
end
  ╠═╡ =#

# ╔═╡ 428afb22-55ab-4f64-805f-7e15ad4cf23f
#=╠═╡
md"""
Solar-to-chemical efficiency as defined above: $(round(STCefficiency(sol,sys,data_embed,grid)[1]*100, sigdigits=2)) %
"""
  ╠═╡ =#

# ╔═╡ 73d1dab0-69dc-40a0-9e38-2326fe9938c2
md"""
## Energy Flow Analysis
"""

# ╔═╡ 89960a27-53a8-40a4-a2af-e338faa0f551
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The different energy flow mechanisms are listed below for analysis.
An energy balance is made around the reactor volume as indicated by the dashed red line.
Energy that enters through the aperture and leaves the balance volume via the indicated mechanisms numbered 1-8.

$(LocalResource("../img/EnergyFlows.png"))
"""

  ╠═╡ =#

# ╔═╡ 604c73ef-3581-4d31-b0c6-564889eb0ed2
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
To derive the relation for irradiation leaving through the quartz window, look at the quartz window in detail:

$(LocalResource("../img/QTopRad.png",:width => 400))

```math
G^{\text{vis}}_{1,\text{top}} = \rho_1^{\text{vis}} G_{\text{HFSS}} + \tau_1^{\text{vis}} (\phi_{12} G_2^{\text{vis}} +\phi_{13} G_3^{\text{vis}})
```
```math
G^{\text{IR}}_{1,\text{top}} = \tau_1^{\text{IR}} (\phi_{12} G_2^{\text{IR}} +\phi_{13} G_3^{\text{IR}}) + \epsilon_1 \sigma T_1^4
```
"""
  ╠═╡ =#

# ╔═╡ e13a0afb-5717-4593-b9d7-e31e34a2e9c0
md"""
Where the surface brightnesses $G_2^{\text{vis/IR}}$ and $G_3^{\text{vis/IR}}$ for the photo-catalyst and uncoated frit surfaces in the visible and IR spectrum respectively are defined as follows:

```math
\begin{align}
G^{\text{IR}}_{2} &= \epsilon_2 \sigma T_2^4 + \rho_2^{\text{IR}} G^{\text{IR}}_{1,\text{bot}}\\
G^{\text{vis}}_{2} &= \rho_2^{\text{vis}} G^{\text{vis}}_{1,\text{bot}}
\end{align}
```

"""

# ╔═╡ 7eb683e4-a183-47b7-8090-3171f420fc7d
md"""
```math
\begin{align}
G^{\text{IR}}_{3} &= \epsilon_3 \sigma T_3^4 + \rho_3^{\text{IR}} G^{\text{IR}}_{1,\text{bot}}\\
G^{\text{vis}}_{3} &= \rho_3^{\text{vis}} G^{\text{vis}}_{1,\text{bot}}
\end{align}
```
"""

# ╔═╡ 50c6434f-03e1-41bf-b1a0-3072af659ea1
md"""
The surface brigthnesses of the catalyst facing side of the quartz window denoted by $G^{\text{IR}}_{1,\text{bot}}$ and $G^{\text{vis}}_{1,\text{bot}}$ are defined as above for irradiation exchange in the top chamber:
"""

# ╔═╡ b842e5b7-fb6a-4225-838a-cf857cf7c28b
md"""
```math
G^{\text{IR}}_{1,\text{bot}} = \frac{\epsilon_1\sigma T_1^4 + \rho_1^{\text{IR}} \left( \phi_{12}\epsilon_2 \sigma T_2^4+\phi_{13}\epsilon_3 \sigma T_3^4\right)}{1-\rho_1^{\text{IR}} \left( \phi_{12} \rho_2^{\text{IR}} + \phi_{13} \rho_3^{\text{IR}} \right)}
```
"""

# ╔═╡ 71a5d5b4-29ac-4dab-b865-e8561365d8ff
md"""
```math
G^{\text{vis}}_{1,\text{bot}} = \frac{\tau_1^{\text{vis}} G_{\text{lamp}}}{1-\rho_1^{\text{vis}} \left( \phi_{12} \rho_2^{\text{vis}} + \phi_{13} \rho_3^{\text{vis}} \right)} \\
```
"""

# ╔═╡ b2eeeb44-d8fe-452b-8905-642dfa043db8
md"""
Due to missing information on the flow field that develops in the chambers, only conductive heat transfer through a stagnant fluid is considered. This is the lower bound for convective heat transfer. The gas thermal conductivity is evaluated at the mean temperature between the two surfaces.

```math
	\begin{align}
	\dot q_{\text{cond}} &= - \lambda \nabla T \\
						 &\approx - \lambda(T_{\text m}) \frac{\Delta T}{\Delta z}
	\end{align}
```
"""

# ╔═╡ Cell order:
# ╠═11ac9b20-6a3c-11ed-0bb6-735d6fbff2d9
# ╠═863c9da7-ef45-49ad-80d0-3594eca4a189
# ╟─94ed0f8b-13c1-4460-a391-5057cff401e0
# ╟─2ed3223e-a604-410e-93d4-016580f49093
# ╟─390c7839-618d-4ade-b9be-ee9ed09a77aa
# ╠═37adb8da-3ad5-4b41-8f08-85da19e15a53
# ╠═e21b9d37-941c-4f2c-9bdf-956964428f90
# ╠═0a911687-aff4-4c77-8def-084293329f35
# ╠═ff58b0b8-2519-430e-8343-af9a5adcb135
# ╟─985718e8-7ed7-4c5a-aa13-29462e52d709
# ╟─2554b2fc-bf5c-4b8f-b5e9-8bc261fe597b
# ╟─f4dcde90-6d8f-4b17-b4ec-367d2372637f
# ╟─3703afb0-93c4-4664-affe-b723758fb56b
# ╟─21d0195b-b170-460d-989e-f9d00b511237
# ╟─8f4843c6-8d2b-4e24-b6f8-4eaf3dfc9bf0
# ╟─66b55f6b-1af5-438d-aaa8-fe4745e85426
# ╟─8528e15f-cce7-44d7-ac17-432f92cc5f53
# ╠═ed7941c4-0485-4d84-ad5b-383eb5cae70a
# ╟─a6afe118-dcbd-4126-8646-c7268acfacf3
# ╠═78cf4646-c373-4688-b1ac-92ed5f922e3c
# ╟─a60ce05e-8d92-4172-b4c1-ac3221c54fe5
# ╟─24374b7a-ce77-45f0-a7a0-c47a224a0b06
# ╟─4865804f-d385-4a1a-9953-5ac66ea50057
# ╟─3bf71cea-4f73-47da-b5ed-2cae3ec3d18b
# ╟─bcaf83fb-f215-428d-9c84-f5b557fe143f
# ╟─7f94d703-2759-4fe1-a8c8-ddf26732a6ca
# ╟─906ad096-4f0c-4640-ad3e-9632261902e3
# ╟─8139166e-42f9-41c3-a360-50d3d4e5ee86
# ╟─44d91c2e-8082-4a90-89cc-81aba783d5ac
# ╟─39e74955-aab6-4bba-a1b8-b2307b45e673
# ╟─6798d5e2-b8c7-4f54-aa71-6ea1ccab78fb
# ╟─ed3609cb-8483-4184-a385-dca307d13f17
# ╟─58d0610b-1739-4260-8d16-5a31ba362d69
# ╠═d0993435-ed0b-4a43-b848-26f5266017a1
# ╟─634d1042-b110-45ef-bfbe-51b827fc922f
# ╠═29f34e55-91ab-4b6d-adb2-a58412af95f6
# ╠═f9ba467a-cefd-4d7d-829d-0889fc6d0f5e
# ╠═ecb8df4c-fb61-472b-ad0b-e55b49905442
# ╠═e459bbff-e065-49e1-b91b-119add7d4a71
# ╠═48b99de4-682d-439b-935e-408510c44e37
# ╠═cfa366bc-f8b9-4219-b210-51b9fc5ff3f6
# ╠═f4ebb596-824a-4124-afb7-c368c1cabb00
# ╟─3aef8203-ce28-4197-a43e-784840f7bc1e
# ╠═2bf61cf4-2f4a-4c48-ad96-5de82403f3a0
# ╠═d9dee38e-6036-46b8-bc06-e545baa06789
# ╠═e58ec04f-023a-4e00-98b8-f9ae85ca506f
# ╟─80d2b5db-792a-42f9-b9c2-91d0e18cfcfb
# ╟─6de46bcb-9aff-4a21-b8d6-f14e64aac95c
# ╟─b1ef2a89-27db-4f21-a2e3-fd6356c394da
# ╠═79e6cb48-53e2-4655-9d77-51f3deb67b94
# ╟─4dae93b2-be63-4ee9-bc1e-871a31ade811
# ╠═9547ed7c-3304-4c63-a6c1-5f84e0001c54
# ╠═9d191a3a-e096-4ad7-aae6-bbd63d478fa2
# ╠═2d51f54b-4cff-4253-b17f-217e3261f36d
# ╠═fdaeafb6-e29f-41f8-8468-9c2b03b9eed7
# ╠═2e631f58-6a4b-4c6d-86ad-b748fd2d463a
# ╟─a395374d-5897-4764-8499-0ebe7f2b4239
# ╟─7db215d7-420e-435b-8889-43c4fff150e0
# ╟─8b485112-6b1a-4b38-af91-deb9b79527e0
# ╟─25edaf7b-4051-4934-b4ad-a4655698a6c7
# ╟─3fe2135d-9866-4367-8faa-56cdb42af7ed
# ╟─652497ee-d07b-45e2-aeaf-87ad5bcc23ad
# ╟─6bd59a54-f059-4646-b053-0fa41ead87fd
# ╟─f5d78670-a98b-46a1-8bf3-3d2599cfdd88
# ╠═162122bc-12ae-4a81-8df6-86498041be40
# ╠═221a1ee4-f7e9-4233-b45c-a715c9edae5f
# ╟─2228bbd4-bc84-4617-a837-2bf9bba76793
# ╟─09d976a7-f6c6-465b-86f4-9bc654ae158c
# ╟─5c9aad13-914a-4f5e-af32-a9c6403c52d0
# ╟─ff764eea-e282-4fed-90b4-5f418ae426f0
# ╟─dfed01d4-8d88-4ce7-afe5-1fd27e2cc746
# ╟─4ebbe06f-0993-4c5c-9af3-76b2b645e592
# ╠═7da59e27-62b9-4b89-b315-d88a4fd34f56
# ╠═40906795-a4dd-4e4a-a62e-91b4639a48fa
# ╠═dec81e61-d1eb-4a6c-8648-a07b98a8a7a9
# ╠═edd9fdd1-a9c4-4f45-8f63-9681717d417f
# ╠═1a8410fa-caac-4023-8c96-a0042f22d88c
# ╠═29d66705-3d9f-40b1-866d-dd3392a1a268
# ╟─02b76cda-ffae-4243-ab40-8d0fe1325776
# ╠═c0de2aff-8f7f-439c-b931-8eb8fbfcd45d
# ╠═51b03d3b-3c7c-4019-8ed8-bf1aaa0b1ddb
# ╠═2191bece-e186-4d8e-8a21-3830441baf11
# ╠═b6381008-0280-404c-a86c-9c9c3c9f82eb
# ╟─44aa5b49-d595-4982-bbc8-100d2f199415
# ╠═333b5c80-259d-47aa-a441-ee7894d6c407
# ╠═aa498412-e970-45f2-8b11-249cc5c2b18d
# ╟─e25e7b7b-47b3-457c-995b-b2ee4a87710a
# ╠═3a35ac76-e1b7-458d-90b7-d59ba4f43367
# ╟─b2df1087-6628-4889-8cd6-c5ee7629cd93
# ╠═3bd80c19-0b49-43f6-9daa-0c87c2ea8093
# ╠═b99610da-ca0c-46e7-b0a7-bc7858a45a35
# ╟─2790b550-3105-4fc0-9070-d142c19678db
# ╠═bea97fb3-9854-411c-8363-15cbef13d033
# ╟─31add356-6854-43c5-9ebd-ef10add6cc3d
# ╠═60b475d1-219e-4037-bd36-e0671cfa1893
# ╠═7198a578-1c74-41b7-aec7-568637d7891b
# ╠═bd7552d2-2c31-4834-97d9-ccdb4652242f
# ╟─e81e803a-d831-4d62-939c-1e4a4fdec74f
# ╟─c4521a0c-c5af-43cd-97bc-a4a7a42d27b1
# ╟─b34c1d1b-a5b8-4de9-bea9-3f2d0503c1c0
# ╠═84373c2e-f673-42fd-9c39-9d1f8a87c262
# ╟─8b1a0902-2542-40ed-9f91-447bffa4290f
# ╠═f39dd714-972c-4d29-bfa8-d2c3795d2eef
# ╟─b06a7955-6c91-444f-9bf3-72cfb4a011ec
# ╟─9f377ad9-3e79-42d2-b340-9020ffce88ca
# ╠═d784528f-1c3f-4522-9a6e-20ce648e2ef7
# ╠═6ae6d894-4923-4408-9b77-1067ba9e2aff
# ╠═7ab38bc2-9ca4-4206-a4c3-5fed673557f1
# ╟─84093807-8ea9-4290-8bf7-bd94c5c28691
# ╠═ed0c6737-1237-40d0-81df-32d6a842cbb0
# ╠═7f2276d5-bb78-487c-868a-0486595d0113
# ╟─a6e61592-7958-4094-8614-e77446eb2223
# ╠═2739bcff-3fb0-4169-8a1a-2b0a14998cec
# ╟─abbe9d3a-a72b-4d8b-9893-a5e17cc817da
# ╠═6737a651-0bd3-44f8-b735-21b8d6c9ed90
# ╟─94292273-ac9a-4894-909f-d078eb61d81e
# ╠═06ad2f54-6ae9-4c32-8391-e6cc8c5cfe68
# ╟─9853479e-cfa9-4078-950e-2539c1e05961
# ╟─d2bc148c-badf-4378-ad61-cbeb52c9ab3f
# ╟─c55a6e07-c848-45ea-8412-5aff5495accb
# ╠═3dfc1ac2-93dc-4b5b-adbc-e66e92bf76b0
# ╟─9952c815-5459-44ff-b1f8-07ab24ce0c53
# ╟─b22387c4-cda5-424d-99b7-d7deda24c678
# ╠═15604034-91fd-4fd4-b09e-e3c5cfe7a265
# ╠═808aed68-7077-4079-be75-1bea962c716d
# ╠═4d9c50ef-3756-492c-90a1-eeb5e51b5515
# ╟─10ded7e6-3b74-43ca-b692-52934c6a95b3
# ╠═58d5952c-375f-4744-8c3e-093afec65f5d
# ╠═68fec8ae-a3ae-4d9e-9a2c-17618d999115
# ╠═a55c5ee7-2274-4447-b0b2-58052f064bc9
# ╠═68e2628a-056a-4ec3-827f-2654f49917d9
# ╠═fec9ca6d-d815-4b50-bec7-f8fb3d8195ba
# ╠═77e8f837-bacb-428a-ad66-bb6dd878a438
# ╠═ffbd1d36-3bd3-4965-9690-db299d9347f0
# ╠═d009a262-fe30-49b1-8f43-463758c23093
# ╠═746e1a39-3c9b-478f-b371-3cb8333e93b1
# ╟─808e5a71-572f-4b0c-aeb3-9513d969dada
# ╟─4cde8752-bbdf-4b83-869e-46b78bb4adb5
# ╟─0f7cce89-3add-4316-bcc2-a924064af884
# ╟─7bf4d925-56e1-4c85-a40a-fef1917f501e
# ╠═b3bf7b7d-eb38-4a32-87f7-5aef098ad03e
# ╟─aebcb161-425c-4a46-aaed-d4b13b3e6654
# ╟─6bbd0496-c275-4cd5-bf48-9ea4e77a9091
# ╟─1f3b2a0e-c76f-462c-a85c-e85815fd7e5a
# ╟─428afb22-55ab-4f64-805f-7e15ad4cf23f
# ╠═47161886-9a5c-41ac-abf5-bbea82096d5a
# ╟─73d1dab0-69dc-40a0-9e38-2326fe9938c2
# ╟─89960a27-53a8-40a4-a2af-e338faa0f551
# ╟─604c73ef-3581-4d31-b0c6-564889eb0ed2
# ╟─e13a0afb-5717-4593-b9d7-e31e34a2e9c0
# ╟─7eb683e4-a183-47b7-8090-3171f420fc7d
# ╟─50c6434f-03e1-41bf-b1a0-3072af659ea1
# ╟─b842e5b7-fb6a-4225-838a-cf857cf7c28b
# ╟─71a5d5b4-29ac-4dab-b865-e8561365d8ff
# ╟─b2eeeb44-d8fe-452b-8905-642dfa043db8