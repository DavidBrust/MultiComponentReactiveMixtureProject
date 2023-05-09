### A Pluto.jl notebook ###
# v0.19.25

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
	using VoronoiFVM
	using ExtendableGrids, GridVisualize
	using NLsolve
	using LinearAlgebra, LinearSolve,LinearSolvePardiso
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots
	using PlutoUI

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 863c9da7-ef45-49ad-80d0-3594eca4a189
PlutoUI.TableOfContents(title="Photo-Catalytic Reactor")

# ╔═╡ 94ed0f8b-13c1-4460-a391-5057cff401e0
md"""
Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=true))
"""

# ╔═╡ 2ed3223e-a604-410e-93d4-016580f49093
md"""
# Domain / Grid
"""

# ╔═╡ 390c7839-618d-4ade-b9be-ee9ed09a77aa
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The modelling domain covers a prism of square shape, which is 1/4 of the symmetric porous foam that is supporting the catalyst layer in the photo-catalytic reactor.
$(LocalResource("../img/Domain.png")) 

"""
  ╠═╡ =#

# ╔═╡ 0a911687-aff4-4c77-8def-084293329f35
begin
	const Γ_side_front = 1 # symmetry bc
	const Γ_side_right = 2 # wall bc
	const Γ_side_back = 3 # wall bc
	const Γ_side_left = 4 # symmetry bc
	const Γ_bottom = 5 # inflow bc
	const Γ_top_frit = 6 # outflow bc, uncoated porous frit 
	const Γ_top_cat = 7 # outflow bc, catalyst coated porous frit 
end;

# ╔═╡ ada45d4d-adfa-484d-9d0e-d3e7febeb3ef
function prism_sq(data; nref=0, w=data.wi, h=data.h, cath=data.cath, catwi=data.catwi)
	
	hw=(w/2.0)/5.0*2.0^(-nref)
	W=collect(0:hw:(w/2.0))
	
	zCL=h-cath
	hhfrit=zCL/5.0
	Hfrit=collect(0:hhfrit:zCL)
	
	hhCL=cath/5.0
	HCL=collect(zCL:hhCL:h)
	H=glue(Hfrit,HCL)
	
	grid=simplexgrid(W,W,H)
	
	# catalyst layer region
	cellmask!(grid,[0.0,0.0,h-cath],[catwi/2,catwi/2,h],2)
	# catalyst layer boundary
	bfacemask!(grid,[0.0,0.0,h],[catwi/2,catwi/2,h],Γ_top_cat)	
end

# ╔═╡ e21b9d37-941c-4f2c-9bdf-956964428f90
const grid_fun = prism_sq

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
- \nabla N_i + R_i &= 0
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

# ╔═╡ 8139166e-42f9-41c3-a360-50d3d4e5ee86
md"""
## Gas species
"""

# ╔═╡ 44d91c2e-8082-4a90-89cc-81aba783d5ac
md"""
Boundary conditions for the transport of gas phase species cover in and outflow boundary conditions at the bottom and top surfaces of the modelling domain with no-flux conditins applied elsewhere. In the catalyst layer, volumetric catalytic reactions take place.
"""

# ╔═╡ 58d0610b-1739-4260-8d16-5a31ba362d69
md"""
## Irradiation Exchange:
Account for the exchange of irradiation between the surfaces in the top and bottom chambers of the reactor. View factors are used to describe the geometrical arrangement of surfaces that exchange irradiation. The view factors are values between 0 and 1, describing the fractions of surface areas that "see each other" (exchange irradiation with each other). In upper chamber, the surfaces are designated as (see figure below):
1. quartz window
2. catalyst layer
3. uncoated frit

For the given geometry, the following view factors are obtained:

```math
\phi_{12} = \frac{A_2}{A_1}, \qquad \phi_{13} = \frac{A_3}{A_1}, \qquad \phi_{21} = \phi_{31} = 1
```
where $A_i$ are the geometric areas of the corresponding surfaces.
"""

# ╔═╡ d9dee38e-6036-46b8-bc06-e545baa06789
md"""
In the following, the net irradiation fluxes through the surfaces of interest marked by the __red__ dashed line (__catalyst layer__, __2__) and __green__ dashed line (__uncoated frit__, __3__) are stated. Per sign convention applied here, __fluxes entering__ the control surfaces are counted as __positive__ while __fluxes exiting__ the control surfaces are counted as __negative__. 
These irradiation fluxes through the surfaces will be implemented in the code as boundary conditions for the thermal energy transport equation.

"""

# ╔═╡ e58ec04f-023a-4e00-98b8-f9ae85ca506f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Top chamber
$(LocalResource("../img/irrad_top_chamber_view_factors.png"))
```math
\begin{align}
F_2 &=-\epsilon_2\sigma T_2^4 + \alpha_2^{\text{vis}} G_1^{\text{vis}} + \alpha_2^{\text{IR}} G_1^{\text{IR}} \\
F_3 &=-\epsilon_3\sigma T_3^4 + \alpha_3^{\text{vis}} G_1^{\text{vis}} + \alpha_3^{\text{IR}} G_1^{\text{IR}}
\end{align}
```
"""
  ╠═╡ =#

# ╔═╡ 80d2b5db-792a-42f9-b9c2-91d0e18cfcfb
md"""
```math
G_1^{\text{vis}} = \frac{\tau_1^{\text{vis}} G_{\text{lamp}}}{1-\rho_1^{\text{vis}} \left( \phi_{12} \rho_2^{\text{vis}} + \phi_{13} \rho_3^{\text{vis}} \right)} \\
```
"""

# ╔═╡ b1ae3b4d-59ca-420f-a1a0-dc698b52e5b0
md"""
```math
G_1^{\text{IR}} = \frac{\epsilon_1\sigma T_1^4 + \rho_1^{\text{IR}} \left( \phi_{12}\epsilon_2 \sigma T_2^4+\phi_{13}\epsilon_3 \sigma T_3^4\right)}{1-\rho_1^{\text{IR}} \left( \phi_{12} \rho_2^{\text{IR}} + \phi_{13} \rho_3^{\text{IR}} \right)}
```
"""

# ╔═╡ b1ef2a89-27db-4f21-a2e3-fd6356c394da
md"""
where $G_1^{\text{vis}}$ and $G_1^{\text{IR}}$ are the surface brightnesses (cumulative irradiation flux resulting from emission, transmission and reflection) of surface 1 in the visible and IR spectral range respectively. The temperature of the quartz window $T_1$ is calculated iteratively after each simulation iteration.
"""

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
	alpha_IR=0.7, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.7, # measurement (ideally at high T) integrated over vis
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
Calculate window temperature that is used in the temperature boundary condition of the top chamber. It can be obtained from an energy balance in thermal equilibrium. The value of window temperature, which is fixed stored in the data struct, is updated after every simulation run.
"""

# ╔═╡ a395374d-5897-4764-8499-0ebe7f2b4239
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
$(LocalResource("../img/WindowEB.png")) 
"""
  ╠═╡ =#

# ╔═╡ 7db215d7-420e-435b-8889-43c4fff150e0
md"""
At steady state, the net transfer of energy to/from the balance region around the quartz window must vanish. This leads to:

```math
0=-\dot q_{\text{conv,top}} - 2 \dot q_{\text{emit}} + \dot q_{\text{cond,bot,2}} + \dot q_{\text{cond,bot,3}} + \dot q_{\text{abs,bot,2}} + \dot q_{\text{abs,bot,3}}
```

with
```math
	\begin{align}
	\dot q_{\text{conv,top}} &= \alpha_{\text{conv}} (T_1-T_{\text{amb}}) \\
	\dot q_{\text{emit}} &= \epsilon_1\sigma T_1^4 \quad(\text{once for top \& bottom surfaces}) \\
	\dot q_{\text{cond,bot,2}} &= -\lambda(\frac{T_1+T_2}{2})\frac{T_1-T_2}{\Delta z} \\
	\dot q_{\text{cond,bot,3}} &= -\lambda(\frac{T_1+T_3}{2})\frac{T_1-T_3}{\Delta z} \\
	\dot q_{\text{abs,bot,2}} &= \alpha_1^{\text{IR}} \phi_{12} G_2^{\text{IR}} + \alpha_1^{\text{vis}} \phi_{12} G_2^{\text{vis}} \\
	\dot q_{\text{abs,bot,3}} &= \alpha_1^{\text{IR}} \phi_{13} G_3^{\text{IR}} + \alpha_1^{\text{vis}} \phi_{13} G_3^{\text{vis}}
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
where the temperature of the bottom Al plate $T_2$ is calculated iteratively outside of the simulation and then set as a parameter in the ModelData structure.
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
0=-\dot q_{\text{conv}} - 2 \dot q_{\text{emit}} - \dot q_{\text{cond}} + \dot q_{\text{abs,1}} 
```

with
```math
	\begin{align}
	\dot q_{\text{conv}} &= \alpha_{\text{conv}} (T_2-T_{\text{amb}}) \\

	\dot q_{\text{emit}} &= \epsilon_2\sigma T_2^4 \quad(\text{once for top \& bottom surfaces}) \\
	\dot q_{\text{cond}} &= -\lambda(\frac{T_1+T_2}{2})\frac{T_1-T_2}{\Delta z} \\
	
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

# ╔═╡ 04058f1d-9622-4b59-bf07-93483bc269f2
md"""
## Conductive Transport
Becuase the precise flow field is not known in the chambers, the transport of heat through the gas volumes in the chambers is estimated with conductive transport through a stagnant layer of gas. This is the lower bound on convective transport. It can be estimated assuming a constant temperature gradient over the chamber height and an thermal conductivity evaluated at an average temperature:

```math
\dot q = -\lambda \nabla T \approx - \lambda(\overline T) \frac{\Delta T}{\Delta z}
```
"""

# ╔═╡ 4ebbe06f-0993-4c5c-9af3-76b2b645e592
md"""
## Implementation
"""

# ╔═╡ 0bc79692-8db4-44a2-9433-5b6ce97b656f
md"""
#### Practical problem
A problem regarding the implementation of above irradiation flux boundary conditions in VoronoiFVM.jl arises: in the physics callback functions that define the boundary conditions, only local values of the solution are available. But since irradiation is not acting locally but acts over distances, values from other part of the modeling domain influence the local solution. These interactions cannot be considered in a streight-forward way.

In the following, a "simplification" is applied: instead of considering the different temperatures $T_2$ and $T_3$, we consider the single, local temperature $T_{2/3}$ in the calculation of the surface brightness. 
Instead of the surface brightness $G_1^{\text{IR}}$ as defined above, the following is implemented in the code:

```math
G_1^{\text{IR}} = \frac{\epsilon_1\sigma T_1^4 + \rho_1^{\text{IR}} \left( \phi_{12}\epsilon_2 \sigma T_{\bf{2/3}}^4+\phi_{13}\epsilon_3 \sigma T_{\bf{2/3}}^4 \right)}{1-\rho_1^{\text{IR}} \left( \phi_{12} \rho_2^{\text{IR}} + \phi_{13} \rho_3^{\text{IR}} \right)}

```
This will lead to an error, that should be estimated.

"""

# ╔═╡ edd9fdd1-a9c4-4f45-8f63-9681717d417f
function side(f,u,bnode,data)
	# side wall boundary condition
	(;iT,α_w,Tamb)=data
	# outer side wall boundaries
	if bnode.region==Γ_side_back || bnode.region==Γ_side_right
		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative
		f[iT] = α_w*(u[iT]-Tamb)
	end
	#iT=data.iT
	#boundary_robin!(f,u,bnode;species=iT,region=[Γ_side_back,Γ_side_right], factor=data.α_w, value=data.Tamb*data.α_w)	
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
#Base.@kwdef mutable struct ModelData <:AbstractModelData

begin
	# !!!ALLOC We need ng as a "static" parameter (known at compile time)
	# !!!ALLOC so we make it a type parameter instead of a struct entry
	Base.@kwdef mutable struct ModelData{NG} <:AbstractModelData
	
	# catalyst / chemistry data
	# kinetic parameters, S3P="simple 3 parameter" kinetics fit to UPV lab scale experimental data
	kinpar::FixedBed.KinData{nreac(S3P)} = S3P
	#kinpar::AbstractKineticsData = XuFroment1989
	
	# number of gas phase species
	#ng::Int64		 		= kinpar.ng
	ip::Int64=NG+1 # index of total pressure variable
	iT::Int64=ip+1 # index of Temperature variable

	# names and fluid indices
	gn::Dict{Int, Symbol} 	= kinpar.gn

	# inverse names and fluid indices
	gni::Dict{Symbol, Int}  = kinpar.gni
	# fluids and respective properties in system
	Fluids::Vector{FluidProps} = kinpar.Fluids
	#Fluids::Vector{AbstractFluidProps} = [N2]
	X0::Vector{Float64} = let
		x=zeros(Float64, NG)
		x[gni[:H2]] = 1.0
		x[gni[:CO2]] = 1.0
		x/sum(x)
	end # inlet composition

	mcat::Float64=3.3*ufac"g" # total catalyst loading
	#mcat::Float64=200.0*ufac"mg" 
	# volume specific cat mass loading, UPV lab scale PC reactor
	lcats::Float64 =1000.0*ufac"kg/m^3"
	#mcats::Float64=80.0*ufac"kg/m^3" # 200 mg cat total loading, 250μm CL 
	isreactive::Bool = 1
	#isreactive::Bool = 0
		
	α_w::Float64=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
	α_nat_conv::Float64=15.0*ufac"W/(m^2*K)" # natural convection heat transfer coefficient	
	
	## porous filter data
	dp::Float64=200.0*ufac"μm" # average pore size, por class 0
	#dp::Float64=100.0*ufac"μm" # average pore size, por class 2

	# frit thickness (applies to 2D & 3D)
	h::Float64=0.5*ufac"cm"
	# catalyst layer thickness (applies to 2D & 3D)
	#cath::Float64 = 500.0*ufac"μm"
	cath::Float64 = 250.0*ufac"μm"

	# upper and lower chamber heights for calculation of conduction b.c.
	uc_h::Float64=17.0*ufac"mm"
	lc_h::Float64=18.0*ufac"mm"

	# prism / 3D
	wi::Float64=12.0*ufac"cm" # prism width/side lenght
	le::Float64=wi # prism width/side lenght
	catwi::Float64=10.0*ufac"cm" # prism width/side lenght	
	
	Ac::Float64=wi*le*ufac"m^2" # cross-sectional area, square

	## irradiation data
	Glamp_target::Float64=100.0*ufac"kW/m^2" # solar simulator irradiation flux
	Glamp::Float64=1.0*ufac"kW/m^2" # solar simulator irradiation flux

	# upper chamber: quartz window eff. optical properties
	uc_window::SurfaceOpticalProps = uc_window
	# upper chamber: catalyst layer eff. optical properties
	uc_cat::SurfaceOpticalProps = uc_cat
	# upper chamber: (uncoated) porous frir eff. optical properties
	uc_frit::SurfaceOpticalProps = uc_frit

	#view factors
	vf_uc_window_cat::Float64 = catwi^2/Ac # A2/A1
	vf_uc_window_frit::Float64 = 1-vf_uc_window_cat # A3/A1 = 1-A2/A1

	# lower chamber: (uncoated) porous frir eff. optical properties
	lc_frit::SurfaceOpticalProps = lc_frit
	# lower chamber:  Al bottom plate eff. optical properties
	lc_plate::SurfaceOpticalProps = lc_plate
	
	Tglass::Float64 = 373.15*ufac"K" # glass temperature on surface facing catalyst
	Tplate::Float64 = 373.15*ufac"K" # bottom Al plate temperature facing frit
	

	# Solid Boro-Silikatglas
	ρs::Float64=2.23e3*ufac"kg/m^3" # density of non-porous Boro-Solikatglas 3.3
	#λs::Float64=1.4*ufac"W/(m*K)" # thermal conductiviy of non-porous SiO2
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


	## Flow data
	#norm conditions
	pn::Float64 = 1.0*ufac"bar"
	Tn::Float64 = 273.15*ufac"K"

	#Qflow::Float64=340.0*ufac"ml/minute" # volumetric feed flow rate (sccm)
	#Qflow::Float64=1480.0*ufac"ml/minute" # volumetric feed flow rate (sccm)
	Qflow::Float64=3400.0*ufac"ml/minute" # volumetric feed flow rate (sccm)
	

	MWin::Float64 = molarweight_mix(Fluids, X0)
	mdotin::Float64=MWin*Qflow*pn/(ph"R"*Tn)*ufac"kg/s"
	
	Tamb::Float64=298.15*ufac"K" # ambient temperature
	p::Float64=1.0*ufac"atm" # reactor pressure

	u0::Float64=Qflow/(Ac)*ufac"m/s" # mean superficial velocity
	ubot::Float64=u0*ufac"m/s" # adjustable parameter to match the outlet mass flow to prescribed inlet mass flow rate
	
end;
	
	# !!!ALLOC Method to be called instead of data.ng
	FixedBed.ngas(::ModelData{NG}) where NG = NG
	
	# !!!ALLOC Additional constructo taking ng as parameter	
	ModelData(;ng=S3P.ng, kwargs...) = ModelData{ng}(;kwargs...)

end;

# ╔═╡ 78cf4646-c373-4688-b1ac-92ed5f922e3c
function reaction(f,u,node,data)
	(;ip,isreactive)=data
	ng=ngas(data)
	
	if node.region == 2 && data.isreactive # catalyst layer
		(;iT,lcats,kinpar,Fluids)=data
		(;rni,nuij)=kinpar
		pi = MVector{ngas(data),eltype(u)}(undef)
		

		for i=1:ng
			pi[i] = u[i]/ufac"bar"
		end
		# negative sign: sign convention of VoronoiFVM: source term < 0 
		# for Xu/Froment 1989 kinetics, ri returns reaction rate in mol/(h gcat)
		unitc=one(eltype(u))
		if kinpar == XuFroment
			unitc *=ufac"mol/(hr*g)"
		end

		#RR = @inline -mcats*ri(data,u[iT],pi)*ufac"mol/(hr*g)"
		RR = @inline -lcats*ri(data,u[iT],pi)*unitc
		for i=1:ng
			f[i] = zero(eltype(u))
			for j=1:nreac(kinpar)
				f[i] += nuij[i,j] * RR[j]
			end			
		end

	end
	
	# ∑xi = 1
	@views f[ip]=u[ip]-sum(u[1:ng])
	
end

# ╔═╡ 7da59e27-62b9-4b89-b315-d88a4fd34f56
function top(f,u,bnode,data)
	# top boundaries (cat layer & frit)
	if bnode.region==Γ_top_frit || bnode.region==Γ_top_cat 
		(;iT,Glamp,Tglass,X0,Fluids,Tamb,uc_window,uc_cat,uc_frit,vf_uc_window_cat,vf_uc_window_frit,uc_h)=data
		ng=ngas(data)
		# convective therm. energy flux
		flux_entahlpy=zero(eltype(u))
		
		# flow velocity is normal to top boundary
		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative
		for i=1:ng
			# f[i] = data.utop*u[i]/(ph"R"*u[iT])
			f[i] = -data.u0*data.X0[i]*data.pn/(ph"R"*data.Tn)
			# use species enthalpy (incl. Δh_formation) for conv. therm. eng. flux 1/2
			# reactant gases enter control volume at Tref
			@inline flux_entahlpy += f[i] * enthalpy_gas(Fluids[i], Tamb)
			#f[i] = -data.u0*data.X0[i]*data.p/(ph"R"*data.Tamb)

		end

		
		#flux_irrad=0.0
		flux_irrad=zero(eltype(u))
		σ=ph"σ"

		# irradiation exchange between quartz window (1), cat surface (2) & frit surface (3)
		# window properties (1)
		tau1_vis=uc_window.tau_vis
		rho1_vis=uc_window.rho_vis
		rho1_IR=uc_window.rho_IR
		eps1=uc_window.eps
		# catalyst layer properties (2)
		rho2_vis=uc_cat.rho_vis
		rho2_IR=uc_cat.rho_IR
		alpha2_vis=uc_cat.alpha_vis
		alpha2_IR=uc_cat.alpha_IR
		eps2=uc_cat.eps
		# uncoated frit properties (3)
		rho3_vis=uc_frit.rho_vis
		rho3_IR=uc_frit.rho_IR
		alpha3_vis=uc_frit.alpha_vis
		alpha3_IR=uc_frit.alpha_IR
		eps3=uc_frit.eps
		#view factors
		ϕ12=vf_uc_window_cat
		ϕ13=vf_uc_window_frit
		
		# surface brigthness of quartz window (1) in vis & IR	
		G1_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))
		# here the simplification is applied: only local value of T (u[iT]) is available, it is used for both surfaces
		G1_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
		
		if bnode.region==Γ_top_cat
			flux_irrad = -eps2*σ*u[iT]^4 + alpha2_vis*G1_vis + alpha2_IR*G1_IR
		else # upper chamber frit
			flux_irrad = -eps3*σ*u[iT]^4 + alpha3_vis*G1_vis + alpha3_IR*G1_IR
		end


		# conductive heat flux through top chamber
		# mean temperature
        Tm=0.5*(u[iT] + Tglass)
        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X0)

        flux_cond = -λf*(Tglass-u[iT])/uc_h
		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative

		# include convective therm. energy flux
		# use species enthalpy (incl. Δh_formation) for conv. therm. eng. flux 2/2
		f[iT] = -flux_irrad + flux_cond + flux_entahlpy
		#f[iT] = -flux_irrad + flux_cond
		
		

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
		(;iT,ip,Fluids,ubot,Tamb,lc_frit,lc_plate,Tplate,lc_h) = data
		ng=ngas(data)

		X=MVector{ngas(data),eltype(u)}(undef)
		@inline mole_frac!(bnode,data,X,u)
		#X=zeros(eltype(u), ng)
		#mole_frac!(bnode,data,X,u)

		# use species enthalpy (incl. Δh_formation) for conv. therm. eng. flux 1/2
		#hmix=enthalpy_mix(Fluids, u[iT], X)
		@inline hmix=enthalpy_mix(Fluids, u[iT], X)
		flux_entahlpy=ubot*u[ip]/(ph"R"*u[iT])*hmix
		#cf=heatcap_mix(Fluids, u[iT], X)		
		#flux_entahlpy=ubot*u[ip]/(ph"R"*u[iT])*cf*(u[iT]-Tamb)

		# irradiation exchange between porous frit (1) and Al bottom plate (2)
		# porous frit properties (1)
		eps1=lc_frit.eps; alpha1_IR=lc_frit.alpha_IR; rho1_IR=lc_frit.rho_IR; 	
		# Al bottom plate properties (2)
		eps2=lc_plate.eps; rho2_IR=lc_plate.rho_IR;
		
		#(;eps1,alpha1_IR,rho1_IR,rho2_IR,eps2,rho2_IR) = lower_chamber
		σ=ph"σ"
		flux_irrad = -eps1*σ*u[iT]^4 + alpha1_IR/(1-rho1_IR*rho2_IR)*(eps2*σ*Tplate^4+rho2_IR*eps1*σ*u[iT]^4)

		# conductive heat flux through top chamber
		# mean temperature
        Tm=0.5*(u[iT] + Tplate)
        # thermal conductivity at Tm and outlet composition X
        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X)

        flux_cond = -λf*(u[iT]-Tplate)/lc_h # positive flux in positive z coord.
		
		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative
		f[iT] = -flux_irrad + flux_entahlpy - flux_cond

		
		
		for i=1:ng
			# specify flux at boundary: flow velocity is normal to bot boundary
			f[i] = ubot*u[i]/(ph"R"*u[iT])			
		end
	end
end

# ╔═╡ 29d66705-3d9f-40b1-866d-dd3392a1a268
function bcond(f,u,bnode,data)
	top(f,u,bnode,data)
	bottom(f,u,bnode,data)
	side(f,u,bnode,data)	
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
		#M[i,i] = -1/DK_eff(data,T,i)
		M[i,i] = 1/DK_eff(data,T,i)
		for j=1:ng
			if j != i
				#M[i,i] -= x[j]/(data.γ_τ*D[i,j])
				M[i,i] += x[j]/(data.γ_τ*D[i,j])
				#M[i,j] = x[i]/(data.γ_τ*D[i,j])
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
	# !!!ALLOC Use MVector with static size information instef of Vector
	F = MVector{ngas(data),eltype(u)}(undef)
	X = MVector{ngas(data),eltype(u)}(undef)

	pk,pl = u[ip,1],u[ip,2]
	δp = pk-pl
	pm = 0.5*(pk+pl)
	Tm=0.5*(u[iT,1]+u[iT,2])
	
	# !!!ALLOC All functions which take MVectors to be stack allocated
	# as argument must be inlined. Here, we use the new "callsite inline"
	# wich is new with Julia 1.8.
        
	@inline mole_frac!(edge,data,X,u)
	@inline cf=heatcap_mix(Fluids, Tm, X)
	@inline μ,λf=dynvisc_thermcond_mix(data, Tm, X)
 	λbed=kbed(data,λf)*λf
	
	#λbed=lambda_eff_AC(data,λf)

	# Darcy flow
	ud=-k/μ * δp
	
	# vh=project(edge,(0,ud)) # 2D
	#vh=project(edge,(0,0,ud)) # 3D
	vh=project(edge,(zero(eltype(u)),zero(eltype(u)),ud)) # 3D

	# !!!ALLOC Use MMatriy with static size information instef of Matrix
	M = MMatrix{ngas(data),ngas(data),eltype(u)}(undef)
	D = MMatrix{ngas(data),ngas(data),eltype(u)}(undef)
	@inline M_matrix!(M,D,Tm,pm,X,data)
   	
		
	for i=1:ng
		DK = DK_eff(data,Tm,i)
		bp,bm=fbernoulli_pm(vh/DK)		
		F[i] = (bm*u[i,1]-bp*u[i,2])/(ph"R"*Tm)
	end
	
	# computation of fluxes J

	# !!!ALLOC All functions which take MMatriy to be stack allocated
	# as argument must be inlined. Here, we use the new "callsite inline"
	# wich is new with Julia 1.8
	# inplace_linsolve  is made available from VoronoiFVM 1.3.2
	@inline inplace_linsolve!(M,F)
	
	@views f[1:ng] .= F	

	#f[iT]= λbed*(Bm*(u[iT,1]-Tamb)-Bp*(u[iT,2]-Tamb))		

		# use species molar fluxes for thermal drift contribution
	conv=zero(eltype(u))
	for i=1:ng
		# use species enthalpy (incl. Δh_formation) for conv. therm. eng. flux 1/3
		#conv += J[i] * enthalpy_gas(Fluids[i], T)
		conv += F[i] * enthalpy_gas(Fluids[i], Tm)
	end
		# use species enthalpy (incl. Δh_formation) for conv. therm. eng. flux 2/3
	Bp,Bm = fbernoulli_pm(conv/λbed/Tm)
		# thermal energy flux
		# use species enthalpy (incl. Δh_formation) for conv. therm. eng. flux 3/3
	f[iT]= λbed*(Bm*u[iT,1]-Bp*u[iT,2])		
		
	#f[ip] via reaction: ∑pi = p
end

# ╔═╡ b2df1087-6628-4889-8cd6-c5ee7629cd93
md"""
## Temperature Plot
"""

# ╔═╡ 2790b550-3105-4fc0-9070-d142c19678db
md"""
## Partial Pressure Plot
"""

# ╔═╡ 31add356-6854-43c5-9ebd-ef10add6cc3d
md"""
### CO molar fraction in Cross-section
"""

# ╔═╡ bd7552d2-2c31-4834-97d9-ccdb4652242f
function SolAlongLine(data,sol)
	ng=ngas(data)		
	grid=grid_fun(data)
	
	mid_x=argmin(abs.(grid[Coordinates][1,:] .-data.wi/4))
	mid_x=grid[Coordinates][1,mid_x]
	mid_y=argmin(abs.(grid[Coordinates][2,:] .-data.le/4))
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
	# bottom - inflow
	Ibot=integrate(sys,bottom,sol; boundary=true)[:,Γ_bottom]
	# top: outflow
	# uncoated outer frit area, inner cat coated area
	Itop=integrate(sys,top,sol; boundary=true)[:,[Γ_top_frit,Γ_top_cat]] 
	Itop=sum(Itop, dims=2)
	Ibot[1:ng],Itop[1:ng]
end

# ╔═╡ ed0c6737-1237-40d0-81df-32d6a842cbb0
function cat_mass(data)
	data.mcat
end

# ╔═╡ 7f2276d5-bb78-487c-868a-0486595d0113
function avg_RR(sol,sys,data)
	# sign convention: source term <0
	PRs=-integrate(sys,reaction,sol;)[:,2] # mol/s
	# catalyst mass specific prodcution rate of CO

	PRs[1] / cat_mass(data) # mol/ (s kgcat)
end

# ╔═╡ a6e61592-7958-4094-8614-e77446eb2223
md"""
##  Average Catalyst Temperature
"""

# ╔═╡ 808aed68-7077-4079-be75-1bea962c716d
function TmaxSide(sol,sys,grid,data)

    (;iT)=data
	sub=subgrid(grid,[Γ_side_back,Γ_side_right],boundary=true)

    Tsides = view(sol[iT, :], sub) .- 273.15
    maximum(Tsides)
		
	
end

# ╔═╡ a55c5ee7-2274-4447-b0b2-58052f064bc9
function areas(sol,sys,grid,data)
	iT = data.iT
	function area(f,u,bnode,data)
		# repurpose temperature index to hold area information
		f[iT] = one(eltype(u))
	end

	integrate(sys,area,sol; boundary=true)[iT,:]
end

# ╔═╡ 333b5c80-259d-47aa-a441-ee7894d6c407
function main(;data=ModelData(),nref=0)

	grid=grid_fun(data;nref=nref)
	(;iT,p,X0,Tamb)=data
	ng=ngas(data)
	
	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)
	enable_species!(sys; species=collect(1:(ng+2))) # gas phase species + p + T
	
	inival=unknowns(sys)
	inival[:,:].=1.0*p
	for i=1:ng
		inival[i,:] .*= X0[i]
	end
	inival[iT,:] .= Tamb

	sol=solve(sys;inival,verbose="an")

	# update view factors as they manifest in the computational grid
	# for coarse grids, the computed view-factors might deviate from the geometric values
	A=areas(sol,sys,grid,data)
	Acat=A[Γ_top_cat]
	Afrit=A[Γ_top_frit]
	data.vf_uc_window_cat = Acat/sum(Acat+Afrit)
	data.vf_uc_window_frit = Afrit/sum(Acat+Afrit)

	function WindowTemperature_(sol,sys,data)
		(;iT,Glamp,X0,α_nat_conv,Tamb,uc_window,uc_cat,uc_frit,vf_uc_window_cat,vf_uc_window_frit,uc_h)=data
		σ=ph"σ"
		
		# irradiation exchange between quartz window (1), cat surface (2) & frit surface (3)
		# window properties (1)
		tau1_vis=uc_window.tau_vis
		rho1_vis=uc_window.rho_vis
		rho1_IR=uc_window.rho_IR
		alpha1_IR=uc_window.alpha_IR
		alpha1_vis=uc_window.alpha_vis
		eps1=uc_window.eps
		# catalyst layer properties (2)
		rho2_vis=uc_cat.rho_vis
		rho2_IR=uc_cat.rho_IR
		alpha2_vis=uc_cat.alpha_vis
		alpha2_IR=uc_cat.alpha_IR
		eps2=uc_cat.eps
		# uncoated frit properties (3)
		rho3_vis=uc_frit.rho_vis
		rho3_IR=uc_frit.rho_IR
		alpha3_vis=uc_frit.alpha_vis
		alpha3_IR=uc_frit.alpha_IR
		eps3=uc_frit.eps
		#view factors
		ϕ12=vf_uc_window_cat
		ϕ13=vf_uc_window_frit
	
		function f!(F, x)		
			Tglass=x[1]
	
			function flux_abs_top_catalyst_layer(f,u,bnode,data)
				if bnode.region==Γ_top_cat
					
					G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))
			
					G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))	        		
			
					G2_vis = rho2_vis*G1_bot_vis		
					
					G2_IR = eps2*σ*u[iT]^4 + rho2_IR*G1_bot_IR
			
					f[iT] = alpha1_vis*G2_vis + alpha1_IR*G2_IR
				end
			end
	
			function flux_abs_top_frit(f,u,bnode,data)
				if bnode.region==Γ_top_frit
					
					G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))
					
					G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
			
					G3_vis = rho3_vis*G1_bot_vis
	
					G3_IR = eps3*σ*u[iT]^4 + rho3_IR*G1_bot_IR
			
					f[iT] = alpha1_vis*G3_vis + alpha1_IR*G3_IR
				end
			end
	
			function flux_conduction_top(f,u,bnode,data)
				if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
			
					# mean temperature
					Tm=0.5*(u[iT] + Tglass)
			
					# thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
					@inline _,λf=dynvisc_thermcond_mix(data, Tm, X0)
			
					# positive flux in positive z-coordinate
					q_cond = -λf*(Tglass-u[iT])/uc_h			
					f[iT] = q_cond			
				end
			end		
			
			Qabs_10 = 4*integrate(sys,flux_abs_top_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
			Qabs_20 = 4*integrate(sys,flux_abs_top_frit,sol; boundary=true)[iT,Γ_top_frit]
	
	
			Qcond_10=4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_cat]
			Qcond_20=4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_frit]
	
			Atop = sum(areas(sol,sys,grid,data)[[Γ_top_cat,Γ_top_frit]])
	
			Qconv0 = 4*Atop*α_nat_conv*(Tglass-Tamb)
			Qemit0 = 4*Atop*uc_window.eps*σ*Tglass^4		
	
			F[1] = -Qconv0 +Qcond_10 +Qcond_20 - 2*Qemit0 +Qabs_10 +Qabs_20
			
		end
	
		sol = nlsolve(f!, [573.15], autodiff=:forward)
		sol.zero[1]
		
	end

	function PlateTemperature_(sol,sys,data)
		(;iT,gni,α_nat_conv,Tamb,lc_frit,lc_plate,lc_h)=data
	
		σ=ph"σ"
		
		# irradiation exchange between frit surface (1) & Al bottom plate (2)
		# frit properties
		eps1=lc_frit.eps
		rho1_IR=lc_frit.rho_IR
		# plate properties
		alpha2_IR=lc_plate.alpha_IR
		eps2=lc_plate.eps
		rho2_IR=lc_plate.rho_IR	
		
		
		function f!(F, x)
			Tplate=x[1]

			function flux_abs_bot(f,u,bnode,data)
			    if bnode.region==Γ_bottom
					
					eps1=lc_frit.eps;  rho1_IR=lc_frit.rho_IR; 	
					# Al bottom plate properties (2)
					eps2=lc_plate.eps; rho2_IR=lc_plate.rho_IR;
					
			
			        G1_IR = (eps1*σ*u[iT]^4 + rho1_IR*eps2*σ*Tplate^4)/(1-rho1_IR*rho2_IR)
			
			        f[iT] = alpha2_IR*G1_IR
			    end
			end

			Qabs_34 = 4*integrate(sys,flux_abs_bot,sol; boundary=true)[iT,Γ_bottom]

			# conductive heat flux through bottom chamber
			function flux_conduction_bottom(f,u,bnode,data)
			    if bnode.region==Γ_bottom
			        			
			        X=zeros(eltype(u), ng)
					mole_frac!(bnode,data,X,u)
			
			        Tm=0.5*(u[iT] + Tplate)
			        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X)
			
			        # positive flux in negative z coord. -> pointing towards bottom plate
			        q_cond = -λf*(Tplate-u[iT])/lc_h 			
			        f[iT] = q_cond			
			    end
			end

			Qcond_34 = 4*integrate(sys,flux_conduction_bottom,sol; boundary=true)[iT,Γ_bottom]
			


			Abot = areas(sol,sys,grid,data)[Γ_bottom]
	
			Qconv4 = 4*Abot*α_nat_conv*(Tplate-Tamb)
			Qemit4 = 4*Abot*lc_plate.eps*σ*Tplate^4		
	
			F[1] = -Qconv4 +Qcond_34 -2*Qemit4 +Qabs_34
		end
	
		sol = nlsolve(f!, [573.15], autodiff=:forward)
		sol.zero[1]
	end

	
	function pre(sol,par)
		ng=ngas(data)
		# set glass emperature
		data.Tglass = WindowTemperature_(sol,sys,data)

		# set bottom plate temperature
		data.Tplate = PlateTemperature_(sol,sys,data)
		
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
		(;catwi,cath,mcat,Glamp_target)=data
		Vcat = catwi^2*cath # m^3
		lcat = mcat/Vcat # kg/m^3"
		lcats=[10.0, lcat]*ufac"kg/m^3"
	 	data.lcats = minimum(lcats) + par*(maximum(lcats)-minimum(lcats))

	 	# irradiation flux density
	 	Glamps=[1.0*ufac"kW/m^2", Glamp_target]
	 	data.Glamp = minimum(Glamps) + par*(maximum(Glamps)-minimum(Glamps))
	 end


	
	 control=SolverControl( ;
	 				  		handle_exceptions=true,
							Δp_min=1.0e-4,					  
	 				  		Δp=0.1,
	 				  		Δp_grow=1.5,
	 				  		Δu_opt=100000.0, # large value, due to unit Pa of pressure?
						verbose="an", # log allocations and newton convergence
						method_linear = KrylovJL_BICGSTAB(),
						precon_linear = A -> factorize(A, UMFPACKFactorization())
	 )
	
	
	embed=[0.0,1.0]
	sol=solve(sys;inival,embed,pre,control)
	sol_end_embed=sol(sol.t[end])
	
	
	sol_old=sol_end_embed
	ΔT = [Inf,Inf]
	
	MAX_ITER=5
	ΔT_MAX=0.1*ufac"K"
	# iteratively solve system and update glass and plate temperatures
	iter=1
	while iter <= MAX_ITER && !all(ΔT .< ΔT_MAX)
		Told=[data.Tglass,data.Tplate]
		sol=solve(sys;inival=sol_old)
		data.Tglass = WindowTemperature_(sol,sys,data)
		data.Tplate = PlateTemperature_(sol,sys,data)
		ΔT = abs.(Told .-[data.Tglass,data.Tplate])
		sol_old=sol
		@show iter
		@show ΔT
		iter += 1
	end	
	
	sol,grid,sys,data
end;

# ╔═╡ aa498412-e970-45f2-8b11-249cc5c2b18d
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	if RunSim
		sol_,grid,sys,data_embed=main(data=ModelData();nref=0);
		if sol_ isa VoronoiFVM.TransientSolution
			sol = copy(sol_(sol_.t[end]))
		else
			sol = copy(sol_)
		end
	end
end;
  ╠═╡ =#

# ╔═╡ 985718e8-7ed7-4c5a-aa13-29462e52d709
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Cutplane at ``z=`` $(@bind zcut PlutoUI.Slider(range(0.0,data_embed.h,length=101),default=data_embed.h,show_value=true)) m
"""
  ╠═╡ =#

# ╔═╡ ff58b0b8-2519-430e-8343-af9a5adcb135
# ╠═╡ skip_as_script = true
#=╠═╡
let
	vis=GridVisualizer(resolution=(600,400),zoom=1.9)
	gridplot!(vis, grid_fun(ModelData{S3P.ng}()), zplane=zcut)
	reveal(vis)
end
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

# ╔═╡ 3bd80c19-0b49-43f6-9daa-0c87c2ea8093
#=╠═╡
let
	iT=data_embed.iT
	#vis=GridVisualizer(resolution=(600,400), title="Temperature °C", Plotter=PyPlot)
    vis=GridVisualizer(resolution=(600,400), title="Temperature °C")
	
	scalarplot!(vis, grid, sol[iT,:].- 273.15, levels= 10, levelalpha=1.0, show=true)

	reveal(vis)
	#save("../img/out/Temeprature.png", reveal(vis), Plotter=PyPlot)
end
  ╠═╡ =#

# ╔═╡ bea97fb3-9854-411c-8363-15cbef13d033
#=╠═╡
let
	data=data_embed
	sol_, grid_,midx,midy = SolAlongLine(data,sol)
	ng=ngas(data)
	(;gn, ip)=data

	cols = distinguishable_colors(ng+1, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
	pcols = map(col -> (red(col), green(col), blue(col)), cols)

	p1=Plots.plot(title="Partial Pressures", size=(450,450), xguide="Height / mm", yguide="Pressure / bar",legend=:outertopright,)
	for i in 1:ng
		Plots.plot!(p1, grid_./ufac"mm", vec(sol_[i])./ufac"bar", label="$(gn[i])", lw=2, ls=:auto)
	end
	Plots.plot!(p1, grid_./ufac"mm", vec(sol_[ip])./ufac"bar", color=cols[ip], label="total p", lw=2)
	lens!(10*[0.45, .5], [0.14, 0.17], inset = (1, bbox(0.1, 0.15, 0.25, 0.25)))
	lens!(10*[0.45, .5], [0.34, 0.37], inset = (1, bbox(0.5, 0.15, 0.25, 0.25)))

	p2 = Plots.plot(xguide="Height / cm", yguide="Pressure / bar",legend=:bottomright)
	Plots.plot!(p2, grid_./ufac"cm", vec(sol_[ip])./ufac"bar", color=cols[ip], label="total p", lw=2)

	p=Plots.plot(p1, size=(500,300))
	#Plots.plot(p1,p2, layout=(1,2), size=(700,450))
	#Plots.savefig(p, "../img/out/pi_flip_lowcat.svg")
	
end
  ╠═╡ =#

# ╔═╡ 0c3eb801-f271-4cf4-a109-25a283a51779
#=╠═╡
data_embed.mdotin/ufac"kg/hr"/4
  ╠═╡ =#

# ╔═╡ 7ab38bc2-9ca4-4206-a4c3-5fed673557f1
#=╠═╡
begin
	ndot_bot, ndot_top = MoleFlows(sol,sys,data_embed)
end;
  ╠═╡ =#

# ╔═╡ 8b1a0902-2542-40ed-9f91-447bffa4290f
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
| $(data_embed.gn[1]) | $(abs(round(ndot_top[1]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_bot[1]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data_embed.gn[2]) | $(abs(round(ndot_top[2]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_bot[2]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data_embed.gn[3]) | $(abs(round(ndot_top[3]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_bot[3]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data_embed.gn[4]) | $(abs(round(ndot_top[4]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_bot[4]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data_embed.gn[5]) | $(abs(round(ndot_top[5]/ufac"mol/hr",sigdigits=2)))  |   $(round(ndot_bot[5]/ufac"mol/hr",sigdigits=2))   |  mol/hr  |
| $(data_embed.gn[6]) | $(abs(round(ndot_top[6]/ufac"mol/hr")))  |   $(round(ndot_bot[6]/ufac"mol/hr"))   |  mol/hr  |


"""
  ╠═╡ =#

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

# ╔═╡ b22387c4-cda5-424d-99b7-d7deda24c678
#=╠═╡
md"""
__Side wall maximum__ temperature: $(round(TmaxSide(sol,sys,grid,data_embed))) °C
"""
  ╠═╡ =#

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

# ╔═╡ 5de9edf7-6059-47d5-b917-e7491068ebdc
function X_avg_bottom(sol,sys,grid,data)

	(;ng)=data
	function X_avg_(f,u,bnode,data)
		if bnode.region==Γ_bottom
			X=zeros(eltype(u), ng)
			mole_frac!(bnode,data,X,u)					
			for i=1:ng
				f[i] = X[i]
			end
		end
	end

	areas_=areas(sol,sys,grid,data)
	
	X_int=integrate(sys,X_avg_,sol; boundary=true)[[1:ng...],Γ_bottom]
	X_int./areas_[Γ_bottom]
end

# ╔═╡ fec9ca6d-d815-4b50-bec7-f8fb3d8195ba
function TopPlane(data,sol)
	ng=ngas(data)
	grid=grid_fun(data;nref=0)
	wi=data.wi/2
	
	bid = maximum(grid[BFaceRegions])+1
	bfacemask!(grid, [0,0,data.h],[wi,wi,data.h],bid)

	# keep x-y coordinates of parent grid
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[2]
	end
	#grid_1D  = subgrid(grid, [bid], boundary=true, transform=_3to1) 
	grid_2D  = subgrid(grid, [bid], boundary=true, transform=_3to2) 

	sol_p = []
	for i=1:(ng+2)
		sol_i = view(sol[i, :], grid_2D)
		push!(sol_p, collect(sol_i))
	end

	sol_p, grid_2D
end

# ╔═╡ 746e1a39-3c9b-478f-b371-3cb8333e93b1
function CutPlane(data,sol)
	ng=ngas(data)
	grid=grid_fun(data;nref=0)
	wi=data.wi/2
	
	bid = maximum(grid[BFaceRegions])+1
	bfacemask!(grid, [0,0.024,0],[wi,0.024,data.h],bid)

	# keep x-z coordinates of parent grid
	function _3to2(a,b)
		a[1]=b[1]
		a[2]=b[3]
	end
	#grid_1D  = subgrid(grid, [bid], boundary=true, transform=_3to1) 
	grid_2D  = subgrid(grid, [bid], boundary=true, transform=_3to2) 

	sol_p = []
	for i=1:(ng+2)
		sol_i = view(sol[i, :], grid_2D)
		push!(sol_p, collect(sol_i))
	end

	sol_p, grid_2D
end

# ╔═╡ 60b475d1-219e-4037-bd36-e0671cfa1893
#=╠═╡
let
	data=data_embed
	(;gni,ip)=data
	sol_xz, grid_xz = CutPlane(data_embed,sol)
	scalarplot(grid_xz, sol_xz[gni[:CO]] ./sol_xz[ip], show=true)
end
  ╠═╡ =#

# ╔═╡ 30393c90-298c-412d-86ce-e36106613d35
#=╠═╡
let
	sol_xy, grid_xy = TopPlane(data_embed,sol)
	sol_xz, grid_xz = CutPlane(data_embed,sol)
	vis=GridVisualizer(layout=(1,2), resolution=(700, 300),)
	scalarplot!(vis[1,1], grid_xy, sol_xy[data_embed.iT] .-273.15, colormap= :inferno, show=true)
	scalarplot!(vis[1,2], grid_xz, sol_xz[data_embed.iT] .-273.15, aspect=4, colormap= :inferno, show=true)
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

	#sol_1D, grid_1D=planeTop(data,sol)

	function Tcat_(f,u,bnode,data)
		iT = data.iT		
		f[iT] = u[iT]		
	end

	catA, catV = CatDims(grid)
	
	
	Tcat_avg_surf=integrate(sys,Tcat_,sol; boundary=true)[data.iT,Γ_top_cat] / catA - 273.15
	Tcat_avg_vol=integrate(sys,Tcat_,sol; )[data.iT,2] / catV - 273.15
	Tcat_avg_surf, Tcat_avg_vol
end

# ╔═╡ 2739bcff-3fb0-4169-8a1a-2b0a14998cec
#=╠═╡
md"""
Average Catalyst layer __surface__ temperature: $(round(Tcatavg(sol,sys,grid,data_embed)[1])) °C
"""
  ╠═╡ =#

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
function STCefficiency(sol,sys,data,grid)
	ndot_bot,ndot_top = MoleFlows(sol,sys,data)
	# R1 = RWGS in S3P kinetic model
	# R2 = RWGS in Xu & Froment kinetic model
	(;gni,kinpar,Glamp,Ac) = data
	DHi = 0.0
	if kinpar == XuFroment1989
		DHi = kinpar.ΔHi[:R2]
	elseif kinpar == S3P
		DHi = kinpar.ΔHi[:R1]
	end
	#STC_Atot = ndot_bot[data.gni[:CO]] * -data.kinpar.ΔHi[:R2] / (data.Glamp*data.Ac/4)
	STC_Atot = ndot_bot[gni[:CO]] * -DHi / (Glamp*Ac/4)
	catA, _ = CatDims(grid)
	#STC_Acat = ndot_bot[data.gni[:CO]] * -data.kinpar.ΔHi[:R2] / (data.Glamp*catA)
	STC_Acat = ndot_bot[gni[:CO]] * -DHi / (Glamp*catA)
	STC_Atot, STC_Acat
end

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
# ╠═ada45d4d-adfa-484d-9d0e-d3e7febeb3ef
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
# ╟─39e74955-aab6-4bba-a1b8-b2307b45e673
# ╟─6798d5e2-b8c7-4f54-aa71-6ea1ccab78fb
# ╟─ed3609cb-8483-4184-a385-dca307d13f17
# ╟─8139166e-42f9-41c3-a360-50d3d4e5ee86
# ╟─44d91c2e-8082-4a90-89cc-81aba783d5ac
# ╟─58d0610b-1739-4260-8d16-5a31ba362d69
# ╟─d9dee38e-6036-46b8-bc06-e545baa06789
# ╟─e58ec04f-023a-4e00-98b8-f9ae85ca506f
# ╟─80d2b5db-792a-42f9-b9c2-91d0e18cfcfb
# ╟─b1ae3b4d-59ca-420f-a1a0-dc698b52e5b0
# ╟─b1ef2a89-27db-4f21-a2e3-fd6356c394da
# ╟─4dae93b2-be63-4ee9-bc1e-871a31ade811
# ╠═9547ed7c-3304-4c63-a6c1-5f84e0001c54
# ╠═9d191a3a-e096-4ad7-aae6-bbd63d478fa2
# ╠═2d51f54b-4cff-4253-b17f-217e3261f36d
# ╠═fdaeafb6-e29f-41f8-8468-9c2b03b9eed7
# ╟─2e631f58-6a4b-4c6d-86ad-b748fd2d463a
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
# ╟─04058f1d-9622-4b59-bf07-93483bc269f2
# ╟─4ebbe06f-0993-4c5c-9af3-76b2b645e592
# ╟─0bc79692-8db4-44a2-9433-5b6ce97b656f
# ╠═7da59e27-62b9-4b89-b315-d88a4fd34f56
# ╠═40906795-a4dd-4e4a-a62e-91b4639a48fa
# ╠═edd9fdd1-a9c4-4f45-8f63-9681717d417f
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
# ╟─2790b550-3105-4fc0-9070-d142c19678db
# ╠═bea97fb3-9854-411c-8363-15cbef13d033
# ╟─31add356-6854-43c5-9ebd-ef10add6cc3d
# ╠═60b475d1-219e-4037-bd36-e0671cfa1893
# ╠═bd7552d2-2c31-4834-97d9-ccdb4652242f
# ╟─e81e803a-d831-4d62-939c-1e4a4fdec74f
# ╟─c4521a0c-c5af-43cd-97bc-a4a7a42d27b1
# ╟─b34c1d1b-a5b8-4de9-bea9-3f2d0503c1c0
# ╠═0c3eb801-f271-4cf4-a109-25a283a51779
# ╟─8b1a0902-2542-40ed-9f91-447bffa4290f
# ╠═f39dd714-972c-4d29-bfa8-d2c3795d2eef
# ╟─b06a7955-6c91-444f-9bf3-72cfb4a011ec
# ╠═6ae6d894-4923-4408-9b77-1067ba9e2aff
# ╠═7ab38bc2-9ca4-4206-a4c3-5fed673557f1
# ╟─84093807-8ea9-4290-8bf7-bd94c5c28691
# ╠═ed0c6737-1237-40d0-81df-32d6a842cbb0
# ╠═7f2276d5-bb78-487c-868a-0486595d0113
# ╟─a6e61592-7958-4094-8614-e77446eb2223
# ╟─2739bcff-3fb0-4169-8a1a-2b0a14998cec
# ╠═30393c90-298c-412d-86ce-e36106613d35
# ╟─9952c815-5459-44ff-b1f8-07ab24ce0c53
# ╟─b22387c4-cda5-424d-99b7-d7deda24c678
# ╠═15604034-91fd-4fd4-b09e-e3c5cfe7a265
# ╠═808aed68-7077-4079-be75-1bea962c716d
# ╠═5de9edf7-6059-47d5-b917-e7491068ebdc
# ╠═a55c5ee7-2274-4447-b0b2-58052f064bc9
# ╠═68e2628a-056a-4ec3-827f-2654f49917d9
# ╠═fec9ca6d-d815-4b50-bec7-f8fb3d8195ba
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
