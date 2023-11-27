### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ c21e1942-628c-11ee-2434-fd4adbdd2b93
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using VoronoiFVM
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve, LinearSolve
	using StaticArrays

	using LessUnitful
	using DataFrames
	using PlutoVista, Plots
	using PlutoUI, Colors
	using FixedBed
	using Printf
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 6da83dc0-3b0c-4737-833c-6ee91552ff5c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=false))
"""
  ╠═╡ =#

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═╡ skip_as_script = true
#=╠═╡
PlutoUI.TableOfContents(title="M-S Transport + Darcy")
  ╠═╡ =#

# ╔═╡ 83fa22fa-451d-4c30-a4b7-834974245996
function grid1D()
	X=(0:0.02:1)*ufac"cm"
	grid=simplexgrid(X)
	# catalyst region
	cellmask!(grid,[0.4]*ufac"cm",[0.6]*ufac"cm",2)	
	#cellmask!(grid,[0.0]*ufac"cm",[0.2]*ufac"cm",2)	
	grid
end

# ╔═╡ 4dae4173-0363-40bc-a9ca-ce5b4d5224cd
function grid2D()
	R=(0:1:7)*ufac"cm"
	Z=(0:0.05:0.5)*ufac"cm"
	grid=simplexgrid(R,Z)
	circular_symmetric!(grid)

	cellmask!(grid,[0.0,0.45].*ufac"cm",[7.0,0.5].*ufac"cm",2) # catalyst region	
	#cellmask!(grid,[0.0,0.4].*ufac"cm",[5.0,0.5].*ufac"cm",2) # catalyst region
	bfacemask!(grid, [0.0,0.5].*ufac"cm",[5.0,0.5].*ufac"cm",5) # top inner
	grid
end

# ╔═╡ 561e96e2-2d48-4eb6-bb9d-ae167a622aeb
function grid3D()
	X=(0:1:14)*ufac"cm"
	Y=(0:0.05:0.5)*ufac"cm"
	Z=(0:1:14)*ufac"cm"
	grid=simplexgrid(X,Y,Z)

	# catalyst region
	cellmask!(grid,[2,0.0,2].*ufac"cm",[12,0.1,12].*ufac"cm",2)
	#bfacemask!(grid, [2,1,2].*ufac"cm",[12,1,12].*ufac"cm",7) # mask outflow
	bfacemask!(grid, [2,0,2].*ufac"cm",[12,0,12].*ufac"cm",7) # mask inflow
	
	grid
end

# ╔═╡ 107a6fa3-60cb-43f0-8b21-50cd1eb5065a
const dim = 2

# ╔═╡ 4e05ab31-7729-4a4b-9c14-145118477715
# ╠═╡ skip_as_script = true
#=╠═╡
if dim == 3
	@bind xcut Slider(linspace(0,14,21)*ufac"cm",show_value=true,default=6.5*ufac"cm")
end
  ╠═╡ =#

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═╡ skip_as_script = true
#=╠═╡
let
	if dim == 1
		gridplot(grid1D(), resolution=(600,200))
	elseif dim == 2
		gridplot(grid2D(), resolution=(660,300), aspect=4.0, zoom=2.8)
	else
		gridplot(grid3D(); xplane=xcut, show=true, outlinealpha=0.0 )
	end
end
  ╠═╡ =#

# ╔═╡ 832f3c15-b75a-4afe-8cc5-75ff3b4704d6
begin
	if dim == 1
		const Γ_left = 1
		const Γ_right = 2
	elseif dim == 2
		const Γ_bottom = 1
		const Γ_outer = 2
		#const Γ_top_mask = 3
		const Γ_top_outer = 3
		const Γ_sym = 4		
		const Γ_top_inner = 5
	else
		# mask inflow
		const Γ_left = 7 
		const Γ_right = 3
		# mask outflow
		# const Γ_left = 1 
		# const Γ_right = 7		
		const Γ_front = 2
		const Γ_back = 4
		const Γ_bottom = 5
		const Γ_top = 6
	end
end;

# ╔═╡ a078e1e1-c9cd-4d34-86d9-df4a052b6b96
md"""
Reactor Simulation of PC reactor for cylindrical symmetrical geometry (2D).
"""

# ╔═╡ 0fadb9d2-1ccf-4d44-b748-b76d911784ca
md"""
## Overall Mass Continuity
Mixture mass flow (overall mass flow) through the pore space of the porous medium. Mixture mass averaged velocity is calculated from Darcy equation. The void fraction (porosity) is given by $\epsilon$.

```math
\begin{align}
	\frac{\partial \epsilon \rho}{\partial t} + \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\end{align}
```
"""

# ╔═╡ b94513c2-c94e-4bcb-9342-47ea48fbfd14
md"""
## Species Mass Continuity and Transport
```math
\begin{align}
	\frac{\partial \epsilon \rho_i}{\partial t} + \nabla \cdot \left( \vec \Phi_i + \rho_i \vec v \right ) - R_i &= 0 ~, \qquad i = 1 ... \nu \\
		\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= -\sum_{j=1 \atop j \neq i}^{\nu} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
		\sum_{i=1}^\nu x_i &= 1
\end{align}
```
"""

# ╔═╡ c886dd12-a90c-40ab-b9d0-32934c17baee
md"""
where $\rho$ is the (total) mixture density, $\vec v$ is the mass-averaged (barycentric)  mixture velocity calculated with the Darcy equation, $x_i$, $w_i$ and $M_i$ are the molar fraction, mass fraction and molar mass of species $i$ respectively, $\vec \Phi_i$ is the mass flux of species $i$ ($\frac{\text{kg}}{\text{m}^2 \text{s}}$) and $R_i$ is the species mass volumetric source/sink ($\frac{\text{kg}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
"""

# ╔═╡ 8f2549f4-b0a6-440f-af94-6880e0814dc2
md"""
## Thermal Energy Transport
"""

# ╔═╡ 78589a1e-2507-4279-ba42-1aaec90d87d0
md"""
```math
\begin{align}
\frac{\partial((1-\epsilon)\rho_{\text s} c_{\text s}+\epsilon c c_p) T}{\partial t} - \nabla \cdot \left(\lambda_{\text{eff}} \nabla T - \sum_i^{\nu} \vec N_i h_i \right)  &= 0
\end{align}
```

where $\lambda_{\text{eff}}$ is the effective thermal conductivity. The heat release from chemical reactions is considered as part of the species enthalpies. The convective heat transport within the porous medium is expressed via the sum over the product of species (mass) fluxes with species mass-specific enthalpies $\vec H_i= \vec N_i h_i$.

As part of the initialisation strategy (see next section) the convective contribution of heat flux is ramped up at the same time with the heat transport boundary conditions after the flow field has been established.
"""

# ╔═╡ 3bb2deff-7816-4749-9f1e-c1e451372b1e
function reaction(f,u,node,data)
	(;m,ip,iT,isreactive)=data
	ng=ngas(data)

	if node.region == 2 && isreactive # catalyst layer
		(;lcat,kinpar,gni)=data
		(;rni,nuij)=kinpar
		
		pi = MVector{ng,eltype(u)}(undef)
		for i=1:ng
            pi[i] = u[ip]*u[i]
		end

        RR = @inline -lcat*ri(data,u[iT],pi)
		
		for i=1:ng
			f[i] = zero(eltype(u))
			for j=1:nreac(kinpar)
				f[i] += nuij[(j-1)*ng+i] * RR[j] * m[i]
			end			
		end
	end
	
	for i=1:ng
		f[ng] += u[i]
	end
	f[ng] = f[ng] - 1.0
end

# ╔═╡ 4af1792c-572e-465c-84bf-b67dd6a7bc93
function storage(f,u,node,data)
	(;ip,iT,m,poros,rhos,cs)=data
	ng=ngas(data)

	c = u[ip]/(ph"R"*u[iT])
	for i=1:ng
		f[i]=c*u[i]*m[i]*poros
	end
	
	# total pressure
	@inline mmix = molarweight_mix(u,data)
	X=MVector{ng,eltype(u)}(undef)
	@inline MoleFrac!(X,u,data)
	@inline cpmix = heatcap_mix(data.Fluids, u[iT], X)
	#cpmix = 0.0
	#for i=1:ng
	#	@inline cpmix += heatcap_gas(data.Fluids[i], u[iT])*u[i]
	#end
	
	f[ip] = mmix*c*poros

	# solid heat capacity is 4 orders of magnitude larger than gas phase heat cap
	#f[iT] = u[iT] * (rhos*cs*(1-poros) + cpmix*c*poros)
	f[iT] = u[iT] * (rhos*cs*(1-poros) + cpmix*c*poros) / 200
	
end

# ╔═╡ f28be8bd-4ccc-473c-b01f-730f2483ac78
md"""
## Boundary Conditions
"""

# ╔═╡ abc28f81-903b-4656-b137-881060ae459c
function radiosity_window(f,u,bnode,data)
    (;iT,iTw,G_lamp,uc_window,uc_cat,uc_mask)=data
    # irrad. exchange between quartz window (1), cat surface (2), masked sruface (3) 
    tau1_vis=uc_window.tau_vis
    rho1_vis=uc_window.rho_vis
    tau1_IR=uc_window.tau_IR
    rho1_IR=uc_window.rho_IR
    eps1=uc_window.eps

    Tglass = u[iTw] # local tempererature of quartz window
    G1_bot_IR = eps1*ph"σ"*Tglass^4
	G1_bot_vis = 0.0
    if bnode.region==Γ_top_inner # catalyst layer (2)
		# flux profile measured behind quarz in plane of cat layer
		G1_bot_vis += G_lamp

    end
    return G1_bot_vis,G1_bot_IR
end

# ╔═╡ a90c6d87-b442-43e8-ab27-827fab25d4f6
function top(f,u,bnode,data)
	(;p,ip,iT,iTw,Tamb,mfluxin,X0,W0,m,mmix0,uc_window,uc_cat,uc_mask,uc_h,Nu,k_nat_conv,dt_hf_enth,dt_hf_irrad)=data
	ng=ngas(data)

	# irrad. exchange between quartz window (1), cat surface (2), masked surface (3)
	alpha1_vis=uc_window.alpha_vis
	alpha1_IR=uc_window.alpha_IR
	
	rho2_vis=uc_cat.rho_vis
	rho2_IR=uc_cat.rho_IR
	alpha2_vis=uc_cat.alpha_vis
	alpha2_IR=uc_cat.alpha_IR
	eps2=uc_cat.eps

	#rho3_vis=uc_mask.rho_vis
	#rho3_IR=uc_mask.rho_IR
	#alpha3_vis=uc_mask.alpha_vis
	#alpha3_IR=uc_mask.alpha_IR
	#eps3=uc_mask.eps

	G1_vis, G1_IR = radiosity_window(f,u,bnode,data)
	
	# top boundaries (inlet and mask cross-sections)
	if bnode.region==Γ_top_inner|| bnode.region==Γ_top_outer

		r_mfluxin = mfluxin*ramp(bnode.time; du=(0.1,1), dt=(0.0,1.0))
		
		f[ip] = -r_mfluxin # total mass flux
		for i=1:(ng-1)
			f[i] = -r_mfluxin*W0[i] # species mass flux
		end		
		# heatflux from enthalpy inflow
		@inline r_hf_enth = mfluxin/mmix0 * enthalpy_mix(data.Fluids, Tamb+100, X0) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_enth)
		
		# heatflux from irradiation
		hflux_irrad = (-eps2*ph"σ"*u[iT]^4 + alpha2_vis*G1_vis + alpha2_IR*G1_IR) 
		
		# heatflux from convection through top chamber
        Tm=0.5*(u[iT] + u[iTw]) # mean temperature
        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X0)
		dh=2*uc_h
		kconv=Nu*λf/dh*ufac"W/(m^2*K)"
		hflux_conv = kconv*(u[iT]-Tm)
		
		f[iT] = -r_hf_enth +(-hflux_irrad + hflux_conv) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)
		
		# calculate local window temperature from (local) flux balance
		hflux_conv_top_w = k_nat_conv*(u[iTw]-Tamb)
		G2_vis = rho2_vis*G1_vis		
		G2_IR = eps2*ph"σ"*u[iT]^4 + rho2_IR*G1_IR
		hflux_abs_w = alpha1_vis*G2_vis + alpha1_IR*G2_IR + alpha1_IR*ph"σ"*Tamb^4
		hflux_emit_w = uc_window.eps*ph"σ"*u[iTw]^4
		f[iTw] = -hflux_conv -hflux_abs_w +hflux_conv_top_w +2*hflux_emit_w
		
	# elseif bnode.region==Γ_top_mask
	#
	#	flux_irrad = -eps3*ph"σ"*u[iT]^4 + alpha3_vis*G1_vis + alpha3_IR*G1_IR
	#	G3_vis = rho3_vis*G1_vis		
	#	G3_IR = eps3*ph"σ"*u[iT]^4 + rho3_IR*G1_IR
	#	flux_abs_w = alpha1_vis*G3_vis + alpha1_IR*G3_IR
	#
	#	f[iT] = -flux_irrad
	end
	
end

# ╔═╡ 7f8955d2-5f0b-4217-b08d-e3f53849dbec
function side(f,u,bnode,data)
	# side wall boundary condition
	(;iT,k_nat_conv,delta_gap,Tamb)=data
	ng=ngas(data)
	
	# all sides for complete domain
	if bnode.region==Γ_outer

		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative
		X=MVector{ng,eltype(u)}(undef)
		@inline MoleFrac!(X,u,data)
		@inline _,λf=dynvisc_thermcond_mix(data, u[iT], X)

		# w/o shell height
		f[iT] = (u[iT]-Tamb)/(delta_gap/λf+1/k_nat_conv)
		
	end
end

# ╔═╡ 389a4798-a9ee-4e9c-8b44-a06201b4c457
function boutflow(f,u,edge,data)
	(;iT,iTp,ip,m,lc_frit,lc_plate,dt_hf_enth)=data
	ng=ngas(data)

	k=outflownode(edge)
	pout = u[ip,k]
	cout = pout/(ph"R"*u[iT,k])
	X = MVector{ng,eltype(u)}(undef)
	
	for i=1:ng
		X[i] = u[i,k]
	end
	@inline mumix, _ = dynvisc_thermcond_mix(data, u[iT,k], X)
	v = DarcyVelo(u,data,mumix)
	
	for i=1:(ng-1)
		f[i] = v*cout*u[i,k]*m[i] # species mass flux at outflow
	end
	# convective heat flux
	@inline r_hf_conv = v *cout * enthalpy_mix(data.Fluids, u[iT,k], X) * ramp(edge.time; du=(0.0,1.0), dt=dt_hf_enth)

	f[iT] = r_hf_conv
end

# ╔═╡ ee797850-f5b5-4178-be07-28192c04252d
function bflux(f,u,bedge,data)
	# window temperature distribution
	if bedge.region == Γ_top_inner || bedge.region == Γ_top_outer
	#if bedge.region == Γ_top_inner
		(;iTw,lambda_window) = data
		f[iTw] = lambda_window * (u[iTw, 1] - u[iTw, 2])
	elseif bedge.region == Γ_bottom
		(;iTp,lambda_Al) = data
		f[iTp] = lambda_Al * (u[iTp, 1] - u[iTp, 2])
	end
end

# ╔═╡ b375a26d-3ee6-4b69-9bd8-c69c3e193dc9
function bstorage(f,u,bnode,data)
	if bnode.region == Γ_top_inner|| bnode.region == Γ_top_outer
		(;iTw) = data
		f[iTw] = u[iTw]
	elseif bnode.region == Γ_bottom
		(;iTp) = data
		f[iTp] = u[iTp]
	end
end

# ╔═╡ 98468f9e-6dee-4b0b-8421-d77ac33012cc
md"""
### Temperature
1) Porous frit + catalyst layer domain
2) Window inner surface
3) Bottom plate
"""

# ╔═╡ c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
md"""
### Molar fractions
1) CO
2) CO2
"""

# ╔═╡ eb9dd385-c4be-42a2-8565-cf3cc9b2a078
md"""
### Flow field
1. Pressure
2. Density
3. Velocity X
4. Velocity Y
"""

# ╔═╡ 68ca72ae-3b24-4c09-ace1-5e340c8be3d4
function len(grid)
	coord = grid[Coordinates]
	L=0.0
	if dim == 1
		L=coord[end]
	elseif dim == 2
		L=coord[1,end]
	else
		L=coord[2,end]
	end
	L*ufac"m"
end

# ╔═╡ db77fca9-4118-4825-b023-262d4073b2dd
md"""
### Peclet Number
```math
\text{Pe}_L= \frac{L \vec v}{D}

```
"""

# ╔═╡ e7497364-75ef-4bd9-87ca-9a8c2d97064c
md"""
### Damköhler Number
```math
\text{Da}= \frac{k_{\text{react}}}{k_{\text{conv}}} = \frac{\dot r}{c_0} \tau = \frac{\dot r L}{c_0 v}

```
"""

# ╔═╡ f6e54602-b63e-4ce5-b8d5-f626fbe5ae7a
md"""
### Reynolds Number
```math
\text{Re} = \frac{\rho v D}{\eta}  
```
where $\rho$ is the density of the fluid, ideal gas in this case, $v$ is the magnitufe of the velocity, $D$ is a representative length scale, here it is the mean pore diameter, and $\eta$ is the dynamic viscosity of the fluid.
"""

# ╔═╡ 0d507c5e-eb0f-4094-a30b-a4a82fd5c302
md"""
### Knudsen Number
```math
\text{Kn} = \frac{\lambda}{l} = \frac{k_{\text B}T}{\sqrt 2 \pi \sigma^2pl}
```
where $\lambda$ is the mean free path length of the fluid, ideal gas in this case and $l$ is the pore diameter of the porous medium. Determine the mean free path length for ideal gas from kinetic gas theory and kinetic collision cross-sections (diameter).

-  $\text{Kn} < 0.01$: Continuum flow
-  $0.01 <\text{Kn} < 0.1$: Slip flow
-  $0.1 < \text{Kn} < 10$: Transitional flow
-  $\text{Kn} > 10$: Free molecular flow
"""

# ╔═╡ 67adda35-6761-4e3c-9d05-81e5908d9dd2
md"""
## Model Data
"""

# ╔═╡ a53197e9-504f-46e0-9c88-37b6f4bad0f6
# optical parameters for quartz window in upper chamber (uc)
const uc_window = SurfaceOpticalProps(	
	alpha_IR=0.3,
	tau_IR=0.7,	
	alpha_vis=0.0,
	tau_vis=0.9,
)

# ╔═╡ fb70db70-33a5-4191-a6c6-0b888aa2dd4a
#  optical parameters for catalyst layer in upper chamber (uc)
const uc_cat = SurfaceOpticalProps(
	alpha_IR=0.45, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.45, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# ╔═╡ be64e957-87f9-42bb-a423-ed6b511060ff
# optical parameters for uncoated frit in upper chamber (uc)
const uc_mask = SurfaceOpticalProps(
	alpha_IR=0.2, # measurement (ideally at high T) integrated over IR
	tau_IR=0.0, # opaque surface
	alpha_vis=0.2, # measurement (ideally at high T) integrated over vis
	tau_vis=0.0 # opaque surface	
)

# ╔═╡ 9131ec72-0c9b-49ee-baf4-e0bb4bc3a581
# optical parameters for uncoated frit in lower chamber (lc) = frit in upper chamber
const lc_frit = uc_mask

# ╔═╡ e215d9f4-4ca8-4d6f-8049-ac387eca9ce3
# optical parameters for Al bottom plate in lower chamber (lc)
const lc_plate = SurfaceOpticalProps(
	# aluminium bottom plate (machined surface)
	alpha_IR=0.1, # see in lit, assume large reflectivity
	tau_IR=0.0, # opaque surface
	alpha_vis=0.1, # see in lit, assume large reflectivity
	tau_vis=0.0 # opaque surface	
)

# ╔═╡ 85439b96-fe21-4b2c-8695-f7088305112f
function bottom(f,u,bnode,data)
	(;iT,iTp,k_nat_conv,Tamb,lc_h,Nu,dt_hf_irrad) = data
	ng=ngas(data)
	if bnode.region == Γ_bottom

		# irradiation heat flux
		rho1_IR=lc_frit.rho_IR
		alpha1_IR=lc_frit.alpha_IR
		eps1=lc_frit.eps
		rho2_IR=lc_plate.rho_IR
		alpha2_IR=lc_plate.alpha_IR
		eps2=lc_plate.eps

		# heatflux from irradiation exchange with bottom plate
		Tplate = u[iTp]	
		hflux_irrad = -eps1*ph"σ"*u[iT]^4 + alpha1_IR/(1-rho1_IR*rho2_IR)*(eps2*ph"σ"*Tplate^4+rho2_IR*eps1*ph"σ"*u[iT]^4)
	
		# heatflux from convective/conductive heat exchange with bottom plate
		Tm=0.5*(u[iT] + Tplate)

		X=MVector{ng,eltype(u)}(undef)
		@inline MoleFrac!(X,u,data)
        # thermal conductivity at Tm and outlet composition X
        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X)

        # flux_cond = -λf*(u[iT]-Tplate)/lc_h # positive flux in positive z coord.
		dh=2*lc_h
		kconv=Nu*λf/dh*ufac"W/(m^2*K)"
		hflux_conv = kconv*(u[iT]-Tm) # positive flux in negative z coord. (towards plate)

		# sign convention: outward pointing fluxes (leaving the domain) as positive, inward pointing fluxes (entering) as negative
		f[iT] = (-hflux_irrad + hflux_conv) * ramp(bnode.time; du=(0.0,1), dt=dt_hf_irrad)
		
		# calculate (local) plate temperature from (local) flux balance
		hflux_conv_bot_p = k_nat_conv*(u[iTp]-Tamb)*2
		hflux_emit_p = lc_plate.eps*ph"σ"*u[iTp]^4
		G1_IR = (eps1*ph"σ"*u[iT]^4 + rho1_IR*eps2*ph"σ"*u[iTp]^4)/(1-rho1_IR*rho2_IR)
		hflux_abs_p = alpha2_IR*G1_IR + alpha2_IR*ph"σ"*Tamb^4
		f[iTp] = -hflux_conv -hflux_abs_p +2*hflux_emit_p +hflux_conv_bot_p		
	end
end

# ╔═╡ 5f88937b-5802-4a4e-81e2-82737514b9e4
function bcond(f,u,bnode,data)
	(;p,ip,iT,Tamb,mfluxin,X0,W0,m,mmix0)=data
	ng=ngas(data)
		
	if dim==2		
		top(f,u,bnode,data)
		side(f,u,bnode,data)
		bottom(f,u,bnode,data)
		
		boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_bottom,value=p)
	else
		boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_left,value=Tamb)
		for i=1:(ng-1)
		boundary_neumann!(f,u,bnode, species=i,region=Γ_left,value=r_mfluxin*W0[i])
		end
		boundary_neumann!(f,u,bnode, species=ip, region=Γ_left, value=r_mfluxin)
		boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_right,value=p)
		boundary_dirichlet!(f,u,bnode, species=iT,region=Γ_right,value=Tamb)
	end
end

# ╔═╡ ab395450-5b1e-40cd-a184-265a2100ca35
ModelData{}()

# ╔═╡ b872acc7-c07c-4d5e-b1bc-32fea369d6a5
typeof(KinData{}())

# ╔═╡ f4dba346-61bc-420d-b419-4ea2d46be6da
function D_matrix!(data, D, T, p)
	(;m,γ_τ)=data
	ng=ngas(data)
	@inbounds for i=1:(ng-1)
		for j=(i+1):ng
			Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
			#Dji *= m[i]*m[j]
			Dji *= m[i]*m[j]*γ_τ # porosity corrected eff. diffusivity
			D[j,i] = Dji
			D[i,j] = Dji
		end
	end
end

# ╔═╡ 5547d7ad-dd58-4b00-8238-6e1abb32874e
function flux(f,u,edge,data)
	(;m,ip,iT,dt_hf_enth)=data
	ng=ngas(data)
		
	F = MVector{ng-1,eltype(u)}(undef)
	X = MVector{ng,eltype(u)}(undef)
	W = MVector{ng,eltype(u)}(undef)
	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)
	D = MMatrix{ng,ng,eltype(u)}(undef)

	pm = 0.5*(u[ip,1]+u[ip,2])
	Tm = 0.5*(u[iT,1]+u[iT,2])
	c = pm/(ph"R"*Tm)
	
	δp = u[ip,1]-u[ip,2]
	
	@inline MoleFrac!(X,u,data)
	@inline mmix = molarweight_mix(X,data)
	@inline MassFrac!(X,W,data)
	
	
	@inline D_matrix!(data, D, Tm, pm)
	@inline mumix, lambdamix = dynvisc_thermcond_mix(data, Tm, X)
	lambda_bed=kbed(data,lambdamix)*lambdamix
	
	rho = c*mmix
	v = DarcyVelo(u,data,mumix)
	
	f[ip] = -rho*v

	@inline M_matrix!(M, W, D, data)
	
	@inbounds for i=1:(ng-1)
		F[i] = ( u[i,1]-u[i,2] + (X[i]-W[i])*δp/pm )*c/mmix
	end				

	@inline inplace_linsolve!(M,F)

	hf_conv = zero(eltype(u))
	@inbounds for i=1:(ng-1)
		f[i] = -(F[i] + c*X[i]*m[i]*v)
	end
	@inline hf_conv = f[ip] * enthalpy_mix(data.Fluids, Tm, X) / mmix * ramp(edge.time; du=(0.0,1), dt=dt_hf_enth) 
	
	Bp,Bm = fbernoulli_pm(hf_conv/lambda_bed/Tm)
	f[iT] = lambda_bed*(Bm*u[iT,1]-Bp*u[iT,2])
end

# ╔═╡ 37b5908c-dd4e-4fb8-9d5b-68402493e10d
function DiffCoeffsMass(ng,m)
	k=0
	D=zeros(Float64, ng,ng)
	for i=1:(ng-1)
		for j=(i+1):ng
			k +=1
			#Dji = k*1.0e-5*ufac"m^2/s"
			Dji = k*1.0e-5*ufac"m^2/s"
			Dji *= m[i]*m[j]
			D[j,i] = Dji
			D[i,j] = Dji
		end			
	end
	D
end

# ╔═╡ 2cbd3e87-289c-47a2-b837-10133974ae82
function DiffCoeffs(data)
	(;m,γ_τ,Tamb,p)=data
	ng=ngas(data)
	D=zeros(Float64, ng,ng)
	for i=1:(ng-1)
		for j=(i+1):ng

			Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], Tamb, p)
			Dji *= γ_τ # porosity corrected eff. diffusivity
			D[j,i] = Dji
			D[i,j] = Dji
		end			
	end
	D
end

# ╔═╡ 1224970e-8a59-48a9-b0ef-76ed776ca15d
function checkinout(sys,sol)	
	tfact=TestFunctionFactory(sys)
	if dim == 2
		tf_in=testfunction(tfact,[Γ_bottom],[Γ_top_inner,Γ_top_outer])
		tf_out=testfunction(tfact,[Γ_top_inner,Γ_top_outer],[Γ_bottom])
	else
		tf_in=testfunction(tfact,[Γ_right],[Γ_left])
		tf_out=testfunction(tfact,[Γ_left],[Γ_right])
	end

	(;in=integrate(sys,tf_in,sol),out=integrate(sys,tf_out,sol) )
end

# ╔═╡ b55e2a48-5a2d-4e29-9f3b-219854461d09
function areas(sol,sys,grid,data)
	iT = data.iT
	function area(f,u,bnode,data)
		f[iT] = one(eltype(u)) # use temperature index to hold area information
	end	
	integrate(sys,area,sol; boundary=true)[iT,:]
end


# ╔═╡ 0a0f0c58-1bca-4f40-93b1-8174892cc4d8
function bareas(bfaceregion,sys,grid)
	area = 0.0
	for ibface =1:num_bfaces(grid)
        for inode=1:num_nodes(grid[BFaceGeometries][1])
			if grid[BFaceRegions][ibface] == bfaceregion
                area +=bfacenodefactors(sys)[inode,ibface]
			end
		end
	end
	area
end

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
#=╠═╡
begin
	
	if dim == 1
		mygrid=grid1D()
		strategy = nothing
		times=[0,10]
	elseif dim == 2
		mygrid=grid2D()
		strategy = nothing
		times=[0,12.0]
	else
		mygrid=grid3D()
		strategy = GMRESIteration(UMFPACKFactorization())
		times=[0,20.0]
	end
	mydata=ModelData()
	(;p,ip,Tamb,iT,iTw,iTp,X0)=mydata
	ng=ngas(mydata)
	
	sys=VoronoiFVM.System( 	mygrid;
							data=mydata,
							flux=flux,
							reaction=reaction,
							storage=storage,
							bcondition=bcond,
							bflux=bflux,
							bstorage=bstorage,
							boutflow=boutflow,
							outflowboundaries=
								[dim == 2 ? Γ_bottom : Γ_right],
							assembly=:edgewise
							)
	
	enable_species!(sys; species=collect(1:(ng+2))) # gas phase species xi, ptotal & T
	enable_boundary_species!(sys, iTw, [Γ_top_inner,Γ_top_outer]) # window temperature as boundary species in upper chamber
	enable_boundary_species!(sys, iTp, [Γ_bottom]) # plate temperature as boundary species in lower chamber
	inival=unknowns(sys)

	inival[ip,:].=p
	inival[[iT,iTw,iTp],:] .= Tamb
	
	for i=1:ng
		inival[i,:] .= X0[i]
	end

	nd_ids = unique(mygrid[CellNodes][:,mygrid[CellRegions] .== 2])
	cat_vol = sum(nodevolumes(sys)[nd_ids])
	mydata.lcat = mydata.mcat/cat_vol
    area_inner = bareas(Γ_top_inner,sys,mygrid)
    area_outer = bareas(Γ_top_outer,sys,mygrid)
	mydata.mfluxin = mydata.mflowin / (area_inner+area_outer)
	
	control = SolverControl(strategy, sys;)
		control.Δt_min=1.0e-6
		control.Δt_max=1.0
		#control.maxiters=200
		control.handle_exceptions=true
		control.Δu_opt=100.0
	function post(sol,oldsol, t, Δt)
		@info "t= "*string(round(t,sigdigits=2))*"\t Δt= "*string(round(Δt,sigdigits=2))		 
	end

	if RunSim
		solt=solve(sys;inival=inival,times,control,post)
	end
end;
  ╠═╡ =#

# ╔═╡ 927dccb1-832b-4e83-a011-0efa1b3e9ffb
#=╠═╡
md"""
## Initialisation and Solve
The simulation is setup as a transient simulation. An initialisation strategy is employed where different physics are enabled step by step once a stationary state is established. Initially, no heat is transported and no chemical reactions take place. 

1. Velocity field (mass flow is ramped up from 1-100 % in T=$(mydata.dt_mf) s)
2. Temperature field through enthalpy flux (ramped up from 0-100 % in T=$(mydata.dt_hf_enth) s)
3. Temperature field through irradiation b.c. (ramped up from 0-100 % in T=$(mydata.dt_hf_irrad) s)

The mass flow boundary condition into the reactor domain is "ramped up" starting from a low value and linearly increasing until the final value is reached. A time delay is given to let the flow stabilize. Once the flow field is established, heat transport is ramped up until a stable temperature field is established. Finally, the reactivity of the catalyst is "ramped up" until its final reactivity value is reached.
"""
  ╠═╡ =#

# ╔═╡ 9a61cf12-9d92-4fdc-9093-79d3fa1f8b90
#=╠═╡
let
	(;lcat,mcat)=mydata
	catvol = mcat/lcat
	1.0*ufac"mol/hr"/catvol
end
  ╠═╡ =#

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

# ╔═╡ b13a76c9-509d-4367-8428-7b5b316ff1ed
#=╠═╡
checkinout(sys,sol)
  ╠═╡ =#

# ╔═╡ 7c7d2f10-d8d2-447e-874b-7be365e0b00c
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;gn,gni,m,nflowin,X0) = mydata
	ng=ngas(mydata)
	in_,out_=checkinout(sys,sol)

	nout(i) = -out_[i]/m[i]
	nin(i) = nflowin*X0[i]
	nout_dry = 0.0
	
	println("Molar species in- & outflows:")
	for i = 1:ng
		@printf "%s\tIN: %2.2f\t OUT: %2.2f mol/hr\n" gn[i] nin(i)/ufac"mol/hr" nout(i)/ufac"mol/hr"
		if i != gni[:H2O] 
			nout_dry += nout(i)
		end
	end


	println("\nDry Product Molar Fractions:")
	for i=1:ng
		if i != gni[:H2O] && i != gni[:N2] 
		@printf "%3s: %2.1f%%\n" gn[i] nout(i)/nout_dry*100
		end
	end
	
	println("\nConversion:")
	for i in(gni[:CO2],gni[:H2])
		@printf "X%4s: %2.2f\n" gn[i] (nin(i)-nout(i))/nin(i)
	end
	println("\nYield & Selectivity (CO2 based):")
	for i in(gni[:CO],gni[:CH4])
		@printf "Y%4s: %2.2f \tS%4s: %2.2f\n" gn[i] nout(i)/nin(gni[:CO2]) gn[i] nout(i)/(nin(gni[:CO2])-nout(gni[:CO2]))
	end
end
  ╠═╡ =#

# ╔═╡ 5d5ac33c-f738-4f9e-bcd2-efc43b638109
#=╠═╡
let
	(;m,ip,gn,poros,mflowin,W0)=mydata
	ng=ngas(mydata)
	vis=GridVisualizer(resolution=(600,300), xlabel="Time / s", ylabel="Molar flow / Total Moles")
	
	tfact=TestFunctionFactory(sys)	
	tf_out=testfunction(tfact,[Γ_top_inner,Γ_top_outer],[Γ_bottom])
		
	inflow_rate=Float64[]
	outflow_rate=Float64[]
	reaction_rate=Float64[]
	stored_amount=Float64[]

	k=5
	for i=2:length(solt)
		m_ = k in 1:ng ? m[k] : 1
		W_ = k in 1:ng ? W0[k] : 1
		#fac = k in 1:ng ? ufac"mol/hr" : ufac"kg/hr"
		# use time dependent version of integrate function for use w testfunction
		ofr=integrate(sys,tf_out,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])

		ifr=mflowin*W_*ramp(solt.t[i]; du=(0.1,1), dt=(0.0,1.0))
		push!(inflow_rate,ifr/m_)		
		push!(outflow_rate,ofr[k]/m_)		
		rr = integrate(sys,reaction,solt[i])[k,2]
		amount = sum(integrate(sys,storage,solt[i]), dims=2)[k]
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
	name = k in 1:ng ?  gn[k] : "Total Mass"
	@printf "%s In: %2.2e \t Out: %2.2e \t React: %2.2e \nIn - Out: %2.4e \nStorage tEnd -t0: %2.4e" name I_in I_out I_reac I_in+I_out-I_reac stored_amount[end]-stored_amount[1]

	
	scalarplot!(vis, solt.t[2:end], inflow_rate, label="Inflow rate")	
	scalarplot!(vis, solt.t[2:end], -outflow_rate, label="Outflow rate", color=:red, clear=false)	
	scalarplot!(vis, solt.t[2:end], -reaction_rate, label="Reaction rate",  color=:blue, clear=false)
	scalarplot!(vis, solt.t[2:end], stored_amount, label="Stored amount", color=:green, clear=false, )
	reveal(vis)	
end

  ╠═╡ =#

# ╔═╡ 99b59260-7651-45d0-b364-4f86db9927f8
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;iT,iTw,iTp)=mydata
	vis=GridVisualizer(layout=(3,1), resolution=(680,900))
	scalarplot!(vis[1,1],mygrid, sol[iT,:] .- 273.15, zoom = 2.8, aspect=4.0)

	# plot temperature of window / boundary species
	function _2to1(a,b)
		a[1]=b[1]
		#a[2]=b[2]
	end
    # window
	bgridw = subgrid(mygrid, [Γ_top_inner]; boundary = true, transform = _2to1)
	bsolw=view(sol[iTw, :], bgridw)
	scalarplot!(vis[2,1],bgridw, bsolw.-273.15, resolution=(680,200))
	# bottom plate
	bgridp = subgrid(mygrid, [Γ_bottom]; boundary = true, transform = _2to1)
	bsolp=view(sol[iTp, :], bgridp)
	scalarplot!(vis[3,1],bgridp, bsolp.-273.15,resolution=(680,200),show=true)	
end
  ╠═╡ =#

# ╔═╡ 111b1b1f-51a5-4069-a365-a713c92b79f4
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;ip,p,gn,gni) = mydata
	ng=ngas(mydata)
	if dim == 1
		cols = distinguishable_colors(ng, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
		pcols = map(col -> (red(col), green(col), blue(col)), cols)
		vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300))
		for i=1:(ng-1)
			scalarplot!(vis, mygrid, sol[i,:], clear=false, color=pcols[i],label=gn[i])
		end
	elseif dim == 2
		vis=GridVisualizer(layout=(3,1), resolution=(680,900))
		scalarplot!(vis[1,1], mygrid, sol[gni[:CO],:], aspect = 4.0,zoom = 2.8) # CO
		scalarplot!(vis[2,1], mygrid, sol[gni[:CO2],:], aspect = 4.0,zoom = 2.8) # CO2

		cols = distinguishable_colors(ng)
		# plot species molar fractions along frit thickness (along y axis)
		function _2to1(a,b)
			a[1]=b[2]
			#a[2]=b[2]
		end
		_grid=grid2D()
		bfacemask!(_grid, [3.0,0.0].*ufac"cm",[3.0,0.5].*ufac"cm",5)
	    grid1D = subgrid(_grid, [5]; boundary = true, transform = _2to1)
		for i=1:ng
			sol1D=view(sol[i, :], grid1D)
			scalarplot!(vis[3,1],grid1D, sol1D, label=gn[i], color=cols[i],clear=false)
		end
		reveal(vis)
	
		
	else
		vis=GridVisualizer(layout=(3,1), resolution=(400,1200), outlinealpha=0.0)
		scalarplot!(vis[1,1], mygrid, sol[1,:], clear=false, label="x1")
		scalarplot!(vis[2,1], mygrid, sol[2,:], clear=false, color=:red, label="x2")
		scalarplot!(vis[3,1], mygrid, sol[3,:], clear=false, color=:blue, label="x3")
	end	
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ de69f808-2618-4add-b092-522a1d7e0bb7
# ╠═╡ skip_as_script = true
#=╠═╡
let
	mydata = ModelData()
	(;p,m,ip,iT,Tamb,mfluxin) = mydata
	ng = ngas(mydata)
	mmix = []
	for j in 1:length(sol[1,:])
		_mmix=0
		for i=1:ng
			_mmix += sol[i,j]*m[i]
		end
		push!(mmix, _mmix)
	end
	
	ps = sol[ip,:]
	Ts = sol[iT,:]
	#rho = @. ps * mmix /(ph"R"*T)
	rho = @. ps * mmix /(ph"R"*Ts)
	
	if dim == 1
		vis=GridVisualizer(legend=:lt, resolution=(600,400), layout = (2, 1))

		scalarplot!(vis[1,1], mygrid, sol[iT,:]/Tamb, clear=false, label="T / Tamb")

		p0 = sol[ip,1]
		rho0 = @. p0 * mmix[1] /(ph"R"*T)
		scalarplot!(vis[2,1], mygrid, rho/rho0, clear=false, label="Rho / Rho0")
		scalarplot!(vis[2,1], mygrid, rho0./rho, clear=false, color=:red, label="v / v0")
		scalarplot!(vis[2,1], mygrid, ps/p0, clear=false, color=:blue, label="p / p0")
	elseif dim == 2
		vis=GridVisualizer(layout=(4,1), resolution=(600,800))
		scalarplot!(vis[1,1], mygrid, ps, aspect=4.0, zoom=3.5) # Total pressure
		scalarplot!(vis[2,1], mygrid, rho, aspect=4.0, zoom=3.5) # Total Density
		nf = nodeflux(sys, sol)
		massflux = nf[:,ip,:]
		scalarplot!(vis[3,1], mygrid, massflux[1,:]./rho, aspect=4.0, zoom=3.5) # Velocity - X
		scalarplot!(vis[4,1], mygrid, massflux[2,:]./rho, aspect=4.0, zoom=3.5) # Velocity - Y
	else
		vis=GridVisualizer(layout=(2,1), resolution=(400,800), outlinealpha=0.0)
		scalarplot!(vis[1,1], mygrid, ps, title="Total Pressure")
		scalarplot!(vis[2,1], mygrid, rho, title="Total Density")
		
	end
	reveal(vis)
end
  ╠═╡ =#

# ╔═╡ ae8c7993-a89f-438a-a72a-d4a0c9a8ce57
# ╠═╡ skip_as_script = true
#=╠═╡
let
	L=maximum(mygrid[Coordinates][dim,:])
	(;mfluxin,mmix0,p,Tamb) = mydata
	ng = ngas(mydata)
	rho0 = p*mmix0/(ph"R"*Tamb)
	v0 = mfluxin / rho0
	D = DiffCoeffs(mydata)
	Pe = L*v0 / minimum(D[D.>0])
end
  ╠═╡ =#

# ╔═╡ e000c100-ee46-454e-b049-c1c29daa9a56
# ╠═╡ skip_as_script = true
#=╠═╡
let
	L=maximum(mygrid[Coordinates][dim,:])
	(;mfluxin,mmix0,lcat,p,X0,gni) = mydata
	T = 650 + 273.15
	c0 = p*X0/(ph"R"*T)
	rho0 = p*mmix0/(ph"R"*T)
	v0 = mfluxin / rho0	
	RR = -lcat*ri(mydata,T,p*X0)
	tau = L/v0
	Da = maximum(RR)/c0[gni[:H2]] * tau
end
  ╠═╡ =#

# ╔═╡ f690c19f-22e3-4428-bc1c-3ed7d1646e71
# ╠═╡ skip_as_script = true
#=╠═╡
let	
	(;mfluxin,mmix0,X0,p,Tamb,m,dp,Fluids) = mydata
	ng = ngas(mydata)
	rho0 = p*mmix0/(ph"R"*Tamb)
	v0 = mfluxin / rho0

	ηf, λf = dynvisc_thermcond_mix(mydata, Tamb, X0)
    
    cf = heatcap_mix(Fluids, Tamb, X0)
	
	Re = v0*rho0*dp/ηf # Reynolds number
	#Pr = cf*ηf/λf # Prandtl number
	#Pe = u0*ρf*cf*d/λf # Peclet
	#Re,Pr,Pe
end
  ╠═╡ =#

# ╔═╡ 5bbe72b2-2f80-4dae-9706-7ddb0b8b6dbe
# ╠═╡ skip_as_script = true
#=╠═╡
let
	(;dp,Tamb,p) = mydata
	σs = Dict(:H2 => 289*ufac"pm", :CH4 => 380*ufac"pm", :H2O => 265*ufac"pm", :N2 => 364*ufac"pm", :CO => 376*ufac"pm", :CO2 => 330*ufac"pm")

	Kn = ph"k_B"*Tamb/(sqrt(2)*pi*σs[:H2O]^2*p*dp)	
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─6da83dc0-3b0c-4737-833c-6ee91552ff5c
# ╠═d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═83fa22fa-451d-4c30-a4b7-834974245996
# ╠═4dae4173-0363-40bc-a9ca-ce5b4d5224cd
# ╠═561e96e2-2d48-4eb6-bb9d-ae167a622aeb
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═4e05ab31-7729-4a4b-9c14-145118477715
# ╠═107a6fa3-60cb-43f0-8b21-50cd1eb5065a
# ╠═832f3c15-b75a-4afe-8cc5-75ff3b4704d6
# ╟─a078e1e1-c9cd-4d34-86d9-df4a052b6b96
# ╟─0fadb9d2-1ccf-4d44-b748-b76d911784ca
# ╟─b94513c2-c94e-4bcb-9342-47ea48fbfd14
# ╟─c886dd12-a90c-40ab-b9d0-32934c17baee
# ╟─8f2549f4-b0a6-440f-af94-6880e0814dc2
# ╟─78589a1e-2507-4279-ba42-1aaec90d87d0
# ╠═5547d7ad-dd58-4b00-8238-6e1abb32874e
# ╠═3bb2deff-7816-4749-9f1e-c1e451372b1e
# ╠═4af1792c-572e-465c-84bf-b67dd6a7bc93
# ╟─f28be8bd-4ccc-473c-b01f-730f2483ac78
# ╠═abc28f81-903b-4656-b137-881060ae459c
# ╠═a90c6d87-b442-43e8-ab27-827fab25d4f6
# ╠═7f8955d2-5f0b-4217-b08d-e3f53849dbec
# ╠═85439b96-fe21-4b2c-8695-f7088305112f
# ╠═5f88937b-5802-4a4e-81e2-82737514b9e4
# ╠═389a4798-a9ee-4e9c-8b44-a06201b4c457
# ╠═ee797850-f5b5-4178-be07-28192c04252d
# ╠═b375a26d-3ee6-4b69-9bd8-c69c3e193dc9
# ╟─927dccb1-832b-4e83-a011-0efa1b3e9ffb
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═5588790a-73d4-435d-950f-515ae2de923c
# ╠═9a61cf12-9d92-4fdc-9093-79d3fa1f8b90
# ╠═b13a76c9-509d-4367-8428-7b5b316ff1ed
# ╠═7c7d2f10-d8d2-447e-874b-7be365e0b00c
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╠═5d5ac33c-f738-4f9e-bcd2-efc43b638109
# ╟─98468f9e-6dee-4b0b-8421-d77ac33012cc
# ╠═99b59260-7651-45d0-b364-4f86db9927f8
# ╟─c9c6ce0b-51f8-4f1f-9c16-1fd92ee78a12
# ╟─111b1b1f-51a5-4069-a365-a713c92b79f4
# ╟─eb9dd385-c4be-42a2-8565-cf3cc9b2a078
# ╠═de69f808-2618-4add-b092-522a1d7e0bb7
# ╟─68ca72ae-3b24-4c09-ace1-5e340c8be3d4
# ╟─db77fca9-4118-4825-b023-262d4073b2dd
# ╠═ae8c7993-a89f-438a-a72a-d4a0c9a8ce57
# ╟─e7497364-75ef-4bd9-87ca-9a8c2d97064c
# ╠═e000c100-ee46-454e-b049-c1c29daa9a56
# ╟─f6e54602-b63e-4ce5-b8d5-f626fbe5ae7a
# ╠═f690c19f-22e3-4428-bc1c-3ed7d1646e71
# ╟─0d507c5e-eb0f-4094-a30b-a4a82fd5c302
# ╠═5bbe72b2-2f80-4dae-9706-7ddb0b8b6dbe
# ╟─67adda35-6761-4e3c-9d05-81e5908d9dd2
# ╠═a53197e9-504f-46e0-9c88-37b6f4bad0f6
# ╠═fb70db70-33a5-4191-a6c6-0b888aa2dd4a
# ╠═be64e957-87f9-42bb-a423-ed6b511060ff
# ╠═9131ec72-0c9b-49ee-baf4-e0bb4bc3a581
# ╠═e215d9f4-4ca8-4d6f-8049-ac387eca9ce3
# ╠═ab395450-5b1e-40cd-a184-265a2100ca35
# ╠═b872acc7-c07c-4d5e-b1bc-32fea369d6a5
# ╠═f4dba346-61bc-420d-b419-4ea2d46be6da
# ╠═37b5908c-dd4e-4fb8-9d5b-68402493e10d
# ╠═2cbd3e87-289c-47a2-b837-10133974ae82
# ╠═1224970e-8a59-48a9-b0ef-76ed776ca15d
# ╠═b55e2a48-5a2d-4e29-9f3b-219854461d09
# ╠═0a0f0c58-1bca-4f40-93b1-8174892cc4d8
