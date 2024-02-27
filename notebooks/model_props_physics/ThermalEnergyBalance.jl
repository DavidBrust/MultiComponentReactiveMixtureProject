### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 5d07bf50-a4ec-11ee-2769-19a4be2d9386
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,"../.."))
	using Revise
	using VoronoiFVM
	using ExtendableGrids, GridVisualize
	using LinearAlgebra
	using LessUnitful
	using PlutoVista
	using PlutoUI
	using MultiComponentReactiveMixtureProject
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ c426878f-67e1-4d2d-924e-73bbba9a06c7
md"""
# TODO
1) List total energy equation 2.2.3 (p.22) $\textcolor{blue}{OK}$
1) write eq 2.2.14 (p.24) $\textcolor{blue}{OK}$
1) for representative case, calculate specific (volumetric) internal, kinetic and gravitational/potential energies, use to justify neglect of kinetic and potential contributions
1) write internal energy equation 2.2.1 after subtraction of kinetic energy eq
1) transform into and list enthalpy equation 2.3.16, using $\rho h=\rho e + p$, this is justified based on availability of thermo data
1) write definition of $\vec Q$ in equation 2.5.3
1) list simplified initial version of $\vec Q$, neglecting cross-effects (Soret and Dufour)
1) calculate representative values for $\vec \Pi:\nabla\vec v$ and $\vec v \cdot\nabla p$ terms to justify their neglect
1) list resulting enthalpy equation that corresponds to the implemented equation
"""

# ╔═╡ 48eb2a52-fabb-4daa-b6da-f8ac9a55668d
md"""
# Energy Conservation and Simplification
"""

# ╔═╡ d906e8af-f315-4355-8b9f-15c3fc7cd744
md"""
Internal energy equation (2.2.21, Giovangigli [1], 3.76, [2]):
```math
\begin{align}
\partial_t(\rho e) + \nabla \cdot (\rho e \vec v) + \nabla \cdot (p \vec v) + \nabla \cdot \vec Q = - (\vec \Pi \:\: \colon \: \nabla \vec v) + \vec v \cdot \nabla p
\end{align}
```
"""

# ╔═╡ a39d4d01-dde9-4a59-85bb-6f590b2bd397
md"""

Total energy equation eq. (2.2.3) following Giovangigli [1]:
```math
\begin{align}
\partial_t(\rho e^{\text{tot}}) + \nabla \cdot ((\rho e^{\text{tot}} + p)\vec v) + \nabla \cdot(\vec Q+ \vec \Pi \cdot \vec v) = \sum_{i=1 \dots N} (\rho_i \vec v + F_i) \cdot \vec b_i
\end{align}
```
where $e^{\text{tot}}$ is the specific total energy and $\vec Q$ is the heat flux.
"""


# ╔═╡ 586fabe6-2f89-4dac-b95c-5149f86afc6f
md"""
Decompose total specific energy of the mixture $e^{\text{tot}}$ into

```math
	e^{\text{tot}} = e + \frac{1}{2}\vec v \cdot \vec v
```
where $e$ is the specific internal energy and $\frac{1}{2}\vec v \cdot \vec v$ is the specific kinetic energy.
"""

# ╔═╡ 3d7ca3da-7446-413c-b905-6d818a778a90
md"""
## Estimation of orders of magnitude
"""

# ╔═╡ 8a40b0bf-a7dd-4200-9ff5-01281272fe23
md"""
### Specific internal energy
```math
\begin{align}
e_i(T) &= e_i^0 + \int_{T_{\text{ref}}}^T c_{v,i}(\widetilde{T}) d\widetilde{T} \\
&\approx e_i^0 + \frac{c_{v,i}(T) + c_{v,i}(T_{\text{ref}})}{2} (T-T_{\text{ref}}) \\
c_{v,i} &= c_{p,i} - \frac{R}{M_i}
\end{align}
```
"""

# ╔═╡ e8749b9f-c394-43aa-9ae1-fa00f49f3b64
md"""
### Specific kinetic energy
```math
k = \frac{1}{2}\vec v \cdot \vec v
```
"""

# ╔═╡ 520848fc-1d79-4cc6-8448-9a5a6bb755bf
function setup_solve(;dim=2, nref=0)
	
	grid, inb, irrb, outb, sb, catr =  PTR_grid_boundaries_regions(dim, nref=nref)
	data = ReactorData(
		dim=dim,
		solve_T_equation=false,
		Tamb = 200 + 273.15,
		X0 = [0,0,0,0,0.0,1.0], # N2
		inlet_boundaries=inb,
		irradiated_boundaries=irrb,
		outlet_boundaries=outb,
		side_boundaries=sb,
		catalyst_regions=catr,
	)
	
	inival,sys = PTR_init_system(dim, grid, data)
	
	control = SolverControl(nothing, sys;)
	control.handle_exceptions=true
	control.Δu_opt=100
	
	times = [0,50]
	solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="a",log=true)
	
	return solt,grid,sys,data
end

# ╔═╡ 82be58bb-8bb3-4769-b7b9-1fbfc5f9d5ea
solt,grid,sys,data = setup_solve()

# ╔═╡ 18756563-680e-4906-8274-f8601c7e6cd7
let
	T = data.Tamb
	Tref = N2.Tref
	cv_T = (heatcap_gas(N2, T)-ph"R")/N2.MW # J / (kg * K)
	cv_Tref = (heatcap_gas(N2, Tref)-ph"R")/N2.MW # J / (kg * K)

	internal_energy = 0.5*(cv_T+cv_Tref) * (T-Tref) # J / (kg)
end

# ╔═╡ b7a0e5f5-3772-419c-8e20-10660275b3aa
sol = solt(solt.t[end]);

# ╔═╡ f53131ee-fafe-4c72-a25d-50c58cbd2484
let
	(;ip,X0,Tamb) = data
	nf = nodeflux(sys, sol)
	massflux = nf[:,ip,:]

	mmix = molarweight_mix(X0,data)
	rho = @. sol[ip,:] * mmix /(ph"R"*Tamb)

	v = massflux./ repeat(reshape(rho,1,:),2)
	v_magnitude = [sqrt(dot(v_,v_)) for v_ in eachcol(v)]
	#@info maximum(v_magnitude)
	k = 0.5* [dot(v_,v_) for v_ in eachcol(v)]
	maximum(k)
end

# ╔═╡ 2af672fe-17b0-4a2d-8209-6caef0d94a5c
md"""
### Flow field
1. Pressure
2. Density
3. Velocity X
4. Velocity Y
"""

# ╔═╡ 3bcaf0df-dbbc-44ff-87c0-8eee495e0ea8
let
	(;p,m,ip,iT,Tamb,mfluxin,solve_T_equation,dim) = data
	ng = ngas(data)
	mmix = []
	for j in 1:length(sol[1,:])
		_mmix=0
		for i=1:ng
			_mmix += sol[i,j]*m[i]
		end
		push!(mmix, _mmix)
	end
	
	ps = sol[ip,:]
	if solve_T_equation
		Ts = sol[iT,:]
	else
		Ts = Tamb
	end
	#rho = @. ps * mmix /(ph"R"*T)
	rho = @. ps * mmix /(ph"R"*Ts)
	
	
	if dim == 2
		vis=GridVisualizer(layout=(4,1), resolution=(600,800))
		scalarplot!(vis[1,1], grid, ps, aspect=4.0, zoom=3.5) # Total pressure
		scalarplot!(vis[2,1], grid, rho, aspect=4.0, zoom=3.5) # Total Density
		nf = nodeflux(sys, sol)
		massflux = nf[:,ip,:]
		scalarplot!(vis[3,1], grid, massflux[1,:]./rho, aspect=4.0, zoom=3.5) # Velocity - X
		scalarplot!(vis[4,1], grid, massflux[2,:]./rho, aspect=4.0, zoom=3.5) # Velocity - Y
	else
		vis=GridVisualizer(layout=(2,1), resolution=(400,800), outlinealpha=0.0)
		scalarplot!(vis[1,1], grid, ps, title="Total Pressure")
		scalarplot!(vis[2,1], grid, rho, title="Total Density")
		
	end
	reveal(vis)
end

# ╔═╡ 826c85a3-9c81-4e7b-b06a-2310f065ef58
md"""
# References
1) Giovangigli, Vincent (1999): Multicomponent flow modeling. Boston, Basel, Berlin: Birkhäuser (Modeling and simulation in science, engineering, and technology).
1) Marquardt, Wolfgang (2014): Modellierung technischer Systeme. Vorlesungsmanuskript, RWTH Aachen University.
"""

# ╔═╡ Cell order:
# ╠═5d07bf50-a4ec-11ee-2769-19a4be2d9386
# ╟─c426878f-67e1-4d2d-924e-73bbba9a06c7
# ╟─48eb2a52-fabb-4daa-b6da-f8ac9a55668d
# ╟─d906e8af-f315-4355-8b9f-15c3fc7cd744
# ╟─a39d4d01-dde9-4a59-85bb-6f590b2bd397
# ╟─586fabe6-2f89-4dac-b95c-5149f86afc6f
# ╟─3d7ca3da-7446-413c-b905-6d818a778a90
# ╟─8a40b0bf-a7dd-4200-9ff5-01281272fe23
# ╠═18756563-680e-4906-8274-f8601c7e6cd7
# ╟─e8749b9f-c394-43aa-9ae1-fa00f49f3b64
# ╠═f53131ee-fafe-4c72-a25d-50c58cbd2484
# ╠═520848fc-1d79-4cc6-8448-9a5a6bb755bf
# ╠═82be58bb-8bb3-4769-b7b9-1fbfc5f9d5ea
# ╠═b7a0e5f5-3772-419c-8e20-10660275b3aa
# ╟─2af672fe-17b0-4a2d-8209-6caef0d94a5c
# ╟─3bcaf0df-dbbc-44ff-87c0-8eee495e0ea8
# ╟─826c85a3-9c81-4e7b-b06a-2310f065ef58
