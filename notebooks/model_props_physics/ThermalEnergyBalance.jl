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
	using ExtendableGrids, GridVisualize,ExtendableSparse
	using NLsolve, LinearSolve,Pardiso, ILUZero
	using StaticArrays

	using LessUnitful
	using DataFrames
	using PlutoVista, Plots
	using PlutoUI, Colors
	using MultiComponentReactiveMixtureProject
	using Printf
	
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

# ╔═╡ a39d4d01-dde9-4a59-85bb-6f590b2bd397
md"""

Total energy equation eq. (2.2.3) following Gionvangigli [1]:
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
let
	T = 200
md"""
## Estimation of orders of magnitude
1) Internal energy / enthalpy of N2 gas at $(T) °C inlet temperature after pre-heating: $(Integer(round(enthalpy_gas(N2, T+273.15)))) J/mol

1) Feed velocity of gas entering
"""
end

# ╔═╡ 18756563-680e-4906-8274-f8601c7e6cd7
let
	enthalpy_gas(N2, 200+273.15)
end

# ╔═╡ 826c85a3-9c81-4e7b-b06a-2310f065ef58
md"""
# References
1) Giovangigli, Vincent (1999): Multicomponent flow modeling. Boston, Basel, Berlin: Birkhäuser (Modeling and simulation in science, engineering, and technology).
"""

# ╔═╡ ef1dfaa4-e070-403f-8854-35072778d655
data=ReactorData()

# ╔═╡ Cell order:
# ╠═5d07bf50-a4ec-11ee-2769-19a4be2d9386
# ╟─c426878f-67e1-4d2d-924e-73bbba9a06c7
# ╟─48eb2a52-fabb-4daa-b6da-f8ac9a55668d
# ╠═a39d4d01-dde9-4a59-85bb-6f590b2bd397
# ╠═586fabe6-2f89-4dac-b95c-5149f86afc6f
# ╠═3d7ca3da-7446-413c-b905-6d818a778a90
# ╠═18756563-680e-4906-8274-f8601c7e6cd7
# ╟─826c85a3-9c81-4e7b-b06a-2310f065ef58
# ╠═ef1dfaa4-e070-403f-8854-35072778d655
