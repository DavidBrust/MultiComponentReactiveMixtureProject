### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 64469510-f2f5-11ed-0dd5-2b60d8d52b40
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using PlutoUI
	using MultiComponentReactiveMixtureProject
end;

# ╔═╡ 7011d5ff-58f2-4ba9-a64c-96fb4df689f4
TableOfContents(aside=false)

# ╔═╡ ceb22984-3af0-4d76-8e27-b5cba9c4e51c
md"""
# The model
Model equations for multi-component transport of reacting (ideal) gas phase mixtures in porous media.
## Overall Mass Continuity
Mixture mass flow (overall mass flow) through the pore space of the porous medium.
Mixture mass averaged velocity is calculated from Darcy equation. The void fraction (porosity) is given by $\epsilon$.

```math
\begin{align}
	\frac{\partial \epsilon \rho}{\partial t} + \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\end{align}
```

## Species Mass Continuity and Transport
```math
\begin{align}
	\frac{\partial \epsilon \rho_i}{\partial t} + \nabla \cdot \left( \vec \Phi_i + \rho_i \vec v \right ) - r_i(\varrho) &= 0 ~, \qquad i = 1 ... n \\
		\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= -\sum_{j=1 \atop j \neq i}^{n} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
		\sum_{i=1}^n x_i &= 1
\end{align}
```

where $\rho$ is the (total) mixture density, $\vec v$ is the mass-averaged (barycentric)  mixture velocity calculated with the Darcy equation,
$x_i$, $w_i$ and $M_i$ are the molar fraction, mass fraction and molar mass of species $i$ respectively,
$\vec \Phi_i$ is the __diffusive__ mass flux of species $i$ ($\frac{\text{kg}}{\text{m}^2 \text{s}}$)
and $r_i$ is the species mass volumetric source/sink ($\frac{\text{kg}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
The __convective__ mass flux of species $i$ is the product of the superficial mean velocity with the partial mass density $\rho_i \vec v$.
The combined __convective-diffusive__ species mass flux $\vec \Psi_i = \vec \Phi_i + \rho_i \vec v$ can be introduced as an auxiliary variable.


"""

# ╔═╡ 00663964-7c47-4ba8-9fc9-e65c63c9f6b8
md"""
## Thermal Energy Transport

Enthalpy equation for __gas phase only__.
Considers convective-diffusive transport of thermal energy,
including "thermal drift" from convective-diffusive species fluxes $\vec \Phi_i + \rho_i \vec v$ carrying enthalpy.

Formulation based on separation of "reference  enthalpy" and "thermal enthalpy".  See section "Derivation of Separated Formulation".
```math
\begin{align}
\frac{\partial (\sum \rho_i h_i^{\text{th}}(T) )}{\partial t} + \nabla \cdot \left( \sum h_i^{\text{th}}(T) \left( \rho_i \vec v + \vec \Phi_i \right) + \vec q \right ) + \sum h_i^0 r_i &= 0
\end{align}
```

```math
\begin{align}
\frac{\partial \overline h}{\partial t} - \nabla \cdot \left(\lambda_{\text{eff}} \nabla T - \sum_i^n \vec \Psi_i h_i / M_i \right) &= 0 \\
\overline h &= (\varepsilon \sum_i^n \rho_i h_i / M_i + [1-\varepsilon] \rho_{\text s} h_{\text s}) \\
h_i(T) &= h_i^0(T_{\text{ref}}) + \int_{T_{\text{ref}}}^{T} c_{p,i} \,dT 
\end{align}
```
where $\overline h$ is the volumetric enthalpy [$\text{J/m}^3$] of the quasi-homogeneous domain comprised of gas and solid volume fractions, the gas phase species enthalpy $h_i$ and $\rho_{\text s}$ and $h_{\text s}$ are the mass density and enthalpy of the solid matrix respectively.

Here $\lambda_{\text{eff}}$ is the effective thermal conductivity in the quasi-homogeneous approach, which is a (complicated) function of the thermal conductivities of the solid and gas phases.

Another process contributing to heat transport within the porous medium is "thermal drift", resulting from gas phase mass flow carrying enthalpy.
It is expressed as the sum  of species mass fluxes and species specific enthalpies $\sum_i^n \vec \Psi_i h_i / M_i$ (division by $M_i$ as $h_i$ is given on a molar basis).
The gas phase species enthalpies consist of a chemical contribution originating from the potential energy stored within the chemical bonds of the species, and a thermal contribution depending on the species' heat capacity and temperature.
During chemical reactions the chemical bonds in the species involved in the reaction are rearranged. If in the process of rearranging of chemical bonds the products have a different chemical potential energy than the reactants, the difference is made up from thermal energy.
In case of an exothermal reaction, the products have less chemical energy thus the difference is released as heat. In contrast for endothermal reactions the products hava more chemical energy than the reactants and thus consume heat during the reaction.
"""

# ╔═╡ Cell order:
# ╠═64469510-f2f5-11ed-0dd5-2b60d8d52b40
# ╠═7011d5ff-58f2-4ba9-a64c-96fb4df689f4
# ╟─ceb22984-3af0-4d76-8e27-b5cba9c4e51c
# ╟─4ac838f9-777e-41d6-a89c-5ed4282b4288
# ╟─37ec2c8a-1711-4e1e-badd-0f5322eef41e
# ╟─52afc3f1-064c-4a45-af0d-942e89e2c524
# ╟─f2668597-d7c1-4200-ad6b-bc6d536068ef
# ╟─83cd18d8-40c8-49c8-8b01-32f0aa8cc17b
# ╟─30ca1e60-ee79-4b6d-9c66-d5e3313ade3a
# ╟─7c3f89da-9b81-4395-8ffc-00c65fdc7529
# ╟─cf75ad38-c5d6-4258-ab39-710ece3b7663
# ╟─70fbdc66-fefc-432e-abeb-05a322b34e00
# ╟─cb64283a-4c74-43ea-be69-a0ee293491fd
# ╟─730348aa-3ff7-43bd-8c43-e9ba6823d318
# ╟─3e53e30d-f601-4d78-9f93-037172e40504
