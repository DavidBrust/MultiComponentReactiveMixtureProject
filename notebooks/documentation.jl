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
"""

# ╔═╡ 4ac838f9-777e-41d6-a89c-5ed4282b4288
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
\vec q = -\lambda \nabla T
\end{align}
```
"""

# ╔═╡ 37ec2c8a-1711-4e1e-badd-0f5322eef41e
md"""
Enthalpy equation for __gas + solid phase__. It is formulated with effective transport parameters that are a consequence of the quasi-homogeneous phase approach. It needs to be checked that this is permissible. 
Formulation based on separation of "reference  enthalpy" and "thermal enthalpy".  See also section "Derivation of Separated Formulation".

```math
\begin{align}
\frac{\partial (\varepsilon \sum (\rho_i h_i^{\text{th}}(T))  + [1-\varepsilon] \rho_{\text s} h_{\text s}) )}{\partial t} + \nabla \cdot \left( \sum h_i^{\text{th}}(T) \left( \rho_i \vec v + \vec \Phi_i \right) + \vec q_{\text{eff}} \right ) + \sum h_i^0 r_i &= 0

\end{align}
```
"""

# ╔═╡ 52afc3f1-064c-4a45-af0d-942e89e2c524
md"""
```math
\begin{align}
h_i &= h_i^0 + \int_{T_{\text{ref}}}^T c_{p,i}(\widetilde{T}) d\widetilde{T} \\
&= h_i^0 + h_i^{\text{th}}(T)
\end{align}
```
"""

# ╔═╡ f2668597-d7c1-4200-ad6b-bc6d536068ef
md"""
```math
\begin{align}
\vec q_{\text{eff}} = -\lambda_{\text{eff}} \nabla T
\end{align}
```
"""

# ╔═╡ 83cd18d8-40c8-49c8-8b01-32f0aa8cc17b
md"""
### Derivation of Separated Formulation
Total species enthalpy $h_i$ consists of two contibutions via
 - a reference enthalpy ("formation enthalpy" at reference conditions, typically 298.15 K, 1 bar) 
 - a thermal enthalpy due to the species temperature unequal the reference temperatrue.

The reference enthalpy must be considered only when species transformations through chemical reactions occur. For non-reacting mixtures it cancels out exactly.

An improvement in numerical convergence was observed after separating the (constant) "formation" enthalpy from the thermal contribution and account for the species transformation via the reaction source term ("reaction enthalpy").

In the following it is shown how to obtain the "separated formulation" from the "combined formulation" of enthalpy balance equation. For simplicity, only the gas phase mixture is considered. 

TODO: Check what happens when also considering the presence of the porous medium.
"""

# ╔═╡ 30ca1e60-ee79-4b6d-9c66-d5e3313ade3a
md"""
```math
\begin{align}
\frac{\partial \rho h}{\partial t} + \nabla \cdot ( \rho h \vec v ) + \nabla \cdot \left(  \sum_i^n h_i \vec \Phi_i \right) + \nabla \cdot \vec q  &= \frac{\partial p}{\partial t} \\
\end{align}
```
"""

# ╔═╡ 7c3f89da-9b81-4395-8ffc-00c65fdc7529
md"""
```math
\begin{align}
\frac{\partial (\sum \rho_i h_i )}{\partial t} + \nabla \cdot \left( \sum \rho_i h_i \vec v \right) + \nabla \cdot \left ( \sum_i^n h_i \vec \Phi_i \right )+ \nabla \cdot \vec q  &= \frac{\partial p}{\partial t} \\
\end{align}
```
"""

# ╔═╡ cf75ad38-c5d6-4258-ab39-710ece3b7663
md"""
```math
h_i = h_i^0 + \int_{T_{\text{ref}}}^T c_{p,i}(\widetilde{T}) d\widetilde{T} = h_i^0 + h_i^{\text{th}}(T)
```
"""

# ╔═╡ 70fbdc66-fefc-432e-abeb-05a322b34e00
md"""
```math
\begin{align}
\frac{\partial (\sum \rho_i h_i^0 )}{\partial t} + \frac{\partial (\sum \rho_i h_i^{\text{th}}(T) )}{\partial t} + \nabla \cdot \left( \sum \rho_i h_i^0 \vec v \right) + \nabla \cdot \left( \sum \rho_i h_i^{\text{th}}(T) \vec v \right) + \dots \\

+ \nabla \cdot \left ( \sum h_i^0 \vec \Phi_i \right ) + \nabla \cdot \left ( \sum h_i^{\text{th}}(T) \vec \Phi_i \right ) + \nabla \cdot \vec q  &= \frac{\partial p}{\partial t} \\
\end{align}
```
"""

# ╔═╡ cb64283a-4c74-43ea-be69-a0ee293491fd
md"""
```math
\begin{align}
\frac{\partial (\sum \rho_i h_i^0 )}{\partial t} = \sum \left[ \frac{\partial (\rho_i h_i^0 )}{\partial t} \right ] \\
\nabla \cdot \left( \sum \rho_i h_i^0 \vec v \right) = \sum \left[ \nabla \cdot \left(  \rho_i h_i^0 \vec v \right) \right] \\
\nabla \cdot \left ( \sum h_i^0 \vec \Phi_i \right ) = \sum \left[ \nabla \cdot \left ( h_i^0 \vec \Phi_i \right ) \right ]
\end{align}
```
"""

# ╔═╡ 730348aa-3ff7-43bd-8c43-e9ba6823d318
md"""
```math
\begin{align}
\sum \left( \frac{\partial (\rho_i h_i^0 )}{\partial t} + \nabla \cdot \left(  \rho_i h_i^0 \vec v \right) + \nabla \cdot \left ( h_i^0 \vec \Phi_i \right ) \right) + \dots\\
 + \frac{\partial (\sum \rho_i h_i^{\text{th}}(T) )}{\partial t} + \nabla \cdot \left( \sum \rho_i h_i^{\text{th}}(T) \vec v \right) 

+ \nabla \cdot \left ( \sum h_i^{\text{th}}(T) \vec \Phi_i \right ) + \nabla \cdot \vec q  &= \frac{\partial p}{\partial t} \\
\end{align}
```
"""

# ╔═╡ 3e53e30d-f601-4d78-9f93-037172e40504
md"""
Multiply the species mass balance with $h_i^0$:
```math
\begin{align}
\frac{\partial \rho_i}{\partial t} + \nabla \cdot \left( \rho_i \vec v \right ) + \nabla \cdot \vec \Phi_i = r_i \\
h_i^0 \frac{\partial \rho_i}{\partial t} + h_i^0 \left( \nabla \cdot \left( \rho_i \vec v \right ) \right) + h_i^0 \left( \nabla \cdot \vec \Phi_i \right )= h_i^0 r_i
\end{align}
```
Because $h_i^0$ are constant:
```math
\begin{align}
\frac{\partial (\rho_i h_i^0 )}{\partial t} + \nabla \cdot \left(  \rho_i h_i^0 \vec v \right) + \nabla \cdot \left ( h_i^0 \vec \Phi_i \right ) = h_i^0 r_i
\end{align}
```
Summing over both sides: 

```math
\begin{align}
\sum \left( \frac{\partial (\rho_i h_i^0 )}{\partial t} + \nabla \cdot \left(  \rho_i h_i^0 \vec v \right) + \nabla \cdot \left ( h_i^0 \vec \Phi_i \right ) \right) = \sum h_i^0 r_i


\end{align}
```
Leading to the finally implemented form of the enthalpy equation:
```math
\begin{align}
\frac{\partial (\sum \rho_i h_i^{\text{th}}(T) )}{\partial t} + \nabla \cdot \left( \sum \rho_i h_i^{\text{th}}(T) \vec v \right) + \nabla \cdot \left ( \sum h_i^{\text{th}}(T) \vec \Phi_i \right ) + \nabla \cdot \vec q + \sum h_i^0 r_i &= \frac{\partial p}{\partial t} \\
\end{align}
```
Or simplifying with $\frac{\partial p}{\partial t} = 0$ for the case of constant pressure:
```math
\begin{align}
\frac{\partial (\sum \rho_i h_i^{\text{th}}(T) )}{\partial t} + \nabla \cdot \left( \sum h_i^{\text{th}}(T) \left( \rho_i \vec v + \vec \Phi_i \right) + \vec q \right ) + \sum h_i^0 r_i &= 0
\end{align}
```

"""

# ╔═╡ f26fdd8c-2d66-49db-9ccd-17aa537925a1
md"""
# Extension to include Soret & Dufour effects
"""

# ╔═╡ b40f02e9-ec79-4122-9156-c469d90152f6
md"""
## Soret effect
The thermo-diffusion or Soret effect leads to a separation of mixture components when a temperature gradient is applied. Separation is strongest for species of strongly varying molecular mass, i.e. it is important when hydrogen ($\text{H}_2$) is part of the mixture.

Including the Soret effect in the model leads to an additional driving force proportional to the temperature gradient in the force-flux formulation for the Maxwell-Stefan species diffusion fluxes [1,2]:
```math
\begin{align}
	\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} + x_i \sum_{j=1}^{n} \frac{x_j}{D_{ij}} (\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}) \frac{\nabla T}{T} \right) &= -\sum_{j=1 \atop j \neq i}^{n} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
\end{align}
```

```math
\begin{align}
	\mathcal{D}_i^{\text T} &= \frac{D_i^{\text T}}{\rho_i} \\
\end{align}
```
"""

# ╔═╡ aa68910a-f86f-4e33-894e-b3d3808473cc
md"""
The notation can be made more compact by considering, that the thermal diffusion coefficients always appear in differences of pairs $\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}$, thus a matrix of so called _Newman–Soret diffusivities_ [2] can be defined that is calculated from binary Soret diffusion coefficients:
```math
\begin{align}
\mathcal{A}_{ij} &= \frac{D_i^{\text T}}{\rho_i} - \frac{D_j^{\text T}}{\rho_j} = \mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T} \\
\end{align}
```

The _thermal diffusion factors_ $\alpha_{ij}$ are commonly used which are related to the _Newman–Soret diffusivities_ via [2]:
```math
\begin{align}
\alpha_{ij} = \frac{\mathcal{A}_{ij} }{D_{ij}}
\end{align}
```
The coefficient of the temperature gradient term can be simplified:
```math
\begin{align}
x_i \sum_{j=1}^{n} \frac{x_j}{D_{ij}} (\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}) = x_i \sum_{j=1}^{n} x_j \alpha_{ij} = x_i \widetilde{\mathcal{X}_i}
\end{align}
```
Where $\widetilde{\mathcal{X}_i}$ correspond to the _rescaled thermal diffusion ratios_ as used in [3]. With this notation we now have, in accordance with [3]:

```math
\begin{align}
	\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} + x_i \widetilde{\mathcal{X}_i} \frac{\nabla T}{T} \right) &= -\sum_{j=1 \atop j \neq i}^{n} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
\end{align}
```

"""

# ╔═╡ 60cf47bf-29d4-4bad-aeb4-88f254c94dce
md"""
## Dufour effect
The Dufor effect is the "counterpart" to the Soret effect and describes the heat flux following from differing diffusive velocities (fluxes) of the different components in the mixture. The total diffusive heat flux consisting of heat conduction and the Dufour effect (summed over all species) can be expressed as [1,2,3]:
```math
\begin{align}
\vec q = -\lambda \nabla T + p \sum_{i=1}^{n} x_i \left( \sum_{j=1}^{n} \frac{x_j}{D_{ij}} (\mathcal{D}_i^{\text T} - \mathcal{D}_j^{\text T}) \right) \vec v_i
\end{align}
```
Using the rescaled thermal diffusion ratios $\widetilde{\mathcal{X}_i}$ and the species diffusice mass flux $\vec \Phi_i$ instead of the diffusive species velocity $\vec v_i$ leads to:
```math
\begin{align}
\vec q &= -\lambda \nabla T + \sum_{i=1}^{n} p x_i \widetilde{\mathcal{X}_i} \vec v_i \\ 
&= -\lambda \nabla T + \sum_{i=1}^{n} RT \frac{\widetilde{\mathcal{X}_i}}{M_i} \vec \Phi_i \\
\end{align}
```
which corresponds to the form in [3].
For the pseudo-homogeneous discription of the enthalpy balance of gas phase in the porous medium using effective thermal conductivity we get:
```math
\begin{align}
\vec q_{\text{eff}} &= -\lambda_{\text{eff}} \nabla T + \sum_{i=1}^{n} RT \frac{\widetilde{\mathcal{X}_i}}{M_i} \vec \Phi_i
\end{align}
```

The enthalpy balance equation including the Dufour effect thus takes the form:
```math
\begin{align}
\frac{\partial (\varepsilon \sum (\rho_i h_i^{\text{th}}(T)) + [1-\varepsilon] \rho_{\text s} h_{\text s}) )}{\partial t} \\ 
+ \nabla \cdot \left( \sum (h_i^{\text{th}}(T) \rho_i \vec v ) + \sum \left( h_i^{\text{th}}(T) + RT \widetilde{\mathcal{X}_i}/M_i \right) \vec \Phi_i + ( -\lambda_{\text{eff}} \nabla T) \right ) \\ 
+ \sum h_i^0 r_i \\
= 0
\end{align}
```
"""

# ╔═╡ f39f687b-1b4b-4ddd-864d-cb33e28d3c78
md"""
# References
1) Kuiken, Gerard D. C. (1994): Thermodynamics of irreversible processes. Applications to diffusion and rheology. New York, NY: Wiley (Wiley tutorial series in theoretical chemistry).
1) Van‐Brunt, Alexander; Farrell, Patrick E.; Monroe, Charles W. (2022): Consolidated theory of fluid thermodiffusion. In: AIChE Journal 68 (5), Artikel e17599. DOI: 10.1002/aic.17599  .
1) Giovangigli, Vincent (2016): Solutions for Models of Chemically Reacting Mixtures. In: Yoshikazu Giga und Antonin Novotny (Hg.): Handbook of Mathematical Analysis in Mechanics of Viscous Fluids. Cham: Springer International Publishing, S. 1–52.
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
# ╟─f26fdd8c-2d66-49db-9ccd-17aa537925a1
# ╟─b40f02e9-ec79-4122-9156-c469d90152f6
# ╟─aa68910a-f86f-4e33-894e-b3d3808473cc
# ╟─60cf47bf-29d4-4bad-aeb4-88f254c94dce
# ╟─f39f687b-1b4b-4ddd-864d-cb33e28d3c78
