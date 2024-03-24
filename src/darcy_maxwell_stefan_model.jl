# abstract type AbstractModelData end

@doc raw"""
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
	\frac{\partial \epsilon \rho_i}{\partial t} + \nabla \cdot \left( \vec J_i + \rho_i \vec v \right ) - r_i(\varrho) &= 0 ~, \qquad i = 1 ... n \\
		\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= -\sum_{j=1 \atop j \neq i}^{n} \frac{w_j \vec J_i-w_i \vec J_j}{D_{ij} M_i M_j} \\
		\sum_{i=1}^n x_i &= 1
\end{align}
```

where $\rho$ is the (total) mixture density, $\vec v$ is the mass-averaged (barycentric)  mixture velocity calculated with the Darcy equation,
$x_i$, $w_i$ and $M_i$ are the molar fraction, mass fraction and molar mass of species $i$ respectively,
$\vec J_i$ is the __diffusive__ mass flux of species $i$ ($\frac{\text{kg}}{\text{m}^2 \text{s}}$)
and $r_i$ is the species mass volumetric source/sink ($\frac{\text{kg}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
The __convective__ mass flux of species $i$ is the product of the superficial mean velocity with the partial mass density $\rho_i \vec v$.
"""
function DMS_Info_isothermal() end

@doc raw"""
## Thermal Energy Transport

Enthalpy equation for __gas phase only__.
Considers convective-diffusive transport of thermal energy,
including "thermal drift" from convective-diffusive species fluxes $\vec J_i + \rho_i \vec v$ carrying enthalpy.

Formulation based on separation of "reference  enthalpy" and "thermal enthalpy".
In the documentation.jl notebook see section "Derivation of Separated Formulation".
```math
\begin{align}
\frac{\partial (\sum \rho_i h_i^{\text{th}}(T) )}{\partial t} + \nabla \cdot \left( \sum h_i^{\text{th}}(T) \left( \rho_i \vec v + \vec J_i \right) + \vec q \right ) + \sum h_i^0 r_i &= 0
\end{align}
```

```math
\begin{align}
\vec q = -\lambda \nabla T
\end{align}
```

Enthalpy equation for __gas + solid phase__.
It is formulated with effective transport parameters
that are a consequence of the quasi-homogeneous phase approach.
Formulation based on separation of "reference  enthalpy" and "thermal enthalpy".
In the documentation.jl notebook see section "Derivation of Separated Formulation".

```math
\begin{align}
\frac{\partial (\varepsilon \sum (\rho_i h_i^{\text{th}}(T))  + [1-\varepsilon] \rho_{\text s} h_{\text s}) )}{\partial t} + \nabla \cdot \left( \sum h_i^{\text{th}}(T) \left( \rho_i \vec v + \vec J_i \right) + \vec q_{\text{eff}} \right ) + \sum h_i^0 r_i &= 0

\end{align}
```

```math
\begin{align}
h_i &= h_i^0 + \int_{T_{\text{ref}}}^T c_{p,i}(\widetilde{T}) d\widetilde{T} \\
&= h_i^0 + h_i^{\text{th}}(T)
\end{align}
```

```math
\begin{align}
\vec q_{\text{eff}} = -\lambda_{\text{eff}} \nabla T
\end{align}
```
"""
function DMS_Info_thermal() end

@doc raw"""
Assemble symmetric Maxwell-Stefan matrix M expressing the flux-force relations in the Maxwell-Stefan
framework for multicomponent diffusion.
"""
# function M_matrix!(M, W, D, data)
# 	ng=ngas(data)
# 	@inbounds for i=1:(ng-1)
# 		for j=1:(ng-1)
# 			M[i,j] = zero(eltype(W))
# 		end
# 		 for j=1:ng
# 			if i != j
# 				M[i,i] -= W[j]/D[i,j]
# 				if j == ng
# 					 for k=1:(ng-1)
# 						M[i,k] -= W[i]/D[i,ng]
# 					end
# 				else
# 					M[i,j] += W[i]/D[i,j]
# 				end
# 			end
# 		end
# 	end			
# end
function M_matrix!(M, W, D, data)
	ng=ngas(data)
	(;m) = data
	@inbounds for i=1:(ng-1)
		for j=1:(ng-1)
			M[i,j] = zero(eltype(W))
		end
		 for j=1:ng
			if i != j
				M[i,i] -= W[j]/(D[i,j]*m[i]*m[j])
				if j == ng
					 for k=1:(ng-1)
						M[i,k] -= W[i]/(D[i,ng]*m[i]*m[ng])
					end
				else
					M[i,j] += W[i]/(D[i,j]*m[i]*m[j])
				end
			end
		end
	end			
end
@doc raw"""
Assemble symmetric Maxwell-Stefan diffusivity matrix D. Account for porous material
through the constriction and tourtuosity factor γ_τ which lowers the diffusivities
 ~1 order of magnitude compared to diffusion through free space.
"""
function D_matrix!(D, T, p, data)
	(;m,γ_τ,constant_properties,constant_binary_diff_coeffs)=data
	ng=ngas(data)
	ind = 1
	@inbounds for i=1:(ng-1) # i: row index
		for j=(i+1):ng # j: col index
            if !constant_properties
			    Dij = binary_diff_coeff_gas(data.Fluids[i], data.Fluids[j], T, p)
			    # Dij *= m[i]*m[j]*γ_τ # eff. diffusivity
            else
                # Dij = 2.0e-5*m[i]*m[j]
				Dij = constant_binary_diff_coeffs[ind]
				# Dij *= m[i]*m[j]
            end
			# Dij *= m[i]*m[j]*γ_τ # eff. diffusivity
			Dij *= γ_τ # eff. diffusivity
			D[i,j] = Dij
			D[j,i] = Dij # M-S diffusion matrix is symmetric
			ind += 1
		end
	end
end

@doc raw"""
Combined function to assemble symmetric Maxwell-Stefan diffusivity and Newman-Soret coefficient matrices D and A.
For M-S diffusivity, account for porous material through the constriction and tourtuosity factor γ_τ which lowers the diffusivities
 ~1 order of magnitude compared to diffusion through free space.
"""
function D_A_matrices!(D, A, T, p, data)
	(;m,γ_τ,constant_properties,constant_binary_diff_coeffs,constant_newman_soret_diff_coeffs)=data
	ng=ngas(data)
	ind = 1
	@inbounds for i=1:(ng-1) # i: row index
		for j=(i+1):ng # j: col index
            # if !constant_properties
			# Dij = binary_diff_coeff_gas(data.Fluids[i], data.Fluids[j], T, p)
			# # Dij *= m[i]*m[j]*γ_τ # eff. diffusivity
			# Aij = constant_newman_soret_diff_coeffs[ind] # in absence of calculation methods use constant values
            # else
                # Dij = 2.0e-5*m[i]*m[j]
				Dij = constant_binary_diff_coeffs[ind]
				# Dij *= m[i]*m[j]
				Aij = constant_newman_soret_diff_coeffs[ind]
            # end
			# Dij *= m[i]*m[j]*γ_τ # eff. diffusivity
			Dij *= γ_τ # eff. diffusivity
			D[i,j] = Dij
			D[j,i] = Dij # M-S diffusion matrix is symmetric
			A[i,j] = Aij
			A[j,i] = -Aij # Newman-Soret coefficient matrix is anti-symmetric
			ind += 1
		end
	end
end


@doc raw"""
Mixture mass flow (bulk convective mass flow) through the pore space of the porous medium. Mixture mass averaged (barycentric) velocity for convective bulk mass flow through the porous medium is calculated from Darcy equation. 

```math
	\vec v  = -\frac{\kappa}{\mu} \vec \nabla p
```
"""

function DarcyVelo(u,data,mu)
	(;ip,perm) = data
	-perm/mu*(u[ip,1]-u[ip,2])	
end

function MoleFrac!(X,u::VoronoiFVM.EdgeUnknowns,data)
	@inbounds for i=1:ngas(data)
		X[i] = 0.5*(u[i,1]+u[i,2])
	end
	nothing
end

function MoleFrac!(X,u::VoronoiFVM.BNodeUnknowns,data)
	@inbounds for i=1:ngas(data)
		X[i] = u[i]
	end
	nothing
end

function MoleFrac!(X,u::VoronoiFVM.NodeUnknowns,data)
	@inbounds for i=1:ngas(data)
		X[i] = u[i]
	end
	nothing
end

function MassFrac!(X,W,data)
	(;m) = data
	@inline mmix = molarweight_mix(X,data)
	@inbounds for i=1:ngas(data)
		W[i] = X[i]*m[i]/mmix
	end
	nothing
end

function ThermalDiffRatio!(TDR,X,A,D,data)
	@inbounds for i=1:ngas(data)
		TDR[i] = zero(eltype(X))		
		for j=1:ngas(data)
			if i != j
				TDR[i] += X[j]*A[i,j]/D[i,j]
			end
		end
	end
	nothing
end

@doc raw"""
Flux function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_flux(f,u,edge,data)
	(;m,ip,iT,dt_hf_enth,solve_T_equation,Tamb,Fluids,include_Soret_Dufour)=data
	ng=ngas(data)
		
	F = MVector{ng-1,eltype(u)}(undef)
	X = MVector{ng,eltype(u)}(undef)
	W = MVector{ng,eltype(u)}(undef)
	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)
	D = MMatrix{ng,ng,eltype(u)}(undef)
	# allocate storage for Thermo-Diffusion parameters
	A = MMatrix{ng,ng,eltype(u)}(undef)
	TDR = MVector{ng,eltype(u)}(undef)

	pm = 0.5*(u[ip,1]+u[ip,2])
    Tm = solve_T_equation ? 0.5*(u[iT,1]+u[iT,2]) : one(eltype(u))*Tamb
	c = pm/(ph"R"*Tm)
	
	δp = u[ip,1]-u[ip,2]
	
	@inline MoleFrac!(X,u,data)
	@inline mmix = molarweight_mix(X,data)
	@inline MassFrac!(X,W,data)
	
	if include_Soret_Dufour
		@inline D_A_matrices!(D, A, Tm, pm, data)
		@inline ThermalDiffRatio!(TDR, X, A, D, data)	
	else
		@inline D_matrix!(D, Tm, pm, data)	
	end
	@inline mumix, lambdamix = dynvisc_thermcond_mix(data, Tm, X)
		
	rho = c*mmix
	v = DarcyVelo(u,data,mumix)
	
	f[ip] = -rho*v

	@inline M_matrix!(M, W, D, data)
	
	@inbounds for i=1:(ng-1)

		F[i] = ( u[i,1]-u[i,2] + (X[i]-W[i])*δp/pm )*c/mmix
		if include_Soret_Dufour
			F[i] += (X[i]*TDR[i]*(u[iT,1]-u[iT,2])/Tm )*c/mmix # Soret effect
		end	
	end				

	@inline inplace_linsolve!(M,F)

	enthalpy_flux = zero(eltype(u))
	mass_flux = zero(eltype(u))
	@inbounds for i=1:(ng-1)

		# convective -diffusive species mass flux
		f[i] = -(F[i] + c*X[i]*m[i]*v)

		if solve_T_equation

			mass_flux += f[i]

			enthalpy_flux += f[i] * enthalpy_gas_thermal(Fluids[i], Tm) / m[i]

			if include_Soret_Dufour
				enthalpy_flux += (-F[i])*ph"R"*Tm*TDR[i]/m[i] # Dufour effect
			end

			
		end
	end	
	
    if solve_T_equation

		enthalpy_flux += (f[ip] - mass_flux) * enthalpy_gas_thermal(Fluids[ng], Tm) / m[ng] # species n

		if include_Soret_Dufour
			enthalpy_flux += (f[ip] - mass_flux) *ph"R"*Tm*TDR[ng]/m[ng] # Dufour effect species n
		end

        lambda_bed=kbed(data,lambdamix)*lambdamix

		f[iT] = lambda_bed*(u[iT,1]-u[iT,2]) + enthalpy_flux * ramp(edge.time; du=(0.0,1), dt=dt_hf_enth)
    end
end

@doc raw"""
Reaction function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_reaction(f,u,node,data)
	(;m,ip,iT,idpdt,is_reactive,solve_T_equation,include_dpdt,catalyst_regions)=data
	ng=ngas(data)

	f[iT] = zero(eltype(u))
	if node.region in catalyst_regions && is_reactive # catalyst layer
		(;lcat,kinpar,Treac,)=data
		(;nuij)=kinpar
		
		pi = MVector{ng,eltype(u)}(undef)
		for i=1:ng
            pi[i] = u[ip]*u[i]
		end

        T = solve_T_equation ? u[iT] : one(eltype(u))*Treac
        RR = @inline -lcat*ri(data,T,pi)
        #RR = @inline -lcat*ri(data,u[iT],pi)
		

		for i=1:ng
			f[i] = zero(eltype(u))
			for j=1:nreac(kinpar)
				f[i] += nuij[(j-1)*ng+i] * RR[j] * m[i]
				if solve_T_equation
					# f[iT] -= nuij[(j-1)*ng+i] * RR[j] * enthalpy_gas(kinpar.Fluids[i], T)
					f[iT] -= nuij[(j-1)*ng+i] * RR[j] * kinpar.Fluids[i].ΔHform 
				end
			end			
		end
	end

	if solve_T_equation && include_dpdt
		f[idpdt] = -u[idpdt]
		f[iT] += u[idpdt] # dpdt source term in enthalpy eq.		
	end

	f[ng] = zero(eltype(u))
	for i=1:ng
		f[ng] += u[i]
	end
	f[ng] = f[ng] - 1.0
end

@doc raw"""
Storage function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_storage(f,u,node,data)
	(;ip,iT,idpdt,Tamb,m,poros,rhos,cs,solve_T_equation,include_dpdt,Fluids)=data
	ng=ngas(data)

    T = solve_T_equation ? u[iT] : one(eltype(u))*Tamb
    c = u[ip]/(ph"R"*T)
    mmix = zero(eltype(u))
	enthalpy_gas = zero(eltype(u))
	for i=1:ng
		f[i]=c*u[i]*m[i]*poros
        mmix += u[i]*m[i]
		if solve_T_equation
			enthalpy_gas += c*u[i]*enthalpy_gas_thermal(Fluids[i], T)
		end
	end
	
	# total density / total pressure
	f[ip] = mmix*c*poros

    if solve_T_equation       

		# f[iT] = u[iT] * rhos*cs*(1-poros) + poros * enthalpy_gas
		f[iT] = (u[iT]-298.15) * rhos*cs*(1-poros) + poros * enthalpy_gas

		if include_dpdt
			f[idpdt] = u[ip]
		end

    end

	
end

@doc raw"""
Outflow boundary function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_boutflow(f,u,edge,data)
	(;iT,ip,m,dt_hf_enth,Tamb,solve_T_equation,Fluids)=data
	ng=ngas(data)

	k=outflownode(edge)
	pout = u[ip,k]
    Tout = solve_T_equation ? u[iT,k] : one(eltype(u))*Tamb
    cout = pout/(ph"R"*Tout)

	X = MVector{ng,eltype(u)}(undef)
	
	for i=1:ng
		X[i] = u[i,k]
	end
	
    @inline mumix, _ = dynvisc_thermcond_mix(data, Tout, X)
	v = DarcyVelo(u,data,mumix)
	
	for i=1:(ng-1)
		f[i] = v*cout*u[i,k]*m[i] # species mass flux at outflow
	end

    if solve_T_equation
		
		hout = zero(eltype(u))
        @inbounds for i=1:ng
            hout += X[i] * enthalpy_gas_thermal(Fluids[i], Tout)
        end

		r_hf_enth = v * cout * hout
        # f[iT] = r_hf_enth
		f[iT] = r_hf_enth * ramp(edge.time; du=(0.0,1.0), dt=dt_hf_enth)
    end
end

