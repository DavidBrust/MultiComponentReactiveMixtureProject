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
function DMS_Info_isothermal() end

@doc raw"""
## Thermal Energy Transport

Enthalpy equation for __gas phase only__.
Considers convective-diffusive transport of thermal energy,
including "thermal drift" from convective-diffusive species fluxes $\vec \Phi_i + \rho_i \vec v$ carrying enthalpy.

Formulation based on separation of "reference  enthalpy" and "thermal enthalpy".
In the documentation.jl notebook see section "Derivation of Separated Formulation".
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

Enthalpy equation for __gas + solid phase__.
It is formulated with effective transport parameters
that are a consequence of the quasi-homogeneous phase approach.
Formulation based on separation of "reference  enthalpy" and "thermal enthalpy".
In the documentation.jl notebook see section "Derivation of Separated Formulation".

```math
\begin{align}
\frac{\partial (\varepsilon \sum (\rho_i h_i^{\text{th}}(T))  + [1-\varepsilon] \rho_{\text s} h_{\text s}) )}{\partial t} + \nabla \cdot \left( \sum h_i^{\text{th}}(T) \left( \rho_i \vec v + \vec \Phi_i \right) + \vec q_{\text{eff}} \right ) + \sum h_i^0 r_i &= 0

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
function M_matrix!(M, W, D, data)
	ng=ngas(data)
	@inbounds for i=1:(ng-1)
		for j=1:(ng-1)
			M[i,j] = zero(eltype(W))
		end
		 for j=1:ng
			if i != j
				M[i,i] -= W[j]/D[i,j]
				if j == ng
					 for k=1:(ng-1)
						M[i,k] -= W[i]/D[i,ng]
					end
				else
					M[i,j] += W[i]/D[i,j]
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
	(;m,γ_τ,constant_properties)=data
	ng=ngas(data)
	@inbounds for i=1:(ng-1)
		for j=(i+1):ng
            if !constant_properties
			    Dji = binary_diff_coeff_gas(data.Fluids[j], data.Fluids[i], T, p)
			    Dji *= m[i]*m[j]*γ_τ # eff. diffusivity
            else
                Dji = 2.0e-5*m[i]*m[j]
            end
			D[j,i] = Dji
			D[i,j] = Dji
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

@doc raw"""
Flux function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_flux(f,u,edge,data)
	(;m,ip,iT,dt_hf_enth,solve_T_equation,Tamb,Fluids)=data
	ng=ngas(data)
		
	F = MVector{ng-1,eltype(u)}(undef)
	X = MVector{ng,eltype(u)}(undef)
	W = MVector{ng,eltype(u)}(undef)
	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)
	D = MMatrix{ng,ng,eltype(u)}(undef)

	pm = 0.5*(u[ip,1]+u[ip,2])
    Tm = solve_T_equation ? 0.5*(u[iT,1]+u[iT,2]) : one(eltype(u))*Tamb
	c = pm/(ph"R"*Tm)
	
	δp = u[ip,1]-u[ip,2]
	
	@inline MoleFrac!(X,u,data)
	@inline mmix = molarweight_mix(X,data)
	@inline MassFrac!(X,W,data)
	
	
	@inline D_matrix!(D, Tm, pm, data)
	@inline mumix, lambdamix = dynvisc_thermcond_mix(data, Tm, X)
		
	rho = c*mmix
	v = DarcyVelo(u,data,mumix)
	
	f[ip] = -rho*v

	@inline M_matrix!(M, W, D, data)
	
	@inbounds for i=1:(ng-1)
		F[i] = ( u[i,1]-u[i,2] + (X[i]-W[i])*δp/pm )*c/mmix
	end				

	@inline inplace_linsolve!(M,F)

	enthalpy_flux = zero(eltype(u))
	mass_flux = zero(eltype(u))
	@inbounds for i=1:(ng-1)
		f[i] = -(F[i] + c*X[i]*m[i]*v)
		if solve_T_equation
			mass_flux += f[i]
			# !!! DEBUG !!!
			# enthalpy_flux += f[i] * enthalpy_gas(Fluids[i], Tm) / m[i]
			# enthalpy_flux += f[i] * heatcap_gas(Fluids[i], Tm) / m[i] * (Tm - Tref)
			enthalpy_flux += f[i] * enthalpy_gas_thermal(Fluids[i], Tm) / m[i]
			# !!! DEBUG !!!
		end
	end	
	
    if solve_T_equation
		# !!! DEBUG !!!
		# enthalpy_flux += (f[ip] - mass_flux) * enthalpy_gas(Fluids[ng], Tm) / m[ng]# species n
		# enthalpy_flux += (f[ip] - mass_flux) * heatcap_gas(Fluids[ng], Tm) / m[ng] * (Tm - Tref)# species n
		enthalpy_flux += (f[ip] - mass_flux) * enthalpy_gas_thermal(Fluids[ng], Tm) / m[ng] # species n
		# !!! DEBUG !!!
        lambda_bed=kbed(data,lambdamix)*lambdamix
        # @inline hf_conv = f[ip] * enthalpy_mix(data, Tm, X) / mmix * ramp(edge.time; du=(0.0,1), dt=dt_hf_enth) 
        # Bp,Bm = fbernoulli_pm(hf_conv/lambda_bed/Tm)
        # f[iT] = lambda_bed*(Bm*u[iT,1]-Bp*u[iT,2])

		# !!! DEBUG !!!
		f[iT] = lambda_bed*(u[iT,1]-u[iT,2]) + enthalpy_flux * ramp(edge.time; du=(0.0,1), dt=dt_hf_enth) 
		# f[iT] = lambda_bed*(u[iT,1]-u[iT,2]) + enthalpy_flux
		# f[iT] = lambda_bed*(u[iT,1]-u[iT,2])
		# !!! DEBUG !!!
    end
end

@doc raw"""
Reaction function definition for use with VoronoiFVM.jl for the Darcy-Maxwell-Stefan
    (DMS) model for Multi-component gas transport in porous media.
"""
function DMS_reaction(f,u,node,data)
	(;m,ip,iT,is_reactive,catalyst_regions)=data
	ng=ngas(data)

	if node.region in catalyst_regions && is_reactive # catalyst layer
		(;lcat,kinpar,iT,Treac,solve_T_equation)=data
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
					f[iT] -= nuij[(j-1)*ng+i] * RR[j] * enthalpy_gas(kinpar.Fluids[i], T)
				end
			end			
		end
	end

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
	(;ip,iT,Tamb,m,poros,rhos,cs,solve_T_equation)=data
	ng=ngas(data)

    T = solve_T_equation ? u[iT] : one(eltype(u))*Tamb
    c = u[ip]/(ph"R"*T)
    mmix = zero(eltype(u))
	for i=1:ng
		f[i]=c*u[i]*m[i]*poros
        mmix += u[i]*m[i]
	end
	
	# total pressure
	f[ip] = mmix*c*poros

    if solve_T_equation
        X=MVector{ng,eltype(u)}(undef)
        @inline MoleFrac!(X,u,data)
        @inline cpmix = heatcap_mix(data, T, X)
		# !!! DEBUG !!!
		f[iT] = u[iT] * (rhos*cs*(1-poros) + cpmix*c*poros)
		# f[iT] = (u[iT]-Tref) * (rhos*cs*(1-poros) + cpmix*c*poros)
		# !!! DEBUG !!!
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
        # @inline r_hf_enth = v *cout * enthalpy_mix(data, Tout, X) * ramp(edge.time; du=(0.0,1.0), dt=dt_hf_enth) # enthalpy heat flux
		# !!! DEBUG !!!
		
		hout = zero(eltype(u))
        @inbounds for i=1:ng
            hout += X[i] * enthalpy_gas_thermal(Fluids[i], Tout)
        end

		#@inline hout = heatcap_mix(data, Tout, X) * (Tout - Tref)
		r_hf_enth = v * cout * hout
        # f[iT] = r_hf_enth
		f[iT] = r_hf_enth * ramp(edge.time; du=(0.0,1.0), dt=dt_hf_enth)
		# f[iT] = zero(eltype(u))
		# !!! DEBUG !!!
    end
end

