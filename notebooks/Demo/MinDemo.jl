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

# ╔═╡ c21e1942-628c-11ee-2434-fd4adbdd2b93
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,"../.."))
	using Revise
	using VoronoiFVM
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve, LinearSolve
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots, Printf
	using PlutoUI, Colors

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 2c4b16ab-8f60-467b-b608-2fea9fbc741c
md"""
# Introduction
Demo Notebook with 3 ideal gas species and constant properties to demonstrate isothermal Multicomponent species transport.

Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=false))
"""

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
PlutoUI.TableOfContents(title="Demo")

# ╔═╡ 83fa22fa-451d-4c30-a4b7-834974245996
function grid1D()
	X=(0:0.02:1)*ufac"cm"
	grid=simplexgrid(X)
	# catalyst region
	cellmask!(grid,[0.4]*ufac"cm",[0.6]*ufac"cm",2)	
	grid
end

# ╔═╡ 4dae4173-0363-40bc-a9ca-ce5b4d5224cd
function grid2D()
	X=(0:0.1:1)*ufac"cm"
	Y=(0:0.1:1)*ufac"cm"
	grid=simplexgrid(X,Y)

	# catalyst region
	cellmask!(grid,[0.4,0.3].*ufac"cm",[0.6,0.7].*ufac"cm",2)
	#bfacemask!(grid, [1,0.2].*ufac"cm",[1,0.8].*ufac"cm",5)
	
	grid
end

# ╔═╡ 561e96e2-2d48-4eb6-bb9d-ae167a622aeb
function grid3D()
	X=(0:0.1:1)*ufac"cm"
	Y=(0:0.1:1)*ufac"cm"
	Z=(0:0.1:1)*ufac"cm"
	grid=simplexgrid(X,Y,Z)

	# catalyst region
	cellmask!(grid,[0.3,0.4,0.3].*ufac"cm",[0.7,0.6,0.7].*ufac"cm",2)
	bfacemask!(grid, [0.2,1,0.2].*ufac"cm",[0.8,1,0.8].*ufac"cm",7)
	
	grid
end

# ╔═╡ 107a6fa3-60cb-43f0-8b21-50cd1eb5065a
const dim = 1

# ╔═╡ 4e05ab31-7729-4a4b-9c14-145118477715
if dim == 3
	@bind xcut Slider(linspace(0,1,21)*ufac"cm",show_value=true,default=0.5*ufac"cm")
end

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
let
	if dim == 1
		gridplot(grid1D(), resolution=(600,200))
	elseif dim == 2
		gridplot(grid2D())
	else
		gridplot(grid3D(); xplane=xcut, show=true, outlinealpha=0.0 )
	end
end

# ╔═╡ 832f3c15-b75a-4afe-8cc5-75ff3b4704d6
begin
	if dim == 1
		const Γ_left = 1
		const Γ_right = 2
	elseif dim == 2
		const Γ_bottom = 1
		const Γ_right = 2
		const Γ_top = 3
		const Γ_left = 4
		#const Γ_right = 5
	else
		const Γ_left = 1
		const Γ_front = 2
		#const Γ_right = 3
		const Γ_right = 7
		const Γ_back = 4
		const Γ_bottom = 5
		const Γ_top = 6
	end
end;

# ╔═╡ 0fadb9d2-1ccf-4d44-b748-b76d911784ca
md"""
## Mass Continuity
```math
\begin{align}
	\frac{\partial \rho}{\partial t} + \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\end{align}
```
"""

# ╔═╡ b94513c2-c94e-4bcb-9342-47ea48fbfd14
md"""
## Species Mass Transport
```math
\begin{align}
	\frac{\partial \rho_i}{\partial t} + \nabla \cdot \left( \vec \Phi_i + \rho_i \vec v \right ) - R_i &= 0 ~, \qquad i = 1 ... \nu \\
		\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= -\sum_{j=1 \atop j \neq i}^{\nu} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
		\sum_{i=1}^\nu x_i &= 1
\end{align}
```
"""

# ╔═╡ c886dd12-a90c-40ab-b9d0-32934c17baee
md"""
where $\rho$ is the (total) mixture density, $\vec v$ is the mass-averaged (barycentric)  mixture velocity calculated with the Darcy equation, $x_i$, $w_i$ and $M_i$ are the molar fraction, mass fraction and molar mass of species $i$ respectively, $\vec \Phi_i$ is the mass flux of species $i$ ($\frac{\text{kg}}{\text{m}^2 \text{s}}$) and $R_i$ is the species mass volumetric source/sink ($\frac{\text{mol}}{\text{m}^3 \text{s}}$) of gas phase species $i$.
"""

# ╔═╡ 3440d4d8-3e03-4ff3-93f1-9afd7aaf9c41
md"""
Reaction leading to increase in moles.

$A \rightarrow 3 B$
"""

# ╔═╡ 3bb2deff-7816-4749-9f1e-c1e451372b1e
function reaction(f,u,node,data)
	(;p,m,ip,Tamb,isreactive)=data
	ng=ngas(data)
	if node.region == 2 && isreactive == 1 # catalyst layer		
		
		# R1: u2 -> 3 u1 
		k1 = 5.0
		c = u[ip]/(ph"R"*Tamb)
		r = k1 * u[2]*c

		f[1] = -3*r*m[1]
		f[2] = r*m[2]		
	end
	
	for i=1:ng
		f[ng] += u[i]
	end
	f[ng] = f[ng] - 1.0
end

# ╔═╡ 4af1792c-572e-465c-84bf-b67dd6a7bc93
function storage(f,u,node,data)
	(;Tamb,ip,m)=data
	ng=ngas(data)
	mmix = zero(eltype(u))
	c = u[ip]/(ph"R"*Tamb)
	
	for i=1:ng
		f[i]=c*u[i]*m[i]
		mmix += u[i]*m[i]
	end
	
	# total pressure
	f[ng+1] = mmix*c
end

# ╔═╡ 5f88937b-5802-4a4e-81e2-82737514b9e4
function bcond(f,u,bnode,data)
	(;p,ip,mfluxin,W0,dt_mf,ng)=data

	r_mfluxin = mfluxin*ramp(bnode.time; du=(0.1,1), dt=dt_mf)
	for i=1:(ng-1)
		boundary_neumann!(f,u,bnode, species=i,region=Γ_left,value=r_mfluxin*W0[i])
	end	
	boundary_neumann!(f,u,bnode, species=ip, region=Γ_left, value=r_mfluxin)
	boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_right,value=p)
end

# ╔═╡ 05949759-2bb9-475b-b2f4-900b32c30e00
md"""
## Molar Fractions
"""

# ╔═╡ d4aaf146-e551-444a-a2f1-e70b05104b53
md"""
## Pressure, Density, Flow
1) Total Pressure
2) Density
3) Velocity - X
4) Velocity - Y
"""

# ╔═╡ 389a4798-a9ee-4e9c-8b44-a06201b4c457
function boutflow(f,u,edge,data)
	(;Tamb,ip,m)=data
	ng=ngas(data)

	k=outflownode(edge)

	pout = u[ip,k]
	cout = pout/(ph"R"*Tamb)
	mumix = 2.0e-5*ufac"Pa*s"
	v = DarcyVelo(u,data,mumix)
	
	for i=1:(ng-1)
		f[i] = v * cout*u[i,k]*m[i]
	end	
end

# ╔═╡ 1f2d79f6-5511-4286-907f-a4cf8c547f8e
function D_matrix!(D,data)
	ng=ngas(data)
	(;m) = data
	@inbounds for i=1:(ng-1)
		for j=(i+1):ng
            Dji = one(eltype(D)) *1.0e-5*ufac"m^2/s" * m[i]*m[j]
			@views D[j,i] = Dji
			@views D[i,j] = Dji
		end
	end
end

# ╔═╡ 5547d7ad-dd58-4b00-8238-6e1abb32874e
function flux(f,u,edge,data)
	(;m,ip,Tamb)=data
	ng=ngas(data)

	F = MVector{ng-1,eltype(u)}(undef)
	X = MVector{ng,eltype(u)}(undef)
	W = MVector{ng,eltype(u)}(undef)
	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)
	D = MMatrix{ng,ng,eltype(u)}(undef)
	
	δp = u[ip,1]-u[ip,2]
	pm = 0.5*(u[ip,1]+u[ip,2])
	c = pm/(ph"R"*Tamb)

	@inline MoleFrac!(X,u,data)
	@inline mmix = molarweight_mix(X,data)
	@inline MassFrac!(X,W,data)
	
	rho = c*mmix
	mumix = 2.0e-5*ufac"Pa*s"
	v = DarcyVelo(u,data,mumix)
	
	
	f[ip] = -rho*v
	
	@inline D_matrix!(D, data)
	@inline M_matrix!(M, W, D, data)


	@inbounds for i=1:(ng-1)
		F[i] = ( u[i,1]-u[i,2] + (X[i]-W[i])*δp/pm )*c/mmix
	end				
	

	@inline inplace_linsolve!(M,F)

	@inbounds for i=1:(ng-1)
		f[i] = -(F[i] + c*X[i]*m[i]*v)
		#f[i] = -F[i]
		#f[i] = u[i,1]-u[i,2] # quadr. conv
	end

end

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
begin
	
	if dim == 1
		mygrid=grid1D()
		strategy = nothing
		times=[0,20.0]
	elseif dim == 2
		mygrid=grid2D()
		strategy = nothing
		times=[0,12.0]
	else
		mygrid=grid3D()
		strategy = GMRESIteration(UMFPACKFactorization())
		times=[0,1.0]
	end
	mydata = ReactorData(;
		ng = 3,
		Tamb = 273.15,
		m = [2.0,6.0,21.0]*ufac"g/mol",
		X0 = [0.2, 0.8, 0.0],
		mfluxin = 0.01*ufac"kg/(m^2*s)"
	)
	
	(;p,ip,X0)=mydata
	ng=ngas(mydata)
	
	sys=VoronoiFVM.System( 	mygrid;
							data=mydata,
							flux=flux,
							reaction=reaction,
							storage=storage,
							bcondition=bcond,
							boutflow=boutflow,
							outflowboundaries=[Γ_right],
							assembly=:edgewise
							)
	
	enable_species!(sys; species=collect(1:(ng+1))) # gas phase species pi & ptotal	
	
	inival=unknowns(sys)

	inival[ip,:].=p
	for i=1:ng
		inival[i,:] .= 1.0/ng
	end

	control = SolverControl(strategy, sys;)
		control.handle_exceptions=true
		control.Δu_opt=1.0e5

	if RunSim
		solt=solve(sys;inival=inival,times,control,verbose="nae")
	end
end;

# ╔═╡ f798e27a-1d7f-40d0-9a36-e8f0f26899b6
@bind t Slider(solt.t,show_value=true,default=solt.t[end])

# ╔═╡ 5588790a-73d4-435d-950f-515ae2de923c
sol = solt(t);

# ╔═╡ 111b1b1f-51a5-4069-a365-a713c92b79f4
let
	(;ip,p) = mydata
	if dim == 1
		vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300))
		scalarplot!(vis, mygrid, sol[1,:], clear=false, label="x1")
		scalarplot!(vis, mygrid, sol[2,:], clear=false, color=:red, label="x2")
		scalarplot!(vis, mygrid, sol[3,:], clear=false, color=:blue, label="x3")
	elseif dim == 2
		vis=GridVisualizer(layout=(1,3), resolution=(700,300))
		scalarplot!(vis[1,1], mygrid, sol[1,:], clear=false, title="x1")
		scalarplot!(vis[1,2], mygrid, sol[2,:], clear=false, title="x2")
		scalarplot!(vis[1,3], mygrid, sol[3,:], clear=false, title="x3")
	else
		vis=GridVisualizer(layout=(3,1), resolution=(400,1200), outlinealpha=0.0)
		scalarplot!(vis[1,1], mygrid, sol[1,:], clear=false, label="x1")
		scalarplot!(vis[2,1], mygrid, sol[2,:], clear=false, color=:red, label="x2")
		scalarplot!(vis[3,1], mygrid, sol[3,:], clear=false, color=:blue, label="x3")
	end
	
	reveal(vis)
end

# ╔═╡ de69f808-2618-4add-b092-522a1d7e0bb7
let
	(;p,m,ip,Tamb, mfluxin) = mydata
	ng = ngas(mydata)
	mmix = []
	for j in 1:length(sol[1,:])
		_mmix=0
		for i=1:ng
			_mmix += sol[i,j]*m[i]
		end
		push!(mmix, _mmix)
	end
	
	w1 = sol[1,:]*m[1] ./mmix
	w2 = sol[2,:]*m[2] ./mmix
	w3 = sol[3,:]*m[3] ./mmix

	ps = sol[ip,:]
	rho = @. ps * mmix /(ph"R"*Tamb)
	
	if dim == 1
		vis=GridVisualizer(legend=:lt, resolution=(600,400), layout = (2, 1))
		scalarplot!(vis[1,1], mygrid, w1, clear=false, label="w1")
		scalarplot!(vis[1,1], mygrid, w2, clear=false, color=:red, label="w2")
		scalarplot!(vis[1,1], mygrid, w3, clear=false, color=:blue, label="w3")

		p0 = sol[ip,1]
		rho0 = @. p0 * mmix[1] /(ph"R"*Tamb)
		scalarplot!(vis[2,1], mygrid, rho/rho0, clear=false, label="Rho / Rho0")
		scalarplot!(vis[2,1], mygrid, rho0./rho, clear=false, color=:red, label="v / v0")
		scalarplot!(vis[2,1], mygrid, ps/p0, clear=false, color=:blue, label="p / p0")
	elseif dim == 2
		vis=GridVisualizer(layout=(2,2), resolution=(660,660))
		scalarplot!(vis[1,1], mygrid, ps)
		scalarplot!(vis[1,2], mygrid, rho)
		nf = nodeflux(sys, sol)
		massflux = nf[:,ip,:]
		scalarplot!(vis[2,1], mygrid, massflux[1,:]./rho,)
		scalarplot!(vis[2,2], mygrid, massflux[2,:]./rho,)
	else
		vis=GridVisualizer(layout=(1,2), resolution=(660,660), outlinealpha=0.0)
		scalarplot!(vis[1,1], mygrid, ps, )
		scalarplot!(vis[1,2], mygrid, rho, )
		
	end
	reveal(vis)
end

# ╔═╡ 1224970e-8a59-48a9-b0ef-76ed776ca15d
function checkinout(sys,sol)	
	tfact=TestFunctionFactory(sys)
	tf_in=testfunction(tfact,[Γ_right],[Γ_left])
	tf_out=testfunction(tfact,[Γ_left],[Γ_right])
	(;in=integrate(sys,tf_in,sol),out=integrate(sys,tf_out,sol) )
end

# ╔═╡ e29848dd-d787-438e-9c32-e9c2136aec4f
checkinout(sys,sol)

# ╔═╡ dac0610e-63d8-409e-961c-caa160d38cca
let
	if RunSim
		(;ip,gn,gni,m,mfluxin,mmix0,X0,ng) = mydata

		in_,out_=checkinout(sys,sol)
	
		nout(i) = -out_[i]/m[i]
		nin(i) = mfluxin/mmix0 *bareas(Γ_left,sys,mygrid)*X0[i]
		nout_dry = 0.0
		RI=sum(integrate(sys,reaction,sol),dims=2) # reaction integral

		println("Total mass inflows and outflows:")
		@printf "IN: %2.6e \t OUT: %2.6e \t REACT: %2.6e kg/hr \n\n" in_[ip]/ufac"kg/hr" out_[ip]/ufac"kg/hr" RI[ip]/ufac"kg/hr"
		
		println("Molar species inflows, outflows and reaction integrals:")
		for i = 1:ng
			@printf "%i\tIN: %2.6e \t OUT: %2.6e \t REACT: %2.6e mol/hr \n" i nin(i)/ufac"mol/hr" nout(i)/ufac"mol/hr" -RI[i]/m[i]/ufac"mol/hr"
			if i != gni[:H2O] 
				nout_dry += nout(i)
			end
		end
	end
end

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╟─2c4b16ab-8f60-467b-b608-2fea9fbc741c
# ╠═d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═83fa22fa-451d-4c30-a4b7-834974245996
# ╠═4dae4173-0363-40bc-a9ca-ce5b4d5224cd
# ╠═561e96e2-2d48-4eb6-bb9d-ae167a622aeb
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═4e05ab31-7729-4a4b-9c14-145118477715
# ╠═107a6fa3-60cb-43f0-8b21-50cd1eb5065a
# ╠═832f3c15-b75a-4afe-8cc5-75ff3b4704d6
# ╟─0fadb9d2-1ccf-4d44-b748-b76d911784ca
# ╟─b94513c2-c94e-4bcb-9342-47ea48fbfd14
# ╟─c886dd12-a90c-40ab-b9d0-32934c17baee
# ╠═5547d7ad-dd58-4b00-8238-6e1abb32874e
# ╟─3440d4d8-3e03-4ff3-93f1-9afd7aaf9c41
# ╠═3bb2deff-7816-4749-9f1e-c1e451372b1e
# ╠═4af1792c-572e-465c-84bf-b67dd6a7bc93
# ╠═5f88937b-5802-4a4e-81e2-82737514b9e4
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═5588790a-73d4-435d-950f-515ae2de923c
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╠═e29848dd-d787-438e-9c32-e9c2136aec4f
# ╟─dac0610e-63d8-409e-961c-caa160d38cca
# ╟─05949759-2bb9-475b-b2f4-900b32c30e00
# ╠═111b1b1f-51a5-4069-a365-a713c92b79f4
# ╟─d4aaf146-e551-444a-a2f1-e70b05104b53
# ╠═de69f808-2618-4add-b092-522a1d7e0bb7
# ╠═389a4798-a9ee-4e9c-8b44-a06201b4c457
# ╠═1f2d79f6-5511-4286-907f-a4cf8c547f8e
# ╠═1224970e-8a59-48a9-b0ef-76ed776ca15d
