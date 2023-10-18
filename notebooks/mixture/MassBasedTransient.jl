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
	Pkg.activate(joinpath(@__DIR__,"..\\.."))
	using Revise
	using VoronoiFVM
	using ExtendableGrids, GridVisualize,ExtendableSparse,SparseArrays
	using NLsolve, LinearSolve
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots, PyPlot
	using PlutoUI, Colors

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 6da83dc0-3b0c-4737-833c-6ee91552ff5c
md"""
Check the box to start the simulation:

__Run Sim__ $(@bind RunSim PlutoUI.CheckBox(default=true))
"""

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
PlutoUI.TableOfContents(title="M-S Transport + Darcy")

# ╔═╡ 83fa22fa-451d-4c30-a4b7-834974245996
function grid1D()
	X=(0:0.05:1)*ufac"cm"
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
	bfacemask!(grid, [1,0.2].*ufac"cm",[1,0.8].*ufac"cm",5)
	
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
const dim = 2

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
		#const Γ_right = 2
		const Γ_top = 3
		const Γ_left = 4
		const Γ_right = 5
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

# ╔═╡ a078e1e1-c9cd-4d34-86d9-df4a052b6b96
md"""
Test notebook with 3 species.
"""

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

# ╔═╡ 220e0a21-9328-4cf3-86cc-58468bb02cf7
md"""
Care must be taken with respect to the sign convention used in VoronoiFVM.jl.
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

# ╔═╡ 3bb2deff-7816-4749-9f1e-c1e451372b1e
function reaction(f,u,node,data)
	(;kp,km,p,m,ip,T, isreactive)=data
	if node.region == 2 && isreactive == 1 # catalyst layer		

		
		# 1 u2 -> 3 u1 
		c = u[ip]/(ph"R"*T)
		c1 = u[1]*c
		c2 = u[2]*c
		c3 = u[3]*c
		
		r = kp* c2

		f[1] = -3*r*m[1]
		f[2] = r*m[2]
		#f[3] = r
		
	end
	
	for i=1:3
		f[3] += u[i]
	end
	f[3] = f[3] - 1.0
end

# ╔═╡ 5f88937b-5802-4a4e-81e2-82737514b9e4
function bcond(f,u,bnode,data)
	(;p,ip,mfluxin,X0)=data

	
	boundary_dirichlet!(f,u,bnode, species=1,region=Γ_left,value=X0[1])
	boundary_dirichlet!(f,u,bnode, species=2,region=Γ_left,value=X0[2])
	
	#boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_left,value=p)
	#boundary_neumann!(f,u,bnode, species=ip, region=Γ_left, value=mfluxin*embedparam(bnode))
	boundary_neumann!(f,u,bnode, species=ip, region=Γ_left, value=mfluxin)


	boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_right,value=p)

end

# ╔═╡ db77fca9-4118-4825-b023-262d4073b2dd
md"""
### Peclet Number
```math
\text{Pe}_L= \frac{L \vec v}{D}

```
"""

# ╔═╡ e183587f-db01-4f04-85fc-8d1a21bc0526
function darcyvelo(u,data)
	(;ip,perm) = data

	μ = 2.0e-5*ufac"Pa*s"
	-perm/μ*(u[ip,1]-u[ip,2])	
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

# ╔═╡ f40a8111-b5cb-40f0-8b12-f57cf59637f1
begin
	
Base.@kwdef mutable struct ModelData{NG}
	#ng::Int64 = 3
	ip::Int64 = NG+1
	p::Float64 = 1.0*ufac"bar"
	T::Float64 = 273.15*ufac"K"
	m::Vector{Float64} = [2.0,6.0,21.0]*ufac"g/mol"
	perm::Float64=1.23e-11*ufac"m^2" # perm. of porous medium, use in Darcy Eq.

	X0::Vector{Float64} = [0.1, 0.7, 0.2]

	Dcoeff::Matrix{Float64} = DiffCoeffsMass(NG,m)
	
	mmix0::Float64 = sum(X0 .* m)
	W0::Vector{Float64} = @. m*X0/mmix0
	
	mfluxin::Float64=0.01*ufac"kg/(m^2*s)"


	#isreactive::Int64 = 0
	isreactive::Int64 = 1
	kp::Float64=5.0e-0
	km::Float64=1.0e-3
end

ModelData(;ng=3, kwargs...) = ModelData{ng}(;kwargs...)

ngas(::ModelData{NG}) where NG = NG
end

# ╔═╡ 5547d7ad-dd58-4b00-8238-6e1abb32874e
function flux(f,u,edge,data)
	(;m,ip,T,Dcoeff)=data
	ng=ngas(data)
		
	F = MVector{ng-1,eltype(u)}(undef)
	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)
	#D = MMatrix{ng,ng,eltype(u)}(undef)

	#@inline D_matrix!(data, D)
		
	pm = 0.5*(u[ip,1]+u[ip,2])
	c = pm/(ph"R"*T)

	
	x1 = 0.5*(u[1,1]+u[1,2])
	x2 = 0.5*(u[2,1]+u[2,2])
	x3 = 0.5*(u[3,1]+u[3,2])
	mmix = x1*m[1]+x2*m[2]+x3*m[3]

	rho = c*mmix
	v = darcyvelo(u,data)
	
	f[ip] = -rho*v
	
	w1 = m[1]*x1/mmix
	w2 = m[2]*x2/mmix
	w3 = m[3]*x3/mmix

	M[1,1] = -w2/Dcoeff[1,2] - (w1+w3)/Dcoeff[1,3]
	M[1,2] = w1/Dcoeff[1,2] - w1/Dcoeff[1,3]
	M[2,1] = w2*(1/Dcoeff[2,1] - 1/Dcoeff[2,3])
	M[2,2] = -w1/Dcoeff[2,1] - (w2+w3)/Dcoeff[2,3]
	
	#M[1,1] = -w2/D[1,2] - (w1+w3)/D[1,3]
	#M[1,2] = w1/D[1,2] - w1/D[1,3]
	#M[2,1] = w2*(1/D[2,1] - 1/D[2,3])
	#M[2,2] = -w1/D[2,1] - (w2+w3)/D[2,3]

	##δx1 = u[1,1]-u[1,2]
	##δx2 = u[2,1]-u[2,2]
	##δp = u[ip,1]-u[ip,2]

	##F[1] = ( δx1 + (x1-w1)*δp/pm )*c/mmix
	##F[2] = ( δx2 + (x2-w2)*δp/pm )*c/mmix
	F[1] = (u[1,1]-u[1,2])*c/mmix
	F[2] = (u[2,1]-u[2,2])*c/mmix

	@inline inplace_linsolve!(M,F)

	## central difference flux
	##f[1] = -(F[1] + rho*w1*v)
	##f[2] = -(F[2] + rho*w2*v)
	f[1] = -(F[1] + c*x1*m[1]*v)
	f[2] = -(F[2] + c*x2*m[2]*v)

end

# ╔═╡ 4af1792c-572e-465c-84bf-b67dd6a7bc93
function storage(f,u,node,data)
	(;T,ip,m)=data
	ng=ngas(data)

	c = u[ip]/(ph"R"*T)
	for i=1:ng
		f[i]=c*u[i]*m[i]
	end
	
	# total pressure
	mmix = u[1]*m[1]+u[2]*m[2]+u[3]*m[3]
	f[ng+1] = mmix*c
end

# ╔═╡ 389a4798-a9ee-4e9c-8b44-a06201b4c457
function boutflow(f,u,edge,data)
	(;T,ip,m)=data
	ng=ngas(data)

	k=outflownode(edge)

	pout = u[ip,k]
	cout = pout/(ph"R"*T)
	
	for i=1:(ng-1)
		# specify flux at boundary
		
		#f[i] = -darcyvelo(u,data) * cout*u[i,k]*m[i]
		f[i] = darcyvelo(u,data) * cout*u[i,k]*m[i]
	end	
end

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
begin
	
	if dim == 1
		mygrid=grid1D()
		strategy = nothing
		times=[0,20]
	elseif dim == 2
		mygrid=grid2D()
		strategy = nothing
		times=[0,10]
	else
		mygrid=grid3D()
		strategy = GMRESIteration(UMFPACKFactorization())
		times=[0,1]
	end
	mydata=ModelData()
	(;p,ip)=mydata
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
	#control = SolverControl(GMRESIteration(UMFPACKFactorization()), sys;)
		control.Δt=1.0e-4
		control.handle_exceptions=true
		control.Δu_opt=1.0e5

	if RunSim
		#times=[0,5.0]
		solt=solve(sys;inival=inival,times,control,verbose="ae")
	end
end;

# ╔═╡ f798e27a-1d7f-40d0-9a36-e8f0f26899b6
@bind t Slider(solt.t,show_value=true,default=solt.t[end])

# ╔═╡ 5588790a-73d4-435d-950f-515ae2de923c
sol = solt(t);

# ╔═╡ 111b1b1f-51a5-4069-a365-a713c92b79f4
let
	(;ip,p) = ModelData()
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
	mydata = ModelData()
	(;p,m,ip,T, mfluxin) = mydata
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
	rho = @. ps * mmix /(ph"R"*T)
	
	if dim == 1
		vis=GridVisualizer(legend=:lt, resolution=(600,400), layout = (2, 1))
		scalarplot!(vis[1,1], mygrid, w1, clear=false, label="w1")
		scalarplot!(vis[1,1], mygrid, w2, clear=false, color=:red, label="w2")
		scalarplot!(vis[1,1], mygrid, w3, clear=false, color=:blue, label="w3")

		p0 = sol[ip,1]
		rho0 = @. p0 * mmix[1] /(ph"R"*T)
		scalarplot!(vis[2,1], mygrid, rho/rho0, clear=false, label="Rho / Rho0")
		scalarplot!(vis[2,1], mygrid, rho0./rho, clear=false, color=:red, label="v / v0")
		scalarplot!(vis[2,1], mygrid, ps/p0, clear=false, color=:blue, label="p / p0")
	elseif dim == 2
		vis=GridVisualizer(layout=(2,2), resolution=(660,660))
		scalarplot!(vis[1,1], mygrid, ps, title="Total Pressure")
		scalarplot!(vis[1,2], mygrid, rho, title="Total Density")
		nf = nodeflux(sys, sol)
		massflux = nf[:,ip,:]
		scalarplot!(vis[2,1], mygrid, massflux[1,:]./rho, title="Velocity - X")
		scalarplot!(vis[2,2], mygrid, massflux[2,:]./rho, title="Velocity - Y")
	else
		vis=GridVisualizer(layout=(2,1), resolution=(400,800), outlinealpha=0.0)
		scalarplot!(vis[1,1], mygrid, ps, title="Total Pressure")
		scalarplot!(vis[2,1], mygrid, rho, title="Total Density")
		
	end
	reveal(vis)
end

# ╔═╡ ca08d9a2-8148-441e-a149-b9d0c5232a6d
function D_matrix!(data, D)
	(;m,Dcoeff)=data
	ng=ngas(data)
	
	for i=1:(ng-1)
		for j=(i+1):ng
			@views D[j,i] = one(eltype(D))
			@views D[i,j] = one(eltype(D))
		end
	end
end

# ╔═╡ 2cbd3e87-289c-47a2-b837-10133974ae82
function DiffCoeffs(ng)
	k=0
	D=zeros(Float64, ng,ng)
	for i=1:(ng-1)
		for j=(i+1):ng
			k +=1
			Dji = k*1.0e-5*ufac"m^2/s"
			D[j,i] = Dji
			D[i,j] = Dji
		end			
	end
	D
end

# ╔═╡ ae8c7993-a89f-438a-a72a-d4a0c9a8ce57
let
	L=mygrid[Coordinates][end]*ufac"m"
	(;mfluxin,mmix0,p,T,m) = mydata
	ng = ngas(mydata)
	rho0 = p*mmix0/(ph"R"*T)
	v0 = mfluxin / rho0
	D = DiffCoeffs(ng)
	Pe = L*v0 / minimum(D[D.>0])
end

# ╔═╡ 1224970e-8a59-48a9-b0ef-76ed776ca15d
function checkinout(sys,sol)
	
	tfact=TestFunctionFactory(sys)
	tf_in=testfunction(tfact,[Γ_right],[Γ_left])
	tf_out=testfunction(tfact,[Γ_left],[Γ_right])
	#tf_in=testfunction(tfact,[Γ_left],[Γ_right])
	#tf_out=testfunction(tfact,[Γ_right],[Γ_left])
	(;in=integrate(sys,tf_in,sol),out=integrate(sys,tf_out,sol) )
end

# ╔═╡ e29848dd-d787-438e-9c32-e9c2136aec4f
checkinout(sys,sol)

# ╔═╡ 05a8556e-8bf4-442c-ad79-171efed51e85
function check_reaction(f,u,node,data)
	(;kp,km,p,m,ip,T)=data
	if node.region == 2 # catalyst layer
		
		# 1 u2 -> 3 u1 
		c = u[ip]/(ph"R"*T)
		c1 = u[1]*c
		c2 = u[2]*c
		c3 = u[3]*c
		
		r = kp* c2

		f[1] = -3*r*m[1]
		f[2] = r*m[2]
		#f[3] = r
	end
	
	for i=1:3
		f[3] += u[i]
	end
	f[3] = f[3] - 1.0
end

# ╔═╡ 370ebbae-9caf-4e56-ad77-020048bf7aa8
integrate(sys,check_reaction,sol)

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
# ╟─220e0a21-9328-4cf3-86cc-58468bb02cf7
# ╟─b94513c2-c94e-4bcb-9342-47ea48fbfd14
# ╟─c886dd12-a90c-40ab-b9d0-32934c17baee
# ╠═5547d7ad-dd58-4b00-8238-6e1abb32874e
# ╠═3bb2deff-7816-4749-9f1e-c1e451372b1e
# ╠═4af1792c-572e-465c-84bf-b67dd6a7bc93
# ╠═5f88937b-5802-4a4e-81e2-82737514b9e4
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═5588790a-73d4-435d-950f-515ae2de923c
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╠═e29848dd-d787-438e-9c32-e9c2136aec4f
# ╠═111b1b1f-51a5-4069-a365-a713c92b79f4
# ╠═de69f808-2618-4add-b092-522a1d7e0bb7
# ╠═370ebbae-9caf-4e56-ad77-020048bf7aa8
# ╟─db77fca9-4118-4825-b023-262d4073b2dd
# ╠═ae8c7993-a89f-438a-a72a-d4a0c9a8ce57
# ╠═f40a8111-b5cb-40f0-8b12-f57cf59637f1
# ╠═e183587f-db01-4f04-85fc-8d1a21bc0526
# ╠═389a4798-a9ee-4e9c-8b44-a06201b4c457
# ╠═ca08d9a2-8148-441e-a149-b9d0c5232a6d
# ╠═37b5908c-dd4e-4fb8-9d5b-68402493e10d
# ╠═2cbd3e87-289c-47a2-b837-10133974ae82
# ╠═1224970e-8a59-48a9-b0ef-76ed776ca15d
# ╠═05a8556e-8bf4-442c-ad79-171efed51e85
