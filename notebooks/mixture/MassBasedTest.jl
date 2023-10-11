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
	using StaticArrays

	using LessUnitful
	
	using PlutoVista, Plots
	using PlutoUI, Colors

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ d3278ac7-db94-4119-8efd-4dd18107e248
PlutoUI.TableOfContents(title="M-S Transport + Darcy")

# ╔═╡ 83fa22fa-451d-4c30-a4b7-834974245996
function grid1D()
	X=(0:0.01:1)*ufac"cm"
	grid=simplexgrid(X)
	# catalyst region
	cellmask!(grid,[0.8]*ufac"cm",[0.9]*ufac"cm",2)	
	grid
end

# ╔═╡ b8819e9f-5262-421a-90ba-74f08d7fc52f
function grid2D()
	X=(0:0.01:1)*ufac"cm"
	Y=(0:0.01:1)*ufac"cm"
	grid=simplexgrid(X,Y)

	# catalyst region
	#cellmask!(grid,[0.7,0.0],[0.9,1.0],2)
	#bfacemask!(grid, [1,0.2],[1,0.8],5)
	
	grid
end

# ╔═╡ a995f83c-6ff7-4b95-a798-ea636ccb1d88
gridplot(grid1D(), resolution=(600,200))
#gridplot(grid2D())

# ╔═╡ 832f3c15-b75a-4afe-8cc5-75ff3b4704d6
begin
	# for 1D domain
	const Γ_left = 1
	const Γ_right = 2
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
	\frac{\partial \rho}{\partial t} - \nabla \cdot \left ( \rho \vec v \right)  &= 0\\
	\vec v  &= -\frac{\kappa}{\mu} \vec \nabla p\\
\end{align}
```
"""

# ╔═╡ 220e0a21-9328-4cf3-86cc-58468bb02cf7
md"""
!!! __check signs in flux function: diffusive and convective flux, respecting VoronoiFVM.jl convention__ !!!
"""

# ╔═╡ b94513c2-c94e-4bcb-9342-47ea48fbfd14
md"""
## Species Mass Transport
```math
\begin{align}
	\frac{\partial \rho_i}{\partial t} - \nabla \cdot \left( \vec \Phi_i - \rho_i \vec v \right ) + R_i &= 0 ~, \qquad i = 1 ... \nu \\
		-\frac{p}{RT}\frac{1}{M_{\text{mix}}} \left( \nabla x_i + (x_i-w_i) \frac{\nabla p}{p} \right) &= \sum_{j=1 \atop j \neq i}^{\nu} \frac{w_j \vec \Phi_i-w_i \vec \Phi_j}{D_{ij} M_i M_j} \\
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
	(;kp,km,p,m0,mf,ip,T, isreactive)=data
	if node.region == 2 && isreactive == 1 # catalyst layer		

		m = @. m0 + (mf-m0)*node.embedparam
		
		# 1 u2 -> 3 u1 
		c = u[ip]/(ph"R"*T)
		c1 = u[1]*c
		c2 = u[2]*c
		c3 = u[3]*c
		
		r = kp* c2
		r *=embedparam(node)

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
	boundary_neumann!(f,u,bnode, species=ip, region=Γ_left, value=mfluxin*embedparam(bnode))


	boundary_dirichlet!(f,u,bnode, species=ip,region=Γ_right,value=p)

end

# ╔═╡ db77fca9-4118-4825-b023-262d4073b2dd
md"""
### Peclet Number
```math
\text{Pe}_L= \frac{L \vec v}{D}

```
"""

# ╔═╡ 83a08a0f-eef6-436e-8cfa-a9af45fc5bb0
md"""
!!! __ask JF why Newton's method is not converging for larger Peclet numbers__ !!!

Idea: the solution does not change for increasing mass fluxes, thus there is no meaningful update step / near singular Jacobi Matrix?
"""

# ╔═╡ f40a8111-b5cb-40f0-8b12-f57cf59637f1
Base.@kwdef mutable struct ModelData
	ng::Int64 = 3
	ip::Int64 = ng+1
	p::Float64 = 1.0*ufac"bar"
	T::Float64 = 273.15*ufac"K"
	mf::Array{Float64,1} = [2.0,6.0,21.0]*ufac"g/mol"
	m0::Array{Float64,1} = [10.0,10.0,10.0]*ufac"g/mol"
	perm::Float64=1.23e-13*ufac"m^2" # perm. of porous medium, use in Darcy Eq.

	X0::Array{Float64,1} = [0.5, 0.3, 0.2]
	mmix0::Float64 = sum(X0 .* mf)
	W0::Array{Float64,1} = @. mf*X0/mmix0
	
	mfluxin::Float64=0.0001*ufac"kg/(m^2*s)"

	#isreactive::Int64 = 0
	isreactive::Int64 = 1
	kp::Float64=5.0
	km::Float64=1.0e-3
end

# ╔═╡ e183587f-db01-4f04-85fc-8d1a21bc0526
function darcyvelo(u,data)
	(;ip,perm) = data

	μ = 2.0e-5*ufac"Pa*s"
	perm/μ*(u[ip,1]-u[ip,2])	
end

# ╔═╡ 389a4798-a9ee-4e9c-8b44-a06201b4c457
function boutflow(f,u,edge,data)
	(;T,ip,m0,mf,ng)=data
	
	m = @. m0 + (mf-m0)*edge.embedparam

	k=outflownode(edge)

	pout = u[ip,k]
	cout = pout/(ph"R"*T)
	
	for i=1:(ng-1)
	#for i=1:ng
		# specify flux at boundary
		
		f[i] = -darcyvelo(u,data) * cout*u[i,k]*m[i]
	end	
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

# ╔═╡ 5547d7ad-dd58-4b00-8238-6e1abb32874e
function flux(f,u,edge,data)
	(;m0,mf,ng,ip,T)=data
		
	F = MVector{ng-1,eltype(u)}(undef)
	#X = MVector{ng,eltype(u)}(undef)
	M = MMatrix{ng-1,ng-1,eltype(u)}(undef)

	m = @. m0 + (mf-m0)*edge.embedparam
	
	D = DiffCoeffsMass(ng,m)
	
	pm = 0.5*(u[ip,1]+u[ip,2])
	c = pm/(ph"R"*T)

	
	x1 = 0.5*(u[1,1]+u[1,2])
	x2 = 0.5*(u[2,1]+u[2,2])
	x3 = 0.5*(u[3,1]+u[3,2])
	mmix = x1*m[1]+x2*m[2]+x3*m[3]

	rho = c*mmix
	v = darcyvelo(u,data)
	# compute total mass flux
	@inline f[ip] =rho*v
	
	w1 = m[1]*x1/mmix
	w2 = m[2]*x2/mmix
	w3 = m[3]*x3/mmix

	
	M[1,1] = -w2/D[1,2] - (w1+w3)/D[1,3]
	M[1,2] = w1/D[1,2] - w1/D[1,3]
	M[2,1] = w2*(1/D[2,1] - 1/D[2,3])
	M[2,2] = -w1/D[2,1] - (w2+w3)/D[2,3]

	#δx1 = u[1,1]-u[1,2]
	#δx2 = u[2,1]-u[2,2]
	#δp = u[ip,1]-u[ip,2]

	#F[1] = ( δx1 + (x1-w1)*δp/pm )*c/mmix
	#F[2] = ( δx2 + (x2-w2)*δp/pm )*c/mmix
	F[1] = (u[1,1]-u[1,2])*c/mmix
	F[2] = (u[2,1]-u[2,2])*c/mmix

	@inline inplace_linsolve!(M,F)

	#!!! check sign convention
	f[1] = -F[1] + rho*w1*v
	f[2] = -F[2] + rho*w2*v
	#!!! check sign convention
	

end

# ╔═╡ 480e4754-c97a-42af-805d-4eac871f4919
begin
	mygrid=grid1D()
	mydata=ModelData()
	(;p,ip)=mydata
	
	sys=VoronoiFVM.System( 	mygrid;
							data=mydata,
							flux=flux,
							reaction=reaction,
							bcondition=bcond,
							boutflow=boutflow,
							outflowboundaries=[Γ_right]
							)
	
	enable_species!(sys; species=collect(1:(3+1))) # gas phase species pi & ptotal	
	
	inival=unknowns(sys)

	inival[ip,:].=p
	for i=1:3
		inival[i,:] .= 1.0/3
	end

	control = SolverControl(nothing, sys;)
		control.Δp=1.0e-1
		#control.Δp_grow=1.2
		control.handle_exceptions=true
		control.Δu_opt=1.0e5
		control.Δp_min=1.0e-5
		#control.maxiters=500

	embed=[0,1]
	solt=solve(sys;inival=inival,embed,control,verbose="ne")	
end;

# ╔═╡ f798e27a-1d7f-40d0-9a36-e8f0f26899b6
@bind t Slider(solt.t,show_value=true,default=solt.t[end])

# ╔═╡ 5588790a-73d4-435d-950f-515ae2de923c
sol = solt(t);

# ╔═╡ 111b1b1f-51a5-4069-a365-a713c92b79f4
let
	(;ip,p) = ModelData()

	vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300))
	scalarplot!(vis, mygrid, sol[1,:], clear=false, label="1")
	scalarplot!(vis, mygrid, sol[2,:], clear=false, color=:red, label="2")
	scalarplot!(vis, mygrid, sol[3,:], clear=false, color=:blue, label="3")
	#scalarplot!(vis, mygrid, sol[ip,:] ./p, clear=false, color=:gray, label="pt/pout")
	reveal(vis)
end

# ╔═╡ de69f808-2618-4add-b092-522a1d7e0bb7
let
	(;p,mf,ng,ip,T, mfluxin) = ModelData()
	mmix = []
	for j in 1:length(sol[1,:])
		_mmix=0
		for i=1:ng
			_mmix += sol[i,j]*mf[i]
		end
		push!(mmix, _mmix)
	end
	
	w1 = sol[1,:]*mf[1] ./mmix
	w2 = sol[2,:]*mf[2] ./mmix
	w3 = sol[3,:]*mf[3] ./mmix
	
	vis=GridVisualizer(legend=:lt, resolution=(600,400), layout = (2, 1))
	#scalarplot!(vis, mygrid, mmix, clear=false, label="mmix")
	scalarplot!(vis[1,1], mygrid, w1, clear=false, label="w1")
	scalarplot!(vis[1,1], mygrid, w2, clear=false, color=:red, label="w2")
	scalarplot!(vis[1,1], mygrid, w3, clear=false, color=:blue, label="w3")

	p0 = sol[ip,1]
	rho0 = @. sol[ip,1] * mmix[1] /(ph"R"*T)
	rho = @. sol[ip,:] * mmix /(ph"R"*T)
	scalarplot!(vis[2,1], mygrid, rho/rho0, clear=false, label="Rho / Rho0")
	scalarplot!(vis[2,1], mygrid, rho0./rho, clear=false, color=:red, label="v / v0")
	scalarplot!(vis[2,1], mygrid, sol[ip,:]/p0, clear=false, color=:blue, label="p / p0")
	reveal(vis)
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
	(;mfluxin,mmix0,p,T,ng,mf) = mydata
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
	(;kp,km,p,m0,mf,ip,T)=data
	if node.region == 2 # catalyst layer
		

		m = @. m0 + (mf-m0)
		
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

# ╔═╡ Cell order:
# ╠═c21e1942-628c-11ee-2434-fd4adbdd2b93
# ╠═d3278ac7-db94-4119-8efd-4dd18107e248
# ╠═83fa22fa-451d-4c30-a4b7-834974245996
# ╠═b8819e9f-5262-421a-90ba-74f08d7fc52f
# ╠═a995f83c-6ff7-4b95-a798-ea636ccb1d88
# ╠═832f3c15-b75a-4afe-8cc5-75ff3b4704d6
# ╟─a078e1e1-c9cd-4d34-86d9-df4a052b6b96
# ╟─0fadb9d2-1ccf-4d44-b748-b76d911784ca
# ╟─220e0a21-9328-4cf3-86cc-58468bb02cf7
# ╟─b94513c2-c94e-4bcb-9342-47ea48fbfd14
# ╟─c886dd12-a90c-40ab-b9d0-32934c17baee
# ╠═5547d7ad-dd58-4b00-8238-6e1abb32874e
# ╠═3bb2deff-7816-4749-9f1e-c1e451372b1e
# ╠═5f88937b-5802-4a4e-81e2-82737514b9e4
# ╠═480e4754-c97a-42af-805d-4eac871f4919
# ╠═5588790a-73d4-435d-950f-515ae2de923c
# ╠═f798e27a-1d7f-40d0-9a36-e8f0f26899b6
# ╠═111b1b1f-51a5-4069-a365-a713c92b79f4
# ╠═de69f808-2618-4add-b092-522a1d7e0bb7
# ╟─db77fca9-4118-4825-b023-262d4073b2dd
# ╟─83a08a0f-eef6-436e-8cfa-a9af45fc5bb0
# ╠═ae8c7993-a89f-438a-a72a-d4a0c9a8ce57
# ╠═e29848dd-d787-438e-9c32-e9c2136aec4f
# ╠═f40a8111-b5cb-40f0-8b12-f57cf59637f1
# ╠═e183587f-db01-4f04-85fc-8d1a21bc0526
# ╠═389a4798-a9ee-4e9c-8b44-a06201b4c457
# ╠═37b5908c-dd4e-4fb8-9d5b-68402493e10d
# ╠═2cbd3e87-289c-47a2-b837-10133974ae82
# ╠═1224970e-8a59-48a9-b0ef-76ed776ca15d
# ╠═05a8556e-8bf4-442c-ad79-171efed51e85
