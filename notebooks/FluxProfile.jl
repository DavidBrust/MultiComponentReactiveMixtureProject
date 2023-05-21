### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 37edeff2-f335-11ed-0a1a-67aa587d3a23
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using LessUnitful
	using CSV
	using DataFrames
	using Statistics

	using VoronoiFVM
	using GridVisualize
	using ExtendableGrids
	using PlutoVista, Plots
	using PlutoUI
	using Interpolations

	using FixedBed

	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 8e4d183d-fc6b-47b5-9c6a-3b18b14c4297
begin
	FluxMap = CSV.read("../data/IrradiationFluxProfiles/IrradFlux.csv", DataFrame, delim=";")
	M=Matrix(FluxMap)
end;

# ╔═╡ c2c77844-c43c-4580-877f-62c238fc68bb
md"""
# Coordinate System
X-Y coordinates, starting from the upper left corner (0,0)

Δx = 0.06579 cm

Δy = 0.06579 cm
"""

# ╔═╡ befaf37d-d1a1-49ae-8a15-6a86e2807e6b
begin
	Dx = 0.06579*ufac"cm"
	Dy = 0.06579*ufac"cm"
end;

# ╔═╡ 6766a4cb-296a-4bbf-8577-a24674243dda
let
	p=Plots.plot(size=(400,400))
	Plots.contourf!(p,M, aspect_ratio=1, levels=maximum(M).*[0.0,0.1,0.3,0.5,0.7,0.9,1.0], lw=0, color=[:black,:blue,:cyan,:lawngreen,:yellow,:orange])
	
end

# ╔═╡ 872522e3-b7e8-4a96-8ce2-2bf3458ca32a
md"""
We only need cut outs of the total flux profile
- 10 cm x 10 cm selection (target: as flat as possible)
- 12 cm x 12 cm selection (entering the aperture)
"""

# ╔═╡ 4379a876-13da-4b1d-a3b5-ec11b1f25a27
md"""
# 10 cm x 10 cm
"""

# ╔═╡ 1c953043-6424-4813-b88b-f8665bd1d317
md"""
From optimization w.r.t. "flatness" of the profile, for 10 cm x 10 cm the selection begins at __col 82__ and __row 80__
"""

# ╔═╡ 60516a56-e663-44d3-ac8e-21b3064a2ab7
function sel10by10(M)
	sr=80
	sc=82

	nx=Integer(round(10.0*ufac"cm"/Dx))
	ny=Integer(round(10.0*ufac"cm"/Dy))
	
	M=Matrix(FluxMap)
	@views M_ = M[sr:(sr+ny),sc:(sc+nx)]

	nx,ny,M_
end

# ╔═╡ 777c6089-7a78-4dfb-8017-cb3bc69f8ab3
begin
	nx10,ny10,M10 = sel10by10(M)
	
	x10 = range(-5.0,5.0,length=nx10+1)
	y10 = range(-5.0,5.0,length=ny10+1)

	itp10 = Interpolations.interpolate((x10,y10), M10, Gridded(Linear()))
end

# ╔═╡ a5a0dd80-0171-4540-9b86-b6be242a733d
md"""
- Mean: __$(round(mean(M10),digits=1))__ $\,\text{kW}\,\text{m}^{-2}$
- Min: __$(round(minimum(M10),digits=1))__ $\,\text{kW}\,\text{m}^{-2}$
- Max: __$(round(maximum(M10),digits=1))__ $\,\text{kW}\,\text{m}^{-2}$
- Std: __$(round(std(M10),digits=1))__ $\,\text{kW}\,\text{m}^{-2}$
- Power: __$(round(sum(M10)*Dx*Dy,digits=1))__ $\,\text{kW}$
"""

# ╔═╡ cd6dfa70-db90-4d38-b265-a89da0abb93c
let
	p=Plots.plot(layout=(1,2))
	Plots.contour!(p[1,1],M10, lw=0, aspect_ratio = 1, fill = true)

	x_ = range(-5.0,5.0,length=20)
	y_ = range(-5.0,5.0,length=20)

	Plots.contour!(p[1,2],x_,y_,itp10(x_,y_), lw=0, aspect_ratio = 1, fill = true)
	p
end

# ╔═╡ 832c53d5-c004-4f71-8947-4146529442f7
md"""
# 12 cm x 12 cm
"""

# ╔═╡ 19951420-f1c2-47ea-9be0-a5db4ceebac0
function sel12by12(M)
	# starting coordinates for optimum 10 cm x 10 cm selection
	sr=80
	sc=82

	# (inverse) resolution of flux measurements, distance between data points
	Dx = 0.06579*ufac"cm"
	Dy = 0.06579*ufac"cm"

	# coordinate offsets, when considering 12 cm x 12 cm selection
	Dsr=Integer(round(1.0*ufac"cm"/Dx))
	Dsc=Integer(round(1.0*ufac"cm"/Dy))

	sr-=Dsr
	sc-=Dsc

	nx=Integer(round(12.0*ufac"cm"/Dx))
	ny=Integer(round(12.0*ufac"cm"/Dy))
	
	M=Matrix(FluxMap)
	@views M_ = M[sr:(sr+ny),sc:(sc+nx)]

	wi=6.0*ufac"cm"
	x = range(-wi,wi,length=nx+1)
	y = range(-wi,wi,length=ny+1)

	itp = Interpolations.interpolate((x,y), M_, Gridded(Linear()))

	M_,itp	
end

# ╔═╡ bfcdd56d-8421-439e-89cb-a84ca6e659e8
md"""
# Irradiation FLux Profile
"""

# ╔═╡ b10b6dcd-cddd-47f0-b6c8-00f79eb7b3fa
let
	M12, itp12 = sel12by12(M)

	wi=12.0*ufac"cm"
	nx,ny=size(M12)
	x_ = range(-wi/2,wi/2,length=nx)
	y_ = range(-wi/2,wi/2,length=ny)
	
	p=Plots.plot(xguide="X Coordinate / cm", yguide="Y Coordinate / cm", colorbar_title="Irradiation Flux / kW m-2", xlim=(-wi/2/ufac"cm", wi/2/ufac"cm"), ylim=(-wi/2/ufac"cm", wi/2/ufac"cm"))
	
	p1=Plots.contour!(p,x_./ufac"cm",y_./ufac"cm",M12, lw=0, aspect_ratio = 1, fill = true)	

	
	x_ = range(-wi/2,wi/2,length=20)
	y_ = range(-wi/2,wi/2,length=20)

	p=Plots.plot(xguide="X Coordinate / cm", yguide="Y Coordinate / cm", colorbar_title="Irradiation Flux / kW m-2", xlim=(-wi/2/ufac"cm", wi/2/ufac"cm"), ylim=(-wi/2/ufac"cm", wi/2/ufac"cm"))
	
	p2=Plots.contour!(p,x_/ufac"cm",y_/ufac"cm",itp12(x_,y_), lw=0, aspect_ratio = 1, fill = true)
	Plots.plot(p1,p2)
end

# ╔═╡ Cell order:
# ╠═37edeff2-f335-11ed-0a1a-67aa587d3a23
# ╠═8e4d183d-fc6b-47b5-9c6a-3b18b14c4297
# ╟─c2c77844-c43c-4580-877f-62c238fc68bb
# ╠═befaf37d-d1a1-49ae-8a15-6a86e2807e6b
# ╠═6766a4cb-296a-4bbf-8577-a24674243dda
# ╠═872522e3-b7e8-4a96-8ce2-2bf3458ca32a
# ╟─4379a876-13da-4b1d-a3b5-ec11b1f25a27
# ╟─1c953043-6424-4813-b88b-f8665bd1d317
# ╟─a5a0dd80-0171-4540-9b86-b6be242a733d
# ╠═60516a56-e663-44d3-ac8e-21b3064a2ab7
# ╠═cd6dfa70-db90-4d38-b265-a89da0abb93c
# ╠═777c6089-7a78-4dfb-8017-cb3bc69f8ab3
# ╟─832c53d5-c004-4f71-8947-4146529442f7
# ╠═19951420-f1c2-47ea-9be0-a5db4ceebac0
# ╟─bfcdd56d-8421-439e-89cb-a84ca6e659e8
# ╠═b10b6dcd-cddd-47f0-b6c8-00f79eb7b3fa
