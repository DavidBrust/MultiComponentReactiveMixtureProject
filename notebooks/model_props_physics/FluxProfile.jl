### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 37edeff2-f335-11ed-0a1a-67aa587d3a23
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,"../.."))
	using Revise
	using LessUnitful
	using CSV, Printf
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

# ╔═╡ b4a07bfc-f30e-4970-b301-bed21f50e943
nom_flux = 40.0

# ╔═╡ 8e4d183d-fc6b-47b5-9c6a-3b18b14c4297
function readFlux(flux)
	path = "../../data/IrradiationFluxProfiles/"
	fluxes = [40.0,60.0,80.0,100.0]
	@assert flux in fluxes
	d_flux_n = Dict(
	40.0=>("FluxMap_Import_20230523_110038_40.csv","Target_coords_20230523_110038_40.csv"),
	60.0=>	("FluxMap_Import_20230523_105203_60.csv","Target_coords_20230523_105203_60.csv"),
	80.0=>("FluxMap_Import_20230523_105025_80.csv","Target_coords_20230523_105025_80.csv"),
	100.0=>("FluxMap_Import_20230523_104820_100.csv","Target_coords_20230523_104820_100.csv")
	)
	fn_flux, fn_coord = d_flux_n[flux]
	FluxMap = CSV.read(path*fn_flux, DataFrame, header=false,delim=";")
	coords = CSV.read(path*fn_coord, DataFrame, header=1,delim=";")
	M=Matrix(FluxMap)
	reverse!(M; dims=1)
	return M, coords
end;

# ╔═╡ 545ad3b6-c2b5-472f-87c4-38a56df9dc78
M, coords = readFlux(nom_flux);

# ╔═╡ 3e816f03-7cb2-4a54-80bd-4ed340c92149
function find_offsets(M,lb,ub)
	inds = 30:1:33
		
	for sc in inds
		for ec in inds
			for sr in inds
				for er in inds
					M_ = @views mean(M[sc:(end-ec),sr:(end-er)])
					if M_ > lb && M_ < ub
						return (sc=sc, ec=ec, sr=sr, er=er)
					end
				end
			end
		end
	end
end

# ╔═╡ 2689148e-2fd5-40cf-95c9-90339c0e1532
D_offsets = Dict(
	40.0=>(39.050003, 39.0500032),
	60.0=>(58.713004, 58.713005),
	80.0=>(79.024620, 79.024621),
	100.0=>(102.494731, 102.494732),
)

# ╔═╡ 2acc8661-8e67-4229-93ff-ebb82ae39d0b
offsets = find_offsets(M,D_offsets[nom_flux]...)

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
Determined row and column offset through brute force search for 10 cm x 10 cm and comparison with average value reported previously.

- Start column: $(offsets.sc)
- End column: $(offsets.ec)
- Start row: $(offsets.sr)
- End row: $(offsets.er)
"""

# ╔═╡ 1a51e389-6fdd-448d-a58f-48a7112af505
begin
	width12 =  coords.X[end] - coords.X[1]
	height12 =  coords.Y[end] - coords.X[1]
end;

# ╔═╡ caca277d-ecc9-47cc-9c88-1dc3d7ee3ce1
md"""
## 12 cm x 12 cm
- Mean: __$(round(mean(M),digits=2))__ $\,\text{kW}\,\text{m}^{-2}$
- Min: __$(round(minimum(M),digits=2))__ $\,\text{kW}\,\text{m}^{-2}$
- Max: __$(round(maximum(M),digits=2))__ $\,\text{kW}\,\text{m}^{-2}$
- Std: __$(round(std(M),digits=2))__ $\,\text{kW}\,\text{m}^{-2}$
- Power: __$(round(sum(M)*width12*height12*ufac"cm"^2/length(M),digits=2))__ $\,\text{kW}$
"""

# ╔═╡ 7031781a-67b6-4ab1-91e1-2e5d1b853068
begin
	(;sc,ec,sr,er) = offsets
	width10 =  coords.X[end-ec] - coords.X[sc]
	height10 =  coords.Y[end-er] - coords.X[sr]
	M10 = M[sc:end-ec,sr:end-er]
end;

# ╔═╡ a5a0dd80-0171-4540-9b86-b6be242a733d
md"""
## 10 cm x 10 cm
- Mean: __$(round(mean(M10),digits=2))__ $\,\text{kW}\,\text{m}^{-2}$
- Min: __$(round(minimum(M10),digits=2))__ $\,\text{kW}\,\text{m}^{-2}$
- Max: __$(round(maximum(M10),digits=2))__ $\,\text{kW}\,\text{m}^{-2}$
- Std: __$(round(std(M10),digits=2))__ $\,\text{kW}\,\text{m}^{-2}$
- Power: __$(round(sum(M10)*width10*height10*ufac"cm"^2/length(M10),digits=2))__ $\,\text{kW}$
"""

# ╔═╡ a10e7372-4abd-4afa-9e29-7bb34a994deb
D_levels = Dict(
	40.0 => [0.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0],
	60.0 => [0.0,10.0,20.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0],
	80.0 => [0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,75.0,80.0,120.0],
	100.0 => [0.0,25.0,50.0,60.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,150.0],
);

# ╔═╡ e2ea9dbb-c2af-439d-bb60-4432df9cd6e1
function plotBorder(p)
	Plots.plot!(p, [1.0;11.0], [1.0;1.0], c=:black, lw=2, label=:none)
	Plots.plot!(p, [1.0;1.0], [1.0;11.0], c=:black, lw=2, label=:none)
	Plots.plot!(p, [1.0;11.0], [11.0;11.0], c=:black, lw=2, label=:none)
	Plots.plot!(p, [11.0;11.0], [1.0;11.0], c=:black, lw=2, label=:none)
end

# ╔═╡ 6766a4cb-296a-4bbf-8577-a24674243dda
let
	p=Plots.plot(size=(400,400))
	Plots.contourf!(p,coords.X,coords.Y,M, aspect_ratio=1, levels=maximum(M).*[0.0,0.1,0.3,0.5,0.7,0.9,1.0], lw=0, color=[:black,:blue,:cyan,:lawngreen,:yellow,:orange])
	plotBorder(p)
end

# ╔═╡ cd6dfa70-db90-4d38-b265-a89da0abb93c
let
	(;sc,ec,sr,er) = offsets
	p = Plots.plot(xlabel="X / cm", ylabel="Y / cm", xlimits=(0,12.0), ylimits=(0,12.0) )

	Plots.contour!(p,coords.X,coords.Y,M, lw=0, aspect_ratio = 1, fill = true, levels=D_levels[nom_flux], clabels=true, colorbar_title="Flux density / kW/m²")
	plotBorder(p)
	#Plots.savefig(p, "./out/img/$(Integer(round(nom_flux)))kw_m2.svg")
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

# ╔═╡ 1d2daf9b-e0df-4e12-bd97-956f975e717b
function flux_interpol(flux)
	nom_flux = 0.0
	if flux >= 40.0 && flux < 50.0
		nom_flux = 40.0
	elseif flux >= 50.0 && flux < 70.0
		nom_flux = 60.0
	elseif flux >= 70.0 && flux < 90.0
		nom_flux = 80.0
	elseif flux >= 90.0 && flux <= 100.0
		nom_flux = 100.0
	end
	M, coords = readFlux(nom_flux);
	M .*= flux/nom_flux*ufac"kW/m^2"
	intp = 		Interpolations.linear_interpolation((coords.X*ufac"cm",coords.Y*ufac"cm"), M,extrapolation_bc=Flat())	
end

# ╔═╡ 76365437-d304-4a4d-be35-057d2ae39927
let
	xs = 0.0:0.01:0.12
	ys = xs
	xys = Vector{Float64}[]
	for x in xs
		for y in ys
			push!(xys, [x,y])
		end
	end
	intp = flux_interpol(70.0)
	intp(0.13,0.13)
end

# ╔═╡ Cell order:
# ╠═37edeff2-f335-11ed-0a1a-67aa587d3a23
# ╠═b4a07bfc-f30e-4970-b301-bed21f50e943
# ╟─a5a0dd80-0171-4540-9b86-b6be242a733d
# ╟─caca277d-ecc9-47cc-9c88-1dc3d7ee3ce1
# ╠═6766a4cb-296a-4bbf-8577-a24674243dda
# ╠═8e4d183d-fc6b-47b5-9c6a-3b18b14c4297
# ╠═545ad3b6-c2b5-472f-87c4-38a56df9dc78
# ╠═3e816f03-7cb2-4a54-80bd-4ed340c92149
# ╠═2689148e-2fd5-40cf-95c9-90339c0e1532
# ╠═2acc8661-8e67-4229-93ff-ebb82ae39d0b
# ╟─872522e3-b7e8-4a96-8ce2-2bf3458ca32a
# ╟─4379a876-13da-4b1d-a3b5-ec11b1f25a27
# ╟─1c953043-6424-4813-b88b-f8665bd1d317
# ╠═1a51e389-6fdd-448d-a58f-48a7112af505
# ╠═7031781a-67b6-4ab1-91e1-2e5d1b853068
# ╠═cd6dfa70-db90-4d38-b265-a89da0abb93c
# ╠═a10e7372-4abd-4afa-9e29-7bb34a994deb
# ╠═e2ea9dbb-c2af-439d-bb60-4432df9cd6e1
# ╟─832c53d5-c004-4f71-8947-4146529442f7
# ╠═19951420-f1c2-47ea-9be0-a5db4ceebac0
# ╠═1d2daf9b-e0df-4e12-bd97-956f975e717b
# ╠═76365437-d304-4a4d-be35-057d2ae39927
