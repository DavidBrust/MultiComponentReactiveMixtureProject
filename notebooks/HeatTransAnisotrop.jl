### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 5c3adaa0-9285-11ed-3ef8-1b57dd870d6f
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	
	using VoronoiFVM
	using ExtendableGrids
	using GridVisualize
	using LessUnitful
	using PlutoVista	
	using PlutoUI

	using FixedBed
	
	GridVisualize.default_plotter!(PlutoVista)
end;

# ╔═╡ 7d8eb6f5-3ba6-46ef-8058-1f24a0938ed1
PlutoUI.TableOfContents(title="Heat Transfer in Fixed Beds")

# ╔═╡ Cell order:
# ╠═7d8eb6f5-3ba6-46ef-8058-1f24a0938ed1
# ╠═5c3adaa0-9285-11ed-3ef8-1b57dd870d6f
