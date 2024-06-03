### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ f7a71f60-03ba-11ef-3a89-dba4bc09d152
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,"../.."))
	using Revise
	
	using LessUnitful
	using PlutoUI, Colors, ColorSchemes, Plots
	using CSV, DataFrames, Tables, DelimitedFiles
	using Dates, Printf
	using Interpolations
	using Random
	using MultiComponentReactiveMixtureProject
end;

# ╔═╡ 4320a8da-21bd-4bfd-9927-79c3a80306e3
begin
	#nom_flux = [40.0, 60.0, 80.0, 100.0] # kW/m^2
	nom_flux = 40.0:5.0:100.0 # kW/m^2
	#n_flow_in = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] # mol/hr
	in_flow_CO2 = 1.0:0.5:7.0 # mol/hr
	in_flow_H2 = 1.0:0.5:7.0 # mol/hr
end;

# ╔═╡ d5f83ff7-1d14-44f8-a8f1-67b7f2bda688
#n_dot_CO(xs,ys) = [log(x*y) for x in xs, y in ys]
n_dot_CO(x1s,x2s,x3s) = [log(x1*x2*x3) for x1 in x1s, x2 in x2s, x3 in x3s]
#n_dot_CO(xs,ys) = [rand()*x/100 + rand()*y for x in xs, y in ys] 

# ╔═╡ bac5b175-6882-4095-b946-3fb245feb7fc
function plot_data_int(itp; x1=nom_flux, x2=n_flow_in)
	
	# orig dat
	n_dot_CO_ = n_dot_CO(x1,x2)
	clims = (minimum(n_dot_CO_), maximum(n_dot_CO_))

	p1 = contour(x1, x2, n_dot_CO_, fill=true, levels=10, clims=clims,)
	
	# test points
	x1t = [42.0, 55.0, 63.0, 72.0, 89.0, 97.0]
	#x1t = x1
	x2t = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5]	
	#x2t = x2
	zt = itp(x1t, x2t)	
	
	#x1t_ = repeat(x1t, outer = length(x2t))
	x1t_ = repeat(x1t, inner = length(x2t))
	#x2t_ = repeat(x2t, inner = length(x1t))	
	x2t_ = repeat(x2t, outer = length(x1t))	

	# colored scatter	
	zt_ = vec(reshape(zt, 1,:))
	#scatter!(xt_, yt_, color=zt_, markersize=20, strokecolor=:black, strokewidth=0)
	
	#p2 = contour(xt, yt, zt, clims=clims, fill=true, levels=10, xlims=xlims(p1), ylims=ylims(p1))
	#contour!(p1, xt, yt, zt, colorrange=cr, levels=10, labels=true)
	scatter!(p1, x1t_, x2t_, marker_z=zt_, markersize=10, legend=:none, markerstrokewidth=0)

end

# ╔═╡ d257408f-b9a7-40e0-90e5-f881871c3975
function make_itp(;x1=nom_flux, x2=in_flow_CO2, x3=in_flow_H2, nx1=2, nx2=2, nx3=2)

	
	#interp_linear = linear_interpolation((x1,x2), n_dot_CO(x1,x2))
	interp_linear = scale(interpolate(n_dot_CO(x1,x2,x3), BSpline(Linear())), (x1,x2,x3))
	
	
	#itp = scale(interpolate(A, BSpline(Quadratic())), (xs, ys))	
	#itp = interpolate(A, BSpline(Linear()))
end

# ╔═╡ 20b37b00-fa9a-4b1b-9598-effe721b232f
itp =  make_itp();

# ╔═╡ d9c9e120-9e99-46bf-a863-0006404afd1b
collect(knots(itp));

# ╔═╡ 509731f2-0601-42d1-a60b-dd4ae568e497
let
	plot_data_int(make_itp())
end

# ╔═╡ Cell order:
# ╠═f7a71f60-03ba-11ef-3a89-dba4bc09d152
# ╠═4320a8da-21bd-4bfd-9927-79c3a80306e3
# ╠═d5f83ff7-1d14-44f8-a8f1-67b7f2bda688
# ╠═20b37b00-fa9a-4b1b-9598-effe721b232f
# ╠═d9c9e120-9e99-46bf-a863-0006404afd1b
# ╠═509731f2-0601-42d1-a60b-dd4ae568e497
# ╠═bac5b175-6882-4095-b946-3fb245feb7fc
# ╠═d257408f-b9a7-40e0-90e5-f881871c3975
