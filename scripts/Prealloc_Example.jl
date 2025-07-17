module Prealloc_Example

using MultiComponentReactiveMixtureProject, VoronoiFVM, LessUnitful
using ExtendableGrids, GridVisualize, CairoMakie, Colors



function run_transient(dim=2;  nref=0)

    grid, inb, irrb, outb, sb, catr =  PTR_grid_boundaries_regions(dim, nref=nref)


	data = ReactorData(
		dim=dim,
		nflowin = 7.4*ufac"mol/hr",
		nom_flux = 100.0*ufac"kW/m^2",		
		dt_hf_irrad = (2.0, 10.0),
		dt_hf_enth = (2.0, 3.0),
		T_gas_in = 273.15 + 25,
		
		X0 = [0,0.5,0,0,0.5,0.0], # H2 / CO2 = 1/1
		inlet_boundaries=inb,
		irradiated_boundaries=irrb,
		outlet_boundaries=outb,
		side_boundaries=sb,
		catalyst_regions=catr,

        is_reactive=true,
        constant_properties=false,

		include_dpdt=true
	)
	
	inival,sys = PTR_init_system(dim, grid, data)

	if dim == 2
		times = [0,1000.0]
		control = SolverControl(nothing, sys;)
	elseif dim == 3
		times = [0,200.0]
		control = SolverControl(nothing,
			# GMRESIteration(MKLPardisoLU(), EquationBlock()),
			sys
		)
	end
	control.handle_exceptions=true
	control.Δu_opt=100
	control.Δt_max=100
		
	solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="aen",log=true)

    return solt,grid,sys,data

end

function run(;nref=0, dim=2)
    println("Running PTR simulation with dimension $(dim)")
    solt,grid,sys,data = run_transient(dim, nref=nref)

    (;gn, gni, ng, nvar) = data
    println("Number of unkonws $(nvar)")

    sol = solt(solt.t[end])

    MaxCO = maximum(sol[gni[:CO],:])
    println("Maximum molefraction CO: $(MaxCO)")

    # vis=GridVisualizer(legend=:lt, title="Molar Fractions", resolution=(600,300), Plotter=CairoMakie)
	# for i=1:ng
	# 	scalarplot!(vis, grid, sol[i,:], clear=false, color=pcols[i], label=String(gn[i]))
    #     # scalarplot!(vis, grid, sol[i,:], clear=false, label=String(gn[i]))
	# end
	# reveal(vis)

end

end