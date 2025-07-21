module PTR3D

using MultiComponentReactiveMixtureProject, VoronoiFVM, LessUnitful, LinearSolve, ExtendableSparse, Test, ExampleJuggler
using GridVisualize, CairoMakie

function run_transient(;dim=3, times=[0,20.0], nref=0, control=control)

    grid, inlet_boundaries, irradiated_boundaries, outlet_boundaries, side_boundaries, catalyst_regions = PTR_grid_boundaries_regions(dim,nref=nref)    

    data=ReactorData(
        dim=dim,
        nflowin = 7.4*ufac"mol/hr",
		nom_flux = 70.0*ufac"kW/m^2",		
        # nom_flux = 100.0*ufac"kW/m^2",		
		dt_hf_irrad = (2.0, 10.0),
		dt_hf_enth = (2.0, 3.0),
		T_gas_in = 273.15 + 25,
		Nu = 0.0,
		X0 = [0,0.5,0,0,0.5,0.0], # H2 / CO2 = 1/1
        inlet_boundaries=inlet_boundaries,
        irradiated_boundaries=irradiated_boundaries,
        outlet_boundaries=outlet_boundaries,
        side_boundaries=side_boundaries,
        catalyst_regions=catalyst_regions,
        rhos=5.0*ufac"kg/m^3" # set solid density to low value to reduce thermal inertia of system
    )

    inival,sys = PTR_init_system(dim, grid, data)
  
    control.handle_exceptions=true
    control.Δu_opt = 40
    #control.Δt_max = 0.5
    control.Δt_max = 2.0
    control.tol_round=1.0e-9
    control.max_round=3
    solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="nae",log=true)

    return solt,grid,sys,data

end

function run_stationary(solt,grid,sys,data;control=control)

    sol_steadystate = VoronoiFVM.solve(
		sys;
		time = solt.t[end],
		inival = solt(solt.t[end]),
        control,
        tol_round=1.0e-9,
        max_round=3,
		verbose="na"
	)

    return sol_steadystate,grid,sys,data
end

function runtests(;dim=3, nref=0)

    if dim == 3
        control = SolverControl(;
			method_linear = LinearSolve.KrylovJL_GMRES(
				precs = ExtendableSparse.ILUZeroPreconBuilder()
			),
   		)
    else
        control = SolverControl(;nothing)
    end
    
    solt,grid,sys,data = run_transient(dim=dim, nref=nref, control=control)
    solss,grid,sys,data = run_stationary(solt, grid, sys, data, control=control)

    # @test isapprox(minimum(solss[data.iT,:]), 410.60941331349306) # 100 suns
    @test isapprox(minimum(solss[data.iT,:]), 382.19131101873813) # 70 suns
end

function run(;dim=3, nref=0)

    if dim == 3
        control = SolverControl(;
			method_linear = LinearSolve.KrylovJL_GMRES(
				precs = ExtendableSparse.ILUZeroPreconBuilder()
			),
   		)
    else
        control = SolverControl(;nothing)
    end

    solt,grid,sys,data = run_transient(dim=dim, nref=nref, control=control)
    # grid = run_transient(dim=dim, nref=nref, control=control)
    # GridVisualize.gridplot(grid, Plotter=CairoMakie, resolution=(1200,900))
    # return 

    solss,grid,sys,data = run_stationary(solt, grid, sys, data, control=control)
    
    # # export grid with scaled z-axis (4x) to improve readability
    # # grid_aspect, inb,irrb,outb,sb,catr = PTR_grid_boundaries_regions(3;nref=nref,H=2,W=16);
    grid_aspect, inb,irrb,outb,sb,catr = PTR_grid_boundaries_regions(3;uf=ufac"m",nref=nref,H=2,W=16);

    # # write 3D solution data to VTK
    WriteSolution3D(solss,grid_aspect,data,desc="nref_$(nref)_aspect_4_cm")

end

end
