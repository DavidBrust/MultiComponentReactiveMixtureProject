module PTR3D

using MultiComponentReactiveMixtureProject, VoronoiFVM, LessUnitful, Pardiso, ExtendableSparse, Test, ExampleJuggler

function run_transient(;dim=3, times=[0,20.0], nref=0)

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

    if dim == 2
        control = SolverControl(nothing, sys)
    else
        control = SolverControl(GMRESIteration(MKLPardisoLU(), EquationBlock()), sys)
    end
    control.handle_exceptions=true
    control.Δu_opt = 100
    #control.Δt_max = 0.5
    control.Δt_max = 2.0
    solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="nae",log=true)

    return solt,grid,sys,data

end

function run_stationary(solt,grid,sys,data)

    (;dim) = data

    if dim == 2
        control = SolverControl(nothing, sys)
    else
        control = SolverControl(GMRESIteration(MKLPardisoLU(), EquationBlock()), sys)
    end

    sol_steadystate = VoronoiFVM.solve(
		sys;
		time = solt.t[end],
		inival = solt(solt.t[end]),
        control,
		verbose="na"
	)

    return sol_steadystate,grid,sys,data
end

function runtests()
    solt,grid,sys,data = run_transient()

    solss,grid,sys,data = run_stationary(solt,grid,sys,data)
    # @test isapprox(minimum(solss[data.iT,:]), 410.60941331349306) # 100 suns
    @test isapprox(minimum(solss[data.iT,:]), 382.19131101873813) # 70 suns
end

function run()
    solt,grid,sys,data = run_transient(nref=1)
    solss,grid,sys,data = run_stationary(solt,grid,sys,data)

    grid_aspect, inb,irrb,outb,sb,catr = PTR_grid_boundaries_regions(3;nref=1,H=0.5*400,W=16*100);
    WriteSolution3D(solss,grid_aspect,data,desc="nref_1_aspect_4_cm")

end

end
