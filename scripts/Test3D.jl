module Test3D

using MultiComponentReactiveMixtureProject, VoronoiFVM, ExtendableSparse, LinearSolve, Krylov, Pardiso, ILUZero, LessUnitful, Test, ExampleJuggler

function run(;dim=3, constant_irradiation_flux_bc=false, times=[0.0,20.0], nref=0)

    grid, inlet_boundaries, irradiated_boundaries, outlet_boundaries, side_boundaries, catalyst_regions = PTR_grid_boundaries_regions(dim,nref=nref)

    data=ReactorData(
        dim=dim,
        nflowin = 7.4*ufac"mol/hr",
		nom_flux = 70.0*ufac"kW/m^2",		
		dt_hf_irrad = (2.0, 10.0),
		dt_hf_enth = (2.0, 3.0),
		T_gas_in = 273.15 + 25,
		Nu = 0.0,
		X0 = [0,0.5,0,0,0.5,0.0], # H2 / CO2 = 1/1
        constant_irradiation_flux_bc=constant_irradiation_flux_bc,
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
        
        # control = SolverControl( DirectSolver(factorization = MKLPardisoFactorize()), sys)
        control = SolverControl(GMRESIteration(MKLPardisoLU(), EquationBlock()), sys)
        # control = SolverControl(GMRESIteration(ILUZeroFactorization()), sys)
        # control = SolverControl(CGIteration(MKLPardisoLU(), EquationBlock()), sys)
        # control = SolverControl(nothing, sys)
    end
    # control = SolverControl(nothing, sys)
    control.handle_exceptions=true
    control.Δu_opt = 100

    solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="nae",log=true)

    # Print_summary_ext(solt,grid,sys,data)

    return solt,grid,sys,data

end

function run_steadystate(solt,sys,grid,data)

    (;dim) = data

    if dim == 2
        control = SolverControl(nothing, sys)
    else
        control = SolverControl(GMRESIteration(MKLPardisoLU(), EquationBlock()), sys)
    end

    sol_steadystate = VoronoiFVM.solve(
		sys;
		time = 100.0,
		inival=solt(solt.t[end]),
        control,
		verbose="na"
	)

    return sol_steadystate,grid,sys,data
end

function runtests()
    solt,grid,sys,data = run()
    sol = solt(solt.t[end])
    # @test isapprox(minimum(sol[data.iT,:]), 410.8628150859246)
    @test trunc(minimum(sol[data.iT,:]), digits=4) == 422.6527
end

end
