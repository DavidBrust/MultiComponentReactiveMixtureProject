module TestSim

using FixedBed, VoronoiFVM, LessUnitful, Pardiso, ExtendableSparse

function run(;dim=3, times=[0,20.0])

    grid, inlet_boundaries, irradiated_boundaries, outlet_boundaries, side_boundaries, catalyst_regions = FixedBed.grid_boundaries_regions(dim)

    data=ReactorData(
        dim=dim,
        nflowin = 7.4*ufac"mol/hr",
		nom_flux = 100.0*ufac"kW/m^2",		
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

    inival,sys = FixedBed.init_system(dim, grid, data)

    if dim == 2
        control = SolverControl(nothing, sys;)
    else
        control = SolverControl(GMRESIteration(MKLPardisoLU(), EquationBlock()), sys)
    end
    control.handle_exceptions=true
    control.Δu_opt=100
    println()
    solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="nae",log=true)

    return solt,grid,sys,data

end

function Summary(solt,grid,sys,data)

    t = solt.t[end]
    sol = solt(t)
    FixedBed.DMS_print_summary_ext(sol,sys,data)
    #return HeatFluxes_EB_I(t,solt,grid,sys,data)

end

    

end