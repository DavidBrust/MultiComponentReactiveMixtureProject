module TestSim

using FixedBed, VoronoiFVM, LessUnitful, Pardiso, ExtendableSparse

    function run(;dim=3, times=[0,20.0])

        grid, inlet_boundaries, irradiated_boundaries, outlet_boundaries, side_boundaries, catalyst_regions = FixedBed.grid_boundaries_regions(dim)

        data=ReactorData(
            dim=dim,
            nflowin=0.8*ufac"mol/hr",
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
    

end