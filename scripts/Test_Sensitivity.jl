module TestSens

using VoronoiFVM, ExtendableGrids, SparseArrays, ExtendableSparse
# using LinearSolve, LinearAlgebra
# using Krylov, Pardiso, ILUZero
using ForwardDiff, DiffResults
using LessUnitful, Test, ExampleJuggler
using GridVisualize, GLMakie
using MultiComponentReactiveMixtureProject


function run(;dim=3, constant_irradiation_flux_bc=false, times=[0,20.0], nref=0, verbose="nae")

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
    solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose=verbose,log=true)

    return solt,grid,sys,data

end

function run_sensititivity(;dim=2)

    solt,grid,sys,data = run(
        dim=dim,
        constant_irradiation_flux_bc=true,
        times=[0,50.0],
        verbose="a")

    @assert data.dim == 2

    # if dim == 2
    #     control = SolverControl(nothing, sys)
    # else
    #     control = SolverControl(GMRESIteration(MKLPardisoLU(), EquationBlock()), sys)
    # end

    sol_steadystate = VoronoiFVM.solve(
		sys;
		time = 100.0,
		inival=solt(solt.t[end]),
        # control,
		verbose="a"
	)

    sys_sens = nothing
    data_sens = nothing


    function g(P)
        Tv = eltype(P)
        if isnothing(sys_sens) || isnothing(data_sens)
            data_sens = ReactorData(
                dim = data.dim,
                nflowin = data.nflowin,
                nom_flux = data.nom_flux,
                dt_hf_irrad = data.dt_hf_irrad,
                dt_hf_enth = data.dt_hf_enth,
                T_gas_in = data.T_gas_in,
                Nu = data.Nu,
                X0 = data.X0,
                constant_irradiation_flux_bc = data.constant_irradiation_flux_bc,
                inlet_boundaries = data.inlet_boundaries,
                irradiated_boundaries = data.irradiated_boundaries,
                outlet_boundaries = data.outlet_boundaries,
                side_boundaries = data.side_boundaries,
                catalyst_regions = data.catalyst_regions,
                rhos = data.rhos, 
                sp=ones(Tv,1)
            )

            _,sys_sens = PTR_init_system(dim, grid, data_sens)
    

        end
        data_sens.sp = P[1]

        inival=unknowns(sys_sens)
        inival .= sol_steadystate

        sol = solve(sys_sens; time=100.0, inival=inival, verbose="na")
        [maximum(sol[data.iT,:])]
        
    end

   
    dresult = DiffResults.JacobianResult(ones(1))

    P = 0.1:0.05:0.3
    G = zeros(0)
    DG = zeros(0)
    @time for p ∈ P
        ForwardDiff.jacobian!(dresult, g, [p])
        push!(G, DiffResults.value(dresult)[1])
        push!(DG, DiffResults.jacobian(dresult)[1])
    end

    vis = GridVisualizer(; Plotter=GLMakie, legend = :lt)
    scalarplot!(vis, P, G; color = :red, label = "g")
    scalarplot!(vis, P, DG; color = :blue, label = "dg", clear = false, show = true)

end


function runtests()
    solt,grid,sys,data = run()
    sol = solt(solt.t[end])
    # @test isapprox(minimum(sol[data.iT,:]), 410.8628150859246)
    @test trunc(minimum(sol[data.iT,:]), digits=4) == 410.8628
end

end
