module TestSens

using VoronoiFVM, ExtendableGrids, SparseArrays, ExtendableSparse
using LinearSolve, LinearAlgebra
using Krylov, Pardiso, ILUZero
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
        control = SolverControl(GMRESIteration(MKLPardisoLU(), EquationBlock()), sys)
        # control = SolverControl(nothing, sys)
    end
    control.handle_exceptions=true
    control.Δu_opt = 100
    solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose=verbose,log=true)

    return solt,grid,sys,data

end

function run_sensititivity(;dim=2,constant_irradiation_flux_bc=true)

    solt,grid,sys,data = run(
        dim=dim,
        constant_irradiation_flux_bc=constant_irradiation_flux_bc,
        times=[0,20.0],
        verbose="nea"
    )


    if dim == 2
        control = SolverControl(nothing, sys)
    else
        control = SolverControl(GMRESIteration(MKLPardisoLU(), EquationBlock()), sys)
    end
    control.handle_exceptions=true
    control.Δu_opt = 100


    sol_steadystate = VoronoiFVM.solve(
		sys;
		time = 100.0,
		inival=solt(solt.t[end]),
        control=control,
		verbose="nea"
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
        data_sens.sp = P

        inival=unknowns(sys_sens)
        inival .= sol_steadystate

        # control = SolverControl(nothing, sys_sens)
        control = SolverControl(
            DirectSolver(;factorization=SparspakFactorization()),
            sys_sens
        )

        # @info control.method_linear
        # @info control.precon_linear
        # control.handle_exceptions=true
        # control.Δu_opt = 100

        sol = VoronoiFVM.solve(sys_sens; time=100.0, inival=inival, control, verbose="nea")

        
        [maximum(sol[data.iT,:])]
        
    end

    function h(P)
        # Tv = eltype(P)
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
                # sp=ones(Tv,1)
            )

            # _,sys_sens = PTR_init_system(dim, grid, data_sens)
            _,sys_sens = PTR_init_system_sens(dim, grid, data_sens)
    

        end
        # data_sens.sp = P 

        inival=unknowns(sys_sens)

        inival .= sol_steadystate

        # control = SolverControl(nothing, sys_sens)
        # control.handle_exceptions=true
        # control.Δu_opt = 100

        sol = VoronoiFVM.solve(sys_sens; time=100.0, tstep=Inf, inival=inival, control, verbose="nea", params = [P[1]])

        return sol
        # sol = VoronoiFVM.solve(sys_sens; inival=inival, verbose="nea", params = [P[1]])


        maxu = [maximum(sol[data.iT,:])]

        dudp = sys_sens.matrix \ vec(sys_sens.dudp[1])

        return maxu, dudp
        
    end

    P = 0.1:0.05:0.3
    # P = 0.0

    # Sensitivity Method H
    ###########################################################
    # h(P[1])
    # H = zeros(0)
    # DH = zeros(0)
    # @time for p ∈ P
    #     maxu, dudp = h(p)
    #     # push!(H, maxu)
    #     # push!(DH, dudp)
    # end
    # return (H, DH)
    ###########################################################

    # Sensitivity Method G
    ###########################################################
    dresult = DiffResults.JacobianResult(ones(1))
    G = zeros(0)
    DG = zeros(0)
    @time for p ∈ P
        ForwardDiff.jacobian!(dresult, g, [p])
        push!(G, DiffResults.value(dresult)[1])
        push!(DG, DiffResults.jacobian(dresult)[1])
    end

    return (G, DG)
    ###########################################################


    # vis = GridVisualizer(; Plotter=GLMakie, legend = :lt)
    # scalarplot!(vis, P, G; color = :red, label = "g")
    # scalarplot!(vis, P, DG; color = :blue, label = "dg", clear = false, show = true)

end


function PTR_init_system_sens(dim, grid, data::ReactorData; assembly=:edgewise, unknown_storage=:dense)

	(;p,ip,Tamb,iT,iTw,iTp,ibf,inlet_boundaries,irradiated_boundaries,outlet_boundaries,catalyst_regions,X0,solve_T_equation,nom_flux,constant_irradiation_flux_bc,sp)=data
	ng=ngas(data)
	Tv = eltype(sp)

	sys=VoronoiFVM.System( 	grid;
							valuetype = Tv,
							data=data,
							flux=DMS_flux,
							reaction=DMS_reaction,
							storage=DMS_storage,
							bcondition=PTR_bcond,
							bflux=PTR_bflux,
							bstorage=PTR_bstorage,
							boutflow=DMS_boutflow,
							outflowboundaries=outlet_boundaries,
							assembly=assembly,
							unknown_storage=unknown_storage,
                            nparams=1
							)

	if solve_T_equation
		enable_species!(sys; species=collect(1:(ng+2))) # gas phase species xi, ptotal & T
		enable_boundary_species!(sys, iTw, irradiated_boundaries) # window temperature as boundary species in upper chamber

		# for 3 dimensional domain, apply measured irradiation flux density as boundary condition
		if dim == 3 && !constant_irradiation_flux_bc
			# boundary flux species, workaround to implement spatially varying irradiation
			enable_boundary_species!(sys, ibf, irradiated_boundaries)
		end
		enable_boundary_species!(sys, iTp, outlet_boundaries) # plate temperature as boundary species in lower chamber
	else
		enable_species!(sys; species=collect(1:(ng+1))) # gas phase species xi, ptotal
	end

	inival=unknowns(sys)
	inival[ip,:].=p

	for i=1:ng
		inival[i,:] .= X0[i]
	end

	if solve_T_equation
		inival[[iT,iTw,iTp],:] .= Tamb
		if dim == 3 && !constant_irradiation_flux_bc

            FluxIntp = flux_interpol(nom_flux)
			function d3tod2(a,b)
				a[1]=b[1]
				a[2]=b[2]
			end
			inival[ibf,:] .= 0.0
			sub=ExtendableGrids.subgrid(grid,irradiated_boundaries,boundary=true, transform=d3tod2 )
				
			for inode in sub[CellNodes]
				c = sub[Coordinates][:,inode]
				inodeip = sub[ExtendableGrids.NodeParents][inode]
				inival[ibf,inodeip] = FluxIntp(c[1]-0.02, c[2]-0.02)
			end
		end
	end

	if dim > 1
		catalyst_nodes = []
		for reg in catalyst_regions
			catalyst_nodes = vcat(catalyst_nodes, unique(grid[CellNodes][:,grid[CellRegions] .== reg]) )
		end
			
		cat_vol = sum(nodevolumes(sys)[unique(catalyst_nodes)])

		data.lcat = data.mcat/cat_vol
		local Ain = 0.0
		for boundary in inlet_boundaries
			Ain += bareas(boundary,sys,grid)
		end
		data.mfluxin = data.mflowin / Ain
	end
	
	return inival,sys
end


function runtests()
    solt,grid,sys,data = run()
    sol = solt(solt.t[end])
    # @test isapprox(minimum(sol[data.iT,:]), 410.8628150859246)
    @test trunc(minimum(sol[data.iT,:]), digits=4) == 410.8628
end

end
