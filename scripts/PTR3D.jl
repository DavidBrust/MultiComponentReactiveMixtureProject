module PTR3D

using MultiComponentReactiveMixtureProject, VoronoiFVM, LessUnitful, Pardiso, ExtendableSparse, Test, ExampleJuggler

# function run_transient(;dim=3, times=[0,20.0], nref=0)
function run_transient(;dim=3, times=[0,20.0], nref=0, ufac=ufac"cm", W_block=3.0ufac, W_window=12.0ufac, H_cat=0.5ufac, mcat=3000ufac"mg")

    data = ReactorData(
		dim=dim,
		kinpar=XuFroment,
		constant_irradiation_flux_bc = false,
		nom_flux = 70.0*ufac"kW/m^2",
		T_gas_in = 273.15 + 25,
		X0 = [0.0,0.5,0.0,0.0,0.5,0.0],
		perm = [0.0,1.0,1.0]*1.23e-10*ufac"m^2",
		poros = [0.0,0.33,0.33],
				
        rhos = 5.0ufac"kg/m^3",

		p = 1.0*ufac"bar",
		nflowin = 7.4*ufac"mol/hr",
        mcat = mcat,

		Nu = 0.0,		
		include_dpdt=true		
	)

	grid = PTR_grid_boundaries_regions!(dim,data; ufac=ufac, nref=nref, W_block=W_block, W_window=W_window, H_cat=H_cat)
	
	
	inival,sys = PTR_init_system(dim, grid, data)
	
	if dim == 2
		times = [0,1000.0]
		control = SolverControl(nothing, sys;)
	elseif dim == 3
		# times = [0,5.0]
		#times = [0,200.0]
		control = SolverControl(
			GMRESIteration(MKLPardisoLU(), EquationBlock()),
			sys
		)
	end
    control.handle_exceptions=true
	control.Δu_opt=100
	control.Δt_max=1

	solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="aen",log=true)	
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

function run(;dim=3,nref=0)
    solt,grid,sys,data = run_transient(dim=dim,nref=nref)
    solss,grid,sys,data = run_stationary(solt,grid,sys,data)

    # grid_aspect, inb,irrb,outb,sb,catr = PTR_grid_boundaries_regions(dim;nref=nref,H=0.5*400,W=16*100);
	grid_aspect = PTR_grid_boundaries_regions!(dim,data;nref=nref,ufac=ufac"m",H=2.0,W=16.0);
	
    WriteSolution3D(solss,grid_aspect,data,desc="nref_0_aspect_4_rescale")
end

end
