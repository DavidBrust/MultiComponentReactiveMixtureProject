module PTR3D

using VoronoiFVM
using LessUnitful
using Pardiso
using ExtendableSparse
# using Test
# using ExampleJuggler
using DataFrames
using CSV
using Dates
using Printf
using MultiComponentReactiveMixtureProject

function RunParSweep(;
	nref=0,
	desc="nref_$(nref)"
	)

	# master_data holds information on chemical species involved,
	# chemical kinetics used, etc.
	master_data = ReactorData()

	df = init_DF(master_data)

	# variation of operating conditions
	NOMFLUX_MIN = 40.0ufac"kW/m^2"
	NOMFLUX_MAX = 80.0ufac"kW/m^2"
	N = 2
	NOMFLUX_STEP = (NOMFLUX_MAX-NOMFLUX_MIN)/N
	# NOMFLUX_STEP = 45.0ufac"kW/m^2"

	# NFLOWIN_H2_MIN = 0.0
	# NFLOWIN_H2_MAX = 10.0ufac"mol/hr"
	NFLOWIN_H2_MIN = 2.0ufac"mol/hr"
	NFLOWIN_H2_MAX = 6.0ufac"mol/hr"
	NFLOWIN_H2_STEP = (NFLOWIN_H2_MAX-NFLOWIN_H2_MIN)/N
	# NFLOWIN_H2_STEP = 5.0ufac"mol/hr"

	NFLOWIN_CO2_MIN = 2.0ufac"mol/hr"
	NFLOWIN_CO2_MAX = 6.0ufac"mol/hr"
	NFLOWIN_CO2_STEP = (NFLOWIN_CO2_MAX-NFLOWIN_CO2_MIN)/N
	# NFLOWIN_CO2_STEP = 5.0ufac"mol/hr"

	for nom_flux = NOMFLUX_MIN:NOMFLUX_STEP:NOMFLUX_MAX
		for nflowin_H2 = NFLOWIN_H2_MIN:NFLOWIN_H2_STEP:NFLOWIN_H2_MAX
			for nflowin_CO2 = NFLOWIN_CO2_MIN:NFLOWIN_CO2_STEP:NFLOWIN_CO2_MAX

				data = SetupInputData(
					master_data,
					nom_flux=nom_flux,
					nflowin_H2=nflowin_H2,
					nflowin_CO2=nflowin_CO2,
				)	
				
				writeInput_ToDF!(data, df)

				# @printf "NOM_FLUX: %e\tNFLOWIN_H2: %e\n" nom_flux/ufac"kW/m^2" nflowin_H2/ufac"mol/hr"
				@printf "NOM_FLUX: %e\tNFLOWIN_H2: %e\tNFLOWIN_CO2: %e\n" nom_flux/ufac"kW/m^2" nflowin_H2/ufac"mol/hr" nflowin_CO2/ufac"mol/hr"
				solt,grid,sys,data = run_transient(data,nref=nref)
				solss,grid,sys,data = run_stationary(solt,grid,sys,data)

				writeOutput_ToDF!(solss, sys, grid, data, df)
			end
		end
	end

	t = now()
    tm = "$(hour(t))_$(minute(t))_$(second(t))"
    desc = isempty(desc) ? desc : "_"*desc
    path = "../data/out/$(Date(t))/$(tm)_datatable$(desc)"
	# path = "../data/out/$(Date(t))/$(tm)/"
    try
        mkpath(path)
    catch e
        println("Directory " * path * " already exists.")
    end
	transform(col,val::Float64) = @sprintf("%e", val)
	df_ = df[:, Cols(:IN_MOLAR_F_CO2, :IN_MOLAR_F_H2, :IN_IRR_FLUX, 10:17)]
	CSV.write("$(path)/$(tm)_PTR_DATA_INTERPOL.csv", df_, transform=transform)
	# CSV.write("$(path)/$(tm)_PTR_DATA_INTERPOL.csv", df, transform=transform)	
end

function init_DF(data)

	(;ng, iT, ip, gn) = data

	df = DataFrame()

	for i=1:ng		
		df[:, "IN_MOLAR_F_"*String(gn[i])] = Float64[] # mol/s
	end
	df[:, "IN_P_GAS"] = Float64[] # ip = ng + 1; Pa
	df[:, "IN_T_GAS"] = Float64[] # iT = ip + 1; °C
	
	df[:, "IN_IRR_FLUX"] = Float64[] # iIrr = iT + 1; W/m2

	for i=1:ng		
		df[:, "OUT_MOLAR_F_"*String(gn[i])] = Float64[] # mol/s
	end
	df[:, "OUT_DELTA_P"] = Float64[] # ip = ng + 1; Pa
	df[:, "OUT_T_GAS"] = Float64[] # iT = ip + 1; °C
	

	return df
end

function writeInput_ToDF!(data, df)
	
	(;ng, ip, iT, nflowin, X0, p, T_gas_in, nom_flux) = data
	push!(df, zeros(Float64, ncol(df)))
	ri = nrow(df)
	for i=1:ng
		df[ri, i] = nflowin * X0[i]
	end
	# in simulation, pressure at outlet is specified via p
	iOut = iT + 1 # corresponds to iIrr (last input variable)
	df[ri, iOut+ip] = p
	# df[ri, ip] = p
	df[ri, iT] = T_gas_in - 273.15 # convert K to °C
	df[ri, iT+1] = nom_flux

	return nothing
end

function writeOutput_ToDF!(sol, sys, grid, data, df)

	(;ng, ip, iT, m) = data
	iOut = iT + 1 # corresponds to iIrr (last input variable)
	_, nflowout = MultiComponentReactiveMixtureProject.BoundaryFluxes(sol,sys,data)

	ri = nrow(df)
	for i=1:ng
		df[ri, iOut+i] = -nflowout[i]/m[i]
	end

	(Avg_inb, Avg_outb) = MultiComponentReactiveMixtureProject.avg_in_out(sol,sys,grid,data)

	# in simulation, pressure at outlet is specified via p and known beforehand
	# the pressure at the inlet boundary results from simulationinpu
	df[ri, ip] = Avg_inb[ip]
	# Dp = Avg_inb[ip] - Avg_outb[ip]
	# df[ri, iOut+ip] = Dp
	

	Tout = Avg_outb[iT] - 273.15 # convert K to °C
	df[ri, iOut+iT] = Tout

	return nothing
end

function SetupInputData(data;
						dim=3,
						nflowin_H2=3.7*ufac"mol/hr",
						nflowin_CO2=3.7*ufac"mol/hr",
						nom_flux=70.0*ufac"kW/m^2",
						T_gas_in=273.15 + 25,
						p=1.0*ufac"bar",
						mcat=3000*ufac"mg",
	)

	(;ng,gni) = data
	nflowin = nflowin_H2 + nflowin_CO2
	X0 = zeros(Float64, ng)
	X0[gni[:H2]] = nflowin_H2/nflowin
	X0[gni[:CO2]] = nflowin_CO2/nflowin


	ReactorData(
		# dt_hf_irrad = (3.0,20.0),
		dim = dim,
		kinpar = XuFroment,
		constant_irradiation_flux_bc = false,
		# is_reactive = false,
		nom_flux = nom_flux,
		T_gas_in = T_gas_in,
		X0 = X0,
		nflowin = nflowin,
		p = p,		
        mcat = mcat,

		perm = [0.0,1.0,1.0]*1.23e-10*ufac"m^2",
		poros = [0.0,0.33,0.33],				
        rhos = 5.0ufac"kg/m^3",
		Nu = 0.0,		
		include_dpdt=true		
	)

end

function run_transient(data; times=[0,15.0], nref=0, ufac=ufac"cm", W_block=3.0ufac, W_window=12.0ufac, H_cat=0.5ufac, mcat=3000ufac"mg", efix=1.0e-2*ufac,)

	(;dim) = data
	grid = PTR_grid_boundaries_regions!(dim,data; ufac=ufac, nref=nref, W_block=W_block, W_window=W_window, H_cat=H_cat, efix=efix)
	
	
	inival,sys = PTR_init_system(dim, grid, data)
	
	if dim == 2
		# times = [0,1000.0]
		control = SolverControl(nothing, sys;)
	elseif dim == 3
		# times = [0,5.0]
		#times = [0,200.0]
		control = SolverControl(
			# GMRESIteration(MKLPardisoLU(), EquationBlock()),
			GMRESIteration(MKLPardisoLU()),
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
	elseif dim == 3
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

# function runtests()
#     solt,grid,sys,data = run_transient()

#     solss,grid,sys,data = run_stationary(solt,grid,sys,data)
#     # @test isapprox(minimum(solss[data.iT,:]), 410.60941331349306) # 100 suns
#     @test isapprox(minimum(solss[data.iT,:]), 382.19131101873813) # 70 suns
# end

function run(data)
    solt,grid,sys,data = run_transient(data)
    solss,grid,sys,data = run_stationary(solt,grid,sys,data)

    # grid_aspect, inb,irrb,outb,sb,catr = PTR_grid_boundaries_regions(dim;nref=nref,H=0.5*400,W=16*100);
	grid_aspect = PTR_grid_boundaries_regions!(dim,data;nref=nref,ufac=ufac"m",H=2.0,W=16.0);
	
    WriteSolution3D(solss,grid_aspect,data,desc="nref_0_aspect_4_rescale")
end

end