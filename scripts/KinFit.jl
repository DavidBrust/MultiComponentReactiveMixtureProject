module KinFit

using VoronoiFVM
using LessUnitful
using Pardiso
using ExtendableSparse
using ILUZero
using LinearSolve
using DataFrames, Tables
using CSV
using Dates
using Printf
using Optim
using Plots
using MultiComponentReactiveMixtureProject

const kinpar = XuFroment

Base.@kwdef struct PCData
    m_cat::Float64 = 500.0ufac"mg"
    W_window::Float64 = 12.0ufac"cm"
    W_block::Float64 = 3.0ufac"cm"
    H_cat::Float64 = 0.5ufac"mm"
end

const PCData23 = PCData(
    m_cat=500.0ufac"mg",
    W_window=12.0ufac"cm",
    W_block=0.0ufac"cm",
    H_cat=0.5ufac"mm"
)

const PCData24 = PCData(
    m_cat=3000.0ufac"mg",
    W_window=12.0ufac"cm",
    W_block=3.0ufac"cm",
    H_cat=5.0ufac"mm"
)

# fit paramters: 
# R1 : CH4 + H2O -> 3 H2 + CO
# k1ref: k_ref(R1)
# Eact1: Eact(R1)
# R2 : CO + H2O -> CO2 + H2
# k2ref: k_ref(R2)
# Eact2: Eact(R2)
# R3 : CH4 + 2 H2O -> 4 H2 + CO2 
# k3ref: k_ref(R3)
# Eact3: Eact(R3)
fitpar_dict = Dict(
    :k1ref => (:ki_ref, :R1),
    :Eact1 => (:Ei, :R1),
    :k2ref => (:ki_ref, :R2),
    :Eact2 => (:Ei, :R2),
    :k3ref => (:ki_ref, :R3),
    :Eact3 => (:Ei, :R3),
)

obs_dict = Dict(
    :CO => :CO_OUT,
    :CH4 => :CH4_OUT,
    :CO2 => :CO2_OUT,
    :H2 => :H2_OUT,
)

fitpar_bak = deepcopy(
    [
        kinpar.ki_ref[:R1],
        kinpar.Ei[:R1],
        kinpar.ki_ref[:R2],
        kinpar.Ei[:R2],
        kinpar.ki_ref[:R3],
        kinpar.Ei[:R3]
    ]
)

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

	df[:, "T_AVG_UC"] = Float64[] # iT_UC = iT + 1; °C
	df[:, "T_AVG_LC"] = Float64[] # iT_LC = iT + 2; °C
	

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

	# temperature in inlet plane (inflow boundary), corresponds to upper chamber
	Tin = Avg_inb[iT] - 273.15 # convert K to °C
	df[ri, iOut+iT+1] = Tin

	df[ri, iOut+iT+2] = Tout

	return nothing
end

function load_data(;
    path="C:\\Users\\brus_dv\\Promotion_FPC\\ReactorDatabase\\Evaluation_FPC_System_0624\\Kinetics_fit\\Exp_data\\PC_Versuche_23\\",
    fn="PC23data.csv"
)
    DataFrame(CSV.File(path*fn))
end


function SetupInputData(kinpar;
    dim=2,
    nflowin_H2=3.7*ufac"mol/hr",
    nflowin_CO2=3.7*ufac"mol/hr",
    nom_flux=70.0*ufac"kW/m^2",
    T_gas_in=273.15 + 25,
    p=1.0*ufac"bar",
    mcat=3000*ufac"mg",
)

    (;ng,gni) = kinpar
    nflowin = nflowin_H2 + nflowin_CO2
    X0 = zeros(Float64, ng)
    X0[gni[:H2]] = nflowin_H2/nflowin
    X0[gni[:CO2]] = nflowin_CO2/nflowin


    ReactorData(
        # dt_hf_irrad = (3.0,15.0),
        dim = dim,
        # kinpar = XuFroment,
        kinpar = kinpar,
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
        include_dpdt=false
    )

end

function run_transient(data; times=[0,15.0], nref=0, ufac=ufac"cm", W_block=3.0ufac, W_window=12.0ufac, H_cat=0.5ufac, efix=1.0e-2*ufac,)

	(;dim) = data
	grid = PTR_grid_boundaries_regions!(dim,data; ufac=ufac, nref=nref, W_block=W_block, W_window=W_window, H_cat=H_cat, efix=efix)
	
	
	inival,sys = PTR_init_system(dim, grid, data)
	
	if dim == 2
		control = SolverControl(nothing, sys;)
        vb = ""
	elseif dim == 3
		control = SolverControl(;
			method_linear = KrylovJL_GMRES(
				#gmres_restart = 10,
				restart = true,
				itmax = 250,
			),
			precon_linear = VoronoiFVM.factorizationstrategy(
				MKLPardisoLU(),
				NoBlock(),
				sys
			)
   		)
        vb = "aen"
	end
    # control.handle_exceptions=true
    control.handle_exceptions=false
	control.Δu_opt=40
    # control.Δu_opt=10
	control.Δt_max=1

	solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose=vb,log=false)
    # solt=VoronoiFVM.solve(sys;inival=inival,times,control,verbose="",log=false)	
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
		# verbose="na"
        verbose=""
	)

    return sol_steadystate,grid,sys,data
end

function simulate_OF(kinpar, PCData, nti, obs_dict, obs_sel; dim=2, writeData=false)


    n_obs_per_row = size(obs_sel,1)
    
    n_data = size(nti)[1] * n_obs_per_row # selected data rows * observations per row
    
    OF_data = Vector{Float64}(undef,n_data)
    OF_calc = Vector{Float64}(undef,n_data)

    for (i,row) in enumerate(nti)
        data = SetupInputData(
            kinpar;
            dim=dim,
            nflowin_H2 = row.H2_IN*ufac"mol/hr",
            nflowin_CO2 = row.CO2_IN*ufac"mol/hr",
            nom_flux = row.IRRFLUX_IN*ufac"kW/m^2",
            T_gas_in=273.15 + 25,
            p = row.PRESSURE_OUT*ufac"bar",
            # p=1.0*ufac"bar",
            mcat=PCData.m_cat,
        )
    
        if writeData
            df = init_DF(data)
	        writeInput_ToDF!(data, df)
        end
        
        # solt,grid,sys,data = run_transient(data, W_block=0.0ufac"cm", W_window=12.0ufac"cm", H_cat=0.05ufac"cm")
        
        solt,grid,sys,data = run_transient(data, W_block=PCData.W_block, W_window=PCData.W_window, H_cat=PCData.H_cat)
        
        solss,grid,sys,data = run_stationary(solt,grid,sys,data)

        if writeData
            writeOutput_ToDF!(solss, sys, grid, data, df)
            t = now()
            tm = "$(hour(t))_$(minute(t))_$(second(t))"
                path = "../data/out/$(Date(t))/$(tm)/"
            try
                mkpath(path)
            catch e
                println("Directory " * path * " already exists.")
            end
            transform(col,val::Float64) = @sprintf("%e", val)
            
            df_ = df[:, Cols(:IN_MOLAR_F_CO2, :IN_MOLAR_F_H2, :IN_IRR_FLUX, 10:17, :T_AVG_UC, :T_AVG_LC)]

            CSV.write("$(path)/$(tm)_PTR_DATA_INTERPOL.csv", df_, transform=transform)

            if dim == 3
                grid_aspect = PTR_grid_boundaries_regions!(dim,data;nref=0,ufac=ufac"m",H=2.0,W=16.0);
                WriteSolution3D(solss,grid_aspect,data,desc="nref_0_aspect_4_rescale")
            end
        end

        (;gni,m) = data       

        _,out_ = MultiComponentReactiveMixtureProject.BoundaryFluxes(solss,sys,data)

        j = (i-1)*n_obs_per_row

        for (k,ob) in enumerate(obs_sel)
            # OF_data[j+k] = row.obs_dict[ob]
            OF_data[j+k] = getproperty(row, obs_dict[ob])            
            OF_calc[j+k] = -out_[gni[ob]] ./ m[gni[ob]] / ufac"mol/hr"
        end

    end
    return OF_data, OF_calc
end

function create_plot(fit_p;
    initial_guess=nothing,
    fitpar_sel=[:k2ref,],
    obs_sel=[:CO, :CH4, :CO2, :H2],
    kinpar=kinpar,
    path="C:\\Users\\brus_dv\\Promotion_FPC\\ReactorDatabase\\Evaluation_FPC_System_0624\\Kinetics_fit\\plots\\",
    fn="ParityPlot",
    fitpar_dict=fitpar_dict,
    fitpar_bak=fitpar_bak,
    obs_dict=obs_dict,
    df=load_data(),
    PCData=PCData23,
    dim=2,
    writeData=false
    )

    # df = load_data()
    
    nti = Tables.namedtupleiterator(df) # whole data set
    # DEBUG: first row only
    # nti = Tables.namedtupleiterator(df[1:1,:]) 

    n_rows = size(nti)[1]

    n_obs_per_row = size(obs_sel,1)


    # Simulate with XuFroment
    # reset_kinpar!()
    OF_data_XF, OF_calc_XF = simulate_OF(kinpar, PCData, nti, obs_dict, obs_sel, dim=dim, writeData=writeData)

    # Simulate with Kinetic fit
    # update_kinpar!(fitpar, kinpar, fitpar_dict, fitpar_sel)
    if !isnothing(initial_guess) # override default values by initial guess if available
        for ig in initial_guess
            # set kinetic paramter, might be fitpar or others
            # so we can set parameters (from a previous fit) for a reaction whose parameters are not fitted in current run 
            prop, Ri = fitpar_dict[ig[1]]
            getproperty(kinpar, prop)[Ri] = ig[2] 
        end
    end

    if !isnothing(fit_p) # override default values by initial guess if available
        for (i,fp) in enumerate(fit_p)
            prop, Ri = fitpar_dict[fitpar_sel[i]]
            getproperty(kinpar, prop)[Ri] = fp 
        end
    end



    OF_data_FIT, OF_calc_FIT = simulate_OF(kinpar, PCData, nti, obs_dict, obs_sel, dim=dim, writeData=writeData)
    SSE = sum(abs2, OF_data_FIT .- OF_calc_FIT) # residual: sum of squared errors

    splots = []

    # Genertate parity plot: CO
    if :CO in obs_sel
        iCO = [(i-1)*n_obs_per_row+1 for i in 1:n_rows]
        p1 = plot(title="CO", xlabel="Exp. flow / mol/hr", ylabel="Calc. flow / mol/hr", aspect_ratio=:equal)
        
        ulim = maximum(vcat(OF_data_XF[iCO],OF_calc_XF[iCO],OF_data_FIT[iCO], OF_calc_FIT[iCO]))
        plot!(p1, [0,ulim], [0,ulim], c=:black, label=:none)
        scatter!(p1, OF_data_XF[iCO], OF_calc_XF[iCO], label="XuFroment", c=1,)
        scatter!(p1, OF_data_FIT[iCO], OF_calc_FIT[iCO], label="Kinetic fit", c=2, m=:cross)

        push!(splots, p1)
    end

    # CH4
    if :CH4 in obs_sel
        iCH4 = [(i-1)*n_obs_per_row+2 for i in 1:n_rows]
        p2 = plot(title="CH4", xlabel="Exp. flow / mol/hr", ylabel="Calc. flow / mol/hr",aspect_ratio=:equal)
        
        ulim = maximum(vcat(OF_data_XF[iCH4],OF_calc_XF[iCH4],OF_data_FIT[iCH4], OF_calc_FIT[iCH4]))
        plot!(p2, [0,ulim], [0,ulim], c=:black, label=:none)
        scatter!(p2, OF_data_XF[iCH4], OF_calc_XF[iCH4], label="XuFroment", c=1)
        scatter!(p2, OF_data_FIT[iCH4], OF_calc_FIT[iCH4], label="Kinetic fit", c=2, m=:cross)

        annotate!(p2, (.65, .4),  Plots.text("SSE: "*string(round(SSE,sigdigits=4)), :left, 8, "Helvetica Bold") )
        for (i,fp) in enumerate(keys(KinFit.fitpar_dict))
            prop, Ri = fitpar_dict[fp]
            val = getproperty(kinpar, prop)[Ri]
            if (!isnothing(fitpar_sel) && fp in fitpar_sel) || (!isnothing(initial_guess) && fp in [ig[1] for ig in initial_guess])
                annotate!(p2, (.65, .3-.05*(i-1)),  Plots.text(string(fp)*": "*string(round(val,sigdigits=4)), :left, 8, "Helvetica Bold") )
            else
                annotate!(p2, (.65, .3-.05*(i-1)),  Plots.text(string(fp)*": "*string(round(val,sigdigits=4)), :left, 8, "Helvetica") )
            end
        end

        push!(splots, p2)
    end
        
    # CO2
    if :CO2 in obs_sel
        iCO2 = [(i-1)*n_obs_per_row+3 for i in 1:n_rows]
        p3 = plot(title="CO2", xlabel="Exp. flow / mol/hr", ylabel="Calc. flow / mol/hr")
        
        ulim = maximum(OF_data_XF[iCO2])
        plot!(p3, [0,ulim], [0,ulim], c=:black, label=:none)
        scatter!(p3, OF_data_XF[iCO2], OF_calc_XF[iCO2], label="XuFroment", c=1)
        scatter!(p3, OF_data_FIT[iCO2], OF_calc_FIT[iCO2], label="Kinetic fit", c=2, m=:cross)
        push!(splots, p3)
    end
    
    # H2
    if :H2 in obs_sel
        iH2 = [(i-1)*n_obs_per_row+4 for i in 1:n_rows]
        p4 = plot(title="H2", xlabel="Exp. flow / mol/hr", ylabel="Calc. flow / mol/hr")
        
        ulim = maximum(OF_data_XF[iH2])
        plot!(p4, [0,ulim], [0,ulim], c=:black, label=:none)
        scatter!(p4, OF_data_XF[iH2], OF_calc_XF[iH2], label="XuFroment", c=1)
        scatter!(p4, OF_data_FIT[iH2], OF_calc_FIT[iH2], label="Kinetic fit", c=2, m=:cross)
        push!(splots, p4)
    end

    

    # p = plot(p1,p2,p3,p4, layout=[2,2], size=(1000,800))
    # p = plot(splots[1], layout=length(splots), aspect_ratio=:equal)

    if length(obs_sel) <= 2
        p = plot(splots..., size=(700,350), aspect_ratio=:equal)
    else
        p = plot(splots..., size=(700,700), layout=4, aspect_ratio=:equal)
    end    
    
    _t = now()
    tm = "$(hour(_t))_$(minute(_t))_$(second(_t))"
    
    # DEBUG: do not export figure
    savefig(p, path*fn*tm*".svg")
    reset_kinpar!()
    return p
end

function update_kinpar!(fit_p, kinpar, fitpar_dict, fitpar_sel)

    for (i,fp) in enumerate(fitpar_sel)
        prop, Ri = fitpar_dict[fp]
        getproperty(kinpar, prop)[Ri] = fit_p[i]
    end
    nothing
end

# residual function, sum of errors squared (SSE): experimental and calculated species outflows (OF)
function RES_OF(fit_p, kinpar, PCData, nti, fitpar_dict, fitpar_sel, obs_dict, obs_sel)

    update_kinpar!(fit_p, kinpar, fitpar_dict, fitpar_sel)
    try
        OF_data, OF_calc = simulate_OF(kinpar, PCData, nti, obs_dict, obs_sel)
        return @views sum(abs2, OF_data .- OF_calc)
    catch err
        return Inf
    end

    
end

function reset_kinpar!(
    kinpar=kinpar,
    fitpar_bak=fitpar_bak
    )
    # reset original state of XuFroment
    kinpar.ki_ref[:R1] = fitpar_bak[1]
    kinpar.Ei[:R1] = fitpar_bak[2]
    kinpar.ki_ref[:R2] = fitpar_bak[3]
    kinpar.Ei[:R2] = fitpar_bak[4]
    kinpar.ki_ref[:R3] = fitpar_bak[5]
    kinpar.Ei[:R3] = fitpar_bak[6]
end

function run_optim(;
        initial_guess=nothing,
        # fitpar_sel=[:k2ref,], # fit only kref of rwgs
        fitpar_sel=[:Eact2,], # fit only kref of rwgs
        obs_sel=[:CO], # consider Outflows of CO & CH4 in fit
        data_row_sel=[1],
        kinpar=kinpar,
        fitpar_dict=fitpar_dict,
        fitpar_bak=fitpar_bak,
        obs_dict=obs_dict,
        df=load_data(),
        PCData=PCData23
    )
    n_fitpar = size(fitpar_sel,1)

    nti = Tables.namedtupleiterator(df[data_row_sel,:])
    
    fit_p0 = Vector{Float64}(undef,n_fitpar)
    for (i,fp) in enumerate(fitpar_sel) # set fitparameter vector to default values
        prop, Ri = fitpar_dict[fp]
        fit_p0[i] = getproperty(kinpar, prop)[Ri]
    end

    # @show fit_p0

    if !isnothing(initial_guess) # override default values by initial guess if available
        for ig in initial_guess
            
            # set kinetic paramter, might be fitpar or others
            # so we can set parameters (from a previous fit) for a reaction whose parameters are not fitted in current run 
            prop, Ri = fitpar_dict[ig[1]]
            getproperty(kinpar, prop)[Ri] = ig[2] 

            idx = findfirst(isequal(ig[1]), fitpar_sel)
            if !isnothing(idx)
                fit_p0[idx] = ig[2] # overwrite initial guess for fitpar
            end
        end
    end

    # @show kinpar.ki_ref
    # @show fit_p0

    res = optimize(fit_p->RES_OF(fit_p, kinpar, PCData, nti, fitpar_dict, fitpar_sel, obs_dict, obs_sel),
                fit_p0,
                Optim.Options(store_trace = false,
                             show_trace = true,
                             show_warnings = true))

    reset_kinpar!()

    return res
end

function run(
    kinpar=kinpar,
    fitpar_dict=fitpar_dict,
    fitpar_bak=fitpar_bak,
    obs_dict=obs_dict
    )

    # data = PCData23
    data = PCData24


    df = load_data(
        # 2023 data
        # path="C:\\Users\\brus_dv\\Promotion_FPC\\ReactorDatabase\\Evaluation_FPC_System_0624\\Kinetics_fit\\Exp_data\\PC_Versuche_23\\",
        # fn="PC23data.csv"

        # 2024 data
        path="C:\\Users\\brus_dv\\Promotion_FPC\\ReactorDatabase\\Evaluation_FPC_System_0624\\Kinetics_fit\\Exp_data\\FPC_System_24\\",
        fn="FPC24data.csv"
    )
    data_row_sel=collect(1:size(df,1)) # whole data set
    # data_row_sel=[2,4,9,11]

    # switch to select what to do
    # :parfit = only fit, no plot, return fit results
    # :plot = only plot, no fit, return parity plot with standard kin. parameters
    # :parfit_plot = fit and plot
    
    # mode = :parfit_plot
    mode = :plot

    dim = 3
    writeData = true

    obs_sel=[:CO, :CH4]
    # obs_sel=[:CO]
    # obs_sel=[:CH4]

    # R1 : CH4 + H2O <-> 3 H2 + CO
    # R2 : CO + H2O <-> CO2 + H2
    # R3 : CH4 + 2 H2O <-> 4 H2 + CO2

    fitpar_sel=[
        :k1ref
    ]

    initial_guess = [
        (:k1ref, -8.1527020),
        (:Eact1, 99931.50063),
        (:k2ref, 3.845787845),
        (:Eact2, 92565.21822),
        (:k3ref, -11.041),
    ]

    # initial_guess = nothing    

    if mode == :parfit || (mode == :parfit_plot)
        println("Branch :parfit || :parfit_plot")
        @show fitpar_sel
        @show initial_guess

        res = run_optim(
            initial_guess=initial_guess,
            fitpar_sel=fitpar_sel, # fit only kref of rwgs
            obs_sel=obs_sel, # consider Outflows of CO & CH4 in fit
            data_row_sel=data_row_sel,
            df=df,
            PCData=data
        )        

        if (mode == :parfit_plot)
            println("Branch :parfit_plot")

            fit_p = Optim.minimizer(res)

            p = create_plot(fit_p,
                    initial_guess=initial_guess,
                    fitpar_sel=fitpar_sel,
                    df=df,
                    PCData=data
                )
        else
            p = nothing
        end
        
    elseif mode == :plot
        println("Branch :plot:\n")
        res = nothing
        fit_p = nothing
        p = create_plot(fit_p,
                initial_guess=initial_guess,
                fitpar_sel=fitpar_sel,
                df=df,
                PCData=data,
                dim=dim,
                writeData=writeData
            )
    
    end

    return res ,p 

end # function run



end # module KinFit
