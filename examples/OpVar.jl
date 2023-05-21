module OpVar

using LessUnitful, GLMakie, CairoMakie

include("../notebooks/PorousCatalystHot3DTopFlowIrrExchange_NonAlloc.jl")


# function runSim(;data=ModelData())

#     sol_,grid,sys,data_=main(;data)
#     if sol_ isa VoronoiFVM.TransientSolution
# 		sol = copy(sol_(sol_.t[end]))
# 	else
# 		sol = copy(sol_)
# 	end

#     sol,grid,sys,data_
# end

function runSim(;data=ModelData())
    function control(sys;kwargs...)
        SolverControl(
#			    direct_umfpack(),
            gmres_umfpack(),
    #			gmres_eqnblock_umfpack(),
    #			gmres_eqnblock_iluzero(),
            sys;
            verbose="na",
            log=true,
            reltol=1.0e-8,
            reltol_linear=1.0e-5,
            kwargs...
        )
    end

    sol_,grid,sys,data_embed=main(;data=data,nref=0,control,
    assembly=:edgewise);
	if sol_ isa VoronoiFVM.TransientSolution
		sol = copy(sol_(sol_.t[end]))
	else
		sol = copy(sol_)
	end

    sol,grid,sys,data_embed
end

function runVar()
    #1
    (sol_,grid_,sys_,data_) = MapTmax(;
    nsteps=6,
    Qflow_min=340.0*ufac"mL/minute",
    Qflow_max=1480.0*ufac"mL/minute",
    Flux_target_min=0.45,
    Flux_target_max=1.0,
    abs_cat_IR=0.7,
    abs_cat_vis=0.7,
    )
    plotHM(sol_,grid_,sys_,data_)
    #2
    (sol_,grid_,sys_,data_) = MapTmax(;
    nsteps=6,
    Qflow_min=3400.0*ufac"mL/minute",
    Qflow_max=14800.0*ufac"mL/minute",
    Flux_target_min=0.45,
    Flux_target_max=1.0,
    abs_cat_IR=0.7,
    abs_cat_vis=0.7,
    )
    plotHM(sol_,grid_,sys_,data_)
    #3
    (sol_,grid_,sys_,data_) = MapTmax(;
    nsteps=6,
    Qflow_min=340.0*ufac"mL/minute",
    Qflow_max=1480.0*ufac"mL/minute",
    Flux_target_min=0.45,
    Flux_target_max=1.0,
    abs_cat_IR=0.5,
    abs_cat_vis=0.5,
    )
    plotHM(sol_,grid_,sys_,data_)
    #4
    (sol_,grid_,sys_,data_) = MapTmax(;
    nsteps=6,
    Qflow_min=3400.0*ufac"mL/minute",
    Qflow_max=14800.0*ufac"mL/minute",
    Flux_target_min=0.45,
    Flux_target_max=1.0,
    abs_cat_IR=0.5,
    abs_cat_vis=0.5,
    )
    plotHM(sol_,grid_,sys_,data_)
    #5
    (sol_,grid_,sys_,data_) = MapTmax(;
    nsteps=6,
    Qflow_min=340.0*ufac"mL/minute",
    Qflow_max=1480.0*ufac"mL/minute",
    Flux_target_min=0.45,
    Flux_target_max=1.0,
    abs_cat_IR=0.3,
    abs_cat_vis=0.3,
    )
    plotHM(sol_,grid_,sys_,data_)
    #6
    (sol_,grid_,sys_,data_) = MapTmax(;
    nsteps=6,
    Qflow_min=3400.0*ufac"mL/minute",
    Qflow_max=14800.0*ufac"mL/minute",
    Flux_target_min=0.45,
    Flux_target_max=1.0,
    abs_cat_IR=0.3,
    abs_cat_vis=0.3,
    )
    plotHM(sol_,grid_,sys_,data_)

end

function MapTmax(;
    nsteps=4,
    Qflow_min=340.0*ufac"mL/minute",
    Qflow_max=1480.0*ufac"mL/minute",
    Flux_target_min=0.45,
    Flux_target_max=1.0,
    # Ieff_min=45.0*ufac"kW/m^2",
    # Ieff_max=100.0*ufac"kW/m^2",
    abs_cat_IR=0.7,
    abs_cat_vis=0.7,
    )
    data=ModelData()
    NG=ngas(data)
    (;kinpar)=data
    (;gni)=kinpar
    Qflows=range(Qflow_min,Qflow_max,length=nsteps)
    #Ieffs=range(Ieff_min,Ieff_max,length=nsteps)
    Ieffs=range(Flux_target_min,Flux_target_max,length=nsteps)

    X0 = begin
		x=zeros(Float64, NG)
		#x[gni[:H2]] = 1.0
		#x[gni[:CO2]] = 1.0
        x[gni[:N2]] = 1.0
		x/sum(x)
	end

    sol_ = []
    data_ = []
    sys_ = []
    grid_ = []
    cnt=1
    for Qflow in Qflows
        for Ieff in Ieffs
            @show (Qflow,Ieff)
            sol,grid,sys,data=runSim(;
                data=ModelData(
                    isreactive=0,
                    Qflow=Qflow,
                    # Glamp_target=Ieff,
                    Flux_target=Ieff,
                    X0=X0,
                    uc_cat=SurfaceOpticalProps(
                        alpha_IR=abs_cat_IR,
                        tau_IR=0.0,
                        alpha_vis=abs_cat_vis,
                        tau_vis=0.0
                    )
                    
                    )
            )
            push!(sol_,sol)
            push!(data_,data)

            push!(sys_,sys)
            push!(grid_,grid)
            
        end
    end

    (sol_,grid_,sys_,data_)

end

function TmaxSide(sol,sys,grid,data)

    (;iT)=data
	if grid_fun == prism_sq_full 
		sub=subgrid(grid,[Γ_side_front,Γ_side_back,Γ_side_right,Γ_side_left],boundary=true)
	else
		sub=subgrid(grid,[Γ_side_back,Γ_side_right],boundary=true)
	end

    Tsides = view(sol[iT, :], sub) .- 273.15
    maximum(Tsides)
		
	
end

function plotHM(sol_,grid_,sys_,data_)
    grid=grid_[1]
    sys=sys_[1]



    len=length(sol_)
    Qflows = zeros(len)
    Ieffs = zeros(len)
    Tmaxs = zeros(len)
    for i=1:len
        Qflows[i]= data_[i].Qflow
        # Ieffs[i]= data_[i].Glamp
        (;wi,Flux_target)=data_[i]
        pwr=sum(PowerIn(sol_[i],sys,grid,data_[i])[[Γ_top_cat,Γ_top_frit]])*Flux_target
        Ieffs[i]= pwr/wi^2/ufac"kW/m^2"
        
        #Tmaxs[i] = Tcatavg(sol_[i],sys,grid,data_[i])[1]
        Tmaxs[i] = TmaxSide(sol_[i],sys,grid,data_[i])
    end
    Qflows=unique(Qflows)./ufac"mL/minute"
    #Ieffs=unique(Ieffs)./ufac"kW/m^2"
    Ieffs=unique(Ieffs)

    Tmaxs=reshape(Tmaxs,size(Qflows,1),:)

    fig = GLMakie.Figure(fontsize=20)
    ax = Axis(
        fig[1,1],
        title= L"\text{Max. Temperature}\,/\,\degree\text{C}",
        xlabel = L"\text{Irradiation}\,/\,\text{kW}\,\text{m}^{-2}",
        ylabel = L"\text{Volumetric Flow}\,/\, \text{mL}\, \text{min}^{-1}",
        xticks = (Ieffs, @. string(Integer(round(Ieffs)))),
        yticks = (Qflows, @. string(Integer(round(Qflows))))
        )
    hm=GLMakie.heatmap!(ax,Ieffs,Qflows,Tmaxs)
    # GLMakie.Colorbar(fig[1, 2], hm, label = L"\text{Max. Temperature}\,/\,\degree\text{C}")
    GLMakie.Colorbar(fig[1, 2], hm,)
    
    for i in eachindex(Ieffs), j in eachindex(Qflows)
        txtcolor = Tmaxs[i, j] < 500.0 ? :white : :black
        GLMakie.text!(ax, "$(Integer(round(Tmaxs[i,j])))", position = (Ieffs[i], Qflows[j]),
            color = txtcolor, align = (:center, :center))
    end

    # add additional info on configuration via textbox an the side
    #tb = GLMakie.Box(fig[1,3], color=:white)
    ax3 = Axis(
        fig[1,3],
        )
    GLMakie.colsize!(fig.layout, 3, Relative(1/4))
    hidedecorations!(ax3)
    
    GLMakie.text!(
        ax3,0.1,1,
        text = "Optical\nProperties",
        space = :relative,
        fontsize = 24,
        offset = (0,-60)
    )
    
    
    abs_cat_IR=data_[end].uc_cat.alpha_IR
    GLMakie.text!(
        ax3,0.1,0.85,
        text = L"\alpha^{\text{IR}}_{\text{cat}} = %$abs_cat_IR",
        space = :relative,
        fontsize = 20,
        offset = (0,-30)
    )
    abs_cat_vis=data_[end].uc_cat.alpha_vis
    GLMakie.text!(
        ax3,0.1,0.8,
        text = L"\alpha^{\text{vis}}_{\text{cat}} = %$abs_cat_vis",
        space = :relative,
        fontsize = 20,
        offset = (0,-30)
    )
    abs_frit_IR=data_[end].uc_frit.alpha_IR    
    GLMakie.text!(
        ax3,0.1,0.75,
        text = L"\alpha^{\text{IR}}_{\text{frit}} = %$abs_frit_IR",
        space = :relative,
        fontsize = 20,
        offset = (0,-30)
    )
    abs_frit_vis=data_[end].uc_frit.alpha_vis
    GLMakie.text!(
        ax3,0.1,0.7,
        text = L"\alpha^{\text{vis}}_{\text{frit}} = %$abs_frit_vis",
        space = :relative,
        fontsize = 20,
        offset = (0,-30)
    )
    tau_window_vis=data_[end].uc_window.tau_vis
    GLMakie.text!(
        ax3,0.1,0.65,
        text = L"\tau^{\text{vis}}_{\text{window}} = %$tau_window_vis",
        space = :relative,
        fontsize = 20,
        offset = (0,-30)
    )
    tau_window_IR=data_[end].uc_window.tau_IR
    GLMakie.text!(
        ax3,0.1,0.6,
        text = L"\tau^{\text{IR}}_{\text{window}} = %$tau_window_IR",
        space = :relative,
        fontsize = 20,
        offset = (0,-30)
    )
    abs_window_IR=data_[end].uc_window.alpha_IR
    GLMakie.text!(
        ax3,0.1,0.55,
        text = L"\alpha^{\text{IR}}_{\text{window}} = %$abs_window_IR",
        space = :relative,
        fontsize = 20,
        offset = (0,-30)
    )

    GLMakie.text!(
        ax3,0.1,0.45,
        text = "Gas Feed",
        space = :relative,
        fontsize = 24,
        offset = (0,-30)
    )

    GLMakie.text!(
        ax3,0.1,0.35,
        text = L"N_2 = 100 %",
        space = :relative,
        fontsize = 20,
        offset = (0,-30)
    )

    # fn="img/out/230517/OpVarHM_AbsCatVis$(abs_cat_vis)_AbsCatIR$(abs_cat_IR)_Qflow$(Integer(round(maximum(Qflows)))).png"
    # GLMakie.save(fn,fig)
    fn="img/out/230517/OpVarHM_AbsCatVis$(abs_cat_vis)_AbsCatIR$(abs_cat_IR)_Qflow$(Integer(round(maximum(Qflows)))).svg"
    CairoMakie.save(fn,fig)
    
    fig

end

end