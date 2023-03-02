module Example2_EnergyBalance

using Plots, DataFrames, CSV


include("../notebooks/PorousCatalystHot3D.jl")

function RunSim(Qflow,C;data=ModelData(Qflow=Qflow))

	grid=prism_sq(data)

	ngas=data.ng
	iT=data.iT
	
	sys=VoronoiFVM.System( 	grid;
							data=data,
							flux=flux,
							reaction=reaction,
							bcondition=bcond
							)
	enable_species!(sys; species=collect(1:(ngas+2))) # gas phase species + p + T
	
	inival=unknowns(sys)
	inival[:,:].=1.0*data.p
	for i=1:ngas
		inival[i,:] .*= data.X0[i]
	end
	inival[iT,:] .= data.Tamb

	sol=solve(sys;inival,)
	
	 function pre(sol,par)
	 	ng=data.ng
	 	iT=data.iT
	 	# iteratively adapt top outflow boundary condition
	 	function Inttop(f,u,bnode,data)
			
	 		X=zeros(eltype(u), ng)
	 		mole_frac!(bnode,data,X,u)
	 		# top boundary(cat/frit)
	 		if bnode.region==Γ_top_frit || bnode.region==Γ_top_cat  
	 			for i=1:ng
	 				f[i] = data.Fluids[i].MW*X[i]
	 			end
	 			f[iT] = u[iT]
	 		end
	 	end
		
	 	MWavg=sum(integrate(sys,Inttop,sol; boundary=true)[1:ng,[Γ_top_frit,Γ_top_cat]])/(data.Ac/4)
	 	ntop=data.mdotin/MWavg
		 
	 	Tavg=sum(integrate(sys,Inttop,sol; boundary=true)[data.iT,[Γ_top_frit,Γ_top_cat]])/(data.Ac/4)
		
	 	utop_calc=ntop*ph"R"*Tavg/(1.0*ufac"bar")/data.Ac
	 	utops=[data.utop, utop_calc]*ufac"m/s"
	 	data.utop = minimum(utops) + par*(maximum(utops)-minimum(utops))
		
		
	 	# specific catalyst loading
	 	mcats=[10.0, 1300.0]*ufac"kg/m^3"
	 	data.mcats= minimum(mcats) + par*(maximum(mcats)-minimum(mcats))

	 	# irradiation flux density
	 	G_lamp=[1.0, C]*ufac"kW/m^2"
	 	data.G_lamp= minimum(G_lamp) + par*(maximum(G_lamp)-minimum(G_lamp))
	 end
	
	 control=SolverControl( ;
	 				  		handle_exceptions=true,
							Δp_min=1.0e-4,					  
	 				  		Δp=0.1,
	 				  		Δp_grow=1.2,
	 				  		Δu_opt=1.0e7, # large value, due to unit Pa of pressure?
	 				  		)
	
	# #sol=solve(sys;inival,)
	
	 sol=solve(sys;inival, embed=[0.0,1.0],pre,control)

	
	sol,grid,sys,data
end;

function ParSweep()
    Qflow =3400.0*ufac"ml/minute" # volumetric feed flow rate (sccm)
    C = 25.0
    #fact = [0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5]
    fact = [2.0,3.0]
    l=length(fact)
    Qflows = Qflow * fact
    Cs = C * fact

    STCs = []
    datas = []
    sols = []
    grid_ = []
    sys_ = []
    for Qflow in Qflows
        for C in Cs
            sol,grid,sys,data = RunSim(Qflow, C)
            solf = sol(sol.t[end])
            push!(STCs, STCefficiency(solf,sys,data,grid)[1])
            push!(datas,data)
            push!(sols,solf)
            if Qflow==Qflows[end] && C ==Cs[end] 
                push!(grid_,grid)
                push!(sys_,sys)
            end
        end
    end
    p=Plots.plot(
            xguide="Feed Flow / sccm",
            yguide="Solar Concentration / kW m-2",
            title="Solar-To-Chemical Efficiency",
            size=(400,300)
            )
    Plots.heatmap!(p,Qflows./ufac"ml/minute",Cs,reshape(STCs,l,:))
    Plots.savefig(p, "./img/out/STC.svg")

    Df=DataFrame((C=repeat(Cs, outer=l),Qflow=repeat(Qflows, inner=l),STCs))
    CSV.write("./data/out/STC.csv",Df)
    Df, datas, sols, grid_, sys_
end

function HeatFluxes(C=[1,10])
    data=ModelData()
    iT=data.iT

    lC=length(C)
    #q_full=zeros(Float64, lC)
    q_in=zeros(Float64, lC)
    q_top_abs=zeros(Float64, lC)
    q_top_refl=zeros(Float64, lC)
    q_top_rerad=zeros(Float64, lC)
    q_top_convec=zeros(Float64, lC)
    q_sides_conv=zeros(Float64, lC)
    q_bot_rerad=zeros(Float64, lC)

    function flux_rerad_top(f,u,bnode,data)
        if bnode.region==6 # top boundary
            f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
        end
    end

    function flux_lamp_abs(f,u,bnode,data)
        if bnode.region==6 # top boundary
            f[iT] = -data.Abs_lamp*data.G_lamp
        end
    end

    function flux_convec_top(f,u,bnode,data)
        if bnode.region==6 # top boundary
            ρf=density_idealgas(data.Fluid, u[iT], data.p)
            cf=heatcap_gas(data.Fluid, u[iT])
            f[iT] = data.u0*ρf*cf*(u[iT]-data.Tamb)
        end
    end

    function sidewalls(f,u,bnode,data)
        # boundary conditions at side walls
            boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
            boundary_robin!(f,u,bnode;species=iT,region=3, factor=data.α_w, value=data.Tamb*data.α_w)
    end

    function bottom(f,u,bnode,data)
		if bnode.region==5 # bottom boundary
            f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
		end
	end

    for (i,C) in enumerate(C)
        data.G_lamp=C*1.0*ufac"kW/m^2"
        sys,sol,data,nref = main3D(data,nref=0)
        # factor 4: due to symmetry, sim domain covers only 1/4 of the frit
        q_top_abs[i] = -4*integrate(sys,flux_lamp_abs,sol; boundary=true)[6]
        q_in[i] = q_top_abs[i]/data.Abs_lamp
        q_top_refl[i] = (1.0-data.Abs_lamp)*q_in[i]
        q_top_rerad[i] = 4*integrate(sys,flux_rerad_top,sol; boundary=true)[6]
        q_top_convec[i] = 4*integrate(sys,flux_convec_top,sol; boundary=true)[6]
        q_sides_conv[i] = 4*sum(integrate(sys,sidewalls,sol; boundary=true)[1,[2,3]])
        q_bot_rerad[i] = 4*integrate(sys,bottom,sol; boundary=true)[5]

    end
    q_in,q_top_abs,q_top_refl,q_top_rerad,q_top_convec,q_sides_conv,q_bot_rerad
    
end


function PlotLosses(C=[1,10,25,50,75,100])
    q_in,q_top_abs,q_top_refl,q_top_rerad,q_top_convec,q_sides_conv,q_bot_rerad = HeatFluxes(C)
    p=plot(xguide="Irradiation / kW m⁻²", yguide="Loss contributions / -", legend=:outertopright, size=(400,250))
	areaplot!(C, [q_top_refl q_top_rerad q_top_convec q_sides_conv q_bot_rerad] ./ (q_in), fillalpha = [0.2 0.2], seriescolor = [1 2 3 4 5], label=["Refl Top" "Rerad Top" "Conv Top" "Conv Sides" "Rerad Bot"])
    savefig("img/loss_contributions.svg")

end

end
