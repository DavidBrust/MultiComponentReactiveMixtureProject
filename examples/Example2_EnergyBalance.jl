module Example2_EnergyBalance

using DataFrames, CSV


include("../notebooks/PorousCatalystHot3DTopFlowIrrExchange.jl")

function RunSim()
    data=ModelData(
        isreactive=0
        )
    sol_embed,grid,sys,data_embed = main(data=data);
    sol = sol_embed(1.0)
    sol,grid,sys,data_embed
end;

function flux_in_aperture(f,u,bnode,data)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
        iT=data.iT
        f[iT] = data.Glamp
    end
end

function flux_radiation_top(f,u,bnode,data)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit

        (;ng,iT,Glamp,Tglass,uc_window,uc_cat,uc_frit,vf_uc_window_cat,vf_uc_window_frit)=data
        σ=ph"σ"

        # irradiation exchange between quartz window (1), cat surface (2) & frit surface (3)
        # window properties (1)
        tau1_vis=uc_window.tau_vis
        rho1_vis=uc_window.rho_vis
        tau1_IR=uc_window.tau_IR
        rho1_IR=uc_window.rho_IR
        eps1=uc_window.eps
        # catalyst layer properties (2)
        rho2_vis=uc_cat.rho_vis
        rho2_IR=uc_cat.rho_IR
        eps2=uc_cat.eps
        # uncoated frit properties (3)
        rho3_vis=uc_frit.rho_vis
        rho3_IR=uc_frit.rho_IR
        eps3=uc_frit.eps
        #view factors
        ϕ12=vf_uc_window_cat
        ϕ13=vf_uc_window_frit


        # surface brigthness of quartz window inwards / towards catalyst facing surface (1) in vis & IR	
        G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))
		# here the simplification is applied: only local value of T (u[iT]) is available, it is used for both surfaces
		G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
                    

        G2_vis = rho2_vis*G1_bot_vis
        G2_IR = eps2*σ*u[iT]^4 + rho2_IR*G1_bot_IR

        G3_vis = rho3_vis*G1_bot_vis
        G3_IR = eps3*σ*u[iT]^4 + rho3_IR*G1_bot_IR

        # surface brigthness of quartz window outward facing surface (1) in vis & IR	
        G1_top_vis = rho1_vis*Glamp + tau1_vis*(ϕ12*G2_vis + ϕ13*G3_vis)
        G1_top_IR = tau1_IR*(ϕ12*G2_IR + ϕ13*G3_IR) + eps1*σ*Tglass^4

    
        iT=data.iT
        f[iT] = G1_top_vis + G1_top_IR
    end
end

function flux_conduction_top(f,u,bnode,data)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
        (;iT,Tglass,X0,uc_h)=data

        # height of upper chamber / distance over which heat conduction occurs
        Tsurf = u[iT]

        # mean temperature
        Tm=0.5*(u[iT] + Tglass)
        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        _,λf=dynvisc_thermcond_mix(data, Tm, X0)

        f[iT] = λf*(Tsurf-Tglass)/uc_h

    end
end


function flux_enthalpy_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
		(;ng,iT,ip,Fluids,ubot,Tamb) = data
		
		X=zeros(eltype(u), ng)
		mole_frac!(bnode,data,X,u)
		cf=heatcap_mix(Fluids, u[iT], X)		
		flux_enth_bot=ubot*u[ip]/(ph"R"*u[iT])*cf*(u[iT]-Tamb)
        f[iT]=flux_enth_bot
    end
end

function flux_radiation_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
		(;iT,lc_frit,lc_plate,Tplate) = data
		
		# irradiation exchange between porous frit (1) and Al bottom plate (2)
		# porous frit properties (1)
		eps1=lc_frit.eps; alpha1_IR=lc_frit.alpha_IR; rho1_IR=lc_frit.rho_IR; 	
		# Al bottom plate properties (2)
		eps2=lc_plate.eps; rho2_IR=lc_plate.rho_IR;
		
		#(;eps1,alpha1_IR,rho1_IR,rho2_IR,eps2,rho2_IR) = lower_chamber
		σ=ph"σ"
		flux_irrad = -eps1*σ*u[iT]^4 + alpha1_IR/(1-rho1_IR*rho2_IR)*(eps2*σ*Tplate^4+rho2_IR*eps1*σ*u[iT]^4)
        f[iT] = -flux_irrad
    end
end

function flux_conduction_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
        (;ng,iT,Tplate,lc_h)=data
        X=zeros(eltype(u), ng)    	

        # height of bottom chamber/ distance over which heat conduction occurs
        Tsurf = u[iT]

        # mean temperature
        Tm=0.5*(u[iT] + Tplate)
        mole_frac!(bnode,data,X,u)
        
        # thermal conductivity at Tm and outlet composition
        _,λf=dynvisc_thermcond_mix(data, Tm, X)

        f[iT] = λf*(Tsurf-Tplate)/lc_h

    end
end

function sidewalls(f,u,bnode,data)
    (;iT,α_w,Tamb) = data
    boundary_robin!(f,u,bnode;species=iT,region=Γ_side_right, factor=α_w, value=Tamb*α_w)
    boundary_robin!(f,u,bnode;species=iT,region=Γ_side_back, factor=α_w, value=Tamb*α_w)
end


function HeatFluxes()
    data=ModelData(
        isreactive=0
    )

    iT=data.iT
    
    sol_embed,grid,sys,data_embed = main(data=data);
    sol = sol_embed(1.0) # select final solution at the end of parameter embedding
    
    q_in = 4*integrate(sys,flux_in_aperture,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    q_in = sum(q_in)

    q_rad_top = 4*integrate(sys,flux_radiation_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    q_rad_top = sum(q_rad_top)

    q_cond_top = 4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    q_cond_top = sum(q_cond_top)

    q_sidewalls = 4*integrate(sys,sidewalls,sol; boundary=true)[iT,[Γ_side_right,Γ_side_back]]        
    q_sidewalls = sum(q_sidewalls)

    q_rad_bottom = 4*integrate(sys,flux_radiation_bottom,sol; boundary=true)[iT,Γ_bottom]
    q_cond_bottom = 4*integrate(sys,flux_conduction_bottom,sol; boundary=true)[iT,Γ_bottom]
    q_enthalpy_bottom = 4*integrate(sys,flux_enthalpy_bottom,sol; boundary=true)[iT,Γ_bottom]



    q_in,q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_enthalpy_bottom

end

function loss_plot(q_in,q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_enthalpy_bottom)

    #q_in,q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_enthalpy_bottom = HeatFluxes()
    x = ["Top Rad", "Top Cond", "Conv Sides", "Bottom Rad", "Bottom Cond", "ΔEnthalpy"]
    total_loss = sum([q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_enthalpy_bottom])
    y = [q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_enthalpy_bottom] ./ q_in
    
    #p=Plots.plot(x, y; st=:bar, texts=round.(y, sigdigits=2), legend=:none)
    p=Plots.plot(x, y; st=:bar, legend=:none)
    Plots.annotate!(p, x, y, round.(y, sigdigits=2), :bottom)
    Plots.annotate!(p, 2.0, 0.3, "Σ flows = "*string(round(total_loss/q_in*100.0,sigdigits=3))*" %", :left)

end

function PlotLosses(C=[1,10,25,50,75,100])
    q_in,q_top_abs,q_top_refl,q_top_rerad,q_top_convec,q_sides_conv,q_bot_rerad = HeatFluxes(C)
    p=plot(xguide="Irradiation / kW m⁻²", yguide="Loss contributions / -", legend=:outertopright, size=(400,250))
	areaplot!(C, [q_top_refl q_top_rerad q_top_convec q_sides_conv q_bot_rerad] ./ (q_in), fillalpha = [0.2 0.2], seriescolor = [1 2 3 4 5], label=["Refl Top" "Rerad Top" "Conv Top" "Conv Sides" "Rerad Bot"])
    savefig("img/loss_contributions.svg")

end

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

end
