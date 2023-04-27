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

function flux_radiation_top(f,u,bnode,data,T2,T3)
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

		#G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
        G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*T2^4+ϕ13*eps3*σ*T3^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))

                    

        G2_vis = rho2_vis*G1_bot_vis
        #G2_IR = eps2*σ*u[iT]^4 + rho2_IR*G1_bot_IR
        G2_IR = eps2*σ*T2^4 + rho2_IR*G1_bot_IR

        G3_vis = rho3_vis*G1_bot_vis
        #G3_IR = eps3*σ*u[iT]^4 + rho3_IR*G1_bot_IR
        G3_IR = eps3*σ*T3^4 + rho3_IR*G1_bot_IR

        # surface brigthness of quartz window outward facing surface (1) in vis & IR	
        G1_top_vis = rho1_vis*Glamp + tau1_vis*(ϕ12*G2_vis + ϕ13*G3_vis)
        G1_top_IR = tau1_IR*(ϕ12*G2_IR + ϕ13*G3_IR) + eps1*σ*Tglass^4

    
        iT=data.iT
        f[iT] = G1_top_vis + G1_top_IR
    end
end

# G1_bot, IR/vis : surface radiosities of glass underside 
function flux_window_underside(f,u,bnode,data,T2,T3)
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


        # surface radiosity of quartz window inwards / towards catalyst facing, vis & IR	
        G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))

		#G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
        G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*T2^4+ϕ13*eps3*σ*T3^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))

        f[iT] = G1_bot_vis + G1_bot_IR

        # G2_vis = rho2_vis*G1_bot_vis
        # #G2_IR = eps2*σ*u[iT]^4 + rho2_IR*G1_bot_IR
        # G2_IR = eps2*σ*T2^4 + rho2_IR*G1_bot_IR

        # G3_vis = rho3_vis*G1_bot_vis
        # #G3_IR = eps3*σ*u[iT]^4 + rho3_IR*G1_bot_IR
        # G3_IR = eps3*σ*T3^4 + rho3_IR*G1_bot_IR

        # # surface brigthness of quartz window outward facing surface (1) in vis & IR	
        # # G1_top_vis = rho1_vis*Glamp + tau1_vis*(ϕ12*G2_vis + ϕ13*G3_vis)
        # # G1_top_IR = tau1_IR*(ϕ12*G2_IR + ϕ13*G3_IR) + eps1*σ*Tglass^4

        # # transmitted radiosity through top window
        # G1_top_vis = tau1_vis * (ϕ12*G2_vis + ϕ13*G3_vis)
        # G1_top_IR = tau1_IR * (ϕ12*G2_IR + ϕ13*G3_IR)
    
        # iT=data.iT
        # f[iT] = G1_top_vis + G1_top_IR
    end
end

# G2_IR/vis : surface radiosities of catalyst layer 
function flux_catalyst_layer(f,u,bnode,data,T2,T3)
    if bnode.region==Γ_top_cat

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


        # surface radiosity of quartz window inwards / towards catalyst facing, vis & IR	
        G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))

        # this is not correct, stricly speaking: using the local values of temperature u[iT] for both, T2 and T3
        # potential remedy: after calculatioin, store averaged values from prev solution in data struct, for use as constants in next
        # simulation?
		G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
        #G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*T2^4+ϕ13*eps3*σ*T3^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))


        G2_vis = rho2_vis*G1_bot_vis

        # consider, wheter to uncomment this: use integrated average value or integrate with actual value
        G2_IR = eps2*σ*u[iT]^4 + rho2_IR*G1_bot_IR
        #G2_IR = eps2*σ*T2^4 + rho2_IR*G1_bot_IR

        f[iT] = G2_vis + G2_IR


    end
end

# G3_IR/vis : surface radiosities of frit layer 
function flux_frit(f,u,bnode,data,T2,T3)
    if bnode.region==Γ_top_frit

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


        # surface radiosity of quartz window inwards / towards catalyst facing, vis & IR	
        G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))

        # this is not correct, stricly speaking: using the local values of temperature u[iT] for both, T2 and T3
        # potential remedy: after calculatioin, store averaged values from prev solution in data struct, for use as constants in next
        # simulation?
		G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
        #G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*T2^4+ϕ13*eps3*σ*T3^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))


        G3_vis = rho3_vis*G1_bot_vis
        # consider, wheter to uncomment this: use integrated average value or integrate with actual value
        G3_IR = eps3*σ*u[iT]^4 + rho3_IR*G1_bot_IR
        #G3_IR = eps3*σ*T3^4 + rho3_IR*G1_bot_IR

        f[iT] = G3_vis + G3_IR

    end
end

function flux_transmission_top(f,u,bnode,data,T2,T3)
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


        # surface radiosity of quartz window inwards / towards catalyst facing, vis & IR	
        G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))

		#G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
        G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*T2^4+ϕ13*eps3*σ*T3^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))

                    

        G2_vis = rho2_vis*G1_bot_vis
        #G2_IR = eps2*σ*u[iT]^4 + rho2_IR*G1_bot_IR
        G2_IR = eps2*σ*T2^4 + rho2_IR*G1_bot_IR

        G3_vis = rho3_vis*G1_bot_vis
        #G3_IR = eps3*σ*u[iT]^4 + rho3_IR*G1_bot_IR
        G3_IR = eps3*σ*T3^4 + rho3_IR*G1_bot_IR

        # surface brigthness of quartz window outward facing surface (1) in vis & IR	
        # G1_top_vis = rho1_vis*Glamp + tau1_vis*(ϕ12*G2_vis + ϕ13*G3_vis)
        # G1_top_IR = tau1_IR*(ϕ12*G2_IR + ϕ13*G3_IR) + eps1*σ*Tglass^4

        # transmitted radiosity through top window
        G1_top_vis = tau1_vis * (ϕ12*G2_vis + ϕ13*G3_vis)
        G1_top_IR = tau1_IR * (ϕ12*G2_IR + ϕ13*G3_IR)
    
        iT=data.iT
        f[iT] = G1_top_vis + G1_top_IR
    end
end



# conductive heat flux through top chamber
function flux_conduction_top(f,u,bnode,data)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
        (;iT,Tglass,X0,uc_h,vf_uc_window_cat,vf_uc_window_frit)=data

		# mean temperature
        Tm=0.5*(u[iT] + Tglass)

        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        _,λf=dynvisc_thermcond_mix(data, Tm, X0)

        # positive flux in positive z-coordinate
        q_cond = -λf*(Tglass-u[iT])/uc_h

        f[iT] = q_cond

    end
end



function areas(sol,sys,grid,data)
	iT = data.iT
	function area(f,u,bnode,data)
		# repurpose temperature index to hold area information
		f[iT] = one(eltype(u))
	end

	integrate(sys,area,sol; boundary=true)[iT,:]
end



function T_avg(sol,sys,grid,data)

	iT = data.iT	
	function T_avg_(f,u,bnode,data)			
		f[iT] = u[iT]		
	end

	areas_=areas(sol,sys,grid,data)
	
	T_int=integrate(sys,T_avg_,sol; boundary=true)[iT,:]
	T_int./areas_	
end




function HeatFluxes()

    data=ModelData()
    
    sol_embed,grid,sys,data_embed = main(data=data);
    sol = sol_embed(1.0) # select final solution at the end of parameter embedding

    (;iT,Tglass,Tamb,α_nat_conv,uc_window)=data_embed
    σ=ph"σ"

    # T2 = Tfrit_avg
    # T3 = Tcat_avg
    Tfrit_avg,Tcat_avg = T_avg(sol,sys,grid,data_embed)[[Γ_top_frit,Γ_top_cat]]

    Atop = sum(areas(sol,sys,grid,data)[[Γ_top_cat,Γ_top_frit]])
    
    # heat flux contributions according to energybalance around upper chamber
    q_in = 4*integrate(sys,flux_in_aperture,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    q_in = sum(q_in)

    # reflection loss from top glass window
    q_RL = uc_window.rho_vis*q_in

    # transmission loss from inside of reactor through window
    q_TL = 4*integrate(
        sys,
        (f,u,bnode,data) -> flux_transmission_top(f,u,bnode,data_embed,Tcat_avg,Tfrit_avg),
        sol;
        boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    q_TL = sum(q_TL)

    # convection from glass top surface
    #q_conv_top = 4*sum(areas_[[Γ_top_cat,Γ_top_frit]])*α_nat_conv*(Tglass-Tamb)
    q_conv_top = 4*Atop*α_nat_conv*(Tglass-Tamb)

    q_emit_top = 4*Atop*uc_window.eps*σ*Tglass^4

    # q_rad_top = 4*integrate(
    #     sys,
    #     (f,u,bnode,data) -> flux_radiation_top(f,u,bnode,data_embed,Tcat_avg,Tfrit_avg),sol; boundary=true
    #     )[iT,[Γ_top_cat,Γ_top_frit]]
    # q_rad_top = sum(q_rad_top)
    
    
    # conduction through upper chamber gas layer
    q_cond = 4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    q_cond = sum(q_cond)


    # radiosity from underside of glass window
    G1_bot_vis_IR = 4*integrate(
        sys,
        (f,u,bnode,data) -> flux_window_underside(f,u,bnode,data_embed,Tcat_avg,Tfrit_avg),
        sol;
        boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    G1_bot_vis_IR = sum(G1_bot_vis_IR)
    

    # radiosity from catalyst layer
    G2_vis_IR = 4*integrate(
        sys,
        (f,u,bnode,data) -> flux_catalyst_layer(f,u,bnode,data_embed,Tcat_avg,Tfrit_avg),
        sol;
        boundary=true)[iT,Γ_top_cat]

    # radiosity from frit
    G3_vis_IR = 4*integrate(
        sys,
        (f,u,bnode,data) -> flux_frit(f,u,bnode,data_embed,Tcat_avg,Tfrit_avg),
        sol;
         boundary=true)[iT,Γ_top_frit]
    
    
    (
        q_in=q_in,
        q_RL=q_RL,
        q_TL=q_TL,
        q_conv_top=q_conv_top,
        q_emit_top=q_emit_top,
        q_cond=q_cond,
        G1_bot_vis_IR=G1_bot_vis_IR,
        G2_vis_IR=G2_vis_IR,
        G3_vis_IR=G3_vis_IR
    )
end

function loss_plot(losses)
    (;q_in,q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_delta_enthalpy,q_delta_enthalpy_reaction,q_rad_bottom_2,q_conv_bottom,q_conv_top) = losses

    q_delta_enthalpy_sensible = q_delta_enthalpy-q_delta_enthalpy_reaction

    lbs = ["Top Rad", "Top Cond", "Conv Sides", "Bottom Rad", "Bottom Cond", "HSensible","HReaction"]
    #lbs = ["Top Rad", "Top Cond", "Top Conv", "Conv Sides", "Bottom Rad", "Bottom Cond", "Bottom Rad 2", "Bottom Conv", "HSensible", "HReaction"]
    x = [string(x) for x in 1:length(lbs)]
    lbsx = [x*": "* l for (l,x) in zip(lbs,x)]
    total_loss = sum([q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_delta_enthalpy])
    y = [q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_delta_enthalpy_sensible,q_delta_enthalpy_reaction] / q_in * 100.0
    #y = [q_rad_top,q_cond_top,q_conv_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_rad_bottom_2,q_conv_bottom,q_delta_enthalpy_sensible,q_delta_enthalpy_reaction] / q_in * 100.0
    
    p=Plots.plot(size=(400,300), yguide="Heat flow contribution / %")
    Plots.plot!(p, permutedims(x), permutedims(y); st=:bar, label=permutedims(lbsx))    
    Plots.annotate!(p, x, y, round.(y, sigdigits=2), :bottom)
    Plots.annotate!(p, 2.0, 30.0, "Σ flows = "*string(round(total_loss/q_in*100.0,sigdigits=3))*" %", :left)
    

end
 function EnergyBalance()
    fluxes=Example2_EnergyBalance.HeatFluxes()
    (;q_in,q_RL,q_TL,q_conv_top,q_emit_top,q_cond,G1_bot_vis_IR,G2_vis_IR,G3_vis_IR) = fluxes

    # energy balance: should sum to 0
    balance = q_in-q_RL-q_TL-q_conv_top-q_emit_top+q_cond-G1_bot_vis_IR+G2_vis_IR+G3_vis_IR

    # relative error of energy balance
    balance/q_in
    
 end


end
