module Example2_EnergyBalance

using DataFrames, CSV


include("../notebooks/PorousCatalystHot3DTopFlowIrrExchange.jl")


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

# G0_bot, IR/vis : surface radiosities of glass underside
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


# G0_bot, IR/vis : surface radiosities of glass underside
function flux_window_underside(f,u,bnode,data)
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

		G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
        #G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*T2^4+ϕ13*eps3*σ*T3^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))

        f[iT] = G1_bot_vis + G1_bot_IR

    end
end


# G2_IR/vis : surface radiosities of catalyst layer 
function flux_catalyst_layer(f,u,bnode,data)
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

# G3_IR/vis : surface radiosities of frit 
function flux_frit(f,u,bnode,data)
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

# tau0*G1_IR/vis : transmitted fraction through window of surface radiosityies of catalyst layer 
function flux_transmission_top_catalyst_layer(f,u,bnode,data)
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

        f[iT] = tau1_vis*G2_vis + tau1_IR*G2_IR


    end
end

# tau0*G2_IR/vis : transmitted fraction through window of surface radiosityies of frit
function flux_transmission_top_frit(f,u,bnode,data)
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

        f[iT] = tau1_vis*G3_vis + tau1_IR*G3_IR

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
        (;iT,Tglass,X0,uc_h,)=data

		# mean temperature
        Tm=0.5*(u[iT] + Tglass)

        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        _,λf=dynvisc_thermcond_mix(data, Tm, X0)

        # positive flux in positive z-coordinate
        q_cond = -λf*(Tglass-u[iT])/uc_h

        f[iT] = q_cond

    end
end

# surface radiosity of underside of frit, facing towards the bottom plate
# solely in IR, consists of emission and reflection
function flux_radiation_frit_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
		(;iT,lc_frit,lc_plate,Tplate) = data
		
		# irradiation exchange between porous frit (1) and Al bottom plate (2)
		# porous frit properties (1)
		eps1=lc_frit.eps;  rho1_IR=lc_frit.rho_IR; 	
		# Al bottom plate properties (2)
		eps2=lc_plate.eps; rho2_IR=lc_plate.rho_IR;
		
		#(;eps1,alpha1_IR,rho1_IR,rho2_IR,eps2,rho2_IR) = lower_chamber
		σ=ph"σ"

        G1_IR = (eps1*σ*u[iT]^4 + rho1_IR*eps2*σ*Tplate^4)/(1-rho1_IR*rho2_IR)

        f[iT] = G1_IR
    end
end

function flux_radiation_plate_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
		(;iT,lc_frit,lc_plate,Tplate) = data
		
		# irradiation exchange between porous frit (1) and Al bottom plate (2)
		# porous frit properties (1)
		eps1=lc_frit.eps;  rho1_IR=lc_frit.rho_IR; 	
		# Al bottom plate properties (2)
		eps2=lc_plate.eps; rho2_IR=lc_plate.rho_IR;
		
		#(;eps1,alpha1_IR,rho1_IR,rho2_IR,eps2,rho2_IR) = lower_chamber
		σ=ph"σ"

        G2_IR = (eps2*σ*Tplate^4 + rho2_IR*eps1*σ*u[iT]^4)/(1-rho1_IR*rho2_IR)

        f[iT] = G2_IR
    end
end



# conductive heat flux through bottom chamber
function flux_conduction_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
        (;ng,iT,Tplate,lc_h)=data

        X=zeros(eltype(u), ng)
		mole_frac!(bnode,data,X,u)

		# mean temperature
        Tm=0.5*(u[iT] + Tplate)
        # thermal conductivity at Tm and outlet composition X
        _,λf=dynvisc_thermcond_mix(data, Tm, X)

        # positive flux in negative z coord. -> pointing towards bottom plate
        q_cond = -λf*(Tplate-u[iT])/lc_h 

        f[iT] = q_cond

    end
end

function flux_abs_top_catalyst_layer(f,u,bnode,data)
    if bnode.region==Γ_top_cat
        (;iT,Tglass,Glamp,uc_window,uc_cat,uc_frit,vf_uc_window_cat,vf_uc_window_frit)=data
        σ=ph"σ"

        # irradiation exchange between quartz window (1), cat surface (2) & frit surface (3)
        # window properties (1)
        tau1_vis=uc_window.tau_vis
        rho1_vis=uc_window.rho_vis
        tau1_IR=uc_window.tau_IR
        rho1_IR=uc_window.rho_IR
        alpha1_vis=uc_window.alpha_vis
        alpha1_IR=uc_window.alpha_IR
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
        
        G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))

        G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))	        		

        G2_vis = rho2_vis*G1_bot_vis		
        
        G2_IR = eps2*σ*u[iT]^4 + rho2_IR*G1_bot_IR

        f[iT] = alpha1_vis*G2_vis + alpha1_IR*G2_IR
    end
end

function flux_abs_top_frit(f,u,bnode,data)
    if bnode.region==Γ_top_cat

        (;iT,Tglass,Glamp,uc_window,uc_cat,uc_frit,vf_uc_window_cat,vf_uc_window_frit)=data
        σ=ph"σ"

        # irradiation exchange between quartz window (1), cat surface (2) & frit surface (3)
        # window properties (1)
        tau1_vis=uc_window.tau_vis
        rho1_vis=uc_window.rho_vis
        tau1_IR=uc_window.tau_IR
        rho1_IR=uc_window.rho_IR
        alpha1_vis=uc_window.alpha_vis
        alpha1_IR=uc_window.alpha_IR
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
        
        G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))
        
        G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))

        G3_vis = rho3_vis*G1_bot_vis

        G3_IR = eps3*σ*u[iT]^4 + rho3_IR*G1_bot_IR

        f[iT] = alpha1_vis*G3_vis + alpha1_IR*G3_IR
    end
end


function flux_enthalpy_top(f,u,bnode,data)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
        
		(;ng,iT,Fluids,u0,X0,pn,Tn,Tamb) = data
		
        flux_enth_top=0.0
        for i=1:ng
             # sign convention: outward pointing surface normal as positive
             # MolarFlow = -u0*X0[i]*pn/(ph"R"*Tn)
             MolarFlow = u0*X0[i]*pn/(ph"R"*Tn)
             flux_enth_top += MolarFlow * enthalpy_gas(Fluids[i], Tamb)
        end

        #hmix=enthalpy_mix(Fluids, Tamb, X0)
		#flux_enth_top=u0*pn/(ph"R"*Tn)*hmix

        f[iT]=flux_enth_top

    end
end

function flux_enthalpy_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
		(;ng,iT,ip,Fluids,ubot,Tamb) = data
		
        X=zeros(eltype(u), ng)
        mole_frac!(bnode,data,X,u)

        hmix=enthalpy_mix(Fluids, u[iT], X)
		flux_enth_bot=ubot*u[ip]/(ph"R"*u[iT])*hmix

		#flux_enth_bot=ubot*u[ip]/(ph"R"*u[iT])*cf*(u[iT]-Tamb)
        f[iT]=flux_enth_bot
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


function runSim()
    #sol_embed,grid,sys,data_embed = main();
    #sol = sol_embed(1.0) # select final solution at the end of parameter embedding

    sol_,grid,sys,data_embed=main();
	if sol_ isa VoronoiFVM.TransientSolution
		sol = copy(sol_(sol_.t[end]))
	else
		sol = copy(sol_)
	end

    sol,grid,sys,data_embed
end

function HeatFluxes_EB_I(sol,grid,sys,data)


    (;iT)=data
    σ=ph"σ"

    Qtop=4*integrate(sys,top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qtop=sum(Qtop)

    Qbot=4*integrate(sys,bottom,sol; boundary=true)[iT,Γ_bottom]

    Hin=4*integrate(sys,flux_enthalpy_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Hin=sum(Hin)

    Qcond_10=4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_cat]
    Qcond_20=4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_frit]

    QG_01=4*integrate(sys,flux_window_underside,sol; boundary=true)[iT,Γ_top_cat]
    QG_02=4*integrate(sys,flux_window_underside,sol; boundary=true)[iT,Γ_top_frit]

    QG_10=4*integrate(sys,flux_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    QG_20=4*integrate(sys,flux_frit,sol; boundary=true)[iT,Γ_top_frit]

    # manually calc Qtop from indiv contributions
    Qtop_calc= Hin +QG_01 +QG_02 -QG_10 -QG_20 -Qcond_10 -Qcond_20
        
    Hout=4*integrate(sys,flux_enthalpy_bottom,sol; boundary=true)[iT,Γ_bottom]

    Qcond_34=4*integrate(sys,flux_conduction_bottom,sol; boundary=true)[iT,Γ_bottom]

    QG_34=4*integrate(sys,flux_radiation_frit_bottom,sol; boundary=true)[iT,Γ_bottom]
    QG_43=4*integrate(sys,flux_radiation_plate_bottom,sol; boundary=true)[iT,Γ_bottom]
    
    Qsides=4*integrate(sys,side,sol; boundary=true)[iT,[Γ_side_right,Γ_side_back]]        
    Qsides=sum(Qsides)

    Qbot_calc= -Hout -Qcond_34 -QG_34 +QG_43 
    
    #Qbot,Qbot_calc
    
    #Qtop,Qbot,Qsides, sum(Qbot+Qtop+Qsides)
    Qtop_calc,Qbot_calc,-Qsides, sum([Qtop_calc,Qbot_calc,-Qsides])
    #Qtop_calc_no_enth =QG_01 +QG_02 -QG_10 -QG_20 -Qcond_10 -Qcond_20

end

function HeatFluxes_EB_II(sol,grid,sys,data)


    (;iT,Tglass,Glamp,Tamb,α_nat_conv,uc_window)=data
    σ=ph"σ"

    Atop = sum(areas(sol,sys,grid,data)[[Γ_top_cat,Γ_top_frit]])

    Qconv0 = 4*Atop*α_nat_conv*(Tglass-Tamb)
    Qemit0 = 4*Atop*uc_window.eps*σ*Tglass^4
    Qin = 4*Atop*Glamp
    Qrefl0 = uc_window.rho_vis*Qin

    Qtrans_10=4*integrate(sys,flux_transmission_top_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    Qtrans_20=4*integrate(sys,flux_transmission_top_frit,sol; boundary=true)[iT,Γ_top_frit]

    ### correct, proven in EB I ###
    Qcond_10=4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_cat]
    Qcond_20=4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_frit]

    QG_01=4*integrate(sys,flux_window_underside,sol; boundary=true)[iT,Γ_top_cat]
    QG_02=4*integrate(sys,flux_window_underside,sol; boundary=true)[iT,Γ_top_frit]

    QG_10=4*integrate(sys,flux_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    QG_20=4*integrate(sys,flux_frit,sol; boundary=true)[iT,Γ_top_frit]
    ### correct, proven in EB I ###

    # calc energy balance II
    EB_bot = -QG_01 +Qcond_10 +QG_10 -QG_02 +Qcond_20 +QG_20
    EB_top = -Qconv0 -Qemit0 +Qin -Qrefl0 -Qtrans_10 -Qtrans_20
    EB = EB_bot + EB_top
    #EB = -Qconv0 -Qemit0 +Qin -Qrefl0 -Qtrans_10 -Qtrans_20 -QG_01 +Qcond_10 +QG_10 -QG_02 +Qcond_20 +QG_20
    
end


function HeatFluxes_EB_III(sol,grid,sys,data)



    (;iT,Tamb,α_nat_conv,Tplate,lc_plate)=data
    σ=ph"σ"

    Abot = areas(sol,sys,grid,data)[Γ_bottom]

    Qconv4 = 4*Abot*α_nat_conv*(Tplate-Tamb)
	Qemit4 = 4*Abot*lc_plate.eps*σ*Tplate^4		

    ### correct, proven in EB I ###
    Qcond_34=4*integrate(sys,flux_conduction_bottom,sol; boundary=true)[iT,Γ_bottom]

    QG_34=4*integrate(sys,flux_radiation_frit_bottom,sol; boundary=true)[iT,Γ_bottom]
    QG_43=4*integrate(sys,flux_radiation_plate_bottom,sol; boundary=true)[iT,Γ_bottom]
    ### correct, proven in EB I ###

    # calc energy balance III
    EB_top = Qcond_34 -QG_43 +QG_34
    EB_bot = -Qconv4 -Qemit4
    
    EB = EB_bot + EB_top
    
    
end

function HeatFluxes_EB_IV_inner(sol,grid,sys,data)
    (;iT,α_nat_conv,Tglass,Tamb,uc_window)=data

    σ=ph"σ"    
        
    Qabs_10 = 4*integrate(sys,flux_abs_top_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    Qabs_20 = 4*integrate(sys,flux_abs_top_frit,sol; boundary=true)[iT,Γ_top_frit]

    Qcond_10=4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_cat]
    Qcond_20=4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_frit]

    Atop = sum(areas(sol,sys,grid,data)[[Γ_top_cat,Γ_top_frit]])

    Qconv0 = 4*Atop*α_nat_conv*(Tglass-Tamb)
    Qemit0 = 4*Atop*uc_window.eps*σ*Tglass^4		

    EB = -Qconv0 +Qcond_10 +Qcond_20 - 2*Qemit0 +Qabs_10 +Qabs_20
        

end

function FluxesSymmetryBC(sol,sys,data)
    (;iT)=data
    tff=TestFunctionFactory(sys)
    Γ_where_T_equal_1=[Γ_side_front,Γ_side_left] # front & left should be symmetry bcs.
    Γ_where_T_equal_0=[Γ_side_right,Γ_side_back,Γ_bottom,Γ_top_cat,Γ_top_frit]
    Tf=testfunction(tff,Γ_where_T_equal_0,Γ_where_T_equal_1)
    integrate(sys,Tf,sol)[iT]
end

end
